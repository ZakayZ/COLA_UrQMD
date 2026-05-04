#!/usr/bin/env python3
"""Sweep edge ON/OFF decision thresholds τ: activate edge iff sigmoid(logit) > τ.

Runs the same pipeline as ``run_reinforce_unsup_knn6.py`` (supervised cache → optional RL with
tuned defaults), then prints metrics vs τ for:

  • Policy after supervised warm-start only (``--skip-rl``), or
  • After REINFORCE with the same RL hyperparameters as the tuned run script defaults.

RL note (threshold mode): ``collect_rollout`` in ``training/reinforce.py`` builds masks with **τ_train = 0.5**
fixed. Training rewards align with that discrete mask; sweeping τ here is **evaluation-only** and
shows sensitivity of physics metrics if you deployed a different cutoff than 0.5.

Bernoulli RL uses sampled binary masks; there is no single τ — see printed note."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import cast

_CLUSTERING_ROOT = Path(__file__).resolve().parent.parent
_TESTS_DIR = Path(__file__).resolve().parent
for _root in (_CLUSTERING_ROOT, _TESTS_DIR):
    _sr = str(_root)
    if _sr not in sys.path:
        sys.path.insert(0, _sr)

import numpy as np
import torch
import torch.optim as optim
from torch.optim.lr_scheduler import CosineAnnealingLR

from models import (
    AffinityGraphConfig,
    AffinityGraphEnv,
    GAT_NODE_IN_DIM,
    GATAffinityPolicy,
    init_policy_all_edges_off,
    load_valid_events_from_pkl,
)
from training.utils import RLActionMode, make_event_sampler
from training.reinforce import train_reinforce
from training.supervised import train_supervised_edges

from run_reinforce_unsup_knn6 import (
    default_supervised_cache_path,
    save_supervised_cache,
    supervised_cache_meta,
    try_load_supervised_cache,
)

# Keep in sync with ``run_reinforce_unsup_knn6.parse_args`` defaults (post-tuning).
RL_BEST = dict(
    updates=12,
    episodes_per_update=24,
    lr=1.5e-6,
    ent_coef=0.001,
    max_grad_norm=1.2,
    policy_coef=0.08,
    anchor_bce_coef=0.006,
    anchor_pos_weight=15.0,
    anchor_focal_gamma=2.0,
    rl_action_mode="threshold",
)


def parse_threshold_list(s: str) -> list[float]:
    parts = [p.strip() for p in s.replace(";", ",").split(",") if p.strip()]
    out: list[float] = []
    for p in parts:
        out.append(float(p))
    return sorted(set(out))


def stats_at_threshold(
    policy: GATAffinityPolicy,
    env: AffinityGraphEnv,
    sampler,
    n_rollouts: int,
    tau: float,
) -> tuple[float, float, float, float]:
    """Mean return (MeV), mean L_pol (MeV), mean gap (MeV), mean fraction of edges ON."""

    policy.eval()
    rets: list[float] = []
    gaps: list[float] = []
    lp: list[float] = []
    frac_on: list[float] = []
    for _ in range(max(int(n_rollouts), 1)):
        pos, mom, isp = sampler()
        obs = env.reset(pos, mom, isp)
        logits = policy(obs)
        prob = torch.sigmoid(logits)
        a = (prob > float(tau)).float()
        loss, _ = env.physics_for_edge_mask(a)
        lb = float(env._baseline_loss)
        Lp = float(loss)
        lp.append(Lp)
        gaps.append(Lp - lb)
        rets.append(float(-Lp))
        frac_on.append(float(a.mean().item()))
    policy.train()
    return (
        float(np.mean(rets)),
        float(np.mean(lp)),
        float(np.mean(gaps)),
        float(np.mean(frac_on)),
    )


def print_table(title: str, taus: list[float], rows: list[tuple[float, float, float, float]]) -> None:
    print(f"\n=== {title} ===")
    print(f"{'tau':>6}  {'mean_G':>10}  {'mean_L_pol':>12}  {'mean_gap':>12}  {'frac_edges_on':>14}")
    print("-" * 62)
    for tau, row in zip(taus, rows):
        g, lpol, gap, fo = row
        print(f"{tau:6.3f}  {g:10.2f}  {lpol:12.2f}  {gap:12.2f}  {fo:14.4f}")


def parse_args() -> argparse.Namespace:
    d_ds = _CLUSTERING_ROOT / "datasets" / "urqmd_nucleons_1k" / "dataset.pkl"
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", type=Path, default=d_ds)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument(
        "--thresholds",
        type=str,
        default=(
            "0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95"
        ),
        help="Comma-separated τ values in (0,1). Edge ON iff sigmoid(logit) > τ.",
    )
    p.add_argument("--eval-rollouts", type=int, default=64)
    p.add_argument(
        "--skip-rl",
        action="store_true",
        help="Only supervised weights; skip REINFORCE (still sweep τ).",
    )
    p.add_argument("--cold-all-off", action="store_true")
    p.add_argument("--sup-steps", type=int, default=80)
    p.add_argument("--sup-events-per-step", type=int, default=8)
    p.add_argument("--sup-lr", type=float, default=1e-3)
    p.add_argument("--sup-pos-weight", type=float, default=20.0)
    p.add_argument("--sup-focal-gamma", type=float, default=2.0)
    p.add_argument("--no-supervised-cache", action="store_true")
    p.add_argument("--supervised-cache", type=Path, default=None)
    p.add_argument("--logit-bias", type=float, default=-10.0)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    taus = parse_threshold_list(args.thresholds)
    if not taus:
        print("No thresholds parsed.", file=sys.stderr)
        return 1
    for tau in taus:
        if not 0.0 < tau < 1.0:
            print(f"Threshold τ={tau} should lie strictly between 0 and 1.", file=sys.stderr)
            return 1

    if not args.dataset.exists():
        print(f"Dataset not found: {args.dataset}", file=sys.stderr)
        return 1

    events = load_valid_events_from_pkl(args.dataset)
    if not events:
        print(f"No events in {args.dataset}", file=sys.stderr)
        return 1

    rng = np.random.default_rng(args.seed)
    cfg = AffinityGraphConfig(k_nn=6, graph_kind="knn")
    env = AffinityGraphEnv(cfg)
    policy = GATAffinityPolicy(
        in_dim=GAT_NODE_IN_DIM,
        hidden=64,
        n_heads=4,
        n_gat_layers=2,
        running_norm=False,
    )
    sampler = make_event_sampler(events=events, rng=rng, fallback_urqmd=None)

    print(
        """
Edge mask rule (deterministic eval): for each edge, ON iff σ(logit) > τ.
  • Larger τ → fewer edges ON (stricter, sparser graphs).
  • Smaller τ → more edges ON (denser graphs).

RL ``threshold`` mode (``training.reinforce.collect_rollout``): training mask uses τ_train = 0.5 only.
  Policy gradient uses Bernoulli log-prob of that threshold action. Tables below vary τ at
  **evaluation** time to show cutoff sensitivity; they do not change how the policy was trained.

RL ``bernoulli`` mode: actions are sampled; no τ — interpret via stochastic rollout_stats instead.
"""
    )

    if args.cold_all_off:
        init_policy_all_edges_off(policy, logit_bias=args.logit_bias)
        print("=== Mode: cold init (all edges off) ===")
    else:
        ns = argparse.Namespace(
            seed=args.seed,
            sup_steps=args.sup_steps,
            sup_events_per_step=args.sup_events_per_step,
            sup_lr=args.sup_lr,
            sup_pos_weight=args.sup_pos_weight,
            sup_focal_gamma=args.sup_focal_gamma,
        )
        meta = supervised_cache_meta(ns, args.dataset)
        cache_path: Path | None = None
        if not args.no_supervised_cache:
            cache_path = (
                args.supervised_cache
                if args.supervised_cache is not None
                else default_supervised_cache_path(args.dataset, _CLUSTERING_ROOT)
            )
        loaded = False
        if cache_path is not None:
            loaded = try_load_supervised_cache(policy, cache_path, meta)
        if loaded:
            print("=== Supervised warm-start: from cache ===")
        else:
            print("=== Supervised warm-start (training) ===")
            sup_opt = optim.Adam(policy.parameters(), lr=args.sup_lr)
            train_supervised_edges(
                policy,
                env,
                sampler,
                steps=max(1, args.sup_steps),
                events_per_step=max(1, args.sup_events_per_step),
                lr=args.sup_lr,
                focal_gamma=args.sup_focal_gamma,
                pos_weight=args.sup_pos_weight,
                optimizer=sup_opt,
                max_grad_norm=0.5,
            )
            if cache_path is not None:
                save_supervised_cache(policy, cache_path, meta)

    rows_sup = [
        stats_at_threshold(policy, env, sampler, args.eval_rollouts, tau) for tau in taus
    ]
    print_table(
        f"After supervised warm-start (eval n={args.eval_rollouts})",
        taus,
        rows_sup,
    )

    if args.skip_rl:
        print("\n(--skip-rl: no REINFORCE.)")
        return 0

    rl_opt = optim.Adam(policy.parameters(), lr=RL_BEST["lr"])
    eta_min = min(1e-5, max(1e-8, float(RL_BEST["lr"]) * 0.05))
    sched = CosineAnnealingLR(
        rl_opt, T_max=max(1, RL_BEST["updates"]), eta_min=eta_min
    )
    print("\n=== REINFORCE (tuned defaults from run_reinforce_unsup_knn6) ===")
    print(RL_BEST)
    train_reinforce(
        policy=policy,
        env=env,
        event_sampler=sampler,
        optimizer=rl_opt,
        lr_scheduler=sched,
        n_updates=max(1, RL_BEST["updates"]),
        episodes_per_update=max(1, RL_BEST["episodes_per_update"]),
        ent_coef=RL_BEST["ent_coef"],
        max_grad_norm=float(RL_BEST["max_grad_norm"]),
        policy_coef=float(RL_BEST["policy_coef"]),
        rl_action_mode=cast(RLActionMode, RL_BEST["rl_action_mode"]),
        anchor_bce_coef=float(RL_BEST["anchor_bce_coef"]),
        anchor_pos_weight=float(RL_BEST["anchor_pos_weight"]),
        anchor_focal_gamma=float(RL_BEST["anchor_focal_gamma"]),
    )

    rows_rl = [
        stats_at_threshold(policy, env, sampler, args.eval_rollouts, tau) for tau in taus
    ]
    print_table(
        f"After REINFORCE + same hyperparameters (eval n={args.eval_rollouts})",
        taus,
        rows_rl,
    )

    i05 = min(range(len(taus)), key=lambda i: abs(taus[i] - 0.5))
    g0, _, _, _ = rows_sup[i05]
    g1, _, _, _ = rows_rl[i05]
    print(
        f"\nAt τ=0.5 (training threshold for RL threshold-mode): "
        f"mean return supervised ≈ {g0:.2f} MeV → after RL ≈ {g1:.2f} MeV "
        f"(Δ {g1 - g0:+.2f} MeV)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
