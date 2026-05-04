#!/usr/bin/env python3
"""kNN k=6: optional supervised warm-start, then REINFORCE; report gains vs post-supervised baseline.

Supervised weights are cached under clustering/.cache/ (see --supervised-cache) so repeated runs
skip edge-BCE training when dataset and supervised hyperparameters match.

RL anchor BCE (train_reinforce) matches logits to the current graph via a fresh env.reset per step.
After warm-start, deterministic quality is often more stable with a modest --anchor-bce-coef,
--policy-coef below 1, low --ent-coef, and a fixed --seed so the cache snapshot stays comparable.

Default CLI hyperparameters are tuned for a short post–warm-start phase (~12 updates) that has
passed the built-in ≥8 MeV deterministic gain check on urqmd_nucleons_1k (seed 0)."""

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
from training.reinforce import collect_rollout, train_reinforce
from training.supervised import train_supervised_edges

SUP_CACHE_FORMAT = "supervised_knn6_v1"


def default_supervised_cache_path(dataset: Path, here: Path) -> Path:
    ds = dataset.resolve()
    try:
        rel = ds.relative_to(here)
        tag = "__".join(rel.with_suffix("").parts)
    except ValueError:
        tag = f"{ds.parent.name}_{ds.stem}"
    return here / ".cache" / f"sup_knn6_{tag}.pt"


def supervised_cache_meta(args: argparse.Namespace, dataset: Path) -> dict:
    return {
        "format": SUP_CACHE_FORMAT,
        "dataset": str(dataset.resolve()),
        "seed": int(args.seed),
        "sup_steps": int(args.sup_steps),
        "sup_events_per_step": int(args.sup_events_per_step),
        "sup_lr": float(args.sup_lr),
        "sup_pos_weight": float(args.sup_pos_weight),
        "sup_focal_gamma": float(args.sup_focal_gamma),
        "policy": {
            "in_dim": GAT_NODE_IN_DIM,
            "hidden": 64,
            "n_heads": 4,
            "n_gat_layers": 2,
            "running_norm": False,
            "k_nn": 6,
        },
    }


def meta_matches_checkpoint(meta: dict, ckpt_meta: dict) -> bool:
    if ckpt_meta.get("format") != SUP_CACHE_FORMAT:
        return False
    keys = (
        "dataset",
        "seed",
        "sup_steps",
        "sup_events_per_step",
        "sup_lr",
        "sup_pos_weight",
        "sup_focal_gamma",
        "policy",
    )
    return all(ckpt_meta.get(k) == meta.get(k) for k in keys)


def try_load_supervised_cache(
    policy: GATAffinityPolicy,
    path: Path,
    meta: dict,
) -> bool:
    if not path.is_file():
        return False
    raw = torch.load(path, map_location="cpu")
    if not isinstance(raw, dict) or "state_dict" not in raw:
        print(f"Supervised cache invalid (missing state_dict): {path}", file=sys.stderr)
        return False
    ckpt_meta = raw.get("meta") or {}
    if not meta_matches_checkpoint(meta, ckpt_meta):
        print(
            f"Supervised cache metadata mismatch — will retrain and overwrite: {path}",
            file=sys.stderr,
        )
        return False
    policy.load_state_dict(raw["state_dict"])
    print(f"Loaded supervised checkpoint from {path}")
    return True


def save_supervised_cache(policy: GATAffinityPolicy, path: Path, meta: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "meta": meta,
        "state_dict": {k: v.detach().cpu() for k, v in policy.state_dict().items()},
    }
    torch.save(payload, path)
    print(f"Saved supervised checkpoint to {path}")


def parse_args() -> argparse.Namespace:
    default_dataset = _CLUSTERING_ROOT / "datasets" / "urqmd_nucleons_1k" / "dataset.pkl"
    p = argparse.ArgumentParser(
        description="REINFORCE after optional supervised warm-start (kNN k=6)."
    )
    p.add_argument("--dataset", type=Path, default=default_dataset)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--cold-all-off", action="store_true", help="Skip supervised; init all edges off.")
    p.add_argument("--sup-steps", type=int, default=80)
    p.add_argument("--sup-events-per-step", type=int, default=8)
    p.add_argument("--sup-lr", type=float, default=1e-3)
    p.add_argument("--sup-pos-weight", type=float, default=20.0)
    p.add_argument("--sup-focal-gamma", type=float, default=2.0)
    p.add_argument(
        "--updates",
        type=int,
        default=12,
        help=(
            "REINFORCE optimizer steps (cosine LR to eta_min). "
            "Short schedules with small LR are easier on deterministic eval after warm-start."
        ),
    )
    p.add_argument(
        "--episodes-per-update",
        type=int,
        default=24,
        help="Rollouts per update (more → lower-variance batch baseline).",
    )
    p.add_argument(
        "--lr",
        type=float,
        default=1.5e-6,
        help="Adam LR for REINFORCE (small after supervised warm-start).",
    )
    p.add_argument(
        "--ent-coef",
        type=float,
        default=0.001,
        help="Entropy bonus coefficient (keep small to preserve sharp supervised logits).",
    )
    p.add_argument(
        "--max-grad-norm",
        type=float,
        default=1.2,
        help="Global grad clip for REINFORCE (avoid clipping every step at 0.5).",
    )
    p.add_argument(
        "--rl-action-mode",
        type=str,
        choices=("bernoulli", "threshold"),
        default="threshold",
        help=(
            "threshold: REINFORCE uses σ(logit)>0.5 masks (matches deterministic eval). "
            "bernoulli: sampled masks (classic REINFORCE)."
        ),
    )
    p.add_argument(
        "--anchor-bce-coef",
        type=float,
        default=0.006,
        help=(
            "Mix supervised BCE toward spatial–momentum baseline edge targets each RL update "
            "(0 disables). Too large a value fights the physics reward; ~0.006–0.01 is a typical stabilizer."
        ),
    )
    p.add_argument(
        "--anchor-pos-weight",
        type=float,
        default=15.0,
        help="Positive-class weight for RL anchor BCE (same role as supervised pos_weight).",
    )
    p.add_argument(
        "--anchor-focal-gamma",
        type=float,
        default=2.0,
        help="Focal gamma for RL anchor BCE.",
    )
    p.add_argument(
        "--policy-coef",
        type=float,
        default=0.08,
        help=(
            "Weight on REINFORCE policy-gradient term (0 = anchor-BCE-only updates; "
            "pair with --ent-coef 0 to avoid entropy-only drift). Values ~0.05–0.15 often "
            "beat pure REINFORCE after supervised warm-start."
        ),
    )
    p.add_argument("--logit-bias", type=float, default=-10.0, help="Only for --cold-all-off.")
    p.add_argument("--eval-rollouts", type=int, default=64)
    p.add_argument("--tail", type=int, default=30)
    p.add_argument(
        "--min-gain-mev",
        type=float,
        default=8.0,
        help=(
            "Exit 0 if deterministic eval (sigmoid>0.5) improves vs post-supervised "
            "snapshot by at least this MeV on mean return and/or gap."
        ),
    )
    p.add_argument(
        "--supervised-cache",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Load/save supervised weights here. Default: clustering/.cache/sup_knn6_<dataset>.pt "
            "Set --no-supervised-cache to always train from scratch."
        ),
    )
    p.add_argument(
        "--no-supervised-cache",
        action="store_true",
        help="Disable supervised checkpoint load/save.",
    )
    p.add_argument(
        "--diag-jsonl",
        type=Path,
        default=None,
        metavar="PATH",
        help=(
            "Append per-update REINFORCE diagnostics (logits, capture stats, grads, "
            "stoch vs det same-event rewards, det eval sweep) as JSON lines."
        ),
    )
    p.add_argument(
        "--diag-every",
        type=int,
        default=1,
        help="Log diagnostics every N updates (with --diag-jsonl).",
    )
    p.add_argument(
        "--diag-det-rollouts",
        type=int,
        default=32,
        help="Deterministic eval episodes per diagnostic row (with --diag-jsonl).",
    )
    return p.parse_args()


def rollout_stats(
    policy: GATAffinityPolicy,
    env: AffinityGraphEnv,
    sampler,
    n: int,
    *,
    deterministic: bool = True,
) -> tuple[float, float, float]:
    """Mean return (MeV), mean L_pol (MeV), mean gap L_pol - L_base (MeV).

    Default ``deterministic=True``: edges on iff sigmoid(logit) > 0.5 (same as notebook
    ``affinity_rollout``). Stochastic sampling can rarely yield huge partition energies
    and distort the pre–REINFORCE baseline mean.
    """

    policy.eval()
    rets: list[float] = []
    gaps: list[float] = []
    lp: list[float] = []
    for _ in range(max(n, 1)):
        pos, mom, isp = sampler()
        obs = env.reset(pos, mom, isp)
        logits = policy(obs)
        if deterministic:
            a = (torch.sigmoid(logits) > 0.5).float()
        else:
            dist = torch.distributions.Bernoulli(logits=logits)
            a = dist.sample()
        loss, _ = env.physics_for_edge_mask(a)
        lb = float(env._baseline_loss)
        Lp = float(loss)
        lp.append(Lp)
        gaps.append(Lp - lb)
        rets.append(float(-Lp))
    policy.train()
    return float(np.mean(rets)), float(np.mean(lp)), float(np.mean(gaps))


def mean_tail(xs: list[float], tail: int) -> float:
    if not xs:
        return float("nan")
    k = min(int(tail), len(xs))
    return float(np.mean(xs[-k:]))


def main() -> int:
    args = parse_args()
    here = _CLUSTERING_ROOT
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

    if args.cold_all_off:
        init_policy_all_edges_off(policy, logit_bias=args.logit_bias)
        print("=== Mode: cold init (all edges off) ===")
    else:
        sup_meta = supervised_cache_meta(args, args.dataset)
        cache_path: Path | None = None
        if not args.no_supervised_cache:
            cache_path = (
                args.supervised_cache
                if args.supervised_cache is not None
                else default_supervised_cache_path(args.dataset, here)
            )
        loaded = False
        if cache_path is not None:
            loaded = try_load_supervised_cache(policy, cache_path, sup_meta)
        if loaded:
            print("=== Supervised warm-start: from cache (skipped training) ===")
        else:
            print("=== Supervised warm-start (baseline edge targets) ===")
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
                save_supervised_cache(policy, cache_path, sup_meta)

    G_sup, L_sup, gap_sup = rollout_stats(
        policy, env, sampler, args.eval_rollouts, deterministic=True
    )
    print(
        f"After {'cold' if args.cold_all_off else 'supervised'} (deterministic @0.5) — "
        f"mean return={G_sup:.2f} MeV  mean L_pol={L_sup:.2f} MeV  "
        f"mean gap={gap_sup:.2f} MeV"
    )

    rl_opt = optim.Adam(policy.parameters(), lr=args.lr)
    n_up = max(1, args.updates)
    # Cosine floor must stay below initial LR; a fixed eta_min=1e-5 inverts the
    # schedule when --lr is smaller (LR would climb toward 1e-5 each update).
    eta_min = min(1e-5, max(1e-8, float(args.lr) * 0.05))
    sched = CosineAnnealingLR(rl_opt, T_max=n_up, eta_min=eta_min)

    history = train_reinforce(
        policy=policy,
        env=env,
        event_sampler=sampler,
        optimizer=rl_opt,
        lr_scheduler=sched,
        n_updates=n_up,
        episodes_per_update=max(1, args.episodes_per_update),
        ent_coef=args.ent_coef,
        max_grad_norm=float(args.max_grad_norm),
        policy_coef=float(args.policy_coef),
        diag_jsonl=args.diag_jsonl,
        diag_every=max(1, args.diag_every),
        diag_det_rollouts=max(1, args.diag_det_rollouts),
        rl_action_mode=cast(RLActionMode, args.rl_action_mode),
        anchor_bce_coef=float(args.anchor_bce_coef),
        anchor_pos_weight=float(args.anchor_pos_weight),
        anchor_focal_gamma=float(args.anchor_focal_gamma),
    )
    if args.diag_jsonl is not None:
        print(f"REINFORCE diagnostics appended to {args.diag_jsonl.resolve()}")

    er = history.get("episode_return", [])
    Lp = history.get("partition_loss", [])
    Lb = history.get("baseline_loss", [])

    G_tail_stoch = mean_tail(er, args.tail)
    L_tail_stoch = mean_tail(Lp, args.tail)
    gap_tail_stoch = float("nan")
    if Lp and Lb and len(Lp) == len(Lb):
        k = min(args.tail, len(Lp))
        gap_tail_stoch = float(
            np.mean(
                np.asarray(Lp[-k:], dtype=np.float64)
                - np.asarray(Lb[-k:], dtype=np.float64)
            )
        )

    G_fin, L_fin, gap_fin = rollout_stats(
        policy, env, sampler, args.eval_rollouts, deterministic=True
    )

    dG_det = G_fin - G_sup
    dgap_det = gap_fin - gap_sup

    print("=== REINFORCE (kNN k=6) ===")
    print(
        f"updates={n_up}  ep/update={args.episodes_per_update}  "
        f"lr={args.lr}  ent_coef={args.ent_coef}  "
        f"rl_action_mode={args.rl_action_mode}  max_grad_norm={args.max_grad_norm}  "
        f"anchor_bce={args.anchor_bce_coef}  policy_coef={args.policy_coef}"
    )
    print(
        f"Training (stochastic): last {min(args.tail, len(er))} updates — "
        f"mean G={G_tail_stoch:.2f} MeV  mean L_pol={L_tail_stoch:.2f} MeV"
    )
    if np.isfinite(gap_tail_stoch):
        print(
            f"Training last-window mean gap (L_pol - L_base): {gap_tail_stoch:.2f} MeV"
        )

    print(
        f"After REINFORCE (deterministic @0.5, n={args.eval_rollouts}) — "
        f"mean return={G_fin:.2f} MeV  mean L_pol={L_fin:.2f} MeV  "
        f"mean gap={gap_fin:.2f} MeV"
    )

    print("--- Deterministic: supervised snapshot → after REINFORCE ---")
    print(f"Δ mean return: {dG_det:+.2f} MeV")
    print(f"Δ mean gap: {dgap_det:+.2f} MeV  (more negative gap is better)")

    improved = dG_det >= float(args.min_gain_mev)
    if not improved:
        improved = dgap_det <= -float(args.min_gain_mev)

    if improved:
        print(
            f"OK: deterministic gain ≥ {args.min_gain_mev:.1f} MeV "
            "on return and/or gap vs supervised snapshot."
        )
        return 0

    print(
        f"WARN: below gain threshold ({args.min_gain_mev:.1f} MeV). "
        "Try --sup-steps, --lr, --ent-coef, --updates, --episodes-per-update."
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
