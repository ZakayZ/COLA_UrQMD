#!/usr/bin/env python3
"""Compare A2C vs PPO vs physics baseline after shared supervised init (urqmd nucleon graphs)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, cast

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

from models import GATAffinityActorCritic
from models import (
    AffinityGraphConfig,
    AffinityGraphEnv,
    GAT_NODE_IN_DIM,
    GATAffinityPolicy,
    load_valid_events_from_pkl,
)
from training.utils import RLActionMode, deterministic_eval_mean, make_event_sampler
from training.a2c import train_actor_critic
from training.ppo import train_ppo
from training.supervised import train_supervised_edges

HERE = _CLUSTERING_ROOT


def clone_state_dict(net: torch.nn.Module) -> dict[str, torch.Tensor]:
    return {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}


def build_policy(
    *,
    hidden: int,
    n_heads: int,
    n_gat_layers: int,
    running_norm: bool,
) -> GATAffinityPolicy:
    return GATAffinityPolicy(
        in_dim=GAT_NODE_IN_DIM,
        hidden=hidden,
        n_heads=n_heads,
        n_gat_layers=n_gat_layers,
        running_norm=running_norm,
    )


def eval_gap_means(
    policy: GATAffinityPolicy,
    env: AffinityGraphEnv,
    sampler,
    n_rollouts: int,
) -> tuple[float, float, float]:
    """Mean return (MeV), mean L_pol (MeV), mean gap L_pol - L_base (MeV)."""
    G, Lp, gap = deterministic_eval_mean(policy, env, sampler, n_rollouts)
    return G, Lp, gap


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--dataset",
        type=Path,
        default=HERE / "datasets" / "urqmd_nucleons_1k" / "dataset.pkl",
    )
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--sup-steps", type=int, default=80)
    p.add_argument("--sup-events-per-step", type=int, default=8)
    p.add_argument("--sup-lr", type=float, default=1e-3)
    p.add_argument("--sup-pos-weight", type=float, default=20.0)
    p.add_argument("--sup-focal-gamma", type=float, default=2.0)
    p.add_argument("--rl-steps", type=int, default=80, help="A2C / PPO outer updates each.")
    p.add_argument("--episodes-per-update", type=int, default=14)
    p.add_argument("--eval-rollouts", type=int, default=96)
    p.add_argument("--hidden", type=int, default=96)
    p.add_argument("--n-heads", type=int, default=4)
    p.add_argument("--n-gat-layers", type=int, default=3)
    p.add_argument("--value-mlp-hidden", type=int, default=112)
    p.add_argument(
        "--running-norm",
        action="store_true",
        help="Enable BatchNorm on nodes/edges (default off like run_reinforce script).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.dataset.is_file():
        print(f"Dataset not found: {args.dataset}", file=sys.stderr)
        return 1

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)

    events = load_valid_events_from_pkl(args.dataset)
    if not events:
        print("No events loaded.", file=sys.stderr)
        return 1

    cfg = AffinityGraphConfig(k_nn=6, graph_kind="knn")
    env = AffinityGraphEnv(cfg)
    sampler = make_event_sampler(events=events, rng=rng, fallback_urqmd=None)

    policy_sup = build_policy(
        hidden=args.hidden,
        n_heads=args.n_heads,
        n_gat_layers=args.n_gat_layers,
        running_norm=args.running_norm,
    )
    print(
        f"=== Supervised warm-start: hidden={args.hidden} layers={args.n_gat_layers} "
        f"heads={args.n_heads} running_norm={args.running_norm} ==="
    )
    sup_opt = optim.Adam(policy_sup.parameters(), lr=float(args.sup_lr))
    train_supervised_edges(
        policy_sup,
        env,
        sampler,
        steps=max(1, args.sup_steps),
        events_per_step=max(1, args.sup_events_per_step),
        lr=float(args.sup_lr),
        focal_gamma=float(args.sup_focal_gamma),
        pos_weight=float(args.sup_pos_weight),
        optimizer=sup_opt,
        max_grad_norm=0.5,
    )
    sup_sd = clone_state_dict(policy_sup)

    G_sup, L_sup, gap_sup = eval_gap_means(
        policy_sup, env, sampler, args.eval_rollouts
    )
    print(
        f"After supervised (deterministic τ=0.5): mean G={G_sup:.3f} MeV  "
        f"L_pol={L_sup:.3f} MeV  gap=L_pol-L_base={gap_sup:.3f} MeV"
    )

    tail = min(15, max(1, args.rl_steps // 5))
    rl_steps = max(50, min(100, int(args.rl_steps)))

    # Tuned for MeV-scale returns: moderate value_coef, small ent, anchor toward baseline edges.
    experiments: list[tuple[str, str, dict[str, Any]]] = [
        (
            "A2C",
            "lr8e-6_vc0.12_ent0.005",
            dict(
                lr=8e-6,
                value_coef=0.12,
                ent_coef=0.005,
                anchor_bce_coef=0.008,
                center_adv=True,
            ),
        ),
        (
            "A2C",
            "lr1.2e-5_vc0.18_ent0.004",
            dict(
                lr=1.2e-5,
                value_coef=0.18,
                ent_coef=0.004,
                anchor_bce_coef=0.006,
                center_adv=True,
            ),
        ),
        (
            "PPO",
            "lr2.5e-5_clip0.15_vc0.12",
            dict(
                lr=2.5e-5,
                clip_range=0.15,
                value_clip_range=0.15,
                value_coef=0.12,
                ent_coef=0.005,
                ppo_epochs=4,
                minibatch_size=min(8, args.episodes_per_update),
                anchor_bce_coef=0.008,
                normalize_advantage=True,
            ),
        ),
        (
            "PPO",
            "lr1.5e-5_clip0.22_vc0.15_ep3",
            dict(
                lr=1.5e-5,
                clip_range=0.22,
                value_clip_range=0.2,
                value_coef=0.15,
                ent_coef=0.004,
                ppo_epochs=3,
                minibatch_size=max(1, args.episodes_per_update // 2),
                anchor_bce_coef=0.006,
                normalize_advantage=True,
            ),
        ),
    ]

    results: list[dict[str, Any]] = []

    for algo, tag, kw in experiments:
        pol = build_policy(
            hidden=args.hidden,
            n_heads=args.n_heads,
            n_gat_layers=args.n_gat_layers,
            running_norm=args.running_norm,
        )
        pol.load_state_dict(sup_sd)
        ac = GATAffinityActorCritic(pol, value_mlp_hidden=int(args.value_mlp_hidden))
        opt = optim.Adam(ac.parameters(), lr=float(kw["lr"]))
        eta_min = min(1e-6, max(1e-9, float(kw["lr"]) * 0.08))
        sched = CosineAnnealingLR(opt, T_max=rl_steps, eta_min=eta_min)

        common_kw = dict(
            optimizer=opt,
            lr_scheduler=sched,
            n_updates=rl_steps,
            episodes_per_update=max(1, args.episodes_per_update),
            max_grad_norm=1.0,
            rl_action_mode=cast(RLActionMode, "threshold"),
            anchor_pos_weight=15.0,
            anchor_focal_gamma=2.0,
        )

        if algo == "A2C":
            hist = train_actor_critic(
                ac,
                env,
                sampler,
                ent_coef=float(kw["ent_coef"]),
                value_coef=float(kw["value_coef"]),
                anchor_bce_coef=float(kw["anchor_bce_coef"]),
                center_adv=bool(kw["center_adv"]),
                **common_kw,
            )
        else:
            hist = train_ppo(
                ac,
                env,
                sampler,
                clip_range=float(kw["clip_range"]),
                value_clip_range=kw.get("value_clip_range"),
                value_coef=float(kw["value_coef"]),
                ent_coef=float(kw["ent_coef"]),
                anchor_bce_coef=float(kw["anchor_bce_coef"]),
                ppo_epochs=int(kw["ppo_epochs"]),
                minibatch_size=int(kw["minibatch_size"]),
                normalize_advantage=bool(kw["normalize_advantage"]),
                **common_kw,
            )

        G_fin, L_fin, gap_fin = eval_gap_means(
            ac.policy, env, sampler, args.eval_rollouts
        )

        lp_tail = hist.get("partition_loss", [])
        lb_tail = hist.get("baseline_loss", [])
        gap_tail = float("nan")
        if lp_tail and lb_tail and len(lp_tail) == len(lb_tail):
            k = min(tail, len(lp_tail))
            gap_tail = float(
                np.mean(
                    np.asarray(lp_tail[-k:], dtype=np.float64)
                    - np.asarray(lb_tail[-k:], dtype=np.float64)
                )
            )

        row = {
            "algo": algo,
            "tag": tag,
            "gap_sup_MeV": gap_sup,
            "gap_fin_MeV": gap_fin,
            "delta_gap_MeV": gap_fin - gap_sup,
            "L_pol_sup_MeV": L_sup,
            "L_pol_fin_MeV": L_fin,
            "G_fin_MeV": G_fin,
            "gap_train_tail_mean_MeV": gap_tail,
            "kw": kw,
        }
        results.append(row)

        print(
            f"\n--- {algo} [{tag}] ---\n"
            f"  deterministic eval (n={args.eval_rollouts}): "
            f"gap {gap_sup:.3f} → {gap_fin:.3f} MeV (Δ {gap_fin - gap_sup:+.3f})\n"
            f"  train tail (~{tail} upd) mean gap (stochastic rollout): {gap_tail:.3f} MeV"
        )

    best = min(results, key=lambda r: r["gap_fin_MeV"])
    print("\n=== Summary (lower gap_fin is better; baseline reference is per-event L_base) ===")
    for r in results:
        print(
            f"  {r['algo']:4} {r['tag']:28}  gap_fin={r['gap_fin_MeV']:+.4f} MeV  "
            f"Δ vs sup={r['delta_gap_MeV']:+.4f} MeV"
        )
    print(
        f"\nBest deterministic gap after RL: {best['algo']} / {best['tag']} "
        f"→ gap_fin={best['gap_fin_MeV']:.4f} MeV (supervised gap was {gap_sup:.4f} MeV)"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
