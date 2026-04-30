"""Standalone runner for edge-collapse debugging in REINFORCE + GAT training."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import torch.optim as optim
from torch.optim.lr_scheduler import LinearLR

from ppo_train import (
    AffinityGraphConfig,
    GAT_NODE_IN_DIM,
    GATAffinityPolicy,
    baseline_edge_targets,
    load_valid_events_from_pkl,
    make_event_sampler,
    train_reinforce,
)
from ppo_train import AffinityGraphEnv


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    default_dataset = here / "datasets" / "urqmd_nucleons_1k" / "dataset.pkl"
    p = argparse.ArgumentParser(description="Run a reproducible REINFORCE edge-debug pass.")
    p.add_argument("--dataset", type=Path, default=default_dataset, help="Path to dataset.pkl")
    p.add_argument("--seed", type=int, default=1234, help="RNG seed")
    p.add_argument("--updates", type=int, default=30, help="PPO update steps")
    p.add_argument("--episodes-per-update", type=int, default=8, help="Episodes per update")
    p.add_argument("--probe-events", type=int, default=64, help="Baseline probe events before training")
    p.add_argument("--k-nn", type=int, default=12, help="k-NN for edge candidates")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.dataset.exists():
        raise FileNotFoundError(f"Dataset not found: {args.dataset}")

    events = load_valid_events_from_pkl(args.dataset)
    if not events:
        raise RuntimeError(f"No valid events found in {args.dataset}")

    rng = np.random.default_rng(args.seed)
    cfg = AffinityGraphConfig(k_nn=args.k_nn)
    env = AffinityGraphEnv(cfg)
    policy = GATAffinityPolicy(
        in_dim=GAT_NODE_IN_DIM,
        hidden=128,
        n_heads=4,
        n_gat_layers=2,
    )
    sampler = make_event_sampler(events=events, rng=rng, fallback_urqmd=None)
    optimizer = optim.Adam(policy.parameters(), lr=3e-3)
    n_up = max(args.updates, 1)
    lr_scheduler = LinearLR(
        optimizer,
        start_factor=1.0,
        end_factor=(1e-4 / 3e-3),
        total_iters=n_up,
    )

    print("Running baseline edge-density probe...")
    for _ in range(max(args.probe_events, 1)):
        pos, mom, isp = sampler()
        env.reset(pos, mom, isp)
        baseline_edge_targets(env)  # Triggers baseline target-ratio instrumentation.

    print("Running REINFORCE debug pass...")
    history = train_reinforce(
        policy=policy,
        env=env,
        event_sampler=sampler,
        optimizer=optimizer,
        lr_scheduler=lr_scheduler,
        n_updates=args.updates,
        episodes_per_update=args.episodes_per_update,
        ent_coef=0.02,
        max_grad_norm=0.5,
        policy_coef=1.0,
    )

    n_tail = min(10, len(history.get("partition_loss", [])))
    if n_tail == 0:
        print("No valid training updates were recorded.")
        return

    l_pol = float(np.mean(history["partition_loss"][-n_tail:]))
    l_base = float(np.mean(history["baseline_loss"][-n_tail:])) if history["baseline_loss"] else float("nan")
    ret = float(np.mean(history["episode_return"][-n_tail:])) if history["episode_return"] else float("nan")
    ent = float(np.mean(history["edge_entropy"][-n_tail:])) if history["edge_entropy"] else float("nan")
    print(f"tail{n_tail}: L_pol={l_pol:.1f} MeV  L_base={l_base:.1f} MeV  G={ret:.1f}  H={ent:.5f}")
    print("Debug run complete.")


if __name__ == "__main__":
    main()
