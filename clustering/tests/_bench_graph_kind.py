#!/usr/bin/env python3
"""Compare graph topology for supervised pretrain with running_norm=False (no BN)."""
from __future__ import annotations

import sys
from pathlib import Path

_CLUSTERING_ROOT = Path(__file__).resolve().parent.parent
_TESTS_DIR = Path(__file__).resolve().parent
for _root in (_CLUSTERING_ROOT, _TESTS_DIR):
    _sr = str(_root)
    if _sr not in sys.path:
        sys.path.insert(0, _sr)

import numpy as np

from models import (
    AffinityGraphConfig,
    AffinityGraphEnv,
    GATAffinityPolicy,
    GAT_NODE_IN_DIM,
    load_valid_events_from_pkl,
)
from training.utils import make_event_sampler
from training.supervised import train_supervised_edges

LR = 1e-3
STEPS = 40
EVENTS_PER_STEP = 8


def run_one(name: str, cfg: AffinityGraphConfig, *, steps: int, events_per_step: int) -> tuple[str, float, float, float]:
    rng = np.random.default_rng(42)
    events = load_valid_events_from_pkl(_CLUSTERING_ROOT / "datasets/urqmd_nucleons_1k/dataset.pkl")
    env = AffinityGraphEnv(cfg)
    policy = GATAffinityPolicy(
        in_dim=GAT_NODE_IN_DIM,
        hidden=64,
        n_heads=4,
        n_gat_layers=2,
        running_norm=False,
    )
    sampler = make_event_sampler(events=events, rng=rng, fallback_urqmd=None)
    hist = train_supervised_edges(
        policy=policy,
        env=env,
        event_sampler=sampler,
        steps=steps,
        events_per_step=events_per_step,
        lr=LR,
        focal_gamma=2.0,
    )
    gap = float(np.mean(hist["pretrain_gap"][-5:]))
    rec = float(np.mean(hist["pretrain_pos_recall_05"][-5:]))
    bce = float(np.mean(hist["supervised_bce"][-5:]))
    return name, gap, rec, bce


def main() -> None:
    # Same schedule for every topology (full graph is O(N²) per event — expect long runtime).
    cfgs: list[tuple[str, AffinityGraphConfig]] = [
        ("kNN k=3", AffinityGraphConfig(k_nn=3, graph_kind="knn")),
        ("kNN k=5", AffinityGraphConfig(k_nn=5, graph_kind="knn")),
        ("kNN k=8", AffinityGraphConfig(k_nn=8, graph_kind="knn")),
        ("radius r=0.9", AffinityGraphConfig(graph_kind="radius", radius_norm=0.9)),
        ("radius r=1.15", AffinityGraphConfig(graph_kind="radius", radius_norm=1.15)),
        ("radius r=1.45", AffinityGraphConfig(graph_kind="radius", radius_norm=1.45)),
        ("full complete", AffinityGraphConfig(graph_kind="full", k_nn=6)),
    ]

    rows: list[tuple[str, float, float, float]] = []
    for name, cfg in cfgs:
        print(
            f"=== {name} (steps={STEPS}, events_per_step={EVENTS_PER_STEP}) ===",
            flush=True,
        )
        rows.append(
            run_one(
                name,
                cfg,
                steps=STEPS,
                events_per_step=EVENTS_PER_STEP,
            )
        )

    print(
        f"\nlr={LR} steps={STEPS} events_per_step={EVENTS_PER_STEP} "
        f"focal=2.0 running_norm=False (metrics = mean last 5 steps)\n"
    )
    hdr = f"{'topology':<18} {'gap5':>10} {'rec5':>8} {'bce5':>8}"
    print(hdr)
    print("-" * len(hdr))
    for name, gap, rec, bce in rows:
        print(f"{name:<18} {gap:>10.1f} {rec:>8.3f} {bce:>8.4f}")
    print("\nIgnore recall≈1 (too dense). Compare manually; gaps fluctuate run-to-run (sampler noise).")


if __name__ == "__main__":
    main()
