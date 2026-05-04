#!/usr/bin/env python3
"""One-off benchmark: supervised edge training variants. Run ``python tests/_bench_supervised_variants.py`` from ``clustering/``."""
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

STEPS = 80
EVENTS_PER_STEP = 8
LR = 1e-3


def main() -> None:
    events = load_valid_events_from_pkl(_CLUSTERING_ROOT / "datasets/urqmd_nucleons_1k/dataset.pkl")

    runs: list[tuple[str, bool, float, AffinityGraphConfig]] = [
        ("BN + focal γ=2 + scales 50/2/5 + cuts", True, 2.0, AffinityGraphConfig(k_nn=6)),
        ("no BatchNorm (Identity)", False, 2.0, AffinityGraphConfig(k_nn=6)),
        ("BN + no focal (γ=0)", True, 0.0, AffinityGraphConfig(k_nn=6)),
        ("raw r,E,k (×1) + BN + focal", True, 2.0, AffinityGraphConfig(
            k_nn=6, feat_scale_r_fm=1.0, feat_scale_e=1000.0, feat_scale_k_fm_inv=1.0)),
        ("raw r,E,k + no BN + focal", False, 2.0, AffinityGraphConfig(
            k_nn=6, feat_scale_r_fm=1.0, feat_scale_e=1000.0, feat_scale_k_fm_inv=1.0)),
    ]

    print(f"dataset events={len(events)}  steps={STEPS}  events/step={EVENTS_PER_STEP}  lr={LR}\n")

    results: list[tuple[str, float, float, float, float, float, float]] = []
    for name, use_bn, focal, cfg in runs:
        # Same event sequence for every row so gaps are comparable.
        rng = np.random.default_rng(42)
        env = AffinityGraphEnv(cfg)
        policy = GATAffinityPolicy(
            in_dim=GAT_NODE_IN_DIM,
            hidden=64,
            n_heads=4,
            n_gat_layers=2,
            running_norm=use_bn,
        )
        sampler = make_event_sampler(events=events, rng=rng, fallback_urqmd=None)

        hist = train_supervised_edges(
            policy=policy,
            env=env,
            event_sampler=sampler,
            steps=STEPS,
            events_per_step=EVENTS_PER_STEP,
            lr=LR,
            focal_gamma=focal,
        )

        gap_last = float(hist["pretrain_gap"][-1])
        gap_avg10 = float(np.mean(hist["pretrain_gap"][-10:]))
        rec_last = float(hist["pretrain_pos_recall_05"][-1])
        rec_avg10 = float(np.mean(hist["pretrain_pos_recall_05"][-10:]))
        bce_last = float(hist["supervised_bce"][-1])
        l_pol_last = float(hist["pretrain_partition_loss"][-1])
        results.append((name, gap_last, gap_avg10, rec_last, rec_avg10, bce_last, l_pol_last))

    hdr = (
        f"{'config':<40} {'gap_last':>9} {'gap_avg10':>10} "
        f"{'rec_last':>9} {'rec_avg10':>10} {'bce':>8} {'L_pol':>10}"
    )
    print(hdr)
    print("-" * len(hdr))
    for row in results:
        name, gl, ga, rl, ra, bce, lp = row
        print(
            f"{name:<40} {gl:>9.1f} {ga:>10.1f} "
            f"{rl:>9.3f} {ra:>10.3f} {bce:>8.4f} {lp:>10.1f}"
        )
    print()
    print()
    print("gap_* = mean(L_pol − L_base) over the micro-batch (MeV); baseline uses spatial–momentum cut.")
    print("gap_avg10 smooths noise over the last 10 steps. rec = positive-edge recall at σ>0.5.")


if __name__ == "__main__":
    main()
