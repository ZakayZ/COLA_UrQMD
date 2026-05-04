"""Spatial/momentum-cut baseline clustering and partition energy."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

from cluster_energy import partition_loss_numpy

R_CUT_FM = 7.0
Q_CUT_GEVC = 0.12
Q_CUT_MEVC = Q_CUT_GEVC * 1000.0

# Same ħc as kNN / GAT node ``k = p / ħc`` (MeV·fm); baseline adjacency uses (r, k) in fm and fm⁻¹.
HBARC_MEV_FM = 197.327
# Momentum cut expressed in wavenumber: equivalent to ``|Δp| < Q_CUT_MEVC`` when ``p`` is in MeV/c.
K_CUT_FM_INV = Q_CUT_MEVC / HBARC_MEV_FM


def baseline_clusters_numpy(
    pos: np.ndarray,
    mom: np.ndarray,
    indices: list[int],
    r_cut_fm: float,
    q_cut_momentum: float,
) -> list[list[int]]:
    """Cut-based clusters: neighbors if ``|Δr| < r_cut`` and ``|Δk| < k_cut``.

    ``pos`` is ``(N, 3)`` ``(x,y,z)`` fm or ``(N, 4)`` ``(t, x, y, z)`` (only spatial columns are used).
    ``mom`` is ``(N, 4)`` with ``(E, px, py, pz)`` in **MeV/c**;
    the cut uses ``k = p_spatial / ħc`` (fm⁻¹), matching the kNN feature space in ``models.env``.
    ``q_cut_momentum`` is the same **MeV/c** scale as ``|Δp|`` (internally converted to ``k_cut``).
    """
    pos = np.asarray(pos, dtype=np.float64)
    if pos.shape[-1] == 4:
        pos = pos[..., 1:4]

    n = len(indices)
    mom = np.asarray(mom, dtype=np.float64)
    k3 = mom[:, 1:4] / HBARC_MEV_FM
    k_cut = float(q_cut_momentum) / HBARC_MEV_FM
    adj = [[] for _ in range(n)]
    for a in range(n):
        ia = indices[a]
        for b in range(a + 1, n):
            ib = indices[b]
            if float(np.linalg.norm(pos[ia] - pos[ib])) < r_cut_fm and float(
                np.linalg.norm(k3[ia] - k3[ib])
            ) < k_cut:
                adj[a].append(b)
                adj[b].append(a)
    used = [False] * n
    comps: list[list[int]] = []
    for s in range(n):
        if used[s]:
            continue
        stack = [s]
        used[s] = True
        comp_local: list[int] = []
        while stack:
            v = stack.pop()
            comp_local.append(indices[v])
            for to in adj[v]:
                if not used[to]:
                    used[to] = True
                    stack.append(to)
        comps.append(sorted(comp_local))
    comps.sort(key=len, reverse=True)
    return comps


@dataclass(frozen=True)
class EventBaseline:
    """Per-event baseline: node-wise cluster ids, partition loss, and cluster lists.

    ``loss`` is the partition energy in **MeV** (same scale as :func:`~cluster_energy.partition_loss_numpy`).
    """

    node_labels: np.ndarray
    loss: float
    partition: list[list[int]]


def compute_event_baseline(
    pos: np.ndarray,
    mom: np.ndarray,
    is_proton: np.ndarray,
    *,
    r_cut_fm: float = R_CUT_FM,
    q_cut_momentum: float = Q_CUT_MEVC,
) -> EventBaseline:
    """Run cut-based clustering (``|Δr|``, ``|Δk|``) and partition loss in **MeV**."""
    n_ev = int(pos.shape[0])
    part_b = baseline_clusters_numpy(pos, mom, list(range(n_ev)), r_cut_fm, q_cut_momentum)
    node_lab = np.empty((n_ev,), dtype=np.int32)
    for ci, c in enumerate(part_b):
        node_lab[np.asarray(c, dtype=np.int64)] = int(ci)
    loss = float(partition_loss_numpy(pos, mom, is_proton, part_b))
    return EventBaseline(node_labels=node_lab, loss=loss, partition=part_b)


def edge_pair_baseline_targets(
    node_labels: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> torch.Tensor:
    """Per-edge binary targets: 1 iff baseline puts both endpoints in the same cluster (CPU)."""
    tgt = (node_labels[edge_i] == node_labels[edge_j]).astype(np.float32)
    return torch.tensor(tgt, dtype=torch.float32)
