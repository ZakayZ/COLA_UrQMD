import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import numpy as np
import torch
from scipy.sparse import csr_matrix
from scipy.sparse import csgraph
import torch.nn as nn
import torch.nn.functional as F
from torch.optim.lr_scheduler import LRScheduler
from torch_geometric.data import Data
import torch_geometric.nn as pyg_nn
from torch_geometric.transforms import KNNGraph
from tqdm.auto import tqdm

from cluster_energy import partition_loss_numpy

from baseline import (
    HBARC_MEV_FM,
    Q_CUT_GEVC,
    Q_CUT_MEVC,
    R_CUT_FM,
    baseline_clusters_numpy,
    compute_event_baseline,
    edge_pair_baseline_targets,
)
from datasets.generate_urqmd_nucleon_dataset import load_dataset_pickle

import colapy

# ``partition_loss_numpy`` returns GeV; training stores **MeV** internally. Use this only for GeV labels/plots.
MEV_PER_GEV = 1000.0

# Per-edge MLP input after ``h_i, h_j``: three copies of the 8-D ``dist`` row (see ``reset``):
# ``Δdist``, ``dist_i``, ``dist_j`` → 8 + 8 + 8 = 24.
EDGE_PHYS_DIM = 24

# GAT node ``x``: r (3 fm), ``t / ħc`` (1, same scale as k = p/ħc), k (3 fm⁻¹), E (1 GeV), is_proton (1).
GAT_NODE_IN_DIM = 9


def extract_nucleons_numpy(particles: list[Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert COLA particles to numpy arrays.

    ``pos`` is ``(N, 4)``: ``(t, x, y, z)`` with ``t`` in fm/c and ``x,y,z`` in fm. ``mom`` is ``(E, px, py, pz)`` in MeV/c.
    """
    pos, mom, is_proton = [], [], []
    for p in particles:
        if p.pdg_code == 2212:
            mom.append([p.momentum.e, p.momentum.x, p.momentum.y, p.momentum.z])
            pos.append([p.position.t, p.position.x, p.position.y, p.position.z])
            is_proton.append(True)
        elif p.pdg_code == 2112:
            mom.append([p.momentum.e, p.momentum.x, p.momentum.y, p.momentum.z])
            pos.append([p.position.t, p.position.x, p.position.y, p.position.z])
            is_proton.append(False)
    if not pos:
        return np.zeros((0, 4), np.float64), np.zeros((0, 4), np.float64), np.zeros((0,), bool)
    return np.asarray(pos, np.float64), np.asarray(mom, np.float64), np.asarray(is_proton, bool)


def try_make_urqmd_event_generator() -> Callable[[], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """COLA / UrQMD: one event → ``(pos, mom, is_proton)`` numpy."""

    class W(colapy.WriterBase):
        events: list[Any] = []

        def __init__(self, **kwargs):
            self.events.clear()

        def __call__(self, event_data):
            self.events.append(event_data)

    config = """
<?xml version="1.0" encoding="UTF-8" ?>
<program>
    <generator name="URQMDGenerator"
        pro="197 79"
        tar="197 79"
        nev="1"
        imp="5."
        elb="100."
        tim="200 200"
        generated_config_file="input_file"/>
    <writer name="PythonWriter" class="W"/>
</program>
"""

    def gen_one() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".xml", delete_on_close=False) as tmp:
            tmp.write(config)
            tmp.close()
            rm = colapy.RunManager().load_module("COLA-Py").load_module("COLA_UrQMD").load_config(tmp.name)
            rm.run(1)
            Path("input_file").unlink(missing_ok=True)
        if not W.events:
            return np.zeros((0, 4)), np.zeros((0, 4)), np.zeros((0,), bool)
        ev = W.events[-1]
        return extract_nucleons_numpy(ev.particles)

    return gen_one


def load_valid_events_from_pkl(pkl: Path) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Load valid events from generated dataset pickle."""
    bundle = load_dataset_pickle(pkl)
    return [(e["pos"], e["mom"], e["is_proton"]) for e in bundle["events"] if e]


def _build_knn_graph(pos: np.ndarray, mom_phys_or_k: np.ndarray, transform: KNNGraph) -> Data:
    """kNN in ``(r, k)`` (6-D): ``r`` fm; second arg is either ``mom`` ``(N,4)`` MeV/c or ``k`` ``(N,3)`` fm⁻¹."""
    m = np.asarray(mom_phys_or_k, dtype=np.float64)
    if m.shape[1] == 4:
        k = m[:, 1:4] / HBARC_MEV_FM
    else:
        k = m
    six = np.concatenate([np.asarray(pos, dtype=np.float64), k], axis=1).astype(np.float32)
    d = transform(Data(pos=torch.tensor(six)))
    src, dst = d.edge_index[0], d.edge_index[1]
    keep = src < dst
    d.edge_pair_i = src[keep]
    d.edge_pair_j = dst[keep]
    d.policy_edge_idx = torch.where(keep)[0].long()
    return d


def labels_to_partition(labels: np.ndarray) -> list[list[int]]:
    by_lbl: dict[int, list[int]] = defaultdict(list)
    for i, lab in enumerate(labels.astype(int).tolist()):
        by_lbl[int(lab)].append(i)
    return [by_lbl[k] for k in sorted(by_lbl.keys())]


def cluster_labels_from_edges(
    n: int,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    edge_on: np.ndarray,
) -> np.ndarray:
    """Connected components on the subgraph of **on** edges (labels ``0 .. n_comp-1``).

    Uses ``scipy.sparse.csgraph.connected_components`` (CSR adjacency, undirected).
    """
    on = np.asarray(edge_on, dtype=bool).reshape(-1)
    ei = np.asarray(edge_i, dtype=np.int64)[on]
    ej = np.asarray(edge_j, dtype=np.int64)[on]

    _, labels = csgraph.connected_components(
        csr_matrix((np.ones(ei.size, dtype=np.int8), (ei, ej)), shape=(n, n)),
        directed=False,
        return_labels=True,
    )
    return labels.astype(np.int32, copy=False)


@dataclass
class AffinityGraphConfig:
    k_nn: int = 5


class AffinityGraphEnv:
    """``reset`` stores ``mom``, ``isp`` as given and normalizes ``pos`` to ``(N, 4)`` ``(t,x,y,z)``.

    If ``pos`` is ``(N, 3)`` spatial fm only, a column of zeros is prepended as coordinate time.

    ``reset`` also fills cut-baseline fields ``_baseline_node_labels`` / ``_baseline_loss`` via
    :func:`baseline.compute_event_baseline`.

    The returned ``Data`` has node features ``x``, ``edge_index`` / ``edge_attr`` from
    ``KNNGraph`` (one ``edge_attr`` row per directed edge), plus ``edge_pair_*`` and
    ``policy_edge_idx`` for the undirected policy edge list (``src < dst``).

    After ``reset``, ``graph`` is the single ``Data`` passed to the policy."""

    def __init__(self, cfg: AffinityGraphConfig) -> None:
        self.cfg = cfg
        self._knn = KNNGraph(k=cfg.k_nn, loop=False, force_undirected=True)
        self.pos = np.zeros((0, 4), dtype=np.float64)
        self.mom = np.zeros((0, 4), dtype=np.float64)
        self.isp = np.zeros((0,), dtype=bool)
        self.graph = Data()
        self._baseline_node_labels = np.zeros((0,), dtype=np.int32)
        self._baseline_loss = 0.0
        self._partition_loss = 0.0

    @property
    def n_real_edges(self) -> int:
        return int(self.graph.edge_pair_i.shape[0])

    def reset(self, pos: np.ndarray, mom: np.ndarray, is_proton: np.ndarray) -> Data:
        pos = np.asarray(pos, dtype=np.float64)
        if pos.shape[1] == 3:
            pos_geo = np.concatenate([np.zeros((pos.shape[0], 1), dtype=np.float64), pos], axis=1)
        elif pos.shape[1] == 4:
            pos_geo = pos
        else:
            raise ValueError(f"pos must be (N, 3) or (N, 4), got shape {pos.shape}")

        k3 = mom[:, 1:] / HBARC_MEV_FM
        t = pos_geo[:, :1] / HBARC_MEV_FM
        r3 = pos_geo[:, 1:4]
        e = mom[:, :1] / 1000.0

        self.pos = np.array(pos_geo, dtype=np.float64, copy=True)
        self.mom = mom
        self.isp = is_proton

        self.graph = _build_knn_graph(r3, k3, self._knn)

        phase_space = np.concatenate(
            [
                t,
                r3,
                e,
                k3,
            ],
            axis=1,
            dtype=np.float32,
        )
        self.graph.x = torch.from_numpy(np.concatenate(
            [
                phase_space,
                is_proton[:, None],
            ],
            axis=1,
            dtype=np.float32,
        ))
        phase_space_t = torch.from_numpy(phase_space)
        i, j = self.graph.edge_index
        self.graph.edge_attr = torch.cat(
            [phase_space_t[i] - phase_space_t[j], phase_space_t[i], phase_space_t[j]],
            dim=1,
        )
        bl = compute_event_baseline(self.pos, self.mom, self.isp)
        self._baseline_node_labels = bl.node_labels.astype(np.int32)
        self._baseline_loss = float(bl.loss)
        return self.graph

    def physics_for_edge_mask(self, edge_on: torch.Tensor) -> tuple[float, np.ndarray]:
        n = self.pos.shape[0]
        on = edge_on.detach().numpy().astype(bool, copy=False).reshape(-1)
        labels = cluster_labels_from_edges(
            n, self.graph.edge_pair_i.numpy(), self.graph.edge_pair_j.numpy(), on
        )
        part = labels_to_partition(labels)
        pl = float(partition_loss_numpy(self.pos, self.mom, self.isp, part)) * MEV_PER_GEV
        self._partition_loss = pl
        return pl, labels


class _MetaEdgeFeatureMLP(nn.Module):
    """``edge_model`` for :class:`torch_geometric.nn.MetaLayer` (see ``torch_geometric.nn.models.meta`` docstring)."""

    def __init__(self, mlp: nn.Module) -> None:
        super().__init__()
        self.mlp = mlp

    def forward(
        self,
        src: torch.Tensor,
        dst: torch.Tensor,
        edge_attr: torch.Tensor,
        u: torch.Tensor | None,
        batch: torch.Tensor | None,
    ) -> torch.Tensor:
        if u is not None or batch is not None:
            raise NotImplementedError("edge readout supports single graphs (u=batch=None) only")
        return self.mlp(torch.cat([src, dst, edge_attr], dim=-1))


class GATAffinityPolicy(nn.Module):
    """Edge-aware stack: ``GATv2Conv`` + ``edge_attr``, ``Sequential`` encoder, ``MetaLayer`` + ``MLP`` edge readout."""

    edge_phys_dim: int = EDGE_PHYS_DIM
    max_value: float = 100.0

    def __init__(
        self,
        in_dim: int,
        hidden: int,
        n_heads: int,
        n_gat_layers: int = 2,
        edge_mlp_depth: int = 3,
        edge_mlp_hidden: int | None = None,
    ) -> None:
        super().__init__()

        head_dim = hidden // n_heads
        h_e = edge_mlp_hidden if edge_mlp_hidden is not None else hidden
        edge_in = 2 * hidden + EDGE_PHYS_DIM
        channels = [edge_in] + [h_e] * edge_mlp_depth + [1]

        enc_layers: list[tuple[nn.Module, str]] = []
        in_ch = in_dim
        for layer_id in range(n_gat_layers):
            enc_layers.append(
                (
                    pyg_nn.GATv2Conv(
                        in_ch,
                        head_dim,
                        heads=n_heads,
                        concat=True,
                        edge_dim=EDGE_PHYS_DIM,
                        add_self_loops=False,
                        residual=(layer_id > 0),
                    ),
                    "x, edge_index, edge_attr -> x",
                )
            )
            enc_layers.append((nn.ELU(), "x -> x"))
            enc_layers.append((pyg_nn.norm.LayerNorm(hidden), "x -> x"))
            in_ch = hidden

        self.encoder = pyg_nn.Sequential("x, edge_index, edge_attr", enc_layers)
        edge_mlp = pyg_nn.models.MLP(
            channels,
            dropout=0.0,
            act="gelu",
            norm="layer_norm",
            plain_last=True,
        )
        self.edge_readout = pyg_nn.MetaLayer(
            edge_model=_MetaEdgeFeatureMLP(edge_mlp),
            node_model=None,
            global_model=None,
        )

    def forward(self, data: Data) -> torch.Tensor:
        h = self.encoder(data.x, data.edge_index, data.edge_attr)
        _, le, _ = self.edge_readout(
            h, data.edge_index, data.edge_attr, u=None, batch=None
        )
        le = le[data.policy_edge_idx]
        return torch.nan_to_num(le.view(-1), nan=0.0, posinf=self.max_value, neginf=-self.max_value).clamp(-self.max_value, self.max_value)


def collect_rollout(
    env: AffinityGraphEnv,
    policy: GATAffinityPolicy,
    pos: np.ndarray,
    mom: np.ndarray,
    isp: np.ndarray,
) -> dict[str, Any]:
    obs = env.reset(pos, mom, isp)
    edge_logits = policy(obs)
    dist = torch.distributions.Bernoulli(logits=edge_logits)
    a = dist.sample()
    ent = dist.entropy().mean()
    loss, labs = env.physics_for_edge_mask(a)
    l_base = float(env._baseline_loss)
    r = -loss
    return {
        "obs": obs,
        "action": a,
        "reward": r,
        "partition_loss": loss,
        "baseline_loss": l_base,
        "edge_entropy": ent,
        "n_clusters": int(len(np.unique(labs))),
    }


def baseline_edge_targets(env: AffinityGraphEnv) -> torch.Tensor:
    """Per-edge binary targets from the spatial–momentum baseline on the current env event.

    Requires ``env.reset`` to have been called for this event so baseline fields are populated.
    """
    return edge_pair_baseline_targets(
        env._baseline_node_labels,
        env.graph.edge_pair_i.numpy(),
        env.graph.edge_pair_j.numpy(),
    )


def weighted_bce_with_logits(
    logits: torch.Tensor,
    targets: torch.Tensor,
    *,
    auto_pos_weight: bool = True,
    pos_weight: float | None = None,
    pos_weight_power: float = 0.5,
    max_pos_weight: float = 300.0,
    focal_gamma: float = 0.0,
) -> tuple[torch.Tensor, float]:
    """BCE over edges with optional positive-class reweighting and focal modulation."""
    if pos_weight is not None:
        pw = float(max(1.0, pos_weight))
    elif auto_pos_weight:
        pos = float(targets.sum().item())
        neg = float((1.0 - targets).sum().item())
        ratio = neg / max(pos, 1.0)
        pw = ratio**float(max(pos_weight_power, 0.0))
    else:
        pw = 1.0
    pw = float(np.clip(pw, 1.0, max_pos_weight))
    weights = 1.0 + (pw - 1.0) * targets
    bce = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
    if focal_gamma > 0.0:
        p = torch.sigmoid(logits)
        pt = torch.where(targets > 0.5, p, 1.0 - p).clamp(1e-6, 1.0 - 1e-6)
        bce = bce * (1.0 - pt) ** float(focal_gamma)
    loss = (bce * weights).sum() / torch.clamp(weights.sum(), min=1.0)
    return loss, pw



type EventSampler = Callable[[], tuple[np.ndarray, np.ndarray, np.ndarray]]


def train_supervised_edges(
    policy: GATAffinityPolicy,
    env: AffinityGraphEnv,
    event_sampler: EventSampler,
    *,
    steps: int,
    events_per_step: int = 8,
    lr: float = 3e-3,
    max_grad_norm: float = 0.5,
    weighted_bce: bool = True,
    pos_weight: float | None = None,
    pos_weight_power: float = 0.5,
    pos_weight_max: float = 300.0,
    focal_gamma: float = 0.0,
    on_update: Callable[[dict[str, list]], None] | None = None,
    optimizer: torch.optim.Optimizer | None = None,
    lr_scheduler: LRScheduler | None = None,
) -> dict[str, list]:
    """BCE on edge logits against spatial–momentum baseline edge targets (warm-start / standalone).

    If ``lr_scheduler`` is set, ``lr_scheduler.step()`` runs after each ``optimizer.step()``.
    """
    opt = optimizer or torch.optim.Adam(policy.parameters(), lr=lr)
    history: dict[str, list] = {
        "supervised_bce": [],
        "pretrain_partition_loss": [],
        "pretrain_baseline_loss": [],
        "pretrain_gap": [],
        "supervised_pos_weight": [],
    }
    policy.train()
    pbar_sup = tqdm(
        range(steps),
        desc="SupEdges",
        miniters=1,
        mininterval=0.0,
        dynamic_ncols=True,
    )
    for _ in pbar_sup:
        opt.zero_grad()
        sup_acc = 0.0
        sup_part: list[float] = []
        sup_base: list[float] = []
        n_eff = 0
        for _ in range(events_per_step):
            pos, mom, isp = event_sampler()
            obs = env.reset(pos, mom, isp)
            logits = policy(obs)
            edge_on = (torch.sigmoid(logits) > 0.5).float()
            l_pol, _ = env.physics_for_edge_mask(edge_on)
            sup_part.append(float(l_pol))
            sup_base.append(float(env._baseline_loss))
            tgt = baseline_edge_targets(env)
            loss_sup, pw = weighted_bce_with_logits(
                logits,
                tgt,
                auto_pos_weight=weighted_bce and pos_weight is None,
                pos_weight=pos_weight,
                pos_weight_power=pos_weight_power,
                max_pos_weight=pos_weight_max,
                focal_gamma=focal_gamma,
            )
            loss_sup.backward()
            sup_acc += float(loss_sup.item())
            history["supervised_pos_weight"].append(float(pw))
            n_eff += 1
        nn.utils.clip_grad_norm_(policy.parameters(), max_grad_norm)
        opt.step()
        if lr_scheduler is not None:
            lr_scheduler.step()
        mean_sup = sup_acc / max(float(n_eff), 1.0)
        history["supervised_bce"].append(mean_sup)
        if sup_part:
            history["pretrain_partition_loss"].append(float(np.mean(sup_part)))
        if sup_base:
            history["pretrain_baseline_loss"].append(float(np.mean(sup_base)))
        if sup_part and sup_base:
            history["pretrain_gap"].append(float(np.mean(sup_part)) - float(np.mean(sup_base)))
        pf_sup: dict[str, float] = {"bce": mean_sup}
        if sup_part:
            pf_sup["L_pol"] = float(np.mean(sup_part)) / MEV_PER_GEV
        if sup_base:
            pf_sup["L_base"] = float(np.mean(sup_base)) / MEV_PER_GEV
        if sup_part and sup_base:
            pf_sup["gap"] = (float(np.mean(sup_part)) - float(np.mean(sup_base))) / MEV_PER_GEV
        pbar_sup.set_postfix(pf_sup, refresh=True)
        if on_update is not None:
            on_update(history)
    return history


def train_reinforce(
    policy: GATAffinityPolicy,
    env: AffinityGraphEnv,
    event_sampler: EventSampler,
    *,
    optimizer: torch.optim.Optimizer,
    n_updates: int = 150,
    episodes_per_update: int = 8,
    lr_scheduler: LRScheduler | None = None,
    ent_coef: float = 0.02,
    max_grad_norm: float = 0.5,
    on_update: Callable[[dict[str, list]], None] | None = None,
    policy_coef: float = 1.0,
) -> dict[str, list]:
    """REINFORCE on edge Bernoulli actions.

    Each update uses batch-centered returns ``R_i - mean(R)`` over episodes collected
    in that update (REINFORCE baseline with no cross-batch state).

    ``lr_scheduler`` is optional; if given, ``step()`` is called after each ``optimizer.step()``.
    """
    opt = optimizer
    history: dict[str, list] = {
        "episode_return": [],
        "partition_loss": [],
        "baseline_loss": [],
        "policy_loss": [],
        "return_baseline": [],
        "edge_entropy": [],
        "n_clusters": [],
        "lr": [],
    }
    pbar = tqdm(range(n_updates), desc="REINFORCE", miniters=1, mininterval=0.0, dynamic_ncols=True)
    for u in pbar:
        ep_returns: list[float] = []
        part_losses: list[float] = []
        base_losses: list[float] = []
        n_clust: list[int] = []
        batch_obs: list[Data] = []
        batch_act: list[torch.Tensor] = []

        for _ in range(episodes_per_update):
            pos, mom, isp = event_sampler()
            ep = collect_rollout(env, policy, pos, mom, isp)
            ep_returns.append(float(ep["reward"]))
            part_losses.append(float(ep["partition_loss"]))
            base_losses.append(float(ep["baseline_loss"]))
            n_clust.append(ep["n_clusters"])
            batch_obs.append(ep["obs"])
            batch_act.append(ep["action"])

        mean_r = float(np.mean(ep_returns)) if ep_returns else 0.0
        history["return_baseline"].append(mean_r)

        policy.train()
        opt.zero_grad()
        pol_acc = ent_acc = 0.0
        n_b = max(len(batch_obs), 1)
        for i in range(len(batch_obs)):
            obs = batch_obs[i]
            logits = policy(obs)
            dist = torch.distributions.Bernoulli(logits=logits)
            a = batch_act[i]
            logp = dist.log_prob(a).mean()
            ent = dist.entropy().mean()
            g = ep_returns[i]
            adv_t = torch.tensor(g - mean_r, dtype=torch.float32)
            pol_t = -(logp * adv_t)
            loss_t = (policy_coef * pol_t - ent_coef * ent) / n_b
            loss_t.backward()
            pol_acc += pol_t.item()
            ent_acc += ent.item()
        nn.utils.clip_grad_norm_(policy.parameters(), max_grad_norm)
        opt.step()
        if lr_scheduler is not None:
            lr_scheduler.step()
        history["lr"].append(float(opt.param_groups[0]["lr"]))

        if ep_returns:
            history["episode_return"].append(float(np.mean(ep_returns)))
        if part_losses:
            history["partition_loss"].append(float(np.mean(part_losses)))
        if base_losses and all(np.isfinite(base_losses)):
            history["baseline_loss"].append(float(np.mean(base_losses)))
        if n_clust:
            history["n_clusters"].append(float(np.mean(n_clust)))
        history["policy_loss"].append(pol_acc / n_b)
        history["edge_entropy"].append(ent_acc / n_b)

        pf: dict[str, float] = {
            "pi": float(pol_acc / n_b),
            "H": float(ent_acc / n_b),
            "lr": float(opt.param_groups[0]["lr"]),
            "Rmean": mean_r / MEV_PER_GEV,
        }
        if ep_returns:
            pf["G"] = float(np.mean(ep_returns)) / MEV_PER_GEV
        if part_losses:
            pf["L_pol"] = float(np.mean(part_losses)) / MEV_PER_GEV
        if base_losses and all(np.isfinite(base_losses)):
            pf["L_base"] = float(np.mean(base_losses)) / MEV_PER_GEV
        if n_clust:
            pf["n_cl"] = float(np.mean(n_clust))
        pbar.set_postfix(pf, refresh=True)
        if on_update is not None:
            on_update(history)
    return history


def make_event_sampler(
    events: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    rng: np.random.Generator,
    fallback_urqmd: Callable[[], tuple[np.ndarray, np.ndarray, np.ndarray]] | None = None,
) -> EventSampler:
    """Build event sampler from preloaded events with optional URQMD fallback."""

    def sample_event() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        if events:
            pos, mom, isp = events[int(rng.integers(0, len(events)))]
            return pos.copy(), mom.copy(), isp.copy()
        pos, mom, isp = fallback_urqmd()
        return pos, mom, isp

    return sample_event


def score_history(history: dict[str, list], tail: int = 10) -> float:
    """Lower is better: mean (L_pol - L_base) over tail updates."""
    l_pol = np.asarray(history.get("partition_loss", []), dtype=np.float64)
    l_base = np.asarray(history.get("baseline_loss", []), dtype=np.float64)
    n = min(len(l_pol), len(l_base))
    if n == 0:
        return float("inf")
    n_tail = min(int(tail), n)
    gap = l_pol[n - n_tail :] - l_base[n - n_tail :]
    return float(np.mean(gap))
