"""
Cluster energy utilities.

1. **UrQMD-inspired model + macroscopic binding prior** on a :class:`NucleonCloud`
   (positions fm, :class:`~pylorentz.Momentum4` in GeV, ``is_proton``). Partition
   score: :func:`partition_loss`.

2. **Excitation from tabulated mass** (MeV): :func:`cluster_excitation_energy_mev`.
"""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

_MASSTABLE_FROM_FILE_PATCHED = False


def _patch_masstable_table_from_file() -> None:
    """
    ``masstable.Table.from_file`` uses ``pd.read_csv(..., delim_whitespace=True)``,
    which pandas deprecates in favor of ``sep=r'\\s+'``.
    """
    global _MASSTABLE_FROM_FILE_PATCHED
    if _MASSTABLE_FROM_FILE_PATCHED:
        return
    import pandas as pd

    import masstable.masstable as msm

    @classmethod
    def from_file(cls, filename: str, name: str = ""):
        df = pd.read_csv(
            filename,
            header=0,
            sep=r"\s+",
            engine="python",
            index_col=[0, 1],
        )["M"]
        df.name = name
        return cls(df=df, name=name)

    msm.Table.from_file = from_file
    _MASSTABLE_FROM_FILE_PATCHED = True
from pylorentz import Momentum4

# -----------------------------------------------------------------------------
# Nucleon cloud (replaces notebook ``Event`` / ``Particle``)
# -----------------------------------------------------------------------------


@dataclass
class NucleonCloud:
    """
    One event worth of nucleons: Cartesian positions, four-momenta, and species.

    ``pos`` shape ``(N, 3)`` in fm. ``four_momentum`` length ``N``: each
    :class:`~pylorentz.Momentum4` in GeV (``e``, ``p_x``, ``p_y``, ``p_z``).
    ``is_proton`` length ``N`` (neutron if False).
    """

    pos: np.ndarray
    four_momentum: tuple[Momentum4, ...]
    is_proton: np.ndarray

    def __post_init__(self) -> None:
        self.pos = np.asarray(self.pos, dtype=np.float64)
        self.is_proton = np.asarray(self.is_proton, dtype=bool)
        n = int(self.pos.shape[0])
        if self.pos.ndim != 2 or self.pos.shape[1] != 3:
            raise ValueError("pos must have shape (N, 3)")
        if len(self.four_momentum) != n:
            raise ValueError("four_momentum length must match pos rows")
        if self.is_proton.shape != (n,):
            raise ValueError("is_proton must have shape (N,)")

    @classmethod
    def from_numpy_mev(
        cls,
        pos: np.ndarray,
        mom_mev: np.ndarray,
        is_proton: np.ndarray,
    ) -> NucleonCloud:
        """
        Build from dataset-style arrays: ``pos`` is ``(N, 3)`` ``(x, y, z)`` fm or ``(N, 4)``
        ``(t, x, y, z)`` with ``t`` in fm/c and ``x,y,z`` in fm (spatial columns are used).
        ``mom_mev`` is ``(N, 4)`` as ``(E, px, py, pz)`` in MeV (stored internally as GeV ``Momentum4``).
        """
        pos = np.asarray(pos, dtype=np.float64)
        mom_mev = np.asarray(mom_mev, dtype=np.float64)
        is_proton = np.asarray(is_proton, dtype=bool)
        if pos.ndim != 2 or pos.shape[1] not in (3, 4):
            raise ValueError("pos must have shape (N, 3) or (N, 4)")
        if mom_mev.ndim != 2 or mom_mev.shape[1] != 4:
            raise ValueError("mom_mev must have shape (N, 4)")
        n = pos.shape[0]
        if mom_mev.shape[0] != n or is_proton.shape[0] != n:
            raise ValueError("pos, mom_mev, is_proton must have the same length")
        pos3 = pos[:, 1:4].copy() if pos.shape[1] == 4 else pos.copy()
        s = 1.0 / 1000.0
        p4: list[Momentum4] = []
        for i in range(n):
            e, px, py, pz = mom_mev[i]
            p4.append(Momentum4(e * s, px * s, py * s, pz * s))
        return cls(pos=pos3, four_momentum=tuple(p4), is_proton=is_proton)

    def summed_momentum(self, indices: Sequence[int]) -> Momentum4:
        """Lorentz sum of ``four_momentum[i]`` for ``i`` in ``indices``."""
        idx = [int(i) for i in indices]
        if not idx:
            return Momentum4(0.0, 0.0, 0.0, 0.0)
        total = self.four_momentum[idx[0]]
        for ii in idx[1:]:
            total = total + self.four_momentum[ii]
        return total

    def four_momentum_numpy_mev(self) -> np.ndarray:
        """Stack ``(E, px, py, pz)`` per row in MeV (for kinematic helpers)."""
        n = len(self.four_momentum)
        out = np.empty((n, 4), dtype=np.float64)
        for i, p in enumerate(self.four_momentum):
            out[i, 0] = p.e * 1000.0
            out[i, 1] = p.p_x * 1000.0
            out[i, 2] = p.p_y * 1000.0
            out[i, 3] = p.p_z * 1000.0
        return out


@dataclass(frozen=True)
class ClusterEnergyResult:
    """UrQMD-style cluster decomposition."""

    A: int
    Z: int
    internal_kinetic: float
    pair_potential: float
    binding_prior: float
    total_energy: float


def _p_vec(cloud: NucleonCloud, i: int) -> np.ndarray:
    p = cloud.four_momentum[i]
    return np.array([p.p_x, p.p_y, p.p_z], dtype=float)


def _mass_gev(cloud: NucleonCloud, i: int) -> float:
    return float(cloud.four_momentum[i].m)


def spatial_distance(cloud: NucleonCloud, i: int, j: int) -> float:
    return float(np.linalg.norm(cloud.pos[i] - cloud.pos[j]))


# -----------------------------------------------------------------------------
# Nuclear binding for the macroscopic prior (MeV → GeV): tables first, SEMF fallback
# -----------------------------------------------------------------------------

M_P = 938.2723
M_N = 939.5656

_MASS_TABLE_NAMES: tuple[str, ...] = ("AME2012all",)


@functools.lru_cache(maxsize=None)
def _table(name: str):
    _patch_masstable_table_from_file()
    from masstable import Table

    return Table(name)


def binding_energy_liquid_drop(a: int, z: int) -> float:
    """SEMF binding energy in MeV (used only when nuclide is absent from mass tables)."""
    if a <= 0 or z < 0 or z > a:
        raise ValueError(f"invalid (A, Z)=({a}, {z})")
    n = a - z
    if a == 1:
        return 0.0
    a_v, a_s, a_c, a_a = 15.75, 17.8, 0.711, 23.7
    vol = a_v * a
    surf = a_s * (a ** (2.0 / 3.0))
    coul = a_c * z * (z - 1.0) / (a ** (1.0 / 3.0))
    asym = a_a * ((n - z) ** 2) / a
    if a % 2 == 1:
        delta = 0.0
    elif z % 2 == 0:
        delta = 12.0 / math.sqrt(a)
    else:
        delta = -12.0 / math.sqrt(a)
    return vol - surf - coul - asym + delta


def binding_energy_from_tables(a: int, z: int) -> float | None:
    """Binding energy in MeV from ``masstable`` if ``(Z, N)`` is in the table."""
    n = a - z
    key = (z, n)
    for name in _MASS_TABLE_NAMES:
        t = _table(name)
        if key in t.binding_energy.df.index:
            return float(t.binding_energy.df.loc[key])
    return None


def get_mass_mev(a: int, z: int) -> float:
    """
    Ground-state mass in MeV: ``Z M_p + N M_n - B`` with ``B`` from tables if present,
    else :func:`binding_energy_liquid_drop`.
    """
    n = a - z
    b = binding_energy_from_tables(a, z)
    if b is None:
        b = binding_energy_liquid_drop(a, z)
    return z * M_P + n * M_N - b


def nuclear_binding_energy_mev(a: int, z: int) -> float:
    """Positive binding ``B`` in MeV: ``Z M_p + N M_n - M_gs``."""
    return z * M_P + (a - z) * M_N - get_mass_mev(a, z)


def binding_prior_gev(a: int, z: int) -> float:
    """
    Binding energy in **GeV** for subtraction from the UrQMD cluster energy (tables
    when available, else SEMF), divided by 1000.
    """
    if a < 1 or z < 0 or z > a:
        return 0.0
    return nuclear_binding_energy_mev(a, z) / 1000.0


# -----------------------------------------------------------------------------
# UrQMD hard EoS parameters (Table 3.1, no Pauli by default)
# -----------------------------------------------------------------------------

URQMD_ALPHA_FM_INV2 = 0.25
URQMD_T1_GEV_FM3 = -7264.04 / 1000.0
URQMD_TGAMMA_GEV_FM6 = 87.65 / 1000.0
URQMD_V0YUK_GEV_FM = -0.498 / 1000.0
URQMD_GAMMA_Y_FM = 1.4
URQMD_E2_GEV_FM = 1.44 / 1000.0

URQMD_P_REL_MAX = 2.0

URQMD_USE_PAULI_APPROX = False
URQMD_V0PAU_GEV = 98.95 / 1000.0
URQMD_Q0_FM = 2.16
URQMD_P0_GEV = 120.0 / 1000.0

URQMD_USE_SK3 = False


def _erf_vec(x: np.ndarray) -> np.ndarray:
    """Vector error function (Abramowitz & Stegun 7.1.26), no SciPy."""
    x = np.asarray(x, dtype=np.float64)
    sign = np.sign(x)
    ax = np.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * ax)
    poly = (
        (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592)
        * t
    )
    y = 1.0 - poly * np.exp(-ax * ax)
    return sign * y


def pair_relative_momentum(cloud: NucleonCloud, i: int, j: int) -> float:
    pi = _p_vec(cloud, i)
    pj = _p_vec(cloud, j)
    return float(np.linalg.norm(pi - pj))


def cluster_internal_kinetic(cloud: NucleonCloud, cluster: list[int]) -> float:
    if len(cluster) <= 1:
        return 0.0
    masses = np.array([_mass_gev(cloud, i) for i in cluster], dtype=float)
    momenta = np.array([_p_vec(cloud, i) for i in cluster], dtype=float)
    total_mass = float(np.sum(masses))
    p_cm = np.sum(momenta, axis=0) / total_mass
    q = momenta - masses[:, None] * p_cm[None, :]
    q2 = np.sum(q * q, axis=1)
    e_rel = np.sqrt(masses * masses + q2) - masses
    return float(np.sum(e_rel))


def urqmd_sk2_pair_energy(r: float) -> float:
    alpha = URQMD_ALPHA_FM_INV2
    return URQMD_T1_GEV_FM3 * (alpha / math.pi) ** 1.5 * math.exp(-alpha * r * r)


def _urqmd_sk2_pair_energy_vec(r: np.ndarray) -> np.ndarray:
    alpha = URQMD_ALPHA_FM_INV2
    return URQMD_T1_GEV_FM3 * (alpha / math.pi) ** 1.5 * np.exp(-alpha * r * r)


def urqmd_yukawa_pair_energy(r: float) -> float:
    if r <= 1.0e-12:
        return 0.0
    alpha = URQMD_ALPHA_FM_INV2
    gamma_y = URQMD_GAMMA_Y_FM
    v0 = URQMD_V0YUK_GEV_FM
    pref = v0 * (1.0 / (2.0 * r)) * math.exp(1.0 / (4.0 * alpha * gamma_y * gamma_y))
    a = 1.0 / (2.0 * gamma_y * math.sqrt(alpha))
    b = math.sqrt(alpha) * r
    term1 = math.exp(-r / gamma_y) * (1.0 - math.erf(a - b))
    term2 = math.exp(+r / gamma_y) * (1.0 - math.erf(a + b))
    return pref * (term1 - term2)


def _urqmd_yukawa_pair_energy_vec(r: np.ndarray) -> np.ndarray:
    out = np.zeros_like(r, dtype=np.float64)
    mask = r > 1.0e-12
    if not np.any(mask):
        return out
    r_m = r[mask]
    alpha = URQMD_ALPHA_FM_INV2
    gamma_y = URQMD_GAMMA_Y_FM
    v0 = URQMD_V0YUK_GEV_FM
    pref = v0 * (1.0 / (2.0 * r_m)) * np.exp(1.0 / (4.0 * alpha * gamma_y * gamma_y))
    a = 1.0 / (2.0 * gamma_y * math.sqrt(alpha))
    b = math.sqrt(alpha) * r_m
    term1 = np.exp(-r_m / gamma_y) * (1.0 - _erf_vec(a - b))
    term2 = np.exp(+r_m / gamma_y) * (1.0 - _erf_vec(a + b))
    out[mask] = pref * (term1 - term2)
    return out


def urqmd_coulomb_pair_energy(cloud: NucleonCloud, i: int, j: int, r: float) -> float:
    if r <= 1.0e-12:
        return 0.0
    zi = 1 if cloud.is_proton[i] else 0
    zj = 1 if cloud.is_proton[j] else 0
    if zi == 0 or zj == 0:
        return 0.0
    alpha = URQMD_ALPHA_FM_INV2
    return (zi * zj * URQMD_E2_GEV_FM / r) * math.erf(math.sqrt(alpha) * r)


def _urqmd_coulomb_pair_energy_vec(z_prod: np.ndarray, r: np.ndarray) -> np.ndarray:
    """``z_prod`` = ``z_i * z_j`` (0 for neutron pairs)."""
    out = np.zeros_like(r, dtype=np.float64)
    mask = (r > 1.0e-12) & (z_prod > 0.0)
    if not np.any(mask):
        return out
    r_m = r[mask]
    z_m = z_prod[mask]
    alpha = URQMD_ALPHA_FM_INV2
    out[mask] = (z_m * URQMD_E2_GEV_FM / r_m) * _erf_vec(np.sqrt(alpha) * r_m)
    return out


def urqmd_pauli_pair_energy_approx(
    cloud: NucleonCloud, i: int, j: int, r: float, p_rel: float
) -> float:
    if not URQMD_USE_PAULI_APPROX:
        return 0.0
    same_isospin = int(bool(cloud.is_proton[i]) == bool(cloud.is_proton[j]))
    if same_isospin == 0:
        return 0.0
    alpha = URQMD_ALPHA_FM_INV2
    q0 = URQMD_Q0_FM
    p0 = URQMD_P0_GEV
    pref = (
        URQMD_V0PAU_GEV
        * (1.0 / (p0 * q0)) ** 3
        * (1.0 + 1.0 / (2.0 * alpha * q0 * q0)) ** (-1.5)
    )
    expo = math.exp(
        -alpha * r * r / (2.0 * alpha * q0 * q0 + 1.0) - (p_rel * p_rel) / (2.0 * p0 * p0)
    )
    return pref * expo * same_isospin


def _urqmd_pauli_pair_energy_vec(
    r: np.ndarray, p_rel: np.ndarray, same_isospin: np.ndarray
) -> np.ndarray:
    if not URQMD_USE_PAULI_APPROX:
        return np.zeros_like(r, dtype=np.float64)
    alpha = URQMD_ALPHA_FM_INV2
    q0 = URQMD_Q0_FM
    p0 = URQMD_P0_GEV
    pref = (
        URQMD_V0PAU_GEV
        * (1.0 / (p0 * q0)) ** 3
        * (1.0 + 1.0 / (2.0 * alpha * q0 * q0)) ** (-1.5)
    )
    expo = np.exp(
        -alpha * r * r / (2.0 * alpha * q0 * q0 + 1.0) - (p_rel * p_rel) / (2.0 * p0 * p0)
    )
    return pref * expo * same_isospin


def urqmd_pair_energy(cloud: NucleonCloud, i: int, j: int) -> float:
    p_rel = pair_relative_momentum(cloud, i, j)
    if p_rel >= URQMD_P_REL_MAX:
        return 0.0
    r = spatial_distance(cloud, i, j)
    e_sk2 = urqmd_sk2_pair_energy(r)
    e_yuk = urqmd_yukawa_pair_energy(r)
    e_coul = urqmd_coulomb_pair_energy(cloud, i, j, r)
    e_pau = urqmd_pauli_pair_energy_approx(cloud, i, j, r, p_rel)
    return e_sk2 + e_yuk + e_coul + e_pau


def urqmd_pair_energies_upper_triangle(
    pos_c: np.ndarray, mom3: np.ndarray, is_proton_c: np.ndarray
) -> np.ndarray:
    """
    Pair energies for all unique pairs in a cluster (upper triangle), shape ``(P,)``
    with ``P = n*(n-1)//2``, same order as ``np.triu_indices(n, k=1)``.
    """
    n = pos_c.shape[0]
    if n <= 1:
        return np.zeros((0,), dtype=np.float64)
    dr = pos_c[:, None, :] - pos_c[None, :, :]
    r = np.linalg.norm(dr, axis=2)
    dp = mom3[:, None, :] - mom3[None, :, :]
    p_rel = np.linalg.norm(dp, axis=2)
    iu, ju = np.triu_indices(n, k=1)
    r_p = r[iu, ju]
    p_rel_p = p_rel[iu, ju]
    zl = is_proton_c.astype(np.float64)
    z_prod = zl[iu] * zl[ju]
    same_iso = (is_proton_c[iu] == is_proton_c[ju]).astype(np.float64)

    e_sk2 = _urqmd_sk2_pair_energy_vec(r_p)
    e_yuk = _urqmd_yukawa_pair_energy_vec(r_p)
    e_coul = _urqmd_coulomb_pair_energy_vec(z_prod, r_p)
    e_pau = _urqmd_pauli_pair_energy_vec(r_p, p_rel_p, same_iso)
    e = e_sk2 + e_yuk + e_coul + e_pau
    return np.where(p_rel_p < URQMD_P_REL_MAX, e, 0.0)


def urqmd_sk3_triplet_energy(cloud: NucleonCloud, i: int, j: int, k: int) -> float:
    if not URQMD_USE_SK3:
        return 0.0
    pij = pair_relative_momentum(cloud, i, j)
    pik = pair_relative_momentum(cloud, i, k)
    pjk = pair_relative_momentum(cloud, j, k)
    if max(pij, pik, pjk) >= URQMD_P_REL_MAX:
        return 0.0
    alpha = URQMD_ALPHA_FM_INV2
    pref = URQMD_TGAMMA_GEV_FM6 * (4.0 * alpha * alpha / (3.0 * math.pi * math.pi)) ** 1.5
    rij = spatial_distance(cloud, i, j)
    rik = spatial_distance(cloud, i, k)
    rjk = spatial_distance(cloud, j, k)
    term_i = math.exp(-alpha * (rij * rij + rik * rik))
    term_j = math.exp(-alpha * (rij * rij + rjk * rjk))
    term_k = math.exp(-alpha * (rik * rik + rjk * rjk))
    return pref * (term_i + term_j + term_k) / 3.0


def cluster_AZ(cloud: NucleonCloud, cluster: list[int]) -> tuple[int, int]:
    a = len(cluster)
    z = sum(1 for i in cluster if cloud.is_proton[i])
    return a, z


def cluster_pair_energy(cloud: NucleonCloud, cluster: list[int]) -> float:
    n = len(cluster)
    if n <= 1:
        return 0.0
    idx = np.asarray(cluster, dtype=int)
    pos_c = cloud.pos[idx]
    mom3 = np.empty((n, 3), dtype=np.float64)
    isp = cloud.is_proton[idx]
    for k in range(n):
        p = cloud.four_momentum[int(idx[k])]
        mom3[k, 0] = p.p_x
        mom3[k, 1] = p.p_y
        mom3[k, 2] = p.p_z
    return float(np.sum(urqmd_pair_energies_upper_triangle(pos_c, mom3, isp)))


def cluster_triplet_energy(cloud: NucleonCloud, cluster: list[int]) -> float:
    if len(cluster) <= 2:
        return 0.0
    e = 0.0
    for a in range(len(cluster)):
        ia = cluster[a]
        for b in range(a + 1, len(cluster)):
            ib = cluster[b]
            for c in range(b + 1, len(cluster)):
                ic = cluster[c]
                e += urqmd_sk3_triplet_energy(cloud, ia, ib, ic)
    return e


def cluster_energy(cloud: NucleonCloud, cluster: list[int]) -> ClusterEnergyResult:
    """UrQMD cluster energy with macroscopic binding :func:`binding_prior_gev`."""
    a, z = cluster_AZ(cloud, cluster)
    tint = cluster_internal_kinetic(cloud, cluster)
    vpair = cluster_pair_energy(cloud, cluster)
    vtrip = cluster_triplet_energy(cloud, cluster)
    urqmd_total = tint + vpair + vtrip
    b_prior = binding_prior_gev(a, z) if a >= 1 else 0.0
    total = urqmd_total - b_prior
    return ClusterEnergyResult(
        A=a,
        Z=z,
        internal_kinetic=tint,
        pair_potential=vpair + vtrip,
        binding_prior=b_prior,
        total_energy=total,
    )


def partition_loss(cloud: NucleonCloud, partition: list[list[int]]) -> float:
    """Sum of :func:`cluster_energy` ``total_energy`` over clusters."""
    total = 0.0
    for c in partition:
        total += cluster_energy(cloud, c).total_energy
    return total


def partition_loss_numpy(
    pos: np.ndarray,
    mom_mev: np.ndarray,
    is_proton: np.ndarray,
    partition: list[list[int]],
) -> float:
    """:func:`partition_loss` on numpy arrays.

    ``pos`` is ``(N, 3)`` fm ``(x,y,z)`` or ``(N, 4)`` ``(t, x, y, z)`` (fm/c, fm). Four-momenta
    rows are MeV/c; **return value is GeV**.
    """
    return partition_loss(NucleonCloud.from_numpy_mev(pos, mom_mev, is_proton), partition)


def _invariant_mass_mev(mom_stack: np.ndarray) -> float:
    total = np.sum(np.asarray(mom_stack, dtype=np.float64), axis=0)
    e, px, py, pz = total[0], total[1], total[2], total[3]
    s = e * e - px * px - py * py - pz * pz
    if s < 0.0 and s > -1e-3 * max(1.0, e * e):
        s = 0.0
    if s < 0.0:
        raise ValueError(f"negative invariant Mandelstam s={s} from summed four-momentum")
    return float(math.sqrt(s))


@dataclass(frozen=True)
class ClusterExcitationResult:
    """Breakdown for :func:`cluster_excitation_energy_mev`."""

    a: int
    z: int
    invariant_mass_mev: float
    ground_state_mass_mev: float
    excitation_energy_mev: float


def cluster_excitation_energy_mev(
    mom_stack: np.ndarray,
    is_proton: np.ndarray,
    *,
    return_details: bool = False,
) -> float | ClusterExcitationResult:
    """Invariant cluster mass minus tabulated ground-state mass (MeV)."""
    mom_stack = np.asarray(mom_stack, dtype=np.float64)
    is_proton = np.asarray(is_proton, dtype=bool)
    if mom_stack.ndim != 2 or mom_stack.shape[1] != 4:
        raise ValueError("mom_stack must have shape (N, 4) with (E, px, py, pz)")
    if is_proton.shape[0] != mom_stack.shape[0]:
        raise ValueError("is_proton length must match number of rows in mom_stack")

    n_nuc = int(mom_stack.shape[0])
    if n_nuc == 0:
        if return_details:
            return ClusterExcitationResult(
                a=0,
                z=0,
                invariant_mass_mev=0.0,
                ground_state_mass_mev=0.0,
                excitation_energy_mev=0.0,
            )
        return 0.0

    if n_nuc == 1:
        z = int(np.sum(is_proton))
        m_inv = _invariant_mass_mev(mom_stack)
        m_gs = get_mass_mev(1, z)
        if return_details:
            return ClusterExcitationResult(
                a=1,
                z=z,
                invariant_mass_mev=m_inv,
                ground_state_mass_mev=m_gs,
                excitation_energy_mev=0.0,
            )
        return 0.0

    a = n_nuc
    z = int(np.sum(is_proton))
    m_inv = _invariant_mass_mev(mom_stack)
    m_gs = get_mass_mev(a, z)
    delta = m_inv - m_gs
    if return_details:
        return ClusterExcitationResult(
            a=a,
            z=z,
            invariant_mass_mev=m_inv,
            ground_state_mass_mev=m_gs,
            excitation_energy_mev=delta,
        )
    return float(delta)


def cluster_excitation_energy_mev_for_indices(
    mom: np.ndarray,
    is_proton: np.ndarray,
    indices: Sequence[int],
    *,
    return_details: bool = False,
) -> float | ClusterExcitationResult:
    """Same as :func:`cluster_excitation_energy_mev` for a subset by index."""
    idx = np.asarray(indices, dtype=int)
    return cluster_excitation_energy_mev(
        mom[idx], is_proton[idx], return_details=return_details
    )


def cluster_excitation_energy_mev_from_cloud(
    cloud: NucleonCloud,
    indices: Sequence[int],
    *,
    return_details: bool = False,
) -> float | ClusterExcitationResult:
    """Excitation energy using Lorentz sum of ``cloud.four_momentum`` (MeV output)."""
    idx = [int(i) for i in indices]
    if len(idx) == 0:
        if return_details:
            return ClusterExcitationResult(
                a=0,
                z=0,
                invariant_mass_mev=0.0,
                ground_state_mass_mev=0.0,
                excitation_energy_mev=0.0,
            )
        return 0.0
    if len(idx) == 1:
        z = int(cloud.is_proton[idx[0]])
        m_inv = float(cloud.four_momentum[idx[0]].m) * 1000.0
        m_gs = get_mass_mev(1, z)
        if return_details:
            return ClusterExcitationResult(
                a=1,
                z=z,
                invariant_mass_mev=m_inv,
                ground_state_mass_mev=m_gs,
                excitation_energy_mev=0.0,
            )
        return 0.0
    tot = cloud.summed_momentum(idx)
    m_inv = float(tot.m) * 1000.0
    isp = cloud.is_proton[np.asarray(idx, dtype=int)]
    a = len(idx)
    z = int(np.sum(isp))
    m_gs = get_mass_mev(a, z)
    delta = m_inv - m_gs
    if return_details:
        return ClusterExcitationResult(
            a=a,
            z=z,
            invariant_mass_mev=m_inv,
            ground_state_mass_mev=m_gs,
            excitation_energy_mev=delta,
        )
    return float(delta)


def partition_excitation_loss_mev(
    mom_mev: np.ndarray,
    is_proton: np.ndarray,
    partition: list[list[int]],
) -> float:
    """Sum of :func:`cluster_excitation_energy_mev` over clusters."""
    total = 0.0
    for c in partition:
        idx = np.asarray(c, dtype=int)
        total += float(cluster_excitation_energy_mev(mom_mev[idx], is_proton[idx]))
    return total


def partition_excitation_loss_mev_from_cloud(
    cloud: NucleonCloud,
    partition: list[list[int]],
) -> float:
    """Sum of :func:`cluster_excitation_energy_mev_from_cloud` over clusters."""
    total = 0.0
    for c in partition:
        total += float(cluster_excitation_energy_mev_from_cloud(cloud, c))
    return total
