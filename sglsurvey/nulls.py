"""Null-ensemble statistics for trajectory stacks (WISE v2 plan §3;
promoted from ``surveys/ztf/scripts/look_elsewhere.py``).

Vocabulary
----------
A *cell* is one (endpoint, role, band) search. Its *trajectories* are
the real (z, mu) grid plus the null constructions:

* spatial offsets — the same cutouts sampled on a ring of offsets from
  the locus (8 *designated* controls define the local threshold T, the
  rest are extra null trials);
* time scrambles — the real-trajectory tensor with the grid-node
  labels permuted independently per epoch, so each epoch keeps its own
  noise, cadence and phase but the trajectory coherence is broken;
* trajectory randomisation — the cutouts sampled on another endpoint's
  (z, mu) sky-motion pattern re-centred on this corridor.

Per trajectory the *grid maximum* S_max is the maximum over grid nodes
of the inverse-variance-weighted stack significance S = A / sqrt(B)
(per-epoch weights capped at WEIGHT_CAP x the per-node median, the v1
"effective-epoch floor"). The *exceedance ratio* is R = S_max / T where
T is the maximum grid maximum over the designated controls (leave one
out when the scored trajectory is itself a designated control). Under
the null, R of the real trajectory is exchangeable with R of every
null trajectory; the survey-wide statistic is max R over the cells of
a family, calibrated by pseudo-experiments that draw one pooled null R
per cell.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

MIN_GOOD_FRAC = 0.7
WEIGHT_CAP = 20.0
MIN_EPOCHS = 5


# -- stacking -----------------------------------------------------------
def stack_weights(v: np.ndarray, g: np.ndarray, epoch_ok=None,
                  f: np.ndarray | None = None, clip_sigma: float | None = None):
    """Capped inverse-variance weights (E, *grid) and the validity mask.

    ``epoch_ok`` (E,) boolean drops whole epochs (quality masks, phase
    or time selections) before the per-node median cap is computed.
    ``clip_sigma`` (with ``f``) drops single-epoch samples with
    |f|/sqrt(v) above it (the ZTF single-epoch clip); WISE v2 does not
    clip.
    """
    v = np.asarray(v, dtype=np.float32)
    g = np.asarray(g, dtype=np.float32)
    valid = np.isfinite(v) & (v > 0) & (g >= MIN_GOOD_FRAC)
    if f is not None:
        valid &= np.isfinite(np.asarray(f, dtype=np.float32))
    if clip_sigma is not None and f is not None:
        with np.errstate(all="ignore"):
            valid &= np.abs(np.asarray(f, dtype=np.float32)) / np.sqrt(
                np.where(v > 0, v, np.inf)) <= clip_sigma
    if epoch_ok is not None:
        valid &= np.asarray(epoch_ok, dtype=bool).reshape(
            (-1,) + (1,) * (valid.ndim - 1))
    w = np.where(valid, 1.0 / np.where(valid, v, 1.0), 0.0).astype(np.float32)
    with np.errstate(all="ignore"):
        wmed = np.nanmedian(np.where(w > 0, w, np.nan), axis=0)
    cap = WEIGHT_CAP * np.nan_to_num(wmed, nan=np.inf, posinf=np.inf)
    w = np.minimum(w, cap[None]).astype(np.float32)
    return w, valid


def stack_S(f: np.ndarray, v: np.ndarray, g: np.ndarray, epoch_ok=None,
            return_parts: bool = False, clip_sigma: float | None = None,
            min_epochs: int = MIN_EPOCHS):
    """Stack significance S over epochs (axis 0) for one trajectory.

    ``f, v, g`` have shape (E, *grid). Nodes with fewer than
    ``min_epochs`` contributing epochs are NaN. With ``return_parts``
    also returns (A, B, n_epochs, w) for downstream diagnostics.
    """
    w, valid = stack_weights(v, g, epoch_ok, f=f, clip_sigma=clip_sigma)
    f = np.where(valid, np.nan_to_num(np.asarray(f, dtype=np.float32)), 0.0)
    A = (f * w).sum(axis=0)
    B = w.sum(axis=0)
    n = (w > 0).sum(axis=0)
    with np.errstate(all="ignore"):
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
    S = np.where(n >= min_epochs, S, np.nan)
    if return_parts:
        return S, A, B, n, w
    return S


def grid_max(S: np.ndarray):
    """(S_max, argmax index tuple) with NaN-safety."""
    if not np.isfinite(S).any():
        return np.nan, None
    i = np.unravel_index(np.nanargmax(S), S.shape)
    return float(S[i]), tuple(int(k) for k in i)


def cell_maxima(f: np.ndarray, v: np.ndarray, g: np.ndarray,
                epoch_ok=None, clip_sigma: float | None = None,
                min_epochs: int = MIN_EPOCHS) -> np.ndarray:
    """Grid maximum per trajectory for a (T, E, *grid) tensor slice."""
    return np.array([grid_max(stack_S(f[t], v[t], g[t], epoch_ok,
                                      clip_sigma=clip_sigma,
                                      min_epochs=min_epochs))[0]
                     for t in range(f.shape[0])])


# -- exceedance ratios --------------------------------------------------
def exceedance_ratios(maxima: np.ndarray, designated: np.ndarray):
    """R for every trajectory given the designated-control indices.

    ``maxima`` (T,) grid maxima; ``designated`` indices into it. The
    threshold for trajectory t is the max over designated controls
    other than t (leave-one-out for the designated controls themselves,
    the full designated max for everything else). Returns (R, T_full).
    """
    maxima = np.asarray(maxima, dtype=float)
    designated = np.asarray(designated, dtype=int)
    dmax = maxima[designated]
    T_full = float(np.nanmax(dmax)) if np.isfinite(dmax).any() else np.nan
    R = np.full(len(maxima), np.nan)
    for t in range(len(maxima)):
        if not np.isfinite(maxima[t]):
            continue
        others = dmax[designated != t]
        if not np.isfinite(others).any():
            continue
        thr = float(np.nanmax(others))
        R[t] = maxima[t] / thr if thr > 0 else np.nan
    return R, T_full


# -- time scrambling ----------------------------------------------------
def phase_parts(f0, v0, g0, phase, epoch_ok=None):
    """Per-parallax-phase stack numerators/denominators (A_p, B_p, n_p)
    of the real trajectory, each (nz, nm, nm); the weight cap is the
    whole-cell one so that sum_p A_p / sqrt(sum_p B_p) == stack_S."""
    w, valid = stack_weights(v0, g0, epoch_ok, f=f0)
    fw = np.where(valid, np.nan_to_num(np.asarray(f0, dtype=np.float32)), 0.0) * w
    phase = np.asarray(phase)
    parts = {}
    for p in sorted(set(phase.tolist())):
        sel = phase == p
        parts[int(p)] = (fw[sel].sum(axis=0), w[sel].sum(axis=0),
                         (w[sel] > 0).sum(axis=0))
    return parts


def phase_scrambled_maxima(parts: dict, n_scramble: int,
                           rng: np.random.Generator) -> np.ndarray:
    """Grid maxima of ``n_scramble`` phase-coherence scrambles.

    Each parallax phase's (A, B, n) grids are relabelled by a common
    random cyclic z-shift and mu-node permutation, independently per
    phase, before the phases are recombined. Within a phase every
    epoch keeps its own values, variance, cadence and static-sky
    contamination (the annual parallax returns the track to the same
    static sources at the same node every year, so single-phase
    contamination is coherent in the real statistic and must stay so in
    the null); only the alignment BETWEEN phases — the coherence a real
    track has and chance contamination lacks — is randomised. A cell
    with one populated phase is invariant under this construction.
    """
    keys = sorted(parts)
    nz, nm, _ = parts[keys[0]][0].shape
    out = np.full(n_scramble, np.nan)
    for k in range(n_scramble):
        A = np.zeros((nz, nm * nm)); B = np.zeros_like(A); n = np.zeros_like(A)
        for p in keys:
            Ap, Bp, npp = (x.reshape(nz, nm * nm) for x in parts[p])
            s = int(rng.integers(0, nz))
            m = rng.permutation(nm * nm)
            A += np.roll(Ap, s, axis=0)[:, m]
            B += np.roll(Bp, s, axis=0)[:, m]
            n += np.roll(npp, s, axis=0)[:, m]
        with np.errstate(all="ignore"):
            S = A / np.sqrt(np.where(B > 0, B, np.inf))
        S = np.where(n >= MIN_EPOCHS, S, np.nan)
        out[k] = np.nanmax(S) if np.isfinite(S).any() else np.nan
    return out


def epoch_scrambled_maxima(f0: np.ndarray, v0: np.ndarray, g0: np.ndarray,
                           n_scramble: int, rng: np.random.Generator,
                           epoch_ok=None) -> np.ndarray:
    """Grid maxima of ``n_scramble`` per-epoch scrambles of the real
    trajectory (E, nz, nm, nm): per epoch the z index is cyclically
    shifted by a random amount and the mu nodes are randomly permuted,
    so every epoch keeps its own values, variance, cadence and phase
    while ALL node coherence — including the static-sky coherence of
    the annual parallax return — is destroyed. This is a noise-only
    floor, not exchangeable with the real statistic in star-contaminated
    cells (measured 2026-08-22 on gj-625/rx/W1: scramble maxima ~12 vs
    spatial controls 26-67); reported as a diagnostic, never pooled.
    """
    w, valid = stack_weights(v0, g0, epoch_ok, f=f0)
    fw = np.where(valid, np.nan_to_num(np.asarray(f0, dtype=np.float32)), 0.0) * w
    E, nz, nm, _ = fw.shape
    fw = fw.reshape(E, nz, nm * nm)
    w = w.reshape(E, nz, nm * nm)
    out = np.full(n_scramble, np.nan)
    zidx = np.arange(nz)
    for k in range(n_scramble):
        shift = rng.integers(0, nz, size=E)
        zperm = (zidx[None, :] + shift[:, None]) % nz           # (E, nz)
        mperm = np.argsort(rng.random((E, nm * nm)), axis=1)     # (E, nm^2)
        e_i = np.arange(E)[:, None, None]
        wsel = w[e_i, zperm[:, :, None], mperm[:, None, :]]
        A = fw[e_i, zperm[:, :, None], mperm[:, None, :]].sum(axis=0)
        B = wsel.sum(axis=0)
        n = (wsel > 0).sum(axis=0)
        with np.errstate(all="ignore"):
            S = A / np.sqrt(np.where(B > 0, B, np.inf))
        S = np.where(n >= MIN_EPOCHS, S, np.nan)
        out[k] = np.nanmax(S) if np.isfinite(S).any() else np.nan
    return out


# -- survey-wide calibration --------------------------------------------
@dataclass
class CellNull:
    key: str
    R_real: float
    pooled: np.ndarray                    # all null R values (finite)
    by_construction: dict = field(default_factory=dict)  # name -> array
    unstable: bool = False
    ks: dict = field(default_factory=dict)

    def scale(self, q: float = 0.95) -> float:
        """Per-cell normaliser: the ``q`` quantile of the pooled null."""
        x = self.pooled[np.isfinite(self.pooled)]
        return float(np.quantile(x, q)) if len(x) else np.nan


def wilson_interval(k: int, n: int, z: float = 1.0) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion k/n."""
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (float(max(0.0, centre - half)), float(min(1.0, centre + half)))


def rank_statement(R_real: float, null: np.ndarray) -> dict:
    """The real R's rank among N exchangeable controls, with a Wilson
    68 % / 95 % interval on the exceedance fraction (review §1 rec. 4)."""
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    n = int(len(null))
    if n == 0 or not np.isfinite(R_real):
        return {"n_controls": n, "n_ge": None, "p_cell": None}
    k = int((null >= R_real).sum())
    return {"n_controls": n, "n_ge": k,
            "p_cell": (k + 1) / (n + 1),
            "frac_ge_ci68": wilson_interval(k, n, 1.0),
            "frac_ge_ci95": wilson_interval(k, n, 1.96)}


def ks_exchangeability(by_construction: dict, alpha: float = 0.01) -> dict:
    """Pairwise two-sample KS tests between null constructions; a cell
    is ``unstable`` if any pair differs at ``alpha`` (plan §3.3)."""
    from itertools import combinations

    from scipy.stats import ks_2samp

    out = {"pairs": {}, "unstable": False, "alpha": alpha}
    names = [k for k, a in by_construction.items()
             if a is not None and np.isfinite(a).sum() >= 5]
    for a, b in combinations(names, 2):
        xa = by_construction[a][np.isfinite(by_construction[a])]
        xb = by_construction[b][np.isfinite(by_construction[b])]
        res = ks_2samp(xa, xb)
        out["pairs"][f"{a}|{b}"] = {"D": float(res.statistic),
                                    "p": float(res.pvalue)}
        if res.pvalue < alpha:
            out["unstable"] = True
    return out


def fwer_threshold(cells: list[CellNull], alpha: float = 0.05,
                   n_draws: int = 10000, seed: int = 0,
                   norm_quantile: float | None = 0.95) -> dict:
    """Survey-wide max null by pseudo-experiments.

    Per draw, one pooled null R per cell (cells treated as independent,
    which makes the threshold conservative for positively correlated
    cells), normalised by the cell's own null ``norm_quantile`` when
    given (R~ = R / q_cell, so that heavily contaminated cells do not
    set the family-wide scale), max over cells. Returns the threshold
    on R~ at family-wise ``alpha``, the draw distribution, the per-cell
    normalisers, and P(max null R~ >= R~_cell) for every cell."""
    rng = np.random.default_rng(seed)
    usable = [c for c in cells if len(c.pooled) > 0 and not c.unstable]
    if not usable:
        return {"R_fwer": np.nan, "n_cells": 0}
    maxes = np.zeros(n_draws)
    scale = {}
    for c in usable:
        q = c.scale(norm_quantile) if norm_quantile else 1.0
        scale[c.key] = q
        draws = rng.choice(c.pooled, size=n_draws, replace=True) / q
        np.maximum(maxes, draws, out=maxes)
    R_fwer = float(np.quantile(maxes, 1 - alpha))
    real_norm = {c.key: (c.R_real / scale[c.key]
                         if c.key in scale and np.isfinite(c.R_real) else None)
                 for c in cells}
    global_p = {k: (float((maxes >= r).mean()) if r is not None else None)
                for k, r in real_norm.items()}
    return {"R_fwer": R_fwer, "alpha": alpha, "n_cells": len(usable),
            "n_draws": n_draws, "seed": seed,
            "norm_quantile": norm_quantile, "scale": scale,
            "R_norm_real": real_norm,
            "max_quantiles": {q: float(np.quantile(maxes, q))
                              for q in (0.5, 0.9, 0.95, 0.99)},
            "global_p": global_p}


def bh_qvalues(pvalues: dict) -> dict:
    """Benjamini-Hochberg q-values for a mapping key -> p (informational)."""
    keys = [k for k, p in pvalues.items() if p is not None and np.isfinite(p)]
    p = np.array([pvalues[k] for k in keys])
    m = len(p)
    if m == 0:
        return {}
    order = np.argsort(p)
    q = np.empty(m)
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, p[i] * m / rank)
        q[i] = running
    return {k: float(q[i]) for i, k in enumerate(keys)}
