"""Archive profile for the v2 engine (project plan §10.1 step 6).

A profile binds the survey-agnostic v2 machinery (tensor build with the
48-offset ring null, trajectory randomisation and phase scramble;
null ensemble and FWER; image-level injections; completeness;
adjudication; geometry check; ledger-driven report tables) to one
archive: where its v1 products are, how to turn one observation into a
FluxMap, what PSF to inject, which catalogue the static test uses,
which observer, grids and zero points. Every numerical parameter of
the decision rule lives on the profile and is hashed into the freeze.

Differences from the WISE v2 scripts (surveys/wise/scripts), which
remain the reference implementation for WISE:

* the effective-epoch weight cap is per FRAME (WEIGHT_CAP x the band's
  median frame weight), so that null trajectories can be accumulated
  streaming without holding 99 x E x grid tensors in memory — needed
  for ZTF's 192 x 25 grids and PS1's 360 x 25;
* the PSF is a hook (empirical PRF grid, Moffat from the frame seeing,
  or the product's PSF cube).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np

from sglsurvey.corridors import CORRIDOR_OF, MEMBERS

RING_ANGLES = tuple(22.5 * k for k in range(16))


def ring_offsets(radii) -> tuple:
    return tuple((round(r * np.cos(np.deg2rad(a)), 6), round(r * np.sin(np.deg2rad(a)), 6))
                 for r in radii for a in RING_ANGLES)


def designated_indices(offsets, radii) -> tuple:
    """(±r1,0) (±r2,0) (±r3,0) (0,±r2): the 8 designated controls."""
    r1, r2, r3 = radii
    return tuple(i for i, (x, y) in enumerate(offsets)
                 if (abs(y) < 1e-6 and abs(x) in (r1, r2, r3))
                 or (abs(x) < 1e-6 and abs(y) == r2))


@dataclass
class ArchiveProfile:
    name: str
    survey_dir: Path
    run_dir: Path                  # runs/<survey>-v2
    v1_run_dir: Path               # runs/<survey>
    registry_path: Path
    hypothesis_version: str
    bands: tuple                   # band names in tensor order (idx 1..)
    z_grid: np.ndarray
    mu_grid: np.ndarray
    t0_mjd: float
    zp_ref: float
    mag_system: str                # "vega" | "ab"
    ring_radii: tuple = (20.0, 30.0, 40.0)
    locus_bin_days: float = 0.05
    locus_tol_arcsec: float = 0.5
    donor_bin_days: float = 1.0
    n_donors: int = 50
    n_scramble: int = 200
    n_epoch_scramble: int = 50
    ks_alpha: float = 0.01
    fwer_alpha: float = 0.05
    n_pseudo: int = 10000
    norm_quantile: float = 0.95
    heavy_tail_ratio: float = 2.5
    min_good_frac: float = 0.7
    weight_cap: float = 20.0
    min_epochs: int = 5
    clip_sigma: float | None = None
    quality_masks: dict = field(default_factory=lambda: {"primary": {}, "strict": {}, "loose": {}})
    holdout_split_mjd: float | None = None
    holdout_min_s_late: float = 3.0
    holdout_min_flux_ratio: float = 0.3
    static_search_arcsec: float = 30.0
    static_factor: float = 2.0
    catalog_min_ndet: int = 3
    n_inj_per_cell: int = 400
    temporal_models: tuple = ("persistent", "flicker", "visit", "block")
    inj_mag_halfwidth: float = 2.0
    visit_gap_days: float = 5.0
    n_z_intervals: int = 8
    recovery_window: tuple = (2, 1)
    stamp_half: int = 16
    default_m90: dict = field(default_factory=dict)
    spectrum: dict = field(default_factory=dict)   # band -> spectrum label
    mc_samples: int = 2000
    mc_seed: int = 20260822
    mc_epochs_jyear: tuple = (2019.0, 2022.0, 2025.0)
    mc_z_au: tuple = (550.0, 10000.0)
    xt_threshold_fwhm: float = 0.5
    psf_fwhm_nominal: dict = field(default_factory=dict)   # band -> arcsec (for the xt threshold)
    split_seed: int = 20260822
    dev_fraction: float = 0.30
    #: corridors forced into the development set before the stratified
    #: draw (e.g. pilot corridors whose data shaped the pipeline);
    #: additive — empty tuple reproduces the original draw exactly.
    forced_dev: tuple = ()
    endpoints: list = field(default_factory=list)          # endpoints with v1 tensors
    extra_params: dict = field(default_factory=dict)
    # hooks (set by the concrete profile)
    geometry_context: Callable[[], Any] = None
    load_inputs: Callable[[], tuple] = None         # -> obs_by_id, usable, manifest
    corridor_centre: Callable[[str, dict, dict], tuple] = None   # (corridor, manifest, obs) -> (ra, dec)
    build_map: Callable[[str, dict, dict, bool, Any], Any] = None  # (oid, obs, row, keep_inputs, inject) -> FluxMap|None
    cutout_files: Callable[[str, dict, dict], list] = None   # (oid, obs, row) -> input file paths (for hashing)
    make_psf: Callable[[str, Any], Any] = None        # (band, fm) -> renderer with .render(x, y, half, element)
    psf_element: Callable[[Any, float, float], Any] = None  # (fm, x, y) -> element key or None
    catalog: Callable[[str, tuple], Any] = None      # (corridor, centre) -> StaticCatalog|None
    v1_m90: Callable[[str, str, str], float | None] = None
    confusion_class: Callable[[str], str] = None
    quality_ok: Callable[[dict, str], bool] = None   # (obs quality_flags, mask) -> bool
    epoch_extra: Callable[[Any, dict], dict] = None  # (fm, obs) -> small dict of per-epoch scalars
    fnu_from_mag: Callable[[str, float], float] = None
    inj_window: Callable[[str, str, str], tuple] = None
    sample_adjust: Callable = None
    stamp_response: Callable = None       # (fm, stamp, ox, oy) -> ResponseWindow; default inject.stamp_response
    linear_wcs: bool = True
    static_epoch_scale: Callable = None   # (tensor d, epoch idx, node) -> per-epoch multiplier on predicted static flux               # False: sample with the full (SIP) WCS instead of the local Jacobian        # (fm, obs, offsets_arcsec(N,2), centre, f, v) -> (f, v): e.g. static-template subtraction   # (endpoint, role, band) -> (lo_mag, hi_mag); default m90_v1 ± halfwidth
    observer_identity: str = ""
    inj_subdir: str = "injections"        # switchable so a joint stage can write a second injection set (same draws, its own window)

    # -- derived -----------------------------------------------------------------
    @property
    def band_idx(self) -> dict:
        return {b: i + 1 for i, b in enumerate(self.bands)}

    @property
    def band_name(self) -> dict:
        return {i + 1: b for i, b in enumerate(self.bands)}

    @property
    def ring(self) -> tuple:
        return ring_offsets(self.ring_radii)

    @property
    def designated(self) -> tuple:
        d = designated_indices(self.ring, self.ring_radii)
        assert len(d) == 8, d
        return d

    @property
    def n_traj(self) -> int:
        return 1 + len(self.ring) + self.n_donors

    @property
    def config_dir(self) -> Path:
        return self.survey_dir / "configs"

    @property
    def results_dir(self) -> Path:
        return self.survey_dir / "results"

    @property
    def freeze_path(self) -> Path:
        return self.config_dir / "v2_freeze.json"

    @property
    def tensor_dir(self) -> Path:
        return self.run_dir / "tensors"

    @property
    def null_dir(self) -> Path:
        return self.run_dir / "nulls"

    @property
    def inj_dir(self) -> Path:
        return self.run_dir / self.inj_subdir

    @property
    def hypotheses_path(self) -> Path:
        return self.survey_dir / "hypotheses.md"

    def corridor_of(self, endpoint: str) -> str:
        return CORRIDOR_OF[endpoint]

    def members(self, corridor: str) -> list:
        return [e for e in MEMBERS[corridor] if e in self.endpoints]

    def corridors(self) -> list:
        return sorted({self.corridor_of(e) for e in self.endpoints})

    def load_freeze(self) -> dict:
        return json.loads(self.freeze_path.read_text())

    def freeze_hash(self) -> str:
        return "sha256:" + hashlib.sha256(self.freeze_path.read_bytes()).hexdigest()

    def frozen_params(self) -> dict:
        return {
            "hypothesis_version": self.hypothesis_version, "bands": list(self.bands),
            "z_grid": {"n": int(len(self.z_grid)), "min_au": float(self.z_grid.max()),
                       "max_au": float(self.z_grid.min()), "spacing": "uniform in 1/z"},
            "mu_grid_arcsec_yr": self.mu_grid.tolist(), "motion_bound_norm": "linf",
            "t0_mjd": self.t0_mjd, "zp_ref": self.zp_ref, "mag_system": self.mag_system,
            "locus_bin_days": self.locus_bin_days, "locus_tol_arcsec": self.locus_tol_arcsec,
            "ring_offsets_arcsec": [list(o) for o in self.ring], "designated_controls": list(self.designated),
            "n_phase_scramble": self.n_scramble, "n_epoch_scramble_diagnostic": self.n_epoch_scramble,
            "n_trajectory_donors": self.n_donors, "ks_alpha": self.ks_alpha,
            "fwer_alpha": self.fwer_alpha, "n_pseudo_experiments": self.n_pseudo,
            "norm_quantile": self.norm_quantile, "heavy_tail_ratio": self.heavy_tail_ratio,
            "pooled_null": "ring", "min_good_frac": self.min_good_frac,
            "weight_cap": {"factor": self.weight_cap, "kind": "per-frame (band median frame weight)"},
            "min_epochs": self.min_epochs, "clip_sigma": self.clip_sigma,
            "quality_masks": self.quality_masks,
            "holdout": {"split_mjd": self.holdout_split_mjd, "min_S_late": self.holdout_min_s_late,
                        "min_flux_ratio": self.holdout_min_flux_ratio, "is_veto": False},
            "static_veto": {"search_arcsec": self.static_search_arcsec, "factor": self.static_factor,
                            "min_ndet": self.catalog_min_ndet},
            "injections": {"n_per_cell": self.n_inj_per_cell, "temporal_models": list(self.temporal_models),
                           "mag_halfwidth": self.inj_mag_halfwidth, "visit_gap_days": self.visit_gap_days,
                           "n_z_intervals": self.n_z_intervals, "recovery_window": list(self.recovery_window),
                           "spectrum": self.spectrum, "stamp_half_pix": self.stamp_half},
            "geometry": {"mc_samples": self.mc_samples, "mc_seed": self.mc_seed,
                         "mc_epochs_jyear": list(self.mc_epochs_jyear), "mc_z_au": list(self.mc_z_au),
                         "xt_threshold_fwhm": self.xt_threshold_fwhm, "psf_fwhm_nominal": self.psf_fwhm_nominal,
                         "observer": self.observer_identity},
            "split": {"seed": self.split_seed, "dev_fraction": self.dev_fraction, "unit": "corridor",
                      "strata": "confusion class", "forced_dev": list(self.forced_dev)},
            "extra": self.extra_params,
        }


def default_quality_ok(q: dict, mask: str, masks: dict) -> bool:
    """Generic mask evaluator: each mask is a dict of ``key: [op, value]``
    with op in {">", ">=", "<", "<=", "==", "!=", "is_false"}; keys
    absent from the quality flags pass (the WISE 'where present' rule)."""
    for key, (op, val) in masks[mask].items():
        x = q.get(key)
        if op == "is_false":
            if x:
                return False
            continue
        if x is None:
            continue
        if op == ">" and not x > val: return False
        if op == ">=" and not x >= val: return False
        if op == "<" and not x < val: return False
        if op == "<=" and not x <= val: return False
        if op == "==" and not x == val: return False
        if op == "!=" and not x != val: return False
    return True
