"""Shared paths, frozen parameters and helpers for the WISE v2 scripts.

Every numerical parameter of the decision rule lives here (and is
hashed into configs/v2_0_freeze.json by freeze_v2.py); scripts import
them rather than redefining them.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
V1_DIR = REPO / "runs" / "wise"
COARSE_DIR = V1_DIR / "coarse_v1"
PRECISE_DIR = V1_DIR / "precise_v1"
CUT_DIR = V1_DIR / "products" / "cut"
MSK_DIR = V1_DIR / "products" / "msk"
V1_CAL = V1_DIR / "calib_v1"
RUN_DIR = REPO / "runs" / "wise-v2"
TENSOR_DIR = RUN_DIR / "tensors"
NULL_DIR = RUN_DIR / "nulls"
INJ_DIR = RUN_DIR / "injections"
PRF_DIR = RUN_DIR / "prf"
SURVEY_DIR = REPO / "surveys" / "wise-v2"
CONFIG_DIR = SURVEY_DIR / "configs"
RESULTS_DIR = SURVEY_DIR / "results"
FREEZE_PATH = CONFIG_DIR / "v2_0_freeze.json"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = SURVEY_DIR / "hypotheses.md"
HYPOTHESIS_VERSION = "wise-hypotheses-v2.1"

# -- grid (unchanged from v1) ------------------------------------------------
Z_GRID = 1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 64)
MU_GRID = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
T0_MJD = 57800.0
ZP_REF = 20.0
LOCUS_BIN_DAYS = 0.5
BAND_IDX = {"W1": 1, "W2": 2, "W3": 3, "W4": 4}
BAND_NAME = {v: k for k, v in BAND_IDX.items()}
PSF_FWHM = {"W1": 6.1, "W2": 6.4, "W3": 6.5, "W4": 12.0}

# -- null ensemble (hypotheses v2.0 §3) ------------------------------------------
RING_RADII = (20.0, 30.0, 40.0)
RING_ANGLES = tuple(22.5 * k for k in range(16))
RING_OFFSETS = tuple(
    (round(r * np.cos(np.deg2rad(a)), 6), round(r * np.sin(np.deg2rad(a)), 6))
    for r in RING_RADII for a in RING_ANGLES)      # 48, (dRA*cos, dDec) arcsec
#: Designated controls: (±20,0) (±30,0) (±40,0) (0,±30) -> indices in RING_OFFSETS
DESIGNATED = tuple(i for i, (x, y) in enumerate(RING_OFFSETS)
                   if (abs(y) < 1e-6 and abs(x) in (20.0, 30.0, 40.0))
                   or (abs(x) < 1e-6 and abs(y) == 30.0))
assert len(DESIGNATED) == 8, DESIGNATED
N_SCRAMBLE = 200
N_EPOCH_SCRAMBLE = 50        # diagnostic only
N_DONORS = 50
KS_ALPHA = 0.01
FWER_ALPHA = 0.05
N_PSEUDO = 10000
NORM_QUANTILE = 0.95
HEAVY_TAIL_RATIO = 2.5       # ring max / q95 above this -> null_heavy_tail, no constraint
POOLED_NULL = "ring"         # the exchangeable ensemble; trajectory / phase -> annotations
MIN_GOOD_FRAC = 0.7
WEIGHT_CAP = 20.0
MIN_EPOCHS = 5

# -- quality masks (hypotheses v2.0 §3.8) ------------------------------------------
QUALITY_MASKS = {
    "primary": {"qual_frame_gt": 0, "qual_scan_ge_if_present": 5, "saa_sep_gt": 0.0},
    "strict": {"qual_frame_eq": 10, "qual_scan_ge_if_present": 5, "saa_sep_gt": 0.0,
               "moon_sep_ge": 30.0},
    "loose": {"qual_frame_ne": 0},
}


def quality_ok(q: dict, mask: str) -> bool:
    m = QUALITY_MASKS[mask]
    qf = q.get("qual_frame")
    if "qual_frame_ne" in m and (qf is None or qf == m["qual_frame_ne"]):
        return False
    if "qual_frame_gt" in m and (qf is None or not qf > m["qual_frame_gt"]):
        return False
    if "qual_frame_eq" in m and qf != m["qual_frame_eq"]:
        return False
    qs = q.get("qual_scan")
    if "qual_scan_ge_if_present" in m and qs is not None and qs < m["qual_scan_ge_if_present"]:
        return False
    saa = q.get("saa_sep")
    if "saa_sep_gt" in m and saa is not None and not saa > m["saa_sep_gt"]:
        return False
    moon = q.get("moon_sep")
    if "moon_sep_ge" in m and (moon is None or moon < m["moon_sep_ge"]):
        return False
    return True


# -- vetoes / hold-out (hypotheses v2.0 §3.7) ----------------------------------------
HOLDOUT_SPLIT_MJD = 59579.0          # 2021-12-31
HOLDOUT_MIN_S_LATE = 3.0
HOLDOUT_MIN_FLUX_RATIO = 0.3
HOLDOUT_IS_VETO = False              # v2.1: annotation only (hypotheses §8)
STATIC_SEARCH_ARCSEC = 30.0
STATIC_FACTOR = 2.0
CATWISE_MIN_NDET = 3

# -- injections (hypotheses v2.0 §4) --------------------------------------------
N_INJ_PER_CELL = 400
TEMPORAL_MODELS = ("persistent", "flicker", "visit", "block")
INJ_MAG_HALFWIDTH = 2.0
VISIT_GAP_DAYS = 5.0
N_Z_INTERVALS = 8
RECOVERY_WINDOW = (2, 1)   # z nodes, mu nodes around the injection
SPECTRUM = {"W1": "flat_fnu", "W2": "flat_fnu",
            "W3": "blackbody_300K", "W4": "blackbody_300K"}
STAMP_HALF = 16

# -- geometry (hypotheses v2.0 §1) ----------------------------------------------
MC_SAMPLES = 2000
MC_SEED = 20260822
MC_EPOCHS_JYEAR = (2010.0, 2017.0, 2024.0)
MC_Z_AU = (550.0, 10000.0)
XT_THRESHOLD_FWHM = 0.5

SPLIT_SEED = 20260822
DEV_FRACTION = 0.30


def sha256_file(path: Path) -> str:
    return "sha256:" + hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_freeze() -> dict:
    return json.loads(FREEZE_PATH.read_text())


def freeze_hash() -> str:
    return sha256_file(FREEZE_PATH)


def frozen_params() -> dict:
    """Every parameter above as a JSON-able dict (hashed into the freeze)."""
    return {
        "hypothesis_version": HYPOTHESIS_VERSION,
        "z_grid": {"n": 64, "min_au": 550.0, "max_au": 10000.0, "spacing": "uniform in 1/z"},
        "mu_grid_arcsec_yr": MU_GRID.tolist(), "motion_bound_norm": "linf",
        "t0_mjd": T0_MJD, "zp_ref": ZP_REF, "locus_bin_days": LOCUS_BIN_DAYS,
        "ring_offsets_arcsec": [list(o) for o in RING_OFFSETS],
        "designated_controls": list(DESIGNATED),
        "n_phase_scramble": N_SCRAMBLE, "n_epoch_scramble_diagnostic": N_EPOCH_SCRAMBLE,
        "n_trajectory_donors": N_DONORS, "ks_alpha": KS_ALPHA,
        "fwer_alpha": FWER_ALPHA, "n_pseudo_experiments": N_PSEUDO,
        "norm_quantile": NORM_QUANTILE, "heavy_tail_ratio": HEAVY_TAIL_RATIO,
        "pooled_null": POOLED_NULL, "min_good_frac": MIN_GOOD_FRAC,
        "weight_cap": WEIGHT_CAP, "min_epochs": MIN_EPOCHS,
        "quality_masks": QUALITY_MASKS,
        "holdout": {"split_mjd": HOLDOUT_SPLIT_MJD, "min_S_late": HOLDOUT_MIN_S_LATE,
                    "min_flux_ratio": HOLDOUT_MIN_FLUX_RATIO, "is_veto": HOLDOUT_IS_VETO},
        "static_veto": {"search_arcsec": STATIC_SEARCH_ARCSEC, "factor": STATIC_FACTOR,
                        "catalog": "CatWISE2020 (VizieR II/365)", "min_ndet": CATWISE_MIN_NDET},
        "injections": {"n_per_cell": N_INJ_PER_CELL, "temporal_models": list(TEMPORAL_MODELS),
                       "mag_halfwidth": INJ_MAG_HALFWIDTH, "visit_gap_days": VISIT_GAP_DAYS,
                       "n_z_intervals": N_Z_INTERVALS, "recovery_window": list(RECOVERY_WINDOW),
                       "spectrum": SPECTRUM, "prf": "IRSA pass2 wise-w?-psf-wpro-09x09",
                       "stamp_half_pix": STAMP_HALF},
        "geometry": {"mc_samples": MC_SAMPLES, "mc_seed": MC_SEED,
                     "mc_epochs_jyear": list(MC_EPOCHS_JYEAR), "mc_z_au": list(MC_Z_AU),
                     "xt_threshold_fwhm": XT_THRESHOLD_FWHM,
                     "observer": "wise-l1b-spacecraft (SUN2SC header table)"},
        "split": {"seed": SPLIT_SEED, "dev_fraction": DEV_FRACTION,
                  "unit": "corridor", "strata": "overlay_v2 confusion class"},
    }
