"""WISE v3 archive profile — the joint-conventions tensor rebuild
(project plan §4 open items; notes/learnings.md §10: "the WISE-joint
(three-archive) stage on the common µ reference epoch, with Vega→AB and
surface-brightness conventions reconciled").

This is NOT a standalone WISE search. The scientific WISE survey remains
the v2.1 script pipeline (surveys/wise/scripts/, configs/v2_0_freeze.json,
report/wise_survey.md, runs/wise/v2/). This profile binds the flattened
engine (sglsurvey/) to the same v1 WISE inputs to rebuild W1/W2 tensors
and injections under the JOINT stage's conventions, consumed only by
surveys/joint/joint.py (freeze v3.0):

  * common µ reference epoch T0 = MJD 59800 (v2 used 57800; trajectories
    at different reference epochs do not correspond on the grid, so a
    rebuild — not a relabel — is required);
  * AB magnitudes on the common flux scale ZP 25: the frame Vega MAGZP
    is shifted by the WISE Vega→AB offsets (W1 +2.699, W2 +3.339, Cutri
    et al., WISE All-Sky Explanatory Supplement §IV.4.h) before the
    engine's zp_ref scaling, so tensor fluxes are directly comparable
    (and summable) with the PS1/ZTF v2 tensors — a flat-Fnu source has
    equal tensor flux in every band of every archive;
  * µ grid 9 nodes at 0.25"/yr: T0 = 59800 sits ~5 yr from the WISE
    mid-baseline, so the optical 0.5"/yr step would quantise the track
    by up to 0.5 x FWHM at the earliest epochs (the PS1 T0 lesson,
    learnings §4); 0.25"/yr keeps it under 0.26 x FWHM. The PS1 5-node
    grid is the [::2] subgrid, so joint nodes map exactly;
  * W1/W2 only: W3/W4 stay threshold-only with the standalone survey
    (scientific-review rule) and take no part in the joint family;
  * no single-epoch clip and the v2 quality masks — the frame-level
    measurement conventions of the frozen WISE v2 survey are preserved.

Inputs reused verbatim: v1 coarse/precise records and cutouts
(runs/wise/{coarse_v1,precise_v1,products}), the v2 spacecraft-observer
table, PRF grids and covariance envelopes (runs/wise/v2/{observer,prf,
geometry}). Outputs: runs/wise/v3/. Operational freeze (split = the
joint = PS1/ZTF corridor split): configs/v3_freeze.json, written by
``write_freeze()`` — NOT by the engine freeze stage, which would draw a
fresh split.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from sglsurvey.adapters.irsa_wise import WiseExactFootprint
from sglsurvey.geometry import register_wise_spacecraft_observer, wise_v2_context
from sglsurvey.inject import PRFGrid
from sglsurvey.photometry import build_flux_map
from sglsurvey.profile import ArchiveProfile
from sglsurvey.records import read_records

REPO = Path(__file__).resolve().parents[2]
V1 = REPO / "runs" / "wise"
CUT_DIR = V1 / "products" / "cut"
MSK_DIR = V1 / "products" / "msk"
V2_RUN = V1 / "v2"                      # observer / prf / geometry reused
RUN_DIR = V1 / "v3"
WISE_V2_FREEZE = REPO / "surveys" / "wise" / "configs" / "v2_0_freeze.json"
JOINT_V2_FREEZE = REPO / "surveys" / "joint" / "configs" / "v2_freeze.json"

#: m_AB = m_Vega + offset (WISE All-Sky Release Explanatory Supplement)
AB_OFFSET = {"W1": 2.699, "W2": 3.339}

# v2 quality-mask semantics (surveys/wise/scripts/v2common.py, hypotheses
# v2.0 §3.8) — the ops differ from profile.default_quality_ok, so the
# evaluator is ported verbatim rather than re-encoded.
QUALITY_MASKS = {
    "primary": {"qual_frame_gt": 0, "qual_scan_ge_if_present": 5, "saa_sep_gt": 0.0},
    "strict": {"qual_frame_eq": 10, "qual_scan_ge_if_present": 5, "saa_sep_gt": 0.0,
               "moon_sep_ge": 30.0},
    "loose": {"qual_frame_ne": 0},
}

_cache: dict = {}


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


def joint_endpoints():
    """The joint = common PS1/ZTF subset (69 endpoints / 62 corridors),
    pinned by the joint v2 freeze rather than re-derived from globs."""
    s = json.loads(JOINT_V2_FREEZE.read_text())["split"]
    return sorted(set(s["development"]["endpoints"]) | set(s["confirmatory"]["endpoints"]))


def load_inputs():
    obs_by_id = {r["observation_id"]: r for r in read_records(V1 / "coarse_v1" / "records" / "observation.jsonl")}
    usable = {}
    for r in read_records(V1 / "precise_v1" / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable.setdefault((r["endpoint_id"], r["role"]), []).append(r["observation_id"])
    manifest = {}
    for line in open(CUT_DIR / "manifest.jsonl"):
        row = json.loads(line)
        if "int" in row.get("files", {}) and "unc" in row.get("files", {}):
            manifest[row["observation_id"]] = row
    return obs_by_id, usable, manifest


def corridor_centre(corridor, manifest, obs_by_id):
    for row in manifest.values():
        cs = row["corridor"] if isinstance(row["corridor"], list) else [row["corridor"]]
        if corridor in cs:
            return (float(row["center_radec"][0]), float(row["center_radec"][1]))
    raise KeyError(corridor)


def _msk_path(obs):
    return MSK_DIR / obs["products"]["msk"]["url"].rsplit("/", 1)[-1]


def cutout_files(oid, obs, row):
    if row is None:
        return []
    return [CUT_DIR / row["files"]["int"], CUT_DIR / row["files"]["unc"], _msk_path(obs)]


def build_map(oid, obs, row, keep_inputs, inject):
    if row is None:
        return None
    int_p = CUT_DIR / row["files"]["int"]
    unc_p = CUT_DIR / row["files"]["unc"]
    msk_p = _msk_path(obs)
    if not (int_p.exists() and unc_p.exists() and msk_p.exists()):
        return None
    band = obs["band"]
    if band not in AB_OFFSET:
        return None
    try:
        fm = build_flux_map(int_p, unc_p, msk_p, band, obs["t_mid_mjd_utc"],
                            WiseExactFootprint.FATAL_MASK,
                            magzp=obs["quality_flags"].get("magzp"),
                            inject=inject, keep_inputs=keep_inputs,
                            pix_scale_from_header=True)
    except Exception:
        return None
    # Vega -> AB before the engine's zp_ref scaling: the tensor flux
    # scale becomes AB ZP 25, common with the PS1/ZTF v2 tensors.
    fm.magzp = fm.magzp + AB_OFFSET[band] if fm.magzp is not None else None
    return fm


def geometry_context():
    tab = np.load(V2_RUN / "observer" / "wise_sc_ephemeris.npz")
    ident = json.loads((V2_RUN / "observer" / "summary.json").read_text())["table_sha256"]
    return wise_v2_context(register_wise_spacecraft_observer(tab["mjd_utc"], tab["xyz_au"], ident))


def _prf(band):
    key = ("prf", band)
    if key not in _cache:
        _cache[key] = PRFGrid.load(band, V2_RUN / "prf" / band.lower())
    return _cache[key]


def make_psf(band, fm):
    return _prf(band)


def psf_element(fm, x, y):
    ox, oy = getattr(fm, "frame_origin", (0, 0))
    return _prf(fm.band).element_of(x + ox, y + oy)


def v1_m90(endpoint, role, band):
    """v1 W1/W2 m90 (Vega), shifted to AB. Informational only — the
    joint stage always overrides the injection window (optical union)."""
    d = _cache.setdefault("m90", np.load(V1 / "calib_v1" / "m90_curves.npz"))
    key = f"{endpoint}__{role}__{band}__0.5"
    if key in d.files and np.isfinite(d[key]).any():
        return float(np.nanmedian(d[key])) + AB_OFFSET[band]
    return None


def confusion_class(corridor):
    cc = _cache.setdefault("cc", json.loads(WISE_V2_FREEZE.read_text())["split"]["confusion_class"])
    return cc[corridor]


@dataclass
class WiseV3Profile(ArchiveProfile):
    """Overrides: the freeze lives at configs/v3_freeze.json (the v2_0
    freeze file is the frozen standalone survey and stays untouched);
    the injection directory is switchable so the joint stage can write
    one per optical band family (identical draws, family-specific
    magnitude windows)."""
    freeze_name: str = "v3_freeze.json"
    inj_subdir: str = "injections"

    @property
    def freeze_path(self) -> Path:
        return self.config_dir / self.freeze_name

    @property
    def inj_dir(self) -> Path:
        return self.run_dir / self.inj_subdir


PROFILE = WiseV3Profile(
    name="wise", survey_dir=REPO / "surveys" / "wise", run_dir=RUN_DIR, v1_run_dir=V1,
    registry_path=REPO / "registries" / "pilot_wise_2026.yaml",
    hypothesis_version="wise-v3.0-joint-conventions",
    bands=("W1", "W2"), z_grid=1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 64),
    mu_grid=np.linspace(-1.0, 1.0, 9), t0_mjd=59800.0, zp_ref=25.0, mag_system="ab",
    locus_bin_days=0.5, locus_tol_arcsec=2.0, donor_bin_days=1.0, clip_sigma=None,
    quality_masks=QUALITY_MASKS, holdout_split_mjd=59579.0,
    default_m90={"W1": 16.7, "W2": 16.3},
    spectrum={"W1": "flat_fnu_ab", "W2": "flat_fnu_ab"},
    mc_epochs_jyear=(2010.0, 2017.0, 2024.0), psf_fwhm_nominal={"W1": 6.1, "W2": 6.4},
    endpoints=joint_endpoints(),
    extra_params={"purpose": "joint v3 colour axis — not a standalone search",
                  "vega_to_ab": AB_OFFSET, "psf": "IRSA pass2 wise-w?-psf-wpro-09x09 (empirical PRF grid)",
                  "mu_grid_note": "9 nodes at 0.25\"/yr; PS1 5-node grid = [::2] subgrid",
                  "reused_v2_products": ["observer", "prf", "geometry/covariance_envelopes.json"]},
    geometry_context=geometry_context, load_inputs=load_inputs, corridor_centre=corridor_centre,
    build_map=build_map, cutout_files=cutout_files, make_psf=make_psf, psf_element=psf_element,
    catalog=None, v1_m90=v1_m90, confusion_class=confusion_class,
    quality_ok=quality_ok, epoch_extra=None,
    fnu_from_mag=lambda band, mag: 3631.0 * 10 ** (-0.4 * mag),
    observer_identity="wise-l1b-spacecraft (SUN2SC header table, runs/wise/v2/observer)",
)


def write_freeze():
    """Operational v3 freeze: the engine parameters of this profile with
    the JOINT split (= the PS1/ZTF corridor split, identical by
    construction) — the scientific freeze is surveys/joint (v3.0). The
    engine freeze stage is bypassed because it would draw a fresh
    stratified split over the 62 joint corridors."""
    from sglseti import canonical_json

    joint = json.loads(JOINT_V2_FREEZE.read_text())
    fr = {
        "survey": "wise", "hypothesis_version": PROFILE.hypothesis_version,
        "purpose": ("joint v3.0 conventions rebuild (T0 59800, AB ZP 25, W1/W2): tensor and "
                    "injection inputs for surveys/joint/joint.py; no standalone WISE search or "
                    "decision rule is defined by this freeze"),
        "hypotheses_hash": "sha256:" + hashlib.sha256((PROFILE.survey_dir / "hypotheses.md").read_bytes()).hexdigest(),
        "standalone_wise_v2_freeze": "sha256:" + hashlib.sha256(WISE_V2_FREEZE.read_bytes()).hexdigest(),
        "joint_v2_freeze": joint["freeze_content_hash"],
        "parameters": PROFILE.frozen_params(),
        "split": joint["split"],
        "frozen_at": "2026-08-24",
    }
    fr["freeze_content_hash"] = "sha256:" + hashlib.sha256(canonical_json(fr).encode()).hexdigest()
    PROFILE.config_dir.mkdir(parents=True, exist_ok=True)
    PROFILE.freeze_path.write_text(json.dumps(fr, indent=1, sort_keys=True) + "\n")
    print("wise v3 freeze", fr["freeze_content_hash"][:23])
    return fr
