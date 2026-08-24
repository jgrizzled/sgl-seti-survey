"""DECam v2 archive profile (project plan §10 step 8 — first survey
built directly on the v2 engine, no v1 exploratory pass).

DECam-specific choices (hypotheses v1.0):
  * per-frame zero point = NSC-star calibration through the identical
    matched filter (`runs/decam/zeropoints.jsonl`,
    scripts/calibrate_zeropoints.py) — never the header MAGZERO
    (measured unreliable, results/flux_scale_check.json);
  * flux maps are per-CCD matched filters pasted onto a TAN canvas
    (photometry.build_flux_map_decam; §8.7 mosaic decision) —
    injections enter each CCD image before the filter;
  * frozen exposure selection (§8.1): EXPTIME >= 30 s, bands grizY,
    obs_type object, applied in load_inputs; per gj-1221's LMC pile-up
    (§8.4) the deepest exposure per night per band is kept everywhere;
  * grid NZ = 360 uniform in 1/z (~1.0" spacing), 5 x 5 mu, T0 = 58500
    (DECam mid-baseline), ZP 25 AB; single-epoch clip 5 sigma;
  * PSF for injection: Moffat(beta = 3) at the frame's fitted FWHM;
  * observer CTIO/W84 (validated <= 1.7 mas,
    results/observer_validation.json); archive still accumulating ->
    epoch hold-out annotation at MJD 60400 (~last observing year);
  * static-test catalogue = NSC DR2 objects from the screen_v1
    snapshots (>= 3 detections);
  * positive control = (60000) Miminko
    (configs/asteroid_control_v1.json).
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from sglsurvey.adapters.noirlab_decam import DQ_FATAL_DEFAULT
from sglsurvey.geometry import GeometryContext
from sglsurvey.inject import MoffatPSF
from sglsurvey.photometry import DECAM_PIX_ARCSEC, build_flux_map_decam
from sglsurvey.records import read_records
from sglsurvey.profile import ArchiveProfile, default_quality_ok
from sglsurvey.vetting import load_nsc_objects

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
from decam_corridors import CORRIDOR_OF as DECAM_CORRIDOR_OF  # noqa: E402
from decam_corridors import PILOT_CORRIDORS  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RUNS = REPO / "runs" / "decam"
CUT_DIR = RUNS / "products" / "cut"
DQ_DIR = RUNS / "products" / "dqmask"
SCREEN = RUNS / "screen_v1"
WISE_FREEZE = REPO / "surveys" / "wise" / "configs" / "v2_0_freeze.json"
MIN_EXPTIME_S = 30.0
BANDS = ("g", "r", "i", "z", "Y")
NIGHT_BIN_DAYS = 1.0
QUALITY_MASKS = {"primary": {},
                 "strict": {"seeing_arcsec": ["<=", 1.2]},
                 "loose": {}}
_cache: dict = {}


def zeropoints():
    if "zp" not in _cache:
        zp = {}
        p = RUNS / "zeropoints.jsonl"
        if p.exists():
            for line in p.read_text().splitlines():
                if line.strip():
                    r = json.loads(line)
                    zp[r["observation_id"]] = r
        _cache["zp"] = zp
    return _cache["zp"]


def manifest_rows():
    if "manifest" not in _cache:
        best = {}
        mp = CUT_DIR / "manifest.jsonl"
        if mp.exists():
            for line in open(mp):
                row = json.loads(line)
                if row.get("schema") == "v1":
                    best[row["observation_id"]] = row
        _cache["manifest"] = {oid: r for oid, r in best.items()
                              if r["files"] and not r["missing"]}
    return _cache["manifest"]


def load_inputs():
    obs_by_id = {r["observation_id"]: r for r in read_records(
        RUNS / "coarse_v1" / "records" / "observation.jsonl")}
    usable = defaultdict(dict)
    for r in read_records(RUNS / "precise_v1" / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])][r["observation_id"]] = \
                r["usable_fraction"] or 0.0
    out = {}
    for pair, frames in usable.items():
        # frozen selection (§8.1) + deepest-per-night-per-band (§8.4)
        nights = defaultdict(list)
        for oid, frac in frames.items():
            o = obs_by_id[oid]
            if (o["band"] not in BANDS
                    or (o["exptime_s"] or 0) < MIN_EXPTIME_S
                    or o["quality_flags"].get("obs_type")
                    not in (None, "object")):
                continue
            key = (o["band"], round(o["t_mid_mjd_utc"] / NIGHT_BIN_DAYS))
            nights[key].append((o["exptime_s"], frac, oid))
        out[pair] = [sorted(g, reverse=True)[0][2]
                     for g in nights.values()]
    for o in obs_by_id.values():
        o["quality_flags"].setdefault("seeing_arcsec", None)
    return obs_by_id, out, manifest_rows()


def corridor_centre(corridor, manifest, obs_by_id):
    for row in manifest.values():
        cs = (row["corridor"] if isinstance(row["corridor"], list)
              else [row["corridor"]])
        if corridor in cs:
            return (float(row["center_ra_deg"]),
                    float(row["center_dec_deg"]))
    ov = _cache.setdefault("overlay", json.loads(
        (REPO / "surveys" / "decam" / "targets"
         / "overlay_v1.json").read_text()))
    for r in ov["rows"]:
        if r["corridor"] == corridor:
            return (float(r["cone"]["ra_deg"]),
                    float(r["cone"]["dec_deg"]))
    raise KeyError(corridor)


def _ccd_files(row):
    return [(CUT_DIR / row["files"][f"image:{c}"],
             CUT_DIR / row["files"][f"wtmap:{c}"]
             if f"wtmap:{c}" in row["files"] else None, c)
            for c in row["ccds"] if f"image:{c}" in row["files"]]


def cutout_files(oid, obs, row):
    if row is None:
        return []
    fs = []
    for img, wt, _c in _ccd_files(row):
        fs.append(img)
        if wt is not None:
            fs.append(wt)
    fs.append(DQ_DIR / row["dqmask"])
    return fs


def build_map(oid, obs, row, keep_inputs, inject):
    if row is None:
        return None
    zp = zeropoints().get(oid)
    if zp is None or zp.get("zp_star") is None \
            or not np.isfinite(zp["zp_star"]):
        return None
    try:
        fm = build_flux_map_decam(
            _ccd_files(row), DQ_DIR / row["dqmask"], DQ_FATAL_DEFAULT,
            obs["band"], obs["t_mid_mjd_utc"],
            (row["center_ra_deg"], row["center_dec_deg"]),
            row["size_pix"], magzp=float(zp["zp_star"]),
            keep_inputs=keep_inputs, inject=inject)
    except Exception:
        return None
    if fm is not None:
        obs["quality_flags"]["seeing_arcsec"] = float(fm.fwhm_arcsec)
    return fm


def make_psf(band, fm):
    return MoffatPSF(fwhm_pix=float(fm.fwhm_arcsec) / DECAM_PIX_ARCSEC,
                     beta=3.0, band=band)


def catalog(corridor, centre):
    key = ("cat", corridor)
    if key not in _cache:
        try:
            _cache[key] = load_nsc_objects(SCREEN, centre[0], centre[1])
        except Exception:
            _cache[key] = None
    return _cache[key]


def confusion_class(corridor):
    cc = _cache.setdefault("cc", json.loads(
        WISE_FREEZE.read_text())["split"]["confusion_class"])
    return cc[corridor]


def epoch_extra(fm, obs):
    return {"seeing": float(fm.fwhm_arcsec),
            "exptime": float(getattr(fm, "exptime", np.nan)),
            "n_ccds": int(getattr(fm, "n_ccds", 1)),
            "var_source": {"wt": 1, "1/wt": 2,
                           "robust": 3}.get(getattr(fm, "var_source",
                                                    ""), 0)}


PROFILE = ArchiveProfile(
    name="decam", survey_dir=REPO / "surveys" / "decam",
    run_dir=RUNS / "v2", v1_run_dir=RUNS,
    registry_path=REPO / "registries" / "pilot_wise_2026.yaml",
    hypothesis_version="decam-hypotheses-v1.0",
    bands=BANDS, z_grid=1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 360),
    mu_grid=np.array([-1.0, -0.5, 0.0, 0.5, 1.0]),
    t0_mjd=58500.0, zp_ref=25.0, mag_system="ab",
    locus_bin_days=0.05, locus_tol_arcsec=0.5, donor_bin_days=0.25,
    clip_sigma=5.0, quality_masks=QUALITY_MASKS,
    holdout_split_mjd=60400.0,
    default_m90={"g": 23.0, "r": 22.8, "i": 22.3, "z": 21.7, "Y": 20.5},
    spectrum={b: "flat_fnu_ab" for b in BANDS},
    mc_epochs_jyear=(2014.0, 2019.0, 2024.0),
    psf_fwhm_nominal={"g": 1.2, "r": 1.1, "i": 1.0, "z": 1.0, "Y": 1.0},
    split_seed=20260824, forced_dev=tuple(PILOT_CORRIDORS),
    endpoints=sorted(DECAM_CORRIDOR_OF),
    extra_params={
        "psf": "moffat beta=3 at fitted frame FWHM",
        "zeropoint": "NSC-star calibration per frame (zp_star); "
                     "header MAGZERO never used",
        "catalog": "nsc-dr2-object (screen_v1 snapshots)",
        "exposure_selection": "EXPTIME>=30 grizY object; deepest per "
                              "night per band (hypotheses v1.0 8.1/8.4)",
        "multi_ccd": "per-CCD matched filter pasted on TAN canvas (8.7)",
        "positive_control": "(60000) Miminko, "
                            "configs/asteroid_control_v1.json",
        "loose_mask_note": "loose == primary (the exposure cut is part "
                           "of the frozen selection, not a mask)"},
    geometry_context=GeometryContext.decam_default,
    load_inputs=load_inputs, corridor_centre=corridor_centre,
    build_map=build_map, cutout_files=cutout_files, make_psf=make_psf,
    psf_element=lambda fm, x, y: None, catalog=catalog,
    v1_m90=None, confusion_class=confusion_class,
    quality_ok=lambda q, m: default_quality_ok(q, m, QUALITY_MASKS),
    epoch_extra=epoch_extra,
    fnu_from_mag=lambda band, mag: 3631.0 * 10 ** (-0.4 * mag),
    observer_identity="ctio-blanco-decam W84 (terrestrial site, "
                      "validated <= 1.7 mas)",
)
PROFILE.frozen_at = "2026-08-24"
