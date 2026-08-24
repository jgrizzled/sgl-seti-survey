"""Pan-STARRS1 v2 archive profile (project plan §10.1 step 6).

v1 discovery, precise pass, screening and the common-T0 (59800) v1
products are reused. Cutouts were purged after the v1 tensor build
(`purge_products.py`), so `scripts/run_batches.sh` re-fetches them per
batch of corridors (fetch_cutouts.py), builds the v2 tensors, runs the
injections, and purges again.

PS1-specific choices:
  * per-warp zero point = the v1 star calibration (DR2 mean-table stars
    through the same matched filter; `calib_t0_59800/zeropoints.jsonl`),
    looked up by observation id — never the header FPA.ZP;
  * skycell duplicates (one exposure on two skycells) collapsed to the
    warp with the larger usable fraction, as in v1;
  * grid NZ = 360 uniform in 1/z (1.0" spacing), 5 x 5 mu, T0 = 59800
    (the joint-stage common epoch), ZP 25 AB; single-epoch clip 5;
  * PSF for injection: Moffat(beta = 3) at the warp's CHIP.SEEING;
  * observer Haleakala; mission over (2009-2014) -> no epoch hold-out
    annotation; static-test catalogue = PS1 DR2 mean objects from the
    v1 screening snapshots (>= 3 detections).
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from sglsurvey.adapters.mast_ps1 import Ps1ExactFootprint
from sglsurvey.geometry import GeometryContext
from sglsurvey.inject import MoffatPSF
from sglsurvey.photometry import PS1_PIX_ARCSEC, build_flux_map_ps1
from sglsurvey.records import read_records
from sglsurvey.profile import ArchiveProfile, default_quality_ok
from sglsurvey.vetting import load_ps1_mean

REPO = Path(__file__).resolve().parents[2]
V1 = REPO / "runs" / "panstarrs"
CUT_DIR = V1 / "products" / "cut"
MSK_DIR = V1 / "products" / "msk"
V1_CAL = V1 / "calib_t0_59800"
WISE_FREEZE = REPO / "surveys" / "wise" / "configs" / "v2_0_freeze.json"
EXPOSURE_DUP_DAYS = 0.001
QUALITY_MASKS = {"primary": {"badflag": ["==", 0]},
                 "strict": {"badflag": ["==", 0], "seeing_arcsec": ["<=", 1.5]},
                 "loose": {}}
_cache: dict = {}


def endpoints_with_tensors():
    """Endpoints of this survey: those with v2 tensors (the v1 tensor
    products were deleted in the v1 retirement, 2026-08-24)."""
    return sorted({p.stem.split("__")[0] for p in (REPO / "runs" / "panstarrs" / "v2" / "tensors").glob("*__[rt]x.npz")})


def zeropoints():
    if "zp" not in _cache:
        zp = {}
        for p in (V1_CAL / "zeropoints.jsonl", V1 / "calib_v1" / "zeropoints.jsonl"):
            if p.exists():
                for line in p.read_text().splitlines():
                    if line.strip():
                        r = json.loads(line)
                        zp.setdefault(r["observation_id"], r)
        _cache["zp"] = zp
    return _cache["zp"]


def load_inputs():
    obs_by_id = {r["observation_id"]: r for r in read_records(V1 / "coarse_v1" / "records" / "observation.jsonl")}
    usable = defaultdict(dict)
    for r in read_records(V1 / "precise_v1" / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])][r["observation_id"]] = r["usable_fraction"] or 0.0
    out = {}
    for pair, frames in usable.items():
        groups = defaultdict(list)
        for oid, frac in frames.items():
            o = obs_by_id[oid]
            groups[(o["band"], round(o["t_start_mjd_utc"] / EXPOSURE_DUP_DAYS))].append((frac, oid))
        out[pair] = [sorted(g, reverse=True)[0][1] for g in groups.values()]
    manifest = {}
    mp = CUT_DIR / "manifest.jsonl"
    if mp.exists():
        for line in open(mp):
            row = json.loads(line)
            if "img" in row.get("files", {}) and row.get("msk"):
                manifest[row["observation_id"]] = row
    for oid, o in obs_by_id.items():
        o["quality_flags"]["seeing_arcsec"] = None
    return obs_by_id, out, manifest


def corridor_centre(corridor, manifest, obs_by_id):
    for row in manifest.values():
        cs = row["corridor"] if isinstance(row["corridor"], list) else [row["corridor"]]
        if corridor in cs:
            return (float(row["center_ra_deg"]), float(row["center_dec_deg"]))
    ov = _cache.setdefault("overlay", json.loads((REPO / "surveys" / "panstarrs" / "targets" / "overlay_v1.json").read_text()))
    for r in ov.get("rows", []):
        if r.get("corridor") == corridor:
            return (float(r["anti_ra"]), float(r["anti_dec"]))
    raise KeyError(corridor)


def cutout_files(oid, obs, row):
    if row is None:
        return []
    fs = [CUT_DIR / row["files"]["img"], MSK_DIR / row["msk"]]
    if "wt" in row["files"]:
        fs.append(CUT_DIR / row["files"]["wt"])
    return fs


def build_map(oid, obs, row, keep_inputs, inject):
    if row is None or not (CUT_DIR / row["files"]["img"]).exists():
        return None
    zp = zeropoints().get(oid)
    if zp is None or zp.get("zp_star") is None or not np.isfinite(zp["zp_star"]):
        return None
    wt = CUT_DIR / row["files"]["wt"] if "wt" in row["files"] else None
    try:
        fm = build_flux_map_ps1(CUT_DIR / row["files"]["img"], wt, MSK_DIR / row["msk"], Ps1ExactFootprint.FATAL_MASK,
                                obs["band"], obs["t_mid_mjd_utc"], magzp=float(zp["zp_star"]),
                                keep_inputs=keep_inputs, inject=inject)
    except Exception:
        return None
    obs["quality_flags"]["seeing_arcsec"] = float(fm.fwhm_arcsec)
    return fm


def make_psf(band, fm):
    return MoffatPSF(fwhm_pix=float(fm.fwhm_arcsec) / PS1_PIX_ARCSEC, beta=3.0, band=band)


def catalog(corridor, centre):
    key = ("cat", corridor)
    if key not in _cache:
        try:
            _cache[key] = load_ps1_mean(V1 / "screen_v1", corridor, band="r")
        except Exception:
            _cache[key] = None
    return _cache[key]


def v1_m90(endpoint, role, band):
    d = _cache.setdefault("m90", np.load(V1 / "calib_v1" / "m90_curves.npz"))
    key = f"{endpoint}__{role}__{band}__0.5"
    if key in d.files and np.isfinite(d[key]).any():
        return float(np.nanmedian(d[key]))
    return None


def confusion_class(corridor):
    cc = _cache.setdefault("cc", json.loads(WISE_FREEZE.read_text())["split"]["confusion_class"])
    return cc[corridor]


def epoch_extra(fm, obs):
    return {"seeing": float(fm.fwhm_arcsec), "exptime": float(getattr(fm, "exptime", np.nan)),
            "var_source": {"wt": 1, "1/wt": 2, "robust": 3}.get(getattr(fm, "var_source", ""), 0)}


PROFILE = ArchiveProfile(
    name="ps1", survey_dir=REPO / "surveys" / "panstarrs", run_dir=REPO / "runs" / "panstarrs" / "v2", v1_run_dir=V1,
    registry_path=REPO / "registries" / "pilot_wise_2026.yaml", hypothesis_version="ps1-hypotheses-v2.0",
    bands=("g", "r", "i", "z", "y"), z_grid=1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 360),
    mu_grid=np.array([-1.0, -0.5, 0.0, 0.5, 1.0]), t0_mjd=59800.0, zp_ref=25.0, mag_system="ab",
    locus_bin_days=0.05, locus_tol_arcsec=0.5, donor_bin_days=0.25, clip_sigma=5.0,
    quality_masks=QUALITY_MASKS, holdout_split_mjd=None,
    default_m90={"g": 21.0, "r": 21.0, "i": 20.7, "z": 19.9, "y": 18.9},
    spectrum={b: "flat_fnu_ab" for b in ("g", "r", "i", "z", "y")},
    mc_epochs_jyear=(2010.0, 2012.0, 2014.0), psf_fwhm_nominal={"g": 1.3, "r": 1.2, "i": 1.1, "z": 1.1, "y": 1.0},
    endpoints=endpoints_with_tensors(),
    extra_params={"psf": "moffat beta=3 at CHIP.SEEING", "zeropoint": "v1 star calibration per warp (zp_star)",
                  "catalog": "ps1-dr2-mean (screen_v1 snapshots)", "duplicates": "skycell duplicates collapsed"},
    geometry_context=GeometryContext.ps1_default, load_inputs=load_inputs, corridor_centre=corridor_centre,
    build_map=build_map, cutout_files=cutout_files, make_psf=make_psf, psf_element=lambda fm, x, y: None,
    catalog=catalog, v1_m90=v1_m90, confusion_class=confusion_class,
    quality_ok=lambda q, m: default_quality_ok(q, m, QUALITY_MASKS), epoch_extra=epoch_extra,
    fnu_from_mag=lambda band, mag: 3631.0 * 10 ** (-0.4 * mag), observer_identity="haleakala-ps1 (terrestrial site)",
)


PARTNER_BAND = {"g": "zg", "r": "zr", "i": "zi"}


def inj_window(endpoint, role, band):
    """Joint-stage magnitude window (union with the ZTF partner band)."""
    m = v1_m90(endpoint, role, band) or PROFILE.default_m90[band]
    ms = [m]
    if band in PARTNER_BAND:
        zd = _cache.setdefault("m90_ztf", np.load(REPO / "runs" / "ztf" / "calib_v1" / "m90_curves.npz"))
        k = f"{endpoint}__{role}__{PARTNER_BAND[band]}__0.5"
        if k in zd.files and np.isfinite(zd[k]).any():
            ms.append(float(np.nanmedian(zd[k])))
    return (min(ms) - PROFILE.inj_mag_halfwidth, max(ms) + PROFILE.inj_mag_halfwidth)


PROFILE.inj_window = inj_window
