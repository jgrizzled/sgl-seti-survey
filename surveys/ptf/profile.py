"""PTF/iPTF v2 archive profile (project plan §4.13 — Pipeline A over
the IRSA PTF level-1 epochal images, 2009-03 → 2015-01; built directly
on the v2 engine like DECam, no v1 exploratory pass).

PTF-specific choices (hypotheses v1.0 §8, frozen 2026-09-10):
  * scie-direct substrate (the archive serves no difference images):
    static sky is in the search image, so the catalogued-static
    flux-consistency test is the calibrated veto at full weight
    (static_epoch_scale = None, PS1 v2 pattern);
  * per-frame zero point = PS1 DR2 mean-star calibration through the
    identical matched filter (`runs/ptf/zeropoints.jsonl`,
    scripts/calibrate_zeropoints.py; g = PS1 g, R = Jordi 2006 Cousins
    R from PS1 r,i + 0.21 mag to AB) — header MAGZPT never used
    (non-photometric-night values, blank on PHTCALEX = 0 frames);
    frames failing the >= 5-star / <= 0.2-mag gate are unusable;
  * grid NZ = 192 uniform in 1/z (ZTF's, same 1.01"/pix scale), 5 x 5
    mu, T0 = MJD 56000 (2012-03, PTF mid-baseline), ZP 25 AB;
    single-epoch clip 5 sigma (layered search);
  * PSF for injection: Moffat(beta = 3) at the header SEEING;
  * observer Palomar P48 = GeometryContext.ztf_default() (PTF and ZTF
    are the same telescope and site); closed archive -> no epoch
    hold-out split (PS1 pattern);
  * static-test catalogue = PS1 DR2 mean objects of the corridor's
    PS1 screen cone (runs/panstarrs/screen_v1), transformed to the PTF
    band, >= 3 detections;
  * full TAN-SIP WCS sampling (linear_wcs = False) unless the recon
    bounds the linear-Jacobian error under the locus tolerance.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from sglsurvey.adapters.irsa_ptf import MASK_FATAL_TEMPLATE
from sglsurvey.geometry import GeometryContext
from sglsurvey.inject import MoffatPSF
from sglsurvey.photometry import PTF_PIX_ARCSEC, build_flux_map_ptf
from sglsurvey.records import read_records
from sglsurvey.profile import ArchiveProfile, default_quality_ok
from sglsurvey.vetting import StaticCatalog

sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
from ptf_corridors import CORRIDOR_OF as PTF_CORRIDOR_OF  # noqa: E402
from ptf_corridors import SEARCHABLE_ENDPOINTS  # noqa: E402
from ptf_calib import mould_r, ps1_mean_rows  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
RUNS = REPO / "runs" / "ptf"
CUT_DIR = RUNS / "products" / "cut"
MSK_DIR = RUNS / "products" / "msk"
WISE_FREEZE = REPO / "surveys" / "wise" / "configs" / "v2_0_freeze.json"
BANDS = ("g", "R")
#: Cousins R Vega -> AB (Blanton & Roweis 2007), applied to the Jordi
#: transform so that both bands are calibrated in AB.
R_AB_OFFSET = 0.21
PILOT_CORRIDORS: tuple = ("ross128", "kapteyn", "ltt1445", "vanmaanen")   # hypotheses v1.0 §8.1
QUALITY_MASKS = {
    "primary": {"seeing": ["<=", 4.0]},
    "strict": {"seeing": ["<=", 2.5], "moonillf_abs": ["<=", 0.8]},
    "loose": {},
}
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


def load_inputs():
    obs_by_id = {r["observation_id"]: r for r in read_records(
        RUNS / "coarse_v1" / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(RUNS / "precise_v1" / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(r["observation_id"])
    manifest = {}
    mp = CUT_DIR / "manifest.jsonl"
    if mp.exists():
        for line in open(mp):
            row = json.loads(line)
            if "sci" in row.get("files", {}) and row.get("msk"):
                manifest[row["observation_id"]] = row
    for o in obs_by_id.values():
        q = o["quality_flags"]
        q.setdefault("seeing", None)
        # moonillf is signed in the IBE metadata (waxing/waning); masks
        # key on its magnitude
        q["moonillf_abs"] = (abs(q["moonillf"]) if q.get("moonillf") is not None
                             else None)
    return obs_by_id, dict(usable), manifest


def corridor_centre(corridor, manifest, obs_by_id):
    for row in manifest.values():
        cs = row["corridor"] if isinstance(row["corridor"], list) else [row["corridor"]]
        if corridor in cs:
            return (float(row["center_ra_deg"]), float(row["center_dec_deg"]))
    ov = _cache.setdefault("overlay", json.loads(
        (REPO / "surveys" / "ptf" / "targets" / "overlay_v1.json").read_text()))
    for r in ov["rows"]:
        if r["corridor"] == corridor:
            return (float(r["antipode_ra_deg"]), float(r["antipode_dec_deg"]))
    raise KeyError(corridor)


def cutout_files(oid, obs, row):
    if row is None:
        return []
    return [CUT_DIR / row["files"]["sci"], MSK_DIR / row["msk"]]


def build_map(oid, obs, row, keep_inputs, inject):
    if row is None:
        return None
    zp = zeropoints().get(oid)
    if zp is None or zp.get("zp_star") is None or not np.isfinite(zp["zp_star"]):
        return None
    try:
        fm = build_flux_map_ptf(CUT_DIR / row["files"]["sci"], MSK_DIR / row["msk"],
                                MASK_FATAL_TEMPLATE, obs["band"], obs["t_mid_mjd_utc"],
                                keep_inputs=keep_inputs, inject=inject)
    except Exception:
        return None
    fm.magzp_header = fm.magzp
    fm.magzp = float(zp["zp_star"])
    return fm


def make_psf(band, fm):
    return MoffatPSF(fwhm_pix=float(fm.fwhm_arcsec) / PTF_PIX_ARCSEC, beta=3.0, band=band)


def catalog(corridor, centre):
    """PS1 DR2 mean objects in the corridor's PS1 screen cone, one
    StaticCatalog per PTF band with the transformed magnitude (objects
    lacking the transform's inputs fall back to the nearest PS1 band)."""
    key = ("cat", corridor)
    if key not in _cache:
        try:
            rows = ps1_mean_rows(corridor)
        except Exception:
            _cache[key] = None
            return None
        ra = np.array([float(r["raMean"]) for r in rows])
        dec = np.array([float(r["decMean"]) for r in rows])
        nd = np.array([int(float(r["nDetections"])) for r in rows])
        g = np.array([float(r["gMeanPSFMag"]) for r in rows])
        rr = np.array([float(r["rMeanPSFMag"]) for r in rows])
        ii = np.array([float(r["iMeanPSFMag"]) for r in rows])
        g = np.where(g < -100, np.nan, g)
        rr = np.where(rr < -100, np.nan, rr)
        ii = np.where(ii < -100, np.nan, ii)
        R = np.where(np.isfinite(ii), mould_r(rr, np.where(np.isfinite(ii), ii, rr)),
                     rr - 0.117) + R_AB_OFFSET
        _cache[key] = {
            "g": StaticCatalog(label="ps1-dr2-mean (g)", ra=ra, dec=dec, mag=g, ndet=nd),
            "R": StaticCatalog(label="ps1-dr2-mean (Jordi R + 0.21 AB)", ra=ra, dec=dec,
                               mag=R, ndet=nd),
        }
    return _cache[key]


def confusion_class(corridor):
    cc = _cache.setdefault("cc", json.loads(WISE_FREEZE.read_text())["split"]["confusion_class"])
    return cc[corridor]


def epoch_extra(fm, obs):
    q = obs["quality_flags"]
    return {"seeing": float(fm.fwhm_arcsec), "exptime": float(getattr(fm, "exptime", np.nan)),
            "photcalflag": float(q.get("photcalflag") or 0),
            "moonillf": float(q.get("moonillf") if q.get("moonillf") is not None else np.nan),
            "airmass": float(q.get("airmass") if q.get("airmass") is not None else np.nan),
            "magzp_header": float(getattr(fm, "magzp_header", None) or np.nan),
            "bg_sigma": float(getattr(fm, "bg_sigma", np.nan))}


PROFILE = ArchiveProfile(
    name="ptf", survey_dir=REPO / "surveys" / "ptf", run_dir=RUNS / "v2", v1_run_dir=RUNS,
    registry_path=REPO / "registries" / "pilot_wise_2026.yaml",
    hypothesis_version="ptf-hypotheses-v1.0",
    bands=BANDS, z_grid=1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 192),
    mu_grid=np.array([-1.0, -0.5, 0.0, 0.5, 1.0]), t0_mjd=56000.0, zp_ref=25.0, mag_system="ab",
    locus_bin_days=0.05, locus_tol_arcsec=0.5, donor_bin_days=0.25, clip_sigma=5.0,
    quality_masks=QUALITY_MASKS, holdout_split_mjd=None,
    default_m90={"g": 21.5, "R": 21.5}, spectrum={b: "flat_fnu_ab" for b in BANDS},
    mc_epochs_jyear=(2009.5, 2012.0, 2014.5), psf_fwhm_nominal={"g": 2.3, "R": 2.1},
    mc_seed=20260910, split_seed=20260910, forced_dev=PILOT_CORRIDORS,
    endpoints=sorted(SEARCHABLE_ENDPOINTS),
    linear_wcs=False,
    extra_params={
        "search_image": "scie-direct (no difference images in the archive)",
        "psf": "moffat beta=3 at header SEEING",
        "zeropoint": "PS1 DR2 mean-star calibration per frame (zp_star; g = PS1 g, "
                     "R = Jordi 2006 Cousins R + 0.21 AB); header MAGZPT never used",
        "zp_gate": ">= 5 calibrators, robust scatter <= 0.2 mag, else unusable",
        "catalog": "ps1-dr2-mean (runs/panstarrs/screen_v1 corridor cones), per-band transformed",
        "wcs_sampling": "full TAN-SIP inverse (linear_wcs=False)",
        "positive_controls": "(798452) 2012 QR36 catalogue regime (configs/asteroid_control_v1.json); "
                             "(388125) 2005 UP482 stack regime (configs/asteroid_control_388125.json)",
    },
    geometry_context=GeometryContext.ztf_default, load_inputs=load_inputs,
    corridor_centre=corridor_centre, build_map=build_map, cutout_files=cutout_files,
    make_psf=make_psf, psf_element=lambda fm, x, y: None, catalog=catalog, v1_m90=None,
    confusion_class=confusion_class,
    quality_ok=lambda q, m: default_quality_ok(q, m, QUALITY_MASKS), epoch_extra=epoch_extra,
    fnu_from_mag=lambda band, mag: 3631.0 * 10 ** (-0.4 * mag),
    observer_identity="palomar-p48 (terrestrial site; PTF = ZTF telescope)",
)
PROFILE.frozen_at = "2026-09-10"
