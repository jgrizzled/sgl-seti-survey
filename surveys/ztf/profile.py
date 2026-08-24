"""ZTF v2 archive profile for the v2 engine (project plan §10.1 step 6;
WISE v2_plan.md §12 / §12.1 amendments).

v1 discovery, precise pass, screening and the retained sci/diff/msk
cutouts (runs/ztf/products) are reused; the tensors, null ensemble,
injections, vetting, geometry check and report are rebuilt.

ZTF-specific choices:
  * search image = v1's pixelwise hybrid (difference image where the
    reference exists, sky-subtracted science image elsewhere); the
    injected source is added to that image before the matched filter
    (the science-image PSF in the difference image is an approximation
    to injecting before differencing — the asteroid control measures
    the chain's throughput on real moving sources);
  * PSF for injection: Moffat(beta = 3) at the frame's SEEING (ZTF ships
    no per-frame PRF in the retained products);
  * grid NZ = 192 uniform in 1/z, 5 x 5 mu, T0 = MJD 59800, ZP 25 AB;
    single-epoch clip |S_e| <= 5 as in v1 (layered search);
  * observer Palomar P48; epoch hold-out split at MJD 60554
    (2024-09-01, the last observing year) as an annotation;
  * static test catalogue: ZTF DR objects snapshotted by v1 screening
    (nearest corridor cone, >= 3 good observations); the predicted
    static flux is scaled per epoch by the science-image fraction
    1 - dfrac at the node — static sources are absent from the
    difference-image regime (found 2026-08-23: without this the veto
    "explained" injected sources with subtracted stars).
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from sglsurvey.adapters.irsa_ztf import ZtfExactFootprint
from sglsurvey.geometry import GeometryContext
from sglsurvey.inject import MoffatPSF
from sglsurvey.photometry import ZTF_PIX_ARCSEC, build_flux_map_ztf
from sglsurvey.records import read_records
from sglsurvey.profile import ArchiveProfile, default_quality_ok
from sglsurvey.vetting import load_ztf_objects

REPO = Path(__file__).resolve().parents[2]
V1 = REPO / "runs" / "ztf"
CUT_DIR = V1 / "products" / "cut"
MSK_DIR = V1 / "products" / "msk"
WISE_FREEZE = REPO / "surveys" / "wise" / "configs" / "v2_0_freeze.json"

QUALITY_MASKS = {
    "primary": {"bad_quality": ["is_false", None], "seeing": ["<=", 4.0]},
    "strict": {"bad_quality": ["is_false", None], "seeing": ["<=", 2.5], "maglimit": [">=", 19.5],
               "moonillf": ["<=", 0.8]},
    "loose": {},
}
_cache: dict = {}


def endpoints_with_tensors():
    """Endpoints of this survey: those with v2 tensors (the v1 tensor
    products were deleted in the v1 retirement, 2026-08-24)."""
    return sorted({p.stem.split("__")[0] for p in (REPO / "runs" / "ztf" / "v2" / "tensors").glob("*__[rt]x.npz")})


def load_inputs():
    obs_by_id = {r["observation_id"]: r for r in read_records(V1 / "coarse_v1" / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(V1 / "precise_v1" / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(r["observation_id"])
    manifest = {}
    for line in open(CUT_DIR / "manifest.jsonl"):
        row = json.loads(line)
        if "sci" in row.get("files", {}) and row.get("msk"):
            manifest[row["observation_id"]] = row
    return obs_by_id, dict(usable), manifest


def corridor_centre(corridor, manifest, obs_by_id):
    for row in manifest.values():
        cs = row["corridor"] if isinstance(row["corridor"], list) else [row["corridor"]]
        if corridor in cs:
            return (float(row["center_ra_deg"]), float(row["center_dec_deg"]))
    raise KeyError(corridor)


def cutout_files(oid, obs, row):
    if row is None:
        return []
    fs = [CUT_DIR / row["files"]["sci"], MSK_DIR / row["msk"]]
    if "diff" in row["files"]:
        fs.append(CUT_DIR / row["files"]["diff"])
    return fs


def build_map(oid, obs, row, keep_inputs, inject):
    if row is None:
        return None
    diff = CUT_DIR / row["files"]["diff"] if "diff" in row["files"] else None
    try:
        return build_flux_map_ztf(CUT_DIR / row["files"]["sci"], diff, MSK_DIR / row["msk"],
                                  ZtfExactFootprint.FATAL_MASK, obs["band"], obs["t_mid_mjd_utc"],
                                  keep_inputs=keep_inputs, inject=inject)
    except Exception:
        return None


def make_psf(band, fm):
    return MoffatPSF(fwhm_pix=float(fm.fwhm_arcsec) / ZTF_PIX_ARCSEC, beta=3.0, band=band)


def catalog(corridor, centre):
    key = ("cat", corridor)
    if key not in _cache:
        try:
            _cache[key] = load_ztf_objects(V1 / "screen_v1", centre[0], centre[1])
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
    return {"seeing": float(fm.fwhm_arcsec), "maglimit": float(obs["quality_flags"].get("maglimit") or np.nan),
            "diff_frac_centre": float(np.nanmedian(fm.aux)) if fm.aux is not None else np.nan}


PROFILE = ArchiveProfile(
    name="ztf", survey_dir=REPO / "surveys" / "ztf", run_dir=REPO / "runs" / "ztf" / "v2", v1_run_dir=V1,
    registry_path=REPO / "registries" / "pilot_wise_2026.yaml", hypothesis_version="ztf-hypotheses-v2.0",
    bands=("zg", "zr", "zi"), z_grid=1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 192),
    mu_grid=np.array([-1.0, -0.5, 0.0, 0.5, 1.0]), t0_mjd=59800.0, zp_ref=25.0, mag_system="ab",
    locus_bin_days=0.05, locus_tol_arcsec=0.5, donor_bin_days=0.25, clip_sigma=5.0,
    quality_masks=QUALITY_MASKS, holdout_split_mjd=60554.0,
    default_m90={"zg": 21.5, "zr": 21.5, "zi": 20.5}, spectrum={b: "flat_fnu_ab" for b in ("zg", "zr", "zi")},
    mc_epochs_jyear=(2018.5, 2022.0, 2025.5), psf_fwhm_nominal={"zg": 2.1, "zr": 2.0, "zi": 1.9},
    endpoints=endpoints_with_tensors(),
    extra_params={"search_image": "hybrid diff/sci (v1)", "psf": "moffat beta=3 at SEEING",
                  "catalog": "ztf_objects_dr24 (screen_v1 snapshots)"},
    geometry_context=GeometryContext.ztf_default, load_inputs=load_inputs, corridor_centre=corridor_centre,
    build_map=build_map, cutout_files=cutout_files, make_psf=make_psf, psf_element=lambda fm, x, y: None,
    catalog=catalog, v1_m90=v1_m90, confusion_class=confusion_class,
    quality_ok=lambda q, m: default_quality_ok(q, m, QUALITY_MASKS), epoch_extra=epoch_extra,
    fnu_from_mag=lambda band, mag: 3631.0 * 10 ** (-0.4 * mag), observer_identity="palomar-p48 (terrestrial site)",
)


PARTNER_BAND = {"zg": "g", "zr": "r", "zi": "i"}


def inj_window(endpoint, role, band):
    """Joint-stage magnitude window: the union of this archive's and the
    PS1 partner band's v1 m90 ± 2 mag, so that the j-th injection of a
    cell is the same physical source in both archives."""
    m = v1_m90(endpoint, role, band) or PROFILE.default_m90[band]
    ms = [m]
    pd = _cache.setdefault("m90_ps1", np.load(REPO / "runs" / "panstarrs" / "calib_v1" / "m90_curves.npz"))
    k = f"{endpoint}__{role}__{PARTNER_BAND[band]}__0.5"
    if k in pd.files and np.isfinite(pd[k]).any():
        ms.append(float(np.nanmedian(pd[k])))
    return (min(ms) - PROFILE.inj_mag_halfwidth, max(ms) + PROFILE.inj_mag_halfwidth)


PROFILE.inj_window = inj_window


def static_epoch_scale(d, eb, node):
    """Static sources are absent from the difference-image regime: the
    predicted catalogue flux at an epoch is scaled by the science-image
    weight fraction 1 - dfrac at the node (tensor ``aux``)."""
    if "aux" not in d.files or d["aux"].size == 0:
        return None
    dfrac = d["aux"][eb][:, node[0], node[1], node[2]].astype(float)
    return 1.0 - np.clip(np.nan_to_num(dfrac, nan=0.0), 0.0, 1.0)


PROFILE.static_epoch_scale = static_epoch_scale
