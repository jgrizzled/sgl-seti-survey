"""SPHEREx v2 archive profile (project plan §10.1 step 6, last in the
order; to be re-run when the next quick release adds a parallax phase).

v1 discovery, precise pass, the retained slim cutouts (IMAGE / FLAGS /
VARIANCE / ZODI / PSF plane / WCS-WAVE, runs/spherex/products/cut) and
the v3 static-sky templates (runs/spherex/calib_v4/templates) are
reused; tensors, nulls, injections, vetting, geometry and report are
rebuilt with the engine.

SPHEREx-specific choices:
  * cells per detector D1–D6 (the v3 six-detector joint stack is not
    carried into v2; a joint cell can be formed from the stored
    accumulators later);
  * search image = IMAGE − ZODI with the exposure's PSF plane binned to
    detector sampling at 2 x 2 sub-pixel phases; the static-sky template
    (per corridor x detector, linear in wavelength) subtracted at
    sampling and the variance scaled by the template's reduced chi2,
    exactly as v3;
  * injection PSF = the exposure's own PSF plane (10x oversampled);
    flat-Fν AB spectrum, so the same µJy amplitude in every detector
    (the plan's "inject an SED": flat Fν is the declared SED);
  * grid NZ = 96 uniform in 1/z (3.7"), µ = (−1, 0, +1)"/yr, T0 = 61000,
    ZP 23.9 (µJy); single-epoch clip |S_e| ≤ 5 (the v1 rule); no epoch hold-out (one quick
    release); observer Earth centre (LEO, ≤ 17 mas);
  * no catalogue flux-consistency veto: the static-sky template IS the
    calibrated treatment of catalogued static sources (a catalogue test
    on template-subtracted fluxes would double count; measured on the
    van Maanen test cell 2026-08-22). The 2MASS / CatWISE loaders are
    kept for annotations. Caveat: the template was fitted without the
    injected sources, so partial absorption of a slow (z ~ 10,000 AU)
    real source into the template is not captured by the injections.
"""

from __future__ import annotations

import json
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning

from sglsurvey.adapters.irsa_spherex import SpherexExactFootprint, wavelength_at
from sglsurvey.geometry import GeometryContext
from sglsurvey.inject import OversampledPSF, stamp_response_multiphase
from sglsurvey.photometry import build_flux_map_spherex
from sglsurvey.records import read_records
from sglsurvey.profile import ArchiveProfile, default_quality_ok
from sglsurvey.vetting import StaticCatalog, load_2mass_vizier, load_catwise_vizier

REPO = Path(__file__).resolve().parents[2]
V1 = REPO / "runs" / "spherex"
V1_CAL = V1 / "calib_v4"
TEMPLATE_DIR = V1_CAL / "templates"
WISE_FREEZE = REPO / "surveys" / "wise" / "configs" / "v2_0_freeze.json"
BANDS = ("D1", "D2", "D3", "D4", "D5", "D6")
QUALITY_MASKS = {"primary": {}, "strict": {"deep": ["is_false", None]}, "loose": {}}
CAT_BAND = {"D1": ("2mass", "J", 0.91), "D2": ("2mass", "J", 0.91), "D3": ("2mass", "H", 1.39),
            "D4": ("2mass", "K", 1.85), "D5": ("catwise", "W1", 2.699), "D6": ("catwise", "W2", 3.339)}
_cache: dict = {}


class StaticTemplate:
    """Bilinear lookup of the per-node (a, b) static-sky fit (v3)."""

    def __init__(self, path):
        d = np.load(path)
        self.ra_c, self.dec_c = float(d["ra_c"]), float(d["dec_c"])
        self.sp = float(d["spacing_arcsec"]); self.lam0 = float(d["lam0"])
        self.a, self.b = np.asarray(d["a"], float), np.asarray(d["b"], float)
        rc = np.asarray(d["rchi2"], float) if "rchi2" in d else np.full_like(self.a, np.nan)
        fin = np.isfinite(rc) & (rc > 0)
        self.vscale_default = float(np.clip(np.median(rc[fin]), 0.05, 20.0)) if fin.any() else 1.0
        self.vscale = np.where(fin, np.clip(rc, 0.05, 20.0), self.vscale_default)
        self.n = self.a.shape[0]

    def __call__(self, ra, dec, lam):
        r0, d0 = np.deg2rad(self.ra_c), np.deg2rad(self.dec_c)
        r, d = np.deg2rad(ra), np.deg2rad(dec)
        cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
        xi = np.rad2deg(np.cos(d) * np.sin(r - r0) / cosc) * 3600
        eta = np.rad2deg((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0)) / cosc) * 3600
        x = xi / self.sp + (self.n - 1) / 2; y = eta / self.sp + (self.n - 1) / 2
        out = np.full(ra.shape, np.nan); vs = np.full(ra.shape, self.vscale_default)
        ok = (x >= 0) & (x <= self.n - 1.001) & (y >= 0) & (y <= self.n - 1.001)
        if ok.any():
            x0 = np.floor(x[ok]).astype(int); y0 = np.floor(y[ok]).astype(int)
            fx, fy = x[ok] - x0, y[ok] - y0

            def bil(arr):
                return (arr[y0, x0] * (1 - fx) * (1 - fy) + arr[y0, x0 + 1] * fx * (1 - fy)
                        + arr[y0 + 1, x0] * (1 - fx) * fy + arr[y0 + 1, x0 + 1] * fx * fy)
            out[ok] = bil(self.a) + bil(self.b) * (lam - self.lam0)
            vs[ok] = bil(self.vscale)
        return out, vs


def templates(corridor):
    key = ("tmpl", corridor)
    if key not in _cache:
        _cache[key] = {tp.stem.split("__")[1]: StaticTemplate(tp) for tp in TEMPLATE_DIR.glob(f"{corridor}__D*.npz")}
    return _cache[key]


def endpoints_with_tensors():
    """Endpoints of this survey: those with v2 tensors (the v1 tensor
    products were deleted in the v1 retirement, 2026-08-24)."""
    return sorted({p.stem.split("__")[0] for p in (REPO / "runs" / "spherex" / "v2" / "tensors").glob("*__[rt]x.npz")})


def load_inputs():
    obs_by_id = {r["observation_id"]: r for r in read_records(V1 / "coarse_v1" / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(V1 / "precise_v1" / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(r["observation_id"])
    manifest = {}
    for line in open(V1 / "precise_v1" / "records" / "cutout_index.jsonl"):
        rec = json.loads(line)
        if rec.get("available"):
            manifest[rec["observation_id"]] = rec
    for o in obs_by_id.values():
        o["quality_flags"]["deep"] = bool(o["collection"].endswith("_deep"))
    return obs_by_id, dict(usable), manifest


def corridor_centre(corridor, manifest, obs_by_id):
    t = templates(corridor)
    if t:
        first = next(iter(t.values()))
        return (first.ra_c, first.dec_c)
    raise KeyError(f"no template for {corridor}")


def cutout_files(oid, obs, row):
    return [REPO / row["path"]] if row else []


def build_map(oid, obs, row, keep_inputs, inject):
    if row is None:
        return None
    path = REPO / row["path"]
    try:
        fm = build_flux_map_spherex(path, SpherexExactFootprint.FATAL_MASK, obs["band"], obs["t_mid_mjd_utc"],
                                    keep_inputs=keep_inputs, inject=inject)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            with fits.open(path) as hdul:
                hdr = hdul["IMAGE"].header
                ny, nx = hdul["IMAGE"].data.shape
                wave, bw = wavelength_at(hdul["WCS-WAVE"].data, hdr, nx / 2.0, ny / 2.0)
        fm.wave_um = float(wave); fm.bandwidth_um = float(bw)
        fm.corridor = (row.get("corridor") or [None])[0] if isinstance(row.get("corridor"), list) else row.get("corridor")
        return fm
    except Exception:
        return None


def sample_adjust(fm, obs, off_arcsec, centre, f, v):
    """Subtract the corridor x detector static template at the sampled
    positions and scale the variance by its reduced chi2 (v3 rule)."""
    cor = _cache.get(("corridor_of_centre", centre))
    if cor is None:
        for c in _cache.get("corridors", []):
            t = templates(c)
            if t and abs(next(iter(t.values())).ra_c - centre[0]) < 1e-6:
                cor = c; break
        _cache[("corridor_of_centre", centre)] = cor
    tm = templates(cor).get(obs["band"]) if cor else None
    if tm is None:
        return f, v
    cosd = np.cos(np.deg2rad(centre[1]))
    ra = centre[0] + off_arcsec[:, 0] / 3600.0 / cosd
    dec = centre[1] + off_arcsec[:, 1] / 3600.0
    t, vs = tm(ra, dec, fm.wave_um)
    has = np.isfinite(t)
    return np.where(has, f - np.nan_to_num(t), f), v * vs


def make_psf(band, fm):
    return OversampledPSF(plane=fm.psf_plane, oversample=10, centre=fm.psf_plane.shape[0] // 2, band=band)


def catalog(corridor, centre):
    key = ("cat", corridor)
    if key not in _cache:
        out = {}
        for band, (src, cb, ab_off) in CAT_BAND.items():
            try:
                if src == "2mass":
                    c = load_2mass_vizier(centre[0], centre[1], 0.12, PROFILE.run_dir / "catalogs", band=cb)
                else:
                    c = load_catwise_vizier(centre[0], centre[1], 0.12, PROFILE.run_dir / "catalogs")
                    if cb == "W2":
                        c = StaticCatalog(label=c.label + "-W2", ra=c.ra, dec=c.dec, mag=c.mag, ndet=c.ndet)
                out[band] = StaticCatalog(label=f"{c.label}-AB", ra=c.ra, dec=c.dec,
                                          mag=(c.mag + ab_off) if c.mag is not None else None, ndet=c.ndet)
            except Exception:
                out[band] = None
        _cache[key] = out
    return _cache[key]


def v1_m90(endpoint, role, band):
    d = _cache.setdefault("m90", np.load(V1_CAL / "m90_curves.npz"))
    key = f"{endpoint}__{role}__{band}__0.5"
    if key in d.files and np.isfinite(d[key]).any():
        return float(np.nanmedian(d[key]))
    return None


def confusion_class(corridor):
    cc = _cache.setdefault("cc", json.loads(WISE_FREEZE.read_text())["split"]["confusion_class"])
    return cc[corridor]


def epoch_extra(fm, obs):
    return {"wave_um": float(fm.wave_um), "bandwidth_um": float(fm.bandwidth_um),
            "var_scale": float(fm.var_scale), "deep": float(obs["quality_flags"].get("deep", False))}


PROFILE = ArchiveProfile(
    name="spherex", survey_dir=REPO / "surveys" / "spherex", run_dir=REPO / "runs" / "spherex" / "v2", v1_run_dir=V1,
    registry_path=REPO / "registries" / "pilot_wise_2026.yaml", hypothesis_version="spherex-hypotheses-v2.0",
    bands=BANDS, z_grid=1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 96), mu_grid=np.array([-1.0, 0.0, 1.0]),
    t0_mjd=61000.0, zp_ref=23.9, mag_system="ab", locus_bin_days=0.05, locus_tol_arcsec=0.5, donor_bin_days=0.25,
    clip_sigma=5.0, quality_masks=QUALITY_MASKS, holdout_split_mjd=None, stamp_half=8,
    default_m90={b: 20.0 for b in BANDS}, spectrum={b: "flat_fnu_ab" for b in BANDS},
    mc_epochs_jyear=(2025.5, 2026.0, 2026.5), psf_fwhm_nominal={b: 6.0 for b in BANDS},
    endpoints=endpoints_with_tensors(),
    extra_params={"v1_calib": "calib_v4", "psf": "exposure PSF plane (10x oversampled)", "template": "v3 static-sky template subtracted at sampling",
                  "catalog": "2MASS J/H/K (D1-D4), CatWISE W1/W2 (D5-D6), AB", "observer_note": "Earth centre (<= 17 mas)"},
    geometry_context=GeometryContext.spherex_default, load_inputs=load_inputs, corridor_centre=corridor_centre,
    build_map=build_map, cutout_files=cutout_files, make_psf=make_psf, psf_element=lambda fm, x, y: None,
    catalog=None, v1_m90=v1_m90, confusion_class=confusion_class,
    quality_ok=lambda q, m: default_quality_ok(q, m, QUALITY_MASKS), epoch_extra=epoch_extra,
    fnu_from_mag=lambda band, mag: 3631.0 * 10 ** (-0.4 * mag), observer_identity="earth-center (SPHEREx LEO, <= 17 mas)",
    sample_adjust=sample_adjust, stamp_response=stamp_response_multiphase, linear_wcs=False,
)
_cache["corridors"] = PROFILE.corridors()
