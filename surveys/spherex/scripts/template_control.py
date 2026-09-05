"""SPHEREx static-sky template: regeneration and the template-absorption
control (plan §5.15 item O3 (ii); notes/learnings.md §10).

Background. The v3 static templates cited by report/spherex_survey.md
lived in runs/spherex/calib_v2/templates; calib_v4/templates was a
symlink into that directory, and the v1 retirement (2026-08-24) deleted
calib_v1–v3, so the templates were lost while the tensors (which embed
the template-subtracted fluxes) and the slim cutouts survived. The fit
is deterministic in the cutouts, so `fit` regenerates them with the
static_template.py algorithm unchanged (same grid, weights, two-pass
3-sigma clip), and `verify` rebuilds one corridor's tensors into a
scratch directory and compares them with the stored v2 tensors.

Control. The template was fitted to data that contain any real source,
so part of a slow real source is absorbed into (a, b) at the nodes it
visits; the v2 injections add a source to the *sampled* fluxes after the
template subtraction and therefore do not model that loss. `absorb`
measures it: for a ladder of injected sources it refits the template
with the injected per-epoch fluxes added (full two-pass clipped fit on
the corridor x detector epoch cache), and evaluates the stack-weighted
fraction of the injected flux that the refitted template removes at the
source's own track.

Usage:
  uv run python surveys/spherex/scripts/template_control.py fit [--corridors c1 c2] [--workers N]
  uv run python surveys/spherex/scripts/template_control.py verify --corridor vanmaanen
  uv run python surveys/spherex/scripts/template_control.py absorb --set dev|confirmatory|all [--workers N]
"""

from __future__ import annotations

import argparse
import json
import sys
import time as _time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from sglsurvey.adapters.irsa_spherex import SpherexExactFootprint, wavelength_at
from sglsurvey.photometry import build_flux_map_spherex
from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spherex_corridors import CORRIDOR_OF, MEMBERS  # noqa: E402
from static_template import (CLIP_SIGMA, MIN_EPOCHS, MIN_EPOCHS_SLOPE, MIN_GOOD_FRAC,  # noqa: E402
                             SPACING_ARCSEC, inverse_gnomonic, tangent_grid)

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "spherex" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "spherex" / "precise_v1"
TEMPLATE_DIR = REPO / "runs" / "spherex" / "calib_v4" / "templates"
V2 = REPO / "runs" / "spherex" / "v2"
OUT = REPO / "runs" / "spherex" / "template_control"
BANDS = ("D1", "D2", "D3", "D4", "D5", "D6")


# -- inputs -------------------------------------------------------------------

def load_groups():
    obs_by_id = {r["observation_id"]: r for r in read_records(COARSE_DIR / "records" / "observation.jsonl")}
    usable_obs = set()
    for r in read_records(PRECISE_DIR / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable_obs.add(r["observation_id"])
    groups = defaultdict(list)
    with (PRECISE_DIR / "records" / "cutout_index.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            if rec["available"] and rec["observation_id"] in usable_obs:
                groups[(CORRIDOR_OF[rec["endpoints"][0]], rec["band"])].append(rec)
    return obs_by_id, dict(groups)


class EpochCache:
    """Per (corridor, band): the template grid and every epoch's node
    samples (lam, f, w) — exactly the arrays static_template.py fits.
    ``nodes`` (global grid indices) restricts the sampling to a subset
    (the absorption control needs only the nodes near the ladder
    tracks; a full 3636 x 3636 grid x epochs does not fit in memory);
    ``load=False`` sets up the grid only (call ``load`` later)."""

    def __init__(self, corridor, band, recs, obs_by_id, keep_maps=False, nodes=None, load=True):
        self.corridor, self.band = corridor, band
        self.recs, self.obs_by_id, self.keep_maps = recs, obs_by_id, keep_maps
        self.ra_c = float(np.median([r["center_ra_deg"] for r in recs]))
        self.dec_c = float(np.median([r["center_dec_deg"] for r in recs]))
        cosd = np.cos(np.deg2rad(self.dec_c))
        dra = np.array([(r["center_ra_deg"] - self.ra_c) * cosd * 3600 for r in recs])
        ddec = np.array([(r["center_dec_deg"] - self.dec_c) * 3600 for r in recs])
        half = float(np.max(np.hypot(dra, ddec))) + 0.5 * 100 * 6.15 + 10
        self.xi, XI, ETA = tangent_grid(self.ra_c, self.dec_c, half, SPACING_ARCSEC)
        self.side = XI.shape[0]
        self.nn = self.side * self.side
        self._XI, self._ETA = XI.ravel(), ETA.ravel()
        self.lam, self.f, self.w, self.oid, self.mjd, self.maps, self.failed = [], [], [], [], [], [], []
        if load:
            self.load(nodes)

    def load(self, nodes=None):
        self.nodes = np.arange(self.nn) if nodes is None else np.unique(np.asarray(nodes, int))
        self.full = nodes is None
        self.local = np.full(self.nn, -1, np.int64); self.local[self.nodes] = np.arange(len(self.nodes))
        self.RA, self.DEC = inverse_gnomonic(self._XI[self.nodes], self._ETA[self.nodes], self.ra_c, self.dec_c)
        for r in self.recs:
            obs = self.obs_by_id[r["observation_id"]]
            path = REPO / r["path"]
            try:
                fm = build_flux_map_spherex(path, SpherexExactFootprint.FATAL_MASK, obs["band"], obs["t_mid_mjd_utc"])
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    with fits.open(path) as hdul:
                        hdr = hdul["IMAGE"].header
                        ny, nx = hdul["IMAGE"].data.shape
                        lam, _ = wavelength_at(hdul["WCS-WAVE"].data, hdr, nx / 2, ny / 2)
            except Exception as exc:
                self.failed.append((r["observation_id"], str(exc)))
                continue
            f, v, g = fm.sample(self.RA, self.DEC)
            ok = np.isfinite(f) & np.isfinite(v) & (v > 0) & (g >= MIN_GOOD_FRAC)
            w = np.where(ok, 1.0 / np.where(v > 0, v, 1.0), 0.0)
            self.lam.append(float(lam)); self.f.append(np.where(ok, f, 0.0).astype(np.float32))
            self.w.append(w.astype(np.float32)); self.oid.append(r["observation_id"])
            self.mjd.append(float(obs["t_mid_mjd_utc"]))
            if self.keep_maps:
                self.maps.append(fm)
        self.lam = np.array(self.lam); self.mjd = np.array(self.mjd)
        self.lam0 = float(np.median(self.lam)) if len(self.lam) else np.nan

    # the static_template.py fit, on optionally perturbed fluxes / node subsets
    def fit(self, delta=None, nodes=None):
        """Two-pass clipped WLS of f = a + b (lam - lam0) per node.
        ``delta``: optional list of per-epoch additive fluxes (E arrays
        over ``nodes``); ``nodes``: node index subset (default all)."""
        idx = np.arange(len(self.nodes)) if nodes is None else self.local[np.asarray(nodes)]
        assert (idx >= 0).all(), "requested nodes were not loaded"
        fs = [f[idx] for f in self.f]; ws = [w[idx] for w in self.w]
        if delta is not None:
            fs = [f + d for f, d in zip(fs, delta)]
        n = len(idx)

        def _fit(mask=None):
            S0 = np.zeros(n); S1 = np.zeros(n); S2 = np.zeros(n)
            F0 = np.zeros(n); F1 = np.zeros(n); N = np.zeros(n, dtype=np.int32)
            for k in range(len(fs)):
                wk = ws[k] if mask is None else ws[k] * mask[k]
                dl = self.lam[k] - self.lam0
                S0 += wk; S1 += wk * dl; S2 += wk * dl * dl
                F0 += wk * fs[k]; F1 += wk * fs[k] * dl; N += (wk > 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                det = S0 * S2 - S1 * S1
                b = np.where(det > 0, (S0 * F1 - S1 * F0) / det, 0.0)
                a_slope = np.where(det > 0, (S2 * F0 - S1 * F1) / det, np.nan)
                a_const = np.where(S0 > 0, F0 / S0, np.nan)
            a = np.where(N >= MIN_EPOCHS_SLOPE, a_slope, a_const)
            b = np.where(N >= MIN_EPOCHS_SLOPE, b, 0.0)
            a[N < MIN_EPOCHS] = np.nan
            b[N < MIN_EPOCHS] = 0.0
            return a, b, N

        a, b, N = _fit()
        mask = []
        for k in range(len(fs)):
            with np.errstate(invalid="ignore"):
                resid = (fs[k] - (a + b * (self.lam[k] - self.lam0))) * np.sqrt(ws[k])
            mask.append((np.abs(np.nan_to_num(resid)) <= CLIP_SIGMA).astype(np.float32))
        a, b, N = _fit(mask)
        chi2 = np.zeros(n); Nc = np.zeros(n)
        for k in range(len(fs)):
            with np.errstate(invalid="ignore"):
                r_ = (fs[k] - (a + b * (self.lam[k] - self.lam0))) ** 2 * ws[k]
            ok = np.isfinite(r_) & (ws[k] > 0)
            chi2 += np.where(ok, r_, 0.0); Nc += ok
        with np.errstate(invalid="ignore", divide="ignore"):
            rchi2 = np.where(Nc > 2, chi2 / np.maximum(Nc - 2, 1), np.nan)
        return a, b, N, rchi2

    def save_template(self, out_dir):
        assert self.full, "save_template needs the full grid"
        a, b, N, rchi2 = self.fit()
        side = self.side
        out_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            out_dir / f"{self.corridor}__{self.band}.npz",
            ra_c=self.ra_c, dec_c=self.dec_c, spacing_arcsec=SPACING_ARCSEC, xi=self.xi,
            lam0=self.lam0, a=a.reshape(side, side).astype(np.float32),
            b=b.reshape(side, side).astype(np.float32),
            n=N.reshape(side, side), rchi2=rchi2.reshape(side, side).astype(np.float32),
            n_epochs=len(self.f))
        fin = np.isfinite(a)
        return (f"[{self.corridor}/{self.band}] {len(self.f)} epochs, grid {side}x{side}, template defined on "
                f"{fin.mean():.2f} of nodes (median N={np.median(N[fin]) if fin.any() else 0:.0f}), "
                f"median reduced chi2 {np.nanmedian(rchi2):.2f}")

    def node_xy(self, xi_arcsec, eta_arcsec):
        """Tangent-plane offsets (arcsec) -> fractional node coordinates."""
        x = xi_arcsec / SPACING_ARCSEC + (self.side - 1) / 2
        y = eta_arcsec / SPACING_ARCSEC + (self.side - 1) / 2
        return x, y


# -- fit ----------------------------------------------------------------------

_G = {}


def _init(obs_by_id, groups):
    _G["obs_by_id"] = obs_by_id; _G["groups"] = groups


def _fit_group(key):
    corridor, band = key
    t0 = _time.monotonic()
    try:
        ec = EpochCache(corridor, band, _G["groups"][key], _G["obs_by_id"])
        if not len(ec.f):
            return f"[{corridor}/{band}] no epochs"
        msg = ec.save_template(TEMPLATE_DIR)
        return msg + f" ({_time.monotonic() - t0:.0f}s; {len(ec.failed)} fluxmap failures)"
    except Exception as exc:
        import traceback
        return f"[{corridor}/{band}] FAILED {exc}\n{traceback.format_exc()}"


def fit(corridors=None, workers=8, only_missing=True):
    if TEMPLATE_DIR.is_symlink() and not TEMPLATE_DIR.exists():
        TEMPLATE_DIR.unlink()          # the dangling calib_v2 link (see module docstring)
    TEMPLATE_DIR.mkdir(parents=True, exist_ok=True)
    obs_by_id, groups = load_groups()
    keys = sorted(k for k in groups if (not corridors or k[0] in corridors))
    if only_missing:
        keys = [k for k in keys if not (TEMPLATE_DIR / f"{k[0]}__{k[1]}.npz").exists()]
    # biggest groups first so the pool tail is short
    keys.sort(key=lambda k: -len(groups[k]))
    print(f"[template fit] {len(keys)} corridor x detector groups, {sum(len(groups[k]) for k in keys)} cutouts", flush=True)
    log = TEMPLATE_DIR.parent / "template_fit.log"
    if workers <= 1:
        _init(obs_by_id, groups)
        for k in keys:
            msg = _fit_group(k); print(msg, flush=True); open(log, "a").write(msg + "\n")
        return
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers, initializer=_init, initargs=(obs_by_id, groups)) as pool:
        for msg in pool.imap_unordered(_fit_group, keys):
            print(msg, flush=True); open(log, "a").write(msg + "\n")


# -- verify -------------------------------------------------------------------

def verify(corridor, scratch):
    """Rebuild the corridor's v2 tensors with the regenerated templates
    into ``scratch`` and compare with runs/spherex/v2/tensors."""
    import dataclasses
    sys.path.insert(0, str(REPO / "surveys" / "spherex"))
    import profile as spx  # noqa: E402
    from sglsurvey import build
    P = dataclasses.replace(spx.PROFILE, run_dir=Path(scratch))
    spx.PROFILE = P
    build.run(P, corridors=[corridor], workers=1, only_missing=False)
    rep = {}
    for tp in sorted(Path(scratch, "tensors").glob("*.npz")):
        ref = V2 / "tensors" / tp.name
        if not ref.exists():
            rep[tp.name] = "no stored tensor"; continue
        with np.load(tp) as a, np.load(ref) as b:
            same_hash = str(a["input_hash"]) == str(b["input_hash"])
            r = {"input_hash_equal": same_hash, "n_epochs": (int(a["mjd"].size), int(b["mjd"].size))}
            if a["f"].shape == b["f"].shape:
                fa, fb = a["f"][0], b["f"][0]
                fin = np.isfinite(fa) & np.isfinite(fb)
                d = np.abs(fa[fin] - fb[fin]); scale = np.abs(fb[fin]) + 1e-6
                r["f_real_traj"] = {"n_finite": int(fin.sum()), "n_nan_mismatch": int((np.isfinite(fa) != np.isfinite(fb)).sum()),
                                    "max_abs_diff": float(d.max()) if d.size else None,
                                    "max_rel_diff": float((d / scale).max()) if d.size else None,
                                    "frac_exact": float((d == 0).mean()) if d.size else None}
                sa, sb = a["summary"][0, 0, :, 0], b["summary"][0, 0, :, 0]
                r["S_max_real_per_band"] = {"rebuilt": [float(x) for x in sa], "stored": [float(x) for x in sb]}
                va, vb = a["v"][0], b["v"][0]
                fin = np.isfinite(va) & np.isfinite(vb)
                r["v_real_traj_max_rel_diff"] = float((np.abs(va[fin] - vb[fin]) / (np.abs(vb[fin]) + 1e-12)).max()) if fin.any() else None
            else:
                r["shape"] = (list(a["f"].shape), list(b["f"].shape))
            rep[tp.name] = r
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"verify_{corridor}.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


# -- validate (end-to-end) ----------------------------------------------------------

def validate(corridor, endpoint, iz=60, dmag=-1.0, scratch=None, band_check=None):
    """Ground truth for `absorb`: add ONE persistent flat-Fnu source to
    the corridor's IMAGES on the endpoint's rx track (z node ``iz``,
    mu = 0, m90_v1(D) + dmag per detector), regenerate the six templates
    from the injected images, rebuild the pair's tensor twice — with the
    injected templates (what a real source suffers) and with the
    original templates (what the v2 injections emulate) — and compare the
    recovered stack signal at the injected node with the module's f_abs."""
    import dataclasses
    from sglsurvey import build
    from sglsurvey.inject import OversampledPSF
    from sglsurvey.photometry import build_flux_map_spherex as _bfm
    sys.path.insert(0, str(REPO / "surveys" / "spherex"))
    import profile as spx  # noqa: E402
    scratch = Path(scratch or (OUT / "validate" / f"{corridor}__{endpoint}"))
    scratch.mkdir(parents=True, exist_ok=True)
    P = spx.PROFILE
    with np.load(P.tensor_dir / f"{endpoint}__rx.npz") as d:
        centre = tuple(float(x) for x in d["centre"]); oids = [str(o) for o in d["oid"]]
        d0 = d["d0"][:, iz, :].astype(float); band_idx = d["band_idx"].astype(int)
        nm = len(P.mu_grid); m0 = nm // 2
        f_orig = d["f"][0, :, iz, m0, m0].astype(float); v_orig = d["v"][0, :, iz, m0, m0].astype(float); g_orig = d["g"][0, :, iz, m0, m0].astype(float)
        cap = d["frame_cap"].astype(float); z_au = float(d["z_grid"][iz])
    cosd = np.cos(np.deg2rad(centre[1]))
    pos = {o: (centre[0] + d0[i, 0] / 3600.0 / cosd, centre[1] + d0[i, 1] / 3600.0) for i, o in enumerate(oids)}
    mags = {b: (spx.v1_m90(endpoint, "rx", b) or 20.0) + dmag for b in BANDS}
    obs_by_id, groups = load_groups()

    def injector(oid, band):
        if oid not in pos:
            return None
        ra, dec = pos[oid]; F_ujy = 10.0 ** (0.4 * (P.zp_ref - mags[band]))

        def inject(resid, hdr, wcs, psf, omega_sr):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                xy = wcs.all_world2pix(np.array([[ra, dec]]), 0, quiet=True)[0]
            x, y = float(xy[0]), float(xy[1])
            ny, nx = resid.shape
            if not (np.isfinite(x) and np.isfinite(y) and -8 <= x < nx + 8 and -8 <= y < ny + 8):
                return resid
            stamp, ox, oy = OversampledPSF(plane=psf, oversample=10, centre=psf.shape[0] // 2, band=band).render(x, y, half=8)
            out = resid.copy(); n = stamp.shape[0]
            y0, y1 = max(0, oy), min(ny, oy + n); x0, x1 = max(0, ox), min(nx, ox + n)
            if y1 > y0 and x1 > x0:
                out[y0:y1, x0:x1] += (F_ujy / (omega_sr * 1e12)) * stamp[y0 - oy:y1 - oy, x0 - ox:x1 - ox]
            return out
        return inject

    # 1. templates from injected images
    tdir = scratch / "templates"; tdir.mkdir(exist_ok=True)
    orig_bfm = build_flux_map_spherex
    me = sys.modules[__name__]          # EpochCache resolves the builder from this module's globals
    for band in BANDS:
        recs = groups[(corridor, band)]

        def patched(path, fatal, b, mjd, *a, **kw):
            oid = next((r["observation_id"] for r in recs if str(REPO / r["path"]) == str(path)), None)
            return _bfm(path, fatal, b, mjd, *a, inject=injector(oid, b), **kw)
        me.build_flux_map_spherex = patched
        try:
            ec = EpochCache(corridor, band, recs, obs_by_id)
            print(ec.save_template(tdir), flush=True)
        finally:
            me.build_flux_map_spherex = orig_bfm
    # 2. tensors: injected images x {injected templates, original templates}
    orig_build_map = spx.build_map

    def build_map_inj(oid, obs, row, keep_inputs, inject):
        if row is None:
            return None
        path = REPO / row["path"]
        try:
            fm = _bfm(path, spx.SpherexExactFootprint.FATAL_MASK, obs["band"], obs["t_mid_mjd_utc"],
                      keep_inputs=keep_inputs, inject=injector(oid, obs["band"]))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", spx.AstropyWarning)
                with fits.open(path) as hdul:
                    hdr = hdul["IMAGE"].header; ny, nx = hdul["IMAGE"].data.shape
                    wave, bw = spx.wavelength_at(hdul["WCS-WAVE"].data, hdr, nx / 2.0, ny / 2.0)
            fm.wave_um = float(wave); fm.bandwidth_um = float(bw)
            return fm
        except Exception:
            return None
    results = {"corridor": corridor, "endpoint": endpoint, "iz": iz, "z_au": z_au, "dmag": dmag, "mags": mags, "bands": {}}
    for label, tpl in (("injected_template", tdir), ("original_template", TEMPLATE_DIR)):
        rd = scratch / label
        spx.TEMPLATE_DIR = tpl; spx._cache.clear()
        Pi = dataclasses.replace(P, run_dir=rd, build_map=build_map_inj, endpoints=[endpoint])
        spx._cache["corridors"] = [corridor]
        spx.PROFILE = Pi
        build.run(Pi, corridors=[corridor], workers=1, only_missing=False)
        with np.load(rd / "tensors" / f"{endpoint}__rx.npz") as d:
            o2 = [str(o) for o in d["oid"]]; idx = [o2.index(o) for o in oids]
            results[label] = {"f": d["f"][0, idx, iz, m0, m0].astype(float), "v": d["v"][0, idx, iz, m0, m0].astype(float),
                              "g": d["g"][0, idx, iz, m0, m0].astype(float)}
    spx.TEMPLATE_DIR = TEMPLATE_DIR; spx._cache.clear(); spx.PROFILE = P; spx._cache["corridors"] = P.corridors()
    # 3. compare the recovered stack signal at the injected node (mu = 0 -> centre mu node)
    for bi, band in enumerate(BANDS):
        eb = band_idx == bi + 1
        if eb.sum() < P.min_epochs:
            continue
        out = {"n_epochs": int(eb.sum()), "mag": mags[band]}
        # stack the injected node with the ORIGINAL weights (v, g of the unperturbed maps) so the
        # comparison isolates the template's effect on the summed flux
        ok = np.isfinite(f_orig[eb]) & np.isfinite(v_orig[eb]) & (v_orig[eb] > 0) & (g_orig[eb] >= P.min_good_frac)
        w = np.where(ok, np.minimum(1.0 / np.where(ok, v_orig[eb], 1.0), cap[eb]), 0.0)
        for label in ("injected_template", "original_template"):
            fi = results[label]["f"][eb]
            okl = ok & np.isfinite(fi)
            A = float(np.nansum(np.where(okl, fi, 0) * w)); B = float(w[okl].sum())
            A0 = float(np.nansum(np.where(okl, f_orig[eb], 0) * w))
            out[label] = {"S": A / np.sqrt(B) if B > 0 else None, "dS_vs_no_source": (A - A0) / np.sqrt(B) if B > 0 else None,
                          "flux_hat_ujy": (A - A0) / B if B > 0 else None, "n_used": int(okl.sum())}
        dS_i, dS_o = out["injected_template"]["dS_vs_no_source"], out["original_template"]["dS_vs_no_source"]
        out["f_abs_end_to_end"] = 1.0 - dS_i / dS_o if (dS_o and dS_o > 0) else None
        out["injected_flux_ujy"] = 10.0 ** (0.4 * (P.zp_ref - mags[band]))
        out["throughput_original_template"] = (out["original_template"]["flux_hat_ujy"] / out["injected_flux_ujy"]) if out["original_template"]["flux_hat_ujy"] else None
        results["bands"][band] = out
        print(f"[validate {endpoint}/rx/{band}] z {z_au:.0f} mag {mags[band]:.2f}: dS original-template {dS_o:.2f}, "
              f"injected-template {dS_i:.2f} -> end-to-end f_abs {out['f_abs_end_to_end']:.3f}; throughput {out['throughput_original_template']:.3f}", flush=True)
    (OUT / "validate").mkdir(parents=True, exist_ok=True)
    (OUT / "validate" / f"{corridor}__{endpoint}__iz{iz}.json").write_text(json.dumps(results, indent=1, default=lambda x: x.tolist() if hasattr(x, "tolist") else float(x)))
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["fit", "verify", "absorb", "validate"])
    ap.add_argument("--endpoint")
    ap.add_argument("--iz", type=int, default=60)
    ap.add_argument("--dmag", type=float, default=-1.0)
    ap.add_argument("--corridors", nargs="*")
    ap.add_argument("--corridor")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--no-only-missing", action="store_true")
    ap.add_argument("--scratch", default=str(OUT / "verify_scratch"))
    ap.add_argument("--set", choices=["dev", "confirmatory", "all"], default="dev")
    a = ap.parse_args()
    if a.mode == "fit":
        fit(a.corridors, a.workers, not a.no_only_missing)
    elif a.mode == "verify":
        verify(a.corridor, a.scratch)
    elif a.mode == "absorb":
        from template_absorb import absorb   # design fixed at the freeze; see that module
        absorb(a.set, a.workers)
    elif a.mode == "validate":
        validate(a.corridor, a.endpoint, iz=a.iz, dmag=a.dmag)


if __name__ == "__main__":
    main()
