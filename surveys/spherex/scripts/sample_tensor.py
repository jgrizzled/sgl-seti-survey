"""Sampling pass for the SPHEREx pilot: per-epoch trajectory sample
tensors for every endpoint x role x detector, real trajectory plus 16
offset controls (same design as the WISE/ZTF sample_tensor scripts).

Writes runs/spherex/calib_v1/tensors/<endpoint>__<role>__<band>.npz:
  mjd (E,), phase (E,) [0/1], wave_um (E,), bandwidth_um (E,),
  psf_fwhm (E,), var_scale (E,), deep (E,) [0/1],
  f/v (9, E, NZ, NM, NM) float32 in uJy / uJy^2, g float16,
  tsub (9, E, NZ, NM, NM) float16 = 1 where the static-sky template
  (static_template.py) was subtracted from f, 0 where undefined.
  v is the per-cutout (MAD) variance multiplied by the template's
  per-node reduced chi2 (v2 variance calibration, 2026-08-21).

SPHEREx-specific choices: NZ = 96 nodes uniform in 1/z (3.7" spacing,
~0.7 PSF FWHM); MU grid 3 nodes (-1, 0, +1 "/yr: residual motion is
unresolved over the 14-month baseline); T0 = MJD 61000 (baseline
midpoint); fluxes from IMAGE - ZODI with the exposure's own PSF as the
matched filter; per-epoch wavelength at the cutout centre from the
WCS-WAVE table.

Usage: uv run python surveys/spherex/scripts/sample_tensor.py [--only-missing]
"""

from __future__ import annotations

import json
import os
import sys
import time as _time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.irsa_spherex import (SpherexExactFootprint,
                                             wavelength_at)
from sglsurvey.geometry import GeometryContext
from sglsurvey.photometry import build_flux_map_spherex
from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spherex_corridors import CORRIDOR_OF  # noqa: E402
from static_template import inverse_gnomonic  # noqa: E402,F401


class StaticTemplate:
    """Bilinear lookup of the per-node (a, b) static-sky fit."""

    def __init__(self, path):
        d = np.load(path)
        self.ra_c, self.dec_c = float(d["ra_c"]), float(d["dec_c"])
        self.sp = float(d["spacing_arcsec"])
        self.lam0 = float(d["lam0"])
        self.a, self.b = np.asarray(d["a"], float), np.asarray(d["b"], float)
        # Variance calibration (v2, 2026-08-21): the per-cutout MAD noise
        # includes static confusion that the template removes, so the
        # true epoch-to-epoch variance at a node is MAD-variance x the
        # template's reduced chi2. Clip to [0.05, 20]; undefined nodes
        # fall back to the per-band median.
        rc = np.asarray(d["rchi2"], float) if "rchi2" in d else np.full_like(self.a, np.nan)
        fin = np.isfinite(rc) & (rc > 0)
        self.vscale_default = float(np.clip(np.median(rc[fin]), 0.05, 20.0)) if fin.any() else 1.0
        self.vscale = np.where(fin, np.clip(rc, 0.05, 20.0), self.vscale_default)
        self.n = self.a.shape[0]
        self.cosd = np.cos(np.deg2rad(self.dec_c))

    def __call__(self, ra, dec, lam, with_vscale=False):
        # forward gnomonic about the template centre
        r0, d0 = np.deg2rad(self.ra_c), np.deg2rad(self.dec_c)
        r, d = np.deg2rad(ra), np.deg2rad(dec)
        cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
        xi = np.rad2deg(np.cos(d) * np.sin(r - r0) / cosc) * 3600
        eta = np.rad2deg((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0)) / cosc) * 3600
        x = xi / self.sp + (self.n - 1) / 2
        y = eta / self.sp + (self.n - 1) / 2
        out = np.full(ra.shape, np.nan)
        vs = np.full(ra.shape, self.vscale_default)
        ok = (x >= 0) & (x <= self.n - 1.001) & (y >= 0) & (y <= self.n - 1.001)
        if ok.any():
            x0 = np.floor(x[ok]).astype(int); y0 = np.floor(y[ok]).astype(int)
            fx, fy = x[ok] - x0, y[ok] - y0

            def bil(arr):
                return (arr[y0, x0] * (1 - fx) * (1 - fy) + arr[y0, x0 + 1] * fx * (1 - fy)
                        + arr[y0 + 1, x0] * (1 - fx) * fy + arr[y0 + 1, x0 + 1] * fx * fy)
            out[ok] = bil(self.a) + bil(self.b) * (lam - self.lam0)
            vs[ok] = bil(self.vscale)
        return (out, vs) if with_vscale else out

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "spherex" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "spherex" / "precise_v1"
OUT_DIR = REPO / "runs" / "spherex" / os.environ.get("SPHEREX_CALIB_RUN", "calib_v1") / "tensors"
TEMPLATE_DIR = REPO / "runs" / "spherex" / os.environ.get("SPHEREX_CALIB_RUN", "calib_v1") / "templates"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"

NZ = 96
Z_GRID = 1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, NZ)
MU_GRID = np.array([-1.0, 0.0, 1.0])
T0_MJD = 61000.0
ZP_REF = 23.9  # AB zero point for uJy
LOCUS_BIN_DAYS = 0.05
# offsets (arcsec): >= 3 PSF FWHM from the track and from each other,
# inside the 615" cutout (locus half-span <= 190"). 16 controls since
# v3 (2026-08-21): FAR per search < 1/16 under the max-of-controls rule.
OFFSETS = [(0.0, 0.0), (40.0, 0.0), (-40.0, 0.0), (60.0, 0.0),
           (-60.0, 0.0), (80.0, 0.0), (-80.0, 0.0), (0.0, 50.0),
           (0.0, -50.0), (50.0, 50.0), (-50.0, 50.0), (50.0, -50.0),
           (-50.0, -50.0), (0.0, 75.0), (0.0, -75.0), (20.0, -60.0),
           (-20.0, 60.0)]


def main() -> None:
    only_missing = "--only-missing" in sys.argv[1:]
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.spherex_default()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(r["observation_id"])
    cut_index = {}
    with (PRECISE_DIR / "records" / "cutout_index.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            if rec["available"]:
                cut_index[rec["observation_id"]] = rec

    corridor_frames = defaultdict(set)
    for (e, role), oids in usable.items():
        corridor_frames[CORRIDOR_OF[e]].update(oids)
    phase_ref = {}
    for c, oids in corridor_frames.items():
        doys = np.array(sorted(obs_by_id[o]["t_mid_mjd_utc"] % 365.25
                               for o in oids))
        ang = doys / 365.25 * 2 * np.pi
        phase_ref[c] = float((np.arctan2(np.sin(ang).mean(),
                                         np.cos(ang).mean())
                              % (2 * np.pi)) / (2 * np.pi) * 365.25)

    def phase_of(corridor, mjd):
        d = (mjd % 365.25) - phase_ref[corridor]
        d = (d + 182.625) % 365.25 - 182.625
        return 0 if abs(d) < 91.3 else 1

    zcache: dict = {}

    def points_at_zgrid(endpoint, role, mjd):
        key = (endpoint, role, round(mjd / LOCUS_BIN_DAYS))
        if key not in zcache:
            al = adaptive_locus(
                target=registry[endpoint], role=Role(role),
                observation_time=Time(key[2] * LOCUS_BIN_DAYS, format="mjd"),
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=0.5, ephemeris=ctx.ephemeris,
                model=ctx.model)
            zs = np.array([p.z_au for p in al.points])
            ra = np.array([p.icrs_ra_deg for p in al.points])
            dec = np.array([p.icrs_dec_deg for p in al.points])
            q = 1.0 / zs
            o = np.argsort(q)
            qg = 1.0 / Z_GRID
            zcache[key] = np.stack([np.interp(qg, q[o], ra[o]),
                                    np.interp(qg, q[o], dec[o])], axis=1)
        return zcache[key]

    pairs_of_frame = defaultdict(list)
    for p, oids in usable.items():
        if only_missing and any(
                (OUT_DIR / f"{p[0]}__{p[1]}__D{d}.npz").exists()
                for d in range(1, 7)):
            continue
        for o in oids:
            pairs_of_frame[o].append(p)
    frames_by_corridor = defaultdict(list)
    for o, pairs in pairs_of_frame.items():
        frames_by_corridor[CORRIDOR_OF[pairs[0][0]]].append(o)
    nz, nm, nt = len(Z_GRID), len(MU_GRID), len(OFFSETS)
    t0 = _time.monotonic()
    n = 0
    for corridor, frames in sorted(frames_by_corridor.items()):
        frames.sort(key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
        rows = defaultdict(list)
        zcache.clear()
        templates = {}
        for tp in TEMPLATE_DIR.glob(f"{corridor}__D*.npz"):
            templates[tp.stem.split("__")[1]] = StaticTemplate(tp)
        print(f"[{corridor}] templates: {sorted(templates)}", flush=True)
        for oid in frames:
            ci = cut_index.get(oid)
            obs = obs_by_id[oid]
            if ci is None:
                continue
            path = REPO / ci["path"]
            try:
                fm = build_flux_map_spherex(
                    path, SpherexExactFootprint.FATAL_MASK, obs["band"],
                    obs["t_mid_mjd_utc"])
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", AstropyWarning)
                    with fits.open(path) as hdul:
                        hdr = hdul["IMAGE"].header
                        ny, nx = hdul["IMAGE"].data.shape
                        wave, bw = wavelength_at(hdul["WCS-WAVE"].data, hdr,
                                                 nx / 2.0, ny / 2.0)
                        psf_fwhm = float(hdr.get("PSF_FWHM", np.nan))
            except Exception as exc:
                print(f"  fluxmap FAIL {oid}: {exc}", flush=True)
                continue
            n += 1
            mjd = obs["t_mid_mjd_utc"]
            dt_yr = (mjd - T0_MJD) / 365.25
            deep = int(obs["collection"].endswith("_deep"))
            for pair in pairs_of_frame[oid]:
                endpoint, role = pair
                ph = phase_of(CORRIDOR_OF[endpoint], mjd)
                base = points_at_zgrid(endpoint, role, mjd)
                cosd = np.cos(np.deg2rad(base[:, 1]))
                dmu = MU_GRID * dt_yr / 3600.0
                ra_all = np.broadcast_to(
                    base[:, 0][:, None, None]
                    + dmu[None, :, None] / cosd[:, None, None], (nz, nm, nm))
                dec_all = np.broadcast_to(
                    base[:, 1][:, None, None] + dmu[None, None, :], (nz, nm, nm))
                F = np.empty((nt, nz, nm, nm), dtype=np.float32)
                V = np.empty_like(F)
                G = np.empty((nt, nz, nm, nm), dtype=np.float16)
                TS = np.empty((nt, nz, nm, nm), dtype=np.float16)
                tmpl = templates.get(obs["band"])
                for ti, (dra, ddec) in enumerate(OFFSETS):
                    ra_q = ra_all + dra / 3600.0 / cosd[:, None, None]
                    dec_q = dec_all + ddec / 3600.0
                    f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
                    if tmpl is not None:
                        t, vs = tmpl(ra_q.ravel(), dec_q.ravel(), wave, with_vscale=True)
                        has = np.isfinite(t)
                        f = np.where(has, f - np.nan_to_num(t), f)
                        v = v * vs
                    else:
                        has = np.zeros(f.shape, dtype=bool)
                    F[ti] = f.reshape(nz, nm, nm).astype(np.float32)
                    V[ti] = v.reshape(nz, nm, nm).astype(np.float32)
                    G[ti] = g.reshape(nz, nm, nm).astype(np.float16)
                    TS[ti] = has.reshape(nz, nm, nm).astype(np.float16)
                rows[(endpoint, role, obs["band"])].append(
                    (mjd, ph, wave, bw, psf_fwhm, fm.var_scale, deep, F, V, G, TS))
            if n % 200 == 0:
                print(f"  {n} maps ({n / (_time.monotonic() - t0):.1f}/s)",
                      flush=True)
        for (endpoint, role, band), rs in rows.items():
            rs.sort(key=lambda r: r[0])
            np.savez_compressed(
                OUT_DIR / f"{endpoint}__{role}__{band}.npz",
                z_grid=Z_GRID, mu_grid=MU_GRID, t0_mjd=T0_MJD, zp_ref=ZP_REF,
                offsets=np.array(OFFSETS),
                mjd=np.array([r[0] for r in rs]),
                phase=np.array([r[1] for r in rs], dtype=np.uint8),
                wave_um=np.array([r[2] for r in rs], dtype=np.float32),
                bandwidth_um=np.array([r[3] for r in rs], dtype=np.float32),
                psf_fwhm=np.array([r[4] for r in rs], dtype=np.float32),
                var_scale=np.array([r[5] for r in rs], dtype=np.float32),
                deep=np.array([r[6] for r in rs], dtype=np.uint8),
                f=np.stack([r[7] for r in rs], axis=1),
                v=np.stack([r[8] for r in rs], axis=1),
                g=np.stack([r[9] for r in rs], axis=1),
                tsub=np.stack([r[10] for r in rs], axis=1))
            print(f"wrote {endpoint}__{role}__{band}: {len(rs)} epochs "
                  f"[{corridor}]", flush=True)
        del rows
    print(f"{n} flux maps built")


if __name__ == "__main__":
    main()
