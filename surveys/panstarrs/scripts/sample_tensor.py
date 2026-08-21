"""Stage-7 sampling pass for the Pan-STARRS1 pilot: per-epoch trajectory
sample tensors for every endpoint x role, real trajectory plus 8 offset
controls (same design as surveys/ztf/scripts/sample_tensor.py).

Writes runs/panstarrs/calib_v1/tensors/<endpoint>__<role>.npz:
  mjd (E,), band_idx (E,) [1 g .. 5 y], phase (E,) [0/1], seeing (E,)
  arcsec, zp_star (E,), zp_hdr (E,), n_cal (E,), var_source (E,),
  f/v (9, E, NZ, 5, 5) float32 at the common ZP_REF = 25, g float16.

PS1-specific choices (hypotheses v1.0 section 9):
  * NZ = 360 nodes uniform in 1/z (~1" spacing, under the seeing);
  * T0 = MJD 56000 (3pi baseline midpoint);
  * flux scale: per-warp empirical zero point from DR2 ``mean``-table
    stars recovered through the same matched filter in the same cutout
    (ZP_star = median(m_DR2 + 2.5 log10 F_mf)); header FPA.ZP recorded
    alongside; frames with < MIN_CAL_STARS calibrators use the filter's
    median (ZP_star - FPA.ZP) offset;
  * duplicate epochs (one exposure warped onto two overlapping
    skycells) collapsed to the warp with the larger usable fraction.

Usage: uv run python surveys/panstarrs/scripts/sample_tensor.py [--only-missing]
"""

from __future__ import annotations

import csv
import io
import json
import sys
import time as _time
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.mast_ps1 import PIX_ARCSEC, Ps1ExactFootprint
from sglsurvey.geometry import GeometryContext
from sglsurvey.photometry import build_flux_map_ps1
from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ps1_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "panstarrs" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "panstarrs" / "precise_v1"
SCREEN_DIR = REPO / "runs" / "panstarrs" / "screen_v1"
CUT_DIR = REPO / "runs" / "panstarrs" / "products" / "cut"
MSK_DIR = REPO / "runs" / "panstarrs" / "products" / "msk"
OUT_DIR = REPO / "runs" / "panstarrs" / "calib_v1" / "tensors"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"

NZ = 360
Z_GRID = 1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, NZ)
MU_GRID = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
T0_MJD = 56000.0
ZP_REF = 25.0
LOCUS_BIN_DAYS = 0.05
EXPOSURE_DUP_DAYS = 0.001
OFFSETS = [(0.0, 0.0), (20.0, 0.0), (-20.0, 0.0), (30.0, 0.0),
           (-30.0, 0.0), (40.0, 0.0), (-40.0, 0.0), (0.0, 25.0),
           (0.0, -25.0)]
BAND_IDX = {"g": 1, "r": 2, "i": 3, "z": 4, "y": 5}
# star calibration
CAL_MAG_RANGE = (15.0, 20.5)
CAL_MIN_NPT = 3
CAL_MAX_ERR = 0.10
CAL_MIN_SNR = 10.0
CAL_MIN_GOOD_FRAC = 0.90
MIN_CAL_STARS = 5
VAR_SOURCE_IDX = {"wt": 1, "1/wt": 2, "robust": 3}


def load_mean_stars(corridor: str) -> dict[str, np.ndarray]:
    stats = json.loads((SCREEN_DIR / "catalog_stats.json").read_text())
    snaps = {s["snapshot_id"]: s for s in read_records(
        SCREEN_DIR / "records" / "query_snapshot.jsonl")}
    rows = []
    for sid in stats[corridor]["mean_snapshots"]:
        p = SCREEN_DIR / snaps[sid]["response_path"]
        rows.extend(csv.DictReader(io.StringIO(p.read_text())))
    out = {"ra": np.array([float(r["raMean"]) for r in rows]),
           "dec": np.array([float(r["decMean"]) for r in rows]),
           "ndet": np.array([int(float(r["nDetections"])) for r in rows])}
    for b in BAND_IDX:
        out[f"{b}_mag"] = np.array([float(r[f"{b}MeanPSFMag"]) for r in rows])
        out[f"{b}_err"] = np.array([float(r[f"{b}MeanPSFMagErr"]) for r in rows])
        out[f"{b}_npt"] = np.array([int(float(r[f"{b}MeanPSFMagNpt"]))
                                    for r in rows])
    return out


def star_zeropoint(fm, stars: dict, band: str):
    """Median zero point from catalogued stars in the cutout; returns
    (zp, n_used, mad)."""
    m = stars[f"{band}_mag"]
    sel = ((m > CAL_MAG_RANGE[0]) & (m < CAL_MAG_RANGE[1])
           & (stars[f"{band}_npt"] >= CAL_MIN_NPT)
           & (stars[f"{band}_err"] < CAL_MAX_ERR) & (stars[f"{band}_err"] > 0))
    if not sel.any():
        return None, 0, None
    f, v, g = fm.sample(stars["ra"][sel], stars["dec"][sel])
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = f / np.sqrt(v)
    ok = (np.isfinite(f) & (f > 0) & (g >= CAL_MIN_GOOD_FRAC)
          & np.isfinite(snr) & (snr > CAL_MIN_SNR))
    if ok.sum() == 0:
        return None, 0, None
    zp = m[sel][ok] + 2.5 * np.log10(f[ok])
    med = float(np.median(zp))
    mad = float(1.4826 * np.median(np.abs(zp - med)))
    # one round of 3-sigma clipping against blends / variables
    keep = np.abs(zp - med) < max(3 * mad, 0.05)
    if keep.sum() >= 3:
        zp = zp[keep]
        med = float(np.median(zp))
        mad = float(1.4826 * np.median(np.abs(zp - med)))
    return med, int(len(zp)), mad


CAL_CUT_DIR = REPO / "runs" / "panstarrs" / "products" / "calcut"
CAL_CUT_PIX = 1600  # 400" box: ~40-100 stars at 15 < r < 20.5 at high latitude
_cal_adapter = None


def calibration_flux_map(oid, obs, manifest_row, corridor):
    """Image-only fitscut cutout at the corridor antipode (overlay centre)
    for the per-warp star zero point when the locus cutout is too sparse.
    Returns a FluxMap on the same warp (robust variance) or None."""
    global _cal_adapter
    from sglsurvey.adapters.base import CutoutSpec
    from sglsurvey.adapters.mast_ps1 import Ps1WarpAdapter
    from sglsurvey.records import Observation
    from ps1_corridors import OVERLAY
    ov = OVERLAY.get(corridor)
    if ov is None:
        return None
    if _cal_adapter is None:
        _cal_adapter = Ps1WarpAdapter()
    o = dict(obs)
    o["corners_icrs_deg"] = tuple(map(tuple, o["corners_icrs_deg"]))
    # centre on the locus cutout (always inside this warp's skycell; the
    # corridor antipode may sit in a neighbouring skycell -> fitscut 400)
    cut = CutoutSpec(manifest_row["center_ra_deg"],
                     manifest_row["center_dec_deg"], CAL_CUT_PIX)
    for attempt in (1, 2, 3):
        try:
            ps = _cal_adapter.fetch(Observation(**o), ["img"], CAL_CUT_DIR,
                                    cutout=cut)
            break
        except FileNotFoundError:
            return None
        except Exception:
            if attempt == 3:
                raise
            _time.sleep(2.0 * attempt)
    return build_flux_map_ps1(ps.products[0].path, None,
                              MSK_DIR / manifest_row["msk"],
                              Ps1ExactFootprint.FATAL_MASK, obs["band"],
                              obs["t_mid_mjd_utc"])


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only-missing", action="store_true")
    ap.add_argument("--corridors", nargs="*", default=None)
    args = ap.parse_args()
    only_missing = args.only_missing
    only = set(args.corridors) if args.corridors else None
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.ps1_default()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    usable = defaultdict(dict)  # pair -> {obs_id: usable_fraction}
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])][r["observation_id"]] = (
                r["usable_fraction"] or 0.0)
    manifest = {}
    with (CUT_DIR / "manifest.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            manifest[rec["observation_id"]] = rec

    # collapse skycell duplicates: same (band, exposure MJD) -> best frac
    dropped = 0
    for pair, frames in usable.items():
        groups = defaultdict(list)
        for oid, frac in frames.items():
            o = obs_by_id[oid]
            groups[(o["band"], round(o["t_start_mjd_utc"] / EXPOSURE_DUP_DAYS))
                   ].append((frac, oid))
        keep = {}
        for g in groups.values():
            g.sort(reverse=True)
            keep[g[0][1]] = g[0][0]
            dropped += len(g) - 1
        usable[pair] = keep
    print(f"{dropped} skycell-duplicate warps collapsed", flush=True)

    corridor_frames = defaultdict(set)
    for (e, role), frames in usable.items():
        corridor_frames[CORRIDOR_OF[e]].update(frames)
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
    for p, frames in usable.items():
        if only is not None and CORRIDOR_OF[p[0]] not in only:
            continue
        if only_missing and (OUT_DIR / f"{p[0]}__{p[1]}.npz").exists():
            continue
        for o in frames:
            pairs_of_frame[o].append(p)
    frames_by_corridor = defaultdict(list)
    for o, pairs in pairs_of_frame.items():
        frames_by_corridor[CORRIDOR_OF[pairs[0][0]]].append(o)
    nz, nm, nt = len(Z_GRID), len(MU_GRID), len(OFFSETS)
    t0 = _time.monotonic()
    n = 0
    zp_log = []
    all_rows = {}
    # merge with earlier batches' zero points (fallback uses everything)
    zp_path = OUT_DIR.parent / "zeropoints.jsonl"
    batch_corridors = set(frames_by_corridor)
    prev_log = []
    if zp_path.exists():
        prev_log = [json.loads(l) for l in zp_path.read_text().splitlines()
                    if l.strip()]
        prev_log = [r for r in prev_log if r["corridor"] not in batch_corridors]
    for corridor, frames in sorted(frames_by_corridor.items()):
        frames.sort(key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
        stars = load_mean_stars(corridor)
        rows = defaultdict(list)
        zcache.clear()
        for oid in frames:
            row = manifest.get(oid)
            obs = obs_by_id[oid]
            if row is None or "img" not in row.get("files", {}) or not row["msk"]:
                continue
            wt = (CUT_DIR / row["files"]["wt"]
                  if "wt" in row["files"] else None)
            try:
                fm = build_flux_map_ps1(
                    CUT_DIR / row["files"]["img"], wt, MSK_DIR / row["msk"],
                    Ps1ExactFootprint.FATAL_MASK, obs["band"],
                    obs["t_mid_mjd_utc"])
            except Exception as exc:
                print(f"  fluxmap FAIL {oid}: {exc}", flush=True)
                continue
            zp_hdr = fm.magzp
            zp_star, n_cal, zp_mad = star_zeropoint(fm, stars, obs["band"])
            cal_source = "locus-cutout"
            if n_cal < MIN_CAL_STARS:
                # dedicated calibration cutout (image only, robust variance)
                # at the corridor antipode; lazily fetched, purged per batch
                try:
                    cfm = calibration_flux_map(oid, obs, row, corridor)
                    if cfm is not None:
                        z2, n2, m2 = star_zeropoint(cfm, stars, obs["band"])
                        if n2 > n_cal:
                            zp_star, n_cal, zp_mad = z2, n2, m2
                            cal_source = "calibration-cutout"
                except Exception as exc:
                    print(f"  calcut FAIL {oid}: {exc}", flush=True)
            zp_log.append({"observation_id": oid, "corridor": corridor,
                           "band": obs["band"], "zp_hdr": zp_hdr,
                           "zp_star": zp_star, "n_cal": n_cal,
                           "zp_mad": zp_mad, "var_source": fm.var_source,
                           "cal_source": cal_source,
                           "fwhm_pix": fm.fwhm_pix, "exptime": fm.exptime})
            n += 1
            # exposure time from the header, mid-exposure epoch
            mjd = fm.mjd_obs + (fm.exptime or 0.0) / 2.0 / 86400.0
            dt_yr = (mjd - T0_MJD) / 365.25
            seeing_arcsec = fm.fwhm_pix * PIX_ARCSEC
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
                for ti, (dra, ddec) in enumerate(OFFSETS):
                    ra_q = ra_all + dra / 3600.0 / cosd[:, None, None]
                    dec_q = dec_all + ddec / 3600.0
                    f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
                    F[ti] = f.reshape(nz, nm, nm).astype(np.float32)
                    V[ti] = v.reshape(nz, nm, nm).astype(np.float32)
                    G[ti] = g.reshape(nz, nm, nm).astype(np.float16)
                # raw DN for now; scaled to ZP_REF once ZP is final
                rows[pair].append([mjd, BAND_IDX[obs["band"]], ph,
                                   seeing_arcsec, zp_hdr, zp_star, n_cal,
                                   VAR_SOURCE_IDX[fm.var_source], F, V, G])
            if n % 50 == 0:
                print(f"  {n} maps ({n / (_time.monotonic() - t0):.2f}/s)",
                      flush=True)
        all_rows.update(rows)
        del rows
    # filter-level fallback offset (ALL corridors) for frames with too few
    # calibrators; sparse high-latitude cutouts (ross128) may have none.
    offs = defaultdict(list)
    for r in prev_log + zp_log:
        if r["n_cal"] >= MIN_CAL_STARS and r["zp_hdr"] is not None:
            offs[r["band"]].append(r["zp_star"] - r["zp_hdr"])
    fallback = {b: float(np.median(v)) for b, v in offs.items()}
    print(f"fallback ZP_star - FPA.ZP by filter: "
          f"{ {b: round(v, 3) for b, v in fallback.items()} }", flush=True)
    # non-photometric guard for fallback frames: header ZP must sit within
    # HDR_ZP_TOL of the filter's median header ZP, else the frame is dropped
    HDR_ZP_TOL = 0.5
    hdr_med = {}
    for b in BAND_IDX.values():
        vals = [r[4] for rs in all_rows.values() for r in rs
                if r[1] == b and r[4] is not None]
        if vals:
            hdr_med[b] = float(np.median(vals))
    n_dropped = 0
    if True:
        for (endpoint, role), rs in all_rows.items():
            corridor = CORRIDOR_OF[endpoint]
            rs.sort(key=lambda r: r[0])
            keep_rows = []
            for r in rs:
                if (r[5] is None or r[6] < MIN_CAL_STARS) and (
                        r[4] is None
                        or abs(r[4] - hdr_med.get(r[1], r[4])) > HDR_ZP_TOL):
                    n_dropped += 1
                    continue
                keep_rows.append(r)
            rs[:] = keep_rows
            for r in rs:
                zp_hdr, zp_star, n_cal = r[4], r[5], r[6]
                if zp_star is None or n_cal < MIN_CAL_STARS:
                    band = {v: k for k, v in BAND_IDX.items()}[r[1]]
                    if zp_hdr is None or band not in fallback:
                        raise RuntimeError(
                            f"no zero point for {endpoint}/{role} band {band}")
                    zp_use = zp_hdr + fallback[band]
                    r[5] = zp_use
                    r[6] = -max(n_cal, 0)  # negative n_cal marks fallback
                scale = 10.0 ** ((ZP_REF - r[5]) / 2.5)
                r[8] = (r[8] * scale).astype(np.float32)
                r[9] = (r[9] * scale * scale).astype(np.float32)
            np.savez_compressed(
                OUT_DIR / f"{endpoint}__{role}.npz",
                z_grid=Z_GRID, mu_grid=MU_GRID, t0_mjd=T0_MJD, zp_ref=ZP_REF,
                offsets=np.array(OFFSETS),
                mjd=np.array([r[0] for r in rs]),
                band_idx=np.array([r[1] for r in rs], dtype=np.uint8),
                phase=np.array([r[2] for r in rs], dtype=np.uint8),
                seeing=np.array([r[3] for r in rs], dtype=np.float32),
                zp_hdr=np.array([r[4] if r[4] is not None else np.nan
                                 for r in rs], dtype=np.float32),
                zp_star=np.array([r[5] for r in rs], dtype=np.float32),
                n_cal=np.array([r[6] for r in rs], dtype=np.int16),
                var_source=np.array([r[7] for r in rs], dtype=np.uint8),
                f=np.stack([r[8] for r in rs], axis=1),
                v=np.stack([r[9] for r in rs], axis=1),
                g=np.stack([r[10] for r in rs], axis=1))
            print(f"wrote {endpoint}__{role}: {len(rs)} epochs "
                  f"[{corridor}]", flush=True)
    print(f"{n_dropped} fallback epochs dropped (header ZP off by > "
          f"{HDR_ZP_TOL} mag)")
    zp_path.write_text(
        "\n".join(json.dumps(r) for r in prev_log + zp_log) + "\n")
    good = [r for r in zp_log if r["n_cal"] >= MIN_CAL_STARS]
    print(f"{n} flux maps built; {len(good)} with >= {MIN_CAL_STARS} "
          f"calibrators")
    for b in BAND_IDX:
        d = [r["zp_star"] - r["zp_hdr"] for r in good
             if r["band"] == b and r["zp_hdr"] is not None]
        if d:
            print(f"  {b}: ZP_star - FPA.ZP median {np.median(d):+.3f} "
                  f"(MAD {1.4826 * np.median(np.abs(d - np.median(d))):.3f}, "
                  f"N={len(d)})")


if __name__ == "__main__":
    main()
