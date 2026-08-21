"""Layer-1 screening for the SPHEREx pilot (plan §6 step 4).

SPHEREx QR2 publishes no source catalogs, so the screening layer is
built from the images: for every usable cutout, find single-epoch
matched-filter peaks with S_e > CLIP within SCREEN_RADIUS of each
endpoint x role locus at that epoch and record them as ScreenMatch
records (distance to track, implied relay distance, flux, wavelength).
Then cluster the matches by sky position: a position recurring across
epochs is a static source (cross-matched against CatWISE2020,
snapshotted); what remains is the track-following residue handed to
the stack layer's candidate census. Catalog absence is neither a
detection nor a null (plan §3.4).

Usage: uv run python surveys/spherex/scripts/track_screen.py
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
import time as _time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.irsa_spherex import (SpherexExactFootprint,
                                             wavelength_at)
from sglsurvey.geometry import GeometryContext
from sglsurvey.photometry import build_flux_map_spherex
from sglsurvey.records import ScreenMatch, append_records, read_records
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spherex_corridors import CORRIDOR_OF, MEMBERS  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "spherex" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "spherex" / "precise_v1"
RUN_DIR = REPO / "runs" / "spherex" / "screen_v1"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "spherex" / "hypotheses.md"
HYPOTHESIS_VERSION = "spherex-hypotheses-v1.0"
TAP_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"

SCREEN_RADIUS_ARCSEC = 10.0
CLIP_SIGMA = 5.0
LOCUS_TOLERANCE_ARCSEC = 1.0
LOCUS_BIN_DAYS = 0.05
CLUSTER_ARCSEC = 6.0   # one pixel
STATIC_MIN_EPOCHS = 3
CATWISE_MATCH_ARCSEC = 6.0


def seg_min_dist_vec(det_ra, det_dec, pts):
    cosd = np.cos(np.deg2rad(det_dec))[:, None]
    x = (pts[None, :, 0] - det_ra[:, None]) * cosd * 3600.0
    y = (pts[None, :, 1] - det_dec[:, None]) * 3600.0
    ax, ay, bx, by = x[:, :-1], y[:, :-1], x[:, 1:], y[:, 1:]
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = np.clip(np.where(L2 > 0, -(ax * dx + ay * dy) / np.where(L2 > 0, L2, 1.0), 0.0), 0.0, 1.0)
    d = np.hypot(ax + t * dx, ay + t * dy)
    i = np.argmin(d, axis=1)
    return d[np.arange(len(d)), i], i


def local_peaks(S, thresh):
    """Indices (y, x) of local maxima above thresh in a 2-D map."""
    from scipy.ndimage import maximum_filter

    Sn = np.where(np.isfinite(S), S, -np.inf)
    mx = maximum_filter(Sn, size=5, mode="nearest")
    ys, xs = np.nonzero((Sn == mx) & (Sn > thresh))
    return ys, xs


def main() -> None:
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.spherex_default()
    store = SnapshotStore(RUN_DIR)
    session = requests.Session()
    hyp_hash = "sha256:" + hashlib.sha256(HYPOTHESES_PATH.read_bytes()).hexdigest()

    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(PRECISE_DIR / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[r["observation_id"]].append((r["endpoint_id"], r["role"]))
    cut_index = {}
    with (PRECISE_DIR / "records" / "cutout_index.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            if rec["available"]:
                cut_index[rec["observation_id"]] = rec

    RUN_DIR.mkdir(parents=True, exist_ok=True)
    run_config = {
        "run_id": "screen_v1", "stage": "track-screen",
        "source_runs": ["coarse_v1", "precise_v1"],
        "screen_radius_arcsec": SCREEN_RADIUS_ARCSEC,
        "single_epoch_sigma": CLIP_SIGMA,
        "cluster_arcsec": CLUSTER_ARCSEC, "static_min_epochs": STATIC_MIN_EPOCHS,
        "context_catalog": "catwise_2020",
        "registry_source_hash": registry.source_hash,
        "hypothesis_version": HYPOTHESIS_VERSION, "hypothesis_hash": hyp_hash,
        "geometry": ctx.identities(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    (RUN_DIR / f"run_config_{run_config['started_utc'][:19]}.json"
     ).write_text(json.dumps(run_config, indent=2))
    match_path = RUN_DIR / "records" / "screen_match.jsonl"
    if match_path.exists():
        match_path.unlink()

    zcache = {}

    def locus(endpoint, role, mjd):
        key = (endpoint, role, round(mjd / LOCUS_BIN_DAYS))
        if key not in zcache:
            al = adaptive_locus(
                target=registry[endpoint], role=Role(role),
                observation_time=Time(key[2] * LOCUS_BIN_DAYS, format="mjd"),
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=LOCUS_TOLERANCE_ARCSEC,
                ephemeris=ctx.ephemeris, model=ctx.model)
            zcache[key] = (np.array([[p.icrs_ra_deg, p.icrs_dec_deg] for p in al.points]),
                           np.array([p.z_au for p in al.points]))
        return zcache[key]

    summary = {}
    for corridor, members in MEMBERS.items():
        oids = sorted((o for o, pairs in usable.items()
                       if o in cut_index and any(e in members for e, _ in pairs)),
                      key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
        if not oids:
            continue
        zcache.clear()
        print(f"[{corridor}] {len(oids)} usable cutouts", flush=True)
        t0 = _time.monotonic()
        matches = []
        n_peaks = 0
        for i, oid in enumerate(oids):
            obs = obs_by_id[oid]
            path = REPO / cut_index[oid]["path"]
            try:
                fm = build_flux_map_spherex(path, SpherexExactFootprint.FATAL_MASK,
                                            obs["band"], obs["t_mid_mjd_utc"])
                from astropy.io import fits
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    with fits.open(path) as hdul:
                        hdr = hdul["IMAGE"].header
                        wave_tab = hdul["WCS-WAVE"].data
            except Exception as exc:
                print(f"  fluxmap FAIL {oid}: {exc}", flush=True)
                continue
            with np.errstate(invalid="ignore", divide="ignore"):
                S = fm.flux / np.sqrt(fm.var)
            S[fm.good_frac < 0.7] = np.nan
            ys, xs = local_peaks(S, CLIP_SIGMA)
            if len(ys) == 0:
                continue
            up = fm.upsample
            px, py = xs / up, ys / up
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                world = fm.wcs.all_pix2world(np.stack([px, py], axis=1), 0)
            det_ra, det_dec = world[:, 0], world[:, 1]
            n_peaks += len(ys)
            for e, r in usable[oid]:
                if e not in members:
                    continue
                pts, zs = locus(e, r, obs["t_mid_mjd_utc"])
                d, seg = seg_min_dist_vec(det_ra, det_dec, pts)
                for k in np.flatnonzero(d <= SCREEN_RADIUS_ARCSEC):
                    lam, bw = wavelength_at(wave_tab, hdr, px[k], py[k])
                    matches.append(ScreenMatch.build(
                        endpoint_id=e, role=r,
                        hypothesis_version=HYPOTHESIS_VERSION,
                        registry_source_hash=registry.source_hash,
                        source_table="spherex-l2-matched-filter-peaks",
                        source_cntr=f"{obs['native_key']['obs_id']}D{obs['native_key']['detector']}:{int(ys[k])}:{int(xs[k])}",
                        source_designation=None,
                        ra_deg=float(det_ra[k]), dec_deg=float(det_dec[k]),
                        sigra_mas=None, sigdec_mas=None,
                        mjd=obs["t_mid_mjd_utc"], scan_id=None, frame_num=None,
                        photometry={"band": obs["band"], "wave_um": round(lam, 4),
                                    "flux_ujy": round(float(fm.flux[ys[k], xs[k]]), 1),
                                    "snr": round(float(S[ys[k], xs[k]]), 1),
                                    "mag_ab": round(23.9 - 2.5 * np.log10(max(fm.flux[ys[k], xs[k]], 1e-3)), 2)},
                        flags={"good_frac": round(float(fm.good_frac[ys[k], xs[k]]), 3)},
                        dist_arcsec=round(float(d[k]), 2),
                        implied_z_au=round(float(zs[seg[k]]), 1),
                        z_segment_au=(round(float(zs[seg[k]]), 1), round(float(zs[seg[k] + 1]), 1)),
                        screen_radius_arcsec=SCREEN_RADIUS_ARCSEC,
                        locus_mjd=obs["t_mid_mjd_utc"], snapshot_id=None,
                        extra={"observation_id": oid, "corridor": corridor}))
            if (i + 1) % 500 == 0:
                print(f"  {i + 1}/{len(oids)} ({(i + 1) / (_time.monotonic() - t0):.1f}/s), "
                      f"{len(matches)} matches", flush=True)
        append_records(match_path, matches)

        # --- recurrence: cluster by position ---------------------------
        clusters = []
        for m in matches:
            for c in clusters:
                if (np.hypot((m.ra_deg - c["ra"]) * np.cos(np.deg2rad(c["dec"])),
                             m.dec_deg - c["dec"]) * 3600 < CLUSTER_ARCSEC):
                    c["members"].append(m)
                    n = len(c["members"])
                    c["ra"] += (m.ra_deg - c["ra"]) / n
                    c["dec"] += (m.dec_deg - c["dec"]) / n
                    break
            else:
                clusters.append({"ra": m.ra_deg, "dec": m.dec_deg, "members": [m]})
        # CatWISE context
        if matches:
            cra = float(np.mean([m.ra_deg for m in matches]))
            cdec = float(np.mean([m.dec_deg for m in matches]))
            q = ("SELECT source_name, ra, dec, w1mpro, w2mpro FROM catwise_2020 WHERE "
                 f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{cra:.6f},{cdec:.6f},0.15))=1")
            request_utc = datetime.now(timezone.utc).isoformat()
            resp = session.get(TAP_SYNC, params={"QUERY": q, "FORMAT": "CSV"}, timeout=900)
            cw = list(csv.DictReader(io.StringIO(resp.text))) if resp.ok else []
            store.store(service_url=TAP_SYNC, query=q, request_utc=request_utc,
                        response_bytes=resp.content, row_count=len(cw),
                        http_status=resp.status_code)
            cwra = np.array([float(r["ra"]) for r in cw]) if cw else np.zeros(0)
            cwdec = np.array([float(r["dec"]) for r in cw]) if cw else np.zeros(0)
        else:
            cw, cwra, cwdec = [], np.zeros(0), np.zeros(0)
        static, moving = [], []
        for c in clusters:
            epochs = sorted({round(m.mjd, 1) for m in c["members"]})
            c["n_epochs"] = len(epochs)
            c["bands"] = sorted({m.photometry["band"] for m in c["members"]})
            c["catwise"] = None
            if len(cwra):
                dd = np.hypot((cwra - c["ra"]) * np.cos(np.deg2rad(c["dec"])), cwdec - c["dec"]) * 3600
                j = int(np.argmin(dd))
                if dd[j] < CATWISE_MATCH_ARCSEC:
                    c["catwise"] = {"name": cw[j]["source_name"], "sep_arcsec": round(float(dd[j]), 2),
                                    "w1mpro": cw[j]["w1mpro"], "w2mpro": cw[j]["w2mpro"]}
            (static if (c["n_epochs"] >= STATIC_MIN_EPOCHS or c["catwise"]) else moving).append(c)
        summary[corridor] = {
            "usable_cutouts": len(oids), "peaks_total": n_peaks,
            "matches_within_radius": len(matches), "clusters": len(clusters),
            "static_clusters": len(static), "unresolved_clusters": len(moving),
            "static": [{"ra": round(c["ra"], 5), "dec": round(c["dec"], 5),
                        "n_matches": len(c["members"]), "n_epochs": c["n_epochs"],
                        "bands": c["bands"], "catwise": c["catwise"],
                        "median_mag_ab": round(float(np.median(
                            [m.photometry["mag_ab"] for m in c["members"]])), 2)}
                       for c in sorted(static, key=lambda c: -len(c["members"]))],
            "unresolved": [{"ra": round(c["ra"], 5), "dec": round(c["dec"], 5),
                            "n_matches": len(c["members"]), "n_epochs": c["n_epochs"],
                            "bands": c["bands"],
                            "snr": [m.photometry["snr"] for m in c["members"]],
                            "implied_z_au": [m.implied_z_au for m in c["members"]],
                            "mjd": [round(m.mjd, 2) for m in c["members"]]}
                           for c in moving],
        }
        print(f"[{corridor}] {n_peaks} peaks, {len(matches)} within {SCREEN_RADIUS_ARCSEC}\", "
              f"{len(clusters)} clusters: {len(static)} static, {len(moving)} unresolved",
              flush=True)
    (RUN_DIR / "summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
