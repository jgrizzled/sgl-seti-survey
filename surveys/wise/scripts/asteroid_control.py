"""Step C: positive control — a numbered main-belt asteroid through the
full WISE v2 chain (hypotheses v2.0 §1.5; v2 plan §2.5).

The JPL Horizons ephemeris (observer = the WISE spacecraft, ``@-163``)
plays the trajectory model. Every NEOWISE/cryo frame containing the
predicted position is discovered through the production adapter,
cut out through IBE, turned into a flux map with the production code
(header pixel scale, kept inputs), sampled on a 5 x 5 residual-offset
grid (±2", the analogue of the mu grid) for the real trajectory, the
48-offset ring and the 8 designated controls, stacked with the v2
statistic, and scored with R = S_max / T and the ring null exactly as
a survey cell. Throughput: the recovered per-frame and stacked fluxes
are compared with the NEOWISE single-exposure source table's own
profile-fit photometry of the asteroid on the same frames (the
pipeline's w?mpro), and with the Horizons V prediction.

Outputs runs/wise/v2/control/<asteroid>/{frames.jsonl, samples.npz,
summary.json}.

Usage: uv run python surveys/wise/scripts/asteroid_control.py
           [--asteroid 14000] [--max-frames 400]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import time as _time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

from sglsurvey import nulls  # noqa: E402
from sglsurvey.adapters.base import ConeRegion, CutoutSpec, MjdRange  # noqa: E402
from sglsurvey.adapters.irsa_wise import (WiseExactFootprint,  # noqa: E402
                                          WiseMergeL1bAdapter)
from sglsurvey.photometry import build_flux_map  # noqa: E402
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
JD_MJD = 2400000.5
CUTOUT_PIX = 120
RESID_GRID = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])   # arcsec
ELONG_TOL = 4.0


def horizons(session, command, **kw):
    params = {"format": "text", "COMMAND": f"'{command};'", "OBJ_DATA": "NO",
              "MAKE_EPHEM": "YES", "EPHEM_TYPE": "OBSERVER", "CENTER": "'@-163'",
              "QUANTITIES": "'1,9,23'", "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES",
              "EXTRA_PREC": "YES"}
    params.update(kw)
    for attempt in range(4):
        try:
            r = session.get(HORIZONS, params=params, timeout=180)
            r.raise_for_status()
            if "$$SOE" in r.text:
                break
        except requests.RequestException:
            pass
        _time.sleep(5 * (attempt + 1))
    else:
        raise RuntimeError("Horizons failed")
    body = r.text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    rows = []
    for line in body:
        c = [x.strip() for x in line.split(",")]
        try:
            rows.append({"date": c[0], "ra": float(c[3]), "dec": float(c[4]),
                         "v": float(c[5]) if c[5] not in ("n.a.", "") else np.nan,
                         "elong": float(c[7]) if len(c) > 7 and c[7] not in ("n.a.", "") else np.nan})
        except (ValueError, IndexError):
            continue
    return rows, r.content, r.url


def _date_to_mjd(date: str) -> float:
    t = datetime.strptime(date[:17].strip(), "%Y-%b-%d %H:%M")
    return (t - datetime(1858, 11, 17)).total_seconds() / 86400.0


def tap_rows(session, query: str):
    for attempt in range(5):
        try:
            r = session.get(TAP, params={"QUERY": query, "FORMAT": "CSV"}, timeout=300)
            if r.status_code == 200:
                return list(csv.DictReader(io.StringIO(r.text))), r.content
        except requests.RequestException:
            pass
        _time.sleep(10 * (attempt + 1))
    return [], b""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--asteroid", default="14000")
    ap.add_argument("--max-frames", type=int, default=400)
    ap.add_argument("--refetch", action="store_true")
    a = ap.parse_args()
    run_dir = C.RUN_DIR / "control" / a.asteroid
    cut_dir = run_dir / "cut"
    run_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    store = SnapshotStore(run_dir)
    adapter = WiseMergeL1bAdapter(session=session)
    if (run_dir / "samples.npz").exists() and not a.refetch:
        recs = [json.loads(l) for l in (run_dir / "frames.jsonl").open()]
        s = np.load(run_dir / "samples.npz")
        analyse(a, run_dir, session, store, recs, s["f"], s["v"], s["g"], s["meta"])
        return

    # 1. daily ephemeris over the mission, WISE as observer
    daily, raw, url = horizons(session, a.asteroid, START_TIME="'2010-01-07'",
                               STOP_TIME="'2024-07-31'", STEP_SIZE="'1d'")
    store.store(service_url=url, query="daily", request_utc=datetime.now(timezone.utc).isoformat(),
                response_bytes=raw, row_count=len(daily), http_status=200)
    for d in daily:
        d["mjd"] = _date_to_mjd(d["date"])
    vis = [d for d in daily if np.isfinite(d["elong"]) and abs(d["elong"] - 90.0) <= ELONG_TOL]
    # group consecutive days into visit windows
    windows = []
    for d in vis:
        if windows and d["mjd"] - windows[-1][-1]["mjd"] <= 1.5:
            windows[-1].append(d)
        else:
            windows.append([d])
    print(f"{len(daily)} daily epochs, {len(vis)} within {ELONG_TOL} deg of quadrature, "
          f"{len(windows)} windows", flush=True)

    # 2. frames per window through the production adapter
    frames = []
    for w in windows:
        mid = w[len(w) // 2]
        ra_span = max(x["ra"] for x in w) - min(x["ra"] for x in w)
        de_span = max(x["dec"] for x in w) - min(x["dec"] for x in w)
        region = ConeRegion(mid["ra"], mid["dec"], 0.5 * max(ra_span, de_span, 0.1) + 0.05)
        rng = MjdRange(w[0]["mjd"] - 0.5, w[-1]["mjd"] + 1.5)
        try:
            for obs in adapter.discover(region, rng, store):
                if obs.band not in ("W1", "W2", "W3", "W4"):
                    continue
                # interpolate the daily ephemeris to the frame time
                m = obs.t_mid_mjd_utc
                mj = np.array([x["mjd"] for x in w])
                if m < mj.min() - 1 or m > mj.max() + 1:
                    continue
                ra = np.interp(m, mj, [x["ra"] for x in w])
                de = np.interp(m, mj, [x["dec"] for x in w])
                if adapter.nominal_footprint(obs, pad_arcsec=-60.0).contains(ra, de):
                    frames.append(obs)
        except Exception as exc:
            print(f"  discover fail window {w[0]['date']}: {exc}", flush=True)
    frames.sort(key=lambda o: o.t_mid_mjd_utc)
    print(f"{len(frames)} frames contain the predicted position", flush=True)
    if len(frames) > a.max_frames:
        idx = np.linspace(0, len(frames) - 1, a.max_frames).astype(int)
        frames = [frames[i] for i in idx]

    # 3. exact ephemeris at frame mid-times
    eph = {}
    for i in range(0, len(frames), 80):
        batch = frames[i:i + 80]
        tlist = " ".join(f"{o.t_mid_mjd_utc + JD_MJD:.6f}" for o in batch)
        rows, raw, url = horizons(session, a.asteroid, TLIST=f"'{tlist}'", TLIST_TYPE="JD")
        store.store(service_url=url, query=f"tlist batch {i}",
                    request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=raw, row_count=len(rows), http_status=200)
        for o, r in zip(batch, rows):
            eph[o.observation_id] = r

    # 4. cutouts, flux maps, samples (real + ring + designated)
    ring = np.array(C.RING_OFFSETS)
    offsets = np.vstack([[0.0, 0.0], ring])          # (49, 2)
    nr, nt = len(RESID_GRID), len(offsets)
    recs, F, V, G, meta = [], [], [], [], []
    t0 = _time.monotonic()
    for k, o in enumerate(frames):
        e = eph.get(o.observation_id)
        if e is None:
            continue
        cut = CutoutSpec(e["ra"], e["dec"], CUTOUT_PIX)
        try:
            ps = adapter.fetch(o, ["int", "unc"], cut_dir, cutout=cut)
            msk_dir = C.MSK_DIR
            msk_name = o.products["msk"]["url"].rsplit("/", 1)[-1]
            msk_path = msk_dir / msk_name
            if not msk_path.exists():
                adapter.fetch(o, ["msk"], cut_dir)
                msk_path = cut_dir / msk_name
        except Exception as exc:
            recs.append({"observation_id": o.observation_id, "status": f"fetch fail: {exc}"})
            continue
        int_p = next(p.path for p in ps.products if p.kind == "int")
        unc_p = next(p.path for p in ps.products if p.kind == "unc")
        try:
            fm = build_flux_map(int_p, unc_p, msk_path, o.band, o.t_mid_mjd_utc,
                                WiseExactFootprint.FATAL_MASK,
                                magzp=o.quality_flags.get("magzp"), keep_inputs=True,
                                pix_scale_from_header=True)
        except Exception as exc:
            recs.append({"observation_id": o.observation_id, "status": f"fluxmap fail: {exc}"})
            continue
        if fm.magzp is None:
            continue
        f0, v0, g0 = fm.sample(e["ra"], e["dec"])
        if not np.isfinite(f0[0]):
            recs.append({"observation_id": o.observation_id, "status": "outside usable cutout"})
            continue
        scale = 10.0 ** ((C.ZP_REF - fm.magzp) / 2.5)
        cosd = np.cos(np.deg2rad(e["dec"]))
        Fk = np.empty((nt, nr, nr), np.float32)
        Vk = np.empty_like(Fk)
        Gk = np.empty_like(Fk)
        for ti, (dx, dy) in enumerate(offsets):
            ra_q = e["ra"] + (RESID_GRID[:, None] + dx) / 3600.0 / cosd + 0 * RESID_GRID[None, :]
            dec_q = e["dec"] + (RESID_GRID[None, :] + dy) / 3600.0 + 0 * RESID_GRID[:, None]
            f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
            Fk[ti] = f.reshape(nr, nr) * scale
            Vk[ti] = np.minimum(v.reshape(nr, nr) * scale * scale, 1e30)
            Gk[ti] = g.reshape(nr, nr)
        F.append(Fk); V.append(Vk); G.append(Gk)
        meta.append((o.t_mid_mjd_utc, C.BAND_IDX[o.band], float(e["v"]), float(fm.magzp),
                     float(f0[0]), float(v0[0])))
        recs.append({"observation_id": o.observation_id, "status": "ok", "band": o.band,
                     "mjd": o.t_mid_mjd_utc, "v_pred": e["v"], "ra": e["ra"], "dec": e["dec"],
                     "single_snr": float(f0[0] / np.sqrt(v0[0])),
                     "single_mag": (float(fm.magzp - 2.5 * np.log10(f0[0])) if f0[0] > 0 else None),
                     "magzp": float(fm.magzp), "native_key": o.native_key,
                     "quality": o.quality_flags})
        if (k + 1) % 20 == 0:
            print(f"  {k + 1}/{len(frames)} ({(k + 1) / (_time.monotonic() - t0):.2f}/s)", flush=True)
    with (run_dir / "frames.jsonl").open("w") as fh:
        for r in recs:
            fh.write(json.dumps(r, default=float) + "\n")
    if not F:
        print("no usable frames"); return
    F, V, G = np.array(F), np.array(V), np.array(G)
    meta = np.array(meta)
    np.savez_compressed(run_dir / "samples.npz", f=F, v=V, g=G, meta=meta,
                        meta_fields=np.array(["mjd", "band", "v_pred", "magzp", "f0", "v0"]),
                        resid_grid=RESID_GRID, offsets=offsets)

    analyse(a, run_dir, session, store, recs, F, V, G, meta)


def analyse(a, run_dir, session, store, recs, F, V, G, meta):
    # 5. NEOWISE pipeline photometry of the asteroid on the same frames
    pipe = {}
    ok_recs = [r for r in recs if r["status"] == "ok"]
    # group the ok frames into visit windows (gaps > 5 d) and query one
    # cone per window and table (IRSA ADQL rejects OR-ed CONTAINS terms)
    ok_sorted = sorted(ok_recs, key=lambda r: r["mjd"])
    groups = []
    for r in ok_sorted:
        if groups and r["mjd"] - groups[-1][-1]["mjd"] <= 5.0:
            groups[-1].append(r)
        else:
            groups.append([r])
    want = {(r["native_key"]["scan_id"], r["native_key"]["frame_num"]): r for r in ok_recs}
    for g in groups:
        ra = np.array([r["ra"] for r in g]); dec = np.array([r["dec"] for r in g])
        cosd = np.cos(np.deg2rad(dec.mean()))
        rad = max(np.hypot((ra - ra.mean()) * cosd, dec - dec.mean()).max() + 0.01, 0.02)
        table = "allsky_4band_p1bs_psd" if g[0]["mjd"] < 55600 else "neowiser_p1bs_psd"
        cols = "scan_id, frame_num, ra, dec, w1mpro, w1sigmpro, w2mpro, w2sigmpro, mjd"
        if table.startswith("allsky"):
            cols += ", w3mpro, w3sigmpro, w4mpro, w4sigmpro"
        q = (f"SELECT {cols} FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),"
             f"CIRCLE('ICRS',{ra.mean():.6f},{dec.mean():.6f},{rad:.5f}))=1 "
             f"AND mjd >= {g[0]['mjd'] - 0.01:.5f} AND mjd <= {g[-1]['mjd'] + 0.01:.5f}")
        rows, raw = tap_rows(session, q)
        store.store(service_url=TAP, query=q, request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=raw, row_count=len(rows), http_status=200)
        if rows and "scan_id" not in rows[0]:
            print(f"  {table}: unexpected response: {raw[:200]!r}", flush=True)
            continue
        for row in rows:
            key = (row["scan_id"], int(row["frame_num"]))
            r = want.get(key)
            if r is None:
                continue
            sep = np.hypot((float(row["ra"]) - r["ra"]) * cosd, float(row["dec"]) - r["dec"]) * 3600
            if sep < 3.0 and (key not in pipe or sep < pipe[key]["_sep"]):
                row["_sep"] = float(sep)
                pipe[key] = row
    (run_dir / "pipeline_photometry.json").write_text(json.dumps(
        {f"{k[0]}/{k[1]}": v for k, v in pipe.items()}, indent=1))

    # 6. stack exactly as the survey does; R vs the ring null
    summary = {"asteroid": a.asteroid, "n_frames_ok": int(len(meta)), "n_pipeline_matches": len(pipe)}
    designated = np.array([1 + i for i in C.DESIGNATED])
    for b in sorted(set(meta[:, 1].astype(int).tolist())):
        sel = meta[:, 1].astype(int) == b
        f, v, g = F[sel].transpose(1, 0, 2, 3), V[sel].transpose(1, 0, 2, 3), G[sel].transpose(1, 0, 2, 3)
        mx = nulls.cell_maxima(f, v, g)
        R, T = nulls.exceedance_ratios(mx, designated)
        S0, A0, B0, n0, _ = nulls.stack_S(f[0], v[0], g[0], return_parts=True)
        s_max, node = nulls.grid_max(S0)
        flux = A0[node] / B0[node]
        mag = C.ZP_REF - 2.5 * np.log10(flux) if flux > 0 else None
        ring_R = R[1:]
        rs = nulls.rank_statement(R[0], ring_R)
        q95 = float(np.nanquantile(ring_R, C.NORM_QUANTILE))
        # throughput vs pipeline photometry, frame by frame
        band = C.BAND_NAME[b]
        diffs = []
        for r in ok_recs:
            if r["band"] != band or r["single_mag"] is None:
                continue
            p = pipe.get((r["native_key"]["scan_id"], r["native_key"]["frame_num"]))
            if p and p.get(f"w{b}mpro") not in (None, "", "null"):
                diffs.append(float(r["single_mag"]) - float(p[f"w{b}mpro"]))
        summary[band] = {
            "n_epochs": int(sel.sum()), "S_max_real": float(s_max), "S_centre": float(S0[2, 2]),
            "T_designated": float(T), "R_real": float(R[0]),
            "ring_R_median": float(np.nanmedian(ring_R)), "ring_R_q95": q95,
            "R_norm_real": float(R[0] / q95) if q95 > 0 else None,
            "rank_statement": rs,
            "peak_offset_arcsec": [float(RESID_GRID[node[0]]), float(RESID_GRID[node[1]])],
            "recovered_stack_mag": mag,
            "predicted_V_median": float(np.median(meta[sel, 2])),
            "single_frame_snr_median": float(np.median(meta[sel, 4] / np.sqrt(meta[sel, 5]))),
            "single_frames_above_5sigma": int((meta[sel, 4] / np.sqrt(meta[sel, 5]) > 5).sum()),
            "throughput_vs_pipeline": {
                "n": len(diffs), "median_dmag": float(np.median(diffs)) if diffs else None,
                "mad_dmag": float(1.4826 * np.median(np.abs(np.array(diffs) - np.median(diffs)))) if diffs else None,
                "note": "forced matched-filter mag (Gaussian kernel) minus pipeline w?mpro on the same frame; positive = our flux lower",
            },
        }
        print(f"[{band}] N={sel.sum()} S_max={s_max:.1f} (centre {S0[2, 2]:.1f}) T={T:.2f} "
              f"R={R[0]:.2f} ring q95 {q95:.2f} -> R~={summary[band]['R_norm_real']} "
              f"rank p={rs['p_cell']:.3f}; stack mag {mag if mag is None else round(mag, 2)} "
              f"vs V {np.median(meta[sel, 2]):.2f}; throughput vs pipeline "
              f"{summary[band]['throughput_vs_pipeline']['median_dmag']} ({len(diffs)} frames)", flush=True)
    (run_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
