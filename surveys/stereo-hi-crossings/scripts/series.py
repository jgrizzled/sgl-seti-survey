"""STEREO-A HI-1 search driver (threshold freeze v1.0).

Stages (role = dev | confirmatory):
  python series.py index <role>    # frame -> task index + patch table
  python series.py run <role>      # fetch (RAL L1) + measure per day, purge
  python series.py reduce <role>   # series -> statistics -> results json

Index: for every covered event of the role's units, every listed frame
within +-13 d of t_ca whose day header centre puts the source inside
the footprint, classified 'arc' (in-beam, 4.05 <= |HPLN| <= 6.0) or
'base' (6 < |HPLN| <= 12). Patches: source (fixed ICRS) + 8 controls
offset by HPLT +-1..+-4 deg at the arc midpoint (fixed ICRS). Tycho-2
VT <= 9 within 3 px drops a patch (S2 source exempt).

Measure: one pass per frame — gate, load, verify WCS, 2-D ZP, top-hat
photometry at every mapped patch (flux normalized to ZP_REF). Every
36th frame dumps calibrators (differential flat, era stability) and a
stamp-response test (5 Gaussian stamps, FWHM 1.7 px, at random empty
positions through the identical chain).

Reduce: per (unit, event, patch): E = F_arc - g*<F_base>; same-frame
differential D_src = E_src - median_k E_k (controls: minus the median of
the other seven); z per event = mean(D_arc)/[MAD-SD(D_arc)/sqrt(n)];
S_stack, S_event, S_pulse; T = max over the 8 control patches.
Planet / Earth / Moon proximity and bright-body in-FOV vetoes and the
frozen recon exclusions are applied here.
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hi_geometry as G
import hi_lib as H

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "stereo-hi-crossings"
CFG_FILE = os.environ.get("HI_FREEZE", "threshold_freeze_v1_2.json")   # v1.1: L2 substrate, 24 px margin; v1.2: 2-deg bright-body proximity veto
CFG = json.loads((SURV / "configs" / CFG_FILE).read_text())
RUNS = REPO / "runs" / "stereo-hi-crossings"
LIST = RUNS / "coverage" / "listings"
HEAD = RUNS / "coverage" / "day_headers"
SER = RUNS / "series"
EVENTS = REPO / "crossings" / "stereoa_v1" / "events.ecsv"
STEM = re.compile(r'href="((\d{8})_(\d{6})_24h1A_br01\.fts)"')

ARC = CFG["geometry"]["arc_abs_hpln_deg"]
BASE = CFG["geometry"]["baseline_abs_hpln_deg"]
R_AU = CFG["geometry"]["in_beam_radius_au"]
OFFS = CFG["controls"]["hplt_offsets_deg"]
HALF_W = G.HI1_HALF_W
MARGIN_PX = CFG["geometry"].get("footprint_margin_px", 8)
MARGIN = MARGIN_PX * 0.01998
PREFER = "ssc" if "L2" in CFG["substrate"]["primary"] else "ral"
DUMP_EVERY = 36
SEED = CFG["seed"]


# ----------------------------------------------------------------- index

def load_listing(day: str):
    p = LIST / f"ssc_{day}.html"
    if not p.exists():
        return []
    out = []
    for fname, d, t in STEM.findall(p.read_text(errors="replace")):
        iso = f"{d[:4]}-{d[4:6]}-{d[6:]}T{t[:2]}:{t[2:4]}:{t[4:]}"
        out.append((f"{d}_{t}", Time(iso, scale="utc").mjd))
    return out


def day_centre(day: str):
    p = HEAD / f"{day}.json"
    if not p.exists():
        return None
    h = json.loads(p.read_text())
    try:
        return (float(h["CRVAL1"]), float(h["CRVAL2"]))
    except (TypeError, ValueError, KeyError):
        return None


def build_index(role: str) -> None:
    t = Table.read(EVENTS)
    by_id = {str(r["event_id"]): r for r in t}
    tyc = np.load(H.TYCHO)
    from scipy.spatial import cKDTree
    sel = tyc["vt"] <= CFG["gates"]["patch"]["tycho2_vt_max"]
    tree = cKDTree(G.radec_to_vec(tyc["ra"][sel], tyc["dec"][sel]))
    lim = 2 * np.sin(np.radians(CFG["gates"]["patch"]["tycho2_proximity_px"] * 0.01998) / 2)

    index: dict[str, list] = defaultdict(list)
    patches, evmeta = {}, {}
    days_needed = set()
    for u in CFG["units"]:
        if u["role"] != role:
            continue
        ukey = f'{u["target_id"]}|{u["unit"]}'
        ch = u["channel"]
        for eid in u["event_ids_covered"]:
            ev = by_id[eid]
            e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
                 "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
            mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
            s_hat = G.source_direction(e, ch)
            # frames within +-13 d
            d0, d1 = int(mjd_ca - 13), int(mjd_ca + 13)
            frames = []
            for d in range(d0, d1 + 1):
                day = Time(d, format="mjd").strftime("%Y%m%d")
                c = day_centre(day)
                if c is None:
                    continue
                for stem, m in load_listing(day):
                    frames.append((stem, m, c))
            if not frames:
                continue
            mj = np.array([f[1] for f in frames])
            hpln, hplt = G.helioprojective(s_hat, mj)
            b_e = G.impact_parameter_au(e, mj)
            kinds = []
            for k, (stem, m, c) in enumerate(frames):
                inside = (abs(hpln[k] - c[0]) <= HALF_W[0] - MARGIN
                          and abs(hplt[k] - c[1]) <= HALF_W[1] - MARGIN)
                if not inside:
                    kinds.append(None)
                elif b_e[k] <= R_AU and ARC[0] <= abs(hpln[k]) <= ARC[1]:
                    kinds.append("arc")
                elif BASE[0] < abs(hpln[k]) <= BASE[1]:
                    kinds.append("base")
                else:
                    kinds.append(None)
            arc_m = [m for (s_, m, c), kd in zip(frames, kinds) if kd == "arc"]
            if len(arc_m) < CFG["gates"]["event"]["min_arc_epochs"]:
                continue
            t_mid = float(np.median(arc_m))
            ln0, lt0 = G.helioprojective(s_hat, t_mid)
            pl = [tuple(map(float, G.vec_to_radec(s_hat)))]
            for off in OFFS:
                v = G.from_helioprojective(ln0[0], lt0[0] + off, t_mid)
                pl.append(tuple(map(float, G.vec_to_radec(v))))
            # Tycho-2 patch mask (S2 source exempt)
            masked = []
            for pi, (ra, dec) in enumerate(pl):
                if pi == 0 and ch == "S2":
                    masked.append(False)
                    continue
                dd, _ = tree.query(G.radec_to_vec(ra, dec), k=1)
                masked.append(bool(dd < lim))
            key = f"{ukey}|{eid}"
            patches[key] = {"radec": pl, "masked": masked, "t_mid": t_mid,
                            "hpln0": float(ln0[0]), "hplt0": float(lt0[0])}
            evmeta[key] = {"mjd_ca": mjd_ca, "side": "pre" if t_mid < mjd_ca else "post",
                           "n_arc_listed": len(arc_m),
                           "n_base_listed": sum(1 for kd in kinds if kd == "base")}
            for (stem, m, c), kd in zip(frames, kinds):
                if kd:
                    index[stem].append({"u": ukey, "e": eid, "k": kd})
                    days_needed.add(stem[:8])
    SER.mkdir(parents=True, exist_ok=True)
    (SER / f"index_{role}_v1.json").write_text(json.dumps(
        {"index": index, "patches": patches, "events": evmeta}, indent=0) + "\n")
    n_arc = sum(1 for v in index.values() if any(t["k"] == "arc" for t in v))
    n_src_masked = sum(1 for p in patches.values() if p["masked"][0])
    print(f"{role}: {len(index)} frames ({n_arc} arc frames) over {len(days_needed)} days; "
          f"{len(patches)} unit-events; source patches Tycho-masked: {n_src_masked}")


# ----------------------------------------------------------------- measure

def _stamp_response(fr: H.Frame, rng: np.random.Generator, n: int = 5) -> list:
    """Inject Gaussian stamps (FWHM 1.7 px) of known ZP-normalized flux
    at random positions through the identical chain (high-pass on the
    raw image + top-hat aperture); returns recovered/injected ratios."""
    out = []
    sig = 1.73 / 2.355
    yy, xx = np.mgrid[-6:7, -6:7]
    for _ in range(n):
        x, y = rng.uniform(40, 980, 2)
        zp = fr.zp_at(x, y)
        if not np.isfinite(zp):
            continue
        f_norm = 10 ** (rng.uniform(-0.6, 0.6))          # V_eq 9.6-12.6
        f_dn = f_norm * 10 ** (0.4 * (zp - H.ZP_REF))
        ix, iy = int(round(x)), int(round(y))
        stamp = f_dn * np.exp(-0.5 * ((xx - (x - ix)) ** 2 + (yy - (y - iy)) ** 2) / sig ** 2)
        stamp /= stamp.sum()
        stamp *= f_dn
        raw = fr.raw.copy()
        raw[iy - 6:iy + 7, ix - 6:ix + 7] += stamp
        # local high-pass recompute on a window (identical filter)
        from scipy import ndimage
        w = 40
        sub = raw[iy - w:iy + w + 1, ix - w:ix + w + 1]
        fill = np.nan_to_num(sub, nan=float(np.nanmedian(sub)))
        hp = sub - ndimage.median_filter(fill, size=H.HP_SIZE)
        f0, _, _ = H.aper_flux(fr.img, x, y)
        f1, _, _ = H.aper_flux(hp, w + (x - ix), w + (y - iy))
        if np.isfinite(f0) and np.isfinite(f1):
            out.append({"x": round(x, 1), "y": round(y, 1), "f_norm": f_norm,
                        "ratio": float((f1 - f0) / f_dn)})
    return out


def measure_frame(args):
    path, tasks, patches, dump = args
    out = {"frame": Path(path).name, "records": []}
    try:
        from astropy.io import fits
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            h = fits.getheader(path)
        reason = H.frame_gate(h)
        if reason:
            out["unusable"] = reason
            return out
        fr = H.load_frame(Path(path))
        ok = H.fit_frame(fr)
        out["fit"] = {"ok": bool(ok), "n_match": fr.n_match, "rms": round(fr.astrom_rms, 3),
                      "zp0": round(float(fr.zp_coef[0]), 4) if fr.zp_coef is not None else None,
                      "zp_mad": round(fr.zp_scatter, 4) if np.isfinite(fr.zp_scatter) else None,
                      "n_calib": fr.n_calib, "mjd": fr.mjd_avg,
                      "centre": [fr.hpln_c, fr.hplt_c], "level": fr.level}
        if not ok or fr.zp_coef is None or fr.zp_scatter > CFG["chain"]["zp"]["max_mad_mag"]:
            out["unusable"] = "astrometry_or_zp"
            return out
        for tk in tasks:
            key = f'{tk["u"]}|{tk["e"]}'
            pt = patches[key]
            for pi, (ra, dec) in enumerate(pt["radec"]):
                if pt["masked"][pi]:
                    continue
                m = H.measure_patch(fr, ra, dec)
                if not H.in_footprint(m["x"], m["y"], MARGIN_PX):
                    continue
                m.update({"u": tk["u"], "e": tk["e"], "k": tk["k"], "p": pi})
                out["records"].append(m)
        if dump:
            import warnings
            from astropy.wcs import WCS
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                wp = WCS(h)
            cal = np.array(fr.calibrators) if fr.calibrators else np.zeros((0, 6))
            if len(cal):
                ln, lt = wp.all_pix2world(cal[:, 0], cal[:, 1], 0)
                zp_at = np.array([fr.zp_at(x, y) for x, y in cal[:, :2]])
                fn = cal[:, 5] * 10 ** (-0.4 * (zp_at - H.ZP_REF))
                out["calibrators"] = [[round(float(a), 1), round(float(b), 1), round(float(c), 3),
                                       round(float(d), 3), round(float(v), 2), round(float(bv), 2),
                                       round(float(f), 4)]
                                      for a, b, c, d, v, bv, f in
                                      zip(cal[:, 0], cal[:, 1], ln, lt, cal[:, 3], cal[:, 4], fn)]
            rng = np.random.default_rng(SEED + int(fr.mjd_avg * 100) % 100000)
            out["stamps"] = _stamp_response(fr, rng)
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
    return out


def run(role: str, workers: int = 20, fetchers: int = 6) -> None:
    idx = json.loads((SER / f"index_{role}_v1.json").read_text())
    patches = idx["patches"]
    out_path = SER / f"measurements_{role}_v1.jsonl"
    done = set()
    if out_path.exists():
        with open(out_path) as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["frame"][:15])
                except Exception:
                    pass
    stems = sorted(s for s in idx["index"] if s not in done)
    by_day = defaultdict(list)
    for s in stems:
        by_day[s[:8]].append(s)
    print(f"{role}: {len(stems)} frames to do over {len(by_day)} days ({len(done)} done)", flush=True)
    fetch_log = open(SER / f"fetch_{role}_v1.log", "a")
    t0 = time.time()
    n_done = 0
    with open(out_path, "a") as fh, ProcessPoolExecutor(workers) as ex, \
            ThreadPoolExecutor(fetchers) as tp:
        days = sorted(by_day)
        # pipeline: fetch day i+1 while measuring day i
        def fetch_day(day):
            res = list(tp.map(lambda s: (s, *H.fetch_frame(day, s, prefer=PREFER, fallback=False)), by_day[day]))
            for s, p, src in res:
                fetch_log.write(f"{s} {src} {p.name if p else '-'}\n")
            fetch_log.flush()
            return res
        pending = None
        for i, day in enumerate(days):
            res = pending if pending is not None else fetch_day(day)
            nxt = None
            if i + 1 < len(days):
                import threading
                holder = {}
                th = threading.Thread(target=lambda: holder.update(r=fetch_day(days[i + 1])))
                th.start()
                nxt = (th, holder)
            args = []
            for j, (s, p, src) in enumerate(res):
                if p is None:
                    fh.write(json.dumps({"frame": s + "_missing", "unusable": "fetch_failed"}) + "\n")
                    continue
                dump = (int(s[9:11]) == 0 and int(s[11:13]) < 40)   # first frame of the day
                args.append((str(p), idx["index"][s], patches, dump))
            for rec in ex.map(measure_frame, args, chunksize=1):
                fh.write(json.dumps(rec) + "\n")
            fh.flush()
            n_done += len(res)
            for s, p, src in res:
                if p is not None:
                    try:
                        p.unlink()
                    except OSError:
                        pass
            try:
                (H.FRAMES / day).rmdir()
            except OSError:
                pass
            if i % 5 == 0:
                el = time.time() - t0
                print(f"  {day}: {n_done}/{len(stems)} frames, {el/60:.1f} min, "
                      f"{n_done/max(el,1)*3600:.0f} frames/h", flush=True)
            if nxt:
                nxt[0].join()
                pending = nxt[1].get("r")
            else:
                pending = None
    fetch_log.close()
    print("done", flush=True)


# ----------------------------------------------------------------- reduce

def _robust_sd(x):
    x = np.asarray(x, float)
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else np.nan


def reduce(role: str, version: str = "v1") -> dict:
    idx = json.loads((SER / f"index_{role}_v1.json").read_text())
    patches, evmeta = idx["patches"], idx["events"]
    excl = CFG["recon_exclusions"]
    excl_mjd = {e: [Time(t).mjd for t in ts] for e, ts in excl.items()}
    gp = CFG["gates"]["epoch"]

    # gather: key (u,e,k,mjd) -> {p: flux}
    epochs = defaultdict(dict)
    frames_ok, frames_bad = 0, defaultdict(int)
    fit_rows, cal_rows, stamp_rows = [], [], []
    with open(SER / f"measurements_{role}_v1.jsonl") as fh:
        for line in fh:
            rec = json.loads(line)
            if "unusable" in rec or "error" in rec:
                frames_bad[rec.get("unusable", "error")] += 1
                continue
            frames_ok += 1
            if "fit" in rec:
                fit_rows.append(rec["fit"])
            cal_rows.extend(rec.get("calibrators", []))
            stamp_rows.extend(rec.get("stamps", []))
            for m in rec["records"]:
                if m["valid"] < gp["aperture_valid_fraction"] or not np.isfinite(m["flux"] or np.nan):
                    continue
                epochs[(m["u"], m["e"], m["k"], round(m["mjd"], 5))][m["p"]] = (m["flux"], m["err"])

    # differential flat g(|HPLN| band, HPLT): ratio of a static star's flux in
    # the arc band to its flux in the baseline band, from the calibrator dump
    g_map = _differential_flat(cal_rows)

    # planet / bright-body vetoes per epoch
    mjds = np.array(sorted({k[3] for k in epochs}))
    pdirs = G.planet_directions(mjds) if len(mjds) else {}
    lim_prox = 2 * np.sin(np.radians(gp["planet_proximity_px"] * 0.01998) / 2)
    n_veto = defaultdict(int)
    keep = {}
    for key, by_p in epochs.items():
        u, e, k, mjd = key
        if e in excl_mjd and any(abs(mjd - x) < 0.02 for x in excl_mjd[e]) and k == "arc":
            n_veto["recon_exclusion"] += 1
            continue
        i = int(np.searchsorted(mjds, mjd))
        pt = patches[f"{u}|{e}"]
        bad = None
        pv = np.array([G.radec_to_vec(ra, dec) for ra, dec in pt["radec"]])
        lim_bright = 2 * np.sin(np.radians(gp.get("bright_body_proximity_deg", 2.0)) / 2)
        for body, dirs in pdirs.items():
            d = dirs[i]
            dist = np.linalg.norm(pv - d, axis=1)
            if np.any(dist < lim_prox):
                bad = f"prox_{body}"
                break
            if body in gp.get("bright_body_list", []) and np.any(dist < lim_bright):
                bad = f"near2deg_{body}"
                break
        if bad:
            n_veto[bad] += 1
            continue
        keep[key] = by_p

    # group usable epochs by (u, e, kind)
    grouped = defaultdict(dict)
    for (u, e, k, mjd), by_p in keep.items():
        grouped[(u, e, k)][mjd] = by_p

    # per (u, e): baseline means per patch, then arc differentials
    results = {}
    series_dump = {}
    units = sorted({k.split("|")[0] + "|" + k.split("|")[1] for k in patches})
    ctrl_idx = list(range(1, 9))
    min_ctrl = CFG["controls"]["min_controls_per_epoch"]
    for ukey in units:
        per_patch = {}   # probe -> {"z": [...], "pulse": [...], "ev": [...]}
        ev_detail = {}
        for key, pt in patches.items():
            if not key.startswith(ukey + "|"):
                continue
            eid = key.split("|")[2]
            arc = grouped.get((ukey, eid, "arc"), {})
            base = grouped.get((ukey, eid, "base"), {})
            if len(arc) < CFG["gates"]["event"]["min_arc_epochs"] or \
                    len(base) < CFG["gates"]["event"]["min_baseline_epochs"]:
                ev_detail[eid] = {"status": "insufficient_epochs", "n_arc": len(arc), "n_base": len(base)}
                continue
            # baseline mean per patch (robust: median), g-scaled
            B = {}
            for pi in range(9):
                vals = [v[pi][0] for v in base.values() if pi in v]
                if len(vals) >= CFG["gates"]["event"]["min_baseline_epochs"]:
                    B[pi] = float(np.median(vals))
            g = g_map.get("ratio", 1.0)
            # per-epoch E and same-frame differential
            D = {pi: [] for pi in range(9)}
            mj_used = []
            for m in sorted(arc):
                v = arc[m]
                E = {pi: v[pi][0] - g * B[pi] for pi in range(9) if pi in v and pi in B}
                cs = [pi for pi in ctrl_idx if pi in E]
                if 0 not in E or len(cs) < min_ctrl:
                    continue
                mj_used.append(m)
                D[0].append(E[0] - float(np.median([E[c] for c in cs])))
                for c in cs:
                    others = [E[o] for o in cs if o != c]
                    D[c].append(E[c] - float(np.median(others)))
            if len(D[0]) < CFG["gates"]["event"]["min_arc_epochs"]:
                ev_detail[eid] = {"status": "insufficient_after_veto", "n_arc": len(D[0])}
                continue
            series_dump[key] = {"mjd": [round(m, 5) for m in mj_used],
                                "D": {str(pi): [round(x, 5) for x in D[pi]] for pi in range(9) if len(D[pi]) == len(mj_used)}}
            for pi in range(9):
                d = np.array(D[pi])
                if len(d) < CFG["gates"]["event"]["min_arc_epochs"]:
                    continue
                sd = _robust_sd(d)
                if not sd > 0:
                    continue
                z = float(np.mean(d) / (sd / np.sqrt(len(d))))
                pulse = float(np.max(d - np.median(d)) / sd)
                pp = per_patch.setdefault(pi, {"z": [], "pulse": [], "ev": []})
                pp["z"].append(z)
                pp["pulse"].append(pulse)
                pp["ev"].append(eid)
            d0 = np.array(D[0])
            ev_detail[eid] = {"status": "included", "n_arc": len(d0), "n_base": len(base),
                              "z": round(float(np.mean(d0) / (_robust_sd(d0) / np.sqrt(len(d0)))), 2),
                              "mean_D": round(float(np.mean(d0)), 4), "sd_D": round(_robust_sd(d0), 4),
                              "side": evmeta[key]["side"], "mjd_ca": round(evmeta[key]["mjd_ca"], 2),
                              "pulse_mjd": round(float(mj_used[int(np.argmax(d0))]), 4),
                              "B_src": round(B.get(0, np.nan), 4)}
        stats = {}
        for pi, pp in per_patch.items():
            if len(pp["z"]) < CFG["gates"]["unit"]["min_events_for_stack"]:
                continue
            z = np.array(pp["z"])
            stats[pi] = {"n_events": int(len(z)), "S_stack": float(np.sum(z) / np.sqrt(len(z))),
                         "S_event": float(np.max(z)), "S_pulse": float(np.max(pp["pulse"])),
                         "z_median": float(np.median(z)), "z_mad": _robust_sd(z)}
        if 0 not in stats:
            results[ukey] = {"status": "constraint_only", "events": ev_detail,
                             "n_included": sum(1 for v in ev_detail.values() if v["status"] == "included")}
            continue
        res = {"status": "searched", "per_patch": {str(k): v for k, v in stats.items()},
               "events": ev_detail}
        src = stats[0]
        ctrl = [v for k, v in stats.items() if k != 0]
        for s in ("S_stack", "S_event", "S_pulse"):
            vals = [c[s] for c in ctrl]
            T = max(vals) if vals else 0.0
            mad = _robust_sd(vals) if len(vals) >= 3 else None
            res[s] = {"S": round(src[s], 3), "T": round(float(T), 3), "n_controls": len(vals),
                      "exceeds": bool(src[s] > max(T, 0.0)),
                      "S_over_ctrl_mad": round(src[s] / mad, 2) if mad else None}
        results[ukey] = res

    n_trials = sum(1 for r in results.values() if r["status"] == "searched" for s in ("S_stack", "S_event", "S_pulse"))
    n_exc = sum(1 for r in results.values() if r["status"] == "searched"
                for s in ("S_stack", "S_event", "S_pulse") if r[s]["exceeds"])
    fit_ok = [f for f in fit_rows if f.get("ok")]
    out = {"role": role, "construction": "threshold_freeze_v1.0",
           "frames": {"usable": frames_ok, "unusable": dict(frames_bad),
                      "median_n_match": float(np.median([f["n_match"] for f in fit_ok])) if fit_ok else None,
                      "median_astrom_rms_px": float(np.median([f["rms"] for f in fit_ok])) if fit_ok else None,
                      "median_zp0": float(np.median([f["zp0"] for f in fit_ok if f["zp0"]])) if fit_ok else None,
                      "median_zp_mad": float(np.median([f["zp_mad"] for f in fit_ok if f["zp_mad"]])) if fit_ok else None},
           "epoch_vetoes": dict(n_veto), "differential_flat": g_map,
           "stamp_response": {"n": len(stamp_rows),
                              "median_ratio": float(np.median([s["ratio"] for s in stamp_rows])) if stamp_rows else None,
                              "mad": _robust_sd([s["ratio"] for s in stamp_rows]) if stamp_rows else None},
           "trials": {"searched": n_trials, "exceedances": n_exc,
                      "expected_control_crossings": round(n_trials / 9, 2)},
           "units": results}
    (SURV / "results" / f"{role}_search_{version}.json").write_text(json.dumps(out, indent=1) + "\n")
    (SER / f"series_{role}_{version}.json").write_text(json.dumps(series_dump) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "units"}, indent=1))
    for u, r in results.items():
        if r["status"] != "searched":
            print(u, "->", r["status"], r.get("n_included"))
            continue
        line = f'{u}: n_ev={r["per_patch"]["0"]["n_events"]} zmed={r["per_patch"]["0"]["z_median"]:+.2f} zmad={r["per_patch"]["0"]["z_mad"]:.2f}'
        for s in ("S_stack", "S_event", "S_pulse"):
            line += f' | {s} S={r[s]["S"]:.2f} T={r[s]["T"]:.2f}{" EXC" if r[s]["exceeds"] else ""}'
        print(line)
    return out


def _differential_flat(cal_rows) -> dict:
    """g = median over stars of flux(arc band)/flux(base band), stars
    identified by rounded catalogue position is unavailable here, so the
    ensemble ratio is formed on the population: median normalized flux
    residual (log f_norm + 0.4 V_inst) per |HPLN| band."""
    if not cal_rows:
        return {"ratio": 1.0, "n": 0, "note": "no calibrator dump"}
    c = np.array(cal_rows, float)   # x, y, hpln, hplt, V, bv, f_norm
    ok = c[:, 6] > 0
    c = c[ok]
    v_inst = c[:, 4] - H.COLOUR_COEFF * (c[:, 5] - H.BV_REF)
    resid = 2.5 * np.log10(c[:, 6]) + v_inst - H.ZP_REF     # 0 if the ZP model is exact
    a = np.abs(c[:, 2])
    arc = (a >= ARC[0]) & (a <= ARC[1])
    base = (a > BASE[0]) & (a <= BASE[1])
    r_arc = float(np.median(resid[arc])) if arc.sum() > 20 else 0.0
    r_base = float(np.median(resid[base])) if base.sum() > 20 else 0.0
    prof = {}
    for lo in np.arange(4, 24, 2.0):
        m = (a >= lo) & (a < lo + 2)
        if m.sum() > 20:
            prof[f"{lo:.0f}-{lo+2:.0f}"] = {"n": int(m.sum()), "resid_mag": round(float(np.median(resid[m])), 4),
                                            "mad": round(_robust_sd(resid[m]), 4)}
    return {"ratio": float(10 ** (0.4 * (r_arc - r_base))), "n": int(len(c)),
            "arc_resid_mag": r_arc, "base_resid_mag": r_base, "profile": prof,
            "note": "ratio = arc-band/base-band ensemble flux residual of static stars after the 2-D ZP"}


if __name__ == "__main__":
    cmd, role = sys.argv[1], sys.argv[2]
    if cmd == "index":
        build_index(role)
    elif cmd == "run":
        run(role, workers=int(sys.argv[3]) if len(sys.argv) > 3 else 20)
    elif cmd == "reduce":
        reduce(role, sys.argv[3] if len(sys.argv) > 3 else "v1")
