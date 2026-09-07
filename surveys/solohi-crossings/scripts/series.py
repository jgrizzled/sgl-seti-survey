"""Solar Orbiter / SoloHI search driver (threshold freeze v1.x).

Stages (role = dev | confirmatory):
  python series.py index <role>    # frame -> task index + patch table (Gaia templates)
  python series.py run <role>      # fetch (NRL + SOAR) + measure per day, purge
  python series.py reduce <role>   # series -> statistics -> results json

Index: for every covered event of the role's units, every inner-tile L2
frame (SOAR inventory snapshot) of the source's nominal tile from
t_ca - 12 d to t_ca + 1 d whose epoch puts the source inside the
planning tile model, classified 'arc' (in-beam, b_e <= 0.1 AU) or
'base' (in tile, off-beam; the nearest `baseline_cap` frames before
the arc). Patches: source (fixed ICRS) + 8 controls at orbit-latitude
offsets (ladder per the config; one-sided toward the tile interior
when the two-sided ladder would cross the seam or the outer edge) at
the arc midpoint, fixed ICRS. Gaia G <= 6 within 4 px drops a patch (S2
source exempt). Each patch carries its Gaia DR3 star list (V, B-V,
sep_px) for the stellar template.

Measure: one pass per frame — gate, load, refine WCS, 2-D ZP, top-hat
photometry at every mapped patch (flux normalised to ZP_REF); the
per-frame footprint test uses the frame's own WCS and DSTART/DSTOP
(12-px margin). Every 12th frame dumps calibrators and a
stamp-response test.

Reduce: per (unit, event, patch): E' = F_arc - median(F_base) (star-
fixed differential, primary; >= min_base_epochs same-tile baseline
epochs, else the event is `no_baseline`); same-frame differential
D_src = E'_src - Q(E'_controls) with Q the quadratic-in-orbit-latitude
interpolation at the source offset (controls: own E' minus the
interpolation of the other seven); z per event =
mean(D_arc)/[MAD-SD(D_arc)/sqrt(n)]; S_stack, S_event, S_pulse
(two-frame persistence); T = max over the 8 control patches. The
Gaia-template route (E = F - (a + b T)) is computed alongside and
reported as the secondary excess. Planet proximity vetoes and the
frozen recon exclusions are applied here.
"""

from __future__ import annotations

import json
import os
import sys
import threading
import time
from collections import defaultdict
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solo_geometry as G
import solohi_lib as W

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "solohi-crossings"
CFG_FILE = os.environ.get("SOLOHI_FREEZE", "threshold_freeze_v1_2.json")
CFG = json.loads((SURV / "configs" / CFG_FILE).read_text())
RUNS = REPO / "runs" / "solohi-crossings"
SER = RUNS / "series"
EVENTS = REPO / "crossings" / "solo_v1" / "events.ecsv"
INV = RUNS / "coverage" / "soar_l2_inventory.json"
COVERAGE = SURV / "results" / "coverage_v1.json"
GAIA_CACHE = RUNS / "coverage" / "gaia_cache.json"

R_AU = CFG["geometry"]["in_beam_radius_au"]
OFFS = CFG["controls"]["orbit_lat_offsets_deg"]
OFFS_ONE = [1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0]
MARGIN_PX = 12
DUMP_EVERY = 12
SEED = CFG["seed"]
BASE_CAP = CFG["volume"]["baseline_cap_frames"]
PRE_D, POST_D = 12.0, 1.0


def unit_events():
    cov = json.loads(COVERAGE.read_text())
    out = defaultdict(list)
    for r in cov["rows"]:
        if r["arc_frames"] >= CFG["gates"]["event"]["min_usable_arc_epochs"]:
            out[f'{r["target_id"]}:{r["channel"]}'].append(r["event_id"])
    return out


def from_ram(lam_deg, bet_deg, mjd):
    u, t, n = G.ram_frame(mjd)
    u, t, n = u[0], t[0], n[0]
    ln, lt = np.radians(lam_deg), np.radians(bet_deg)
    return np.cos(lt) * np.cos(ln) * u + np.cos(lt) * np.sin(ln) * t + np.sin(lt) * n


def inventory_frames():
    inv = json.loads(INV.read_text())
    ci = {c: i for i, c in enumerate(inv["columns"])}
    out = {"1": [], "2": []}
    for r in inv["rows"]:
        d = r[ci["descriptor"]]
        if d in ("solohi-1ft", "solohi-2ft"):
            out[d[7]].append((r[ci["filename"]], Time(r[ci["begin_time"]]).mjd, r[ci["begin_time"]][:10].replace("-", "")))
    for k in out:
        out[k].sort(key=lambda x: x[1])
    return out


def ladder_for(tile: str, bet0: float) -> list:
    b0, b1 = G.TILES[tile]["bet_deg"]
    if bet0 - 6.0 >= b0 + 0.3 and bet0 + 6.0 <= b1 - 0.3:
        return list(OFFS)
    # one-sided toward the tile interior
    centre = 0.5 * (b0 + b1)
    sign = 1.0 if bet0 < centre else -1.0
    return [sign * o for o in OFFS_ONE]


def build_index(role: str) -> None:
    t = Table.read(EVENTS)
    by_id = {str(r["event_id"]): r for r in t}
    frames = inventory_frames()
    fm = {k: np.array([f[1] for f in v]) for k, v in frames.items()}
    orbits = G.perihelia()
    cache = json.loads(GAIA_CACHE.read_text()) if GAIA_CACHE.exists() else {}
    star_tab = CFG.get("s2_stellar_template", {}).get("targets_v_bv", {})
    gmask = CFG["gates"]["patch"]["gaia_g_mask"]
    mask_px = CFG["gates"]["patch"]["mask_radius_px"]

    ue = unit_events()
    units = CFG["split"][role]
    index: dict[str, list] = defaultdict(list)
    patches, evmeta = {}, {}
    todo = []
    for ukey in units:
        tid, ch = ukey.split(":")
        for eid in ue.get(ukey, []):
            ev = by_id[eid]
            e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
                 "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
            mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
            s_hat = G.source_direction(e, ch)
            # nominal tile at the arc midpoint (from the census geometry)
            mj_probe = mjd_ca + np.arange(-PRE_D, POST_D, 1 / 24)
            b_p, _, _ = G.axis_geometry(e, mj_probe)
            _, lam_p, bet_p = G.solohi_coords(s_hat, mj_probe)
            tl_p = G.tile_of(lam_p, bet_p)
            arc_p = (b_p <= R_AU) & np.isin(tl_p, ["1", "2"])
            if not arc_p.any():
                continue
            tile = str(np.unique(tl_p[arc_p])[0])
            lo, hi = np.searchsorted(fm[tile], mjd_ca - PRE_D), np.searchsorted(fm[tile], mjd_ca + POST_D)
            if hi - lo == 0:
                continue
            mj = fm[tile][lo:hi]
            b_e, along, rh = G.axis_geometry(e, mj)
            eps, lam, bet = G.solohi_coords(s_hat, mj)
            intile = G.in_tile(lam, bet, G.TILES[tile])
            kinds = np.where(~intile, None, np.where(b_e <= R_AU, "arc", "base"))
            arc_idx = [k for k in range(hi - lo) if kinds[k] == "arc"]
            if len(arc_idx) < CFG["gates"]["event"]["min_usable_arc_epochs"]:
                continue
            base_idx = [k for k in range(hi - lo) if kinds[k] == "base" and k < arc_idx[0]]
            base_idx = base_idx[-BASE_CAP:]
            t_mid = float(np.median(mj[arc_idx]))
            _, lam0, bet0 = G.solohi_coords(s_hat, t_mid)
            lam0, bet0 = float(lam0[0]), float(bet0[0])
            ladder = ladder_for(tile, bet0)
            pl = [tuple(map(float, G.vec_to_radec(s_hat)))]
            for off in ladder:
                pl.append(tuple(map(float, G.vec_to_radec(from_ram(lam0, bet0 + off, t_mid)))))
            key = f"{ukey}|{eid}"
            patches[key] = {"radec": pl, "t_mid": t_mid, "ladder": ladder, "tile": tile,
                            "lam0": lam0, "bet0": bet0}
            todo.append((key, tid, ch))
            evmeta[key] = {"mjd_ca": mjd_ca, "tile": tile, "n_arc_listed": len(arc_idx), "n_base_listed": len(base_idx),
                           "orbit": next((x["orbit"] for x in orbits if x["mjd_start"] <= mjd_ca <= x["mjd_end"]), None),
                           "arc_rel_tca_h": [round(float((mj[arc_idx[0]] - mjd_ca) * 24), 2), round(float((mj[arc_idx[-1]] - mjd_ca) * 24), 2)],
                           "arc_b_e_rsun": [round(float(b_e[arc_idx].min() / G.RSUN_AU), 2), round(float(b_e[arc_idx].max() / G.RSUN_AU), 2)],
                           "arc_eps_deg": [round(float(eps[arc_idx].min()), 2), round(float(eps[arc_idx].max()), 2)],
                           "arc_r_au": [round(float(rh[arc_idx].min()), 4), round(float(rh[arc_idx].max()), 4)]}
            for k in arc_idx + base_idx:
                fname, m, day = frames[tile][lo + k]
                index[fname].append({"u": ukey, "e": eid, "k": str(kinds[k]), "day": day})
    # Gaia templates (cached cone queries, radius 0.15 deg ~ 7 px)
    rad = 0.15
    lim_mask = mask_px * W.PIX_ARCSEC / 3600

    def fetch(key_ra_dec):
        ra, dec = key_ra_dec
        ck = f"{ra:.4f},{dec:.4f}"
        if ck not in cache:
            cache[ck] = W.gaia_cone(ra, dec, rad)
        return ck
    keys = sorted({(round(ra, 4), round(dec, 4)) for p in patches.values() for ra, dec in p["radec"]})
    with ThreadPoolExecutor(6) as ex:
        list(ex.map(fetch, keys))
    GAIA_CACHE.write_text(json.dumps(cache))
    for key, tid, ch in todo:
        pt = patches[key]
        masked, stars = [], []
        for pi, (ra, dec) in enumerate(pt["radec"]):
            g = cache[f"{round(ra, 4):.4f},{round(dec, 4):.4f}"]
            v0 = G.radec_to_vec(ra, dec)
            bright_near = False
            for row in g:
                if row[4] <= gmask:
                    sv = G.radec_to_vec(row[0], row[1])
                    if np.degrees(np.arccos(np.clip(sv @ v0, -1, 1))) < lim_mask:
                        bright_near = True
            masked.append(bright_near and not (pi == 0 and ch == "S2"))
            if pi == 0 and ch == "S2" and tid in star_tab:
                vv, bv = star_tab[tid]
                stars.append(W.template_stars_gaia(ra, dec, pt["t_mid"], g, extra=[(vv, bv, 0.0)], exclude_sep_arcsec=20.0))
            else:
                stars.append(W.template_stars_gaia(ra, dec, pt["t_mid"], g))
        pt["masked"], pt["stars"], pt["n_gaia"] = masked, stars, [len(cache[f"{round(ra, 4):.4f},{round(dec, 4):.4f}"]) for ra, dec in pt["radec"]]
    SER.mkdir(parents=True, exist_ok=True)
    (SER / f"index_{role}_v1.json").write_text(json.dumps(
        {"index": index, "patches": patches, "events": evmeta}, indent=0) + "\n")
    n_arc = sum(1 for v in index.values() if any(t["k"] == "arc" for t in v))
    n_src_masked = sum(1 for p in patches.values() if p["masked"][0])
    days = {v[0]["day"] for v in index.values()}
    print(f"{role}: {len(index)} frames ({n_arc} arc frames) over {len(days)} days; "
          f"{len(patches)} unit-events; source patches Gaia-masked: {n_src_masked}; "
          f"one-sided ladders: {sum(1 for p in patches.values() if p['ladder'][0] * p['ladder'][1] > 0)}")


# ----------------------------------------------------------------- measure

def _stamp_response(fr: W.Frame, rng: np.random.Generator, n: int = 5) -> list:
    from scipy import ndimage
    out = []
    sig = W.PSF_FWHM_PX / 2.355
    yy, xx = np.mgrid[-6:7, -6:7]
    for _ in range(n):
        x, y = rng.uniform(40, fr.nx - 40), rng.uniform(40, fr.ny - 40)
        zp = fr.zp_at(x, y)
        if not np.isfinite(zp):
            continue
        f_norm = 10 ** (rng.uniform(-0.6, 0.6))
        f_unit = f_norm * 10 ** (0.4 * (zp - W.ZP_REF))
        ix, iy = int(round(x)), int(round(y))
        stamp = np.exp(-0.5 * ((xx - (x - ix)) ** 2 + (yy - (y - iy)) ** 2) / sig ** 2)
        stamp *= f_unit / stamp.sum()
        raw = fr.raw.copy()
        raw[iy - 6:iy + 7, ix - 6:ix + 7] += stamp
        w = 30
        sub = raw[iy - w:iy + w + 1, ix - w:ix + w + 1]
        fill = np.nan_to_num(sub, nan=float(np.nanmedian(sub)))
        hp = sub - ndimage.median_filter(fill, size=W.HP_SIZE)
        f0, _, _ = W.aper_flux(fr.img, x, y)
        f1, _, _ = W.aper_flux(hp, w + (x - ix), w + (y - iy))
        if np.isfinite(f0) and np.isfinite(f1):
            out.append({"x": round(x, 1), "y": round(y, 1), "f_norm": f_norm,
                        "ratio": float((f1 - f0) / f_unit)})
    return out


def measure_frame(args):
    path, tasks, patches, dump = args
    out = {"frame": Path(path).name, "records": []}
    try:
        import warnings
        from astropy.io import fits
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            h = fits.getheader(path)
        reason = W.frame_gate(h)
        if reason:
            out["unusable"] = reason
            return out
        fr = W.load_frame(Path(path))
        ok = W.fit_frame(fr)
        out["fit"] = {"ok": bool(ok), "n_match": fr.n_match, "rms": round(fr.astrom_rms, 3),
                      "zp0": round(float(fr.zp_coef[0]), 4) if fr.zp_coef is not None else None,
                      "zp_mad": round(fr.zp_scatter, 4) if np.isfinite(fr.zp_scatter) else None,
                      "n_calib": fr.n_calib, "mjd": fr.mjd_avg, "r_au": round(fr.r_au, 5), "zp_ref": W.ZP_REF,
                      "xposure": fr.xposure, "gain": fr.gain, "tile": fr.tile, "mode": fr.mode, "roll": round(fr.roll, 2),
                      "sat_masked": round(fr.sat_frac_masked, 4)}
        if not ok or fr.zp_coef is None or fr.zp_scatter / np.sqrt(max(fr.n_calib, 1)) > W.MAX_ZP_UNC:
            out["unusable"] = "astrometry_or_zp"
            return out
        for tk in tasks:
            key = f'{tk["u"]}|{tk["e"]}'
            pt = patches[key]
            for pi, (ra, dec) in enumerate(pt["radec"]):
                if pt["masked"][pi]:
                    continue
                m = W.measure_patch(fr, ra, dec)
                if not W.in_footprint(fr, m["x"], m["y"], MARGIN_PX):
                    continue
                m.update({"u": tk["u"], "e": tk["e"], "k": tk["k"], "p": pi})
                out["records"].append(m)
        if dump:
            cal = np.array(fr.calibrators) if fr.calibrators else np.zeros((0, 6))
            if len(cal):
                zp_at = np.array([fr.zp_at(x, y) for x, y in cal[:, :2]])
                fn = cal[:, 5] * 10 ** (-0.4 * (zp_at - W.ZP_REF))
                out["calibrators"] = [[round(float(a), 1), round(float(b), 1), round(float(v), 2),
                                       round(float(bv), 2), round(float(f), 5)]
                                      for a, b, v, bv, f in zip(cal[:, 0], cal[:, 1], cal[:, 3], cal[:, 4], fn)]
            rng = np.random.default_rng(SEED + int(fr.mjd_avg * 100) % 100000)
            out["stamps"] = _stamp_response(fr, rng)
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
    return out


def run(role: str, workers: int = 12, fetchers: int = 16) -> None:
    idx = json.loads((SER / f"index_{role}_v1.json").read_text())
    patches = idx["patches"]
    out_path = SER / f"measurements_{role}_v1.jsonl"
    done = set()
    if out_path.exists():
        with open(out_path) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                    if rec["frame"].endswith("_missing"):
                        continue
                    done.add(rec["frame"])
                except Exception:
                    pass
    names = sorted(n for n in idx["index"] if n not in done)
    by_day = defaultdict(list)
    for n in names:
        by_day[idx["index"][n][0]["day"]].append(n)
    print(f"{role}: {len(names)} frames to do over {len(by_day)} days ({len(done)} done)", flush=True)
    fetch_log = open(SER / f"fetch_{role}_v1.log", "a")
    t0 = time.time()
    n_done = 0
    with open(out_path, "a") as fh, ProcessPoolExecutor(workers, mp_context=mp.get_context("forkserver")) as ex, \
            ThreadPoolExecutor(fetchers) as tp:
        days = sorted(by_day)

        def fetch_day(day):
            names_d = by_day[day]
            res = list(tp.map(lambda a: (a[1], *W.fetch_frame(day, a[1], source="nrl" if a[0] % 2 == 0 else "soar")),
                              list(enumerate(names_d))))
            for n, p, src in res:
                fetch_log.write(f"{n} {src}\n")
            fetch_log.flush()
            return res
        pending = None
        for i, day in enumerate(days):
            res = pending if pending is not None else fetch_day(day)
            nxt = None
            if i + 1 < len(days):
                holder = {}
                th = threading.Thread(target=lambda: holder.update(r=fetch_day(days[i + 1])))
                th.start()
                nxt = (th, holder)
            args = []
            for j, (n, p, src) in enumerate(res):
                if p is None:
                    fh.write(json.dumps({"frame": n + "_missing", "unusable": "fetch_" + src}) + "\n")
                    continue
                args.append((str(p), idx["index"][n], patches, j % DUMP_EVERY == 0))
            for rec in ex.map(measure_frame, args, chunksize=1):
                fh.write(json.dumps(rec) + "\n")
            fh.flush()
            n_done += len(res)
            for n, p, src in res:
                if p is not None:
                    try:
                        p.unlink()
                    except OSError:
                        pass
            try:
                (W.FRAMES / day).rmdir()
            except OSError:
                pass
            if i % 5 == 0:
                el = time.time() - t0
                print(f"  {day}: {n_done}/{len(names)} frames, {el/60:.1f} min, "
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


INTERP = os.environ.get("SOLOHI_INTERP", "quad")   # 'quad' = frozen v1.x; 'median' = post-blind diagnostic only


def _quad_interp(offs, vals, at):
    """Quadratic-in-latitude interpolation of control values at `at`;
    linear when < 4 points, median when < 3. SOLOHI_INTERP=median gives
    the robust median of the controls (diagnostic, not the frozen rule)."""
    offs, vals = np.asarray(offs, float), np.asarray(vals, float)
    if INTERP == "median":
        return float(np.median(vals)) if len(vals) else np.nan
    if len(offs) < 3:
        return float(np.median(vals)) if len(vals) else np.nan
    deg = 2 if len(offs) >= 4 else 1
    c = np.polyfit(offs, vals, deg)
    return float(np.polyval(c, at))


def _template_response(rows: np.ndarray) -> dict:
    T, F = rows[:, 0], rows[:, 1]
    edges = [0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 1e9]
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (T >= lo) & (T < hi)
        if m.sum() >= 15:
            bins.append((float(np.median(T[m])), float(np.median(F[m])), int(m.sum()), float(_robust_sd(F[m]))))
    if len(bins) < 3:
        return {"a": 0.0, "b": 1.0, "n": int(len(T)), "bins": bins, "note": "too few bins; identity"}
    B = np.array(bins)
    w = np.sqrt(B[:, 2]) / (B[:, 3] + 1e-9)
    A = np.c_[np.ones(len(B)), B[:, 0]]
    coef, *_ = np.linalg.lstsq(A * w[:, None], B[:, 1] * w, rcond=None)
    return {"a": float(coef[0]), "b": float(coef[1]), "n": int(len(T)),
            "bins": [{"T_median": round(t, 4), "F_median": round(f, 4), "n": n, "F_mad": round(sd, 4)} for t, f, n, sd in bins],
            "note": "E = F - (a + b T); a = blank-patch aperture offset, b = template scale"}


def _colour_system(cal_rows) -> dict:
    if len(cal_rows) < 200:
        return {"coeff": W.COLOUR_COEFF, "n": len(cal_rows), "note": "frozen value (dump too small)"}
    c = np.array(cal_rows, float)
    ok = (c[:, 4] > 0) & (c[:, 3] > -0.3) & (c[:, 3] < 1.6) & (c[:, 2] >= CFG["chain"].get("colour_calibrator_vmin", 0.0))
    c = c[ok]
    resid = 2.5 * np.log10(c[:, 4]) + (c[:, 2] - W.COLOUR_COEFF * (c[:, 3] - W.BV_REF)) - W.ZP_REF
    A = np.c_[np.ones(len(c)), c[:, 3] - W.BV_REF]
    coef, *_ = np.linalg.lstsq(A, resid, rcond=None)
    r2 = resid - A @ coef
    return {"coeff": float(W.COLOUR_COEFF + coef[1]), "delta": float(coef[1]), "n": int(len(c)),
            "scatter_mad": _robust_sd(r2),
            "by_v": {f"{lo}-{lo+1}": round(float(np.median(r2[(c[:, 2] >= lo) & (c[:, 2] < lo + 1)])), 3)
                     for lo in range(3, 9) if ((c[:, 2] >= lo) & (c[:, 2] < lo + 1)).sum() >= 10}}


def reduce(role: str, version: str = "v1") -> dict:
    idx = json.loads((SER / f"index_{role}_v1.json").read_text())
    patches, evmeta = idx["patches"], idx["events"]
    excl = {e: [Time(x["utc"]).mjd for x in xs] for e, xs in CFG.get("recon_exclusions", {}).items()}
    gp = CFG["gates"]["epoch"]

    epochs = defaultdict(dict)
    frames_ok, frames_bad = 0, defaultdict(int)
    fit_rows, cal_rows, stamp_rows = [], [], []
    last, n_corrupt = {}, 0
    n_noise_gated = defaultdict(int)
    with open(SER / f"measurements_{role}_v1.jsonl") as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                n_corrupt += 1
                continue
            last[rec["frame"]] = rec
    for name in [n for n in last if n.endswith("_missing") and n[:-8] in last]:
        del last[name]
    for rec in last.values():
        if "unusable" in rec or "error" in rec:
            frames_bad[rec.get("unusable", "error")] += 1
            continue
        frames_ok += 1
        if "fit" in rec:
            fit_rows.append(rec["fit"])
        cal_rows.extend(rec.get("calibrators", []))
        stamp_rows.extend(rec.get("stamps", []))
        for m in rec["records"]:
            if m["valid"] < gp["aperture_finite_frac_min"] or not np.isfinite(m["flux"] or np.nan):
                continue
            if m["err"] is None or not np.isfinite(m["err"]) or m["err"] > gp.get("max_epoch_err_units", np.inf):
                n_noise_gated[m["k"]] += 1           # v1.2 structure-noise gate (the sunward-edge band)
                continue
            epochs[(m["u"], m["e"], m["k"], round(m["mjd"], 5))][m["p"]] = (m["flux"], m["err"])

    colour = _colour_system(cal_rows)
    ccoef = colour["coeff"]

    mjds = np.array(sorted({k[3] for k in epochs}))
    pdirs = G.planet_directions(mjds) if len(mjds) else {}
    lim_prox = 2 * np.sin(np.radians(gp["body_proximity_px"] * W.PIX_ARCSEC / 3600) / 2)
    bright_prox = {b: np.radians(r) for b, r in gp.get("bright_body_proximity_deg", {}).items()}
    n_veto = defaultdict(int)
    keep = {}
    diag = defaultdict(list)
    for key, by_p in epochs.items():
        u, e, k, mjd = key
        i = min(int(np.searchsorted(mjds, mjd)), len(mjds) - 1)
        if e in excl and any(abs(mjd - x) < 0.003 for x in excl[e]):
            n_veto["recon_exclusion"] += 1
            continue
        pt = patches[f"{u}|{e}"]
        pv = np.array([G.radec_to_vec(ra, dec) for ra, dec in pt["radec"]])
        bad, sep_body = None, None
        for body, dirs in pdirs.items():
            if np.any(np.linalg.norm(pv - dirs[i], axis=1) < lim_prox):
                bad = f"prox_{body}"
                break
        if not bad:
            for body, rad in bright_prox.items():
                if body in pdirs:
                    sep = np.arccos(np.clip(pv @ pdirs[body][i], -1, 1))
                    if np.any(sep < rad):
                        bad = f"near_{body}"
                        sep_body = float(np.degrees(sep.min()))
                        break
        cls = bad if bad else "clean"
        if sep_body is not None:
            cls = f"{bad}_{'0-2' if sep_body < 2 else '2-5' if sep_body < 5 else '5-10' if sep_body < 10 else '10+'}deg"
        if k == "arc" and len(by_p) >= 6:
            cs = [pi for pi in range(1, 9) if pi in by_p]
            if len(cs) >= 5:
                E = {c: by_p[c][0] for c in cs}
                for c in cs:
                    diag[cls].append(E[c] - float(np.median([E[o] for o in cs if o != c])))
        if bad:
            n_veto[bad] += 1
            continue
        keep[key] = by_p
    bright_diag = {c: {"n": len(v), "ctrl_diff_mad": _robust_sd(v)} for c, v in diag.items()}

    grouped = defaultdict(dict)
    for (u, e, k, mjd), by_p in keep.items():
        grouped[(u, e, k)][mjd] = by_p

    results, series_dump = {}, {}
    units = sorted({k.rsplit("|", 1)[0] for k in patches})
    min_arc = CFG["gates"]["event"]["min_usable_arc_epochs"]
    min_base = CFG["gates"]["event"]["min_base_epochs"]
    ctrl_idx = list(range(1, 9))
    min_ctrl = CFG["controls"]["min_controls_per_epoch"]
    bs = CFG.get("bright_star_class", {})
    star_tab = CFG.get("s2_stellar_template", {}).get("targets_v_bv", {})
    bright_units, unsearchable = set(), set()
    for ukey in units:
        tid, ch = ukey.split(":")
        if ch == "S2" and tid in star_tab and bs:
            if star_tab[tid][0] < bs["v_unsearchable"]:
                unsearchable.add(ukey)
            elif star_tab[tid][0] < bs["v_self_calibrated"]:
                bright_units.add(ukey)
    # template response on the control ensemble (secondary route)
    resp_rows = []
    for key, pt in patches.items():
        ukey_, eid = key.rsplit("|", 1)
        arc = grouped.get((ukey_, eid, "arc"), {})
        if len(arc) < min_arc:
            continue
        for pi in range(1, 9):
            f = [v[pi][0] for v in arc.values() if pi in v]
            if len(f) >= min_arc:
                resp_rows.append((W.template_flux(pt["stars"][pi], ccoef), float(np.mean(f))))
    template_response = _template_response(np.array(resp_rows)) if resp_rows else {"a": 0.0, "b": 1.0, "n": 0}
    ra_, rb_ = template_response["a"], template_response["b"]

    for ukey in units:
        per_patch, ev_detail = {}, {}
        if ukey in unsearchable:
            results[ukey] = {"status": "bright_star_unsearchable", "events": {}, "n_included": 0}
            continue
        for key, pt in patches.items():
            if not key.startswith(ukey + "|"):
                continue
            eid = key.split("|")[1]
            arc = grouped.get((ukey, eid, "arc"), {})
            base = grouped.get((ukey, eid, "base"), {})
            if len(arc) < min_arc:
                ev_detail[eid] = {"status": "insufficient_epochs", "n_arc": len(arc), "n_base": len(base)}
                continue
            offs = {pi: (0.0 if pi == 0 else pt["ladder"][pi - 1]) for pi in range(9)}
            # star-fixed baseline per patch
            Fb = {}
            for pi in range(9):
                fb = [v[pi][0] for v in base.values() if pi in v]
                if len(fb) >= min_base:
                    Fb[pi] = float(np.median(fb))
            route = "star_fixed" if 0 in Fb and sum(1 for c in ctrl_idx if c in Fb) >= min_ctrl else None
            T = {pi: ra_ + rb_ * W.template_flux(pt["stars"][pi], ccoef) for pi in range(9)}
            if route is None:
                ev_detail[eid] = {"status": "no_baseline", "n_arc": len(arc), "n_base": len(base),
                                  "n_base_src": len([1 for v in base.values() if 0 in v])}
                continue
            D, Dt = {pi: [] for pi in range(9)}, {pi: [] for pi in range(9)}
            mj_used = []
            for m in sorted(arc):
                v = arc[m]
                E = {pi: v[pi][0] - Fb[pi] for pi in range(9) if pi in v and pi in Fb}
                Et = {pi: v[pi][0] - T[pi] for pi in range(9) if pi in v}
                cs = [pi for pi in ctrl_idx if pi in E]
                if 0 not in E or len(cs) < min_ctrl:
                    continue
                mj_used.append(m)
                D[0].append(E[0] - _quad_interp([offs[c] for c in cs], [E[c] for c in cs], 0.0))
                for c in cs:
                    oth = [o for o in cs if o != c]
                    D[c].append(E[c] - _quad_interp([offs[o] for o in oth], [E[o] for o in oth], offs[c]))
                cst = [pi for pi in ctrl_idx if pi in Et]
                if 0 in Et and len(cst) >= min_ctrl:
                    Dt[0].append(Et[0] - _quad_interp([offs[c] for c in cst], [Et[c] for c in cst], 0.0))
            if len(D[0]) < min_arc:
                ev_detail[eid] = {"status": "insufficient_after_veto", "n_arc": len(D[0])}
                continue
            series_dump[key] = {"mjd": [round(m, 5) for m in mj_used],
                                "D": {str(pi): [round(x, 5) for x in D[pi]] for pi in range(9) if len(D[pi]) == len(mj_used)},
                                "D_template_src": [round(x, 5) for x in Dt[0]]}
            for pi in range(9):
                d = np.array(D[pi])
                if len(d) < min_arc:
                    continue
                sd = _robust_sd(d)
                if not sd > 0:
                    continue
                pp = per_patch.setdefault(pi, {"z": [], "pulse": [], "ev": []})
                pp["z"].append(float(np.mean(d) / (sd / np.sqrt(len(d)))))
                x = (d - np.median(d)) / sd
                pp["pulse"].append(float(np.max(np.minimum(x[:-1], x[1:]))))
                pp["ev"].append(eid)
            d0 = np.array(D[0])
            det = {"status": "included", "route": route, "n_arc": len(d0), "n_base": len(base),
                   "z": round(float(np.mean(d0) / (_robust_sd(d0) / np.sqrt(len(d0)))), 2),
                   "mean_D": round(float(np.mean(d0)), 5), "sd_D": round(_robust_sd(d0), 5),
                   "baseline_src": round(Fb[0], 5), "template_src": round(T[0], 5),
                   "tile": pt["tile"], "orbit": evmeta[key]["orbit"], "mjd_ca": round(evmeta[key]["mjd_ca"], 3),
                   "pulse_mjd": round(float(mj_used[int(np.argmax(np.minimum(d0[:-1], d0[1:])))]), 4)}
            if len(Dt[0]) >= min_arc:
                dt = np.array(Dt[0])
                det["z_template"] = round(float(np.mean(dt) / (_robust_sd(dt) / np.sqrt(len(dt)) + 1e-30)), 2)
                det["mean_D_template"] = round(float(np.mean(dt)), 5)
            ev_detail[eid] = det
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
        res = {"status": "searched", "per_patch": {str(k): v for k, v in stats.items()}, "events": ev_detail,
               "n_included": stats[0]["n_events"]}
        src = stats[0]
        ctrl = [v for k, v in stats.items() if k != 0]
        for s in ("S_stack", "S_event", "S_pulse"):
            vals = [c[s] for c in ctrl]
            Tm = max(vals) if vals else 0.0
            mad = _robust_sd(vals) if len(vals) >= 3 else None
            res[s] = {"S": round(src[s], 3), "T": round(float(Tm), 3), "n_controls": len(vals),
                      "exceeds": bool(src[s] > max(Tm, 0.0)),
                      "S_over_ctrl_mad": round(src[s] / mad, 2) if mad else None}
        results[ukey] = res

    n_trials = sum(3 for r in results.values() if r["status"] == "searched")
    n_exc = sum(1 for r in results.values() if r["status"] == "searched"
                for s in ("S_stack", "S_event", "S_pulse") if r[s]["exceeds"])
    fit_ok = [f for f in fit_rows if f.get("ok")]
    out = {"role": role, "construction": CFG["version"],
           "frames": {"usable": frames_ok, "unusable": dict(frames_bad), "corrupt_lines_skipped": n_corrupt,
                      "median_n_match": float(np.median([f["n_match"] for f in fit_ok])) if fit_ok else None,
                      "median_astrom_rms_px": float(np.median([f["rms"] for f in fit_ok])) if fit_ok else None,
                      "median_zp0": float(np.median([f["zp0"] for f in fit_ok if f["zp0"]])) if fit_ok else None,
                      "median_zp_mad": float(np.median([f["zp_mad"] for f in fit_ok if f["zp_mad"]])) if fit_ok else None},
           "epoch_vetoes": dict(n_veto), "noise_gated_records": dict(n_noise_gated),
           "bright_body_diagnostic": bright_diag, "colour_system": colour,
           "bright_star_class": {"self_calibrated_units": sorted(bright_units), "unsearchable": sorted(unsearchable),
                                 "note": "star-fixed route: the star's own baseline is the template; no scale needed"},
           "template_response": template_response,
           "stamp_response": {"n": len(stamp_rows),
                              "median_ratio": float(np.median([s["ratio"] for s in stamp_rows])) if stamp_rows else None,
                              "mad": _robust_sd([s["ratio"] for s in stamp_rows]) if stamp_rows else None},
           "trials": {"searched": n_trials, "exceedances": n_exc, "expected_control_crossings": round(n_trials / 9, 2)},
           "units": results}
    (SURV / "results" / f"{role}_search_{version}.json").write_text(json.dumps(out, indent=1) + "\n")
    (SER / f"series_{role}_{version}.json").write_text(json.dumps(series_dump) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("units", "bright_body_diagnostic")}, indent=1))
    for u, r in results.items():
        if r["status"] != "searched":
            print(u, "->", r["status"], r.get("n_included"), {k: v["status"] for k, v in r.get("events", {}).items()})
            continue
        line = f'{u}: n_ev={r["per_patch"]["0"]["n_events"]} zmed={r["per_patch"]["0"]["z_median"]:+.2f} zmad={r["per_patch"]["0"]["z_mad"]:.2f}'
        for s in ("S_stack", "S_event", "S_pulse"):
            line += f' | {s} S={r[s]["S"]:.2f} T={r[s]["T"]:.2f}{" EXC" if r[s]["exceeds"] else ""}'
        print(line)
    return out


if __name__ == "__main__":
    cmd, role = sys.argv[1], sys.argv[2]
    if cmd == "index":
        build_index(role)
    elif cmd == "run":
        run(role, workers=int(sys.argv[3]) if len(sys.argv) > 3 else 12)
    elif cmd == "reduce":
        reduce(role, sys.argv[3] if len(sys.argv) > 3 else "v1")
