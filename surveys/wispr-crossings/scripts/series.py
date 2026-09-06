"""PSP/WISPR-I search driver (threshold freeze v1.x).

Stages (role = dev | confirmatory):
  python series.py index <role>    # frame -> task index + patch table
  python series.py run <role>      # fetch (NRL L3) + measure per day, purge
  python series.py reduce <role>   # series -> statistics -> results json

Index: for every covered event of the role's units, every synoptic
WISPR-I L3 frame (inventory snapshot) within +-9 d of t_ca whose epoch
puts the source inside the planning field model, classified 'arc'
(in-beam, b_e <= 0.1 AU) or 'base' (in field, out of beam). Patches:
source (fixed ICRS) + 8 controls at orbit-latitude offsets (ladder per
the config) at the arc midpoint, fixed ICRS. Tycho-2 VT <= 7 within 4 px
drops a patch (S2 source exempt). Each patch carries its catalogue
star list (V, B-V, sep_px) for the stellar template.

Measure: one pass per frame — gate, load (L3 rescaled), refine WCS,
2-D ZP, top-hat photometry at every mapped patch (flux normalised to
ZP_REF); the per-frame footprint test uses the frame's own WCS (12-px
margin). Every 12th frame dumps calibrators (colour system, era
stability) and a stamp-response test.

Reduce: per (unit, event, patch): E = F - F_template (colour
coefficient re-measured from the calibrator dump); same-frame
differential D_src = E_src - median_k E_k; z per event =
mean(D_arc)/[MAD-SD(D_arc)/sqrt(n)]; S_stack, S_event, S_pulse; T = max
over the 8 control patches. Planet/comet proximity, bright-body
in-FOV veto, encounter/extended-campaign gate and the frozen recon
exclusions are applied here. Star-fixed differential reported where
a baseline exists.
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
import psp_geometry as G
import wispr_lib as W

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "wispr-crossings"
CFG_FILE = os.environ.get("WISPR_FREEZE", "threshold_freeze_v1_5.json")
CFG = json.loads((SURV / "configs" / CFG_FILE).read_text())
RUNS = REPO / "runs" / "wispr-crossings"
SER = RUNS / "series"
EVENTS = REPO / "crossings" / "psp_v1" / "events.ecsv"
INV_L3 = RUNS / "coverage" / "l3_inventory.json"
COVERAGE = SURV / "results" / "coverage_v1.json"

R_AU = CFG["geometry"]["in_beam_radius_au"]
OFFS = CFG["controls"]["orbit_lat_offsets_deg"]
OFFS_NORTH = [-1.5, -3.0, -4.5, -6.0, -7.5, -9.0, -10.5, -12.0]
MARGIN_PX = 12
DUMP_EVERY = 12
SEED = CFG["seed"]
STEP = 10 / 1440.0
FOV_I = G.WISPR_I


def stamp_mjd(s):
    return Time(f"{s[:4]}-{s[4:6]}-{s[6:8]}T{s[9:11]}:{s[11:13]}:{s[13:15]}").mjd


def unit_events():
    """{unit_key: [event_id, ...]} from the coverage rows (WISPR-I arc rows)."""
    cov = json.loads(COVERAGE.read_text())
    out = defaultdict(list)
    for r in cov["rows"]:
        if r["wispr_i"]["arc_frames"] >= 10:
            out[f'{r["target_id"]}:{r["channel"]}'].append(r["event_id"])
    return out


def from_ram(lam_deg, bet_deg, mjd):
    """ICRS unit vector of the direction at (ram longitude, orbit latitude)."""
    u, t, n = G.ram_frame(mjd)
    u, t, n = u[0], t[0], n[0]
    ln, lt = np.radians(lam_deg), np.radians(bet_deg)
    return np.cos(lt) * np.cos(ln) * u + np.cos(lt) * np.sin(ln) * t + np.sin(lt) * n


def build_index(role: str) -> None:
    t = Table.read(EVENTS)
    by_id = {str(r["event_id"]): r for r in t}
    inv = json.loads(INV_L3.read_text())
    frames = []
    for od, files in inv["days"].items():
        for f in files or []:
            if f["wxyz"][0] == "1" and f["wxyz"][3] in "12":
                frames.append((f["file"], stamp_mjd(f["stamp"]), od))
    frames.sort(key=lambda x: x[1])
    fm = np.array([f[1] for f in frames])
    enc = G.encounters()
    tyc = np.load(W.TYCHO)
    tv = G.radec_to_vec(tyc["ra"], tyc["dec"])
    tyc_tree = cKDTree(tv)
    sel = tyc["vt"] <= CFG["gates"]["patch"]["tycho_vt_mask"]
    mask_tree = cKDTree(tv[sel])
    lim = 2 * np.sin(np.radians(CFG["gates"]["patch"]["mask_radius_px"] * W.PIX_ARCSEC / 3600) / 2)
    hra, hdec = W.catalog_at(60000.0)
    hip_tree = cKDTree(G.radec_to_vec(hra, hdec))

    ue = unit_events()
    units = CFG["split"][role]
    index: dict[str, list] = defaultdict(list)
    patches, evmeta = {}, {}
    for ukey in units:
        tid, ch = ukey.split(":")
        for eid in ue.get(ukey, []):
            ev = by_id[eid]
            e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
                 "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
            mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
            s_hat = G.source_direction(e, ch)
            lo, hi = np.searchsorted(fm, mjd_ca - 9), np.searchsorted(fm, mjd_ca + 9)
            if hi - lo == 0:
                continue
            mj = fm[lo:hi]
            b_e, along, rh = G.axis_geometry(e, mj)
            eps, lam, bet = G.wispr_coords(s_hat, mj)
            infov = G.in_wispr(lam, bet, FOV_I)
            kinds = np.where(~infov, None, np.where(b_e <= R_AU, "arc", "base"))
            arc_idx = [k for k in range(hi - lo) if kinds[k] == "arc"]
            if len(arc_idx) < CFG["gates"]["event"]["min_usable_arc_epochs"]:
                continue
            t_mid = float(np.median(mj[arc_idx]))
            _, lam0, bet0 = G.wispr_coords(s_hat, t_mid)
            lam0, bet0 = float(lam0[0]), float(bet0[0])
            ladder = OFFS_NORTH if bet0 > 9.0 else OFFS
            pl = [tuple(map(float, G.vec_to_radec(s_hat)))]
            for off in ladder:
                pl.append(tuple(map(float, G.vec_to_radec(from_ram(lam0, bet0 + off, t_mid)))))
            masked, stars = [], []
            star_tab = CFG.get("s2_stellar_template", {}).get("targets_v_bv", {})
            for pi, (ra, dec) in enumerate(pl):
                v = G.radec_to_vec(ra, dec)
                dd, _ = mask_tree.query(v, k=1)
                masked.append(bool(dd < lim) and not (pi == 0 and ch == "S2"))
                if pi == 0 and ch == "S2" and tid in star_tab:
                    vv, bv = star_tab[tid]
                    stars.append(W.template_stars(ra, dec, t_mid, tyc, tyc_tree, hip_tree,
                                                  extra=[(vv, bv, 0.0)], exclude_sep_arcsec=20.0))
                else:
                    stars.append(W.template_stars(ra, dec, t_mid, tyc, tyc_tree, hip_tree))
            key = f"{ukey}|{eid}"
            patches[key] = {"radec": pl, "masked": masked, "t_mid": t_mid, "ladder": ladder,
                            "lam0": lam0, "bet0": bet0, "stars": stars}
            evmeta[key] = {"mjd_ca": mjd_ca, "n_arc_listed": len(arc_idx),
                           "n_base_listed": int(np.sum(kinds == "base")),
                           "encounter": next((x["encounter"] for x in enc if x["mjd_start"] - 9 <= mjd_ca <= x["mjd_end"] + 9), None),
                           "arc_rel_tca_h": [round(float((mj[arc_idx[0]] - mjd_ca) * 24), 2), round(float((mj[arc_idx[-1]] - mjd_ca) * 24), 2)],
                           "arc_b_e_rsun": [round(float(b_e[arc_idx].min() / G.RSUN_AU), 2), round(float(b_e[arc_idx].max() / G.RSUN_AU), 2)],
                           "arc_r_au": [round(float(rh[arc_idx].min()), 4), round(float(rh[arc_idx].max()), 4)]}
            for k in range(hi - lo):
                if kinds[k]:
                    fname, m, od = frames[lo + k]
                    index[fname].append({"u": ukey, "e": eid, "k": str(kinds[k]), "od": od})
    SER.mkdir(parents=True, exist_ok=True)
    (SER / f"index_{role}_v1.json").write_text(json.dumps(
        {"index": index, "patches": patches, "events": evmeta}, indent=0) + "\n")
    n_arc = sum(1 for v in index.values() if any(t["k"] == "arc" for t in v))
    n_src_masked = sum(1 for p in patches.values() if p["masked"][0])
    days = {v[0]["od"] for v in index.values()}
    print(f"{role}: {len(index)} frames ({n_arc} arc frames) over {len(days)} orbit-days; "
          f"{len(patches)} unit-events; source patches Tycho-masked: {n_src_masked}")


# ----------------------------------------------------------------- measure

def _stamp_response(fr: W.Frame, rng: np.random.Generator, n: int = 5) -> list:
    from scipy import ndimage
    out = []
    sig = W.PSF_FWHM_PX / 2.355
    yy, xx = np.mgrid[-6:7, -6:7]
    for _ in range(n):
        x, y = rng.uniform(40, W.NX - 40), rng.uniform(40, W.NY - 40)
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
                      "xposure": fr.xposure, "gain": fr.gain, "code": fr.code}
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
                if not W.in_footprint(m["x"], m["y"], MARGIN_PX):
                    continue
                m.update({"u": tk["u"], "e": tk["e"], "k": tk["k"], "p": pi})
                out["records"].append(m)
        if dump:
            cal = np.array(fr.calibrators) if fr.calibrators else np.zeros((0, 6))
            if len(cal):
                zp_at = np.array([fr.zp_at(x, y) for x, y in cal[:, :2]])
                fn = cal[:, 5] * 10 ** (-0.4 * (zp_at - W.ZP_REF))
                out["calibrators"] = [[round(float(a), 1), round(float(b), 1), round(float(v), 2),
                                       round(float(bv), 2), round(float(f), 5)]   # f in 1e-12 units
                                      for a, b, v, bv, f in zip(cal[:, 0], cal[:, 1], cal[:, 3], cal[:, 4], fn)]
            rng = np.random.default_rng(SEED + int(fr.mjd_avg * 100) % 100000)
            out["stamps"] = _stamp_response(fr, rng)
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
    return out


def run(role: str, workers: int = 16, fetchers: int = 4) -> None:
    idx = json.loads((SER / f"index_{role}_v1.json").read_text())
    patches = idx["patches"]
    out_path = SER / f"measurements_{role}_v1.jsonl"
    done = set()
    if out_path.exists():
        with open(out_path) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                    if rec["frame"].endswith("_missing") or rec.get("unusable") == "astrometry_or_zp":
                        continue                                  # failed fetches / gate rejects are retried on a rerun
                    done.add(rec["frame"])
                except Exception:
                    pass
    names = sorted(n for n in idx["index"] if n not in done)
    by_day = defaultdict(list)
    for n in names:
        by_day[idx["index"][n][0]["od"]].append(n)
    print(f"{role}: {len(names)} frames to do over {len(by_day)} orbit-days ({len(done)} done)", flush=True)
    fetch_log = open(SER / f"fetch_{role}_v1.log", "a")
    t0 = time.time()
    n_done = 0
    with open(out_path, "a") as fh, ProcessPoolExecutor(workers, mp_context=mp.get_context("forkserver")) as ex, \
            ThreadPoolExecutor(fetchers) as tp:
        days = sorted(by_day)

        def fetch_day(od):
            res = list(tp.map(lambda n: (n, *W.fetch_frame(od, n)), by_day[od]))
            for n, p, src in res:
                fetch_log.write(f"{n} {src}\n")
            fetch_log.flush()
            return res
        pending = None
        for i, od in enumerate(days):
            res = pending if pending is not None else fetch_day(od)
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
                (W.FRAMES / od.split("/")[1]).rmdir()
            except OSError:
                pass
            if i % 3 == 0:
                el = time.time() - t0
                print(f"  {od}: {n_done}/{len(names)} frames, {el/60:.1f} min, "
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


def _template_response(rows: np.ndarray) -> dict:
    """<F_arc> = a + b T on the control-patch ensemble: medians of <F> in
    template bins (T in flux units), then a weighted least-squares line
    through the bin medians (bins with >= 15 patches)."""
    T, F = rows[:, 0], rows[:, 1]
    edges = [0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 1e9]
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (T >= lo) & (T < hi)
        if m.sum() >= 15:
            bins.append((float(np.median(T[m])), float(np.median(F[m])), int(m.sum()),
                         float(_robust_sd(F[m]))))
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
    """Ensemble colour coefficient from the calibrator dump: regress the
    ZP residual (2.5 log f_norm + V - ZP_REF, formed with the frozen
    coefficient) on (B-V - 0.65) — the correction to COLOUR_COEFF."""
    if len(cal_rows) < 200:
        return {"coeff": W.COLOUR_COEFF, "n": len(cal_rows), "note": "frozen value (dump too small)"}
    c = np.array(cal_rows, float)   # x, y, V, bv, f_norm
    ok = (c[:, 4] > 0) & (c[:, 3] > -0.3) & (c[:, 3] < 1.6)
    c = c[ok]
    resid = 2.5 * np.log10(c[:, 4]) + (c[:, 2] - W.COLOUR_COEFF * (c[:, 3] - W.BV_REF)) - W.ZP_REF
    A = np.c_[np.ones(len(c)), c[:, 3] - W.BV_REF]
    coef, *_ = np.linalg.lstsq(A, resid, rcond=None)
    # resid = a + d (B-V - 0.65) => the total coefficient is COLOUR_COEFF + d
    r2 = resid - A @ coef
    return {"coeff": float(W.COLOUR_COEFF + coef[1]), "delta": float(coef[1]), "n": int(len(c)),
            "scatter_mad": _robust_sd(r2),
            "by_v": {f"{lo}-{lo+1}": round(float(np.median(r2[(c[:, 2] >= lo) & (c[:, 2] < lo + 1)])), 3)
                     for lo in range(3, 8) if ((c[:, 2] >= lo) & (c[:, 2] < lo + 1)).sum() >= 10}}


def reduce(role: str, version: str = "v1") -> dict:
    idx = json.loads((SER / f"index_{role}_v1.json").read_text())
    patches, evmeta = idx["patches"], idx["events"]
    excl = {e: [Time(x["utc"]).mjd for x in xs] for e, xs in CFG["recon_exclusions"].items()}
    gp = CFG["gates"]["epoch"]
    enc = G.encounters()

    epochs = defaultdict(dict)
    frames_ok, frames_bad = 0, defaultdict(int)
    fit_rows, cal_rows, stamp_rows = [], [], []
    # last record per frame wins (re-measured frames after a gate amendment)
    last, n_corrupt = {}, 0
    with open(SER / f"measurements_{role}_v1.jsonl") as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                n_corrupt += 1               # truncated line (power loss mid-write)
                continue
            last[rec["frame"]] = rec
    for name in [n for n in last if n.endswith("_missing") and n[:-8] in last]:
        del last[name]                           # superseded fetch-failure records
    for rec in last.values():
        if True:
            if "unusable" in rec or "error" in rec:
                frames_bad[rec.get("unusable", "error")] += 1
                continue
            frames_ok += 1
            # legacy-unit records (dev run launched under ZP_REF -21.35): flux x 1e12
            fscale = 1.0 if rec.get("fit", {}).get("zp_ref") == W.ZP_REF else 10 ** (0.4 * (W.ZP_REF + 21.35))
            if "fit" in rec:
                fit_rows.append(rec["fit"])
            if fscale == 1.0:
                cal_rows.extend(rec.get("calibrators", []))
            stamp_rows.extend(rec.get("stamps", []))
            for m in rec["records"]:
                if m["valid"] < gp["aperture_finite_frac_min"] or not np.isfinite(m["flux"] or np.nan):
                    continue
                epochs[(m["u"], m["e"], m["k"], round(m["mjd"], 5))][m["p"]] = (m["flux"] * fscale, m["err"] * fscale)

    colour = _colour_system(cal_rows)
    ccoef = colour["coeff"]

    # vetoes per epoch: proximity (12 px) to any body, bright-body in-FOV,
    # encounter gate, recon exclusion
    mjds = np.array(sorted({k[3] for k in epochs}))
    pdirs = G.planet_directions(mjds) if len(mjds) else {}
    lim_prox = 2 * np.sin(np.radians(gp["body_proximity_px"] * W.PIX_ARCSEC / 3600) / 2)
    in_enc = np.zeros(len(mjds), bool)
    for x in enc:
        in_enc |= (mjds >= x["mjd_start"]) & (mjds <= x["mjd_end"])
    bright_prox = {b: np.radians(r) for b, r in gp.get("bright_body_proximity_deg", {}).items()}
    bright_in_fov = {}
    if not bright_prox:                      # v1.0-v1.2 frame-wide veto
        for body in gp.get("bright_body_in_fov_veto", []):
            if body in pdirs:
                lam = np.array([G.wispr_coords(pdirs[body][i], mjds[i])[1][0] for i in range(len(mjds))])
                bet = np.array([G.wispr_coords(pdirs[body][i], mjds[i])[2][0] for i in range(len(mjds))])
                bright_in_fov[body] = G.in_wispr(lam, bet, FOV_I)
    n_veto = defaultdict(int)
    keep = {}
    diag = defaultdict(list)          # veto class -> control same-frame differentials (arc epochs)
    for key, by_p in epochs.items():
        u, e, k, mjd = key
        i = int(np.searchsorted(mjds, mjd))
        i = min(i, len(mjds) - 1)
        if e in excl and any(abs(mjd - x) < 0.003 for x in excl[e]):
            n_veto["recon_exclusion"] += 1
            continue
        pt = patches[f"{u}|{e}"]
        pv = np.array([G.radec_to_vec(ra, dec) for ra, dec in pt["radec"]])
        bad = None
        for body, dirs in pdirs.items():
            if np.any(np.linalg.norm(pv - dirs[i], axis=1) < lim_prox):
                bad = f"prox_{body}"
                break
        sep_body = None
        if not bad:
            for body, rad in bright_prox.items():           # v1.3 proximity veto
                if body in pdirs:
                    sep = np.arccos(np.clip(pv @ pdirs[body][i], -1, 1))
                    if np.any(sep < rad):
                        bad = f"near_{body}"
                        sep_body = float(np.degrees(sep.min()))
                        break
            for body, flags in bright_in_fov.items():
                if flags[i]:
                    bad = f"infov_{body}"
                    sep_body = float(np.degrees(np.arccos(np.clip(pv[0] @ pdirs[body][i], -1, 1))))
                    break
        cls = bad if bad else "clean"
        if sep_body is not None:
            cls = f"{bad}_{'0-2' if sep_body < 2 else '2-5' if sep_body < 5 else '5-10' if sep_body < 10 else '10+'}deg"
        if k == "arc" and len(by_p) >= 6:
            cs = [pi for pi in range(1, 9) if pi in by_p]
            if len(cs) >= 5:
                E = {c: by_p[c][0] for c in cs}
                encn = evmeta[f"{u}|{e}"]["encounter"]
                for c in cs:
                    diag[cls].append(E[c] - float(np.median([E[o] for o in cs if o != c])))
                    diag[f"{encn}|{cls.split('_')[0] if cls != 'clean' else 'clean'}"].append(diag[cls][-1])
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
    # v1.3 bright-star S2 class: self-calibrated template scale k_unit =
    # median over the unit's events of <F_src>/T_src; V < 3 not searchable
    bs = CFG.get("bright_star_class", {})
    star_tab = CFG.get("s2_stellar_template", {}).get("targets_v_bv", {})
    bright_units, unsearchable, k_unit = set(), set(), {}
    for ukey in units:
        tid, ch = ukey.split(":")
        if ch == "S2" and tid in star_tab and bs:
            if star_tab[tid][0] < bs["v_unsearchable"]:
                unsearchable.add(ukey)
            elif star_tab[tid][0] < bs["v_self_calibrated"]:
                bright_units.add(ukey)
    # v1.4 template response: robust linear fit <F_arc> = a + b T over the
    # CONTROL patches of every unit-event (source patches excluded), in
    # template-brightness bins; frozen from the ensemble at reduce
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
    qgroup = {}
    for grp, q in CFG["era"]["q_groups_rsun"].items():
        a_, b_ = grp.split("-")
        for n_ in range(int(a_[1:]), int(b_[1:]) + 1):
            qgroup[f"E{n_:02d}"] = grp
    k_group = {}
    for ukey in bright_units:
        ratios, by_grp = [], defaultdict(list)
        for key, pt in patches.items():
            if not key.startswith(ukey + "|"):
                continue
            eid = key.split("|")[1]
            arc = grouped.get((ukey, eid, "arc"), {})
            if len(arc) < min_arc:
                continue
            T0 = ra_ + rb_ * W.template_flux(pt["stars"][0], ccoef)
            f = [v[0][0] for v in arc.values() if 0 in v]
            if T0 > 0 and len(f) >= min_arc:
                ratios.append(float(np.mean(f)) / T0)
                by_grp[qgroup.get(evmeta[key]["encounter"])].append(ratios[-1])
        k_unit[ukey] = float(np.median(ratios)) if ratios else 1.0
        # per q-group scale when the group has >= 3 events (the ratio drifts
        # with the exposure regime), else the unit value
        k_group[ukey] = {g: (float(np.median(v)) if len(v) >= 3 else k_unit[ukey]) for g, v in by_grp.items()}
    ctrl_idx = list(range(1, 9))
    min_ctrl = CFG["controls"]["min_controls_per_epoch"]
    min_arc = CFG["gates"]["event"]["min_usable_arc_epochs"]
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
            T = {pi: ra_ + rb_ * W.template_flux(pt["stars"][pi], ccoef) for pi in range(9)}
            if ukey in bright_units:                    # v1.3/1.4 self-calibrated template, per q-group
                T[0] = T[0] * k_group.get(ukey, {}).get(qgroup.get(evmeta[key]["encounter"]), k_unit.get(ukey, 1.0))
            D = {pi: [] for pi in range(9)}
            mj_used = []
            for m in sorted(arc):
                v = arc[m]
                E = {pi: v[pi][0] - T[pi] for pi in range(9) if pi in v}
                cs = [pi for pi in ctrl_idx if pi in E]
                if 0 not in E or len(cs) < min_ctrl:
                    continue
                mj_used.append(m)
                D[0].append(E[0] - float(np.median([E[c] for c in cs])))
                for c in cs:
                    D[c].append(E[c] - float(np.median([E[o] for o in cs if o != c])))
            if len(D[0]) < min_arc:
                ev_detail[eid] = {"status": "insufficient_after_veto", "n_arc": len(D[0])}
                continue
            series_dump[key] = {"mjd": [round(m, 5) for m in mj_used],
                                "D": {str(pi): [round(x, 5) for x in D[pi]] for pi in range(9) if len(D[pi]) == len(mj_used)}}
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
                pp["pulse"].append(float(np.max(np.minimum(x[:-1], x[1:]))))   # v1.5: two consecutive frames
                pp["ev"].append(eid)
            d0 = np.array(D[0])
            det = {"status": "included", "n_arc": len(d0), "n_base": len(base),
                   "z": round(float(np.mean(d0) / (_robust_sd(d0) / np.sqrt(len(d0)))), 2),
                   "mean_D": round(float(np.mean(d0)), 5), "sd_D": round(_robust_sd(d0), 5),
                   "template_src": round(T[0], 5), "encounter": evmeta[key]["encounter"],
                   "mjd_ca": round(evmeta[key]["mjd_ca"], 3),
                   "pulse_mjd": round(float(mj_used[int(np.argmax(np.minimum(d0[:-1], d0[1:])))]), 4)}
            # secondary: star-fixed differential where a baseline exists
            if len(base) >= 30:
                fb = [b[0][0] for b in base.values() if 0 in b]
                fa = [arc[m][0][0] for m in mj_used]
                if len(fb) >= 30:
                    det["star_fixed_excess"] = round(float(np.mean(fa) - np.median(fb)), 5)
                    det["star_fixed_z"] = round(float((np.mean(fa) - np.median(fb)) / (_robust_sd(fa) / np.sqrt(len(fa)) + 1e-30)), 2)
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
        res = {"status": "searched", "per_patch": {str(k): v for k, v in stats.items()}, "events": ev_detail}
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
           "epoch_vetoes": dict(n_veto), "bright_body_diagnostic": bright_diag, "colour_system": colour,
           "bright_star_class": {"self_calibrated": {u: round(k, 4) for u, k in k_unit.items()},
                                 "per_q_group": {u: {g: round(k, 4) for g, k in d_.items()} for u, d_ in k_group.items()},
                                 "unsearchable": sorted(unsearchable)},
           "template_response": template_response,
           "stamp_response": {"n": len(stamp_rows),
                              "median_ratio": float(np.median([s["ratio"] for s in stamp_rows])) if stamp_rows else None,
                              "mad": _robust_sd([s["ratio"] for s in stamp_rows]) if stamp_rows else None},
           "trials": {"searched": n_trials, "exceedances": n_exc, "expected_control_crossings": round(n_trials / 9, 2)},
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


if __name__ == "__main__":
    cmd, role = sys.argv[1], sys.argv[2]
    if cmd == "index":
        build_index(role)
    elif cmd == "run":
        run(role, workers=int(sys.argv[3]) if len(sys.argv) > 3 else 16)
    elif cmd == "reduce":
        reduce(role, sys.argv[3] if len(sys.argv) > 3 else "v1")
