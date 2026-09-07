"""Positive control (freeze D5): a bright moving body through the SoloHI
inner tiles, measured by the identical patch chain along its Horizons
track (CENTER='@-144').

`scan`: over the covered era, query Horizons (2-h steps) for the bright
numbered asteroids and record inner-tile passages with V <= V_MAX
(planning tile model). `measure <n> <label> [every]`: fetch the tile's
L2 frames of the passage (every k-th), run the chain, measure the
aperture flux at the interpolated position and at 8 orbit-latitude
control positions, compare V_eq (B-V 0.75) with the Horizons V.
Output: results/moving_body_scan_v1.json, results/moving_body_control_v1.json
"""
from __future__ import annotations
import json, sys, time, urllib.parse, urllib.request
from datetime import datetime
from pathlib import Path
import numpy as np
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import solo_geometry as G
import solohi_lib as W

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "solohi-crossings"
RUNS = REPO / "runs" / "solohi-crossings"
ASTEROIDS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 27, 29, 39, 40, 41, 42, 43, 44, 51, 52, 63, 64, 68, 88, 89, 115, 192, 216, 324, 349, 354, 387, 471, 511, 532, 704]
V_MAX = 9.0
BV_AST = 0.75
ERA = ("2022-01-01", "2026-04-10")


def horizons_track(body: str, t0: str, t1: str, step: str = "2h"):
    q = {"format": "text", "COMMAND": f"'{body}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'@-144'", "START_TIME": f"'{t0}'",
         "STOP_TIME": f"'{t1}'", "STEP_SIZE": f"'{step}'", "QUANTITIES": "'1,9'",
         "CSV_FORMAT": "'YES'", "ANG_FORMAT": "'DEG'", "TIME_DIGITS": "'MINUTES'"}
    for attempt in range(3):
        try:
            txt = urllib.request.urlopen("https://ssd.jpl.nasa.gov/api/horizons.api?"
                                         + urllib.parse.urlencode(q), timeout=180).read().decode()
            if "$$SOE" not in txt:
                return None
            rows = []
            for line in txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
                p = [x.strip() for x in line.split(",")]
                try:
                    rows.append((Time(datetime.strptime(p[0][:17], "%Y-%b-%d %H:%M")).mjd,
                                 float(p[3]), float(p[4]), float(p[5]) if p[5] not in ("n.a.", "") else np.nan))
                except Exception:
                    pass
            return np.array(rows)
        except Exception:
            time.sleep(5)
    return None


def scan():
    out = []
    for n in ASTEROIDS:
        tr = horizons_track(f"{n};", ERA[0], ERA[1])
        if tr is None or len(tr) == 0:
            continue
        m, ra, dec, v = tr.T
        if np.nanmin(v) > V_MAX:
            continue
        s = G.radec_to_vec(ra, dec)
        _, lam, bet = G.solohi_coords(s, m) if False else (None, None, None)
        u, t, nrm = G.ram_frame(m)
        lam = np.degrees(np.arctan2(np.sum(s * t, -1), np.sum(s * u, -1)))
        bet = np.degrees(np.arcsin(np.clip(np.sum(s * nrm, -1), -1, 1)))
        tiles = G.tile_of(lam, bet)
        inf = np.isin(tiles, ["1", "2"]) & (v <= V_MAX)
        if not inf.any():
            continue
        # split into passages (gaps > 2 d)
        idx = np.nonzero(inf)[0]
        breaks = np.nonzero(np.diff(m[idx]) > 2.0)[0]
        for seg in np.split(idx, breaks + 1):
            out.append({"asteroid": n, "label": Time(m[seg[0]], format="mjd").iso[:10],
                        "tile": str(tiles[seg[0]]), "hours_in_field": int(len(seg) * 2),
                        "V_range": [round(float(np.nanmin(v[seg])), 2), round(float(np.nanmax(v[seg])), 2)],
                        "mjd_range": [float(m[seg].min()), float(m[seg].max())],
                        "utc_range": [Time(m[seg].min(), format="mjd").iso[:16], Time(m[seg].max(), format="mjd").iso[:16]],
                        "eps_range": [round(float(np.degrees(np.arccos(np.sum(s * u, -1)))[seg].min()), 1),
                                      round(float(np.degrees(np.arccos(np.sum(s * u, -1)))[seg].max()), 1)]})
            print(out[-1], flush=True)
    (SURV / "results" / "moving_body_scan_v1.json").write_text(json.dumps(out, indent=1) + "\n")
    print("wrote scan:", len(out), "passages")


def measure(n: int, label: str, every: int = 4):
    scan_rows = json.loads((SURV / "results" / "moving_body_scan_v1.json").read_text())
    row = next(r for r in scan_rows if r["asteroid"] == n and r["label"] == label)
    m0, m1 = row["mjd_range"]
    tr = horizons_track(f"{n};", Time(m0 - 0.1, format="mjd").iso[:16], Time(m1 + 0.1, format="mjd").iso[:16], "10m")
    inv = json.loads((RUNS / "coverage" / "soar_l2_inventory.json").read_text())
    ci = {c: i for i, c in enumerate(inv["columns"])}
    frames = []
    for r in inv["rows"]:
        if r[ci["descriptor"]] in ("solohi-1ft", "solohi-2ft"):
            mj = Time(r[ci["begin_time"]]).mjd
            if m0 <= mj <= m1:
                frames.append((mj, r[ci["begin_time"]][:10].replace("-", ""), r[ci["filename"]]))
    frames.sort()
    frames = frames[::every]
    print(f"asteroid {n} {label}: {len(frames)} frames to measure", flush=True)
    recs = []
    for mj, day, fname in frames:
        p, src = W.fetch_frame(day, fname)
        if p is None:
            continue
        try:
            fr = W.load_frame(p)
            ok = W.fit_frame(fr)
            if not ok or fr.zp_coef is None or fr.zp_scatter / np.sqrt(max(fr.n_calib, 1)) > W.MAX_ZP_UNC:
                continue
            ra = np.interp(fr.mjd_avg, tr[:, 0], tr[:, 1])
            dec = np.interp(fr.mjd_avg, tr[:, 0], tr[:, 2])
            v = np.interp(fr.mjd_avg, tr[:, 0], tr[:, 3])
            m = W.measure_patch(fr, ra, dec)
            if not W.in_footprint(fr, m["x"], m["y"], 12) or not np.isfinite(m["flux"]) or m["flux"] <= 0:
                continue
            s_hat = G.radec_to_vec(ra, dec)
            _, lam, bet = G.solohi_coords(s_hat, fr.mjd_avg)
            ctrl = []
            for off in (1.5, -1.5, 3.0, -3.0, 4.5, -4.5, 6.0, -6.0):
                u, t, nrm = G.ram_frame(fr.mjd_avg)
                ln, lt = np.radians(lam[0]), np.radians(bet[0] + off)
                vv = np.cos(lt) * np.cos(ln) * u[0] + np.cos(lt) * np.sin(ln) * t[0] + np.sin(lt) * nrm[0]
                cra, cdec = G.vec_to_radec(vv)
                c = W.measure_patch(fr, float(cra), float(cdec))
                if W.in_footprint(fr, c["x"], c["y"], 12) and np.isfinite(c["flux"]):
                    ctrl.append(c["flux"])
            v_eq = W.ZP_REF - 2.5 * np.log10(m["flux"])
            v_pred = v - W.COLOUR_COEFF * (BV_AST - W.BV_REF)
            recs.append({"mjd": fr.mjd_avg, "utc": Time(fr.mjd_avg, format="mjd").iso[:16], "x": m["x"], "y": m["y"],
                         "tile": fr.tile, "flux": m["flux"], "err": m["err"], "V_horizons": round(float(v), 2),
                         "V_pred_inst": round(float(v_pred), 2), "V_eq": round(float(v_eq), 2),
                         "delta_mag": round(float(v_eq - v_pred), 3),
                         "ctrl_flux_median": float(np.median(ctrl)) if ctrl else None, "n_ctrl": len(ctrl),
                         "r_au": fr.r_au, "xposure": fr.xposure, "zp0": float(fr.zp_coef[0])})
        finally:
            try:
                p.unlink()
            except OSError:
                pass
    d = np.array([r["delta_mag"] for r in recs])
    summ = {"asteroid": n, "label": label, "n_frames": len(recs),
            "delta_mag_median": round(float(np.median(d)), 3) if len(d) else None,
            "delta_mag_mad": round(float(1.4826 * np.median(np.abs(d - np.median(d)))), 3) if len(d) else None,
            "records": recs}
    out_p = SURV / "results" / "moving_body_control_v1.json"
    allr = json.loads(out_p.read_text()) if out_p.exists() else []
    allr = [r for r in allr if not (r["asteroid"] == n and r["label"] == label)] + [summ]
    out_p.write_text(json.dumps(allr, indent=1) + "\n")
    print({k: v for k, v in summ.items() if k != "records"})


if __name__ == "__main__":
    if sys.argv[1] == "scan":
        scan()
    else:
        measure(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 4)
