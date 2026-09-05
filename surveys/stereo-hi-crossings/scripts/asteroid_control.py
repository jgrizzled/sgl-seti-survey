"""Bright-asteroid arc-passage positive control (freeze D5).

Stage `find`: Horizons OBSERVER ephemerides of (4) Vesta, (1) Ceres,
(2) Pallas with STEREO-A as the observer (CENTER='@-234'), daily
2007-2026 -> helioprojective track -> days inside the frozen arc band
(4.05 <= |HPLN| <= 6.0 deg, inside the era footprint). Picks up to
three passages (the faintest predicted V <= 10 preferred, unsaturated)
and snapshots the responses.

Stage `run`: fetches the passage frames (RAL L1), interpolates the
Horizons track (10-min ephemeris) to each frame's DATE-AVG, measures
the moving source through the identical chain (verified WCS, 2-D ZP,
top-hat aperture, ZP_REF normalization), and compares V_eq =
ZP_REF - 2.5 log10(F_norm) against the Horizons V corrected by the
frozen colour term (V_inst = V - 0.595 (B-V - 0.65); asteroid B-V:
Vesta 0.78, Ceres 0.72, Pallas 0.66). Also measures the same frames'
8 parallel-track control patches around the asteroid position for a
window-statistic recovery (the D construction with a moving source).

Output: results/asteroid_control_v1.json; frames purged after use.
"""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hi_geometry as G
import hi_lib as H
import series as S

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "stereo-hi-crossings"
OUT = SURV / "results" / "asteroid_control_v1.json"
RUNS = REPO / "runs" / "stereo-hi-crossings" / "asteroid_control"
BODIES = {"vesta": ("4;", 0.78), "ceres": ("1;", 0.72), "pallas": ("2;", 0.66)}


def horizons(command: str, start: str, stop: str, step: str) -> str:
    q = {"format": "text", "COMMAND": f"'{command}'", "OBJ_DATA": "'NO'", "MAKE_EPHEM": "'YES'",
         "EPHEM_TYPE": "'OBSERVER'", "CENTER": "'@-234'", "START_TIME": f"'{start}'",
         "STOP_TIME": f"'{stop}'", "STEP_SIZE": f"'{step}'", "QUANTITIES": "'1,9,20'",
         "CSV_FORMAT": "'YES'", "ANG_FORMAT": "'DEG'"}
    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" + urllib.parse.urlencode(q)
    return urllib.request.urlopen(url, timeout=300).read().decode()


def parse(txt: str):
    rows = []
    for line in txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
        p = [x.strip() for x in line.split(",")]
        if len(p) < 7:
            continue
        try:
            t = _ptime(p[0])
            rows.append((t.mjd, float(p[3]), float(p[4]), float(p[5])))
        except (ValueError, KeyError):
            continue
    return np.array(rows)


def _ptime(s: str) -> Time:
    mon = {m: i + 1 for i, m in enumerate(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                                            "Aug", "Sep", "Oct", "Nov", "Dec"])}
    d, hm = s.split()
    y, m, dd = d.split("-")
    return Time(f"{y}-{mon[m]:02d}-{dd}T{hm}:00", scale="utc")


def find() -> dict:
    RUNS.mkdir(parents=True, exist_ok=True)
    cands = []
    for name, (cmd, bv) in BODIES.items():
        txt = horizons(cmd, "2007-01-01", "2026-09-01", "1d")
        (RUNS / f"{name}_daily_horizons.txt").write_text(txt)
        eph = parse(txt)
        v = G.radec_to_vec(eph[:, 1], eph[:, 2])
        ln, lt = zip(*[G.helioprojective(vv, m) for vv, m in zip(v, eph[:, 0])])
        ln, lt = np.array(ln).ravel(), np.array(lt).ravel()
        infov = G.in_fov_era(ln, lt, eph[:, 0])
        arc = infov & (np.abs(ln) >= S.ARC[0]) & (np.abs(ln) <= S.ARC[1])
        # group consecutive days
        idx = np.nonzero(arc)[0]
        groups = np.split(idx, np.nonzero(np.diff(idx) > 2)[0] + 1) if len(idx) else []
        for g in groups:
            if len(g) == 0:
                continue
            cands.append({"body": name, "bv": bv, "start": Time(eph[g[0], 0] - 0.5, format="mjd").iso[:10],
                          "stop": Time(eph[g[-1], 0] + 1.5, format="mjd").iso[:10],
                          "v_pred_median": float(np.median(eph[g, 3])), "days": int(len(g)),
                          "hplt": float(np.median(lt[g]))})
    # prefer faint (unsaturated) passages V >= 6.5; take up to three, spread over eras
    ok = sorted([c for c in cands if c["v_pred_median"] >= 6.5], key=lambda c: -c["v_pred_median"])
    chosen = ok[:3]
    res = {"candidates": cands, "chosen": chosen}
    (RUNS / "passages.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))
    return res


def run(workers: int = 12) -> None:
    ps = json.loads((RUNS / "passages.json").read_text())["chosen"]
    results = []
    for c in ps:
        cmd, bv = BODIES[c["body"]]
        txt = horizons(cmd, c["start"], c["stop"], "10m")
        (RUNS / f'{c["body"]}_{c["start"]}_10m.txt').write_text(txt)
        eph = parse(txt)
        # frames in the span
        d0, d1 = int(Time(c["start"]).mjd), int(Time(c["stop"]).mjd)
        stems = []
        for d in range(d0, d1 + 1):
            day = Time(d, format="mjd").strftime("%Y%m%d")
            stems += [s for s, m in S.load_listing(day)]
        rows = []
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(6) as tp:
            fetched = list(tp.map(lambda s: (s, *H.fetch_frame(s[:8], s, dest_dir=RUNS / "frames", prefer=S.PREFER, fallback=False)), stems))
        for s, p, src in fetched:
            if p is None:
                continue
            try:
                from astropy.io import fits
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    h = fits.getheader(p)
                if H.frame_gate(h):
                    continue
                fr = H.load_frame(p)
                if not H.fit_frame(fr) or fr.zp_coef is None or fr.zp_scatter > 0.12:
                    continue
                ra = np.interp(fr.mjd_avg, eph[:, 0], eph[:, 1])
                dec = np.interp(fr.mjd_avg, eph[:, 0], eph[:, 2])
                vmag = np.interp(fr.mjd_avg, eph[:, 0], eph[:, 3])
                s_hat = G.radec_to_vec(ra, dec)
                ln, lt = G.helioprojective(s_hat, fr.mjd_avg)
                m = H.measure_patch(fr, ra, dec)
                if not H.in_footprint(m["x"], m["y"], S.MARGIN_PX) or m["valid"] < 0.9 or not np.isfinite(m["flux"]):
                    continue
                # 8 parallel-track controls at the same epoch
                ctrl = []
                for off in S.OFFS:
                    vv = G.from_helioprojective(ln[0], lt[0] + off, fr.mjd_avg)
                    cra, cdec = G.vec_to_radec(vv)
                    mc = H.measure_patch(fr, float(cra), float(cdec))
                    if H.in_footprint(mc["x"], mc["y"], S.MARGIN_PX) and np.isfinite(mc["flux"]):
                        ctrl.append(mc["flux"])
                v_inst = vmag - H.COLOUR_COEFF * (bv - H.BV_REF)
                v_eq = H.ZP_REF - 2.5 * np.log10(m["flux"]) if m["flux"] > 0 else np.nan
                rows.append({"mjd": fr.mjd_avg, "hpln": float(ln[0]), "hplt": float(lt[0]),
                             "v_horizons": float(vmag), "v_inst_pred": float(v_inst),
                             "v_eq_measured": float(v_eq), "flux": m["flux"], "err": m["err"],
                             "ctrl_median": float(np.median(ctrl)) if ctrl else None,
                             "in_arc": bool(S.ARC[0] <= abs(ln[0]) <= S.ARC[1])})
            except Exception as exc:
                rows.append({"error": str(exc)[:100]})
            finally:
                try:
                    p.unlink()
                except Exception:
                    pass
        good = [r for r in rows if "error" not in r and np.isfinite(r["v_eq_measured"])]
        arc = [r for r in good if r["in_arc"]]
        d = np.array([r["v_eq_measured"] - r["v_inst_pred"] for r in good])
        da = np.array([r["v_eq_measured"] - r["v_inst_pred"] for r in arc])
        results.append({**c, "n_frames": len(rows), "n_measured": len(good), "n_in_arc": len(arc),
                        "offset_mag_all": {"median": float(np.median(d)) if len(d) else None,
                                           "mad": S._robust_sd(d) if len(d) else None},
                        "offset_mag_arc": {"median": float(np.median(da)) if len(da) else None,
                                           "mad": S._robust_sd(da) if len(da) else None},
                        "flux_over_ctrl_median_arc": float(np.median([r["flux"] / max(abs(r["ctrl_median"]), 1e-3)
                                                                       for r in arc if r["ctrl_median"] is not None])) if arc else None,
                        "rows": rows})
        print(c["body"], c["start"], "measured", len(good), "arc", len(arc),
              "offset all", results[-1]["offset_mag_all"], "arc", results[-1]["offset_mag_arc"])
    OUT.write_text(json.dumps({"stage": "asteroid_control_v1", "colour_coeff": H.COLOUR_COEFF,
                               "gate_mag": 0.2, "passages": results}, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    if sys.argv[1] == "find":
        find()
    else:
        run()
