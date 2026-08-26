"""Dev search driver (threshold freeze v1.0): series + statistics for
the 5 dev units.

Design: one pass over fetched frames. A task index maps each frame to
the (unit, event, kind) photometry it serves; a worker loads + fits
the frame once and measures every mapped position (source + the 8
frozen PA-ring controls). Records stream to a jsonl; the reduce step
assembles per-(unit, event, dpa) series, applies the frozen detrend
(30 d baseline running median + variance rescale k), and computes
S_stack / S_event / S_pulse with T = max over ring controls.

Positions are solar-frame: in-window epochs at the event's
source_track (r(t), PA(t)); baseline epochs at the event's median
visible (r, PA). Calibrator records (V, B-V, r, flux, exptime) are
dumped from every fitted frame for the L2 colour-coefficient
measurement. C2 frames that fail the 5-star gate are still measured
when >= 2 stars match (finding L1); every epoch carries fit_quality
so the L1/L2 amendment can gate after the fact.

Stages:
  python dev_series.py index    # build the frame->task index
  python dev_series.py measure  # per-frame pass (resumable, parallel)
  python dev_series.py reduce   # series -> statistics -> results
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lasco_lib as L
import tracks

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "heliospheric-crossings"
CFG = json.loads((SURV / "configs" / "threshold_freeze_v1.json").read_text())
OCC = json.loads((SURV / "results" / "occulter_radii_v1.json").read_text())["adopted"]
FRAMES = REPO / "runs" / "heliospheric-crossings" / "frames"
DEVDIR = REPO / "runs" / "heliospheric-crossings" / "dev"
EVENTS = REPO / "crossings" / "soho_v1" / "events.ecsv"

RSUN_AU = 0.00465047
AU_KM = 1.495978707e8
BASELINE_D = 110
RING_DPA = CFG["controls"]["ring_dpa_deg"]
ANN = {"c2": OCC["c2_rsun"], "c3": OCC["c3_rsun"]}
RUNG_R = {"S1_2.5Rsun": 2.5 * RSUN_AU, "S1_0.1AU": 0.1, "S2_0.1AU": 0.1}
COMBOS = {"S1": ("outbound", "target"), "S2": ("inbound", "anti_target")}


def dev_events(role: str = "dev"):
    """{unit_key: [event dicts]} for the given role's units."""
    t = Table.read(EVENTS)
    out = {}
    for u in CFG["units"]:
        if u["role"] != role:
            continue
        ch = u["unit"].split("_")[0]
        r_au = RUNG_R[u["unit"]]
        direction, side = COMBOS[ch]
        m = ((np.asarray(t["link_direction"]) == direction)
             & (np.asarray(t["side"]) == side)
             & (np.asarray(t["target_id"]) == u["target_id"])
             & (np.asarray(t["b_min_au"]) <= r_au))
        evs = []
        for ev in t[m]:
            b, v = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
            half = float(np.sqrt(max(r_au**2 - b**2, 0)) * AU_KM / (v * 86400.0))
            evs.append({
                "event_id": str(ev["event_id"]),
                "star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
                "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"]),
                "mjd_ca": Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd,
                "half_days": half, "b_min_au": b,
            })
        out[f'{u["target_id"]}|{u["unit"]}'] = {"camera": u["camera"], "events": evs}
    return out


def day_str(mjd: float) -> str:
    return Time(float(int(mjd)), format="mjd").strftime("%y%m%d")


def build_index(role: str = "dev") -> dict:
    """frame relpath -> list of tasks {unit, event_id, kind}."""
    units = dev_events(role)
    index: dict[str, list] = defaultdict(list)
    meta = {}
    for ukey, u in units.items():
        cam = u["camera"]
        for ev in u["events"]:
            # median visible (r, PA) for the baseline position
            mjds = ev["mjd_ca"] + np.linspace(-ev["half_days"], ev["half_days"], 25)
            rs, pas = [], []
            for m in mjds:
                s = tracks.source_track(ev, m)
                if ANN[cam][0] <= s["r_rsun"] <= ANN[cam][1]:
                    rs.append(s["r_rsun"])
            if not rs:
                meta[f'{ukey}|{ev["event_id"]}'] = {"visible": False}
                continue
            meta[f'{ukey}|{ev["event_id"]}'] = {"visible": True,
                                                "r_med": float(np.median(rs))}
            for d in range(int(ev["mjd_ca"] - ev["half_days"]),
                           int(ev["mjd_ca"] + ev["half_days"]) + 1):
                dday = FRAMES / cam / day_str(d)
                if dday.is_dir():
                    for f in sorted(dday.iterdir()):
                        index[str(f.relative_to(REPO))].append(
                            {"u": ukey, "e": ev["event_id"], "k": "win"})
            for d in range(int(ev["mjd_ca"] - BASELINE_D),
                           int(ev["mjd_ca"] + BASELINE_D) + 1):
                dday = FRAMES / cam / day_str(d)
                if dday.is_dir():
                    files = sorted(dday.iterdir())
                    if files:
                        index[str(files[0].relative_to(REPO))].append(
                            {"u": ukey, "e": ev["event_id"], "k": "base"})
    DEVDIR.mkdir(parents=True, exist_ok=True)
    (DEVDIR / f"task_index_{role}_v1.json").write_text(json.dumps(
        {"index": index, "event_meta": meta,
         "units": {k: v for k, v in units.items()}}, indent=0) + "\n")
    n_win = sum(1 for v in index.values() if any(t["k"] == "win" for t in v))
    print(f"{len(index)} frames indexed ({n_win} window frames); "
          f"{len(meta)} unit-events")
    return index


def measure_frame(args):
    relpath, tasks, units_events = args
    out = {"frame": relpath, "records": [], "calibrators": []}
    try:
        reason = L.frame_usable(REPO / relpath)
        if reason:
            out["unusable"] = reason
            return out
        fr = L.load_frame(REPO / relpath)
        fit_ok = L.fit_frame(fr)
        out["fit"] = {"ok": bool(fit_ok), "n_match": fr.n_match,
                      "zp_scatter": float(fr.zp_scatter),
                      "rot": fr.rot_deg, "mjd": fr.mjd_mid}
        if fr.rot_deg is None:
            return out  # no astrometry at all -> epoch lost (recorded)
        # calibrator dump for the L2 ensemble
        if fr.zp_r is not None:
            out["calibrators"] = [{"r": float(r), "zp": float(z)}
                                  for r, z in zip(fr.zp_r, fr.zp_v)]
        for tk in tasks:
            ev = units_events[tk["u"]]["events_by_id"][tk["e"]]
            emeta = units_events["_meta"].get(f'{tk["u"]}|{tk["e"]}', {})
            if not emeta.get("visible"):
                continue
            for dpa in [0.0] + [float(x) for x in RING_DPA]:
                if tk["k"] == "win":
                    s = tracks.source_track(ev, fr.mjd_mid, dpa_deg=dpa)
                    r_rsun = s["r_rsun"]
                    cam = units_events[tk["u"]]["camera"]
                    if not (ANN[cam][0] <= r_rsun <= ANN[cam][1]):
                        continue
                    ra, dec = s["ra"], s["dec"]
                else:
                    # baseline: the event's median visible solar-frame
                    # radius at the *current* transverse azimuth + dpa —
                    # samples the same coronal radius year-round
                    s = tracks.source_track(ev, fr.mjd_mid, dpa_deg=dpa)
                    scale = emeta["r_med"] / max(s["r_rsun"], 1e-9)
                    # rescale the offset to the median radius
                    ra = fr.sun_ra + (s["ra"] - fr.sun_ra) * scale
                    dec = fr.sun_dec + (s["dec"] - fr.sun_dec) * scale
                p = L.forced_photometry(fr, ra, dec)
                if fr.exptime > 0:   # exposure-normalized flux units
                    p["flux"] /= fr.exptime
                    p["err"] /= fr.exptime
                p.update({"u": tk["u"], "e": tk["e"], "k": tk["k"], "dpa": dpa})
                for drop in ("x", "y", "zp_scatter", "n_match"):
                    p.pop(drop, None)
                out["records"].append(p)
    except Exception as exc:
        out["error"] = str(exc)[:200]
    return out


def run_measure(workers: int = 6, role: str = "dev") -> None:
    idx = json.loads((DEVDIR / f"task_index_{role}_v1.json").read_text())
    units = idx["units"]
    for u in units.values():
        u["events_by_id"] = {e["event_id"]: e for e in u["events"]}
    units["_meta"] = idx["event_meta"]
    done = set()
    out_path = DEVDIR / f"measurements_{role}_v1.jsonl"
    if out_path.exists():
        with open(out_path) as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["frame"])
                except Exception:
                    pass
    todo = [(rp, tasks, units) for rp, tasks in idx["index"].items()
            if rp not in done and (REPO / rp).exists()]
    print(f"{len(todo)} frames to measure ({len(done)} already done)")
    with open(out_path, "a") as fh, ProcessPoolExecutor(workers) as ex:
        for i, rec in enumerate(ex.map(measure_frame, todo, chunksize=8)):
            fh.write(json.dumps(rec) + "\n")
            if (i + 1) % 500 == 0:
                fh.flush()
                print(f"  {i+1}/{len(todo)}", flush=True)


def reduce_series() -> None:
    idx = json.loads((DEVDIR / "task_index_v1.json").read_text())
    series = defaultdict(list)   # (u, e, dpa) -> [(mjd, flux, err, k(ind))]
    fit_stats = defaultdict(list)
    with open(DEVDIR / "measurements_v1.jsonl") as fh:
        for line in fh:
            rec = json.loads(line)
            if "fit" in rec:
                cam = "c2" if "/c2/" in rec["frame"] else "c3"
                fit_stats[cam].append((rec["fit"]["n_match"], rec["fit"]["zp_scatter"],
                                       rec["fit"]["ok"]))
            for p in rec.get("records", []):
                if p.get("valid_fraction", 0) < 0.9 or not np.isfinite(p.get("flux", np.nan)):
                    continue
                # degenerate error = photometry annulus in a constant
                # fill region (occulted core plateau) — physically
                # invalid epoch, and 1/err^2 weights explode
                if not np.isfinite(p.get("err", np.nan)) or p["err"] <= 1e-6:
                    continue
                series[(p["u"], p["e"], p["dpa"])].append(
                    (p["mjd"], p["flux"], p["err"], p["k"]))

    results = {}
    for ukey in idx["units"]:
        cam = idx["units"][ukey]["camera"]
        stats = {}
        for dpa in [0.0] + [float(x) for x in RING_DPA]:
            zs, pulses = [], []
            for ev in idx["units"][ukey]["events"]:
                rows = sorted(series.get((ukey, ev["event_id"], dpa), []))
                win = [r for r in rows if r[3] == "win"]
                base = [r for r in rows if r[3] == "base"]
                if len(win) < CFG["gates"]["min_visible_epochs_per_event"] \
                        or len(base) < 30:
                    continue
                # relative error floor: epochs with errors far below the
                # series norm are residual degenerate-annulus cases
                med_err = np.median([r[2] for r in win + base])
                win = [r for r in win if r[2] >= 0.05 * med_err]
                base = [r for r in base if r[2] >= 0.05 * med_err]
                if len(win) < CFG["gates"]["min_visible_epochs_per_event"] \
                        or len(base) < 30:
                    continue
                bm, bf = np.array([r[0] for r in base]), np.array([r[1] for r in base])
                # frozen detrend: 30 d running median of the baseline
                def detr(mjds, fluxes):
                    med = np.array([np.median(bf[np.abs(bm - m) <= 15]) if
                                    (np.abs(bm - m) <= 15).any() else np.median(bf)
                                    for m in mjds])
                    return fluxes - med
                wr = detr(np.array([r[0] for r in win]), np.array([r[1] for r in win]))
                we = np.array([r[2] for r in win])
                br = detr(bm, bf)
                be = np.array([r[2] for r in base])
                k = max(1.0, float(np.median((br / be) ** 2) / 0.4549))
                we_k = we * np.sqrt(k)
                x = float(np.sum(wr / we_k**2) / np.sum(1 / we_k**2))
                sx = float(1 / np.sqrt(np.sum(1 / we_k**2)))
                zs.append((x / sx, sx))
                pulses.append(float(np.max(wr / we_k)))
            if len(zs) >= CFG["gates"]["min_events_for_stack"]:
                z = np.array([a for a, _ in zs])
                stats[f"dpa{dpa:+.0f}"] = {
                    "n_events": len(zs),
                    "S_stack": float(np.sum(z) / np.sqrt(len(z))),
                    "S_event": float(np.max(z)),
                    "S_pulse": float(np.max(pulses)) if cam == "c2" else None,
                }
        if not stats:
            results[ukey] = {"status": "constraint_only_or_pending"}
            continue
        src = stats.get("dpa+0", {})
        ctrl = [v for kk, v in stats.items() if kk != "dpa+0"]
        res = {"status": "searched", "per_dpa": stats}
        for s in ("S_stack", "S_event", "S_pulse"):
            if src.get(s) is None:
                continue
            ring_vals = [c[s] for c in ctrl if c.get(s) is not None]
            T = max(ring_vals, default=0.0)
            ring_mad = float(np.median(np.abs(np.array(ring_vals)
                             - np.median(ring_vals))) * 1.4826) if len(ring_vals) >= 3 else None
            res[s] = {"S": src[s], "T": float(T),
                      "exceeds": bool(src[s] > max(T, 0.0)),
                      "S_ring_norm": round(src[s] / ring_mad, 2) if ring_mad else None}
        results[ukey] = res

    summary = {"fit_stats": {cam: {
        "n_frames": len(v),
        "fit_ok_frac": float(np.mean([r[2] for r in v])),
        "median_n_match": float(np.median([r[0] for r in v])),
        "median_zp_scatter": float(np.median([r[1] for r in v if np.isfinite(r[1])])) if v else None,
    } for cam, v in fit_stats.items()}, "units": results}
    (SURV / "results" / "dev_search_v1.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


def _planet_tables(mjds: np.ndarray) -> dict:
    """ICRS (ra, dec) of solar-system contaminants as seen from SOHO,
    vectorized over the unique frame epochs."""
    from astropy.coordinates import get_body_barycentric

    t = Time(mjds, format="mjd", scale="utc")
    soho = np.stack([np.interp(mjds, L._soho["mjd_utc"], L._soho["xyz_au"][:, i])
                     for i in range(3)], axis=1)
    out = {}
    for body in ("mercury", "venus", "mars", "jupiter", "saturn", "moon"):
        xyz = get_body_barycentric(body, t).xyz.to_value("AU").T - soho
        r = np.linalg.norm(xyz, axis=1)
        out[body] = (np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0])) % 360.0,
                     np.degrees(np.arcsin(xyz[:, 2] / r)))
    return out


def reduce_v11(role: str = "dev") -> None:
    """Amendment-v1.1 reduce: (A1) the frozen planet-proximity epoch
    mask, implemented; (A2) differential statistic — source minus the
    mean of the ring controls on the same frame (controls: ring_i minus
    the mean of the other rings); (A3) empirical night-level sigma from
    the differential baseline replacing formal-error z-scores."""
    idx = json.loads((DEVDIR / f"task_index_{role}_v1.json").read_text())
    units = idx["units"]
    for u in units.values():
        u["events_by_id"] = {e["event_id"]: e for e in u["events"]}

    # gather epoch records grouped by (u, e, kind, mjd) across dpa
    epochs = defaultdict(dict)   # key -> {dpa: flux}
    zp_by_key = {}
    fit_stats = defaultdict(lambda: [0, 0])
    with open(DEVDIR / f"measurements_{role}_v1.jsonl") as fh:
        for line in fh:
            rec = json.loads(line)
            if "fit" in rec:
                cam = "c2" if "/c2/" in rec["frame"] else "c3"
                fit_stats[cam][0] += 1
                fit_stats[cam][1] += rec["fit"]["rot"] is not None
            for p in rec.get("records", []):
                if (p.get("valid_fraction", 0) < 0.9
                        or not np.isfinite(p.get("flux", np.nan))
                        or not np.isfinite(p.get("err", np.nan)) or p["err"] <= 1e-6):
                    continue
                key = (p["u"], p["e"], p["k"], round(p["mjd"], 6))
                epochs[key][p["dpa"]] = p["flux"]
                if p["dpa"] == 0.0:
                    zp_by_key[key] = p.get("zp", float("nan"))

    # A1: planet mask per epoch at the source position
    all_keys = list(epochs)
    mjds = np.array(sorted({k[3] for k in all_keys}))
    __import__("tracks").warm_geometry_cache(mjds)
    ptab = _planet_tables(mjds)
    # v1.2 rule 8: bright-planet in-FOV veto (frame-wide stray light)
    sun_radec = {}
    for m in mjds:
        ra, dec, _, _ = L.sun_from_soho(float(m))
        sun_radec[float(m)] = (ra, dec)
    fov_deg = {"c2": 1.7, "c3": 8.5}
    infov = {}
    for i, m in enumerate(mjds):
        sra, sdec = sun_radec[float(m)]
        bad = False
        for body in ("venus", "jupiter"):
            pra, pdec = ptab[body][0][i], ptab[body][1][i]
            dra = (pra - sra + 180) % 360 - 180
            sep = np.hypot(dra * np.cos(np.radians(sdec)), pdec - sdec)
            if sep < 8.5:
                infov[float(m)] = sep
                bad = True
        # store max-fov sep; per-camera cut applied at use
    n_planet_fov = [0]
    n_planet_masked = 0
    keep = {}
    for key in all_keys:
        ukey, eid, kind, mjd = key
        ev = units[ukey]["events_by_id"][eid]
        cam = units[ukey]["camera"]
        s = __import__("tracks").source_track(ev, mjd)
        lim = 10 * (11.9 if cam == "c2" else 56.0) / 3600.0   # 10 px in deg
        i = int(np.searchsorted(mjds, mjd))
        bad = False
        for body, (pra, pdec) in ptab.items():
            dra = (pra[i] - s["ra"] + 180) % 360 - 180
            if np.hypot(dra * np.cos(np.radians(s["dec"])), pdec[i] - s["dec"]) < lim:
                bad = True
                break
        if bad:
            n_planet_masked += 1
        else:
            keep[key] = epochs[key]

    # bright-star proximity mask (3 px): positions recomputed per
    # (epoch, dpa) exactly as measure did; V<=8 Hipparcos via KDTree
    from scipy.spatial import cKDTree
    tyc = np.load(REPO / "runs" / "heliospheric-crossings" / "starcat" / "tyc2_v11.npz")
    bra, bdec = np.radians(tyc["ra"]), np.radians(tyc["dec"])
    bvec = np.stack([np.cos(bdec) * np.cos(bra), np.cos(bdec) * np.sin(bra),
                     np.sin(bdec)], 1)
    star_tree = cKDTree(bvec)

    def star_hit(ra_deg, dec_deg, cam):
        lim = np.radians(3 * (11.9 if cam == "c2" else 56.0) / 3600.0)
        ra, dec = np.radians(ra_deg), np.radians(dec_deg)
        v = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
        d, _ = star_tree.query(v, k=1)
        return d < 2 * np.sin(lim / 2) + 1e-12

    # L3: subtract the target star's expected flux from S2 in-window
    # epochs when the star is bright enough to appear at stack depth
    STAR_V = {"gj-908": (8.98, 1.49), "ross-154": (10.44, 1.70),
              "ross-128": (11.15, 1.75), "van-maanen": (12.37, 0.55),
              "wolf-359": (13.53, 2.03), "teegarden": (15.1, 2.0),
              "gj-1276": (16.0, 2.0)}
    TEMPLATE_VMAX = 12.5
    COLOR_COEFF = json.loads((SURV / "results" / "color_ensemble_v1.json").read_text())         if (SURV / "results" / "color_ensemble_v1.json").exists() else {}
    expmap = json.loads((DEVDIR / "exptime_map_v1.json").read_text())
    exp_by_mjd = {}
    for rp, v in expmap.items():
        if v:
            exp_by_mjd[("c2" if "/c2/" in rp else "c3", round(v[1], 6))] = v[0]

    dpas = [float(x) for x in RING_DPA]
    results = {}
    n_template = [0, 0]   # applied, dropped-for-missing-zp
    for ukey in units:
        cam = units[ukey]["camera"]
        per_dpa_stats = {}
        n_star_masked = [0]
        for probe in [0.0] + dpas:
            zs, pulses, n_pl, ev_ids = [], [], 0, []
            for ev in units[ukey]["events"]:
                win_d, base_d = defaultdict(list), defaultdict(list)
                for (uk, eid, kind, mjd), by_dpa in keep.items():
                    if uk != ukey or eid != ev["event_id"]:
                        continue
                    refs = [by_dpa[d] for d in dpas if d != probe and d in by_dpa]
                    if probe not in by_dpa or len(refs) < 5:
                        continue
                    s_pos = __import__("tracks").source_track(ev, mjd, dpa_deg=probe)
                    # S2 source-star exemption: the in-window source IS
                    # the target star (a Tycho entry for bright targets);
                    # the stellar template models it — masking it would
                    # delete the channel (found on ross-154, v1.3)
                    sep_fov = infov.get(mjd)
                    if sep_fov is not None and sep_fov < fov_deg[cam]:
                        if probe == 0.0:
                            n_planet_fov[0] += 1
                        continue
                    exempt = (probe == 0.0 and kind == "win" and "S2" in ukey)
                    if not exempt and star_hit(s_pos["ra"], s_pos["dec"], cam):
                        n_star_masked[0] += 1
                        continue
                    f_probe = by_dpa[probe]
                    if (probe == 0.0 and kind == "win" and "S2" in ukey
                            and STAR_V[ukey.split("|")[0]][0] <= TEMPLATE_VMAX):
                        vstar, bvstar = STAR_V[ukey.split("|")[0]]
                        zp = zp_by_key.get((uk, eid, kind, mjd), float("nan"))
                        expt = exp_by_mjd.get((cam, mjd), 0.0)
                        if not np.isfinite(zp) or expt <= 0:
                            n_template[1] += 1
                            continue   # cannot predict the star: drop epoch
                        c = COLOR_COEFF.get(cam, {}).get("color_coeff_mag_per_bv", 0.0)
                        m_inst = vstar + c * (bvstar - 0.65)
                        f_probe -= 10 ** ((zp - m_inst) / 2.5) / expt
                        n_template[0] += 1
                    diff = f_probe - float(np.median(refs))
                    (win_d if kind == "win" else base_d)[int(mjd)].append((mjd, diff))
                # night-level series (medians: star transits through
                # the aperture make single epochs heavy-tailed)
                wn = [np.median([d for _, d in v]) for v in win_d.values()]
                bn = [np.median([d for _, d in v]) for v in base_d.values()]
                if len(wn) < 1 or len(bn) < 20:
                    continue
                sd = float(np.std(bn, ddof=1))
                if sd <= 0:
                    continue
                # centre the window differential on the baseline
                # differential (cancels the stationary F-corona PA bias)
                b0 = float(np.median(bn))
                zs.append(float((np.mean(wn) - b0) / (sd / np.sqrt(len(wn)))))
                ev_ids.append((ev["event_id"], ev["mjd_ca"]))
                flat = [d for v in win_d.values() for _, d in v]
                bflat = [d for v in base_d.values() for _, d in v]
                sd_e = float(np.std(bflat, ddof=1))
                pulses.append(float((np.max(flat) - np.median(bflat)) / sd_e)
                              if sd_e > 0 else np.nan)
            if len(zs) >= CFG["gates"]["min_events_for_stack"]:
                z = np.array(zs)
                per_dpa_stats[f"dpa{probe:+.0f}"] = {
                    "n_events": len(zs),
                    "z_by_event": {eid: round(zz, 2) for (eid, _), zz
                                   in zip(ev_ids, zs)} if probe == 0.0 else None,
                    "S_stack": float(np.sum(z) / np.sqrt(len(z))),
                    "S_event": float(np.max(z)),
                    "S_pulse": float(np.nanmax(pulses)) if cam == "c2" and pulses else None,
                    "z_median": float(np.median(z)), "z_mad": float(
                        np.median(np.abs(z - np.median(z))) * 1.4826),
                }
        if not per_dpa_stats or "dpa+0" not in per_dpa_stats:
            results[ukey] = {"status": "constraint_only_or_pending"}
            continue
        src = per_dpa_stats["dpa+0"]
        ctrl = [v for kk, v in per_dpa_stats.items() if kk != "dpa+0"]
        res = {"status": "searched", "per_dpa": per_dpa_stats}
        for s in ("S_stack", "S_event", "S_pulse"):
            if src.get(s) is None:
                continue
            ring_vals = [c[s] for c in ctrl if c.get(s) is not None]
            T = max(ring_vals, default=0.0)
            ring_mad = float(np.median(np.abs(np.array(ring_vals)
                             - np.median(ring_vals))) * 1.4826) if len(ring_vals) >= 3 else None
            res[s] = {"S": src[s], "T": float(T),
                      "exceeds": bool(src[s] > max(T, 0.0)),
                      "S_ring_norm": round(src[s] / ring_mad, 2) if ring_mad else None}
        results[ukey] = res

    n_trials = sum(1 for r in results.values() if r.get("status") == "searched"
                   for st in ("S_stack", "S_event", "S_pulse") if st in r)
    n_exc = sum(1 for r in results.values() if r.get("status") == "searched"
                for st in ("S_stack", "S_event", "S_pulse")
                if st in r and r[st]["exceeds"])
    out = {"construction": "v1.3 amended (A1 planet mask; A2 ring-differential; "
                           "A3 empirical night sigma; Tycho-2 V<=11 star mask; "
                           "L3 stellar template on bright-S2; ring-normalized reporting)",
           "trials": {"searched": n_trials, "exceedances": n_exc,
                      "expected_control_crossings": round(n_trials / 9, 2)},
           "stellar_template": {"applied_epochs": n_template[0],
                                "dropped_missing_zp": n_template[1]},
           "bright_planet_infov_dropped": n_planet_fov[0],
           "n_planet_masked_epochs": n_planet_masked,
           "construction_note": "v1.2 iteration: baseline-median centring, robust ring median, V<=8 star mask 3px, night medians",
           "fit_stats_astrom": {cam: {"n": v[0], "astrom_frac": v[1] / max(v[0], 1)}
                                for cam, v in fit_stats.items()},
           "units": results}
    name = "dev_search_v11.json" if role == "dev" else f"{role}_search_v1.json"
    (SURV / "results" / name).write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "units"}, indent=1))
    for u, r in out["units"].items():
        if r.get("status") != "searched":
            print(u, "->", r.get("status"))
            continue
        line = f'{u}: n_ev={r["per_dpa"]["dpa+0"]["n_events"]} zmed={r["per_dpa"]["dpa+0"]["z_median"]:+.2f} zmad={r["per_dpa"]["dpa+0"]["z_mad"]:.2f}'
        for s in ("S_stack", "S_event", "S_pulse"):
            if s in r:
                line += f' | {s} S={r[s]["S"]:.2f} T={r[s]["T"]:.2f} exc={r[s]["exceeds"]}'
        print(line)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "index"
    role = sys.argv[2] if len(sys.argv) > 2 else "dev"
    if cmd == "index":
        build_index(role)
    elif cmd == "measure":
        run_measure(role=role)
    elif cmd == "reduce":
        reduce_series()
    elif cmd == "reduce11":
        reduce_v11(role)
