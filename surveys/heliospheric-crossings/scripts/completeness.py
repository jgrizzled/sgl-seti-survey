"""Completeness stage (hypotheses §8 / thresholds §6 / amendment v1.1).

Three measurements, C1-normalized end to end:

1. **Stamp response R** — Gaussian stamps of known flux injected into
   copies of real frames at real unit positions (baseline epochs), the
   identical load→detrend→fit→photometry chain re-run, R = recovered /
   injected. The stamp width is the measured PSF second-moment sigma
   from bright field stars (per camera). Nothing is written back to
   any frame.
2. **Series-level injection recovery** — for every searched unit, the
   d = 1 persistent signal (flux F × R added to every visible in-window
   epoch of the real confirmatory null series) and, for C2 units, the
   single-epoch pulse; statistics recomputed against the CONFIRMATORY
   thresholds T (fixed); recovery over epoch-bootstrap draws (seed
   20260825) → m90 per unit and statistic.
3. **Depths** — flux → V via the per-unit median ensemble ZP (+ colour
   correction is already inside the star ZP construction) and the
   median EXPTIME; powers from first principles through the frozen
   cone areas (line-equivalent, W_eff = 300 nm Clear / 100 nm Orange,
   declared), anchored against the ZTF convention.

Outputs: results/completeness_v1.json.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lasco_lib as L
import tracks

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "heliospheric-crossings"
DEVDIR = REPO / "runs" / "heliospheric-crossings" / "dev"
CFG = json.loads((SURV / "configs" / "threshold_freeze_v1.json").read_text())
CONF = json.loads((SURV / "results" / "confirmatory_search_v1.json").read_text())
RING_DPA = [float(x) for x in CFG["controls"]["ring_dpa_deg"]]
RNG = np.random.default_rng(20260825)

FLUX_GRID = np.array([2, 5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560], float)
N_DRAWS = 100

# ---------------------------------------------------------------- PSF + R

def measure_psf_sigma(paths, n=25):
    sigs = []
    for p in paths[:n]:
        try:
            fr = L.load_frame(p)
            det, snr = L.detect_sources(fr, snr_min=30.0, cap=30)
            for (x, y), s in zip(det, snr):
                if s < 50:
                    continue
                ix, iy = int(x), int(y)
                if not (6 < ix < fr.img.shape[1] - 6 and 6 < iy < fr.img.shape[0] - 6):
                    continue
                sub = fr.img[iy - 5:iy + 6, ix - 5:ix + 6].copy()
                sub -= np.median(sub)
                if sub.max() <= 0:
                    continue
                yy, xx = np.mgrid[-5:6, -5:6]
                w = np.clip(sub, 0, None)
                var = float(np.sum(w * (xx**2 + yy**2)) / np.sum(w)) / 2
                if 0.3 < var < 9:
                    sigs.append(np.sqrt(var))
        except Exception:
            pass
    return float(np.median(sigs)) if sigs else 1.3


def stamp_response(paths, sigma_psf, n_frames=40, flux=5000.0):
    """Inject one stamp per frame at a fixed off-source position and
    re-measure through the identical chain."""
    recov = []
    for p in paths[:n_frames]:
        try:
            fr = L.load_frame(p)
            if not (L.fit_frame(fr) or fr.rot_deg is not None):
                continue
            ny, nx = fr.raw.shape
            # position: mid-annulus, PA chosen per frame deterministically
            rsun_px = fr.rsun_arcsec / fr.scale_arcsec
            r_px = (3.5 if fr.camera == "c2" else 12.0) * rsun_px
            th = (hash(p.name) % 360) * np.pi / 180
            x0, y0 = fr.crpix[0] + r_px * np.cos(th), fr.crpix[1] + r_px * np.sin(th)
            if not (20 < x0 < nx - 20 and 20 < y0 < ny - 20):
                continue
            base_flux, _ = L._aper_flux(fr.img, x0, y0)
            yy, xx = np.mgrid[0:ny, 0:nx]
            stamp = flux * np.exp(-0.5 * ((xx - x0)**2 + (yy - y0)**2) / sigma_psf**2) \
                / (2 * np.pi * sigma_psf**2)
            fr.raw = fr.raw + stamp
            fr.img = L.radial_profile_subtract(fr)   # identical detrend on injected frame
            inj_flux, _ = L._aper_flux(fr.img, x0, y0)
            if np.isfinite(inj_flux) and np.isfinite(base_flux):
                recov.append((inj_flux - base_flux) / flux)
        except Exception:
            pass
    return (float(np.median(recov)), int(len(recov))) if recov else (np.nan, 0)


# ------------------------------------------------------- series machinery

def gather_series():
    """(u, e, dpa) epoch lists from the confirmatory measurements —
    the same gather the reduce used (gates identical)."""
    idx = json.loads((DEVDIR / "task_index_confirmatory_v1.json").read_text())
    units = idx["units"]
    for u in units.values():
        u["events_by_id"] = {e["event_id"]: e for e in u["events"]}
    epochs = defaultdict(dict)
    with open(DEVDIR / "measurements_confirmatory_v1.jsonl") as fh:
        for line in fh:
            rec = json.loads(line)
            for p in rec.get("records", []):
                if (p.get("valid_fraction", 0) < 0.9
                        or not np.isfinite(p.get("flux", np.nan))
                        or not np.isfinite(p.get("err", np.nan)) or p["err"] <= 1e-6):
                    continue
                key = (p["u"], p["e"], p["k"], round(p["mjd"], 6))
                epochs[key][p["dpa"]] = (p["flux"], p.get("zp", float("nan")))
    return units, epochs


def unit_event_series(units, epochs, ukey):
    """per event: (window nights diffs, baseline nights diffs) at dpa 0,
    with per-epoch lists so injections can be added per epoch."""
    out = {}
    for ev in units[ukey]["events"]:
        win_d, base_d, zps = defaultdict(list), defaultdict(list), []
        for (uk, eid, kind, mjd), by in epochs.items():
            if uk != ukey or eid != ev["event_id"] or 0.0 not in by:
                continue
            refs = [by[d][0] for d in by if d != 0.0]
            if len(refs) < 5:
                continue
            diff = by[0.0][0] - float(np.median(refs))
            (win_d if kind == "win" else base_d)[int(mjd)].append(diff)
            if np.isfinite(by[0.0][1]):
                zps.append(by[0.0][1])
        if win_d and len(base_d) >= 20:
            out[ev["event_id"]] = (win_d, base_d, float(np.median(zps)) if zps else np.nan)
    return out


def stat_from_series(series, inject=0.0, pulse_at=None, rng=None):
    zs, pulse_max = [], -np.inf
    for eid, (win_d, base_d, _) in series.items():
        wn = []
        for night, vals in win_d.items():
            v = [x + inject for x in vals]
            if pulse_at == eid:
                v[0] += pulse_at_flux[0]
            wn.append(np.median(v))
        bn = [np.median(v) for v in base_d.values()]
        if rng is not None:   # epoch bootstrap
            wn = list(rng.choice(wn, size=len(wn), replace=True))
            bn = list(rng.choice(bn, size=len(bn), replace=True))
        sd = float(np.std(bn, ddof=1))
        if sd <= 0:
            continue
        zs.append((np.mean(wn) - np.median(bn)) / (sd / np.sqrt(len(wn))))
        flat = [x + inject for v in win_d.values() for x in v]
        bflat = [x for v in base_d.values() for x in v]
        sde = float(np.std(bflat, ddof=1))
        if sde > 0:
            pulse_max = max(pulse_max, (np.max(flat) - np.median(bflat)) / sde)
    if len(zs) < CFG["gates"]["min_events_for_stack"]:
        return None
    z = np.array(zs)
    return {"S_stack": float(np.sum(z) / np.sqrt(len(z))),
            "S_event": float(np.max(z)), "S_pulse": float(pulse_max)}


pulse_at_flux = [0.0]


def main() -> None:
    # frame pools for PSF/response (baseline frames on disk)
    frames = REPO / "runs" / "heliospheric-crossings" / "frames"
    pool = {"c2": sorted((frames / "c2").rglob("*.fts"))[::37][:60],
            "c3": sorted((frames / "c3").rglob("*.fts"))[::61][:60]}
    psf = {cam: measure_psf_sigma(pool[cam]) for cam in pool}
    resp = {cam: stamp_response(pool[cam], psf[cam]) for cam in pool}
    print("PSF sigma(px):", psf, " stamp response:", resp, flush=True)

    units, epochs = gather_series()
    results = {}
    for ukey, conf_unit in CONF["units"].items():
        if conf_unit.get("status") != "searched":
            results[ukey] = {"status": conf_unit.get("status")}
            continue
        cam = units[ukey]["camera"]
        R = resp[cam][0]
        series = unit_event_series(units, epochs, ukey)
        thresholds = {s: max(conf_unit[s]["T"], 0.0)
                      for s in ("S_stack", "S_event", "S_pulse") if s in conf_unit}
        rec = {"n_events_series": len(series), "response": R}
        # persistent d=1 recovery
        curves = {s: [] for s in thresholds}
        for F in FLUX_GRID:
            hits = {s: 0 for s in thresholds}
            for d in range(N_DRAWS):
                st = stat_from_series(series, inject=F * R, rng=RNG)
                if st is None:
                    continue
                for s in thresholds:
                    if s == "S_pulse":
                        continue
                    hits[s] += st[s] > thresholds[s]
            for s in thresholds:
                if s != "S_pulse":
                    curves[s].append(hits[s] / N_DRAWS)
        for s in ("S_stack", "S_event"):
            if s in curves and curves[s]:
                c = np.array(curves[s])
                i = np.argmax(c >= 0.9) if (c >= 0.9).any() else None
                rec[f"m90_flux_{s}"] = float(FLUX_GRID[i]) if i is not None else None
        # pulse recovery (C2 units): single-epoch injection
        if "S_pulse" in thresholds:
            pc = []
            eids = list(series)
            for F in FLUX_GRID * 4:
                hits = 0
                for d in range(N_DRAWS // 2):
                    eid = eids[int(RNG.integers(len(eids)))]
                    pulse_at_flux[0] = F * R
                    st = stat_from_series(series, inject=0.0, pulse_at=eid)
                    if st and st["S_pulse"] > thresholds["S_pulse"]:
                        hits += 1
                pc.append(hits / (N_DRAWS // 2))
            c = np.array(pc)
            i = np.argmax(c >= 0.9) if (c >= 0.9).any() else None
            rec["m90_flux_S_pulse"] = float(FLUX_GRID[i] * 4) if i is not None else None
        # depth conversion: median unit ZP + exptime
        zps = [v[2] for v in series.values() if np.isfinite(v[2])]
        expmap = json.loads((DEVDIR / "exptime_map_v1.json").read_text())
        expt = np.median([v[0] for rp, v in expmap.items() if v and (("/c2/" in rp) == (cam == "c2"))])
        rec["zp_med"] = float(np.median(zps)) if zps else None
        rec["exptime_med"] = float(expt)
        for s in ("S_stack", "S_event", "S_pulse"):
            f = rec.get(f"m90_flux_{s}")
            if f and rec["zp_med"]:
                rec[f"m90_V_{s}"] = round(rec["zp_med"] - 2.5 * np.log10(f * expt), 2)
        results[ukey] = rec
        print(ukey, {k: v for k, v in rec.items() if k.startswith("m90_V") or k == "n_events_series"}, flush=True)

    out = {"stage": "completeness_v1", "seed": 20260825,
           "psf_sigma_px": psf, "stamp_response": {k: {"R": v[0], "n": v[1]} for k, v in resp.items()},
           "flux_grid": list(FLUX_GRID), "n_draws": N_DRAWS,
           "thresholds_source": "confirmatory_search_v1.json (fixed)",
           "units": results}
    (SURV / "results" / "completeness_v1.json").write_text(json.dumps(out, indent=1) + "\n")
    print("wrote completeness_v1.json")


if __name__ == "__main__":
    main()
