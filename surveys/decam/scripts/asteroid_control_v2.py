"""DECam v2 positive control (hypotheses v1.0 §6/§8.6): recover
(60000) Miminko through the full v2 chain — adapter fetch, per-CCD
matched filter + canvas paste, in-frame NSC star calibration, and the
v2 ring-48 rule (sglsurvey.control.rescore).

For each exposure in configs/asteroid_control_v1.json: fetch the
dqmask + the CCD containing the Horizons position, build the flux map
centred there, star-calibrate it against a snapshotted NSC object cone
at the field, and re-score. Predicted magnitudes = Horizons V + solar
colours. Output: runs/decam/v2/control/60000{,_noclip}/summary.json.

Usage: uv run python surveys/decam/scripts/asteroid_control_v2.py
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from profile import PROFILE as P  # noqa: E402

from sglsurvey.adapters.base import CutoutSpec  # noqa: E402
from sglsurvey.adapters.noirlab_decam import (DQ_FATAL_DEFAULT,  # noqa: E402
                                              RETRIEVE_URL,
                                              DecamInstcalAdapter)
from sglsurvey.photometry import build_flux_map_decam  # noqa: E402
from sglsurvey.records import Observation  # noqa: E402
from sglsurvey.snapshots import SnapshotStore  # noqa: E402
from sglsurvey.control import rescore  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
CFG = REPO / "surveys" / "decam" / "configs" / "asteroid_control_v1.json"
WORK = REPO / "runs" / "decam" / "control_v1"
CANVAS_PIX = 2000   # +-4.4 arcmin: must contain the star-ZP cone
STAR_MAG = (16.0, 21.0)
MIN_STARS = 6
#: Solar colours, band = V + colour (v1 convention, SDSS-derived).
COLOUR = {"g": 0.25, "r": -0.19, "i": -0.33, "z": -0.37, "Y": -0.39}


def observation_from(e: dict) -> Observation:
    prods = {k: {"md5": m, "url": RETRIEVE_URL.format(md5=m)}
             for k, m in e["products"].items()}
    return Observation.build(
        archive_id="noirlab-decam", collection="decam/instcal",
        release="control", native_key={"expnum": e["expnum"],
                                       "image_md5": e["products"]["image"]},
        band=e["band"], wavelength_um=None,
        t_start_mjd_utc=e["mjd_mid"], t_mid_mjd_utc=e["mjd_mid"],
        t_stop_mjd_utc=e["mjd_mid"], exptime_s=e["exptime_s"],
        corners_icrs_deg=(), wcs={}, products=prods)


def star_zp(fm, ra, dec, band, store) -> tuple[float | None, int]:
    from dl import queryClient as qc
    col = {"g": "gmag", "r": "rmag", "i": "imag", "z": "zmag",
           "Y": "ymag"}.get(band)
    if col is None:
        return None, 0
    sql = (f"SELECT ra,dec,{col},class_star,ndet FROM nsc_dr2.object "
           f"WHERE q3c_radial_query(ra,dec,{ra:.6f},{dec:.6f},0.055) "
           f"AND {col} BETWEEN {STAR_MAG[0]} AND {STAR_MAG[1]} "
           f"AND class_star > 0.7 AND ndet >= 5")
    request_utc = datetime.now(timezone.utc).isoformat()
    text = qc.query(sql=sql, fmt="csv", timeout=120)
    store.store(service_url="datalab:queryClient/query", query=sql,
                request_utc=request_utc, response_bytes=text.encode(),
                row_count=text.count("\n") - 1, http_status=200)
    st = []
    for line in text.strip().splitlines()[1:]:
        p = line.split(",")
        try:
            st.append((float(p[0]), float(p[1]), float(p[2])))
        except ValueError:
            continue
    if not st:
        return None, 0
    st = np.array(st)
    f, v, g = fm.sample(st[:, 0], st[:, 1])
    ok = (np.isfinite(f) & (f > 0) & (g > 0.9) & np.isfinite(v)
          & (f / np.sqrt(np.maximum(v, 1e-30)) > 7))
    if ok.sum() < MIN_STARS:
        return None, int(ok.sum())
    zps = st[ok, 2] + 2.5 * np.log10(f[ok])
    med = float(np.median(zps))
    mad = 1.4826 * float(np.median(np.abs(zps - med)))
    keep = np.abs(zps - med) < 3 * max(mad, 0.02)
    return float(np.median(zps[keep])), int(keep.sum())


def main() -> None:
    cfg = json.loads(CFG.read_text())
    store = SnapshotStore(WORK)
    adapter = DecamInstcalAdapter()
    frames, maps, positions, vpred = [], [], [], {b: [] for b in P.bands}
    for e in cfg["exposures"]:
        hz = e.get("horizons") or {}
        if "vmag" not in hz:
            continue
        obs = observation_from(e)
        try:
            adapter.fetch(obs, ["dqmask"], WORK / "products")
            cut = CutoutSpec(ra_deg=hz["ra_deg"], dec_deg=hz["dec_deg"],
                             size_pix=CANVAS_PIX)
            ps = adapter.fetch(obs, ["image", "wtmap"], WORK / "products",
                               cutout=cut, dqmask_dir=WORK / "products")
        except Exception as exc:
            print(f"  EXPNUM {e['expnum']}: fetch failed "
                  f"({type(exc).__name__} {str(exc)[:80]})", flush=True)
            continue
        img = next(p for p in ps.products if p.kind == "image")
        wt = next((p for p in ps.products if p.kind == "wtmap"), None)
        extname = img.path.stem.rsplit("_", 1)[-1]
        fm = build_flux_map_decam(
            [(img.path, wt.path if wt else None, extname)],
            WORK / "products" / f"exp{e['expnum']}_dqmask.fits.fz",
            DQ_FATAL_DEFAULT, e["band"], e["mjd_mid"],
            (hz["ra_deg"], hz["dec_deg"]), CANVAS_PIX)
        if fm is None:
            print(f"  EXPNUM {e['expnum']}: no flux map", flush=True)
            continue
        zp, n_st = star_zp(fm, hz["ra_deg"], hz["dec_deg"], e["band"],
                           store)
        if zp is None:
            print(f"  EXPNUM {e['expnum']}: star ZP failed "
                  f"({n_st} stars)", flush=True)
            continue
        fm.magzp = zp
        frames.append({"band": e["band"], "mjd": e["mjd_mid"],
                       "v_pred": hz["vmag"]})
        maps.append(fm)
        positions.append((hz["ra_deg"], hz["dec_deg"]))
        vpred[e["band"]].append(hz["vmag"])
        print(f"  EXPNUM {e['expnum']} {e['band']} "
              f"{e['exptime_s']:.0f}s: zp* {zp:.3f} ({n_st} stars), "
              f"V_pred {hz['vmag']:.2f}", flush=True)
    pred = {b: float(-2.5 * np.log10(
        np.mean(10 ** (-0.4 * np.array(v)))) + COLOUR[b])
        for b, v in vpred.items() if v}
    rescore(P, frames, maps, positions, P.run_dir / "control" / "60000",
            "60000", predicted_mag=pred,
            extra={"source": "asteroid_control_v1.json exposures, "
                             "fetched via DecamInstcalAdapter",
                   "colours": COLOUR})
    rescore(P, frames, maps, positions,
            P.run_dir / "control" / "60000_noclip",
            "60000 (no single-epoch clip)", predicted_mag=pred,
            extra={"source": "asteroid_control_v1.json exposures",
                   "colours": COLOUR,
                   "note": "the frozen clip removes frames where a "
                           "bright mover is individually detected"},
            clip_sigma=None)


if __name__ == "__main__":
    main()
