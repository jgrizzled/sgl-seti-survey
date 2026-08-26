"""L2 colour-coefficient ensemble (amendment v1.1 input).

Refits a random sample of measured frames and records every matched
calibrator's (V, B-V, r_rsun, aperture flux, exptime, camera). Fits
the global Clear/Orange-vs-V colour coefficient c in
    zp_star = V + 2.5 log10(f) = ZP(r) + c (B-V - 0.65)
after removing the per-frame radial ZP, and reports the achievable
per-frame scatter after colour correction. Seeded (20260825).
"""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lasco_lib as L

REPO = Path(__file__).resolve().parents[3]
DEVDIR = REPO / "runs" / "heliospheric-crossings" / "dev"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "color_ensemble_v1.json"

N_PER_CAM = 400


def frame_calibrators(path: Path) -> list[dict]:
    fr = L.load_frame(path)
    L.fit_frame(fr)
    if fr.rot_deg is None:
        return []
    dyr = (fr.mjd_mid - 48348.5625) / 365.25
    ra = L._cat["ra"] + L._cat["pmra_masyr"] * dyr / 3.6e6 / np.cos(np.radians(L._cat["dec"]))
    dec = L._cat["dec"] + L._cat["pmde_masyr"] * dyr / 3.6e6
    sep = np.hypot((ra - fr.sun_ra) * np.cos(np.radians(dec)), dec - fr.sun_dec)
    m = sep < fr.raw.shape[0] / 2 * fr.scale_arcsec / 3600.0 * 1.5
    ra, dec, v, bv = ra[m], dec[m], L._cat["vmag"][m], L._cat["bv"][m]
    px, py = L._sky_to_pix(fr, ra, dec, fr.rot_deg)
    px, py = px + fr.translation[0], py + fr.translation[1]
    rsun_px = fr.rsun_arcsec / fr.scale_arcsec
    r_lo, r_hi = (2.2, 6.0) if fr.camera == "c2" else (4.4, 29.0)
    out = []
    for i in range(len(ra)):
        if v[i] < 4.5 or not (-0.5 <= bv[i] <= 2.0):
            continue
        if not (10 < px[i] < fr.raw.shape[1] - 10 and 10 < py[i] < fr.raw.shape[0] - 10):
            continue
        r_star = float(np.hypot(px[i] - fr.crpix[0] - fr.translation[0],
                                py[i] - fr.crpix[1] - fr.translation[1]) / rsun_px)
        if not (r_lo <= r_star <= r_hi):
            continue
        fx, er = L._aper_flux(fr.img, px[i], py[i])
        if np.isfinite(fx) and er > 0 and fx / er >= 10.0:
            out.append({"V": float(v[i]), "bv": float(bv[i]), "r": r_star,
                        "zp": float(v[i] + 2.5 * np.log10(fx)),
                        "cam": fr.camera, "mjd": fr.mjd_mid})
    return out


def main() -> None:
    frames = {"c2": [], "c3": []}
    with open(DEVDIR / "measurements_v1.jsonl") as fh:
        for line in fh:
            rec = json.loads(line)
            if "fit" in rec and rec["fit"]["rot"] is not None:
                cam = "c2" if "/c2/" in rec["frame"] else "c3"
                frames[cam].append(rec["frame"])
    rng = random.Random(20260825)
    cal = []
    for cam in ("c2", "c3"):
        pick = rng.sample(frames[cam], min(N_PER_CAM, len(frames[cam])))
        for i, rp in enumerate(pick):
            try:
                cal.extend(frame_calibrators(REPO / rp))
            except Exception:
                pass
            if (i + 1) % 100 == 0:
                print(f"{cam} {i+1}", flush=True)

    result = {"n_records": len(cal)}
    for cam in ("c2", "c3"):
        rows = [c for c in cal if c["cam"] == cam]
        if len(rows) < 50:
            result[cam] = {"n": len(rows), "note": "too few"}
            continue
        bv = np.array([c["bv"] for c in rows])
        zp = np.array([c["zp"] for c in rows])
        r = np.array([c["r"] for c in rows])
        # remove the radial profile (median per 3-Rsun bin over ensemble)
        binw = 1.0 if cam == "c2" else 3.0
        prof = {}
        for lo in np.arange(0, r.max() + binw, binw):
            mm = (r >= lo) & (r < lo + binw)
            if mm.sum() >= 10:
                prof[lo] = np.median(zp[mm])
        keys = sorted(prof)
        resid = zp - np.interp(r, keys, [prof[k] for k in keys])
        # robust linear fit resid vs (bv - 0.65): iterate one clip
        use = np.abs(resid - np.median(resid)) < 3 * np.median(np.abs(resid - np.median(resid))) * 1.4826 + 1e-9
        for _ in range(2):
            A = np.stack([np.ones(use.sum()), bv[use] - 0.65], 1)
            coef, *_ = np.linalg.lstsq(A, resid[use], rcond=None)
            fitted = coef[0] + coef[1] * (bv - 0.65)
            sc = np.abs(resid - fitted)
            use = sc < 3 * np.median(sc[use]) * 1.4826 + 1e-9
        after = resid - fitted
        result[cam] = {
            "n": len(rows), "n_used": int(use.sum()),
            "color_coeff_mag_per_bv": round(float(coef[1]), 4),
            "scatter_before_mad": round(float(np.median(np.abs(resid - np.median(resid))) * 1.4826), 3),
            "scatter_after_mad": round(float(np.median(np.abs(after[use] - np.median(after[use]))) * 1.4826), 3),
        }
    np.savez_compressed(DEVDIR / "color_ensemble_records.npz",
                        **{k: np.array([c[k] for c in cal]) for k in ("V", "bv", "r", "zp", "mjd")},
                        cam=np.array([c["cam"] for c in cal]))
    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
