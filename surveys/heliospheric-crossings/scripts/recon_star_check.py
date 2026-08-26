"""LASCO recon: end-to-end astrometry + star-visibility check on one
level-1 C3 frame (2010-06-15 00:18 UT, 35227512.fts.gz, fetched to
runs/heliospheric-crossings/recon/).

Level-1 headers are helioprojective (Sun-center CRPIX, solar-north-up
after derolling, CROTA residual roll) with no celestial WCS. The
adapter plan is: build a synthetic celestial WCS from the Sun's
apparent position + solar P-angle + CROTA, then refine per-frame on
catalog stars. This script tests that plan: predict pixels for three
bright Taurus stars in the field, search a window, report hits, the
common offset (SOHO-vs-geocenter Sun parallax + pointing error), and
detection S/N. The orientation convention (rotation sign, X flip) is
determined empirically by which combo lands all stars at one common
offset — the result documents the convention for the adapter.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import get_sun
from astropy.io import fits
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
FRAME = REPO / "runs" / "heliospheric-crossings" / "recon" / "c3_L1_35227512.fts.gz"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "recon_star_check_v1.json"

# (name, ICRS ra deg, dec deg, V)
STARS = [
    ("bet-tau", 81.5730, 28.6074, 1.65),
    ("zet-tau", 84.4112, 21.1425, 3.01),
    ("119-tau", 83.0529, 18.5942, 4.32),
]

SEARCH_PX = 30  # window half-size: absorbs SOHO parallax + pointing error


def solar_p_angle_deg(t: Time) -> float:
    """Position angle of solar north vs celestial north (Meeus ch. 29)."""
    jd = t.tt.jd
    theta = np.nan  # unused; kept minimal
    d = jd - 2451545.0
    eps = np.radians(23.4393 - 3.563e-7 * d)
    incl = np.radians(7.25)
    omega = np.radians(73.6667 + 1.3958333 * (jd - 2396758.0) / 36525.0)
    # apparent solar longitude (low precision is fine at 0.01 deg here)
    g = np.radians((357.529 + 0.98560028 * d) % 360.0)
    q = 280.459 + 0.98564736 * d
    lam = np.radians((q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g)) % 360.0)
    x = np.arctan(-np.cos(lam) * np.tan(eps))
    y = np.arctan(-np.cos(lam - omega) * np.tan(incl))
    return float(np.degrees(x + y))


def main() -> None:
    with gzip.open(FRAME) as g:
        hdu = fits.open(g)[0]
        hdr, img = hdu.header, np.asarray(hdu.data, dtype=float)

    t = Time(hdr["DATE-OBS"], scale="utc")
    sun = get_sun(t)  # geocentric apparent; SOHO offset absorbed in search window
    p_ang = solar_p_angle_deg(t)
    scale = hdr["CDELT1"]  # arcsec/px
    crpix = np.array([hdr["CRPIX1"], hdr["CRPIX2"]])  # 1-based, Sun center
    crota = float(hdr["CROTA"])

    sun_ra, sun_dec = sun.ra.deg, sun.dec.deg

    # point-source extraction: high-pass against a running median, then
    # the brightest local maxima (the F-corona gradient is the noise here)
    from scipy.ndimage import maximum_filter, median_filter

    hp = img - median_filter(img, size=11)
    sigma = np.median(np.abs(hp)) * 1.4826 + 1e-30
    snr = hp / sigma
    is_peak = (maximum_filter(snr, size=7) == snr) & (snr > 8)
    is_peak[:20, :] = is_peak[-20:, :] = is_peak[:, :20] = is_peak[:, -20:] = False
    pys, pxs = np.nonzero(is_peak)
    order = np.argsort(snr[pys, pxs])[::-1][:300]
    det = np.stack([pxs[order], pys[order]], axis=1).astype(float)
    det_snr = snr[pys[order], pxs[order]]

    # pattern match per orientation: one common translation must land
    # every star on a detection (translation absorbs SOHO parallax +
    # pointing error; the inter-star vectors pin rotation/flip)
    results = {}
    for rot_sign in (+1, -1):
        for xflip in (+1, -1):
            key = f"rot{rot_sign:+d}_x{xflip:+d}"
            rot = np.radians(rot_sign * p_ang + crota)
            pred = []
            for name, ra, dec, vmag in STARS:
                d_e = (ra - sun_ra) * np.cos(np.radians(dec)) * 3600.0  # East +
                d_n = (dec - sun_dec) * 3600.0  # North +
                e = d_e * np.cos(rot) - d_n * np.sin(rot)
                n = d_e * np.sin(rot) + d_n * np.cos(rot)
                # celestial East appears LEFT when North is up -> -x; xflip tests it
                px = crpix[0] - 1 + xflip * (-e / scale)
                py = crpix[1] - 1 + n / scale
                pred.append((name, vmag, px, py))
            # candidate translations: anchor = brightest star to any detection
            best = None
            _, _, ax, ay = pred[0]
            for dx, dy in det - np.array([ax, ay]):
                if abs(dx) > SEARCH_PX or abs(dy) > SEARCH_PX:
                    continue
                matches, resid = [], 0.0
                for name, vmag, px, py in pred:
                    d = np.hypot(det[:, 0] - (px + dx), det[:, 1] - (py + dy))
                    j = int(np.argmin(d))
                    matches.append((name, vmag, float(d[j]), float(det_snr[j]), det[j].tolist()))
                    resid += min(float(d[j]), 10.0)
                n_hit = sum(1 for m in matches if m[2] < 4.0)
                score = (n_hit, -resid)
                if best is None or score > best[0]:
                    best = (score, [dx, dy], matches)
            if best is None:
                results[key] = {"n_matched": 0}
                continue
            (n_hit, _), t, matches = best
            results[key] = {
                "n_matched_lt_4px": n_hit,
                "translation_px": [round(t[0], 1), round(t[1], 1)],
                "stars": [
                    {"star": m[0], "V": m[1], "resid_px": round(m[2], 1), "det_snr": round(m[3], 1)}
                    for m in matches
                ],
            }

    out = {
        "frame": FRAME.name,
        "date_obs": hdr["DATE-OBS"],
        "sun_radec_geocentric": [round(sun_ra, 4), round(sun_dec, 4)],
        "p_angle_deg": round(p_ang, 3),
        "crota_deg": crota,
        "platescale_arcsec": scale,
        "search_halfwidth_px": SEARCH_PX,
        "orientations": results,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
