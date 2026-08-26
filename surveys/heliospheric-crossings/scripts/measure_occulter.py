"""Measure the usable occulter inner radius (and outer edge) of LASCO
C2/C3 from level-1 frames (freeze D7 gate input).

Fetches a small cross-era sample of NRL level-1 frames, computes the
radial valid-pixel fraction about the Sun-center CRPIX, and reports
the inner radius where the valid fraction first exceeds 0.5 and the
outer radius where it last exceeds 0.5, in R_sun (header RSUN).
Writes results/occulter_radii_v1.json; frames cached under
runs/heliospheric-crossings/recon/.
"""

from __future__ import annotations

import gzip
import io
import json
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO = Path(__file__).resolve().parents[3]
CACHE = REPO / "runs" / "heliospheric-crossings" / "recon"
OUT = REPO / "surveys" / "heliospheric-crossings" / "results" / "occulter_radii_v1.json"

SAMPLES = [  # (camera, yymmdd) — spread across eras with level-1 coverage
    ("c2", "970601"), ("c2", "040601"), ("c2", "100615"), ("c2", "150101"),
    ("c3", "970601"), ("c3", "100615"), ("c3", "150101"),
]


def first_l1_frame(cam: str, day: str) -> tuple[str, bytes]:
    base = f"https://lasco-www.nrl.navy.mil/lz/level_1/{day}/{cam}/"
    html = urllib.request.urlopen(base, timeout=30).read().decode()
    import re

    names = re.findall(r'href="(\d+\.fts(?:\.gz)?)"', html)
    if not names:
        raise RuntimeError(f"no L1 frames {day}/{cam}")
    name = names[len(names) // 2]  # midday-ish frame
    dest = CACHE / f"{cam}_L1_{day}_{name}"
    if not dest.exists():
        data = urllib.request.urlopen(base + name, timeout=120).read()
        if len(data) < 100_000:
            raise RuntimeError(f"truncated fetch {name}: {len(data)} B")
        dest.write_bytes(data)
    return name, dest.read_bytes()


def radial_valid_profile(hdr, img) -> dict:
    """Occulter geometry from a level-1 frame.

    C3's occulted core is zero-filled, so a valid-fraction threshold
    finds the inner edge. C2's is not — but its radial median profile
    peaks sharply at the occulter-edge diffraction/straylight ring and
    drops inward, so the ring radius marks the edge and usable corona
    starts just outside it. Outer edge capped at the inscribed circle.
    """
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r_px = np.hypot(xx - (hdr["CRPIX1"] - 1), yy - (hdr["CRPIX2"] - 1))
    rsun_px = hdr["RSUN"] / hdr["CDELT1"]  # solar radius in pixels
    r = r_px / rsun_px
    valid = np.isfinite(img) & (img != 0)
    bins = np.arange(0.0, 34.0, 0.1)
    frac, med = [], []
    for lo in bins:
        m = (r >= lo) & (r < lo + 0.1)
        frac.append(valid[m].mean() if m.any() else 0.0)
        med.append(float(np.median(img[m])) if m.any() else 0.0)
    frac, med = np.array(frac), np.array(med)
    ok = frac > 0.5
    inner_valid = float(bins[np.argmax(ok)]) if ok.any() else None
    inscribed = (min(nx, ny) / 2 - 4) / rsun_px
    # diffraction-ring radius: brightest annulus within r < 6 Rsun
    core = bins < 6.0
    ring = float(bins[core][np.argmax(med[core])])
    return {"inner_rsun_valid_gt50": inner_valid,
            "ring_peak_rsun": ring,
            "inner_rsun_usable": max(inner_valid or 0.0, ring + 0.2),
            "outer_rsun_usable": float(min(
                bins[len(ok) - 1 - np.argmax(ok[::-1])] if ok.any() else 0.0,
                inscribed)),
            "rsun_px": float(rsun_px)}


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    results = {}
    for cam, day in SAMPLES:
        try:
            name, blob = first_l1_frame(cam, day)
            if name.endswith(".gz"):
                blob = gzip.decompress(blob)
            hdu = fits.open(io.BytesIO(blob))[0]
            prof = radial_valid_profile(hdu.header, np.asarray(hdu.data, float))
            prof["frame"] = name
            results[f"{cam}_{day}"] = prof
            print(cam, day, name, prof["inner_rsun_usable"], "->", prof["outer_rsun_usable"],
                  f"(ring {prof['ring_peak_rsun']})")
        except Exception as exc:  # keep going; report per-sample failures
            results[f"{cam}_{day}"] = {"error": str(exc)}
            print(cam, day, "ERROR", exc)

    def agg(cam):
        vals = [v for k, v in results.items() if k.startswith(cam) and "inner_rsun_usable" in v]
        return {
            "inner_rsun_adopted_max": max(v["inner_rsun_usable"] for v in vals) if vals else None,
            "outer_rsun_adopted_min": min(v["outer_rsun_usable"] for v in vals) if vals else None,
            "n_samples": len(vals),
        }

    # Adopted gate values (freeze D7), by inspection of the per-era
    # profiles: C2 inner = 2.2 (diffraction ring at 1.9-2.0, usable
    # past it; consistent 1997-2015), outer = 6.3 (inscribed edge;
    # the 1997 C2 sample is a masked/binned subframe — excluded).
    # C3: occulted core is a constant fill plateau (MAD 0) to r~4.2
    # in 2010/2015 and zero-masked to 4.4 in 1997 → inner = 4.4
    # conservative; outer = 29 (inscribed).
    adopted = {"c2_rsun": [2.2, 6.3], "c3_rsun": [4.4, 29.0]}
    out = {"samples": results, "c2": agg("c2"), "c3": agg("c3"),
           "adopted": adopted,
           "note": "valid = finite & nonzero on L1 frames; automated "
                   "metrics are era-heterogeneous (fill conventions "
                   "differ) — the 'adopted' block is the freeze gate "
                   "input, rationale in comments here"}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"c2": out["c2"], "c3": out["c3"]}, indent=1))


if __name__ == "__main__":
    main()
