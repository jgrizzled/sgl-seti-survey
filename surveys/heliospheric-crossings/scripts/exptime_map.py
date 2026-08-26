"""Header-only EXPTIME (+ MID-time) map over the fetched frames —
joins the L3 stellar-flux template (flux = 10^((zp - V*)/2.5)/exptime)
to the exptime-normalized series without re-measuring."""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from astropy.io import fits

REPO = Path(__file__).resolve().parents[3]
FRAMES = REPO / "runs" / "heliospheric-crossings" / "frames"
OUT = REPO / "runs" / "heliospheric-crossings" / "dev" / "exptime_map_v1.json"


def one(p: Path):
    try:
        h = fits.getheader(p)
        mjd = float(h.get("MID_DATE", 0)) + float(h.get("MID_TIME", 0)) / 86400.0
        return str(p.relative_to(REPO)), [float(h.get("EXPTIME", 0)), round(mjd, 6)]
    except Exception:
        return str(p.relative_to(REPO)), None


def main() -> None:
    paths = sorted(FRAMES.rglob("*.fts"))
    with ThreadPoolExecutor(8) as ex:
        out = dict(ex.map(one, paths))
    OUT.write_text(json.dumps(out) + "\n")
    ok = sum(1 for v in out.values() if v)
    print(f"{ok}/{len(out)} headers mapped")


if __name__ == "__main__":
    main()
