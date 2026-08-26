"""One-time fetch of the bright-star catalog for per-frame WCS/ZP fits.

Hipparcos (VizieR I/239/hip_main) V <= 9.5 with proper motions, via the
VizieR TAP sync endpoint. Raw response snapshotted; npz written to
runs/heliospheric-crossings/starcat/hip_v95.npz (ra/dec deg ICRS epoch
1991.25, pm mas/yr, Vmag).
"""

from __future__ import annotations

import hashlib
import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
DEST = REPO / "runs" / "heliospheric-crossings" / "starcat"

URL = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/239/hip_main"
       "&-out=RAICRS,DEICRS,pmRA,pmDE,Vmag,B-V&Vmag=%3C=9.5&-out.max=200000")


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(URL, headers={"User-Agent": "sgl-seti-survey/1.0"})
    raw = urllib.request.urlopen(req, timeout=300).read()
    (DEST / "hip_v95_response.tsv").write_bytes(raw)
    rows = []
    for line in raw.decode("utf-8", "replace").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 6:
            continue
        try:
            ra, de, v = float(parts[0]), float(parts[1]), float(parts[4])
        except ValueError:
            continue  # column-name / unit / separator rows
        def f(x):
            try:
                return float(x)
            except ValueError:
                return 0.0
        rows.append((ra, de, f(parts[2]), f(parts[3]), v, f(parts[5]) if parts[5].strip() else 99.0))
    arr = np.array(rows)
    np.savez_compressed(DEST / "hip_v95.npz", ra=arr[:, 0], dec=arr[:, 1],
                        pmra_masyr=arr[:, 2], pmde_masyr=arr[:, 3], vmag=arr[:, 4], bv=arr[:, 5],
                        epoch_jyear=np.array(1991.25))
    meta = {"url": URL, "n_stars": int(len(arr)),
            "raw_sha256": hashlib.sha256(raw).hexdigest(),
            "npz_sha256": hashlib.sha256((DEST / "hip_v95.npz").read_bytes()).hexdigest()}
    (DEST / "hip_v95_meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
