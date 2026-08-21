"""Delete regenerable image products (locus cutouts, calibration cutouts,
full masks) for corridors whose sample tensors exist, keeping manifests,
cutout index, checksums and tensors (plan section 9 data retention:
content hashes + durable archive identifiers are kept; large images are
re-downloadable with fetch_cutouts.py / precise_pass.py).

Usage: uv run python surveys/panstarrs/scripts/purge_products.py <corridor ...>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ps1_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
P = REPO / "runs" / "panstarrs"


def main(corridors: list[str]) -> None:
    corridors = set(corridors)
    tens = {p.stem for p in (P / "calib_v1" / "tensors").glob("*.npz")}
    ready = set()
    for c in corridors:
        eps = [e for e, cc in CORRIDOR_OF.items() if cc == c]
        if eps and all(f"{e}__{r}" in tens for e in eps for r in ("rx", "tx")):
            ready.add(c)
    if corridors - ready:
        print(f"skipping (tensors missing): {sorted(corridors - ready)}")
    if not ready:
        return
    # observations of the ready corridors (via precise evaluations)
    obs_ids = set()
    for r in read_records(P / "precise_v1" / "records"
                          / "intersection_evaluation.jsonl"):
        if CORRIDOR_OF.get(r["endpoint_id"]) in ready:
            obs_ids.add(r["observation_id"])
    # but keep observations still needed by a not-ready corridor
    keep = set()
    for r in read_records(P / "precise_v1" / "records"
                          / "intersection_evaluation.jsonl"):
        c = CORRIDOR_OF.get(r["endpoint_id"])
        if c is not None and c not in ready and r["observation_id"] in obs_ids:
            keep.add(r["observation_id"])
    obs_ids -= keep
    n = freed = 0
    with (P / "products" / "cut" / "manifest.jsonl").open() as fh:
        for line in fh:
            m = json.loads(line)
            if m["observation_id"] not in obs_ids:
                continue
            paths = [P / "products" / "cut" / f for f in m["files"].values()]
            if m.get("msk"):
                paths.append(P / "products" / "msk" / m["msk"])
                stem = m["msk"].replace(".mask.fits", ".fits")
                paths += list((P / "products" / "calcut").glob(f"cut*-{stem}"))
            for q in paths:
                if q.exists():
                    freed += q.stat().st_size
                    q.unlink()
                    n += 1
    print(f"purged {n} files ({freed / 1e9:.1f} GB) for {sorted(ready)}")


if __name__ == "__main__":
    main(sys.argv[1:])
