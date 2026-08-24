"""Step A: hypothesis v2.0 freeze with the stratified corridor hold-out
(hypotheses v2.0 §5, §7). Writes configs/v2_0_freeze.json.

The split unit is the corridor (both roles and all components of a
system share frames), stratified by the v1 overlay confusion class,
drawn with the recorded seed. Re-running with the same inputs
reproduces the same file byte for byte; any change of inputs changes
the freeze hash, which every v2 AnalysisRun pins.

Usage: uv run python surveys/wise/scripts/freeze_v2.py
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402
from sglseti import canonical_json, load_target_registry  # noqa: E402
from sglsurvey.corridors import CORRIDOR_OF, MEMBERS  # noqa: E402

OVERLAY = C.REPO / "surveys" / "wise" / "targets" / "overlay_v2.json"
# Overlay system name -> corridor where the normalised names differ.
MANUAL = {"lalande": "Lalande 21185", "proxima": "Proxima Cen",
          "alphacen": "Proxima Cen", "gj54": "YZ Cet",
          "wise0855": "WISE 0855-0714", "gj667": "GJ 667 ABC",
          "ltt1445": "LTT 1445 ABC"}


def _norm(s: str) -> str:
    s = re.sub(r"\(.*?\)", "", s.lower())
    s = s.replace("'s star", "").replace(" star", "").replace(" ab", "")
    return re.sub(r"[^a-z0-9]", "", s)


def confusion_classes() -> dict[str, str]:
    systems = json.loads(OVERLAY.read_text())["systems"]
    by_norm = {_norm(s["system"]): s["confusion"] for s in systems}
    by_name = {s["system"]: s["confusion"] for s in systems}
    alias = {"cyg61": "61cyg", "grb34": "groombridge34"}
    out = {}
    for c in MEMBERS:
        if c in MANUAL:
            out[c] = by_name[MANUAL[c]]
        else:
            key = alias.get(c, _norm(c))
            if key not in by_norm:
                raise KeyError(f"no overlay confusion class for corridor {c}")
            out[c] = by_norm[key]
    return out


def main() -> None:
    registry = load_target_registry(C.REGISTRY_PATH)
    classes = confusion_classes()
    strata = defaultdict(list)
    for c in sorted(MEMBERS):
        strata[classes[c]].append(c)
    rng = np.random.default_rng(C.SPLIT_SEED)
    dev, conf = [], []
    per_stratum = {}
    for cls in sorted(strata):
        cs = sorted(strata[cls])
        n_dev = max(1, int(round(C.DEV_FRACTION * len(cs))))
        pick = sorted(rng.choice(cs, size=n_dev, replace=False).tolist())
        dev.extend(pick)
        conf.extend(c for c in cs if c not in pick)
        per_stratum[cls] = {"n_corridors": len(cs), "n_dev": n_dev, "dev": pick}
    dev, conf = sorted(dev), sorted(conf)
    dev_endpoints = sorted(e for c in dev for e in MEMBERS[c])
    conf_endpoints = sorted(e for c in conf for e in MEMBERS[c])
    assert set(dev_endpoints) | set(conf_endpoints) == set(CORRIDOR_OF)
    assert not set(dev) & set(conf)

    freeze = {
        "hypothesis_version": C.HYPOTHESIS_VERSION,
        "hypotheses_hash": C.sha256_file(C.HYPOTHESES_PATH),
        "v1_hypotheses_hash": C.sha256_file(C.REPO / "surveys" / "wise" / "hypotheses_v1.md"),
        "registry_source_hash": registry.source_hash,
        "n_endpoints": len(CORRIDOR_OF), "n_corridors": len(MEMBERS),
        "frozen_at": "2026-08-22",
        "parameters": C.frozen_params(),
        "split": {
            "seed": C.SPLIT_SEED, "dev_fraction": C.DEV_FRACTION,
            "strata": per_stratum, "confusion_class": classes,
            "development": {"corridors": dev, "endpoints": dev_endpoints,
                            "n_corridors": len(dev), "n_endpoints": len(dev_endpoints)},
            "confirmatory": {"corridors": conf, "endpoints": conf_endpoints,
                             "n_corridors": len(conf), "n_endpoints": len(conf_endpoints)},
        },
    }
    import hashlib
    freeze["freeze_content_hash"] = "sha256:" + hashlib.sha256(
        canonical_json(freeze).encode()).hexdigest()
    C.CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    C.FREEZE_PATH.write_text(json.dumps(freeze, indent=1, sort_keys=True) + "\n")
    print(f"freeze written: {C.FREEZE_PATH}")
    print(f"  content hash {freeze['freeze_content_hash']}")
    print(f"  development: {len(dev)} corridors / {len(dev_endpoints)} endpoints")
    for cls, s in per_stratum.items():
        print(f"    {cls:18s} {s['n_dev']}/{s['n_corridors']}: {', '.join(s['dev'])}")
    print(f"  confirmatory: {len(conf)} corridors / {len(conf_endpoints)} endpoints")


if __name__ == "__main__":
    main()
