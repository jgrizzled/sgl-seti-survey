"""v2 engine: hypothesis freeze with the stratified corridor hold-out
(same construction as surveys/wise/scripts/freeze_v2.py). The split
is drawn per survey with the profile's seed over the corridors that
have v1 tensors, stratified by the WISE confusion class (sky density
is archive-independent), so that every survey's confirmatory set is
declared before its v2 scripts touch it."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict

import numpy as np

from sglseti import canonical_json, load_target_registry


def run(P) -> dict:
    registry = load_target_registry(P.registry_path)
    corridors = P.corridors()
    strata = defaultdict(list)
    for c in corridors:
        strata[P.confusion_class(c)].append(c)
    rng = np.random.default_rng(P.split_seed)
    forced = set(getattr(P, "forced_dev", ()) or ())
    dev, conf, per = [], [], {}
    for cls in sorted(strata):
        cs = sorted(strata[cls])
        n_dev = max(1, int(round(P.dev_fraction * len(cs)))) if len(cs) > 1 else 0
        pre = sorted(c for c in cs if c in forced)
        pool = [c for c in cs if c not in forced]
        n_draw = max(n_dev - len(pre), 0)
        drawn = (sorted(rng.choice(pool, size=n_draw, replace=False).tolist())
                 if n_draw and pool else [])
        pick = sorted(pre + drawn)
        dev.extend(pick); conf.extend(c for c in cs if c not in pick)
        per[cls] = {"n_corridors": len(cs), "n_dev": len(pick), "dev": pick,
                    **({"forced": pre} if pre else {})}
    dev, conf = sorted(dev), sorted(conf)
    dev_e = sorted(e for c in dev for e in P.members(c))
    conf_e = sorted(e for c in conf for e in P.members(c))
    freeze = {
        "survey": P.name, "hypothesis_version": P.hypothesis_version,
        "hypotheses_hash": ("sha256:" + hashlib.sha256(P.hypotheses_path.read_bytes()).hexdigest()
                            if P.hypotheses_path.exists() else None),
        "registry_source_hash": registry.source_hash,
        "n_endpoints": len(P.endpoints), "n_corridors": len(corridors),
        "frozen_at": getattr(P, "frozen_at", "2026-08-22"),
        "parameters": P.frozen_params(),
        "split": {"seed": P.split_seed, "dev_fraction": P.dev_fraction, "strata": per,
                  "confusion_class": {c: P.confusion_class(c) for c in corridors},
                  "development": {"corridors": dev, "endpoints": dev_e, "n_corridors": len(dev), "n_endpoints": len(dev_e)},
                  "confirmatory": {"corridors": conf, "endpoints": conf_e, "n_corridors": len(conf), "n_endpoints": len(conf_e)}},
    }
    freeze["freeze_content_hash"] = "sha256:" + hashlib.sha256(canonical_json(freeze).encode()).hexdigest()
    P.config_dir.mkdir(parents=True, exist_ok=True)
    P.freeze_path.write_text(json.dumps(freeze, indent=1, sort_keys=True) + "\n")
    print(f"[{P.name}] freeze {freeze['freeze_content_hash'][:23]}…: dev {len(dev)} corridors / {len(dev_e)} "
          f"endpoints; confirmatory {len(conf)} / {len(conf_e)}")
    for cls, s in per.items():
        print(f"   {cls:18s} {s['n_dev']}/{s['n_corridors']}: {', '.join(s['dev'])}")
    return freeze
