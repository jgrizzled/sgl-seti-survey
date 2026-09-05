"""Fetch every spectrum named in results/units_v1.json (in-window + null
draws) into runs/spectral-archives/v1/data/<instrument>/, with a
sha256 + size record per file (runs/.../fetch_manifest.jsonl).

usage: fetch.py [dev|confirmatory|all]
CARMENES: one zip per star (whole DR1 VIS set), members extracted on demand.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
RUN = REPO / "runs" / "spectral-archives" / "v1"
DATA = RUN / "data"
MANIFEST = RUN / "fetch_manifest.jsonl"
units = json.load(open(HERE / "results" / "units_v1.json"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def record(**kw):
    kw["utc"] = datetime.now(timezone.utc).isoformat()
    with open(MANIFEST, "a") as f:
        f.write(json.dumps(kw) + "\n")


def download(url, dest, tries=4, timeout=900):
    if dest.exists() and dest.stat().st_size > 1000:
        return "cached"
    tmp = dest.with_suffix(dest.suffix + ".part")
    for k in range(tries):
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            tmp.rename(dest)
            return "fetched"
        except Exception as exc:
            print(f"  retry {k+1} {url}: {exc}", file=sys.stderr, flush=True)
            time.sleep(15 * (k + 1))
    return "failed"


def main():
    fam = sys.argv[1] if len(sys.argv) > 1 else "all"
    DATA.mkdir(parents=True, exist_ok=True)
    wanted = {}   # file -> (instrument, url/zip info)
    combos_needed = set()
    for u in units["units"]:
        if fam != "all" and u["family"] != fam:
            continue
        combos_needed.add(f"{u['target']}|{u['instrument']}")
        for s in u["spectra"]:
            wanted[(u["instrument"], s["file"])] = s
    for key in combos_needed:
        c = units["combos"][key]
        for s in c["null_draw"]:
            wanted[(c["instrument"], s["file"])] = s
    print(f"{len(wanted)} files wanted for family={fam}", flush=True)

    # CARMENES zips first (one per star)
    zips = {}
    for (inst, fn), s in wanted.items():
        if "zip" in s:
            zips[s["zip"]] = s["target"]
    for zurl, tid in zips.items():
        zdest = DATA / "CARMENES-VIS" / f"{tid}_VIS.zip"
        zdest.parent.mkdir(parents=True, exist_ok=True)
        st = download(zurl, zdest, timeout=3600)
        print(f"  zip {tid}: {st} {zdest.stat().st_size if zdest.exists() else 0} bytes", flush=True)
        record(kind="zip", target=tid, url=zurl, path=str(zdest.relative_to(RUN)), status=st,
               size=zdest.stat().st_size if zdest.exists() else None, sha256=sha256(zdest) if zdest.exists() else None)
        if zdest.exists():
            with zipfile.ZipFile(zdest) as z:
                names = {Path(n).name: n for n in z.namelist()}
                for (inst, fn), s in wanted.items():
                    if s.get("zip") == zurl:
                        dest = DATA / inst / fn
                        if dest.exists():
                            continue
                        member = names.get(fn)
                        if member is None:
                            record(kind="member", target=tid, file=fn, status="missing_in_zip")
                            continue
                        dest.write_bytes(z.read(member))
                        record(kind="member", target=tid, file=fn, path=str(dest.relative_to(RUN)), status="extracted",
                               size=dest.stat().st_size, sha256=sha256(dest))

    # HTTP files in parallel
    jobs = [(inst, fn, s) for (inst, fn), s in wanted.items() if "url" in s]
    def job(inst, fn, s):
        dest = DATA / inst / fn
        dest.parent.mkdir(parents=True, exist_ok=True)
        st = download(s["url"], dest)
        ok = dest.exists()
        record(kind="file", target=s["target"], instrument=inst, file=fn, url=s["url"], status=st,
               path=str(dest.relative_to(RUN)) if ok else None, size=dest.stat().st_size if ok else None,
               sha256=sha256(dest) if ok and st == "fetched" else None)
        return fn, st
    n = {"fetched": 0, "cached": 0, "failed": 0}
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(job, *j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            fn, st = f.result()
            n[st] += 1
            if i % 25 == 0 or st == "failed":
                print(f"  [{i}/{len(jobs)}] {fn} {st}  totals {n}", flush=True)
    print("done", n, flush=True)


if __name__ == "__main__":
    main()
