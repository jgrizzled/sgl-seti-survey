"""Select the positive-control asteroid exposure set (freeze decision
#6): (60000) Miminko — the PS1 control, for cross-survey continuity —
through DECam instcal exposures found by CADC SSOIS.

Steps: parse the SSOIS listing (snapshotted), pick the deepest
instcal exposure per night, resolve EXPNUM + pointing from the file's
primary header (``?hdus=0``), re-discover the exposure triplet through
the survey adapter's own box search (EXPNUM must match — validates the
adapter path end-to-end), and pull the predicted V magnitude and
position from JPL Horizons at each exposure epoch. Writes
surveys/decam/configs/asteroid_control_v1.json.

Usage: uv run python surveys/decam/scripts/asteroid_control_recon.py
"""

from __future__ import annotations

import io
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import requests

from sglsurvey.adapters.base import ConeRegion, MjdRange
from sglsurvey.adapters.noirlab_decam import DecamInstcalAdapter
from sglsurvey.snapshots import SnapshotStore

REPO = Path(__file__).resolve().parents[3]
RUN_DIR = REPO / "runs" / "decam" / "asteroid_control_v1"
OUT = (REPO / "surveys" / "decam" / "configs"
       / "asteroid_control_v1.json")
ASTEROID = "60000"          # (60000) Miminko, the PS1 control
SSOIS_URL = ("https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/"
             "ssos/ssosclf.pl")
HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"
MIN_EXPTIME_S = 30.0
MAX_NIGHTS = 99   # all qualifying exposures: rescore needs >= 5 epochs per band


def fetch_ssois(store: SnapshotStore) -> list[list[str]]:
    params = ("lang=en;object=" + ASTEROID + ";search=bynameall;"
              "epoch1=2013+01+01;epoch2=2024+12+31;eunits=none;"
              "extres=no;xyres=no;format=tsv")
    url = SSOIS_URL + "?" + params
    request_utc = datetime.now(timezone.utc).isoformat()
    resp = requests.get(url, timeout=180)
    resp.raise_for_status()
    store.store(service_url=SSOIS_URL, query=params,
                request_utc=request_utc, response_bytes=resp.content,
                row_count=resp.text.count("\n"), http_status=200)
    return [l.split("\t") for l in resp.text.splitlines()[1:]]


def horizons_vmag(mjds: list[float], store: SnapshotStore) -> dict:
    tlist = " ".join(f"'{m + 2400000.5:.6f}'" for m in mjds)
    params = {"format": "json", "COMMAND": f"'{ASTEROID};'",
              "EPHEM_TYPE": "OBSERVER", "CENTER": "'W84'",
              "TLIST": tlist, "TLIST_TYPE": "JD",
              "QUANTITIES": "'1,9'", "CSV_FORMAT": "YES"}
    request_utc = datetime.now(timezone.utc).isoformat()
    resp = requests.get(HORIZONS_URL, params=params, timeout=120)
    resp.raise_for_status()
    store.store(service_url=HORIZONS_URL,
                query=urllib.parse.urlencode(params),
                request_utc=request_utc, response_bytes=resp.content,
                row_count=len(mjds), http_status=200)
    text = resp.json()["result"]
    lines = text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    out = {}
    for mjd, line in zip(mjds, lines):
        # CSV: date, sun-flag, moon-flag, RA " H M S", Dec "sD M S",
        # APmag, S-brt (QUANTITIES 1,9)
        parts = [p.strip() for p in line.split(",")]
        try:
            h, m, sec = (float(v) for v in parts[3].split())
            dd = parts[4].split()
            sign = -1.0 if dd[0].startswith("-") else 1.0
            dec = sign * (abs(float(dd[0])) + float(dd[1]) / 60
                          + float(dd[2]) / 3600)
            out[round(mjd, 5)] = {
                "ra_deg": round(15.0 * (h + m / 60 + sec / 3600), 6),
                "dec_deg": round(dec, 6),
                "vmag": float(parts[5])}
        except (IndexError, ValueError):
            out[round(mjd, 5)] = {"raw": line}
    return out


def main() -> None:
    from astropy.io import fits

    store = SnapshotStore(RUN_DIR)
    adapter = DecamInstcalAdapter()
    rows = fetch_ssois(store)
    dec = [r for r in rows if "DECam" in r[7] and "_ooi_" in r[0]
           and float(r[3]) >= MIN_EXPTIME_S]
    nights = defaultdict(list)
    for r in dec:
        nights[int(float(r[1]))].append(r)
    print(f"{len(dec)} instcal rows over {len(nights)} nights")

    picks = []
    for n in sorted(nights):
        for r in sorted(nights[n], key=lambda r: -float(r[3])):
            md5 = r[0].split()[-1]
            hdr = None
            for attempt in (1, 2, 3):
                try:
                    raw = urllib.request.urlopen(
                        f"https://astroarchive.noirlab.edu/api/retrieve/"
                        f"{md5}/?hdus=0", timeout=180).read()
                    hdr = fits.open(io.BytesIO(raw))[0].header
                    break
                except Exception as exc:
                    print(f"  night {n}: header attempt {attempt} "
                          f"failed ({type(exc).__name__})", flush=True)
                    time.sleep(3 * attempt)
            if hdr is None:
                continue
            expnum = int(hdr["EXPNUM"])
            # box around the SSOIS object position (inside the FOV by
            # construction) — header pointing keywords vary by CP era
            mjd = float(r[1])
            obs = [o for o in adapter.discover(
                ConeRegion(ra_deg=float(r[4]), dec_deg=float(r[5]),
                           radius_deg=0.05),
                MjdRange(mjd - 0.2, mjd + 0.2), store)
                if o.native_key["expnum"] == expnum]
            if not obs:
                print(f"  night {n}: EXPNUM {expnum} triplet not found "
                      f"via adapter", flush=True)
                continue
            o = obs[0]
            picks.append({
                "night_mjd": n, "expnum": expnum,
                "observation_id": o.observation_id, "band": o.band,
                "exptime_s": o.exptime_s, "mjd_mid": o.t_mid_mjd_utc,
                "object_ra_deg": float(r[4]),
                "object_dec_deg": float(r[5]),
                "proposal": o.extra.get("proposal"),
                "products": {k: v["md5"] for k, v in o.products.items()}})
            print(f"  night {n}: EXPNUM {expnum} {o.band} "
                  f"{o.exptime_s:.0f}s OK ({o.observation_id})",
                  flush=True)
        if len(picks) >= MAX_NIGHTS:
            break

    eph = horizons_vmag([p["mjd_mid"] for p in picks], store)
    for p in picks:
        p["horizons"] = eph.get(round(p["mjd_mid"], 5))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "version": "asteroid_control_v1",
        "asteroid": "(60000) Miminko",
        "why": ("PS1 positive control reused for cross-survey "
                "continuity; SSOIS-listed DECam instcal coverage over "
                "2013-2023 in griz"),
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "selection": ("deepest instcal exposure per night, EXPTIME >= "
                      f"{MIN_EXPTIME_S:.0f} s, triplet re-discovered "
                      "through DecamInstcalAdapter"),
        "exposures": picks}, indent=2))
    bands = defaultdict(int)
    for p in picks:
        bands[p["band"]] += 1
    vs = [p["horizons"]["vmag"] for p in picks
          if p.get("horizons") and "vmag" in p["horizons"]]
    print(f"\n{len(picks)} control exposures over {len(set(p['night_mjd'] for p in picks))} "
          f"nights, bands {dict(bands)}, V {min(vs):.1f}-{max(vs):.1f}"
          if vs else f"\n{len(picks)} control exposures (no V mags)")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
