"""Dev-stage data acquisition (threshold freeze v1.0 §7 order of work).

For the 5 dev units (ross-154/gj-908 S1+S2 0.1 AU C3; ross-128 S1
2.5 R_sun C2): completes the daily-listing snapshots for baseline days,
builds the frame fetch manifest per the frozen subsampling (D6), and
fetches + verifies the frames.

Frozen selection rules (config `subsampling`):
- C2 in-window days: every frame (full cadence).
- C3 in-window days: 8 frames/day — deterministic uniform pick over
  the day's time-sorted listing (implementation of "3-h bins": file
  counters are time-ordered; exact times recorded from headers at the
  search stage). Implementation note recorded in the manifest.
- Baseline days (t_ca ± 110 d, both cameras): first frame of day.

Every frame is verified (size > 500 kB, FITS magic) before acceptance;
failed fetches retry on the alternate archive. Resumable: verified
files on disk are skipped. Frames land under
runs/heliospheric-crossings/frames/{cam}/{yymmdd}/; the manifest with
per-frame source + sha256 under runs/heliospheric-crossings/dev/.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "heliospheric-crossings"
CFG = json.loads((SURV / "configs" / "threshold_freeze_v1.json").read_text())
EVENTS = REPO / "crossings" / "soho_v1" / "events.ecsv"
SNAP = REPO / "runs" / "heliospheric-crossings" / "coverage" / "listings"
FRAMES = REPO / "runs" / "heliospheric-crossings" / "frames"
DEVDIR = REPO / "runs" / "heliospheric-crossings" / "dev"
OUT = SURV / "results" / "dev_fetch_v1.json"

RSUN_AU = 0.00465047
AU_KM = 1.495978707e8
SDAC = "https://umbra.nascom.nasa.gov/pub/lasco_level05/"
NRL = "https://lasco-www.nrl.navy.mil/lz/level_05/"
SDAC_END_MJD = Time("2025-02-28").mjd
BASELINE_D = 110
C3_PER_DAY = 8

RUNG_R = {"S1_2.5Rsun": 2.5 * RSUN_AU, "S1_0.1AU": 0.1, "S2_0.1AU": 0.1}
COMBOS = {"S1": ("outbound", "target"), "S2": ("inbound", "anti_target")}


_local = __import__("threading").local()


def get(url: str, timeout: int = 60) -> bytes:
    """Keep-alive fetch: one requests.Session per worker thread —
    per-request TCP/TLS handshakes serialized the 6-worker pool at
    ~1 connection's bandwidth on the first run."""
    if not hasattr(_local, "s"):
        import requests

        _local.s = requests.Session()
        _local.s.headers["User-Agent"] = "sgl-seti-survey/1.0"
    r = _local.s.get(url, timeout=timeout)
    r.raise_for_status()
    return r.content


def listing(cam: str, day: str, mjd: float) -> list[str]:
    """Time-ordered frame names for a day (snapshot-backed)."""
    snap = SNAP / f"{cam}_{day}.html"
    if not snap.exists():
        order = [("sdac", SDAC), ("nrl", NRL)]
        if mjd > SDAC_END_MJD:
            order = order[::-1]
        html, src = "", "empty"
        for s, base in order:
            for attempt in range(3):
                try:
                    h = get(f"{base}{day}/{cam}/", 30).decode("utf-8", "replace")
                    if re.search(r'href="\d+\.fts"', h):
                        html, src = h, s
                    break
                except urllib.error.HTTPError as e:
                    if e.code == 404:
                        break
                    time.sleep(1 + attempt)
                except Exception:
                    time.sleep(1 + attempt)
            if html:
                break
        snap.write_text(f"<!--{src}-->\n" + html)
    txt = snap.read_text()
    return sorted(set(re.findall(r'href="(\d+\.fts)"', txt)))


def snap_source(cam: str, day: str) -> str:
    with open(SNAP / f"{cam}_{day}.html") as fh:
        head = fh.read(100)
    m = re.match(r"<!--(\w+)-->", head)
    return m.group(1) if m else "sdac"


def fetch_frame(cam: str, day: str, name: str) -> dict:
    dest = FRAMES / cam / day / name
    if dest.exists() and dest.stat().st_size > 500_000:
        return {"path": str(dest.relative_to(REPO)), "cached": True}
    dest.parent.mkdir(parents=True, exist_ok=True)
    src = snap_source(cam, day)
    bases = [SDAC, NRL] if src == "sdac" else [NRL, SDAC]
    err = "unknown"
    for base in bases:
        for attempt in range(3):
            try:
                data = get(f"{base}{day}/{cam}/{name}", 120)
                if len(data) > 500_000 and data[:6] == b"SIMPLE":
                    dest.write_bytes(data)
                    return {"path": str(dest.relative_to(REPO)),
                            "sha256": hashlib.sha256(data).hexdigest(),
                            "source": base[8:20], "bytes": len(data)}
                err = f"bad content {len(data)}B"
            except Exception as exc:
                err = str(exc)[:80]
            time.sleep(1 + attempt)
    return {"path": str(dest.relative_to(REPO)), "error": err}


def main() -> None:
    role = sys.argv[1] if len(sys.argv) > 1 else "dev"
    t = Table.read(EVENTS)
    dev_units = [u for u in CFG["units"] if u["role"] == role]
    # (cam, day) -> mode: "full" | "sub" | "base"  (higher mode wins)
    rank = {"base": 0, "sub": 1, "full": 2}
    daymode: dict[tuple[str, str], str] = {}
    daymjd: dict[tuple[str, str], float] = {}

    for u in dev_units:
        ch = u["unit"].split("_")[0]
        r_au = RUNG_R[u["unit"]]
        direction, side = COMBOS[ch]
        m = ((np.asarray(t["link_direction"]) == direction)
             & (np.asarray(t["side"]) == side)
             & (np.asarray(t["target_id"]) == u["target_id"])
             & (np.asarray(t["b_min_au"]) <= r_au))
        for ev in t[m]:
            b, v = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
            half = float(np.sqrt(max(r_au**2 - b**2, 0)) * AU_KM / (v * 86400.0))
            mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
            win_mode = "full" if u["camera"] == "c2" else "sub"
            for d in range(int(np.floor(mjd_ca - half)), int(np.floor(mjd_ca + half)) + 1):
                key = (u["camera"], Time(float(d), format="mjd").strftime("%y%m%d"))
                if rank[win_mode] > rank.get(daymode.get(key, "base"), -1) or key not in daymode:
                    daymode[key] = win_mode
                daymjd[key] = float(d)
            for d in range(int(np.floor(mjd_ca - BASELINE_D)), int(np.floor(mjd_ca + BASELINE_D)) + 1):
                key = (u["camera"], Time(float(d), format="mjd").strftime("%y%m%d"))
                daymode.setdefault(key, "base")
                daymjd.setdefault(key, float(d))

    print(f"{len(dev_units)} dev units -> {len(daymode)} (cam,day) entries "
          f"({sum(1 for v in daymode.values() if v == 'full')} full / "
          f"{sum(1 for v in daymode.values() if v == 'sub')} sub / "
          f"{sum(1 for v in daymode.values() if v == 'base')} base)", flush=True)

    # ensure listings, build frame list
    keys = sorted(daymode)
    with ThreadPoolExecutor(max_workers=10) as ex:
        listings = dict(zip(keys, ex.map(
            lambda k: listing(k[0], k[1], daymjd[k]), keys)))

    frames: list[tuple[str, str, str]] = []
    for key in keys:
        cam, day = key
        names = listings[key]
        if not names:
            continue
        mode = daymode[key]
        if mode == "full":
            pick = names
        elif mode == "sub":
            idx = np.unique(np.linspace(0, len(names) - 1, min(C3_PER_DAY, len(names))).astype(int))
            pick = [names[i] for i in idx]
        else:
            pick = [names[0]]
        frames.extend((cam, day, n) for n in pick)
    print(f"{len(frames)} frames to fetch (~{len(frames)*2/1024:.0f} GB)", flush=True)

    DEVDIR.mkdir(parents=True, exist_ok=True)
    results, n_err = [], 0
    with ThreadPoolExecutor(max_workers=10) as ex:
        for i, rec in enumerate(ex.map(lambda f: fetch_frame(*f), frames)):
            results.append(rec)
            n_err += "error" in rec
            if (i + 1) % 1000 == 0:
                print(f"  {i+1}/{len(frames)} ({n_err} errors)", flush=True)

    manifest = {"config_sha": (SURV / "results" / "threshold_freeze_v1_sha.txt").read_text().strip(),
                "subsampling_note": "C3 in-window = deterministic uniform 8/day pick over the "
                                    "time-ordered listing (3-h-bin intent; exact times from headers "
                                    "at search stage)",
                "day_modes": {f"{c}_{d}": m for (c, d), m in daymode.items()},
                "frames": results}
    (DEVDIR / f"fetch_manifest_{role}_v1.json").write_text(json.dumps(manifest, indent=0) + "\n")
    summary = {
        "n_days": len(daymode), "n_frames_planned": len(frames),
        "n_fetched_ok": sum(1 for r in results if "error" not in r),
        "n_cached": sum(1 for r in results if r.get("cached")),
        "n_errors": n_err,
        "errors_sample": [r for r in results if "error" in r][:20],
        "bytes_total": sum(r.get("bytes", 0) for r in results),
    }
    out_path = OUT if role == "dev" else OUT.with_name(f"{role}_fetch_v1.json")
    out_path.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "errors_sample"}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
