"""The pre-registered Swift/XRT look (hypotheses.md §9/D11): wolf-359 A
0.1 AU, obsid 00010119025 (2018-03-06, +2.9 d from t_ca), PC-mode cleaned
event list; 0.3-10 keV light curve in 100-s bins inside a 47" aperture
(20 XRT pixels) with a 100-200" background annulus, against the star's
2017 campaign (obsids 00010119001-024, off-window) as the null. A look,
not a trial: no threshold, no constraint. Writes results/xrt_look_v1.json/.md."""
from __future__ import annotations

import gzip
import io
import json
import re

import numpy as np
import requests
from astropy.io import fits
from astropy.time import Time

import hlib as H
import recon_scan as R

XRT_DIR = H.RUN / "xrt"
STAR = (164.1057, 7.0045)          # era-mean channel-A position of wolf-359 (recon)
R_SRC = 47.0 / 3600.0; R_BG0 = 100.0 / 3600.0; R_BG1 = 200.0 / 3600.0
BIN_S = 100.0
IN_WINDOW = "00010119025"
CAMPAIGN = [f"000101190{k:02d}" for k in range(1, 25)]
MJDREF_SWIFT = 51910.0 + 7.428703703703703e-4


def fetch_evt(obsid, ym):
    XRT_DIR.mkdir(parents=True, exist_ok=True)
    base = f"https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{ym}/{obsid}/xrt/event/"
    idx = requests.get(base, timeout=120).text
    names = sorted(set(re.findall(r'href="(sw\d+xpc\w+po_cl\.evt\.gz)"', idx)))
    out = []
    for n in names:
        p = XRT_DIR / n
        if not p.exists():
            r = requests.get(base + n, timeout=600)
            if r.status_code == 200:
                p.write_bytes(r.content)
        if p.exists():
            out.append(p)
    return out


def lightcurve(path):
    with fits.open(path) as h:
        d = h["EVENTS"].data
        gti = h["GTI"].data if "GTI" in h else h["STDGTI"].data
        hd = h["EVENTS"].header
        ra = np.asarray(d["RA"]) if "RA" in d.columns.names else None
        t = np.asarray(d["TIME"]); pi = np.asarray(d["PI"])
        # sky coordinates from X/Y if RA/DEC absent
        if ra is None:
            w = {k: hd.get(k) for k in ("TCRVL2", "TCRVL3", "TCDLT2", "TCDLT3", "TCRPX2", "TCRPX3")}
            x = np.asarray(d["X"]); y = np.asarray(d["Y"])
            dec = w["TCRVL3"] + (y - w["TCRPX3"]) * w["TCDLT3"]
            ra = w["TCRVL2"] + (x - w["TCRPX2"]) * w["TCDLT2"] / np.cos(np.radians(dec))
        else:
            dec = np.asarray(d["DEC"])
    sep = H.sep_deg_vec(ra, dec, *STAR)
    band = (pi >= 30) & (pi < 1000)        # 0.3-10 keV (PI = 10 eV)
    src = band & (sep <= R_SRC); bg = band & (sep >= R_BG0) & (sep <= R_BG1)
    area_ratio = R_SRC ** 2 / (R_BG1 ** 2 - R_BG0 ** 2)
    bins = []
    for g0, g1 in zip(gti["START"], gti["STOP"]):
        t0 = g0
        while t0 < g1:
            t1 = min(t0 + BIN_S, g1)
            if t1 - t0 >= 50.0:
                ns = int(((t[src] >= t0) & (t[src] < t1)).sum()); nb = int(((t[bg] >= t0) & (t[bg] < t1)).sum())
                bins.append({"met0": float(t0), "mjd0": float(MJDREF_SWIFT + t0 / 86400.0), "dt": float(t1 - t0),
                             "n_src": ns, "n_bg": nb, "rate": (ns - nb * area_ratio) / (t1 - t0)})
            t0 = t1
    return {"file": path.name, "n_src": int(src.sum()), "n_bg": int(bg.sum()), "expo_s": float((gti["STOP"] - gti["START"]).sum()),
            "rate": float((src.sum() - bg.sum() * area_ratio) / max((gti["STOP"] - gti["START"]).sum(), 1)), "bins": bins}


def main():
    inwin = [lightcurve(p) for p in fetch_evt(IN_WINDOW, "2018_03")]
    camp = []
    for o in CAMPAIGN:
        ym = "2017_05" if o.endswith(("01", "02")) else ("2017_06" if int(o[-2:]) <= 16 else "2017_07")
        try:
            for p in fetch_evt(o, ym):
                camp.append(lightcurve(p))
        except Exception as exc:
            print(f"[xrt] {o}: {exc}")
    # second null: the 2021-12 monitoring campaign (obsids 00014969001-084) and the 2022-06 pair
    camp2 = []
    for o in [f"000149690{k:02d}" for k in range(1, 85)]:
        try:
            for p in fetch_evt(o, "2021_12"):
                camp2.append(lightcurve(p))
        except Exception as exc:
            print(f"[xrt] {o}: {exc}")
    win = [w for w in H.unit_windows("A", "wolf-359", "0.1AU") if w["t_ca_utc"].startswith("2018-03-03")][0]
    rates_in = np.array([b["rate"] for lc in inwin for b in lc["bins"]])
    rates_c = np.array([b["rate"] for lc in camp for b in lc["bins"]])
    rates_c2 = np.array([b["rate"] for lc in camp2 for b in lc["bins"]])
    out = {"utc": Time.now().isot, "in_window": inwin, "campaign_n_files": len(camp),
           "campaign2_summary": {"n_files": len(camp2), "n_bins": int(len(rates_c2)), "median_rate": float(np.median(rates_c2)) if len(rates_c2) else None,
                                 "p95_rate": float(np.percentile(rates_c2, 95)) if len(rates_c2) else None, "max_rate": float(rates_c2.max()) if len(rates_c2) else None,
                                 "obs_rates": [(lc["file"], lc["rate"], lc["expo_s"]) for lc in camp2]},
           "campaign_summary": {"n_bins": int(len(rates_c)), "median_rate": float(np.median(rates_c)) if len(rates_c) else None,
                                "p95_rate": float(np.percentile(rates_c, 95)) if len(rates_c) else None,
                                "max_rate": float(rates_c.max()) if len(rates_c) else None,
                                "obs_rates": [(lc["file"], lc["rate"], lc["expo_s"]) for lc in camp]},
           "in_window_summary": {"n_bins": int(len(rates_in)), "median_rate": float(np.median(rates_in)) if len(rates_in) else None,
                                 "max_rate": float(rates_in.max()) if len(rates_in) else None,
                                 "bins_above_campaign_p95": int((rates_in > np.percentile(rates_c, 95)).sum()) if len(rates_c) and len(rates_in) else None,
                                 "window_utc": [Time(win["mjd0"], format="mjd").isot, Time(win["mjd1"], format="mjd").isot]}}
    # flare morphology: contiguous bins above 3x the in-window median
    if len(rates_in):
        med = np.median(rates_in); hi = rates_in > 3 * med
        runs = []
        k = 0
        while k < len(hi):
            if hi[k]:
                j = k
                while j < len(hi) and hi[j]:
                    j += 1
                runs.append({"n_bins": j - k, "peak_rate": float(rates_in[k:j].max()), "peak_over_median": float(rates_in[k:j].max() / max(med, 1e-9))})
                k = j
            else:
                k += 1
        out["in_window_summary"]["excursions_3x_median"] = runs
    (H.RES / "xrt_look_v1.json").write_text(json.dumps(H.clean(out), indent=1))
    s = out["in_window_summary"]; c = out["campaign_summary"]
    md = ["# Swift/XRT look v1 — wolf-359 A 0.1 AU, obsid 00010119025 (2018-03-06)", "",
          f"Window {s['window_utc'][0][:16]} → {s['window_utc'][1][:16]} UTC; PC-mode cleaned events; 47″ aperture, 100–200″ annulus, 0.3–10 keV, 100-s bins.", "",
          f"- In-window: {len(inwin)} event file(s), {sum(l['n_src'] for l in inwin)} source counts in {sum(l['expo_s'] for l in inwin):.0f} s "
          f"(mean net rate {np.mean([l['rate'] for l in inwin]) if inwin else float('nan'):.3f} ct/s); {s['n_bins']} bins, median {s['median_rate']:.3f}, max {s['max_rate']:.3f} ct/s; "
          f"{s['bins_above_campaign_p95']} bins above the campaign's 95th percentile; excursions > 3× median: {s.get('excursions_3x_median')}",
          f"- 2017 campaign null: {c['n_bins']} bins over {len(camp)} files; median {c['median_rate']:.3f}, p95 {c['p95_rate']:.3f}, max {c['max_rate']:.3f} ct/s.",
          f"- 2021-12 campaign null: {out['campaign2_summary']['n_bins']} bins over {len(camp2)} files; median {out['campaign2_summary']['median_rate']:.3f}, "
          f"p95 {out['campaign2_summary']['p95_rate']:.3f}, max {out['campaign2_summary']['max_rate']:.3f} ct/s; per-visit mean rates: "
          + ", ".join(f"{r[1]:.3f}" for r in sorted(out['campaign2_summary']['obs_rates'])[:84]),
          f"- 2017 per-visit mean rates: " + ", ".join(f"{r[1]:.3f}" for r in sorted(c['obs_rates'])),
          "", "A look, not a trial (hypotheses §9): no threshold, no constraint claimed."]
    (H.RES / "xrt_look_v1.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
