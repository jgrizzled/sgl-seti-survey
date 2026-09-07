"""Summarise the recon_scan.py stage outputs into results/recon_summary_v0.md
(tables quoted by notes/highenergy_recon_2026-09-07.md)."""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.time import Time

RES = Path(__file__).resolve().parents[1] / "results"
TODAY_MJD = Time("2026-09-07", format="iso", scale="utc").mjd
RUNGS = ("1.2Rsun", "2.5Rsun", "0.1AU")


def load(name):
    p = RES / name
    return json.loads(p.read_text()) if p.exists() else None


def past(rec):
    return Time(rec["t_ca_utc"], format="isot", scale="utc").mjd <= TODAY_MJD


def main():
    out = ["# High-energy recon — summary tables (v0, 2026-09-07)", ""]

    m = load("recon_masters_v0.json")
    if m:
        out += ["## Pointed X-ray missions (HEASARC master tables, era-mean channel positions)", "",
                "| Channel | Target | Chandra (25′) | XMM pointed (20′) | Swift (15′) |", "|---|---|---|---|---|"]
        byk = defaultdict(dict)
        for c in m["cones"]:
            byk[(c["channel"], c["target"])][c["table"]] = c["n_rows"]
        for (ch, tid), d in sorted(byk.items()):
            out.append(f"| {ch} | {tid} | {d.get('chanmaster')} | {d.get('xmmmaster')} | {d.get('swiftmastr')} |")
        out += ["", f"In-window pointed units: **{len(m['units'])}**", ""]
        for u in m["units"]:
            out.append(f"- {u['channel']} {u['target']} {u['rung']} t_ca {u['t_ca_utc'][:16]} (b {u['b_min_rsun']:.2f} R☉): "
                       f"{u['table']} obsid {u['obsid']} ({u['name']}) {u['obs_start_utc'][:16]} → {u['obs_stop_utc'][11:16]}, "
                       f"{u['offset_days']:+.1f} d from t_ca, XRT {u['exposure_s']:.0f} s, pointing offset {u['pointing_offset_arcmin']:.1f}′")
        out.append("")

    s = load("recon_slew_v0.json")
    if s:
        out += ["## XMM-Newton slew survey (XSA `v_slew_exposure` footprint CONTAINS)", "",
                "| Channel | Target | Slew exposures over the position | Nearest to any t_ca (d) |", "|---|---|---|---|"]
        for p in s["positions"]:
            nd = p["nearest_days"]
            out.append(f"| {p['channel']} | {p['target']} | {p['n_rows']} | {'—' if nd is None else f'{nd:+.0f}'} |")
        out += ["", f"In-window slew units: **{len(s['units'])}**", ""]

    b = load("recon_bat_v0.json")
    if b:
        out += ["## Swift/BAT coded-mask coverage of the windows (swiftmastr per-window queries, survey+event exposure, pointing offset ≤ 30° / 20°)", "",
                "| Channel | Rung | Windows (t_ca ≤ today) | ≥ 1 ks within 30° | ≥ 1 ks within 20° | median s (30°) | p90 s (30°) |",
                "|---|---|---|---|---|---|---|"]
        for ch in "AB":
            for rung in RUNGS:
                xs = [e["rungs"][rung] for e in b["events"] if e["channel"] == ch and rung in e["rungs"] and past(e)]
                if not xs:
                    continue
                s30 = np.array([x["bat_s_30deg"] for x in xs]); s20 = np.array([x["bat_s_20deg"] for x in xs])
                out.append(f"| {ch} | {rung} | {len(xs)} | {(s30 >= 1000).sum()} | {(s20 >= 1000).sum()} | {np.median(s30):.0f} | {np.percentile(s30, 90):.0f} |")
        # per-target grazing detail
        out += ["", "Grazing-family (≤ 2.5 R☉) windows with ≥ 1 ks BAT within 20° (the best-coded subset), per target-channel:", ""]
        cnt = Counter(); tot = Counter()
        for e in b["events"]:
            if not past(e) or "2.5Rsun" not in e["rungs"]:
                continue
            tot[(e["channel"], e["target"])] += 1
            if e["rungs"]["2.5Rsun"]["bat_s_20deg"] >= 1000:
                cnt[(e["channel"], e["target"])] += 1
        for k in sorted(tot):
            out.append(f"- {k[0]} {k[1]}: {cnt[k]}/{tot[k]}")
        out.append("")

    l = load("recon_lat_v0.json")
    if l:
        out += [f"## Fermi-LAT in-FoV livetime in the windows (weekly spacecraft files; θ ≤ {l['theta_max_deg']:.0f}°, zenith ≤ {l['zenith_max_deg']:.0f}°, DATA_QUAL > 0, LAT_CONFIG = 1)", "",
                "| Channel | Rung | Windows | livetime min / median / max (ks) | duty median | windows with < 1 ks |",
                "|---|---|---|---|---|---|"]
        for ch in "AB":
            for rung in RUNGS:
                xs = [e["rungs"][rung] for e in l["events"] if e["channel"] == ch and rung in e["rungs"] and past(e)]
                if not xs:
                    continue
                lt = np.array([x["livetime_s"] for x in xs]) / 1e3
                duty = np.array([x["duty"] for x in xs])
                out.append(f"| {ch} | {rung} | {len(xs)} | {lt.min():.1f} / {np.median(lt):.1f} / {lt.max():.1f} | {np.median(duty):.2f} | {(lt < 1).sum()} |")
        out += ["", "Split at the 2018-03-16 solar-array-drive anomaly (the post-anomaly survey profile):", "",
                "| Era | Channel | Rung | Windows | < 1 ks | median ks |", "|---|---|---|---|---|---|"]
        T = Time("2018-03-16", format="iso", scale="utc").mjd
        for era, (lo, hi) in (("pre", (0, T)), ("post", (T, TODAY_MJD))):
            for ch in "AB":
                for rung in RUNGS:
                    xs = [e["rungs"][rung] for e in l["events"] if e["channel"] == ch and rung in e["rungs"]
                          and lo <= Time(e["t_ca_utc"], format="isot", scale="utc").mjd <= hi]
                    if not xs:
                        continue
                    lt = np.array([x["livetime_s"] for x in xs]) / 1e3
                    out.append(f"| {era} | {ch} | {rung} | {len(xs)} | {(lt < 1).sum()} | {np.median(lt):.1f} |")
        zero = [e for e in l["events"] if past(e) and any(v["livetime_s"] < 1000 for v in e["rungs"].values())]
        if zero:
            out += ["", "Windows with < 1 ks in-FoV livetime at some rung:", ""]
            for e in zero:
                out.append(f"- {e['channel']} {e['target']} {e['t_ca_utc'][:10]}: " + ", ".join(f"{k} {v['livetime_s']:.0f} s ({v['n_weeks_found']}/{len(v['weeks'])} weeks)" for k, v in e["rungs"].items()))
        out.append("")

    l1 = load("recon_lat01_v0.json")
    if l1:
        out += ["## Fermi-LAT in-FoV livetime in the 0.1 AU windows, all seven targets (`recon_lat01_v0.json`)", "",
                "| Era | Channel | Windows | < 1 ks | < 10 ks | min / median / max (ks) |", "|---|---|---|---|---|---|"]
        T = Time("2018-03-16", format="iso", scale="utc").mjd
        for era, (lo, hi) in (("all", (0, TODAY_MJD)), ("pre-2018-03", (0, T)), ("post-2018-03", (T, TODAY_MJD))):
            for ch in "AB":
                xs = [e for e in l1["events"] if e["channel"] == ch and lo <= Time(e["t_ca_utc"], format="isot", scale="utc").mjd <= hi]
                lt = np.array([e["rungs"]["0.1AU"]["livetime_s"] for e in xs]) / 1e3
                out.append(f"| {era} | {ch} | {len(xs)} | {(lt < 1).sum()} | {(lt < 10).sum()} | {lt.min():.1f} / {np.median(lt):.1f} / {lt.max():.1f} |")
        low = [e for e in l1["events"] if past(e) and e["rungs"]["0.1AU"]["livetime_s"] < 10000]
        out += ["", "0.1 AU windows below 10 ks: " + "; ".join(f"{e['channel']} {e['target']} {e['t_ca_utc'][:10]} ({e['rungs']['0.1AU']['livetime_s']:.0f} s)" for e in low), ""]

    c = load("recon_catalogs_v0.json")
    if c:
        out += ["## Catalogue screens on the 88 anti-star corridors (7′) and the 7 deep-family stars (2′; 4FGL 30′, BAT 12′)", "",
                "| Catalogue | Corridors with ≥ 1 source | Sources total (corridors) | Stars detected | Errors |", "|---|---|---|---|---|"]
        tabs = sorted(set(x["table"] for x in c["cones"]))
        for t in tabs:
            cs = [x for x in c["cones"] if x["table"] == t and x["kind"] == "corridor"]
            ss = [x for x in c["cones"] if x["table"] == t and x["kind"] == "star"]
            ncor = sum(1 for x in cs if x["n_rows"]); nsrc = sum(x["n_rows"] or 0 for x in cs)
            stars = [x["target"] for x in ss if x["n_rows"]]
            errs = sum(1 for x in cs + ss if x["error"])
            out.append(f"| {t} | {ncor}/{len(cs)} | {nsrc} | {', '.join(stars) or '—'} | {errs} |")
        cs = [x for x in c["csc"] if x["kind"] == "corridor"]; ss = [x for x in c["csc"] if x["kind"] == "star"]
        out.append(f"| CSC 2.1 (csc21tap) | {sum(1 for x in cs if x['n_rows'])}/{len(cs)} | {sum(x['n_rows'] or 0 for x in cs)} | {', '.join(x['target'] for x in ss if x['n_rows']) or '—'} | {sum(1 for x in cs + ss if x['error'])} |")
        er = load("recon_erodat_rows_v0.json")
        escs = er["cones"] if er else c["erosita_scs"]
        for cat in ("DR1_Main", "DR2_Main"):
            cs = [x for x in escs if x["cat"] == cat and x["kind"] == "corridor"]
            ss = [x for x in escs if x["cat"] == cat and x["kind"] == "star"]
            out.append(f"| eRODat {cat} | {sum(1 for x in cs if x['n_rows'])}/{len(cs)} | {sum(x['n_rows'] or 0 for x in cs)} | {', '.join(x['target'] for x in ss if x['n_rows']) or '—'} | {sum(1 for x in cs + ss if x.get('error'))} |")
        out += ["", "Corridor matches (any catalogue except 4FGL/BAT wide cones), with separation from the corridor centre:", ""]
        for x in c["cones"] + c["csc"]:
            if x["kind"] != "corridor" or not x["n_rows"] or x.get("table") in ("fermilpsc", "swbat157m"):
                continue
            for r in x["rows"]:
                nm = r.get("name") or r.get("srcid") or r.get("detid") or "?"
                out.append(f"- {x['target']} {x.get('table', 'csc21')}: {nm} at {r.get('_sep_arcmin', float('nan')) if '_sep_arcmin' in r else '?'}′")
        if er:
            for x in er["cones"]:
                if x["kind"] == "corridor" and x["n_rows"]:
                    for r in x["rows"]:
                        out.append(f"- {x['target']} eRODat {x['cat']}: {r.get('iauname')} at {r['_sep_arcmin']:.1f}′ (det_like {r.get('det_like_0', '?')})")
        for x in c["cones"]:
            if x["kind"] == "corridor" and x["n_rows"] and x["table"] in ("fermilpsc", "swbat157m"):
                for r in x["rows"]:
                    out.append(f"- {x['target']} {x['table']} (wide): {r.get('name')} at {r['_sep_arcmin']:.1f}′")
        ul = c.get("erosita_ul") or {}
        rows = ul.get("rows") or []
        if rows:
            out += ["", "### eRODat upper limits (multi API)", ""]
            for drs, band in (("DR1_eRASS1", "024"), ("DR1_eRASS1", "021"), ("DR2_eRASSc3", "024")):
                rr = [r for r in rows if r.get("dr_survey") == drs and r.get("band") == band and r["t"]["kind"] == "corridor"]
                de = [r for r in rr if r.get("de_sky")]
                ulb = np.array([r["UL_B"] for r in de if r.get("UL_B") is not None])
                ex = np.array([r["Exposure"] for r in de if r.get("Exposure") is not None])
                cnt = np.array([r["Count"] for r in de if r.get("Count") is not None])
                cnt_s = f"; counts median {np.median(cnt):.0f}, max {cnt.max():.0f}" if len(cnt) else " (no counts field)"
                out.append(f"- {drs} band {band}: {len(de)}/{len(rr)} corridors in the DE sky; exposure median {np.median(ex):.0f} s "
                           f"(min {ex.min():.0f}, max {ex.max():.0f}); UL_B median {np.median(ulb):.2e} erg/cm²/s (min {ulb.min():.1e}, max {ulb.max():.1e})" + cnt_s)
            de_list = sorted(set(r["t"]["target"] for r in rows if r["t"]["kind"] == "corridor" and r.get("de_sky")))
            out.append(f"- DE-sky corridors ({len(de_list)}): {', '.join(de_list)}")
            st = [r for r in rows if r["t"]["kind"] == "star" and r.get("dr_survey") == "DR1_eRASS1" and r.get("band") == "024"]
            out.append("- Deep-family stars (DR1 eRASS1 band 024): " + "; ".join(
                f"{r['t']['target']} {'DE' if r.get('de_sky') else 'RU'}" + (f" {r['Count']} ct / {r['Exposure']:.0f} s" if r.get("de_sky") else "") for r in st))
        out.append("")

    r = load("recon_routes_v0.json")
    if r:
        out += ["## Data-product routes", "", "| Probe | HTTP | Content |", "|---|---|---|"]
        for p in r["probes"]:
            out.append(f"| {p['label']} | {p.get('http', p.get('error'))} | {p.get('content_type', '')} {p.get('content_length', '')} |")
        out.append("")

    (RES / "recon_summary_v0.md").write_text("\n".join(out))
    print("\n".join(out))


if __name__ == "__main__":
    main()
