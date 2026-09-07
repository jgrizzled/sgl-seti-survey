"""Corridor catalogue screen v1 (corridor_screen.md): hits within 7' of the
88 anti-star corridor centres from the recon catalogue pulls, the SGL
track test per hit (two-epoch rule with 5XMM detections and the
eRASS1/eRASS:3 pair; stack rule with LSXPS spans and eRASS:3 compactness),
and the eRODat upper-limit ledger. Writes results/corridor_screen_v1.json/.md."""
from __future__ import annotations

import json
from collections import Counter, defaultdict

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

import hlib as H
import recon_scan as R

Z_MIN_AU = 10000.0
PAR_MIN_ARCSEC = 206265.0 / Z_MIN_AU
MATCH_ARCSEC = 10.0
DISP_MIN_ARCSEC = 15.0
ERASS1_MID = Time("2020-03-15").mjd
ERASS_SCAN2_MID = Time("2020-09-15").mjd
ERASS_SCAN3_MID = Time("2021-03-15").mjd
EROSITA_HEW = 30.0
XRT_HEW = 18.0


def earth_au(mjd):
    t = Time(mjd, format="mjd", scale="utc")
    e = get_body_barycentric("earth", t).xyz.to_value("AU")
    s = get_body_barycentric("sun", t).xyz.to_value("AU")
    return np.asarray(e - s)


def relay_displacement_arcsec(ra, dec, mjd1, mjd2, z_au=Z_MIN_AU):
    """Predicted apparent displacement of a relay at z_au on the anti-star
    axis toward (ra, dec) between two epochs (perpendicular projection)."""
    u = np.array([np.cos(np.radians(dec)) * np.cos(np.radians(ra)), np.cos(np.radians(dec)) * np.sin(np.radians(ra)), np.sin(np.radians(dec))])
    d = earth_au(mjd1) - earth_au(mjd2)
    d_perp = d - u * (d @ u)
    return 206265.0 / z_au * float(np.linalg.norm(d_perp))


def epoch_supplement(cat):
    """Supplementary HEASARC pull (snapshotted): first/last detection epochs
    for every 5XMM (xmmssc time/end_time), LSXPS and 2SXPS hit — the recon
    stored only one of the two. Cached under runs/.../v1/catalogs."""
    cache = H.RUN / "catalogs" / "epoch_supplement.json"
    if cache.exists():
        return json.loads(cache.read_text())
    from sglsurvey.snapshots import SnapshotStore
    store = SnapshotStore(H.RUN / "catalogs")
    out = {}
    for table, idcol, cols in (("xmmssc", "srcid", 'srcid, "time", end_time'),
                               ("swiftlsxps", "name", 'name, "time", end_time, exposure'),
                               ("swift2sxps", "name", 'name, "time", end_time, exposure')):
        ids = sorted(set(str(r.get(idcol)) for x in cat["cones"] if x["table"] == table and x["kind"] == "corridor" and x["n_rows"] for r in x["rows"]))
        for i in range(0, len(ids), 40):
            chunk = ids[i:i + 40]
            adql = f"SELECT {cols} FROM {table} WHERE " + " OR ".join(f"{idcol}='{v}'" for v in chunk)
            rows, err = R.heasarc(store, adql, f"epochs {table} {i}")
            for r in rows or []:
                out[f"{table}:{r[idcol]}"] = {"time": R.ffloat(r.get("time")), "end_time": R.ffloat(r.get("end_time"))}
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(out, indent=1))
    return out


def main():
    cat = json.loads((H.RES / "recon_catalogs_v0.json").read_text())
    epochs_sup = epoch_supplement(cat)
    ero = json.loads((H.RES / "recon_erodat_rows_v0.json").read_text())
    ul = cat["erosita_ul"]["rows"]
    corridors = sorted(set(x["target"] for x in cat["cones"] if x["kind"] == "corridor"))
    pos = {x["target"]: (x["ra"], x["dec"]) for x in cat["cones"] if x["kind"] == "corridor"}
    out = {"utc": Time.now().isot, "rule": {"z_min_au": Z_MIN_AU, "match_arcsec": MATCH_ARCSEC, "disp_min_arcsec": DISP_MIN_ARCSEC}, "corridors": []}
    census = Counter(); n_hits = 0
    for tid in corridors:
        ra, de = pos[tid]
        rec = {"target": tid, "ra": ra, "dec": de, "hits": [], "upper_limits": []}
        # upper limits
        for r in ul:
            if r["t"]["kind"] == "corridor" and r["t"]["target"] == tid:
                rec["upper_limits"].append({k: r.get(k) for k in ("dr_survey", "band", "de_sky", "Exposure", "UL_B", "Count", "Bkg_counts")})
        rec["de_sky"] = any(r.get("de_sky") for r in rec["upper_limits"])
        # --- collect hits (position, catalogue, epoch info)
        hits = []
        for x in cat["cones"]:
            if x["kind"] != "corridor" or x["target"] != tid or not x["n_rows"] or x["table"] in ("fermilpsc", "swbat157m"):
                continue
            for r in x["rows"]:
                sup = epochs_sup.get(f"{x['table']}:{r.get('name') or r.get('srcid')}", {})
                hits.append({"cat": x["table"], "name": r.get("name") or r.get("srcid"), "ra": R.ffloat(r["ra"]), "dec": R.ffloat(r["dec"]),
                             "sep_arcmin": r["_sep_arcmin"], "epoch_mjd": sup.get("time") or R.ffloat(r.get("time")),
                             "end_mjd": sup.get("end_time") or R.ffloat(r.get("end_time")),
                             "err_arcsec": R.ffloat(r.get("error_radius"))})
        for x in cat["csc"]:
            if x["kind"] == "corridor" and x["target"] == tid and x["n_rows"]:
                for r in x["rows"]:
                    hits.append({"cat": "csc21", "name": r["name"], "ra": float(r["ra"]), "dec": float(r["dec"]), "sep_arcmin": r["_sep_arcmin"],
                                 "epoch_mjd": None, "end_mjd": None, "err_arcsec": R.ffloat(r.get("err_ellipse_r0"))})
        for x in ero["cones"]:
            if x["kind"] == "corridor" and x["target"] == tid and x["n_rows"]:
                for r in x["rows"]:
                    hits.append({"cat": "erodat_" + x["cat"], "name": r["iauname"], "ra": float(r["ra"]), "dec": float(r["dec"]), "sep_arcmin": r["_sep_arcmin"],
                                 "epoch_mjd": ERASS1_MID if x["cat"] == "DR1_Main" else None, "end_mjd": None,
                                 "err_arcsec": R.ffloat(r.get("radec_err")), "ext": R.ffloat(r.get("ext")), "det_like": R.ffloat(r.get("det_like_0"))})
        # --- group hits into sources (10" clustering) and apply the track test
        groups = []
        for h in sorted(hits, key=lambda h: h["sep_arcmin"]):
            for g in groups:
                if R.sep_deg(g["ra"], g["dec"], h["ra"], h["dec"]) * 3600 <= MATCH_ARCSEC:
                    g["members"].append(h); break
            else:
                groups.append({"ra": h["ra"], "dec": h["dec"], "members": [h]})
        for g in groups:
            m = g["members"]
            cats = sorted(set(x["cat"] for x in m))
            disp = None; reason = None
            # rule 1: two epochs with a decisive predicted displacement
            epochs = [(x["epoch_mjd"], x["ra"], x["dec"]) for x in m if x["epoch_mjd"]]
            # merged multi-epoch catalogue sources (5XMM, LSXPS, 2SXPS): the
            # catalogue merged the first and last detections into one position
            for x in m:
                if x["cat"] in ("xmmssc", "swiftlsxps", "swift2sxps") and x["epoch_mjd"] and x["end_mjd"] and x["end_mjd"] > x["epoch_mjd"] + 1:
                    epochs.append((x["end_mjd"], x["ra"], x["dec"]))
            # eRASS:3 counterpart supplies the scan-2/3 epochs for an eRASS1 detection
            if any(x["cat"] == "erodat_DR1_Main" for x in m) and any(x["cat"] == "erodat_DR2_Main" for x in m):
                epochs.append((ERASS_SCAN2_MID, *[(x["ra"], x["dec"]) for x in m if x["cat"] == "erodat_DR2_Main"][0]))
            pairs = []
            for i in range(len(epochs)):
                for j in range(i + 1, len(epochs)):
                    d_pred = relay_displacement_arcsec(ra, de, epochs[i][0], epochs[j][0])
                    d_obs = R.sep_deg(epochs[i][1], epochs[i][2], epochs[j][1], epochs[j][2]) * 3600
                    pairs.append((d_pred, d_obs, epochs[i][0], epochs[j][0]))
            decisive = [p for p in pairs if p[0] >= DISP_MIN_ARCSEC and p[1] < MATCH_ARCSEC]
            if decisive:
                disp = "fixed_sky"; reason = f"two-epoch: {len(decisive)} decisive pairs, e.g. predicted {decisive[0][0]:.0f}\" vs observed {decisive[0][1]:.1f}\""
            else:
                # rule 2: stack compactness
                ls = [x for x in m if x["cat"] == "swiftlsxps" and x["epoch_mjd"] and x["end_mjd"] and (x["end_mjd"] - x["epoch_mjd"]) >= 180]
                e3 = [x for x in m if x["cat"] == "erodat_DR2_Main" and (x.get("ext") or 0) == 0]
                if ls:
                    disp = "fixed_sky"; reason = f"stack: LSXPS compact over {ls[0]['end_mjd'] - ls[0]['epoch_mjd']:.0f} d"
                elif e3:
                    disp = "fixed_sky"; reason = "stack: eRASS:3 compact (ext = 0) over three scans"
                else:
                    disp = "single_epoch_unresolved"; reason = "no second epoch / no multi-scan stack"
            g["catalogues"] = cats; g["disposition"] = disp; g["reason"] = reason
            g["sep_arcmin"] = min(x["sep_arcmin"] for x in m); g["n_pairs"] = len(pairs)
            g["names"] = sorted(set(str(x["name"]) for x in m))[:4]
            census[disp] += 1; n_hits += 1
        rec["hits"] = [{k: v for k, v in g.items() if k != "members"} | {"n_members": len(g["members"])} for g in groups]
        rec["n_hits"] = len(groups)
        rec["status"] = "no_source_within_7arcmin" if not groups else ("all_fixed_sky" if all(g["disposition"] == "fixed_sky" for g in groups) else "has_unresolved")
        out["corridors"].append(rec)
    out["census"] = {"corridors": len(corridors), "corridors_with_hits": sum(1 for c in out["corridors"] if c["n_hits"]),
                     "sources": n_hits, "dispositions": dict(census),
                     "de_sky_corridors": sum(1 for c in out["corridors"] if c["de_sky"]),
                     "status": dict(Counter(c["status"] for c in out["corridors"]))}
    (H.RES / "corridor_screen_v1.json").write_text(json.dumps(H.clean(out), indent=1))
    md = ["# Corridor catalogue screen v1 — 2026-09-07", "", f"Census: {json.dumps(out['census'])}", "",
          "| Corridor | DE sky | eRASS1 UL 024 (erg/cm²/s) | eRASS:3 UL 024 | sources ≤ 7′ | fixed_sky | unresolved | closest (′) |", "|---|---|---|---|---|---|---|---|"]
    for c in out["corridors"]:
        u1 = [u for u in c["upper_limits"] if u["dr_survey"] == "DR1_eRASS1" and u["band"] == "024" and u.get("UL_B")]
        u3 = [u for u in c["upper_limits"] if u["dr_survey"] == "DR2_eRASSc3" and u.get("UL_B")]
        fx = sum(1 for g in c["hits"] if g["disposition"] == "fixed_sky"); un = c["n_hits"] - fx
        cl = min([g["sep_arcmin"] for g in c["hits"]], default=None)
        md.append(f"| {c['target']} | {'yes' if c['de_sky'] else 'no'} | {u1[0]['UL_B']:.1e} | {u3[0]['UL_B']:.1e} | {c['n_hits']} | {fx} | {un} | {'' if cl is None else f'{cl:.1f}'} |"
                  if u1 and u3 else f"| {c['target']} | {'yes' if c['de_sky'] else 'no'} | — | — | {c['n_hits']} | {fx} | {un} | {'' if cl is None else f'{cl:.1f}'} |")
    md += ["", "Unresolved hits (single epoch, no multi-scan stack):", ""]
    for c in out["corridors"]:
        for g in c["hits"]:
            if g["disposition"] != "fixed_sky":
                md.append(f"- {c['target']}: {g['names']} ({', '.join(g['catalogues'])}) at {g['sep_arcmin']:.1f}′ — {g['reason']}")
    md += ["", "Closest fixed-sky sources (≤ 2′):", ""]
    for c in out["corridors"]:
        for g in c["hits"]:
            if g["disposition"] == "fixed_sky" and g["sep_arcmin"] <= 2.0:
                md.append(f"- {c['target']}: {g['names']} ({', '.join(g['catalogues'])}) at {g['sep_arcmin']:.1f}′ — {g['reason']}")
    (H.RES / "corridor_screen_v1.md").write_text("\n".join(md))
    print(json.dumps(out["census"], indent=1))


if __name__ == "__main__":
    main()
