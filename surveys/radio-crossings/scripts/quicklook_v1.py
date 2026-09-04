"""Radio quick-look v1 (plan §5.15 item O1; hand-offs (a)+(b) of
report/radio_crossings_ext.md §4 and the VLASS six of
report/radio_crossings.md).

A *look*, not a search: the radio decision (§5.5, geometry-only) is
unchanged. For every 0.1 AU in-window epoch found by the coverage
intersections (ASKAP v2 in-footprint rows, the LoTSS v2 in-beam row,
the VLASS v1 rows) this script asks one question per epoch — is there
a catalogued or visible broadband continuum source at the relay
(channel B, antipode) or star (channel A) position in *that* epoch —
and records the answer with the local noise.

Stages
  1. Catalogue cone-search (anonymous TAP):
       * VAST full-survey epochs: the per-SBID component catalogue
         (CASDA ``AS207.vast_extragal_dr1_<field>_sb<sbid>_components_v01``).
       * RACS epochs: the matching RACS release catalogue
         (``AS110.racs_low2_v01`` for the 2022 low re-epoch incl. its
         Stokes V table, ``AS110.racs_high_components_v01``,
         ``AS110.racs_mid_components_v01`` with its ``sbid`` column).
       * LoTSS: ``lotss_dr3.main_sources`` (ASTRON VO) — a mosaic of
         every run on the pointing, so not epoch-resolved.
       * VLASS: CADC CAOM2 lists the quick-look tiles containing the
         position with exact time bounds (sharper than the daily
         tile-list dates the v1 arm used).
     The CASDA level-5 Selavy catalogues of the VAST *pilot* and FLASH
     SBIDs are not TAP tables; those epochs are recorded as
     ``catalogue: not_exposed`` (needs an OPAL login — stage 2 item).
  2. Image cutouts: VLASS quick-look tiles via CADC SODA (anonymous) —
     peak within one beam of the position, local rms from the annulus.
     ASKAP: CASDA DataLink needs an OPAL login (HTTP 401 anonymous;
     ``OPAL_USERNAME``/``OPAL_PASSWORD`` from the environment or the
     repo ``.env``). With it, the DataLink VOTable (basic auth) yields
     an authenticated token per product; the SODA async job at
     ``casda_data_access/data/async`` then takes ``ID=<token>`` with
     no further auth (the astroquery.casda flow). Access is per
     product: some (the VAST-pilot v1 products, a REJECTED FLASH
     image) return "no permission" and are recorded as such.
  1b. For epochs whose catalogue is not a TAP table (VAST pilot,
     FLASH, the 2024 RACS-mid re-observation) the per-SBID Selavy
     component catalogue is downloaded through the same DataLink and
     cone-searched locally.

Match radius: a component within MATCH_ARCSEC (= one ASKAP/VLASS
restoring beam or so) counts as coincident; the nearest component
within NEAR_ARCSEC is reported either way. Positions are the
per-event ICRS positions from crossings/universal_v1 (star: PM-
propagated to t_ca; relay: the axis antipode at t_ca; the relay's
apparent motion over a ±6 d window is < 1 arcmin).

Raw responses (TAP CSVs, CAOM rows, cutout FITS) are stored under
runs/radio-crossings/quicklook/.
"""

from __future__ import annotations

import json
import os
import re
import sys
import time as _time
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.snapshots import SnapshotStore  # noqa: E402

CASDA_TAP = "https://casda.csiro.au/casda_vo_tools/tap/sync"
ASTRON_TAP = "https://vo.astron.nl/__system__/tap/run/tap/sync"
CADC_TAP = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync"
CADC_FILES = "https://ws-cadc.canfar.net/minoc/files/"
CASDA_DATALINK = "https://data.csiro.au/casda_vo_proxy/vo/datalink/links?ID="
CASDA_SODA = "https://casda.csiro.au/casda_data_access/data/async"
ASKAP_ANN_ARCSEC = (30.0, 150.0)    # rms annulus in ASKAP cutouts


def opal_auth():
    """(user, password) from the environment or the repo .env; None if absent."""
    env = dict(os.environ)
    envfile = REPO / ".env"
    if envfile.exists():
        for line in envfile.read_text().splitlines():
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                env.setdefault(k.strip(), v.strip())
    u, p = env.get("OPAL_USERNAME"), env.get("OPAL_PASSWORD")
    return (u, p) if u and p else None

RES = REPO / "surveys" / "radio-crossings" / "results"
RUN = REPO / "runs" / "radio-crossings" / "quicklook"
MATCH_ARCSEC = 20.0
NEAR_ARCSEC = 180.0
RMS_PROXY_DEG = 0.5          # radius over which catalogue rms is medianed
CUT_RADIUS_DEG = 0.05        # VLASS cutout radius (3 arcmin)
VLASS_BEAM_ARCSEC = 2.5


# ----------------------------------------------------------------- helpers
def sep_arcsec(ra1, de1, ra2, de2):
    c1, c2 = np.radians(de1), np.radians(de2)
    dra = np.radians(np.asarray(ra1) - ra2)
    x = np.sin(c1) * np.sin(c2) + np.cos(c1) * np.cos(c2) * np.cos(dra)
    return np.degrees(np.arccos(np.clip(x, -1, 1))) * 3600.0


def tap_csv(session, store, url, query, label=""):
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
              "QUERY": query, "MAXREC": "200000"}
    for attempt in range(4):
        try:
            r = session.get(url, params=params, timeout=600)
            r.raise_for_status()
            text = r.text
            break
        except Exception as exc:
            print(f"  [{label}] attempt {attempt + 1}: {exc}", flush=True)
            _time.sleep(15 * (attempt + 1))
    else:
        return None
    if text.lstrip().startswith("<"):
        print(f"  [{label}] non-CSV response: {text[:200]}", flush=True)
        store.store(service_url=url, query=query, request_utc=Time.now().isot,
                    response_bytes=r.content, row_count=None,
                    http_status=r.status_code)
        return None
    tab = Table.read(text, format="ascii.csv") if text.count("\n") > 1 \
        else Table(names=text.strip().split(","))
    store.store(service_url=url, query=query, request_utc=Time.now().isot,
                response_bytes=r.content, row_count=len(tab),
                http_status=r.status_code)
    return tab


def fnum(v):
    try:
        if v is None or np.ma.is_masked(v):
            return None
        return float(v)
    except Exception:
        return None


# ---------------------------------------------------------------- epochs
def load_epochs():
    ev = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    by_id = {str(r["event_id"]): r for r in ev}

    def pos(ch, eid):
        r = by_id[eid]
        if ch == "A":
            return float(r["star_icrs_ra_deg"]), float(r["star_icrs_dec_deg"])
        return float(r["relay_icrs_ra_deg"]), float(r["relay_icrs_dec_deg"])

    epochs = []
    a = Table.read(RES / "askap_inwindow_v2.ecsv")
    a = a[(a["radius_au"] == 0.1) & a["in_footprint"]]
    for r in a:
        ra, de = pos(str(r["channel"]), str(r["event_id"]))
        epochs.append({
            "arm": "ASKAP", "channel": str(r["channel"]),
            "target_id": str(r["target_id"]), "event_id": str(r["event_id"]),
            "b_rsun": float(r["b_min_au"]) * 215.032,
            "t_ca_mjd": float(r["t_ca_mjd"]),
            "half_window_days": float(r["window_days"]) / 2,
            "collection": str(r["collection"]), "sbid": str(r["sbid"]),
            "field": str(r["field"]), "obs_start_mjd": float(r["obs_start_mjd"]),
            "obs_end_mjd": float(r["obs_end_mjd"]),
            "offset_days": float(r["offset_days"]),
            "freq_mhz": float(r["freq_mhz"]), "sep_field_deg": float(r["sep_deg"]),
            "quality": str(r["quality"]), "ra": ra, "dec": de,
        })
    lo = Table.read(RES / "lotss_inwindow_v2.ecsv")
    lo = lo[(lo["radius_au"] == 0.1) & lo["in_beam"]]
    for r in lo:
        ra, de = pos(str(r["channel"]), str(r["event_id"]))
        epochs.append({
            "arm": "LoTSS", "channel": str(r["channel"]),
            "target_id": str(r["target_id"]), "event_id": str(r["event_id"]),
            "b_rsun": float(r["b_min_au"]) * 215.032,
            "t_ca_mjd": float(r["t_ca_mjd"]),
            "half_window_days": float(r["window_days"]) / 2,
            "collection": "LoTSS DR3", "sbid": "", "field": str(r["pointing"]),
            "obs_start_mjd": float(r["obs_mid_mjd"]) - 4 / 24,
            "obs_end_mjd": float(r["obs_mid_mjd"]) + 4 / 24,
            "offset_days": float(r["offset_days"]), "freq_mhz": 144.0,
            "sep_field_deg": float(r["sep_deg"]), "quality": "",
            "ra": ra, "dec": de,
        })
    v = Table.read(RES / "vlass_inwindow_v1.ecsv")
    v = v[v["radius_au"] == 0.1]
    for r in v:
        ra, de = pos(str(r["channel"]), str(r["event_id"]))
        epochs.append({
            "arm": "VLASS", "channel": str(r["channel"]),
            "target_id": str(r["target_id"]), "event_id": str(r["event_id"]),
            "b_rsun": float(by_id[str(r["event_id"])]["b_min_solar_radii"]),
            "t_ca_mjd": float(r["t_ca_mjd"]),
            "half_window_days": float(r["window_days"]) / 2,
            "collection": str(r["vlass_epoch"]), "sbid": "", "field": "",
            "obs_start_mjd": float(r["tile_obs_mjd"]),
            "obs_end_mjd": float(r["tile_obs_mjd"]) + 1,
            "offset_days": float(r["offset_days"]), "freq_mhz": 3000.0,
            "sep_field_deg": float("nan"), "quality": "",
            "ra": ra, "dec": de,
        })
    return epochs


# ------------------------------------------------------------ stage 1
def casda_tables(session, store):
    t = tap_csv(session, store, CASDA_TAP,
                "SELECT table_name FROM TAP_SCHEMA.tables", "casda-tables")
    return set(str(x) for x in t["table_name"])


def cone_generic(session, store, url, table, racol, deccol, cols, ra, de,
                 radius_deg, extra="", label=""):
    q = (f"SELECT {cols} FROM {table} WHERE "
         f"1=CONTAINS(POINT('ICRS', {racol}, {deccol}), "
         f"CIRCLE('ICRS', {ra:.6f}, {de:.6f}, {radius_deg:.4f})){extra}")
    return tap_csv(session, store, url, q, label)


def summarise_cone(tab, ra, de, racol, deccol, peakcol, rmscol, epeakcol=None):
    """Nearest component, coincident flag, local-rms proxy."""
    out = {"n_cone": 0, "nearest_arcsec": None, "nearest_peak_mjy": None,
           "nearest_peak_err_mjy": None, "coincident": False,
           "rms_proxy_mjy": None}
    if tab is None:
        out["n_cone"] = None
        return out
    out["n_cone"] = len(tab)
    if len(tab) == 0:
        return out
    s = sep_arcsec(np.asarray(tab[racol], float), np.asarray(tab[deccol], float),
                   ra, de)
    i = int(np.argmin(s))
    out["nearest_arcsec"] = float(s[i])
    out["nearest_peak_mjy"] = fnum(tab[peakcol][i])
    if epeakcol and epeakcol in tab.colnames:
        out["nearest_peak_err_mjy"] = fnum(tab[epeakcol][i])
    out["coincident"] = bool(s[i] <= MATCH_ARCSEC)
    if rmscol in tab.colnames:
        rms = np.asarray(tab[rmscol], float)
        rms = rms[np.isfinite(rms) & (rms > 0)]
        if rms.size:
            out["rms_proxy_mjy"] = float(np.median(rms))
    return out


def stage1_askap(ep, session, store, tables):
    ra, de = ep["ra"], ep["dec"]
    coll = ep["collection"]
    sb = ep["sbid"].replace("ASKAP-", "")
    field = ep["field"].lower().replace("-", "m").replace("+", "p")
    res = {"catalogue": None, "table": None}
    if "VAST" in coll.upper() or ep["field"].startswith("VAST"):
        cand = [t for t in tables
                if t.startswith("AS207.") and f"_sb{sb}_" in t
                and t.endswith("components_v01")]
        if not cand:
            res["catalogue"] = "not_exposed"
            res["note"] = ("per-SBID Selavy catalogue not a TAP table "
                           "(VAST pilot); CASDA level-5 file needs OPAL login")
            return res
        table = cand[0]
        res["table"] = table
        cone = cone_generic(session, store, CASDA_TAP, table,
                            "ra_deg_cont", "dec_deg_cont",
                            "component_id, ra_deg_cont, dec_deg_cont, flux_peak, "
                            "flux_peak_err, flux_int, rms_image, maj_axis, min_axis",
                            ra, de, RMS_PROXY_DEG, label=table)
        res.update(summarise_cone(cone, ra, de, "ra_deg_cont", "dec_deg_cont",
                                  "flux_peak", "rms_image", "flux_peak_err"))
        res["catalogue"] = "epoch_resolved"
        return res
    if "Rapid ASKAP" in coll or ep["field"].startswith("RACS"):
        f = ep["freq_mhz"]
        if f < 1000:
            table = "AS110.racs_low2_v01" if ep["obs_start_mjd"] > 59500 \
                else "AS110.racs_dr1_sources_galacticcut_v2021_08_v02"
            vtable = "AS110.racs_low2_v_v01" if table.endswith("low2_v01") else None
        elif f < 1500:
            table, vtable = "AS110.racs_mid_components_v01", None
        else:
            table, vtable = "AS110.racs_high_components_v01", None
        res["table"] = table
        cols = "id, ra, dec, peak_flux, e_peak_flux, total_flux, noise, field_id"
        if table.endswith("mid_components_v01"):
            cols += ", sbid, scan_start_mjd"
        cone = cone_generic(session, store, CASDA_TAP, table, "ra", "dec", cols,
                            ra, de, RMS_PROXY_DEG, label=table)
        res.update(summarise_cone(cone, ra, de, "ra", "dec", "peak_flux",
                                  "noise", "e_peak_flux"))
        # which SBID/field fed the catalogue at this position?
        if cone is not None and len(cone):
            fields = sorted(set(str(x) for x in cone["field_id"]))
            res["catalogue_fields"] = fields
            if "sbid" in cone.colnames:
                sbids = sorted(set(str(x) for x in cone["sbid"]))
                res["catalogue_sbids"] = sbids
                res["catalogue"] = "epoch_resolved" if sb in sbids \
                    else "different_epoch"
            else:
                res["catalogue"] = "release_epoch"   # single-epoch survey
        else:
            res["catalogue"] = "empty_cone"
        if vtable:
            vc = cone_generic(session, store, CASDA_TAP, vtable, "ra", "dec",
                              "id, ra, dec, peak_flux, peak_flux_v, noise_v, "
                              "peak_vi, sbid, field_id",
                              ra, de, RMS_PROXY_DEG, label=vtable)
            vs = summarise_cone(vc, ra, de, "ra", "dec", "peak_flux_v", "noise_v")
            res["stokes_v"] = {"table": vtable, **vs}
            if vc is not None and len(vc) and "sbid" in vc.colnames:
                res["stokes_v"]["sbids"] = sorted(set(str(x) for x in vc["sbid"]))
        return res
    res["catalogue"] = "not_exposed"
    res["note"] = (f"{coll}: per-SBID Selavy catalogue not a TAP table; "
                   "CASDA level-5 file needs OPAL login")
    return res


def stage1_lotss(ep, session, store):
    ra, de = ep["ra"], ep["dec"]
    table = "lotss_dr3.main_sources"
    cone = cone_generic(session, store, ASTRON_TAP, table, "ra", "dec",
                        "source_name, ra, dec, peak_flux, e_peak_flux, "
                        "total_flux, isl_rms, mosaic_url",
                        ra, de, RMS_PROXY_DEG, label=table)
    res = {"catalogue": "mosaic_all_runs", "table": table}
    res.update(summarise_cone(cone, ra, de, "ra", "dec", "peak_flux",
                              "isl_rms", "e_peak_flux"))
    return res


# ------------------------------------------------------ CASDA (auth)
def casda_datalink(session, store, auth, product_id):
    """Parse the DataLink VOTable of one CASDA product: sync download URL,
    cutout-service token, or a permission error."""
    url = CASDA_DATALINK + product_id
    r = session.get(url, auth=auth, timeout=300)
    store.store(service_url=url, query=product_id, request_utc=Time.now().isot,
                response_bytes=r.content, row_count=None, http_status=r.status_code)
    if r.status_code != 200:
        return {"error": f"HTTP {r.status_code}"}
    x = r.text
    if "do not have permission" in x:
        return {"error": "no_permission"}
    out = {}
    m = re.search(r'(https://casda\.csiro\.au/casda_data_access/data/sync\?id=[A-Za-z0-9_-]+)', x)
    if m:
        out["sync_url"] = m.group(1)
    # rows: ID, access_url, service_def, ... ; the cutout row has
    # service_def = cutout_service and an authenticated token cell
    for row in re.findall(r"<TR>(.*?)</TR>", x, re.S):
        cells = re.findall(r"<TD>(.*?)</TD>", row, re.S)
        if "cutout_service" in cells:
            toks = [c for c in cells if re.fullmatch(r"[A-Za-z0-9_-]{40,}", c)]
            if toks:
                out["cutout_token"] = toks[-1]
    return out


def casda_soda_cutout(session, token, ra, de, radius_deg, dest, label=""):
    """Run one SODA async cutout job; returns the result URL or None."""
    r = session.post(CASDA_SODA, params={"ID": token}, timeout=300,
                     allow_redirects=False)
    job = r.headers.get("Location")
    if not job:
        print(f"  [{label}] SODA create: HTTP {r.status_code}", flush=True)
        return None
    session.post(job + "/parameters",
                 data={"CIRCLE": f"{ra:.6f} {de:.6f} {radius_deg:.4f}"},
                 timeout=300, allow_redirects=False)
    session.post(job + "/phase", data={"PHASE": "RUN"}, timeout=300,
                 allow_redirects=False)
    for _ in range(60):
        _time.sleep(10)
        ph = session.get(job + "/phase", timeout=300).text.strip()
        if ph in ("COMPLETED", "ERROR", "ABORTED"):
            break
    jx = session.get(job, timeout=300).text
    if ph != "COMPLETED":
        print(f"  [{label}] SODA job {ph}: {jx[:300]}", flush=True)
        return None
    hrefs = re.findall(r'xlink:href="([^"]+\.fits)"', jx)
    if not hrefs:
        return None
    rr = session.get(hrefs[0], timeout=600)
    if rr.status_code == 200 and rr.content[:6] == b"SIMPLE":
        dest.write_bytes(rr.content)
        return hrefs[0]
    print(f"  [{label}] download HTTP {rr.status_code}", flush=True)
    return None


def askap_cutout_photometry(path, ra, de):
    with fits.open(path) as h:
        data = np.squeeze(h[0].data).astype(float)
        hdr = h[0].header
    w = WCS(hdr).celestial
    bmaj = float(hdr["BMAJ"]) * 3600
    scale = abs(float(hdr.get("CDELT2", hdr.get("CD2_2")))) * 3600
    x, y = w.all_world2pix(ra, de, 0)
    ny, nx = data.shape
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - x, yy - y) * scale
    inbeam = data[(r <= bmaj) & np.isfinite(data)]
    ann = data[(r >= ASKAP_ANN_ARCSEC[0]) & (r <= ASKAP_ANN_ARCSEC[1]) & np.isfinite(data)]
    rms = 1.4826 * np.median(np.abs(ann - np.median(ann))) if ann.size else np.nan
    xi, yi = int(round(float(x))), int(round(float(y)))
    at = float(data[yi, xi]) if (0 <= xi < nx and 0 <= yi < ny) else np.nan
    peak = float(np.nanmax(inbeam)) if inbeam.size else np.nan
    return {"value_at_position_mjy": at * 1e3, "peak_in_beam_mjy": peak * 1e3,
            "rms_annulus_mjy": float(rms) * 1e3,
            "snr_value_at_position": float(at / rms) if rms > 0 else None,
            "snr_peak_in_beam": float(peak / rms) if rms > 0 else None,
            "bmaj_arcsec": bmaj, "pixel_arcsec": scale,
            "date_obs": hdr.get("DATE-OBS"), "n_annulus_px": int(ann.size),
            "inside_image": bool(inbeam.size > 0)}


def stage1b_selavy(ep, session, store, auth):
    """Download the per-SBID Selavy component catalogue(s) and cone-search."""
    ra, de = ep["ra"], ep["dec"]
    q = ("SELECT filename, access_url FROM ivoa.obscore WHERE "
         f"obs_id='{ep['sbid']}' AND dataproduct_subtype='catalogue.continuum.component'")
    t = tap_csv(session, store, CASDA_TAP, q, "selavy-" + ep["sbid"])
    out = {"catalogue": "selavy_download", "files": []}
    if t is None or len(t) == 0:
        out["catalogue"] = "no_selavy_product"
        return out
    best = None
    for row in t:
        pid = str(row["access_url"]).split("ID=")[-1]
        fn = str(row["filename"])
        dl = casda_datalink(session, store, auth, pid)
        rec = {"product": pid, "filename": fn, **{k: v for k, v in dl.items()
                                                  if k != "cutout_token"}}
        if "sync_url" in dl:
            dest = RUN / "selavy" / f"{ep['sbid']}_{pid}.xml"
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists():
                rr = session.get(dl["sync_url"], timeout=900)
                if rr.status_code == 200:
                    dest.write_bytes(rr.content)
                else:
                    rec["download"] = f"HTTP {rr.status_code}"
            if dest.exists():
                cat = Table.read(dest, format="votable")
                cat.rename_columns(cat.colnames, [c.replace("col_", "") for c in cat.colnames])
                sel = cat[sep_arcsec(np.asarray(cat["ra_deg_cont"], float),
                                     np.asarray(cat["dec_deg_cont"], float),
                                     ra, de) <= RMS_PROXY_DEG * 3600]
                summ = summarise_cone(sel, ra, de, "ra_deg_cont", "dec_deg_cont",
                                      "flux_peak", "rms_image", "flux_peak_err")
                rec.update({"n_catalogue": len(cat), "file": str(dest.relative_to(REPO)),
                            **summ})
                if best is None or (summ["n_cone"] or 0) > (best.get("n_cone") or 0):
                    best = rec
        out["files"].append(rec)
    if best:
        out.update({k: best[k] for k in ("n_cone", "nearest_arcsec", "nearest_peak_mjy",
                                         "nearest_peak_err_mjy", "coincident",
                                         "rms_proxy_mjy")})
        out["table"] = best["filename"]
    else:
        out["catalogue"] = "selavy_not_accessible"
    return out


def stage2_askap(ep, session, store, auth, obscore_rows):
    ra, de = ep["ra"], ep["dec"]
    rows = []
    for fn, url in obscore_rows:
        pid = url.split("ID=")[-1]
        rec = {"product": pid, "filename": fn}
        dl = casda_datalink(session, store, auth, pid)
        if "cutout_token" not in dl:
            rec["status"] = dl.get("error", "no_cutout_service")
            rows.append(rec)
            continue
        dest = RUN / "cutouts" / f"{ep['target_id']}_{ep['channel']}_{ep['event_id']}_{pid}.fits"
        dest.parent.mkdir(parents=True, exist_ok=True)
        if not dest.exists():
            href = casda_soda_cutout(session, dl["cutout_token"], ra, de,
                                     CUT_RADIUS_DEG, dest, label=pid)
            rec["soda_result"] = href
        if dest.exists():
            rec["status"] = "ok"
            rec["cutout"] = str(dest.relative_to(REPO))
            try:
                rec.update(askap_cutout_photometry(dest, ra, de))
            except Exception as exc:
                rec["photometry_error"] = str(exc)
        else:
            rec["status"] = "cutout_failed"
        rows.append(rec)
    return {"status": "ok" if any(r.get("status") == "ok" for r in rows)
            else "no_accessible_product", "cutouts": rows}


# ------------------------------------------------------------ VLASS
def vlass_tiles(ep, session, store):
    ra, de = ep["ra"], ep["dec"]
    q = ("SELECT o.observationID, o.proposal_id, p.productID, "
         "p.time_bounds_lower, p.time_bounds_upper, a.uri "
         "FROM caom2.Observation o JOIN caom2.Plane p ON o.obsID=p.obsID "
         "JOIN caom2.Artifact a ON p.planeID=a.planeID "
         f"WHERE o.collection='VLASS' AND CONTAINS(POINT('ICRS',{ra:.6f},{de:.6f}),"
         "p.position_bounds)=1 AND a.productType='science' "
         "AND p.productID LIKE '%quicklook%' AND a.uri LIKE '%tt0.subim.fits'")
    return tap_csv(session, store, CADC_TAP, q, "cadc-vlass")


def cutout_photometry(path, ra, de):
    with fits.open(path) as h:
        data = np.squeeze(h[0].data).astype(float)
        hdr = h[0].header
        w = WCS(hdr).celestial
        bmaj = float(hdr.get("BMAJ", VLASS_BEAM_ARCSEC / 3600)) * 3600
        date = hdr.get("DATE-OBS")
    x, y = w.all_world2pix(ra, de, 0)
    ny, nx = data.shape
    yy, xx = np.mgrid[:ny, :nx]
    scale = abs(float(hdr.get("CDELT2", hdr.get("CD2_2", 1 / 3600)))) * 3600
    r = np.hypot(xx - x, yy - y) * scale
    inbeam = data[(r <= bmaj) & np.isfinite(data)]
    ann = data[(r >= 20) & (r <= 120) & np.isfinite(data)]
    rms = 1.4826 * np.median(np.abs(ann - np.median(ann))) if ann.size else np.nan
    peak = float(np.nanmax(inbeam)) if inbeam.size else np.nan
    xi, yi = int(round(float(x))), int(round(float(y)))
    at = float(data[yi, xi]) if (0 <= xi < nx and 0 <= yi < ny) else np.nan
    return {"peak_in_beam_mjy": peak * 1e3, "value_at_position_mjy": at * 1e3,
            "snr_value_at_position": float(at / rms) if np.isfinite(rms) and rms > 0
            else None,
            "rms_annulus_mjy": float(rms) * 1e3,
            "snr_at_position": float(peak / rms) if np.isfinite(rms) and rms > 0
            else None, "bmaj_arcsec": bmaj, "date_obs": date,
            "inside_image": bool(0 <= x < nx and 0 <= y < ny and inbeam.size > 0)}


def stage2_vlass(ep, session, store, tiles):
    ra, de = ep["ra"], ep["dec"]
    lo = ep["t_ca_mjd"] - ep["half_window_days"]
    hi = ep["t_ca_mjd"] + ep["half_window_days"]
    rows = []
    if tiles is None:
        return {"tiles": None}
    for t in tiles:
        t0, t1 = fnum(t["time_bounds_lower"]), fnum(t["time_bounds_upper"])
        uri = str(t["uri"])
        in_win = (t0 is not None) and (t1 >= lo) and (t0 <= hi)
        row = {"observationID": str(t["observationID"]),
               "epoch": str(t["proposal_id"]), "uri": uri,
               "t_start_mjd": t0, "t_end_mjd": t1,
               "offset_days": (None if t0 is None
                               else round(0.5 * (t0 + t1) - ep["t_ca_mjd"], 3)),
               "in_window": bool(in_win),
               "in_window_v1_tolerance": bool(t0 is not None and
                                             abs(0.5 * (t0 + t1) - ep["t_ca_mjd"])
                                             <= ep["half_window_days"] + 3)}
        cut = RUN / "cutouts" / f"{ep['target_id']}_{ep['channel']}_{ep['event_id']}_{row['observationID']}.fits"
        cut.parent.mkdir(parents=True, exist_ok=True)
        if not cut.exists():
            url = f"{CADC_FILES}{uri}?CIRCLE={ra:.6f}+{de:.6f}+{CUT_RADIUS_DEG}"
            for attempt in range(3):
                try:
                    r = session.get(url, timeout=600)
                    if r.status_code == 200 and r.content[:6] == b"SIMPLE":
                        cut.write_bytes(r.content)
                        break
                    print(f"  cutout {row['observationID']}: HTTP {r.status_code}",
                          flush=True)
                except Exception as exc:
                    print(f"  cutout attempt {attempt + 1}: {exc}", flush=True)
                _time.sleep(10)
        if cut.exists():
            row["cutout"] = str(cut.relative_to(REPO))
            try:
                row.update(cutout_photometry(cut, ra, de))
            except Exception as exc:
                row["photometry_error"] = str(exc)
            if t0 is None and row.get("date_obs"):
                # CAOM time bounds absent (some VLASS4.1 planes): use the
                # tile header DATE-OBS (+ ~30 min) instead
                with fits.open(cut) as h:
                    d0 = h[0].header.get("DATE-OBS")
                    d1 = h[0].header.get("DATE-END")
                m0 = Time(d0).mjd
                m1 = Time(d1).mjd if d1 else m0 + 0.02
                row.update({"t_start_mjd": m0, "t_end_mjd": m1,
                            "time_source": "cutout_header",
                            "offset_days": round(0.5 * (m0 + m1) - ep["t_ca_mjd"], 3),
                            "in_window": bool(m1 >= lo and m0 <= hi),
                            "in_window_v1_tolerance": bool(
                                abs(0.5 * (m0 + m1) - ep["t_ca_mjd"])
                                <= ep["half_window_days"] + 3)})
            else:
                row["time_source"] = "caom2"
        rows.append(row)
    rows.sort(key=lambda r: (r["t_start_mjd"] or 0))
    return {"tiles": rows}


# ------------------------------------------------------------------ main
def main():
    RUN.mkdir(parents=True, exist_ok=True)
    store = SnapshotStore(RUN)
    session = requests.Session()
    session.headers["User-Agent"] = "sgl-seti-survey radio quicklook v1"
    epochs = load_epochs()
    print(f"{len(epochs)} epochs", flush=True)
    tables = casda_tables(session, store)
    auth = opal_auth()
    print("OPAL credentials:", "found" if auth else "absent", flush=True)
    out = []
    for ep in epochs:
        tag = f"{ep['arm']} {ep['channel']} {ep['target_id']} {ep['collection']} {ep['sbid']}"
        print(f"== {tag}  ({ep['ra']:.4f}, {ep['dec']:.4f})", flush=True)
        rec = dict(ep)
        if ep["arm"] == "ASKAP":
            rec["stage1"] = stage1_askap(ep, session, store, tables)
            # DataLink IDs of the restored image + noise map for stage 2
            q = ("SELECT dataproduct_subtype, filename, access_url, quality_level "
                 f"FROM ivoa.obscore WHERE obs_id='{ep['sbid']}' AND "
                 "filename LIKE 'image.i.%restored%' AND "
                 "dataproduct_subtype='cont.restored.t0'")
            t = tap_csv(session, store, CASDA_TAP, q, "obscore-" + ep["sbid"])
            obs_rows = [(str(x["filename"]), str(x["access_url"])) for x in t] \
                if t is not None else []
            if rec["stage1"].get("catalogue") in ("not_exposed", "different_epoch") \
                    and auth:
                rec["stage1_release"] = rec["stage1"]
                rec["stage1"] = stage1b_selavy(ep, session, store, auth)
            if auth:
                rec["stage2"] = stage2_askap(ep, session, store, auth, obs_rows)
            else:
                rec["stage2"] = {"status": "blocked_auth",
                                 "note": "CASDA DataLink returns HTTP 401 anonymously; "
                                         "needs an OPAL login",
                                 "datalink": [u for _, u in obs_rows],
                                 "filenames": [f for f, _ in obs_rows]}
        elif ep["arm"] == "LoTSS":
            rec["stage1"] = stage1_lotss(ep, session, store)
            rec["stage2"] = {"status": "not_attempted",
                             "note": "DR3 mosaic cutouts are not epoch-resolved"}
        else:
            tiles = vlass_tiles(ep, session, store)
            rec["stage1"] = {"catalogue": "caom2_tiles",
                             "n_tiles_all_epochs": None if tiles is None else len(tiles)}
            rec["stage2"] = stage2_vlass(ep, session, store, tiles)
        s1 = rec["stage1"]
        print(f"   stage1: {s1.get('catalogue')} n={s1.get('n_cone')} "
              f"nearest={s1.get('nearest_arcsec')} peak={s1.get('nearest_peak_mjy')} "
              f"rms~{s1.get('rms_proxy_mjy')} coincident={s1.get('coincident')}",
              flush=True)
        if ep["arm"] == "ASKAP" and rec["stage2"].get("cutouts"):
            for r in rec["stage2"]["cutouts"]:
                print(f"   {r['product']:12s} {r['status']:10s} at={r.get('value_at_position_mjy')} "
                      f"peak={r.get('peak_in_beam_mjy')} rms={r.get('rms_annulus_mjy')} "
                      f"{r.get('filename')}", flush=True)
        if ep["arm"] == "VLASS" and rec["stage2"].get("tiles"):
            for r in rec["stage2"]["tiles"]:
                print(f"   {r['epoch']:9s} {r['observationID']} off={r['offset_days']} "
                      f"win={r['in_window']} peak={r.get('peak_in_beam_mjy')} "
                      f"rms={r.get('rms_annulus_mjy')}", flush=True)
        out.append(rec)
    (RES / "quicklook_v1.json").write_text(json.dumps(out, indent=1, default=str))
    print("wrote", RES / "quicklook_v1.json")


if __name__ == "__main__":
    main()
