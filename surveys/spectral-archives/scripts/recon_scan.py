"""Spectral-archive family reachability recon (plan §5.15 ledger item O4;
2026-09-05).

Channel A only (on-star: the star's uplink to the relay passes Earth
while Earth sits on the Sun-star axis, star at opposition). A spectrum
of the target star taken inside a crossing window is the substrate for
a laser-line inspection. Channel B (antipode) has no spectral substrate
(nobody takes spectra of empty sky at arbitrary positions; multi-object
sky fibres are shown to be geometrically hopeless in the note).

For each of the 7 deep-family targets (channel A b <= 0.1 AU:
gj-1276, gj-908, ross-128, ross-154, teegarden, van-maanen, wolf-359)
this script pulls every science spectrum in the public metadata of

  KOA (Keck: HIRES, NIRSPEC, KPF, NIRES, ESI, MOSFIRE, LRIS, DEIMOS, KCWI)
  ESO (dbo.raw all instruments + ivoa.ObsCore phase-3 spectra)
  CADC (CAOM: Gemini incl. MAROON-X/GNIRS/IGRINS/GRACES, HST, CFHT
        ESPaDOnS/SPIRou, DAO, JCMT, ...)
  MAST (CAOM spectra: HST, IUE, FUSE, GALEX grism ...)
  SDSS DR17 (apogeeVisit, specObjAll) + DR19 SDSS-V (mwm_apogee_allvisit,
        spAll_allepoch) via SkyServer
  NEID archive (neidl2)
  Breakthrough Listen Open Data (APF)
  CARMENES GTO DR1 (per-target SERVAL time series = per-spectrum list)
  SOPHIE archive (OHP; per-object HTML table)
  Data Lab DESI DR1 (zpix)
  LAMOST DR10 v1.0 (combined spectra, form POST)

and intersects the exposure mid-times with the universal channel-A
windows at the 1.2 Rsun, 2.5 Rsun and 0.1 AU rungs (flat-chord windows,
as in every crossings survey). Responses are snapshotted verbatim
under runs/spectral-archives/recon/.

No signal statistic is formed; nothing here is a search.
"""
from __future__ import annotations

import csv
import html as htmlmod
import io
import json
import re
import sys
import time
import urllib.parse
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

RUN = REPO / "runs" / "spectral-archives" / "recon"
OUT_DIR = Path(__file__).resolve().parents[1] / "results"
EVENTS = REPO / "crossings" / "universal_v1" / "events.ecsv"

KM_PER_AU = 1.495978707e8
RSUN_AU = 0.00465047
RUNGS = [("1.2Rsun", 1.2 * RSUN_AU), ("2.5Rsun", 2.5 * RSUN_AU), ("0.1AU", 0.1)]
CONE_DEG = 5.0 / 60.0          # 5 arcmin metadata cone
ONSTAR_ARCSEC = 180.0          # rows within 3' of the era-mean position count as on-star
DEEP = ["gj-1276", "gj-908", "ross-128", "ross-154", "teegarden", "van-maanen", "wolf-359"]

# archive-side names (SOPHIE object search; CARMENES Karmn ids)
SOPHIE_NAMES = {
    "wolf-359": ["Gl 406", "GJ406", "Wolf 359", "CN Leo"],
    "ross-128": ["Gl 447", "GJ447", "Ross 128"],
    "ross-154": ["Gl 729", "GJ729", "Ross 154"],
    "gj-908": ["Gl 908", "GJ908", "BR Psc"],
    "van-maanen": ["Gl 35", "GJ35", "van Maanen", "Wolf 28"],
    "teegarden": ["Teegarden", "GAT 1370", "SO0253"],
    "gj-1276": ["GJ 1276", "GJ1276"],
}
CARMENES_KARMN = {"teegarden": "J02530+168", "wolf-359": "J10564+070",
                  "ross-128": "J11477+008", "ross-154": "J18498-238"}

# wavelength coverage of the instruments met (nm), for the line-cell table
INSTRUMENT_BAND_NM = {
    "HIRES": (300, 1000), "KPF": (445, 870), "NIRSPEC": (950, 5500), "NIRES": (940, 2450),
    "ESI": (390, 1100), "MOSFIRE": (970, 2450), "LRIS": (320, 1000), "DEIMOS": (410, 1100),
    "HARPS": (378, 691), "ESPRESSO": (378, 789), "UVES": (300, 1100), "XSHOOTER": (300, 2480),
    "CRIRES": (950, 5300), "NIRPS": (980, 1800), "FEROS": (350, 920), "MAROON-X": (500, 920),
    "GNIRS": (800, 2500), "IGRINS-2": (1450, 2500), "IGRINS": (1450, 2500), "GRACES": (400, 1050),
    "ESPaDOnS": (370, 1050), "SPIRou": (980, 2350), "PHOENIX": (1000, 5000),
    "APOGEE": (1510, 1700), "BOSS": (360, 1040), "SDSS": (380, 920), "APF": (374, 970),
    "NEID": (380, 930), "CARMENES-VIS": (520, 960), "CARMENES-NIR": (960, 1710),
    "SOPHIE": (387, 694), "DESI": (360, 980), "STIS": (115, 1030), "COS": (90, 320),
}
LINES_NM = {"266": 266.0, "355": 355.0, "532": 532.0, "1064": 1064.0, "1550": 1550.0}


# ----------------------------------------------------------------- utils
def now_utc():
    return datetime.now(timezone.utc).isoformat()


def store_text(store, url, query, text, nrows, status=200):
    store.store(service_url=url, query=query, request_utc=now_utc(),
                response_bytes=text.encode("utf-8", "replace"), row_count=nrows,
                http_status=status)


def get(store, url, params=None, label="", timeout=300, tries=3, sleep=10):
    last = None
    for k in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            text = r.text
            store.store(service_url=url, query=urllib.parse.urlencode(params or {}),
                        request_utc=now_utc(), response_bytes=r.content,
                        row_count=None, http_status=r.status_code)
            if r.status_code >= 500:
                raise RuntimeError(f"HTTP {r.status_code}")
            return r.status_code, text
        except Exception as exc:
            last = exc
            print(f"  [{label}] attempt {k+1}: {exc}", file=sys.stderr, flush=True)
            time.sleep(sleep)
    return None, f"ERROR {last}"


def tap_csv(store, base, adql, label, timeout=300):
    status, text = get(store, base, {"REQUEST": "doQuery", "LANG": "ADQL",
                                     "FORMAT": "csv", "QUERY": adql}, label, timeout)
    if status is None or text.lstrip().startswith("<") or status >= 400:
        return None, text[:400]
    rows = list(csv.DictReader(io.StringIO(text)))
    return rows, None


def skyserver_csv(store, release, sql, label):
    base = f"https://skyserver.sdss.org/{release}/SkyServerWS/SearchTools/SqlSearch"
    status, text = get(store, base, {"cmd": sql, "format": "csv"}, label, 300)
    if status is None or text.lstrip().startswith("{"):
        return None, text[:400]
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    rows = list(csv.DictReader(io.StringIO("\n".join(lines))))
    return rows, None


def sep_arcsec(ra1, de1, ra2, de2):
    c1, c2 = np.radians(de1), np.radians(de2)
    x = np.sin(c1) * np.sin(c2) + np.cos(c1) * np.cos(c2) * np.cos(np.radians(ra1 - ra2))
    return float(np.degrees(np.arccos(np.clip(x, -1, 1))) * 3600)


def to_mjd(s):
    """Parse an ISO-ish UTC string to MJD; None on failure."""
    if s is None:
        return None
    s = str(s).strip().replace("Z", "")
    if not s:
        return None
    for fmt in ("isot", "iso"):
        try:
            return float(Time(s, format=fmt, scale="utc").mjd)
        except Exception:
            pass
    try:
        return float(Time(s.replace(" ", "T"), format="isot", scale="utc").mjd)
    except Exception:
        return None


def ffloat(x, default=None):
    try:
        v = float(x)
        return v if np.isfinite(v) else default
    except Exception:
        return default


# -------------------------------------------------------------- geometry
def load_events():
    t = Table.read(EVENTS)
    side = np.asarray(t["axis_distance_au"])
    A = t[(np.asarray(t["link_direction"]) == "inbound") & (side > 0)]
    A["mjd"] = Time(list(A["t_ca_utc"]), format="isot", scale="utc").mjd
    return A[np.asarray(A["b_min_au"]) <= 0.1]


def windows_for(ev):
    b = float(ev["b_min_au"]); vp = float(ev["v_perp_km_s"]); tc = float(ev["mjd"])
    out = {}
    for name, lim in RUNGS:
        if b >= lim:
            out[name] = None
        else:
            hw = np.sqrt(lim * lim - b * b) * KM_PER_AU / vp / 86400.0
            out[name] = (tc - hw, tc + hw, hw)
    return out


# ------------------------------------------------------------- archives
# every archive function returns a list of dict rows with at least:
#   archive, instrument, ident, target_name, ra, dec, t_start_mjd, exptime_s,
#   t_mid_mjd, band_nm (lo, hi) or None, level, public (bool|None), extra

def koa(store, ra, dec):
    base = "https://koa.ipac.caltech.edu/TAP/sync"
    tables = {
        "koa_hires": ("HIRES", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint,iodin,deckname,waveblue,wavered"),
        "koa_nirspec": ("NIRSPEC", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint,filter,waveblue,wavered"),
        "koa_kpf": ("KPF", "koaid,targname,ra,dec,date_beg,date_end,elaptime,progid,koaimtyp,propint"),
        "koa_nires": ("NIRES", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint"),
        "koa_esi": ("ESI", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint"),
        "koa_mosfire": ("MOSFIRE", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint"),
        "koa_lris": ("LRIS", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint"),
        "koa_deimos": ("DEIMOS", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint"),
        "koa_kcwi": ("KCWI", "koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint"),
    }
    rows, errs = [], {}
    for tb, (inst, cols) in tables.items():
        adql = (f"SELECT {cols} FROM {tb} WHERE CONTAINS(POINT('ICRS',ra,dec),"
                f"CIRCLE('ICRS',{ra:.6f},{dec:.6f},{CONE_DEG:.5f}))=1")
        got, err = tap_csv(store, base, adql, f"koa:{tb}")
        if got is None:
            # retry without the optional columns
            adql = (f"SELECT koaid,targname,ra,dec,date_obs,ut,elaptime,progid,koaimtyp,propint FROM {tb} "
                    f"WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.6f},{dec:.6f},{CONE_DEG:.5f}))=1")
            got, err = tap_csv(store, base, adql, f"koa:{tb}:min")
        if got is None:
            errs[tb] = err
            continue
        for r in got:
            if inst == "KPF":
                t0 = to_mjd(r.get("date_beg"))
            else:
                d = (r.get("date_obs") or "")[:10]
                t0 = to_mjd(f"{d}T{r.get('ut') or '00:00:00'}")
            exptime = ffloat(r.get("elaptime"), 0.0)
            band = INSTRUMENT_BAND_NM.get(inst)
            wb, wr = ffloat(r.get("waveblue")), ffloat(r.get("wavered"))
            if wb and wr and wr > wb:
                band = (wb / 10.0, wr / 10.0)   # KOA wave* are Angstrom
            rows.append(dict(archive="KOA", instrument=inst, ident=r.get("koaid"),
                             target_name=r.get("targname"), ra=ffloat(r.get("ra")),
                             dec=ffloat(r.get("dec")), t_start_mjd=t0, exptime_s=exptime,
                             band_nm=band, level="raw", public=(str(r.get("propint", "")).strip() not in ("", "0") or True),
                             imtype=r.get("koaimtyp"), progid=r.get("progid"),
                             extra={k: r.get(k) for k in ("iodin", "deckname", "filter") if k in r}))
    return rows, errs


def eso(store, ra, dec):
    base = "http://archive.eso.org/tap_obs/sync"
    rows, errs = [], {}
    adql = ("SELECT dp_id,instrument,object,ra,dec,exp_start,exposure,dp_tech,dp_type,prog_id,release_date "
            f"FROM dbo.raw WHERE dp_cat='SCIENCE' AND INTERSECTS(CIRCLE('ICRS',{ra:.6f},{dec:.6f},{CONE_DEG:.5f}), s_region)=1")
    got, err = tap_csv(store, base, adql, "eso:raw", timeout=600)
    if got is None:
        errs["dbo.raw"] = err
    else:
        for r in got:
            tech = (r.get("dp_tech") or "").upper()
            if not any(k in tech for k in ("SPECTRUM", "ECHELLE", "MOS", "IFU")):
                continue
            inst = (r.get("instrument") or "").split("/")[0]
            t0 = to_mjd(r.get("exp_start"))
            rel = to_mjd(r.get("release_date"))
            rows.append(dict(archive="ESO-raw", instrument=inst, ident=r.get("dp_id"),
                             target_name=r.get("object"), ra=ffloat(r.get("ra")), dec=ffloat(r.get("dec")),
                             t_start_mjd=t0, exptime_s=ffloat(r.get("exposure"), 0.0),
                             band_nm=INSTRUMENT_BAND_NM.get(inst), level="raw",
                             public=(rel is not None and rel <= Time.now().mjd),
                             imtype=tech, progid=r.get("prog_id"), extra={"dp_type": r.get("dp_type")}))
    adql = ("SELECT dp_id,instrument_name,target_name,s_ra,s_dec,t_min,t_max,t_exptime,em_min,em_max,em_res_power,"
            "obs_release_date,access_url,obs_collection FROM ivoa.ObsCore WHERE dataproduct_type='spectrum' AND "
            f"CONTAINS(POINT('ICRS',s_ra,s_dec),CIRCLE('ICRS',{ra:.6f},{dec:.6f},{CONE_DEG:.5f}))=1")
    got, err = tap_csv(store, base, adql, "eso:obscore", timeout=600)
    if got is None:
        errs["ivoa.ObsCore"] = err
    else:
        for r in got:
            inst = (r.get("instrument_name") or "").split("/")[0]
            emn, emx = ffloat(r.get("em_min")), ffloat(r.get("em_max"))
            band = (emn * 1e9, emx * 1e9) if emn and emx else INSTRUMENT_BAND_NM.get(inst)
            rel = to_mjd(r.get("obs_release_date"))
            rows.append(dict(archive="ESO-phase3", instrument=inst, ident=r.get("dp_id"),
                             target_name=r.get("target_name"), ra=ffloat(r.get("s_ra")), dec=ffloat(r.get("s_dec")),
                             t_start_mjd=ffloat(r.get("t_min")), exptime_s=ffloat(r.get("t_exptime"), 0.0),
                             band_nm=band, level="reduced-1D", public=(rel is not None and rel <= Time.now().mjd),
                             imtype="spectrum", progid=r.get("obs_collection"),
                             extra={"res_power": r.get("em_res_power"), "access_url": r.get("access_url")}))
    return rows, errs


def cadc(store, ra, dec):
    base = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync"
    adql = ("SELECT o.collection,o.instrument_name,o.target_name,o.observationID,p.productID,p.dataProductType,"
            "p.calibrationLevel,p.time_bounds_lower,p.time_bounds_upper,p.time_exposure,p.energy_bounds_lower,"
            "p.energy_bounds_upper,p.energy_resolvingPower,p.dataRelease,o.targetPosition_coordinates_cval1,o.targetPosition_coordinates_cval2 "
            "FROM caom2.Observation o JOIN caom2.Plane p ON o.obsID=p.obsID "
            f"WHERE INTERSECTS(CIRCLE('ICRS',{ra:.6f},{dec:.6f},{CONE_DEG:.5f}), p.position_bounds)=1 "
            "AND (p.dataProductType='spectrum' OR p.dataProductType LIKE '%spectrum')")
    got, err = tap_csv(store, base, adql, "cadc", timeout=600)
    if got is None:
        return [], {"caom2": err}
    rows = []
    for r in got:
        inst = (r.get("instrument_name") or "").split("/")[0]
        elo, ehi = ffloat(r.get("energy_bounds_lower")), ffloat(r.get("energy_bounds_upper"))
        band = (elo * 1e9, ehi * 1e9) if elo and ehi else INSTRUMENT_BAND_NM.get(inst)
        cra = ffloat(r.get("targetPosition_coordinates_cval1"), ra)
        cdec = ffloat(r.get("targetPosition_coordinates_cval2"), dec)
        rel = to_mjd(r.get("dataRelease"))
        t0 = ffloat(r.get("time_bounds_lower"))
        rows.append(dict(archive=f"CADC/{r.get('collection')}", instrument=inst, ident=f"{r.get('observationID')}/{r.get('productID')}",
                         target_name=r.get("target_name"), ra=cra, dec=cdec, t_start_mjd=t0,
                         exptime_s=ffloat(r.get("time_exposure"), 0.0), band_nm=band,
                         level=f"cal{r.get('calibrationLevel')}", public=(rel is not None and rel <= Time.now().mjd),
                         imtype=r.get("dataProductType"), progid=r.get("collection"),
                         extra={"res_power": r.get("energy_resolvingPower")}))
    return rows, {}


def mast(store, ra, dec):
    """MAST CAOM via the Mashup cone service (the TAP cone on dbo.obspointing
    takes > 3 min per position; the Mashup cone answers in ~1 min)."""
    api = "https://mast.stsci.edu/api/v0/invoke"
    req = {"service": "Mast.Caom.Cone", "params": {"ra": ra, "dec": dec, "radius": CONE_DEG},
           "format": "json", "pagesize": 50000, "page": 1}
    status, text = get(store, api, {"request": json.dumps(req)}, "mast:cone", 600)
    if status is None:
        return [], {"Mast.Caom.Cone": text[:200]}
    try:
        doc = json.loads(text)
    except Exception as exc:
        return [], {"Mast.Caom.Cone": f"json: {exc}"}
    rows = []
    for r in doc.get("data", []):
        if r.get("dataproduct_type") != "spectrum":
            continue
        inst = (r.get("instrument_name") or "").split("/")[0]
        emn, emx = ffloat(r.get("em_min")), ffloat(r.get("em_max"))
        band = (emn, emx) if emn and emx else INSTRUMENT_BAND_NM.get(inst)   # MAST em_* in nm
        rows.append(dict(archive=f"MAST/{r.get('obs_collection')}", instrument=inst, ident=r.get("obs_id"),
                         target_name=r.get("target_name"), ra=ffloat(r.get("s_ra")), dec=ffloat(r.get("s_dec")),
                         t_start_mjd=ffloat(r.get("t_min")), exptime_s=ffloat(r.get("t_exptime"), 0.0),
                         band_nm=band, level=f"cal{r.get('calib_level')}", public=(str(r.get("dataRights")).upper() == "PUBLIC"),
                         imtype="spectrum", progid=r.get("proposal_id"), extra={"filters": r.get("filters")}))
    return rows, {"status": doc.get("status")} if doc.get("status") != "COMPLETE" else {}


def lamost(store, ra, dec):
    """LAMOST DR10 v1.0 combined-spectra cone search (form POST; pipe-separated csv)."""
    url = "https://www.lamost.org/dr10/v1.0/q"
    fields = {"pos.type": "cone", "pos.racenter": f"{ra:.6f}", "pos.deccenter": f"{dec:.6f}",
              "pos.radius": f"{CONE_DEG*3600:.0f}", "output.fmt": "csv"}
    for c in ("obsid", "uid", "designation", "obsdate", "lmjd", "mjd", "planid", "spid", "fiberid",
              "ra_obs", "dec_obs", "snrr", "class", "subclass", "fibertype"):
        fields[f"output.combined.{c}"] = "on"
    last = None
    text = ""
    for k in range(3):
        try:
            r = requests.post(url, files={k2: (None, v) for k2, v in fields.items()}, timeout=300)
            store.store(service_url=url, query=urllib.parse.urlencode(fields), request_utc=now_utc(),
                        response_bytes=r.content, row_count=None, http_status=r.status_code)
            if r.status_code >= 400:
                raise RuntimeError(f"HTTP {r.status_code}")
            text = r.text
            break
        except Exception as exc:
            last = exc
            time.sleep(10)
    else:
        return [], {"dr10": f"{last}"}
    first = text.splitlines()[0] if text.strip() else ""
    if "|" not in first:
        return [], {"dr10": text[:200]}
    rows = []
    for r in csv.DictReader(io.StringIO(text), delimiter="|"):
        plan = r.get("combined_planid") or ""
        rows.append(dict(archive="LAMOST-DR10", instrument=("LAMOST-MRS" if re.search(r"M\d\d$", plan) else "LAMOST-LRS"),
                         ident=r.get("combined_obsid"), target_name=f"{r.get('combined_designation')} {r.get('combined_class')}/{r.get('combined_subclass')}",
                         ra=ffloat(r.get("combined_ra_obs")), dec=ffloat(r.get("combined_dec_obs")),
                         t_start_mjd=ffloat(r.get("combined_mjd")), exptime_s=0.0,
                         band_nm=(370, 900), level="combined", public=True, imtype=r.get("combined_fibertype"),
                         progid=plan, extra={"snrr": r.get("combined_snrr"), "time_resolution": "integer MJD (exposure UT in FITS header)"}))
    return rows, {}


def sdss(store, ra, dec):
    rows, errs = [], {}
    d = CONE_DEG
    dra = d / max(np.cos(np.radians(dec)), 0.1)
    box = f"ra BETWEEN {ra-dra:.5f} AND {ra+dra:.5f} AND dec BETWEEN {dec-d:.5f} AND {dec+d:.5f}"
    # DR17 APOGEE visits (SDSS-III/IV)
    got, err = skyserver_csv(store, "dr17", f"SELECT apogee_id,ra,dec,mjd,telescope,field,plate,fiberid,jd,dateobs FROM apogeeVisit WHERE {box}", "sdss17:apogeeVisit")
    if got is None:
        got, err = skyserver_csv(store, "dr17", f"SELECT apogee_id,ra,dec,mjd,telescope,field,plate,fiberid FROM apogeeVisit WHERE {box}", "sdss17:apogeeVisit:min")
    if got is None:
        errs["dr17.apogeeVisit"] = err
    else:
        for r in got:
            t0 = to_mjd(r.get("dateobs")) if r.get("dateobs") else ffloat(r.get("mjd"))
            rows.append(dict(archive="SDSS-DR17", instrument="APOGEE", ident=f"{r.get('plate')}-{r.get('mjd')}-{r.get('fiberid')}",
                             target_name=r.get("apogee_id"), ra=ffloat(r.get("ra")), dec=ffloat(r.get("dec")),
                             t_start_mjd=t0, exptime_s=ffloat(r.get("exptime"), 0.0), band_nm=INSTRUMENT_BAND_NM["APOGEE"],
                             level="visit", public=True, imtype="fiber", progid=r.get("field"),
                             extra={"telescope": r.get("telescope"), "time_resolution": "dateobs" if r.get("dateobs") else "integer MJD"}))
    # DR17 optical (SDSS/BOSS/eBOSS) spectra
    got, err = skyserver_csv(store, "dr17", f"SELECT specObjID,plate,mjd,fiberID,ra,dec,class,subClass,survey,instrument,programname FROM specObjAll WHERE {box}", "sdss17:specObjAll")
    if got is None:
        errs["dr17.specObjAll"] = err
    else:
        for r in got:
            rows.append(dict(archive="SDSS-DR17", instrument=(r.get("instrument") or "SDSS/BOSS"), ident=f"{r.get('plate')}-{r.get('mjd')}-{r.get('fiberID')}",
                             target_name=f"{r.get('class')}/{r.get('subClass')}", ra=ffloat(r.get("ra")), dec=ffloat(r.get("dec")),
                             t_start_mjd=ffloat(r.get("mjd")), exptime_s=0.0, band_nm=INSTRUMENT_BAND_NM["BOSS"],
                             level="coadd", public=True, imtype="fiber", progid=r.get("programname"),
                             extra={"survey": r.get("survey"), "time_resolution": "integer MJD (per-exposure TAI in spec file)"}))
    # DR19 SDSS-V APOGEE visits
    got, err = skyserver_csv(store, "dr19", f"SELECT sdss_id,ra,dec,mjd,telescope,field,fiber,exptime,date_obs FROM mwm_apogee_allvisit WHERE {box}", "sdss19:mwm_apogee_allvisit")
    if got is None:
        got, err = skyserver_csv(store, "dr19", f"SELECT sdss_id,ra,dec,mjd,telescope FROM mwm_apogee_allvisit WHERE {box}", "sdss19:mwm_apogee_allvisit:min")
    if got is None:
        errs["dr19.mwm_apogee_allvisit"] = err
    else:
        for r in got:
            t0 = to_mjd(r.get("date_obs")) if r.get("date_obs") else ffloat(r.get("mjd"))
            rows.append(dict(archive="SDSS-DR19", instrument="APOGEE", ident=f"{r.get('sdss_id')}@{r.get('mjd')}",
                             target_name=str(r.get("sdss_id")), ra=ffloat(r.get("ra")), dec=ffloat(r.get("dec")),
                             t_start_mjd=t0, exptime_s=ffloat(r.get("exptime"), 0.0), band_nm=INSTRUMENT_BAND_NM["APOGEE"],
                             level="visit", public=True, imtype="fiber", progid=r.get("field"),
                             extra={"telescope": r.get("telescope")}))
    # DR19 SDSS-V BOSS per-epoch
    boxc = box.replace("ra BETWEEN", "racat BETWEEN").replace("dec BETWEEN", "deccat BETWEEN")
    got, err = skyserver_csv(store, "dr19", f"SELECT sdss_id,racat,deccat,mjd,mjd_final,tai_list,nexp,exptime,firstcarton,programname,survey,obs FROM spAll_allepoch WHERE {boxc}", "sdss19:spAll_allepoch")
    if got is None:
        errs["dr19.spAll_allepoch"] = err
    else:
        for r in got:
            rows.append(dict(archive="SDSS-DR19", instrument="BOSS", ident=f"{r.get('sdss_id')}@{r.get('mjd')}",
                             target_name=str(r.get("sdss_id")), ra=ffloat(r.get("racat")), dec=ffloat(r.get("deccat")),
                             t_start_mjd=ffloat(r.get("mjd_final")) or ffloat(r.get("mjd")), exptime_s=ffloat(r.get("exptime"), 0.0),
                             band_nm=INSTRUMENT_BAND_NM["BOSS"], level="epoch", public=True, imtype="fiber",
                             progid=r.get("firstcarton"), extra={"tai_list": r.get("tai_list"), "nexp": r.get("nexp"), "obs": r.get("obs")}))
    return rows, errs


def neid(store, ra, dec):
    base = "https://neid.ipac.caltech.edu/TAP/sync"
    rows, errs = [], {}
    for tb in ("neidl2", "neidl1", "neidl0"):
        adql = (f"SELECT object,qrad,qdecd,obsmjd,obsdate,exptime,program,obsmode,obstype FROM {tb} "
                f"WHERE CONTAINS(POINT('ICRS',qrad,qdecd),CIRCLE('ICRS',{ra:.6f},{dec:.6f},{CONE_DEG:.5f}))=1")
        got, err = tap_csv(store, base, adql, f"neid:{tb}")
        if got is None:
            errs[tb] = err
            continue
        for r in got:
            rows.append(dict(archive="NEID", instrument="NEID", ident=f"{tb}:{r.get('obsdate')}",
                             target_name=r.get("object"), ra=ffloat(r.get("qrad")), dec=ffloat(r.get("qdecd")),
                             t_start_mjd=ffloat(r.get("obsmjd")), exptime_s=ffloat(r.get("exptime"), 0.0),
                             band_nm=INSTRUMENT_BAND_NM["NEID"], level=tb, public=True, imtype=r.get("obstype"),
                             progid=r.get("program"), extra={"obsmode": r.get("obsmode")}))
        if got:
            break   # l2 present -> no need for l1/l0 duplicates
    return rows, errs


def bl_apf(store, ra, dec):
    api = "http://seti.berkeley.edu/opendata/api/query-files"
    status, text = get(store, api, {"target": "", "pos-ra": f"{ra:.5f}", "pos-dec": f"{dec:.5f}",
                                    "pos-rad": "0.5", "limit": "100000"}, "bl", 300)
    if status is None:
        return [], {"query-files": text}
    try:
        doc = json.loads(text)
    except Exception as exc:
        return [], {"query-files": f"json: {exc}"}
    rows = []
    for r in doc.get("data", []):
        if r.get("telescope") != "APF":
            continue
        t0 = None
        try:
            t0 = float(Time(datetime.strptime(r["utc"], "%a, %d %b %Y %H:%M:%S GMT"), scale="utc").mjd)
        except Exception:
            pass
        rows.append(dict(archive="BL-OpenData", instrument="APF", ident=r.get("url", "").rsplit("/", 1)[-1],
                         target_name=r.get("target"), ra=ffloat(r.get("ra")), dec=ffloat(r.get("decl")),
                         t_start_mjd=t0, exptime_s=0.0, band_nm=INSTRUMENT_BAND_NM["APF"], level="raw-2D",
                         public=True, imtype="spectrum", progid="BL", extra={"url": r.get("url"), "size": r.get("size")}))
    return rows, {}


def carmenes(store, tid):
    karmn = CARMENES_KARMN.get(tid)
    if not karmn:
        return [], {"dr1": "target not in CARMENES GTO DR1 (362-star list checked 2026-09-05)"}
    url = "http://carmenes.cab.inta-csic.es/gto/getDR1DataPublic.action"
    status, text = get(store, url, {"id": f"{karmn}_SERVAL+RACOON.csv"}, f"carmenes:{karmn}", 120)
    if status is None or status >= 400 or not text.startswith("BJD"):
        return [], {"dr1": f"HTTP {status}: {text[:120]}"}
    rows = []
    for r in csv.DictReader(io.StringIO(text)):
        bjd = ffloat(r.get("BJD"))
        tim = r.get("TIMEID") or ""
        m = re.search(r"car-(\d{8})T(\d\d)h(\d\d)m(\d\d)s", tim)
        t0 = to_mjd(f"{m.group(1)[:4]}-{m.group(1)[4:6]}-{m.group(1)[6:]}T{m.group(2)}:{m.group(3)}:{m.group(4)}") if m else (bjd - 2400000.5 if bjd else None)
        rows.append(dict(archive="CARMENES-DR1", instrument="CARMENES-VIS", ident=tim, target_name=karmn,
                         ra=None, dec=None, t_start_mjd=t0, exptime_s=ffloat(r.get("EXPTIME"), 0.0),
                         band_nm=INSTRUMENT_BAND_NM["CARMENES-VIS"], level="reduced-1D(caracal)", public=True,
                         imtype="spectrum", progid="GTO", extra={"bjd_mid": bjd, "snref": r.get("SNREF"), "flag": r.get("FLAG")}))
    return rows, {}


def sophie(store, tid):
    base = "http://atlas.obs-hp.fr/sophie/sophie.cgi"
    rows, errs, seen = [], {}, set()
    for name in SOPHIE_NAMES.get(tid, []):
        status, text = get(store, base, {"n": "sophies", "c": "o", "of": "1,leda,simbad", "nra": "l", "nsr": "1", "o": name}, f"sophie:{name}", 120)
        if status is None:
            errs[name] = text[:120]
            continue
        txt = re.sub(r"<[^>]+>", " ", htmlmod.unescape(text))
        # table rows: objname  J-coord  S E seq slen date mode fiber_b exptime sn26 ...
        for m in re.finditer(r"\s(\S+)\s+J(\d{6}\.\d[+-]\d{6})\s+S\s+E\s+(\d+)\s+(\d+)\s+(\d{4}-\d\d-\d\d)\s+(H[ER])\s+(\S+)\s+([\d.]+)\s+([\d.]+)", txt):
            seq = m.group(3)
            if seq in seen:
                continue
            seen.add(seq)
            rows.append(dict(archive="SOPHIE", instrument="SOPHIE", ident=f"seq{seq}", target_name=m.group(1),
                             ra=None, dec=None, t_start_mjd=to_mjd(m.group(5)), exptime_s=ffloat(m.group(8), 0.0),
                             band_nm=INSTRUMENT_BAND_NM["SOPHIE"], level="reduced (e2ds/s1d)", public=None,
                             imtype="spectrum", progid=m.group(6), extra={"jcoord": m.group(2), "fiber_b": m.group(7), "sn26": m.group(9),
                                                                          "time_resolution": "date only in the listing (UT in header)"}))
    return rows, errs


def desi(store, ra, dec):
    base = "https://datalab.noirlab.edu/tap/sync"
    rows, errs = [], {}
    d = CONE_DEG
    dra = d / max(np.cos(np.radians(dec)), 0.1)
    adql = (f"SELECT targetid,mean_fiber_ra,mean_fiber_dec,survey,program,healpix,z,spectype,coadd_numexp,coadd_exptime,"
            f"min_mjd,max_mjd,mean_mjd FROM desi_dr1.zpix WHERE mean_fiber_ra BETWEEN {ra-dra:.6f} AND {ra+dra:.6f} "
            f"AND mean_fiber_dec BETWEEN {dec-d:.6f} AND {dec+d:.6f}")
    got, err = tap_csv(store, base, adql, "desi:zpix", timeout=300)
    if got is None:
        errs["zpix"] = err
        return rows, errs
    for r in got:
        rows.append(dict(archive="DESI-DR1", instrument="DESI", ident=str(r.get("targetid")),
                         target_name=f"{r.get('spectype')}/{r.get('survey')}/{r.get('program')}",
                         ra=ffloat(r.get("mean_fiber_ra")), dec=ffloat(r.get("mean_fiber_dec")),
                         t_start_mjd=ffloat(r.get("min_mjd")), exptime_s=ffloat(r.get("coadd_exptime"), 0.0),
                         band_nm=INSTRUMENT_BAND_NM["DESI"], level="coadd", public=True, imtype="fiber",
                         progid=r.get("program"), extra={"numexp": r.get("coadd_numexp"), "max_mjd": r.get("max_mjd"), "mean_mjd": r.get("mean_mjd"), "z": r.get("z")}))
    return rows, errs


# ------------------------------------------------------------------ main
def main():
    RUN.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    store = SnapshotStore(RUN)
    ev = load_events()
    tids = sorted(set(str(t) for t in ev["target_id"]))
    assert set(tids) == set(DEEP), tids

    report = {"generated_utc": now_utc(), "cone_deg": CONE_DEG, "onstar_arcsec": ONSTAR_ARCSEC, "targets": {}}
    all_rows = []
    for tid in DEEP:
        e = ev[np.asarray(ev["target_id"]) == tid]
        ra_m = float(np.mean(e["star_icrs_ra_deg"])); dec_m = float(np.mean(e["star_icrs_dec_deg"]))
        drift = max(sep_arcsec(ra_m, dec_m, float(r["star_icrs_ra_deg"]), float(r["star_icrs_dec_deg"])) for r in e)
        print(f"\n=== {tid}: era-mean position {ra_m:.5f} {dec_m:.5f}, max PM drift from mean {drift:.0f}\"", flush=True)
        events = []
        for r in e:
            w = windows_for(r)
            events.append(dict(event_id=str(r["event_id"]), t_ca_utc=str(r["t_ca_utc"]), t_ca_mjd=float(r["mjd"]),
                               b_rsun=float(r["b_min_au"]) / RSUN_AU, v_perp=float(r["v_perp_km_s"]),
                               windows={k: (None if v is None else {"lo": v[0], "hi": v[1], "half_days": v[2]}) for k, v in w.items()},
                               star_ra=float(r["star_icrs_ra_deg"]), star_dec=float(r["star_icrs_dec_deg"])))

        rows, errs = [], {}
        for name, fn in (("KOA", lambda: koa(store, ra_m, dec_m)), ("ESO", lambda: eso(store, ra_m, dec_m)),
                         ("CADC", lambda: cadc(store, ra_m, dec_m)), ("MAST", lambda: mast(store, ra_m, dec_m)),
                         ("SDSS", lambda: sdss(store, ra_m, dec_m)), ("NEID", lambda: neid(store, ra_m, dec_m)),
                         ("BL", lambda: bl_apf(store, ra_m, dec_m)), ("CARMENES", lambda: carmenes(store, tid)),
                         ("SOPHIE", lambda: sophie(store, tid)), ("DESI", lambda: desi(store, ra_m, dec_m)),
                         ("LAMOST", lambda: lamost(store, ra_m, dec_m))):
            try:
                rr, ee = fn()
            except Exception as exc:       # keep going; record
                rr, ee = [], {"exception": repr(exc)}
            for r in rr:
                r["target_id"] = tid
                if r.get("ra") is not None and r.get("dec") is not None:
                    r["sep_arcsec"] = sep_arcsec(ra_m, dec_m, r["ra"], r["dec"])
                else:
                    r["sep_arcsec"] = 0.0   # name-resolved archives (CARMENES, SOPHIE)
                r["t_mid_mjd"] = (r["t_start_mjd"] + 0.5 * (r.get("exptime_s") or 0.0) / 86400.0) if r.get("t_start_mjd") is not None else None
            print(f"  {name:9s} rows={len(rr):5d}  errors={ee if ee else '-'}", flush=True)
            rows += rr
            if ee:
                errs[name] = ee

        # ---- intersect ----
        onstar = [r for r in rows if r["sep_arcsec"] <= ONSTAR_ARCSEC and r.get("t_mid_mjd") is not None]
        for r in onstar:
            r["in_window"] = {}
            r["nearest_event_dt_days"] = None
            best = None
            for evd in events:
                dt = r["t_mid_mjd"] - evd["t_ca_mjd"]
                if best is None or abs(dt) < abs(best[0]):
                    best = (dt, evd["event_id"], evd["t_ca_utc"][:10])
                for rung, w in evd["windows"].items():
                    if w and w["lo"] <= r["t_mid_mjd"] <= w["hi"]:
                        r["in_window"][rung] = evd["event_id"]
            r["nearest_event_dt_days"] = best[0] if best else None
            r["nearest_event"] = best[1] if best else None
            r["nearest_event_date"] = best[2] if best else None

        # per instrument summary
        summ = defaultdict(lambda: {"n": 0, "n_onstar": 0, "t0": None, "t1": None, "in_window": defaultdict(int),
                                    "in_window_events": defaultdict(set), "band_nm": None, "nearest_dt_days": None})
        for r in rows:
            key = f"{r['archive']}:{r['instrument']}"
            s = summ[key]
            s["n"] += 1
            if r["sep_arcsec"] <= ONSTAR_ARCSEC:
                s["n_onstar"] += 1
                t = r.get("t_mid_mjd")
                if t is not None:
                    s["t0"] = t if s["t0"] is None else min(s["t0"], t)
                    s["t1"] = t if s["t1"] is None else max(s["t1"], t)
                    d = r.get("nearest_event_dt_days")
                    if d is not None and (s["nearest_dt_days"] is None or abs(d) < abs(s["nearest_dt_days"])):
                        s["nearest_dt_days"] = d
                for rung, eid in r.get("in_window", {}).items():
                    s["in_window"][rung] += 1
                    s["in_window_events"][rung].add(eid)
            if r.get("band_nm"):
                s["band_nm"] = list(r["band_nm"])
        for s in summ.values():
            s["in_window"] = dict(s["in_window"])
            s["in_window_events"] = {k: sorted(v) for k, v in s["in_window_events"].items()}
            for k in ("t0", "t1"):
                if s[k] is not None:
                    s[k] = Time(s[k], format="mjd").utc.isot[:10]
        report["targets"][tid] = {"era_mean_ra": ra_m, "era_mean_dec": dec_m, "max_pm_drift_arcsec": drift,
                                  "events": events, "archive_errors": errs, "per_instrument": dict(sorted(summ.items())),
                                  "in_window_spectra": [r for r in onstar if r["in_window"]]}
        all_rows += rows
        hits = [r for r in onstar if r["in_window"]]
        print(f"  --> {len(rows)} rows, {len(onstar)} on-star with times, {len(hits)} in-window")
        for r in sorted(hits, key=lambda r: r["t_mid_mjd"]):
            print(f"      {Time(r['t_mid_mjd'], format='mjd').utc.isot[:16]} {r['archive']}:{r['instrument']} "
                  f"{r['in_window']} dt={r['nearest_event_dt_days']:+.2f} d exp={r['exptime_s']} '{r['target_name']}'")

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple, set)):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        return o
    (OUT_DIR / "recon_scan_v0.json").write_text(json.dumps(clean(report), indent=1, default=str))
    (OUT_DIR / "recon_rows_v0.json").write_text(json.dumps(clean(all_rows), indent=0, default=str))
    print(f"\nwrote {OUT_DIR / 'recon_scan_v0.json'} and recon_rows_v0.json ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
