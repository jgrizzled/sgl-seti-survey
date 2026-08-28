"""Shared chain for the Rubin DP2 crossings survey.

Implements the frozen constructions (hypotheses.md v1.0 +
threshold_freeze_v1.json) once, used identically by the dev stage and
the blind confirmatory: z-grid apparent positions, dedup, the
footprint gate, DiaSource gates, S_det, and the 8-control rule.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

Z_GRID_AU = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
R_ASSOC_ARCSEC = 1.0
STATIC_EXCL_ARCSEC = 2.0
CONE_ARCSEC = 60.0            # one discovery cone covers unit + controls
N_CONTROLS = 8
CONTROL_STEP_DEG = 40.0
CONTROL_RESOLVE_DEG = 5.0
LOCUS_AVOID_ARCSEC = 2.5      # amendment v1.1 (was 10.0 in v1.0)
KM_PER_AU = 1.495978707e8
RSUN_KM = 695_700.0

#: freeze section 7 DiaSource exclusion flags (column names)
GATE_FLAGS = (
    "pixelFlags_bad", "pixelFlags_saturatedCenter", "pixelFlags_crCenter",
    "pixelFlags_edge", "pixelFlags_nodataCenter",
    "pixelFlags_interpolatedCenter", "pixelFlags_suspectCenter",
    "pixelFlags_streakCenter", "pixelFlags_injected",
    "pixelFlags_injectedCenter", "pixelFlags_injected_template",
    "pixelFlags_injected_templateCenter", "psfFlux_flag", "centroid_flag",
)
DIA_COLS = (("diaSourceId", "diaObjectId", "ssObjectId", "ra", "dec",
             "band", "visit", "detector", "midpointMjdTai", "psfFlux",
             "psfFluxErr", "reliability", "extendedness")
            + GATE_FLAGS)


def star_unit_vector(ev):
    ra = np.radians(float(ev["star_icrs_ra_deg"]))
    de = np.radians(float(ev["star_icrs_dec_deg"]))
    return np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra),
                     np.sin(de)])


def _to_radec(v):
    v = v / np.linalg.norm(v)
    return (float(np.degrees(np.arctan2(v[1], v[0])) % 360.0),
            float(np.degrees(np.arcsin(np.clip(v[2], -1, 1)))))


def relay_apparent(ev, z_au, t_mjd):
    """Apparent ICRS ra/dec of a relay at z_au on the anti-star axis."""
    t = Time(t_mjd, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", t).xyz.to_value("AU")
    earth = get_body_barycentric("earth", t).xyz.to_value("AU")
    return _to_radec((sun - z_au * star_unit_vector(ev)) - earth)


def axis_point(ev):
    """The z→∞ axis direction = the star antipode (pattern pivot)."""
    return _to_radec(-star_unit_vector(ev))


def sep_arcsec(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    c = (np.sin(d1) * np.sin(d2)
         + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))) * 3600.0


def offset_polar(pivot, pos):
    """(offset arcsec, position angle deg E of N) of pos from pivot."""
    dra = ((pos[0] - pivot[0] + 180.0) % 360.0 - 180.0) * np.cos(
        np.radians(pivot[1])) * 3600.0
    dde = (pos[1] - pivot[1]) * 3600.0
    return float(np.hypot(dra, dde)), float(
        np.degrees(np.arctan2(dra, dde)) % 360.0)


def polar_offset(pivot, r_arcsec, pa_deg):
    dra = r_arcsec * np.sin(np.radians(pa_deg)) / 3600.0
    dde = r_arcsec * np.cos(np.radians(pa_deg)) / 3600.0
    return (pivot[0] + dra / np.cos(np.radians(pivot[1])),
            pivot[1] + dde)


def z_positions(ev, t_mjd):
    """Deduplicated z-grid apparent positions (threshold freeze)."""
    pos = [relay_apparent(ev, z, t_mjd) for z in Z_GRID_AU]
    kept = []
    for p in pos:
        if all(sep_arcsec(*p, *q) >= R_ASSOC_ARCSEC for q in kept):
            kept.append(p)
    return kept


def point_in_quad(ra, dec, corners):
    cosd = np.cos(np.radians(dec))
    dra = [((c[0] - ra + 180.0) % 360.0 - 180.0) * cosd for c in corners]
    dde = [c[1] - dec for c in corners]
    inside = False
    n = len(corners)
    for i in range(n):
        x1, y1 = dra[i], dde[i]
        x2, y2 = dra[(i + 1) % n], dde[(i + 1) % n]
        if (y1 > 0) != (y2 > 0):
            if x1 + (0 - y1) * (x2 - x1) / (y2 - y1) > 0:
                inside = not inside
    return inside


def fetch_visit_detectors(client, store, visit_id):
    _, rows = client.query(
        "SELECT detector, llcra, llcdec, lrcra, lrcdec, urcra, urcdec, "
        "ulcra, ulcdec, magLim FROM dp2.VisitDetector WHERE "
        f"visitId = {visit_id}", store)
    return [{"detector": int(r[0]),
             "corners": [(float(r[1]), float(r[2])),
                         (float(r[3]), float(r[4])),
                         (float(r[5]), float(r[6])),
                         (float(r[7]), float(r[8]))],
             "magLim": float(r[9]) if r[9] not in ("", None) else None}
            for r in rows if all(x not in ("", None) for x in r[1:9])]


def on_detector(dets, pos):
    for d in dets:
        if d["magLim"] is not None and point_in_quad(pos[0], pos[1],
                                                     d["corners"]):
            return d
    return None


def fetch_diasources(client, store, visit_id, center, cone_arcsec):
    cols = ", ".join(DIA_COLS)
    _, rows = client.query(
        f"SELECT {cols} FROM dp2.DiaSource WHERE visit = {visit_id} "
        f"AND CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', "
        f"{center[0]:.7f}, {center[1]:.7f}, {cone_arcsec / 3600.0:.7f}))"
        " = 1", store)
    out = []
    for r in rows:
        rec = dict(zip(DIA_COLS, r))
        for k in ("ra", "dec", "psfFlux", "psfFluxErr"):
            rec[k] = float(rec[k]) if rec[k] not in ("", None) else np.nan
        rec["reliability"] = (float(rec["reliability"])
                              if rec["reliability"] not in ("", None)
                              else np.nan)
        for k in GATE_FLAGS:
            rec[k] = rec[k] in ("True", "1", "true", "T")
        out.append(rec)
    return out


def fetch_static_objects(client, store, center, cone_arcsec):
    mags = [f"{b}_psfMag" for b in "ugrizy"]
    _, rows = client.query(
        "SELECT objectId, coord_ra, coord_dec, " + ", ".join(mags) +
        " FROM dp2.Object WHERE CONTAINS(POINT('ICRS', coord_ra, "
        f"coord_dec), CIRCLE('ICRS', {center[0]:.7f}, {center[1]:.7f}, "
        f"{cone_arcsec / 3600.0:.7f})) = 1", store)
    out = []
    for r in rows:
        vals = [float(x) for x in r[3:9] if x not in ("", None)]
        out.append({"objectId": r[0], "ra": float(r[1]),
                    "dec": float(r[2]),
                    "min_psfMag": min(vals) if vals else np.nan})
    return out


def static_excluded(pos, objects, maglim):
    for o in objects:
        if (np.isfinite(o["min_psfMag"])
                and o["min_psfMag"] < maglim + 0.5
                and sep_arcsec(pos[0], pos[1], o["ra"], o["dec"])
                < STATIC_EXCL_ARCSEC):
            return True
    return False


def gate_pass(d):
    return (np.isfinite(d["psfFlux"]) and np.isfinite(d["psfFluxErr"])
            and d["psfFluxErr"] > 0
            and not any(d[k] for k in GATE_FLAGS))


def s_det(positions, diasources):
    """Frozen statistic + its best association (None if S_det = 0)."""
    best, best_d = 0.0, None
    for d in diasources:
        if not gate_pass(d):
            continue
        if any(sep_arcsec(p[0], p[1], d["ra"], d["dec"]) < R_ASSOC_ARCSEC
               for p in positions):
            snr = d["psfFlux"] / d["psfFluxErr"]
            if snr > best:
                best, best_d = snr, d
    return best, best_d


def control_patterns(ev, t_mjd, dets, objects, maglim):
    """The 8 frozen control patterns with the +5 deg resolution rule.

    Returns (patterns, census) where each pattern is a list of
    positions and census records the rotations applied.
    """
    pivot = axis_point(ev)
    base = [offset_polar(pivot, p) for p in z_positions(ev, t_mjd)]
    patterns, census = [], []
    for k in range(1, N_CONTROLS + 1):
        extra = 0.0
        while True:
            pat = [polar_offset(pivot, r,
                                pa + k * CONTROL_STEP_DEG + extra)
                   for r, pa in base]
            bad = any(on_detector(dets, p) is None
                      or static_excluded(p, objects, maglim)
                      for p in pat)
            near_locus = any(
                sep_arcsec(*p, *q) < LOCUS_AVOID_ARCSEC
                for p in pat for q in
                [polar_offset(pivot, r, pa) for r, pa in base])
            if not bad and not near_locus:
                break
            extra += CONTROL_RESOLVE_DEG
            if extra >= 360.0:
                pat = None
                break
        patterns.append(pat)
        census.append({"k": k, "extra_rotation_deg": extra,
                       "valid": pat is not None})
    return patterns, census


def run_chain(ev, visit_id, t_mjd, maglim_band, client, store,
              objects=None):
    """The full frozen chain for one (event, visit): returns a record
    with S_det, threshold T, exceedance, and the census."""
    dets = fetch_visit_detectors(client, store, visit_id)
    pos = z_positions(ev, t_mjd)
    hits = [on_detector(dets, p) for p in pos]
    usable = [(p, h) for p, h in zip(pos, hits) if h is not None]
    if not usable:
        return {"visit": visit_id, "status": "off_detector"}
    maglim = min(h["magLim"] for _, h in usable)
    pivot = axis_point(ev)
    if objects is None:
        objects = fetch_static_objects(client, store, pivot,
                                       CONE_ARCSEC + 30.0)
    live = [p for p, h in usable
            if not static_excluded(p, objects, maglim)]
    if not live:
        return {"visit": visit_id, "status": "not_constrainable",
                "reason": "all positions static-excluded"}
    dias = fetch_diasources(client, store, visit_id, pivot, CONE_ARCSEC)
    s, best = s_det(live, dias)
    patterns, census = control_patterns(ev, t_mjd, dets, objects, maglim)
    t_vals = []
    for pat in patterns:
        if pat is None:
            t_vals.append(np.nan)
            continue
        live_c = [p for p in pat
                  if not static_excluded(p, objects, maglim)]
        sc, _ = s_det(live_c, dias)
        t_vals.append(sc)
    t_arr = [v for v in t_vals if np.isfinite(v)]
    threshold = max(t_arr) if t_arr else np.nan
    return {"visit": visit_id, "status": "searched",
            "t_mjd": t_mjd, "maglim": maglim,
            "n_positions": len(pos), "n_on_detector": len(usable),
            "n_live": len(live), "n_diasources_cone": len(dias),
            "s_det": s,
            "association": ({k: (best[k] if not isinstance(best[k], float)
                                 or np.isfinite(best[k]) else None)
                             for k in ("diaSourceId", "diaObjectId",
                                       "ssObjectId", "ra", "dec",
                                       "psfFlux", "psfFluxErr",
                                       "reliability", "extendedness")}
                            if best else None),
            "control_s": t_vals, "threshold": threshold,
            "exceedance": bool(np.isfinite(threshold)
                               and s > max(threshold, 0.0)),
            "control_census": census}
