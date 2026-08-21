"""Build the universal (survey-agnostic) SGL target list.

Implements the simplified staged method from
notes/sgl_seti_star_ranking_methods.md:

  Stage 0  CNS5 census within the horizon, consolidated into systems.
  Stage 1  Network-side scores: distance rank; adjacency votes across
           seven sparse-graph families (Delaunay, Gabriel, RNG, MST,
           mutual-kNN k=1..3); linear-motion closest approach.
  Stage 2  Categorical host flags (multiplicity class, white dwarf,
           substellar, orbit availability). Labels only — nothing is
           removed from the distance core.
  Stage 2b Engineering columns (v2): mass/radius (TIC, FLAME, or
           Mann+2015/2019 M_Ks relations), lens compactness M/R^2,
           companion-acceleration class (PMa, HGCA chi2, RUWE),
           activity class (log R'HK, CARMENES H-alpha, ROSAT Lx/Lbol
           or its upper limit), evolutionary class. Each with a
           measured/limit/unknown provenance, per the ranking note
           §3.3. Feeds the *desirability* levels used by the
           selective re-vote — which is a network-taste filter,
           distinct from the solution gate (our solvability).
  Stage 3  Universal searchability columns: antipode coordinates,
           |b_gal|, |beta_ecl|, corridor drift rate, solution gate,
           already-searched status.
  Stage 4  Basket union -> labeled system list -> component tracks,
           trimmed to the track budget (never from the distance core).
           v2 baskets: DISTANCE_CORE, ROBUST_GEOMETRIC_NEIGHBOR,
           SELECTIVE_NETWORK_NEIGHBOR (re-vote on desirable systems,
           at each desirability level), ENGINEERING_BACKBONE (rank-sum
           of compactness/quiet/clean/lifetime among desirable
           dwarfs), SCIENCE_INTEREST (curated science_interest.yaml,
           tier 1 auto, tier 2 budget-permitting), HISTORICAL_NEIGHBOR,
           COMPACT_LENS_WILDCARD.

The network prior and archive searchability are kept as separate
columns and never merged into one number. Per-survey overlays (e.g.
surveys/wise/targets/) consume the output.

Usage: uv run python targets/scripts/build_universal_list.py
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import yaml
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import Delaunay
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.stats import rankdata

HERE = Path(__file__).resolve().parent
TDIR = HERE.parent
CFG = yaml.safe_load((TDIR / "config.yaml").read_text())

ALREADY_SEARCHED = {  # registry v1.1 systems (WISE shakedown complete)
    "alpha Cen": "wise",  # includes Proxima? Proxima is its own corridor
}
SEARCHED_SYSTEM_KEYS = ["alpha Cen", "Proxima", "Barnard", "Sirius",
                        "Ross 154", "Ross 248", "Lalande 21185",
                        "Wolf 359", "GJ 65", "CN Leo", "HH And",
                        # batch 2 (2026-08-19)
                        "Ross 128", "eps Ind", "tau Cet",
                        "YZ Cet",  # registry id gj-54 (v1 name truncation)
                        "Teegarden", "Lacaille 8760", "van Maanen",
                        "GJ 908", "GJ 784",
                        # batch 3 (2026-08-19)
                        "eps Eri", "Lacaille 9352", "GJ 1061",
                        "GJ 12724", "Wolf 1061",
                        # batch 4 + unblocked deferred (2026-08-20)
                        "61 Cyg", "Struve 2398", "Groombridge 34",
                        "GJ 1111", "Luyten", "Kapteyn", "LP 145-141",
                        "GJ 1221", "GJ 9193", "GJ 783", "GJ 11068",
                        "WISE 0855", "EZ Aqr", "Luhman 16", "Procyon",
                        # batch 5 (2026-08-21): universal v2 expansion
                        "GJ 876", "GJ 1002", "GJ 832", "GJ 526", "GJ 581",
                        "GJ 514", "Fomalhaut", "Wolf 437", "GJ 915",
                        "GJ 518", "GJ 1276", "61 Vir", "GJ 2012",
                        "GJ 11547", "LHS 1723", "82 Eri", "GJ 338",
                        "GJ 625", "GJ 293", "GJ 3306", "GJ 3112", "GJ 687",
                        "GJ 674", "GJ 682", "GJ 251", "sigma Dra",
                        "HD 219134", "GJ 1087", "GJ 318", "Wolf 1069",
                        "GJ 588", "GJ 3512", "GJ 13157", "GJ 2066",
                        "GJ 367", "GJ 229", "GJ 667", "LTT 1445", "GJ 66"]

KM_S_TO_PC_MYR = 1.0227
SCIENCE = yaml.safe_load((TDIR / "science_interest.yaml").read_text())


def fnum(v):
    try:
        return float(v) if v not in (None, "") else None
    except ValueError:
        return None


def load_engineering():
    path = TDIR / CFG["census"].get("engineering", "")
    if not path.is_file():
        return {}
    with open(path) as fh:
        return {r["CNS5"].strip(): r for r in csv.DictReader(fh)}


ENG = load_engineering()

GJ_NAMES = {
    "551": "Proxima Cen", "559": "alpha Cen AB", "699": "Barnard's Star",
    "406": "Wolf 359", "411": "Lalande 21185", "244": "Sirius AB",
    "65": "GJ 65 AB", "729": "Ross 154", "905": "Ross 248",
    "144": "eps Eri", "887": "Lacaille 9352", "447": "Ross 128",
    "866": "EZ Aqr", "280": "Procyon AB", "820": "61 Cyg AB",
    "725": "Struve 2398 AB", "15": "Groombridge 34 AB", "845": "eps Ind",
    "71": "tau Cet", "1061": "GJ 1061", "273": "Luyten's Star",
    "191": "Kapteyn's Star", "825": "Lacaille 8760", "860": "Kruger 60 AB",
    "234": "Ross 614 AB", "628": "Wolf 1061", "35": "van Maanen's Star",
    "1": "GJ 1", "473": "Wolf 424 AB", "687": "GJ 687", "876": "GJ 876",
    "166": "40 Eri ABC", "702": "70 Oph AB", "34": "eta Cas AB",
    "388": "AD Leo", "440": "LP 145-141 (WD)", "380": "Groombridge 1618",
    "832": "GJ 832", "682": "GJ 682", "674": "GJ 674",
    "54.1": "YZ Cet", "139": "82 Eri", "892": "HD 219134",
    "764": "sigma Dra", "663": "36 Oph AB", "664": "36 Oph C",
    "506": "61 Vir", "881": "Fomalhaut", "721": "Vega", "768": "Altair",
    "3193": "LTT 1445 ABC", "1253": "Wolf 1069", "486": "Wolf 437",
    "3323": "LHS 1723", "667": "GJ 667 ABC", "570": "GJ 570 ABCD",
    "229": "GJ 229 AB", "338": "GJ 338 AB", "412": "GJ 412 AB",
    "1245": "GJ 1245 ABC", "752": "GJ 752 AB", "1": "GJ 1",
    "191": "Kapteyn's Star", "205": "GJ 205", "588": "GJ 588",
    "433": "GJ 433", "436": "GJ 436", "357": "GJ 357", "367": "GJ 367",
    "3512": "GJ 3512", "251": "GJ 251", "581": "GJ 581", "1002": "GJ 1002",
    "809": "GJ 809", "445": "GJ 445", "526": "GJ 526", "1151": "GJ 1151",
    "625": "GJ 625", "686": "GJ 686", "849": "GJ 849", "176": "GJ 176",
    "514": "GJ 514", "3618": "GJ 3618", "83.1": "GJ 83.1 (TZ Ari)",
    "1116": "GJ 1116 AB", "1005": "GJ 1005 AB", "169.1": "Stein 2051 AB",
    "293": "GJ 293 (WD)", "518": "GJ 518 (WD)", "1087": "GJ 1087 (WD)",
}
CNS5_NAMES = {"2653": "Luhman 16 AB", "2194": "WISE 0855-0714",
              "723": "Teegarden's Star"}


def resolve_name(member):
    gj_key = member["gj"].strip() if member["gj"] else ""
    if gj_key.endswith(".0"):   # CNS5 writes e.g. "144.0" for GJ 144
        gj_key = gj_key[:-2]
    if member["cns5"] in CNS5_NAMES:
        return CNS5_NAMES[member["cns5"]]
    if gj_key and gj_key in GJ_NAMES:
        return GJ_NAMES[gj_key]
    if gj_key:
        return f"GJ {gj_key}"
    return f"CNS5 {member['cns5']}"


def comp_count(member):
    letters = [c for c in (member["comp"] or "") if c.isalpha()]
    return max(1, len(letters))



def load_census():
    rows = []
    with open(TDIR / CFG["census"]["snapshot"]) as fh:
        for r in csv.DictReader(fh):
            def f(k):
                v = r.get(k, "").strip()
                return float(v) if v else None
            plx = f("plx")
            rows.append({
                "cns5": r["CNS5"].strip(), "gj": r["GJ"].strip(),
                "comp": r["Comp"].strip(), "ncomp": f("NComp"),
                "gaia": r["GaiaDR3"].strip(), "hip": r["HIP"].strip(),
                "ra": f("RAJ2000"), "dec": f("DEJ2000"),
                "plx": plx, "e_plx": f("e_plx"),
                "pmra": f("pmRA") or 0.0, "pmde": f("pmDE") or 0.0,
                "rv": f("RV"), "g": f("Gmag"), "rp": f("RPmag"),
                "k": f("Ksmag"), "w1": f("W1mag"),
                "name": r["SimbadName"].strip(),
                "eng": ENG.get(r["CNS5"].strip(), {}),
            })
    rows = [r for r in rows if r["plx"] and r["ra"] is not None
            and 1000.0 / r["plx"] <= CFG["census"]["horizon_pc"]]
    for r in rows:
        d = 1000.0 / r["plx"]
        sc = SkyCoord(ra=r["ra"] * u.deg, dec=r["dec"] * u.deg,
                      distance=d * u.pc)
        c = sc.cartesian
        r["xyz"] = np.array([c.x.value, c.y.value, c.z.value])
        r["dist_pc"] = d
    return rows


def group_systems(rows):
    n = len(rows)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    thr = CFG["census"]["system_grouping_pc"]
    # Group on *projected* separation plus parallax consistency: 3D
    # separation splits genuine systems whose components carry
    # parallaxes from different catalogs (GJ 667 AB Hipparcos vs
    # C Gaia differ by 6%, i.e. 0.4 pc radially).
    for i in range(n):
        for j in range(i + 1, n):
            a, b = rows[i], rows[j]
            cosd = np.cos(np.deg2rad(a["dec"]))
            dth = np.deg2rad(np.hypot((a["ra"] - b["ra"]) * cosd,
                                      a["dec"] - b["dec"]))
            proj = dth * (a["dist_pc"] + b["dist_pc"]) / 2
            dplx = abs(a["plx"] - b["plx"]) / max(a["plx"], b["plx"])
            if proj < thr and dplx < CFG["census"]["system_plx_tol"]:
                parent[find(i)] = find(j)
    systems = defaultdict(list)
    for i in range(n):
        systems[find(i)].append(rows[i])
    out = []
    for members in systems.values():
        members.sort(key=lambda r: (r["g"] if r["g"] is not None else 99))
        prim = members[0]
        # Projected (sky-plane) separation: 3D separations are unusable
        # for close pairs because parallax errors inflate the radial
        # dimension (1 mas at 2.7 pc ~ 2000 AU).
        seps_au = []
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                cosd = np.cos(np.deg2rad(a["dec"]))
                dtheta = np.hypot((a["ra"] - b["ra"]) * cosd,
                                  a["dec"] - b["dec"]) * 3600.0
                seps_au.append(dtheta * (a["dist_pc"] + b["dist_pc"]) / 2)
        rvs = [m["rv"] for m in members if m["rv"] is not None]
        names = [resolve_name(m) for m in members]
        sysname = next((n for n in names if not n.startswith(("GJ ", "CNS5 "))),
                       names[0])
        out.append({
            "members": members, "primary": prim, "name": sysname,
            "xyz": prim["xyz"], "dist_pc": prim["dist_pc"],
            "n_comp": len(members),
            "min_sep_au": min(seps_au) if seps_au else None,
            "rv": float(np.mean(rvs)) if rvs else None,
            "pmra": prim["pmra"], "pmde": prim["pmde"],
            "ra": prim["ra"], "dec": prim["dec"], "plx": prim["plx"],
        })
    out.sort(key=lambda s: s["dist_pc"])
    return out


def simbad_short(name):
    return " ".join(name.split()) if name else None


# ---- Stage 1: graphs -------------------------------------------------

def adjacency_votes(systems):
    pts = np.vstack([[0.0, 0.0, 0.0]] + [s["xyz"] for s in systems])
    n = len(pts)
    dist = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
    votes = np.zeros(n, dtype=int)
    graphs = {}

    tri = Delaunay(pts)
    dedges = set()
    for simplex in tri.simplices:
        for a in simplex:
            for b in simplex:
                if a < b:
                    dedges.add((a, b))
    graphs["delaunay"] = {(a, b) for a, b in dedges}

    gabriel, rng = set(), set()
    for a in range(n):
        for b in range(a + 1, n):
            mid = (pts[a] + pts[b]) / 2
            r = dist[a, b] / 2
            others = [k for k in range(n) if k not in (a, b)]
            od = np.linalg.norm(pts[others] - mid, axis=1)
            if np.all(od >= r - 1e-12):
                gabriel.add((a, b))
            if not any(max(dist[a, k], dist[b, k]) < dist[a, b]
                       for k in others):
                rng.add((a, b))
    graphs["gabriel"] = gabriel
    graphs["rng"] = rng

    mst = minimum_spanning_tree(dist).tocoo()
    graphs["mst"] = {(min(i, j), max(i, j))
                     for i, j in zip(mst.row, mst.col)}

    order = np.argsort(dist, axis=1)
    for k in CFG["stage1"]["knn_ks"]:
        edges = set()
        for a in range(n):
            for b in order[a][1:k + 1]:
                if a in order[b][1:k + 1]:
                    edges.add((min(a, int(b)), max(a, int(b))))
        graphs[f"mknn{k}"] = edges

    sun_votes = {}
    for gi, (gname, edges) in enumerate(graphs.items()):
        for (a, b) in edges:
            if a == 0:
                sun_votes.setdefault(b, []).append(gname)
    return {i - 1: v for i, v in sun_votes.items()}, graphs


def closest_approach(sys):
    if sys["rv"] is None:
        return None
    sc = SkyCoord(ra=sys["ra"] * u.deg, dec=sys["dec"] * u.deg,
                  distance=(1000.0 / sys["plx"]) * u.pc,
                  pm_ra_cosdec=sys["pmra"] * u.mas / u.yr,
                  pm_dec=sys["pmde"] * u.mas / u.yr,
                  radial_velocity=sys["rv"] * u.km / u.s)
    v = sc.velocity.d_xyz.to_value(u.km / u.s) * KM_S_TO_PC_MYR
    r0 = sys["xyz"]
    W = CFG["stage1"]["historical_window_myr"]
    v2 = float(v @ v)
    if v2 == 0:
        return float(np.linalg.norm(r0))
    tstar = float(np.clip(-(r0 @ v) / v2, -W, W))
    return float(np.linalg.norm(r0 + v * tstar))


# ---- Stage 2/3 -------------------------------------------------------

def host_class(sys):
    merged_pair = any(comp_count(m) > 1 for m in sys["members"])
    if sys["n_comp"] > 1 or merged_pair:
        ms = sys["min_sep_au"]
        if merged_pair or (ms is not None
                           and ms < CFG["stage2"]["close_multiple_au"]):
            cls = "close"
        elif ms is not None and ms < CFG["stage2"]["wide_multiple_au"]:
            cls = "intermediate"
        else:
            cls = "wide"
    else:
        cls = "single"
    prim = sys["primary"]
    flags = []
    if prim["g"] is None and (prim["k"] is None or prim["k"] > 8.0) \
            and not prim["hip"]:  # HIP stars lacking G/Ks are saturated
        # No Gaia G and faint in Ks: substellar/ultracool. Bright stars
        # missing Gaia G (Sirius, Procyon, alpha Cen merged rows) are
        # saturated, not substellar.
        flags.append("substellar_or_faint")
    elif prim["g"] is not None:
        mg = prim["g"] + 5 * np.log10(prim["plx"] / 100.0)
        grp = (prim["g"] - prim["rp"]) if prim["rp"] is not None else None
        if mg > 10 and grp is not None and grp < 1.0:
            flags.append("white_dwarf")
    member_names = [resolve_name(m) for m in sys["members"]]
    orbit = any(any(k.lower() in n.lower() for n in
                    member_names + [sys["name"]])
                for k in CFG["stage2"]["orbit_available"])
    if cls == "close" and not orbit:
        flags.append("close_no_orbit")
    if "substellar_or_faint" in flags and not sys["primary"]["gaia"]:
        flags.append("substellar_no_gaia")
    return cls, flags, orbit


MANN19_A = [-0.642, -0.208, -8.43e-4, 7.87e-3, 1.42e-4, -2.13e-4]
# Last-resort main-sequence (mass, radius, L) by spectral type
# (Pecaut & Mamajek 2013 table), for Gaia-saturated stars lacking Ks/TIC.
SPT_TABLE = {"F0": (1.6, 1.7, 6.5), "F5": (1.33, 1.47, 3.3),
             "G0": (1.06, 1.1, 1.35), "G5": (0.98, 0.98, 0.83),
             "G8": (0.94, 0.91, 0.66), "K0": (0.9, 0.85, 0.46),
             "K2": (0.8, 0.78, 0.3), "K5": (0.7, 0.7, 0.16),
             "K7": (0.63, 0.63, 0.1), "M0": (0.57, 0.58, 0.07),
             "M2": (0.44, 0.43, 0.03), "M4": (0.2, 0.25, 0.007)}


def spt_lookup(spt):
    if not spt:
        return None
    import re
    m = re.match(r"([FGKM])(\d)", spt.strip())
    if not m:
        return None
    letter, num = m.group(1), int(m.group(2))
    keys = [k for k in SPT_TABLE if k[0] == letter]
    best = min(keys, key=lambda k: abs(int(k[1]) - num))
    return SPT_TABLE[best]


def mass_radius(prim, flags):
    """Mass, radius (solar) and provenance for the primary."""
    e = prim["eng"]
    if "white_dwarf" in flags:
        return 0.6, 0.0125, 3e-4, "wd_nominal"
    if "substellar_or_faint" in flags:
        return None, None, None, "unknown_substellar"
    mk = None
    if prim["k"] is not None and prim["plx"]:
        mk = prim["k"] + 5 * np.log10(prim["plx"] / 100.0)
    tic_m, tic_r = fnum(e.get("tic_mass")), fnum(e.get("tic_rad"))
    fl_m, fl_r = fnum(e.get("flame_mass")), fnum(e.get("flame_rad"))
    lum = fnum(e.get("tic_lum")) or None   # TIC writes 0.0 for missing

    def ml(m):  # crude MS mass-luminosity fallback
        return 0.23 * m ** 2.3 if m < 0.43 else m ** 4
    if mk is not None and 4.5 < mk < 10.5:
        # Mann+2019 mass, Mann+2015 radius (radius extrapolated mildly
        # beyond 9.8); preferred over TIC for uniformity on M dwarfs.
        z = mk - 7.5
        m = 10 ** sum(a * z ** i for i, a in enumerate(MANN19_A))
        r = 1.9515 - 0.3520 * mk + 0.01680 * mk ** 2
        if lum is None:  # BC_K ~ 2.6 for M dwarfs (Mann+2015)
            lum = 10 ** (-0.4 * (mk + 2.6 - 4.74))
        return m, r, lum, "mann_mks"
    if tic_m and tic_r:
        return tic_m, tic_r, lum or ml(tic_m), "tic"
    if fl_m and fl_r:
        return fl_m, fl_r, lum or ml(fl_m), "flame"
    if tic_r and fl_m:
        return fl_m, tic_r, lum or ml(fl_m), "flame+tic"
    t = spt_lookup(e.get("sptype"))
    if t and "V" in (e.get("sptype") or ""):
        return t[0], t[1], lum or t[2], "sptype_table"
    return None, None, lum, "unknown"


def activity_class(prim, mass, lum):
    """(class, evidence). Precedence: any measured active > any
    measured quiet > X-ray upper limit > unknown."""
    e = prim["eng"]
    c = CFG["stage2b"]["activity"]
    active, quiet = [], []
    lxlb = None
    rhk = fnum(e.get("logRpHK"))
    if rhk is not None:
        (active if rhk > c["logRpHK_active"] else quiet).append(
            f"logR'HK={rhk:.2f}")
    ha = fnum(e.get("carm_pEWHa"))
    if ha is not None:
        (active if ha < c["pEWHa_active_A"] else quiet).append(
            f"pEWHa={ha:+.2f}A")
    cr = fnum(e.get("rosat_crate"))
    lbol = None
    if lum:
        lbol = lum * 3.828e33
    elif mass:
        lbol = 3.828e33 * mass ** 4  # crude MS scaling, flagged below
    d_cm = 1000.0 / prim["plx"] * 3.086e18
    if cr is not None and lbol:
        lx = 4 * np.pi * d_cm ** 2 * cr * c["rosat_ecf_erg_cm2_ct"]
        ratio = np.log10(lx / lbol)
        lxlb = float(ratio)
        (active if ratio > c["log_lx_lbol_active"] else quiet).append(
            f"logLx/Lbol={ratio:.1f}")
    elif cr is None and lbol:
        lx_lim = 4 * np.pi * d_cm ** 2 * c["rosat_limit_ct_s"] * \
            c["rosat_ecf_erg_cm2_ct"]
        lim = np.log10(lx_lim / lbol)
        lxlb = float(lim)
        if lim < c["log_lx_lbol_active"]:
            quiet.append(f"logLx/Lbol<{lim:.1f}(RASS limit)")
    if active:
        return "active", active + quiet, lxlb
    if quiet:
        return "quiet", quiet, lxlb
    return "unknown", [], lxlb


def accel_class(prim):
    e = prim["eng"]
    c = CFG["stage2b"]["accel"]
    pma, chi2, ruwe = (fnum(e.get("pma_snr")), fnum(e.get("hgca_chi2")),
                       fnum(e.get("ruwe")))
    ev = []
    det = False
    if pma is not None:
        ev.append(f"PMa S/N={pma:.1f}")
        det |= pma > c["pma_snr"]
    if chi2 is not None:
        ev.append(f"HGCA chi2={chi2:.1f}" + (" (PMa decides)" if pma is not None else ""))
        det |= chi2 > c["hgca_chi2"] and pma is None
    bright = prim["g"] is not None and prim["g"] < c["ruwe_bright_g"]
    if ruwe is not None:
        ev.append(f"RUWE={ruwe:.2f}" + (" (bright, ignored)" if bright else ""))
        det |= ruwe > c["ruwe"] and not bright
    if det:
        return "accel_detected", ev
    if pma is not None or chi2 is not None:
        return "clean", ev
    if ruwe is not None:
        return "clean_ruwe_only", ev
    return "unknown", ev


def evol_class(sys, prim, mass):
    e = prim["eng"]
    c = CFG["stage2b"]
    names = [resolve_name(m) for m in sys["members"]] + [sys["name"]]
    if any(k.lower() in n.lower() for n in names for k in c["known_evolved"]):
        return "evolved", ["curated"]
    ev = fnum(e.get("flame_evol"))
    if ev is not None and ev >= c["flame_evol_subgiant"]:
        return "evolved", [f"FLAME evolstage={ev:.0f}"]
    lg, te = fnum(e.get("gspphot_logg")), fnum(e.get("gspphot_teff"))
    if lg is not None and te is not None and te < 7000 and lg < 4.0:
        return "evolved", [f"logg={lg:.2f}"]
    if mass is not None and mass > c["massive_short_lived_msun"]:
        return "massive_short_lived", [f"M={mass:.2f}"]
    return "main_sequence", []


def engineering(sys, cls, flags):
    prim = sys["primary"]
    mass, rad, lum, src = mass_radius(prim, flags)
    act, act_ev, lxlb = activity_class(prim, mass, lum)
    acc, acc_ev = accel_class(prim)
    evo, evo_ev = evol_class(sys, prim, mass)
    names = [resolve_name(m) for m in sys["members"]] + [sys["name"]]
    if any(k.lower() in n.lower() for n in names
           for k in CFG["stage2b"]["known_rapid_rotators"]):
        acc_ev = acc_ev + ["rapid rotator (curated)"]
    comp = (mass / rad ** 2) if (mass and rad) else None
    age = fnum(prim["eng"].get("flame_age"))
    # Desirability levels (network taste; NOT the solution gate).
    lv = CFG["stage2b"]["desirability"]
    undesirable1 = (cls == "close" or "substellar_or_faint" in flags
                    or evo != "main_sequence" and "white_dwarf" not in flags)
    undesirable2 = undesirable1 or cls == "intermediate" or \
        acc == "accel_detected" or act == "active" or \
        any(k.lower() in n.lower() for n in names
            for k in CFG["stage2b"]["known_rapid_rotators"])
    if undesirable1:
        level = 0
    elif undesirable2:
        level = 1
    else:
        level = 2
    return {"mass_msun": round(mass, 3) if mass else None,
            "radius_rsun": round(rad, 4) if rad else None,
            "lum_lsun": round(lum, 5) if lum else None,
            "mass_radius_source": src,
            "compactness_m_r2": round(comp, 2) if comp else None,
            "activity": act, "activity_evidence": act_ev,
            "log_lx_lbol": round(lxlb, 2) if lxlb is not None else None,
            "acceleration": acc, "acceleration_evidence": acc_ev,
            "evolution": evo, "evolution_evidence": evo_ev,
            "flame_age_gyr": age,
            "sptype": prim["eng"].get("sptype") or None,
            "desirability_level": level}


def engineering_score(systems):
    """Rank-sum over desirable (level 2) main-sequence dwarfs: higher
    is better. Components: compactness, quietness, cleanliness,
    lifetime, distance (attachment convenience)."""
    pool = [s for s in systems
            if s["eng"]["desirability_level"] == 2
            and "white_dwarf" not in s["flags"]
            and s["eng"]["compactness_m_r2"] is not None
            and (s["eng"]["activity"] == "quiet"
                 or not CFG["stage4"]["engineering_require_measured_quiet"])]
    if not pool:
        return {}
    n = len(pool)
    comp = np.array([s["eng"]["compactness_m_r2"] for s in pool])
    dist = np.array([s["dist_pc"] for s in pool])
    # Continuous quietness: log Lx/Lbol (measured or RASS limit); stars
    # classed quiet by R'HK/H-alpha alone sit at a neutral -4.5.
    act = np.array([s["eng"]["log_lx_lbol"] if s["eng"]["log_lx_lbol"]
                    is not None else -4.5 for s in pool])
    act += np.array([{"quiet": 0.0, "unknown": 1.0}.get(
        s["eng"]["activity"], 2.0) for s in pool])
    acc = np.array([{"clean": 0.0, "clean_ruwe_only": 0.25,
                     "unknown": 0.5}.get(s["eng"]["acceleration"], 1.0)
                    for s in pool])
    mass = np.array([s["eng"]["mass_msun"] for s in pool])
    life = np.clip((mass - 1.0) / 0.3, 0, 1)  # >1.0 Msun: shorter MS life
    power = np.array([np.log10(s["eng"]["lum_lsun"] or 1e-4) for s in pool])

    def rank_hi(x):  # 1.0 = best (largest x), ~0 = worst; ties share rank
        return rankdata(x, method="average") / n

    w = CFG["stage4"]["engineering_weights"]
    score = (w["compactness"] * rank_hi(comp)
             + w["activity"] * rank_hi(-act)
             + w["acceleration"] * rank_hi(-acc)
             + w["lifetime"] * rank_hi(-life)
             + w["power"] * rank_hi(power)
             + w["distance"] * rank_hi(-dist))
    return {s["name"]: round(float(sc), 3) for s, sc in zip(pool, score)}


def searchability(sys, cls, flags, orbit):
    anti = SkyCoord(ra=(sys["ra"] + 180.0) % 360.0 * u.deg,
                    dec=-sys["dec"] * u.deg)
    b = float(abs(anti.galactic.b.deg))
    beta = float(abs(anti.barycentricmeanecliptic.lat.deg))
    drift = float(np.hypot(sys["pmra"], sys["pmde"]) / 1000.0)  # "/yr
    if orbit:
        gate = "clean"
    elif "substellar_no_gaia" in flags:
        gate = "hard"
    elif cls == "close":
        gate = "orbit_needed"
    else:
        gate = "clean"
    member_names = [resolve_name(m) for m in sys["members"]]
    searched = any(any(k.lower() in n.lower() for n in
                       member_names + [sys["name"]])
                   for k in SEARCHED_SYSTEM_KEYS)
    return {"anti_ra": round(float(anti.ra.deg), 3),
            "anti_dec": round(float(anti.dec.deg), 3),
            "abs_gal_b": round(b, 1), "abs_ecl_beta": round(beta, 1),
            "drift_arcsec_yr": round(drift, 2), "solution_gate": gate,
            "already_searched_wise": searched}


# ---- main ------------------------------------------------------------

def main():
    rows = load_census()
    systems = group_systems(rows)
    print(f"{len(rows)} census objects -> {len(systems)} systems "
          f"within {CFG['census']['horizon_pc']} pc")

    votes, _ = adjacency_votes(systems)
    for i, s in enumerate(systems):
        s["votes"] = votes.get(i, [])
        s["d_min_pc"] = closest_approach(s)
        s["cls"], s["flags"], s["orbit"] = host_class(s)
        s["eng"] = engineering(s, s["cls"], s["flags"])
        s["search"] = searchability(s, s["cls"], s["flags"], s["orbit"])
        s["gj"] = s["primary"]["gj"].strip().removesuffix(".0")

    # Basket C: selective re-vote on the *desirable* subgraph at each
    # desirability level (the network's taste, not our solvability).
    c4 = CFG["stage4"]
    sel_maps = {}
    for level in CFG["stage2b"]["desirability"]["revote_levels"]:
        keep = [s for s in systems if s["eng"]["desirability_level"] >= level]
        sel_votes, _ = adjacency_votes(keep)
        keep_names = [s["name"] for s in keep]
        sel_maps[level] = {keep_names[i]: v for i, v in sel_votes.items()}
        print(f"  level>={level}: {len(keep)} desirable systems, "
              f"{sum(1 for v in sel_votes.values() if len(v) >= c4['adjacency_vote_min'])}"
              f" Sun-neighbors with >= {c4['adjacency_vote_min']} votes")

    eng_scores = engineering_score(systems)
    eng_rank = sorted(eng_scores, key=lambda k: -eng_scores[k])
    # Basket size counts additions beyond the distance core; core
    # systems scoring above the cutoff carry the label too.
    core_names = {s["name"] for s in systems[:c4["n_distance_core"]]}
    noncore = [k for k in eng_rank if k not in core_names]
    eng_cut = (eng_scores[noncore[c4["n_engineering"] - 1]]
               if len(noncore) >= c4["n_engineering"] else -1)
    sci1 = {str(k): v for k, v in SCIENCE["tier1"].items()}
    sci2 = {str(k): v for k, v in SCIENCE["tier2"].items()}

    baskets = defaultdict(list)
    for rank, s in enumerate(systems):
        if rank < c4["n_distance_core"]:
            baskets[s["name"]].append("DISTANCE_CORE")
        if len(s["votes"]) >= c4["adjacency_vote_min"]:
            baskets[s["name"]].append("ROBUST_GEOMETRIC_NEIGHBOR")
        s["selective_votes"] = {lv: len(m.get(s["name"], []))
                                for lv, m in sel_maps.items()}
        if any(n >= c4["adjacency_vote_min"]
               for n in s["selective_votes"].values()) and \
                len(s["votes"]) < c4["adjacency_vote_min"]:
            baskets[s["name"]].append("SELECTIVE_NETWORK_NEIGHBOR")
        if "white_dwarf" in s["flags"] and s["n_comp"] == 1:
            baskets[s["name"]].append("COMPACT_LENS_WILDCARD")
        if s["name"] in eng_scores and eng_scores[s["name"]] >= eng_cut:
            baskets[s["name"]].append("ENGINEERING_BACKBONE")
        s["science"] = None
        if s["gj"] in sci1:
            baskets[s["name"]].append("SCIENCE_INTEREST")
            s["science"] = {"tier": 1, **sci1[s["gj"]]}
        elif s["gj"] in sci2:
            baskets[s["name"]].append("SCIENCE_INTEREST_T2")
            s["science"] = {"tier": 2, **sci2[s["gj"]]}
        s["eng_score"] = eng_scores.get(s["name"])
    hist = sorted([s for s in systems if s["d_min_pc"] is not None],
                  key=lambda s: s["d_min_pc"])
    n_hist = 0
    for s in hist:
        if s["name"] not in baskets and n_hist < c4["n_historical"]:
            baskets[s["name"]].append("HISTORICAL_NEIGHBOR")
            n_hist += 1
        elif s["name"] in baskets and s["d_min_pc"] < 0.75 * s["dist_pc"]:
            baskets[s["name"]].append("HISTORICAL_NEIGHBOR")

    selected = [s for s in systems if s["name"] in baskets]

    # Component tracks: one per stellar member with astrometry; members
    # of a close pair share a corridor (marginal cost only).
    for s in selected:
        s["n_tracks"] = sum(comp_count(m) for m in s["members"])
    # Trim to budget: never from the distance core. Tier-2 science
    # systems enter last; otherwise more independent baskets first.
    def trim_key(s):
        b = set(baskets[s["name"]])
        core = "DISTANCE_CORE" in b
        t2only = b == {"SCIENCE_INTEREST_T2"}
        strong = len(b - {"SCIENCE_INTEREST_T2"})
        return (0 if core else 1, 1 if t2only else 0, -strong,
                -len(s["votes"]), s["dist_pc"])
    selected.sort(key=trim_key)
    kept, n_tracks = [], 0
    for s in selected:
        nt = s["n_tracks"]
        if n_tracks + nt > c4["track_budget"] and \
                "DISTANCE_CORE" not in baskets[s["name"]]:
            baskets[s["name"]].append("DEFERRED_BUDGET")
            continue
        kept.append(s)
        n_tracks += nt
    kept.sort(key=lambda s: s["dist_pc"])

    out = []
    for s in kept:
        out.append({
            "system": s["name"], "gj": s["gj"],
            "dist_pc": round(s["dist_pc"], 3),
            "n_comp": s["n_comp"],
            "min_sep_au": (round(s["min_sep_au"], 1)
                           if s["min_sep_au"] else None),
            "multiplicity": s["cls"], "flags": s["flags"],
            "adjacency_votes": len(s["votes"]),
            "vote_graphs": s["votes"],
            "selective_votes": s["selective_votes"],
            "engineering": s["eng"],
            "engineering_score": s["eng_score"],
            "science": s["science"],
            "dmin_1myr_pc": (round(s["d_min_pc"], 2)
                             if s["d_min_pc"] is not None else None),
            "baskets": sorted(set(baskets[s["name"]])),
            **s["search"],
            "n_tracks": s["n_tracks"],
            "components": [
                {"comp": m["comp"] or "", "name": resolve_name(m),
                 "gaia_dr3": m["gaia"], "gmag": m["g"],
                 "tracks": comp_count(m)}
                for m in s["members"]],
        })
    result = {"config": CFG, "n_systems": len(kept),
              "n_tracks": sum(s["n_tracks"] for s in kept),
              "systems": out}
    ver = CFG["version"]
    (TDIR / f"universal_{ver}.json").write_text(json.dumps(result, indent=1))
    write_markdown(result, kept, baskets, TDIR / f"universal_{ver}.md")
    print(f"selected {len(kept)} systems, "
          f"{result['n_tracks']} component tracks")
    for s in out:
        print(f"  {s['system']:24s} {s['dist_pc']:6.2f} pc "
              f"votes={s['adjacency_votes']} gate={s['solution_gate']:12s} "
              f"{'+'.join(s['baskets'])}")


BASKET_LABELS = {
    "DISTANCE_CORE": "Distance Core",
    "ROBUST_GEOMETRIC_NEIGHBOR": "Robust Geometric Neighbor",
    "SELECTIVE_NETWORK_NEIGHBOR": "Selective Network Neighbor",
    "ENGINEERING_BACKBONE": "Engineering Backbone",
    "SCIENCE_INTEREST": "Science Interest",
    "SCIENCE_INTEREST_T2": "Science Interest (tier 2)",
    "HISTORICAL_NEIGHBOR": "Historical Neighbor",
    "COMPACT_LENS_WILDCARD": "Compact Lens Wildcard",
    "DEFERRED_BUDGET": "Deferred (budget)",
}


def write_markdown(result, kept, baskets, path):
    import datetime as dt
    ver = CFG["version"]
    rows = result["systems"]
    n_new = sum(1 for r in rows if not r["already_searched_wise"])
    L = [f"---", f'title: "Universal SGL target list — {ver}"',
         f"date: {dt.date.today().isoformat()}",
         f'status: "{result["n_systems"]} systems / {result["n_tracks"]} '
         f'component tracks; {n_new} systems not yet searched in WISE"',
         "---", "", f"# Universal target list {ver}", "",
         "Survey-agnostic target portfolio built by "
         "`scripts/build_universal_list.py` from `config.yaml` and the "
         "curated `science_interest.yaml`. Census: CNS5 within "
         f"{CFG['census']['horizon_pc']:g} pc (`{CFG['census']['snapshot']}`), "
         f"engineering joins in `{CFG['census'].get('engineering', '')}`. "
         "The network prior, the engineering/desirability columns, and "
         "archive searchability are separate columns, never merged.", "",
         "Baskets: **Distance Core** (nearest "
         f"{CFG['stage4']['n_distance_core']} systems, protected), "
         "**Robust Geometric Neighbor** (≥"
         f"{CFG['stage4']['adjacency_vote_min']} of 7 sparse-graph families "
         "among all systems), **Selective Network Neighbor** (same vote on "
         "the *desirable* subgraph — level 1 drops close multiples, "
         "substellar objects and evolved stars; level 2 also drops "
         "intermediate multiples, detected accelerations, active and "
         "rapidly rotating stars), **Engineering Backbone** (top "
         f"{CFG['stage4']['n_engineering']} rank-sum of compactness M/R², "
         "quietness, dynamical cleanliness, lifetime and distance among "
         "level-2 dwarfs), **Science Interest** (curated tier 1 / tier 2), "
         "**Historical Neighbor** (closest approach within ±"
         f"{CFG['stage1']['historical_window_myr']:g} Myr), "
         "**Compact Lens Wildcard** (isolated white dwarfs).", "",
         "Desirability is the network's presumed taste and is distinct "
         "from the solution gate (our ability to predict the track).", "",
         "| System | d (pc) | Tracks | Mult. | SpT | M/R² | Activity | Accel. | Votes all / L1 / L2 | Eng. score | Gate | WISE | Baskets |",
         "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        e = r["engineering"]
        sv = r["selective_votes"]
        lv = sorted(sv)
        votes = f"{r['adjacency_votes']} / " + " / ".join(str(sv[k]) for k in lv)
        L.append("| " + " | ".join([
            r["system"], f"{r['dist_pc']:.2f}", str(r["n_tracks"]),
            r["multiplicity"], e["sptype"] or "—",
            f"{e['compactness_m_r2']:.1f}" if e["compactness_m_r2"] else "—",
            e["activity"], e["acceleration"].replace("_", " "), votes,
            f"{r['engineering_score']:.2f}" if r["engineering_score"] else "—",
            r["solution_gate"], "done" if r["already_searched_wise"] else "queue",
            ", ".join(BASKET_LABELS.get(b, b) for b in r["baskets"])]) + " |")
    L += ["", "## Science-interest rationale", ""]
    for r in rows:
        if r["science"]:
            L.append(f"- **{r['system']}** (tier {r['science']['tier']}): "
                     f"{r['science']['reason']} — {r['science']['ref']}")
    deferred = [(n, b) for n, b in baskets.items() if "DEFERRED_BUDGET" in b]
    if deferred:
        L += ["", "## Deferred by track budget", ""]
        for n, b in deferred:
            L.append(f"- {n}: " + ", ".join(BASKET_LABELS.get(x, x)
                                            for x in b if x != "DEFERRED_BUDGET"))
    prev = TDIR / "universal_v1.json"
    if prev.is_file() and ver != "v1":
        old = {x["system"] for x in json.load(open(prev))["systems"]}
        old |= {"YZ Cet" if n == "GJ 54" else n for n in old}
        new = {r["system"] for r in rows}
        added = [r for r in rows if r["system"] not in old]
        dropped = sorted(n for n in old if n not in new and n != "GJ 54")
        L += ["", "## Changes from v1", "",
              f"Added ({len(added)}):", ""]
        for r in added:
            L.append(f"- {r['system']} ({r['dist_pc']:.2f} pc): "
                     + ", ".join(BASKET_LABELS.get(b, b) for b in r["baskets"]))
        L += ["", f"Dropped ({len(dropped)}): " + (", ".join(dropped) or "none"),
              "", "Renamed: GJ 54 → YZ Cet (GJ 54.1; v1 truncated decimal GJ ids)."]
    L += ["", "## Notes", "",
          "- Mass/radius provenance per system is in the JSON "
          "(`mass_radius_source`: mann_mks / tic / flame / wd_nominal / unknown).",
          "- Activity evidence strings record which indicator decided the "
          "class; `unknown` is not `quiet`. ROSAT non-detections are "
          "recorded as Lx/Lbol upper limits and count as quiet only when "
          "the limit is below the active threshold.",
          "- Adjacency votes are computed against all systems in the census "
          "(all / level-1 / level-2 subgraphs), so only ~15 systems can hold "
          "Sun-Delaunay edges at once; a vote of ≥3 is meaningful."]
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
