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
  Stage 3  Universal searchability columns: antipode coordinates,
           |b_gal|, |beta_ecl|, corridor drift rate, solution gate,
           already-searched status.
  Stage 4  Basket union -> labeled system list -> component tracks,
           trimmed to the track budget (never from the distance core).

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
                        "Ross 128", "eps Ind", "tau Cet", "GJ 54",
                        "Teegarden", "Lacaille 8760", "van Maanen",
                        "GJ 908", "GJ 784",
                        # batch 3 (2026-08-19)
                        "eps Eri", "Lacaille 9352", "GJ 1061",
                        "GJ 12724", "Wolf 1061",
                        # batch 4 + unblocked deferred (2026-08-20)
                        "61 Cyg", "Struve 2398", "Groombridge 34",
                        "GJ 1111", "Luyten", "Kapteyn", "LP 145-141",
                        "GJ 1221", "GJ 9193", "GJ 783", "GJ 11068",
                        "WISE 0855", "EZ Aqr", "Luhman 16", "Procyon"]

KM_S_TO_PC_MYR = 1.0227

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
}
CNS5_NAMES = {"2653": "Luhman 16 AB", "2194": "WISE 0855-0714",
              "723": "Teegarden's Star"}


def resolve_name(member):
    gj = member["gj"].split(".")[0].strip() if member["gj"] else ""
    gj_key = gj.strip()
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
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(rows[i]["xyz"] - rows[j]["xyz"]) < thr:
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
    if prim["g"] is None and (prim["k"] is None or prim["k"] > 8.0):
        # No Gaia G and faint in Ks: substellar/ultracool. Bright stars
        # missing Gaia G (Sirius, Procyon, alpha Cen merged rows) are
        # saturated, not substellar.
        flags.append("substellar_or_faint")
    else:
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
        s["search"] = searchability(s, s["cls"], s["flags"], s["orbit"])

    # Basket C: selective re-vote after dropping flagged hosts.
    drop = set(CFG["stage4"]["selective_drop_classes"])
    keep = [s for s in systems if not (set(s["flags"]) & drop)]
    sel_votes, _ = adjacency_votes(keep)
    keep_names = [s["name"] for s in keep]
    sel_map = {keep_names[i]: v for i, v in sel_votes.items()}

    c4 = CFG["stage4"]
    baskets = defaultdict(list)
    for rank, s in enumerate(systems):
        if rank < c4["n_distance_core"]:
            baskets[s["name"]].append("DISTANCE_CORE")
        if len(s["votes"]) >= c4["adjacency_vote_min"]:
            baskets[s["name"]].append("ROBUST_GEOMETRIC_NEIGHBOR")
        sv = sel_map.get(s["name"], [])
        if len(sv) >= c4["adjacency_vote_min"] and \
                len(s["votes"]) < c4["adjacency_vote_min"]:
            baskets[s["name"]].append("SELECTIVE_NETWORK_NEIGHBOR")
        if "white_dwarf" in s["flags"] and s["n_comp"] == 1:
            baskets[s["name"]].append("COMPACT_LENS_WILDCARD")
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
    # Trim to budget: never from the distance core.
    def trim_key(s):
        core = "DISTANCE_CORE" in baskets[s["name"]]
        return (0 if core else 1, -len(s["votes"]), s["dist_pc"])
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
            "system": s["name"], "dist_pc": round(s["dist_pc"], 3),
            "n_comp": s["n_comp"],
            "min_sep_au": (round(s["min_sep_au"], 1)
                           if s["min_sep_au"] else None),
            "multiplicity": s["cls"], "flags": s["flags"],
            "adjacency_votes": len(s["votes"]),
            "vote_graphs": s["votes"],
            "selective_votes": len(sel_map.get(s["name"], [])),
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
    (TDIR / "universal_v1.json").write_text(json.dumps(result, indent=1))
    print(f"selected {len(kept)} systems, "
          f"{result['n_tracks']} component tracks")
    for s in out:
        print(f"  {s['system']:24s} {s['dist_pc']:6.2f} pc "
              f"votes={s['adjacency_votes']} gate={s['solution_gate']:12s} "
              f"{'+'.join(s['baskets'])}")


if __name__ == "__main__":
    main()
