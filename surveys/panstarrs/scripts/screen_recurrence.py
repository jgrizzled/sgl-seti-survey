"""(PS1 port of the ZTF/WISE script.) First-cut recurrence analysis over ScreenMatch records (plan §3.4:
"sources that recur along the predicted parallax, anti-target-motion,
and role-dependent track").

A persistent relay at fixed relay distance z produces ~one match per
visit whose implied z is consistent across visits; chance background
matches scatter in z and epoch. For each endpoint x role we bin matches
in log-z (each match's z_segment marks every bin it overlaps, once per
visit) and compare each bin's visit-coverage count against a binomial
background estimated from all bins of that endpoint x role.

Output: per endpoint x role, the most-recurrent z bins with rough
binomial tail probabilities. These are NOT trials-corrected and NOT
candidates — plan §3.4 requires motion-model comparison and held-out
epoch prediction before anything is retained. This is a screening
triage that says where (if anywhere) stage-2 forced photometry should
look first.

Usage: uv run python surveys/panstarrs/scripts/screen_recurrence.py
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
SCREEN_DIR = REPO / "runs" / "panstarrs" / "screen_v1"
VISIT_GAP_DAYS = 3.0  # TTI pairs + same-night visits collapse; revisits are weeks apart
N_ZBINS = 24
Z_MIN, Z_MAX = 550.0, 10000.0


def binom_tail(n: int, k: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    from math import comb
    p = min(max(p, 1e-12), 1 - 1e-12)
    return float(sum(comb(n, j) * p**j * (1 - p)**(n - j)
                     for j in range(k, n + 1)))


def gap_cluster_ids(mjds: np.ndarray, gap: float) -> np.ndarray:
    order = np.argsort(mjds)
    ids = np.zeros(len(mjds), dtype=int)
    cur = 0
    last = None
    for idx in order:
        if last is not None and mjds[idx] - last > gap:
            cur += 1
        ids[idx] = cur
        last = mjds[idx]
    return ids


def main() -> None:
    matches = [json.loads(l) for l in
               open(SCREEN_DIR / "records" / "screen_match.jsonl")]
    print(f"{len(matches)} screen matches")
    edges = np.geomspace(Z_MIN, Z_MAX, N_ZBINS + 1)

    by_pair = defaultdict(list)
    for m in matches:
        by_pair[(m["endpoint_id"], m["role"])].append(m)

    report = {}
    for (endpoint, role), ms in sorted(by_pair.items()):
        mjds = np.array([m["mjd"] for m in ms])
        visit_ids = gap_cluster_ids(mjds, VISIT_GAP_DAYS)
        n_visits = int(visit_ids.max()) + 1
        # bin -> set of visits with >= 1 match overlapping that z bin
        bin_visits = defaultdict(set)
        for m, v in zip(ms, visit_ids):
            z0, z1 = sorted(m["z_segment_au"])
            lo = int(np.searchsorted(edges, z0, side="right") - 1)
            hi = int(np.searchsorted(edges, z1, side="left"))
            for b in range(max(lo, 0), min(hi, N_ZBINS - 1) + 1):
                bin_visits[b].add(int(v))
        counts = np.array([len(bin_visits.get(b, ()))
                           for b in range(N_ZBINS)])
        p_bg = float(counts.mean()) / max(n_visits, 1)
        rows = []
        for b in np.argsort(counts)[::-1][:3]:
            k = int(counts[b])
            rows.append({
                "z_bin_au": [round(float(edges[b]), 1),
                             round(float(edges[b + 1]), 1)],
                "visits_with_match": k, "n_visits": n_visits,
                "background_p_per_visit": round(p_bg, 3),
                "binom_tail_uncorrected": float(
                    f"{binom_tail(n_visits, k, p_bg):.3g}"),
            })
        report[f"{endpoint}/{role}"] = {
            "matches": len(ms), "visits": n_visits,
            "mean_bin_coverage": round(p_bg, 3), "top_bins": rows}
        top = rows[0]
        print(f"{endpoint:16s} {role}: {len(ms):6d} matches, "
              f"{n_visits:3d} visits, mean bin coverage {p_bg:.2f}; "
              f"top z-bin {top['z_bin_au']} in "
              f"{top['visits_with_match']}/{n_visits} visits "
              f"(P>= {top['binom_tail_uncorrected']:.3g})")

    # ---- Stage 1.5: fixed-z point-matched filter --------------------
    # Bin occupancy saturates at this background density (a 10" strip
    # along the whole arc catches chance matches in every visit). A
    # relay at fixed z predicts a POINT per epoch; requiring matches
    # within R_TIGHT of the interpolated locus point at each trial z
    # shrinks the chance area by ~two orders of magnitude.
    R_TIGHT_ARCSEC = 1.5  # ~1 seeing FWHM at PS1
    Z_GRID = 1.0 / np.linspace(1.0 / Z_MIN, 1.0 / Z_MAX, 240)

    from astropy.time import Time

    from sglseti import Role, adaptive_locus, load_target_registry
    from sglsurvey.geometry import GeometryContext

    registry = load_target_registry(REPO / "registries"
                                    / "pilot_wise_2026.yaml")
    ctx = GeometryContext.ps1_default()
    interp_cache: dict = {}

    def points_at_zgrid(endpoint, role, locus_mjd):
        key = (endpoint, role, locus_mjd)
        if key not in interp_cache:
            al = adaptive_locus(
                target=registry[endpoint], role=Role(role),
                observation_time=Time(locus_mjd, format="mjd"),
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=0.5, ephemeris=ctx.ephemeris,
                model=ctx.model)
            zs = np.array([p.z_au for p in al.points])
            ra = np.array([p.icrs_ra_deg for p in al.points])
            dec = np.array([p.icrs_dec_deg for p in al.points])
            # Interpolate linearly in q = 1/z (the sampler's variable).
            q = 1.0 / zs
            order = np.argsort(q)
            qg = 1.0 / Z_GRID
            ra_g = np.interp(qg, q[order], ra[order])
            dec_g = np.interp(qg, q[order], dec[order])
            interp_cache[key] = np.stack([ra_g, dec_g], axis=1)
        return interp_cache[key]

    def _supports(m, endpoint, role, z, r, patz, b):
        pts = patz(endpoint, role, m["locus_mjd"])
        cosd = np.cos(np.deg2rad(m["dec_deg"]))
        d = np.hypot((pts[b, 0] - m["ra_deg"]) * cosd,
                     pts[b, 1] - m["dec_deg"]) * 3600.0
        return d <= r

    print(f"\n=== fixed-z point filter (R={R_TIGHT_ARCSEC}\") ===")
    point_report = {}
    for (endpoint, role), ms in sorted(by_pair.items()):
        mjds = np.array([m["mjd"] for m in ms])
        visit_ids = gap_cluster_ids(mjds, VISIT_GAP_DAYS)
        n_visits = int(visit_ids.max()) + 1
        support = np.zeros((len(Z_GRID), n_visits), dtype=bool)
        for m, v in zip(ms, visit_ids):
            pts = points_at_zgrid(endpoint, role, m["locus_mjd"])
            cosd = np.cos(np.deg2rad(m["dec_deg"]))
            d = np.hypot((pts[:, 0] - m["ra_deg"]) * cosd,
                         pts[:, 1] - m["dec_deg"]) * 3600.0
            support[d <= R_TIGHT_ARCSEC, int(v)] = True
        k = support.sum(axis=1)
        p_bg = float(k.mean()) / max(n_visits, 1)
        b = int(np.argmax(k))
        tail = binom_tail(n_visits, int(k[b]), max(p_bg, 1e-6))

        # Static-star discriminator: supporting matches of a genuine
        # track-follower move with the predicted parallax+drift path;
        # a static background star clumps at one fixed sky position.
        sup = [(m, int(v)) for m, v in zip(ms, visit_ids)
               if _supports(m, endpoint, role, float(Z_GRID[b]),
                            R_TIGHT_ARCSEC, points_at_zgrid, b)]
        ras = np.array([m["ra_deg"] for m, _ in sup])
        decs = np.array([m["dec_deg"] for m, _ in sup])
        cosd = np.cos(np.deg2rad(decs.mean()))
        sky_rms = float(np.hypot(
            (ras - ras.mean()) * cosd, decs - decs.mean()
        ).std() * 3600.0) if len(sup) > 1 else None
        # Track motion scale over the supporting epochs: predicted point
        # spread at best z (parallax + secular drift).
        pred = np.array([points_at_zgrid(endpoint, role,
                                         m["locus_mjd"])[b]
                         for m, _ in sup])
        track_rms = float(np.hypot(
            (pred[:, 0] - pred[:, 0].mean()) * cosd,
            pred[:, 1] - pred[:, 1].mean()).std() * 3600.0) \
            if len(sup) > 1 else None
        verdict = "indeterminate"
        if sky_rms is not None and track_rms is not None:
            verdict = ("static-background-star"
                       if sky_rms < max(1.5, 0.3 * track_rms)
                       else "track-consistent-or-mixed")

        point_report[f"{endpoint}/{role}"] = {
            "r_tight_arcsec": R_TIGHT_ARCSEC,
            "mean_visits_per_z": round(float(k.mean()), 2),
            "n_visits": n_visits,
            "best_z_au": round(float(Z_GRID[b]), 1),
            "best_visits_with_match": int(k[b]),
            "binom_tail_uncorrected": float(f"{tail:.3g}"),
            "support_sky_rms_arcsec": (round(sky_rms, 2)
                                       if sky_rms is not None else None),
            "predicted_track_rms_arcsec": (round(track_rms, 2)
                                           if track_rms is not None
                                           else None),
            "verdict": verdict,
        }
        print(f"{endpoint:16s} {role}: mean {k.mean():5.2f}/{n_visits} "
              f"visits per trial z; best z={Z_GRID[b]:7.1f} AU with "
              f"{k[b]}/{n_visits} (P>= {tail:.3g}); "
              f"sky rms {sky_rms if sky_rms is None else round(sky_rms, 1)}\" "
              f"vs track rms {track_rms if track_rms is None else round(track_rms, 1)}\" "
              f"-> {verdict}")

    out = SCREEN_DIR / "recurrence_report.json"
    out.write_text(json.dumps(
        {"bin_occupancy": report, "fixed_z_point_filter": point_report},
        indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
