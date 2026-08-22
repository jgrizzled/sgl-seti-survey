"""Shared candidate-vetting helpers (plan §3.4 layer-1 screening applied
to stack-layer exceedances).

The catalogued-static-source test: a background star sitting within a
PSF width of the predicted track at the epochs that dominate a stack
reproduces a high real-track S with a chance minor-phase S, and passed
the automatic phase rules in 2 of 6 PS1 and 2 of 2 joint PS1+ZTF
retained cells (2026-08-21). It is cheap (catalog cones are already
snapshotted by each adapter's screening stage), deterministic, and
belongs in every adapter's automatic candidate rules rather than in a
manual stage-7 pass.

Catalog loaders read the screening snapshots offline where they exist
(PS1 DR2 ``mean``, ZTF DR24 ``objects``) and otherwise query CatWISE2020
through VizieR (CDS, not IRSA), caching the response under the run dir.
"""

from __future__ import annotations

import csv
import io
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from astropy.time import Time

from sglseti import Role, adaptive_locus

#: Default test radius per archive, ~1 PSF FWHM.
STATIC_RADIUS_ARCSEC = {"ps1": 2.0, "ztf": 2.5, "wise": 6.0, "spherex": 6.0}
TEST_QUANTILES = (10, 50, 90)
#: A catalog entry counts as a static source only with this many
#: detections: single-detection "objects" are transients, artefacts or
#: the very excess under test.
MIN_CATALOG_DETECTIONS = 3


@dataclass
class StaticCatalog:
    label: str
    ra: np.ndarray
    dec: np.ndarray
    mag: np.ndarray | None = None  # any single magnitude, for the report
    ndet: np.ndarray | None = None

    def nearest(self, ra_deg: float, dec_deg: float,
                min_ndet: int = MIN_CATALOG_DETECTIONS):
        if len(self.ra) == 0:
            return None
        cosd = np.cos(np.deg2rad(dec_deg))
        sep = np.hypot((self.ra - ra_deg) * cosd, self.dec - dec_deg) * 3600.0
        if self.ndet is not None:
            sep = np.where(self.ndet >= min_ndet, sep, np.inf)
            if not np.isfinite(sep).any():
                return None
        k = int(np.argmin(sep))
        return {"sep_arcsec": round(float(sep[k]), 2),
                "mag": (None if self.mag is None or not np.isfinite(self.mag[k])
                        or self.mag[k] <= 0 else round(float(self.mag[k]), 2)),
                "ndet": None if self.ndet is None else int(self.ndet[k]),
                "catalog": self.label}


def track_position(ctx, target, role: str, z_au: float, mu: Sequence[float],
                   t0_mjd: float, mjd: float) -> tuple[float, float]:
    """ICRS position (deg) of the (z, mu, T0) trajectory at ``mjd``."""
    al = adaptive_locus(target=target, role=Role(role),
                        observation_time=Time(mjd, format="mjd"),
                        observer=ctx.observer, relay_range=ctx.relay_range,
                        tolerance_arcsec=0.5, ephemeris=ctx.ephemeris,
                        model=ctx.model)
    zs = np.array([p.z_au for p in al.points])
    ra = np.array([p.icrs_ra_deg for p in al.points])
    dec = np.array([p.icrs_dec_deg for p in al.points])
    q = 1.0 / zs
    o = np.argsort(q)
    ra0 = float(np.interp(1.0 / z_au, q[o], ra[o]))
    dec0 = float(np.interp(1.0 / z_au, q[o], dec[o]))
    dt = (mjd - t0_mjd) / 365.25
    return (ra0 + mu[0] * dt / 3600.0 / np.cos(np.deg2rad(dec0)),
            dec0 + mu[1] * dt / 3600.0)


def static_source_test(ctx, target, role: str, z_au: float, mu: Sequence[float],
                       t0_mjd: float, mjd: np.ndarray, phase: np.ndarray,
                       catalog: StaticCatalog, radius_arcsec: float,
                       quantiles: Sequence[int] = TEST_QUANTILES) -> dict:
    """Nearest catalogued source to the track at the epoch quantiles of
    each populated parallax phase. ``static`` is True when, at the phase
    that carries the most epochs, a source lies within ``radius_arcsec``
    at a majority of the test epochs — the signature of a background
    star driving the stack."""
    mjd, phase = np.asarray(mjd, float), np.asarray(phase)
    per_phase = {}
    for p in sorted(set(phase.tolist())):
        sel = mjd[phase == p]
        if len(sel) == 0:
            continue
        tests = []
        for q in quantiles:
            t = float(np.percentile(sel, q))
            ra, dec = track_position(ctx, target, role, z_au, mu, t0_mjd, t)
            near = catalog.nearest(ra, dec)
            tests.append({"mjd": round(t, 2), "ra": ra, "dec": dec, "nearest": near})
        per_phase[str(p)] = {"n_epochs": int(len(sel)), "tests": tests}
    if not per_phase:
        return {"static": False, "per_phase": {}, "radius_arcsec": radius_arcsec}
    major = max(per_phase, key=lambda k: per_phase[k]["n_epochs"])
    hits = [t["nearest"] for t in per_phase[major]["tests"]
            if t["nearest"] and t["nearest"]["sep_arcsec"] < radius_arcsec]
    static = len(hits) * 2 > len(per_phase[major]["tests"])
    closest = min((t["nearest"]["sep_arcsec"] for t in per_phase[major]["tests"]
                   if t["nearest"]), default=None)
    return {"static": bool(static), "major_phase": major,
            "closest_major_phase_arcsec": closest,
            "hit": hits[0] if hits else None,
            "radius_arcsec": radius_arcsec, "per_phase": per_phase}


def static_reason(res: dict) -> str:
    h = res["hit"] or {}
    return (f"catalogued static source ({h.get('catalog')}) {h.get('sep_arcsec')}\" "
            f"from the track at the major parallax phase"
            + (f", mag {h['mag']}" if h.get("mag") is not None else "")
            + (f", {h['ndet']} detections" if h.get("ndet") is not None else "")
            + f" (< {res['radius_arcsec']}\" test radius)")


# -- catalog loaders -----------------------------------------------------
def _snapshot_rows(run_dir: Path, snapshot_ids: Sequence[str]) -> list[dict]:
    from sglsurvey.records import read_records

    snaps = {s["snapshot_id"]: s for s in read_records(
        run_dir / "records" / "query_snapshot.jsonl")}
    rows = []
    for sid in snapshot_ids:
        p = run_dir / snaps[sid]["response_path"]
        rows.extend(csv.DictReader(io.StringIO(p.read_text())))
    return rows


def load_ps1_mean(screen_dir: Path, corridor: str, band: str = "r") -> StaticCatalog:
    """PS1 DR2 ``mean`` objects snapshotted by surveys/panstarrs catalog_screen."""
    stats = json.loads((screen_dir / "catalog_stats.json").read_text())
    rows = _snapshot_rows(screen_dir, stats[corridor]["mean_snapshots"])
    return StaticCatalog(
        label="ps1-dr2-mean",
        ra=np.array([float(r["raMean"]) for r in rows]),
        dec=np.array([float(r["decMean"]) for r in rows]),
        mag=np.array([float(r[f"{band}MeanPSFMag"]) for r in rows]),
        ndet=np.array([int(float(r["nDetections"])) for r in rows]))


def load_ztf_objects(screen_dir: Path, ra_deg: float, dec_deg: float,
                     table: str = "ztf_objects_dr24") -> StaticCatalog | None:
    """ZTF DR objects table snapshotted per corridor by surveys/ztf
    catalog_screen (one 0.2-deg cone per corridor); picks the snapshot
    whose cone centre is nearest (ra, dec)."""
    import re

    from sglsurvey.records import read_records

    best, bdist = None, 1e9
    for s in read_records(screen_dir / "records" / "query_snapshot.jsonl"):
        if table not in s["query"]:
            continue
        m = re.search(r"CIRCLE\('ICRS',([-\d.]+),([-\d.]+)", s["query"])
        if not m:
            continue
        d = np.hypot((float(m.group(1)) - ra_deg) * np.cos(np.deg2rad(dec_deg)),
                     float(m.group(2)) - dec_deg)
        if d < bdist:
            best, bdist = s, d
    if best is None or bdist > 0.3:
        return None
    rows = list(csv.DictReader(io.StringIO(
        (screen_dir / best["response_path"]).read_text())))
    return StaticCatalog(
        label=table,
        ra=np.array([float(r["ra"]) for r in rows]),
        dec=np.array([float(r["dec"]) for r in rows]),
        mag=np.array([float(r["medianmag"] or "nan") for r in rows]),
        ndet=np.array([int(float(r["ngoodobs"] or 0)) for r in rows]))


def load_catwise_vizier(ra_deg: float, dec_deg: float, radius_deg: float,
                        cache_dir: Path, timeout_s: float = 120.0
                        ) -> StaticCatalog:
    """CatWISE2020 (VizieR II/365/catwise, CDS — not IRSA) cone, cached
    under ``cache_dir`` as a verbatim TSV with its query."""
    import hashlib

    import requests

    params = {"-source": "II/365/catwise", "-c": f"{ra_deg:.6f} {dec_deg:.6f}",
              "-c.rs": f"{radius_deg * 3600:.1f}", "-out.max": "unlimited",
              "-out": "RA_ICRS,DE_ICRS,W1mproPM,W2mproPM,nW1"}
    key = hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()[:16]
    cache_dir.mkdir(parents=True, exist_ok=True)
    p = cache_dir / f"catwise_{key}.tsv"
    if not p.exists():
        r = requests.get("https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
                         params=params, timeout=timeout_s)
        r.raise_for_status()
        p.write_text("# " + json.dumps(params) + "\n" + r.text)
    lines = [l for l in p.read_text().splitlines() if l and not l.startswith("#")]
    rows = []
    if lines:
        hdr = lines[0].split("\t")
        for l in lines[3:]:  # header, units, dashes
            v = l.split("\t")
            if len(v) == len(hdr):
                rows.append(dict(zip(hdr, [x.strip() for x in v])))

    def f(r, k):
        try:
            return float(r[k])
        except (KeyError, ValueError):
            return np.nan

    return StaticCatalog(
        label="catwise2020-vizier",
        ra=np.array([f(r, "RA_ICRS") for r in rows]),
        dec=np.array([f(r, "DE_ICRS") for r in rows]),
        mag=np.array([f(r, "W1mproPM") for r in rows]),
        ndet=np.array([int(f(r, "nW1")) if np.isfinite(f(r, "nW1")) else 0 for r in rows]))


# -- v2 additions (wise v2_plan §6, §8.5) ---------------------------------
def track_position_xt(ctx, target, role: str, z_au: float, mu: Sequence[float],
                      t0_mjd: float, mjd: float, xt_arcsec: float = 0.0
                      ) -> tuple[float, float]:
    """As :func:`track_position` with a cross-track offset ``xt_arcsec``
    applied perpendicular to the local corridor tangent (the direction
    of increasing z at this epoch), for the §2.3 cross-track cells."""
    ra0, dec0 = track_position(ctx, target, role, z_au, mu, t0_mjd, mjd)
    if xt_arcsec == 0.0:
        return ra0, dec0
    q = 1.0 / z_au
    ra1, dec1 = track_position(ctx, target, role, 1.0 / (q * 1.02), mu, t0_mjd, mjd)
    cosd = np.cos(np.deg2rad(dec0))
    tx, ty = (ra1 - ra0) * cosd, dec1 - dec0
    n = np.hypot(tx, ty) or 1.0
    nx, ny = -ty / n, tx / n
    return (ra0 + xt_arcsec * nx / 3600.0 / cosd, dec0 + xt_arcsec * ny / 3600.0)


def parallax_phase_test(S_by_phase: Mapping[int, float],
                        n_by_phase: Mapping[int, int] | None = None,
                        single_phase_ratio: float = 0.25) -> dict:
    """Annotation (never a veto): the stack significance split by
    parallax phase. A persistent in-scope source contributes at both
    phases; a static star near one phase position gives high S at one
    phase only; an intermittent in-scope source can legitimately be
    single-phase, which is why this is an annotation (v2 plan §6)."""
    vals = {int(k): float(v) for k, v in S_by_phase.items() if np.isfinite(v)}
    if len(vals) < 2:
        return {"n_phases_populated": len(vals), "S_by_phase": vals,
                "signature": "single-phase-data", "ratio": None}
    hi = max(vals.values())
    lo = min(vals.values())
    ratio = (lo / hi) if hi > 0 else None
    sig = ("single-phase" if (ratio is not None and ratio < single_phase_ratio)
           or (lo <= 0 < hi) else "both-phases")
    return {"n_phases_populated": len(vals), "S_by_phase": vals,
            "n_by_phase": ({int(k): int(v) for k, v in n_by_phase.items()}
                           if n_by_phase else None),
            "ratio": None if ratio is None else round(ratio, 3),
            "signature": sig}


def radial_response_table(template: np.ndarray, kernel: np.ndarray,
                          oversample: int = 8, max_pix: float = 16.0):
    """Matched-filter response (flux units per unit source flux) of a
    unit-sum Gaussian ``kernel`` (detector pixels) to a PRF ``template``
    (oversampled by ``oversample``) as a function of radial offset in
    detector pixels: r_pix (N,), response (N,). Used by the flux-
    consistency test for both the core (static star on the track) and
    the PSF wings (bright-star halo)."""
    from scipy.signal import fftconvolve

    n = template.shape[0]
    c = (n - 1) // 2
    # bin the template to detector sampling at zero sub-pixel phase
    half = int(max_pix) + kernel.shape[0]
    ii = np.arange(-half, half + 1)
    u = c + oversample * ii
    grid = template[np.ix_(np.clip(u, 0, n - 1), np.clip(u, 0, n - 1))] * oversample ** 2
    grid = np.where((np.abs(u) < n)[:, None] & (np.abs(u) < n)[None, :], grid, 0.0)
    num = fftconvolve(grid, kernel[::-1, ::-1], mode="same")
    resp = num / (kernel ** 2).sum()
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    r = np.hypot(xx, yy).ravel()
    o = np.argsort(r)
    r, v = r[o], resp.ravel()[o]
    # radial average in 0.25-pixel bins
    bins = np.arange(0, max_pix + 0.25, 0.25)
    idx = np.digitize(r, bins) - 1
    out_r, out_v = [], []
    for k in range(len(bins) - 1):
        sel = idx == k
        if sel.any():
            out_r.append(0.5 * (bins[k] + bins[k + 1]))
            out_v.append(float(np.mean(v[sel])))
    return np.array(out_r), np.array(out_v)


def flux_consistency(catalog: StaticCatalog, track_ra: np.ndarray,
                     track_dec: np.ndarray, w: np.ndarray, phase: np.ndarray,
                     zp_ref: float, resp_r_arcsec: np.ndarray,
                     resp_v: np.ndarray, S_measured: float,
                     S_by_phase: Mapping[int, float] | None = None,
                     search_arcsec: float = 30.0, factor: float = 2.0,
                     min_ndet: int = MIN_CATALOG_DETECTIONS) -> dict:
    """Can the catalogued static sources near the track account for the
    measured stack S? For every catalogue source within ``search_arcsec``
    of any test epoch's track position, the predicted per-epoch matched-
    filter flux is 10^((zp_ref - mag)/2.5) x response(separation), and
    the predicted S is sum_e w_e f_e / sqrt(sum_e w_e) over the same
    capped weights as the real stack. ``consistent`` is True when the
    prediction accounts for the measurement within ``factor`` (and
    likewise in the dominant phase when ``S_by_phase`` is given). This,
    together with proximity, is what makes the static/halo rejection a
    calibrated veto (v2 plan §6); proximity alone is an annotation."""
    track_ra = np.asarray(track_ra, float)
    track_dec = np.asarray(track_dec, float)
    w = np.asarray(w, float)
    phase = np.asarray(phase)
    ok = np.isfinite(track_ra) & np.isfinite(track_dec) & (w > 0)
    out = {"n_sources_considered": 0, "S_pred": 0.0, "S_measured": float(S_measured),
           "consistent": False, "sources": [], "by_phase": {}}
    if not ok.any() or len(catalog.ra) == 0:
        return out
    dec0 = float(np.nanmean(track_dec[ok]))
    cosd = np.cos(np.deg2rad(dec0))
    # candidate sources: within search radius of the track's bounding box
    ra_min, ra_max = track_ra[ok].min(), track_ra[ok].max()
    de_min, de_max = track_dec[ok].min(), track_dec[ok].max()
    pad = search_arcsec / 3600.0
    sel = ((catalog.ra >= ra_min - pad / cosd) & (catalog.ra <= ra_max + pad / cosd)
           & (catalog.dec >= de_min - pad) & (catalog.dec <= de_max + pad))
    if catalog.ndet is not None:
        sel &= catalog.ndet >= min_ndet
    if catalog.mag is not None:
        sel &= np.isfinite(catalog.mag)
    idx = np.where(sel)[0]
    if len(idx) == 0:
        return out
    fpred = np.zeros(len(track_ra))
    contrib = []
    for k in idx:
        d = np.hypot((catalog.ra[k] - track_ra) * cosd,
                     catalog.dec[k] - track_dec) * 3600.0
        resp = np.interp(d, resp_r_arcsec, resp_v, right=0.0)
        f = 10 ** (0.4 * (zp_ref - float(catalog.mag[k]))) * resp
        f = np.where(ok, f, 0.0)
        fpred += f
        Sk = float((w * f).sum() / np.sqrt(w[ok].sum()))
        if Sk > 0.05 * max(abs(S_measured), 1.0):
            contrib.append({"mag": round(float(catalog.mag[k]), 2),
                            "min_sep_arcsec": round(float(np.nanmin(d[ok])), 2),
                            "S_pred": round(Sk, 2),
                            "ndet": None if catalog.ndet is None else int(catalog.ndet[k])})
    Bsum = w[ok].sum()
    S_pred = float((w * fpred).sum() / np.sqrt(Bsum)) if Bsum > 0 else 0.0
    out.update({"n_sources_considered": int(len(idx)), "S_pred": round(S_pred, 2),
                "sources": sorted(contrib, key=lambda s: -s["S_pred"])[:5]})
    consistent = S_pred >= abs(S_measured) / factor and S_measured > 0
    if S_by_phase:
        byp = {}
        for p, Sm in S_by_phase.items():
            selp = ok & (phase == int(p))
            Bp = w[selp].sum()
            Sp = float((w[selp] * fpred[selp]).sum() / np.sqrt(Bp)) if Bp > 0 else 0.0
            byp[int(p)] = {"S_pred": round(Sp, 2), "S_measured": round(float(Sm), 2)}
        out["by_phase"] = byp
        if byp:
            dom = max(byp, key=lambda p: abs(byp[p]["S_measured"]))
            consistent &= byp[dom]["S_pred"] >= abs(byp[dom]["S_measured"]) / factor
    out["consistent"] = bool(consistent)
    return out


def holdout_prediction_test(S_early_node: float, f_early: float,
                            S_late: float, f_late: float, n_late: int,
                            min_S_late: float = 3.0, min_flux_ratio: float = 0.3
                            ) -> dict:
    """Calibrated post-hoc held-out-epoch prediction test (v2 plan §1.8).

    Given the refit on the early epochs (argmax node, mean flux
    ``f_early``) and the forced photometry of the late epochs at that
    node, the test PASSES when the late stack is significant and the
    late flux is at least ``min_flux_ratio`` of the early flux. Its
    false-pass rate on null exceedances and true-pass rate on
    injections are measured, not assumed, and quoted with any result.
    """
    if n_late < MIN_EPOCHS_HOLDOUT or not np.isfinite(S_late):
        return {"applicable": False, "n_late": int(n_late)}
    ratio = (f_late / f_early) if f_early > 0 else np.nan
    passed = bool(S_late >= min_S_late and np.isfinite(ratio) and ratio >= min_flux_ratio)
    return {"applicable": True, "pass": passed, "S_early": round(float(S_early_node), 2),
            "S_late": round(float(S_late), 2), "f_early": float(f_early),
            "f_late": float(f_late), "flux_ratio": (None if not np.isfinite(ratio)
                                                   else round(float(ratio), 3)),
            "n_late": int(n_late), "rule": f"S_late >= {min_S_late} and f_late/f_early >= {min_flux_ratio}"}


MIN_EPOCHS_HOLDOUT = 5
