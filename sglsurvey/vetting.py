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
from typing import Sequence

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
