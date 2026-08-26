"""MAST gPhoton photon-database adapter (GALEX crossings survey).

The GALEX substrate is not exposure products but a photon/aspect
database behind the MAST Mashup SQL service (recon
``surveys/galex-crossings/notes/galex_recon_2026-08-26.md``): free-form
SQL against the ``aspect`` table (1 row per second per band, boresight
ra0/dec0, flag), the ``NUVPhotonsV``/``FUVPhotonsV`` photon-event views
(5 ms-tick times, aspect-corrected ra/dec, flag), and the GR6+7 MCAT.
Endpoint facts pinned at recon: the service lives on ``mastcomp``
(the ``mast.`` host 404s); ``fGetTimeRanges`` does not exist; the
``band`` column is ``'FUV/NUV'`` when both detectors are on
(substring-match). Every query is snapshotted verbatim through the
caller's SnapshotStore.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import requests

MASHUP_URL = ("https://mastcomp.stsci.edu/portal/Mashup/MashupQuery.asmx/"
              "GalexPhotonListQueryTest")
#: "GALEX time" = UNIX - offset, seconds; DB time columns are ms.
GALEX_UNIX_OFFSET_S = 315_964_800
#: visits reconstruct as contiguous aspect-second runs; a gap larger
#: than this splits (aspect cadence is 1 s).
VISIT_GAP_MS = 2_000

MJD_UNIX_EPOCH = 40587.0  # MJD of 1970-01-01T00:00:00 UTC


def mjd_to_galex_ms(mjd_utc: float) -> float:
    return ((mjd_utc - MJD_UNIX_EPOCH) * 86400.0
            - GALEX_UNIX_OFFSET_S) * 1000.0


def galex_ms_to_mjd(t_ms: float) -> float:
    return (t_ms / 1000.0 + GALEX_UNIX_OFFSET_S) / 86400.0 + MJD_UNIX_EPOCH


def galex_ms_to_iso(t_ms: float) -> str:
    return datetime.fromtimestamp(
        t_ms / 1000.0 + GALEX_UNIX_OFFSET_S,
        tz=timezone.utc).isoformat()


class GphotonError(RuntimeError):
    pass


class GphotonClient:
    """Snapshot-disciplined SQL access to the gPhoton database."""

    def __init__(self, timeout_s: float = 600.0):
        self.timeout_s = timeout_s
        self.session = requests.Session()

    def query(self, sql: str, store=None):
        """Run one Mashup SQL query; return (columns, rows).

        Snapshots the verbatim response when a SnapshotStore is given.
        Raises GphotonError on a non-COMPLETE service status.
        """
        request_utc = datetime.now(timezone.utc).isoformat()
        resp = self.session.get(MASHUP_URL,
                                params={"query": sql, "format": "extjs"},
                                timeout=self.timeout_s)
        payload = json.loads(resp.text)
        status = payload.get("status")
        tables = (payload.get("data") or {}).get("Tables") or []
        rows = tables[0]["Rows"] if tables else []
        cols = ([c["text"] for c in tables[0]["Columns"]]
                if tables else [])
        if store is not None:
            store.store(service_url=MASHUP_URL, query=sql,
                        request_utc=request_utc,
                        response_bytes=resp.content,
                        row_count=len(rows) if status == "COMPLETE" else None,
                        http_status=resp.status_code)
        if resp.status_code != 200 or status != "COMPLETE":
            raise GphotonError(
                f"gPhoton query failed (HTTP {resp.status_code}, "
                f"status {status!r}): {payload.get('msg')!r} :: {sql[:200]}")
        return cols, rows

    # -- aspect ---------------------------------------------------------

    def aspect_near(self, ra_deg: float, dec_deg: float,
                    radius_arcmin: float, t0_ms: float, t1_ms: float,
                    store=None):
        """(time_ms, band, distance_arcmin) aspect rows near a position."""
        ra_deg, dec_deg = float(ra_deg), float(dec_deg)
        radius_arcmin = float(radius_arcmin)
        sql = (f"select time, band, distance from fGetNearbyAspectEq("
               f"{ra_deg!r},{dec_deg!r},{radius_arcmin!r},"
               f"{int(t0_ms)},{int(t1_ms)}) "
               f"where time >= {int(t0_ms)} and time < {int(t1_ms)}")
        _, rows = self.query(sql, store)
        return [(int(r[0]), str(r[1]), float(r[2])) for r in rows]

    def aspect_range(self, t0_ms: float, t1_ms: float, store=None):
        """Full aspect rows (time, ra0, dec0, band, flag) for a range.

        Use only on visit-scale ranges (~1 row/s/band): the aspect
        table is global, not per-field.
        """
        sql = (f"select time, ra0, dec0, band, flag from aspect "
               f"where time >= {int(t0_ms)} and time < {int(t1_ms)}")
        _, rows = self.query(sql, store)
        return [(int(r[0]), float(r[1]), float(r[2]), str(r[3]), int(r[4]))
                for r in rows]

    # -- photons --------------------------------------------------------

    def photons_box(self, band: str, ra_deg: float, dec_deg: float,
                    half_deg: float, t0_ms: float, t1_ms: float,
                    store=None):
        """Photon events (time, ra, dec, flag) in a time-range + box."""
        view = {"NUV": "NUVPhotonsV", "FUV": "FUVPhotonsV"}[band]
        ra_deg, dec_deg, half_deg = (float(ra_deg), float(dec_deg),
                                     float(half_deg))
        sql = (f"select time, ra, dec, flag from {view} "
               f"where time >= {int(t0_ms)} and time < {int(t1_ms)} "
               f"and ra between {ra_deg - half_deg!r} and "
               f"{ra_deg + half_deg!r} "
               f"and dec between {dec_deg - half_deg!r} and "
               f"{dec_deg + half_deg!r}")
        _, rows = self.query(sql, store)
        return [(int(r[0]), float(r[1]), float(r[2]), int(r[3]))
                for r in rows]

    # -- mcat -----------------------------------------------------------

    def mcat_box(self, ra_deg: float, dec_deg: float, half_deg: float,
                 store=None, limit: int = 500):
        """photoobjall rows (objid, nuv_mag, fuv_mag, ra, dec) in a box."""
        ra_deg, dec_deg, half_deg = (float(ra_deg), float(dec_deg),
                                     float(half_deg))
        sql = (f"select top {limit} objid, nuv_mag, fuv_mag, ra, dec "
               f"from GR6Plus7.dbo.photoobjall "
               f"where ra between {ra_deg - half_deg!r} and "
               f"{ra_deg + half_deg!r} "
               f"and dec between {dec_deg - half_deg!r} and "
               f"{dec_deg + half_deg!r}")
        _, rows = self.query(sql, store)
        return rows


def group_visits(times_ms):
    """Contiguous (start_ms, stop_ms) runs from sorted aspect seconds."""
    times = sorted(set(int(t) for t in times_ms))
    if not times:
        return []
    runs = []
    start = prev = times[0]
    for t in times[1:]:
        if t - prev > VISIT_GAP_MS:
            runs.append((start, prev))
            start = t
        prev = t
    runs.append((start, prev))
    return runs
