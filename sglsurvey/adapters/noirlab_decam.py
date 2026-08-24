"""NOIRLab Astro Archive DECam instcal adapter (fifth adapter; first
NOIRLab one — plan §7 TODO southern corridors / §10 step 8).

Discovery uses the Astro Archive advanced-search API (JSON POST,
snapshotted verbatim); image/dqmask/wtmap files of one exposure are
joined on the EXPNUM aux field. Products are fetched from
``api/retrieve/<md5>/``; single-CCD subsets use the ``?hdus=0,N``
parameter (plain HTTP Range is ignored by the service). Facts verified
before writing this module are in notes/decam_recon_2026-08-24.md.

Key archive facts this module relies on (all verified 2026-08-24):

- every instcal exposure has image (~321 MB), dqmask (~6 MB fpack) and
  wtmap (~154 MB) siblings sharing EXPNUM and MJD-OBS;
- all CCD HDUs of one exposure share CRVAL = the exposure pointing
  centre, with full TPV distortion per CCD — so a static tangent-plane
  focal-plane layout gives a nominal footprint from the pointing centre
  alone, and the dqmask alone gives exact astrometry + usable pixels;
- image/dqmask/wtmap of one exposure share HDU ordering, but the CCD
  set varies across eras (dead CCDs), so the EXTNAME->HDU-index map is
  derived per exposure from its dqmask and asserted after every
  single-HDU fetch;
- the Data Lab cutout service strips TPV to TAN — never use it for
  astrometry; ``?hdus=`` CCD fetches are the cutout mechanism here.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
import requests

from sglsurvey.adapters.base import (ConeRegion, CutoutSpec, LocalProduct,
                                     MjdRange, ProductSet)
from sglsurvey.records import Observation
from sglsurvey.snapshots import SnapshotStore

SEARCH_URL = "https://astroarchive.noirlab.edu/api/adv_search/find/"
RETRIEVE_URL = "https://astroarchive.noirlab.edu/api/retrieve/{md5}/"
COLLECTION = "decam/instcal"

#: DECam focal plane: 61-62 science CCDs inside a ~2.2 deg-wide hexagon.
#: Maximum angular distance from the pointing centre to any live pixel.
FOV_RADIUS_DEG = 1.1
PIX_ARCSEC = 0.2637

FILTER_WAVELENGTH_UM = {"u": 0.355, "g": 0.473, "r": 0.642, "i": 0.784,
                        "z": 0.926, "Y": 1.009, "VR": 0.626}

#: CP dqmask codes (ood product). Bits per the NOIRLab community-pipeline
#: instcal documentation; 4 (interpolated) is not fatal by default —
#: the hypothesis freeze owns the final convention.
DQ_BAD_DETECTOR = 1
DQ_SATURATED = 2
DQ_INTERPOLATED = 4
DQ_COSMIC_RAY = 16
DQ_BLEED = 64
DQ_TRANSIENT = 128
DQ_FATAL_DEFAULT = 0xFFFF & ~DQ_INTERPOLATED

_OUTFIELDS = ["md5sum", "archive_filename", "original_filename",
              "proc_type", "prod_type", "obs_type", "proposal", "caldat",
              "ifilter", "exposure", "ra_center", "dec_center",
              "dateobs_min", "EXPNUM", "MJD-OBS", "filesize", "url"]


def _band(ifilter: str | None) -> str:
    return (ifilter or "?").split()[0] if (ifilter or "").strip() else "?"


def _tangent_offsets(ra_deg, dec_deg, ra0_deg: float, dec0_deg: float):
    """Gnomonic (TAN) projection of points about (ra0, dec0), in deg.

    xi grows with RA (east), eta with Dec (north) — the convention the
    focal-plane layout is built in.
    """
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    ra0, dec0 = np.deg2rad(ra0_deg), np.deg2rad(dec0_deg)
    cosc = (np.sin(dec0) * np.sin(dec)
            + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0))
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = np.cos(dec) * np.sin(ra - ra0) / cosc
        eta = (np.cos(dec0) * np.sin(dec)
               - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    bad = cosc <= 0  # opposite hemisphere
    xi = np.where(bad, np.inf, xi)
    eta = np.where(bad, np.inf, eta)
    return np.rad2deg(xi), np.rad2deg(eta)


def build_focal_plane_layout(dqmask_path: Path) -> dict:
    """Static per-CCD tangent-plane bounding boxes, derived once from a
    reference dqmask (all CCDs share CRVAL = pointing centre; DECam has
    no field rotator, so the layout is fixed on the sky)."""
    from astropy.io import fits
    from astropy.wcs import WCS

    layout = {}
    with fits.open(dqmask_path) as hdul:
        ra0 = float(hdul[1].header["CRVAL1"])
        dec0 = float(hdul[1].header["CRVAL2"])
        for hdu in hdul[1:]:
            h = hdu.header
            nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
            wcs = WCS(h)
            corners = np.array([[0.5, 0.5], [nx + 0.5, 0.5],
                                [nx + 0.5, ny + 0.5], [0.5, ny + 0.5]])
            sky = wcs.all_pix2world(corners, 1)
            xi, eta = _tangent_offsets(sky[:, 0], sky[:, 1], ra0, dec0)
            layout[h["EXTNAME"]] = [float(xi.min()), float(xi.max()),
                                    float(eta.min()), float(eta.max())]
    return layout


class DecamNominalFootprint:
    """Containment against the static focal-plane layout placed at the
    exposure pointing centre. Coarse stage only: real dead pixels, era
    CCD outages and the TPV distortion residual are invisible here (the
    latter is < 1.5" and covered by the pad)."""

    def __init__(self, layout: dict, ra_center_deg: float,
                 dec_center_deg: float, pad_arcsec: float = 0.0):
        self._boxes = np.array(list(layout.values()))  # (n, 4)
        self._ra0, self._dec0 = ra_center_deg, dec_center_deg
        self._pad_deg = pad_arcsec / 3600.0

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        return self.contains_any(np.array([[ra_deg, dec_deg]]))

    def contains_any(self, radec_deg) -> bool:
        radec_deg = np.asarray(radec_deg, dtype=float)
        xi, eta = _tangent_offsets(radec_deg[:, 0], radec_deg[:, 1],
                                   self._ra0, self._dec0)
        p = self._pad_deg
        b = self._boxes
        inside = ((xi[:, None] >= b[None, :, 0] - p)
                  & (xi[:, None] <= b[None, :, 1] + p)
                  & (eta[:, None] >= b[None, :, 2] - p)
                  & (eta[:, None] <= b[None, :, 3] + p))
        return bool(np.any(inside))


class DecamExactFootprint:
    """Exact usable-pixel test over a full dqmask file: per-CCD TPV WCS
    + data-quality codes. A pixel is usable when it lies on a CCD and
    its dqmask value has no fatal bit.

    Mask data is decompressed lazily per CCD (a full exposure would be
    ~2 GB int32; a corridor locus touches 1-2 CCDs). Call ``close()``
    (or use as a context manager) when done."""

    def __init__(self, dqmask_path: Path,
                 fatal_mask: int = DQ_FATAL_DEFAULT):
        from astropy.io import fits
        from astropy.wcs import WCS

        self.fatal_mask = fatal_mask
        self._hdul = fits.open(dqmask_path)
        self.primary = dict(self._hdul[0].header)
        self._ccds = []  # (extname, hdu_index, wcs, nx, ny)
        self.hdu_index = {}
        self._masks: dict[str, np.ndarray] = {}
        pad = 5.0 / 3600.0  # bbox prefilter pad vs corner projection
        for i, hdu in enumerate(self._hdul[1:], start=1):
            h = hdu.header
            name = h["EXTNAME"]
            self.hdu_index[name] = i
            wcs = WCS(h)
            nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
            corners = wcs.all_pix2world(
                np.array([[0.5, 0.5], [nx + 0.5, 0.5],
                          [nx + 0.5, ny + 0.5], [0.5, ny + 0.5]]), 1)
            cosd = max(np.cos(np.deg2rad(corners[:, 1].mean())), 1e-3)
            bbox = (corners[:, 0].min() - pad / cosd,
                    corners[:, 0].max() + pad / cosd,
                    corners[:, 1].min() - pad,
                    corners[:, 1].max() + pad)
            self._ccds.append((name, i, wcs, nx, ny, bbox))

    def close(self) -> None:
        self._hdul.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def _mask(self, name: str) -> np.ndarray:
        if name not in self._masks:
            self._masks[name] = np.asarray(
                self._hdul[self.hdu_index[name]].data)
        return self._masks[name]

    def locate(self, radec_deg):
        """Per point: (extname or None, x, y, usable)."""
        radec_deg = np.asarray(radec_deg, dtype=float)
        out = [(None, -1, -1, False)] * len(radec_deg)
        todo = set(range(len(radec_deg)))
        for name, _i, wcs, nx, ny, bbox in self._ccds:
            if not todo:
                break
            idx = [j for j in sorted(todo)
                   if bbox[0] <= radec_deg[j, 0] <= bbox[1]
                   and bbox[2] <= radec_deg[j, 1] <= bbox[3]]
            if not idx:
                continue
            pix = wcs.all_world2pix(radec_deg[idx], 0)
            with np.errstate(invalid="ignore"):
                x = np.rint(pix[:, 0]).astype(int)
                y = np.rint(pix[:, 1]).astype(int)
                inb = (x >= 0) & (x < nx) & (y >= 0) & (y < ny)
            if not inb.any():
                continue
            mask = self._mask(name)
            for k, j in enumerate(idx):
                if inb[k]:
                    usable = (int(mask[y[k], x[k]]) & self.fatal_mask) == 0
                    out[j] = (name, int(x[k]), int(y[k]), usable)
                    todo.discard(j)
        return out

    def status(self, radec_deg):
        """Per-point (on_ccd, usable) boolean arrays for (N,2) deg —
        the FootprintTest-adjacent bulk form the precise pass consumes."""
        loc = self.locate(radec_deg)
        onccd = np.array([l[0] is not None for l in loc])
        usable = np.array([l[3] for l in loc])
        return onccd, usable

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        return bool(self.locate([[ra_deg, dec_deg]])[0][3])

    def ccd_of(self, ra_deg: float, dec_deg: float):
        """(extname, hdu_index) of the CCD containing the point, or
        (None, -1) if the point falls off every CCD."""
        name = self.locate([[ra_deg, dec_deg]])[0][0]
        return (name, self.hdu_index[name]) if name else (None, -1)


class DecamInstcalAdapter:
    archive_id = "noirlab-decam"
    collection = COLLECTION

    def __init__(self, layout: dict | None = None,
                 layout_path: Path | None = None,
                 session: requests.Session | None = None,
                 timeout_s: float = 600.0, page_limit: int = 50000):
        if layout is None and layout_path is not None:
            layout = json.loads(Path(layout_path).read_text())["ccds"]
        self._layout = layout
        self._session = session or requests.Session()
        self._timeout = timeout_s
        self._page_limit = page_limit

    # -- discovery ---------------------------------------------------------
    def search(self, region: ConeRegion, store: SnapshotStore,
               extra_search: list | None = None) -> list[dict]:
        """Paged advanced-search over a box that covers every exposure
        whose focal plane can touch the cone; raw responses snapshotted.
        Returns one row per FILE (3 prod_types per exposure)."""
        r_deg = region.radius_deg + FOV_RADIUS_DEG
        cosd = max(np.cos(np.deg2rad(region.dec_deg)), 1e-3)
        search = [["instrument", "decam"], ["proc_type", "instcal"],
                  ["ra_center", (region.ra_deg - r_deg / cosd) % 360.0,
                   (region.ra_deg + r_deg / cosd) % 360.0],
                  ["dec_center", max(region.dec_deg - r_deg, -90.0),
                   min(region.dec_deg + r_deg, 90.0)]]
        # NOTE a box crossing RA 0/360 would need two queries; assert
        # instead of silently returning nothing.
        assert (region.ra_deg - r_deg / cosd) >= 0.0 and \
               (region.ra_deg + r_deg / cosd) <= 360.0, \
            "RA wrap-around box not implemented"
        search.extend(extra_search or [])
        payload = {"outfields": _OUTFIELDS, "search": search}
        body = json.dumps(payload).encode()
        rows, offset = [], 0
        while True:
            url = (f"{SEARCH_URL}?limit={self._page_limit}"
                   + (f"&offset={offset}" if offset else ""))
            request_utc = datetime.now(timezone.utc).isoformat()
            resp = self._session.post(
                url, data=body, headers={"Content-Type": "application/json"},
                timeout=self._timeout)
            resp.raise_for_status()
            batch = resp.json()
            got = batch[1:]
            store.store(service_url=url, query=body.decode(),
                        request_utc=request_utc, response_bytes=resp.content,
                        row_count=len(got), http_status=resp.status_code)
            rows.extend(got)
            if len(got) < self._page_limit:
                break
            offset += len(got)
        return rows

    def discover(self, region: ConeRegion, time_range: MjdRange,
                 store: SnapshotStore) -> Iterator[Observation]:
        rows = self.search(region, store)
        by_exp: dict[int, dict[str, dict]] = {}
        for row in rows:
            expnum = row.get("EXPNUM")
            ptype = row.get("prod_type")
            if expnum is None or ptype not in ("image", "dqmask", "wtmap"):
                continue
            prods = by_exp.setdefault(int(expnum), {})
            # Several CP versions of one file may coexist; keep the
            # lexicographically latest archive_filename (…_v2 > …_v1).
            old = prods.get(ptype)
            if old is None or (row.get("archive_filename") or "") > \
                    (old.get("archive_filename") or ""):
                prods[ptype] = row
        for expnum in sorted(by_exp):
            prods = by_exp[expnum]
            if set(prods) != {"image", "dqmask", "wtmap"}:
                continue  # unusable without all three siblings
            img = prods["image"]
            mjd = img.get("MJD-OBS")
            if mjd is None:
                continue
            if not (time_range.start_mjd_utc <= mjd
                    <= time_range.stop_mjd_utc):
                continue
            yield self._observation(expnum, prods)

    @staticmethod
    def _cp_version(row: dict) -> str:
        name = (row.get("archive_filename") or "").rsplit("/", 1)[-1]
        stem = name.split(".fits")[0]
        return stem.rsplit("_", 1)[-1] if "_" in stem else "?"

    def _observation(self, expnum: int, prods: dict) -> Observation:
        img = prods["image"]
        band = _band(img.get("ifilter"))
        exptime = float(img.get("exposure") or 0.0)
        t0 = float(img["MJD-OBS"])
        half = exptime / 2.0 / 86400.0
        products = {
            ptype: {"md5": row["md5sum"],
                    "url": RETRIEVE_URL.format(md5=row["md5sum"]),
                    "filesize": row.get("filesize"),
                    "archive_filename": row.get("archive_filename")}
            for ptype, row in prods.items()}
        return Observation.build(
            archive_id=self.archive_id, collection=self.collection,
            release=self._cp_version(img),
            native_key={"expnum": expnum,
                        "image_md5": img["md5sum"]},
            band=band, wavelength_um=FILTER_WAVELENGTH_UM.get(band),
            t_start_mjd_utc=t0, t_mid_mjd_utc=t0 + half,
            t_stop_mjd_utc=t0 + 2 * half, exptime_s=exptime,
            corners_icrs_deg=(),
            wcs={"ra_center": float(img["ra_center"]),
                 "dec_center": float(img["dec_center"]),
                 "layout": "decam_focal_plane_v1"},
            quality_flags={"obs_type": img.get("obs_type")},
            products=products,
            extra={"proposal": img.get("proposal"),
                   "caldat": img.get("caldat"),
                   "ifilter": img.get("ifilter"),
                   "original_filename": img.get("original_filename"),
                   "dateobs_min": img.get("dateobs_min")},
            snapshot_id=None,  # search() snapshots are shared, not per-row
            discovered_utc=None,
        )

    # -- footprints --------------------------------------------------------
    def nominal_footprint(self, obs: Observation,
                          pad_arcsec: float = 0.0) -> DecamNominalFootprint:
        if self._layout is None:
            raise RuntimeError("focal-plane layout not loaded")
        return DecamNominalFootprint(self._layout, obs.wcs["ra_center"],
                                     obs.wcs["dec_center"], pad_arcsec)

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> DecamExactFootprint:
        dq = next(p for p in products.products if p.kind == "dqmask")
        return DecamExactFootprint(dq.path)

    # -- fetch -------------------------------------------------------------
    def fetch(self, obs: Observation, kinds: Sequence[str], dest: Path,
              *, cutout: CutoutSpec | None = None,
              dqmask_dir: Path | None = None,
              extname: str | None = None) -> ProductSet:
        """Download products. ``dqmask`` is always the full file (needed
        for the exact footprint and the per-exposure HDU map) and its
        archive md5 is verified. For ``image``/``wtmap`` a single CCD
        HDU is fetched via ``?hdus=``, selected either by ``extname``
        directly or by ``cutout`` (the CCD containing (ra, dec);
        size_pix is ignored — the unit of retrieval is the CCD; the
        centre may sit in a chip gap, so callers that know the covered
        locus should pass ``extname`` instead). The EXTNAME of the
        received HDU is asserted, and the checksum is a local sha256
        (the archive md5 covers only the full file)."""
        dest.mkdir(parents=True, exist_ok=True)
        expnum = obs.native_key["expnum"]
        products = []
        exact = None
        for kind in list(kinds):
            info = obs.products[kind]
            if kind == "dqmask" or (cutout is None and extname is None):
                path = dest / f"exp{expnum}_{kind}.fits.fz"
                if not path.exists():
                    resp = self._session.get(info["url"],
                                             timeout=self._timeout)
                    if resp.status_code == 404:
                        raise FileNotFoundError(info["url"])
                    resp.raise_for_status()
                    md5 = hashlib.md5(resp.content).hexdigest()
                    if md5 != info["md5"]:
                        raise IOError(
                            f"md5 mismatch for {info['url']}: "
                            f"got {md5}, archive says {info['md5']}")
                    path.write_bytes(resp.content)
                digest = hashlib.md5(path.read_bytes()).hexdigest()
                if digest != info["md5"]:
                    raise IOError(f"stale local file {path}")
                products.append(LocalProduct(
                    kind=kind, path=path, checksum=f"md5:{digest}"))
                continue
            # single-CCD fetch: needs the exposure's dqmask HDU map
            if exact is None:
                dq_local = [p for p in products if p.kind == "dqmask"]
                if dq_local:
                    exact = DecamExactFootprint(dq_local[0].path)
                else:
                    dq_set = self.fetch(obs, ["dqmask"],
                                        dqmask_dir or dest)
                    exact = DecamExactFootprint(dq_set.products[0].path)
            if extname is not None:
                if extname not in exact.hdu_index:
                    raise FileNotFoundError(
                        f"no CCD {extname} in exposure {expnum}")
                ccd_name, hdu = extname, exact.hdu_index[extname]
            else:
                ccd_name, hdu = exact.ccd_of(cutout.ra_deg,
                                             cutout.dec_deg)
                if ccd_name is None:
                    raise FileNotFoundError(
                        f"({cutout.ra_deg}, {cutout.dec_deg}) is off "
                        f"every CCD of exposure {expnum}")
            path = dest / f"exp{expnum}_{kind}_{ccd_name}.fits"
            if not path.exists():
                content = self._fetch_hdu(info, hdu, ccd_name, expnum)
                path.write_bytes(content)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            products.append(LocalProduct(
                kind=kind, path=path, checksum=f"sha256:{digest}"))
        if exact is not None:
            exact.close()
        return ProductSet(observation_id=obs.observation_id,
                          products=tuple(products), cutout=cutout)

    def _fetch_hdu(self, info: dict, hdu: int, ccd_name: str,
                   expnum: int) -> bytes:
        """Fetch one CCD HDU with two recovery paths for known archive
        failure modes (~1% of files, recon addendum):

        - the target file's HDU order can differ from its dqmask
          sibling's (seen on a corrupt wtmap): on EXTNAME mismatch the
          per-file order is read from the api/header page and the
          fetch retried at the right index;
        - ``?hdus=`` can 500 for specific files: fall back to the full
          file (md5-verified) and extract the HDU locally.
        """
        url = info["url"]
        resp = self._session.get(url, params={"hdus": f"0,{hdu}"},
                                 timeout=self._timeout)
        if resp.status_code == 404:
            raise FileNotFoundError(url)
        if resp.status_code >= 500:
            return self._fetch_hdu_via_full_file(info, ccd_name, expnum)
        resp.raise_for_status()
        if not resp.content.startswith(b"SIMPLE"):
            raise IOError(f"non-FITS response for {url}: "
                          f"{resp.content[:80]!r}")
        got = self._extname_of(resp.content)
        if got == ccd_name:
            return resp.content
        order = self._hdu_order_from_header_page(info["md5"])
        if order and ccd_name in order:
            idx = order.index(ccd_name) + 1  # order excludes primary
            resp = self._session.get(url, params={"hdus": f"0,{idx}"},
                                     timeout=self._timeout)
            resp.raise_for_status()
            if resp.content.startswith(b"SIMPLE") and \
                    self._extname_of(resp.content) == ccd_name:
                return resp.content
        return self._fetch_hdu_via_full_file(info, ccd_name, expnum)

    def _fetch_hdu_via_full_file(self, info: dict, ccd_name: str,
                                 expnum: int) -> bytes:
        """Fallback full-file fetch. The served bytes are md5-verified
        against the archive metadata; on mismatch the file is still
        accepted IF it parses as FITS and its EXPNUM matches — the
        archive serves repaired versions of some damaged files under
        the original md5sum (recon addendum; 1 wtmap of 214 pilot
        exposures). Unparsable served bytes raise (1 image of 214 is
        corrupt server-side with no recovery)."""
        import io
        import tempfile

        from astropy.io import fits
        resp = self._session.get(info["url"], timeout=self._timeout)
        resp.raise_for_status()
        md5 = hashlib.md5(resp.content).hexdigest()
        md5_mismatch = md5 != info["md5"]
        with tempfile.NamedTemporaryFile(suffix=".fits.fz") as tmp:
            tmp.write(resp.content)
            tmp.flush()
            try:
                hdul = fits.open(tmp.name, output_verify="silentfix")
            except Exception as exc:
                raise IOError(
                    f"unusable archive file {info['url']}: "
                    f"{'md5 mismatch and ' if md5_mismatch else ''}"
                    f"unparsable ({exc})") from exc
            with hdul:
                if md5_mismatch:
                    got_exp = hdul[0].header.get("EXPNUM")
                    if got_exp != expnum:
                        raise IOError(
                            f"md5 mismatch AND EXPNUM {got_exp} != "
                            f"{expnum} for {info['url']}")
                    print(f"  note: served bytes for {info['md5']} "
                          f"(EXPNUM {expnum}) differ from archive md5 "
                          f"(served {md5}); content-verified, accepted",
                          flush=True)
                target = None
                for h in hdul[1:]:
                    if h.header.get("EXTNAME") == ccd_name:
                        target = h
                        break
                if target is None:
                    raise FileNotFoundError(
                        f"no CCD {ccd_name} in {info['url']}")
                buf = io.BytesIO()
                fits.HDUList([fits.PrimaryHDU(header=hdul[0].header),
                              target]).writeto(buf,
                                               output_verify="silentfix")
                return buf.getvalue()

    def _hdu_order_from_header_page(self, md5: str) -> list[str] | None:
        """EXTNAME order from the api/header HTML page; None when the
        page fails (e.g. BADFFILE for corrupt files)."""
        import re
        try:
            resp = self._session.get(
                f"https://astroarchive.noirlab.edu/api/header/{md5}/",
                timeout=self._timeout)
            if resp.status_code != 200:
                return None
            names = re.findall(
                r"EXTNAME</b></th>\s*<td>([A-Za-z0-9]+)</td>", resp.text)
            return names or None
        except Exception:
            return None

    @staticmethod
    def _extname_of(content: bytes) -> str | None:
        import io

        from astropy.io import fits
        with fits.open(io.BytesIO(content)) as hdul:
            return hdul[1].header.get("EXTNAME") if len(hdul) > 1 else None

    @staticmethod
    def _assert_extname(content: bytes, expected: str, url: str) -> None:
        import io

        from astropy.io import fits
        with fits.open(io.BytesIO(content)) as hdul:
            got = hdul[1].header.get("EXTNAME")
        if got != expected:
            raise IOError(f"HDU-order mismatch for {url}: asked for "
                          f"{expected}, got {got}")
