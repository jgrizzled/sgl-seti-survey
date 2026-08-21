"""IRSA SPHEREx Level-2 spectral-image archive adapter.

Discovery runs against the IRSA TAP service (``spherex.obscore`` joined
to the CAOM ``plane``/``artifact`` views, which carry the product URI
and the published MD5); products are read by HTTP byte-range requests
from the public S3 mirror ``nasa-irsa-spherex`` (IBE serves the same
files and its cutout service is the fallback). Facts verified before
writing this module: surveys/spherex/notes/irsa_recon.md.

A Level-2 file is a multi-extension FITS: IMAGE (MJy/sr, SIP WCS and a
per-pixel wavelength lookup ``WAVE-TAB``), FLAGS (bitmask), VARIANCE
(MJy^2/sr^2), ZODI (model), PSF (11 x 11 detector zones, 10x
oversampled) and WCS-WAVE (9 x 9 wavelength/bandwidth table). ``fetch``
assembles a *slim cutout* of all of these around a sky position: the
sub-arrays of the four image planes with CRPIX shifted, the PSF plane
of the zone containing the cutout centre, and the wavelength table.
"""

from __future__ import annotations

import csv
import hashlib
import io
import re
import warnings
from pathlib import Path
from typing import Iterator, Sequence

import numpy as np
import requests

from sglsurvey.adapters.base import (ConeRegion, CutoutSpec, LocalProduct,
                                     MjdRange, ProductSet)
from sglsurvey.records import Observation
from sglsurvey.snapshots import SnapshotStore

TAP_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"
IRSA_ROOT = "https://irsa.ipac.caltech.edu/"
S3_ROOT = "https://nasa-irsa-spherex.s3.us-east-1.amazonaws.com/"
COLLECTIONS = ("spherex_qr2", "spherex_qr2_deep")

NX = NY = 2040
PIX_ARCSEC = 6.15
#: Half-diagonal of the 3.5 deg detector field, for the centre-distance
#: discovery superset (ADQL INTERSECTS on s_region is not supported).
HALF_DIAGONAL_DEG = 2.6

#: FLAGS bits (extension header MP_* cards, verified 2026-08-20).
FLAG_BITS = {
    "TRANSIENT": 0, "OVERFLOW": 1, "SUR_ERROR": 2, "PHANTOM": 4,
    "REFERENCE": 5, "NONFUNC": 6, "DICHROIC": 7, "MISSING_DATA": 9,
    "HOT": 10, "COLD": 11, "FULLSAMPLE": 12, "PHANMISS": 14,
    "NONLINEAR": 15, "PERSIST": 17, "OUTLIER": 19, "SOURCE": 21,
}
#: Every bit except FULLSAMPLE (12, informational) and SOURCE (21, a
#: known source is present — exactly what a search must not exclude).
FATAL_FLAG_TEMPLATE = sum(1 << b for n, b in FLAG_BITS.items()
                          if n not in ("FULLSAMPLE", "SOURCE"))

OBSCORE_COLS = [
    "o.obs_id", "o.s_ra", "o.s_dec", "o.energy_bandpassname", "o.em_min",
    "o.em_max", "o.t_min", "o.t_max", "o.t_exptime", "o.obs_collection",
    "o.obs_release_date", "o.data_rights", "o.s_region", "o.s_resolution",
    "o.s_xel1", "o.s_xel2", "o.s_pixel_scale", "o.access_url",
    "o.obs_publisher_did", "p.planeid", "p.provenance_version",
    "p.quality_flag", "a.uri", "a.contentchecksum", "a.contentlength",
]


def _f(row: dict, key: str) -> float | None:
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return None


def parse_region(s_region: str) -> tuple[tuple[float, float], ...]:
    """``POLYGON ICRS ra dec ra dec ...`` -> vertices (closing vertex
    dropped if repeated)."""
    nums = [float(x) for x in re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?",
                                         s_region.split("ICRS", 1)[-1])]
    pts = [(nums[i], nums[i + 1]) for i in range(0, len(nums) - 1, 2)]
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    return tuple(pts)


def _gnomonic(ra_deg, dec_deg, ra0_deg, dec0_deg):
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    ra0, dec0 = np.deg2rad(ra0_deg), np.deg2rad(dec0_deg)
    cosc = (np.sin(dec0) * np.sin(dec)
            + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0))
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.cos(dec) * np.sin(ra - ra0) / cosc
        y = (np.cos(dec0) * np.sin(dec)
             - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    x = np.where(cosc > 0, x, np.nan)
    y = np.where(cosc > 0, y, np.nan)
    return np.rad2deg(x), np.rad2deg(y)


class SpherexNominalFootprint:
    """Containment test from the ObsCore ``s_region`` polygon projected
    gnomonically about the field centre — coarse stage only (the true
    SIP footprint differs by up to a few arcsec at the edges)."""

    def __init__(self, obs: Observation, pad_arcsec: float = 0.0):
        self._ra0 = float(obs.wcs["s_ra"])
        self._dec0 = float(obs.wcs["s_dec"])
        vx, vy = _gnomonic([p[0] for p in obs.corners_icrs_deg],
                           [p[1] for p in obs.corners_icrs_deg],
                           self._ra0, self._dec0)
        self._vx, self._vy = np.asarray(vx), np.asarray(vy)
        self._pad_deg = pad_arcsec / 3600.0

    def _inside_and_dist(self, x, y):
        vx, vy = self._vx, self._vy
        n = len(vx)
        inside = np.zeros(x.shape, dtype=bool)
        dmin = np.full(x.shape, np.inf)
        for i in range(n):
            ax, ay = vx[i], vy[i]
            bx, by = vx[(i + 1) % n], vy[(i + 1) % n]
            cond = ((ay > y) != (by > y))
            with np.errstate(divide="ignore", invalid="ignore"):
                xint = ax + (y - ay) * (bx - ax) / (by - ay)
            inside ^= cond & (x < xint)
            dx, dy = bx - ax, by - ay
            L2 = dx * dx + dy * dy
            t = np.clip(((x - ax) * dx + (y - ay) * dy) / L2, 0.0, 1.0)
            d = np.hypot(x - (ax + t * dx), y - (ay + t * dy))
            dmin = np.minimum(dmin, d)
        return inside, dmin

    def contains_any(self, radec_deg) -> bool:
        radec_deg = np.asarray(radec_deg, dtype=float)
        x, y = _gnomonic(radec_deg[:, 0], radec_deg[:, 1],
                         self._ra0, self._dec0)
        ok = np.isfinite(x) & np.isfinite(y)
        if not ok.any():
            return False
        inside, dmin = self._inside_and_dist(x[ok], y[ok])
        return bool(np.any(inside | (dmin <= self._pad_deg)))

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        return self.contains_any(np.array([[ra_deg, dec_deg]]))


class SpherexExactFootprint:
    """Exact usable-pixel test on a slim cutout: SIP WCS from the FLAGS
    extension plus the FATAL_FLAG_TEMPLATE bitplanes. Pixels outside the
    cutout are out-of-bounds (not unusable) — size cutouts to enclose
    the locus and treat out-of-bounds as "not evaluated"."""

    FATAL_MASK = FATAL_FLAG_TEMPLATE

    def __init__(self, path: Path):
        from astropy.io import fits
        from astropy.utils.exceptions import AstropyWarning
        from astropy.wcs import WCS

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            with fits.open(path) as hdul:
                hdu = hdul["FLAGS"]
                self._mask = np.asarray(hdu.data, dtype=np.int64)
                self._wcs = WCS(hdu.header)
                self.header = dict(hdul["IMAGE"].header)
        self._ny, self._nx = self._mask.shape

    def status(self, radec_deg):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pix = self._wcs.all_world2pix(np.asarray(radec_deg, float), 0,
                                          quiet=True)
        with np.errstate(invalid="ignore"):
            x = np.rint(pix[:, 0]).astype(int)
            y = np.rint(pix[:, 1]).astype(int)
            inb = (x >= 0) & (x < self._nx) & (y >= 0) & (y < self._ny)
        usable = np.zeros(len(radec_deg), dtype=bool)
        if inb.any():
            usable[inb] = (self._mask[y[inb], x[inb]] & self.FATAL_MASK) == 0
        return inb, usable

    def contains(self, ra_deg: float, dec_deg: float) -> bool:
        _, usable = self.status(np.array([[ra_deg, dec_deg]]))
        return bool(usable[0])


class HttpRangeFile:
    """Minimal seekable read-only file over HTTP byte-range requests with
    a block cache, good enough for astropy's lazy HDU loading and
    ``ImageHDU.section`` slicing."""

    def __init__(self, url: str, session: requests.Session,
                 block_size: int = 1 << 20, timeout_s: float = 120.0,
                 size: int | None = None):
        self.url, self._s, self._bs, self._to = url, session, block_size, timeout_s
        self._pos = 0
        self._blocks: dict[int, bytes] = {}
        self.bytes_fetched = 0
        self.requests = 0
        self.etag = None
        if size is None:
            r = session.head(url, timeout=timeout_s)
            r.raise_for_status()
            size = int(r.headers["Content-Length"])
            self.etag = r.headers.get("ETag")
        self._size = size
        self.mode = "rb"
        self.closed = False

    # -- file protocol ---------------------------------------------------
    def seekable(self):
        return True

    def readable(self):
        return True

    def writable(self):
        return False

    def tell(self):
        return self._pos

    def seek(self, off, whence=0):
        if whence == 0:
            self._pos = off
        elif whence == 1:
            self._pos += off
        else:
            self._pos = self._size + off
        return self._pos

    def _block(self, i: int) -> bytes:
        b = self._blocks.get(i)
        if b is None:
            lo = i * self._bs
            hi = min(self._size, lo + self._bs) - 1
            for attempt in (1, 2, 3):
                try:
                    r = self._s.get(self.url, headers={"Range": f"bytes={lo}-{hi}"},
                                    timeout=self._to)
                    if r.status_code not in (200, 206):
                        r.raise_for_status()
                    b = r.content
                    break
                except Exception:
                    if attempt == 3:
                        raise
            self.requests += 1
            self.bytes_fetched += len(b)
            self._blocks[i] = b
        return b

    def read(self, n=-1):
        if n is None or n < 0:
            n = self._size - self._pos
        n = max(0, min(n, self._size - self._pos))
        out = bytearray()
        pos = self._pos
        while len(out) < n:
            i, off = divmod(pos, self._bs)
            chunk = self._block(i)[off:off + (n - len(out))]
            if not chunk:
                break
            out += chunk
            pos += len(chunk)
        self._pos = pos
        return bytes(out)

    def readinto(self, buf):
        data = self.read(len(buf))
        buf[:len(data)] = data
        return len(data)

    def close(self):
        self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


def psf_zone_index(header, x_pix: float, y_pix: float) -> int:
    """0-based plane index of the PSF zone containing detector pixel
    (x, y) (0-based, full-frame coordinates)."""
    best, bestd = 0, np.inf
    n = int(header["NAXIS3"])
    for k in range(1, n + 1):
        d = np.hypot(header[f"XCTR_{k}"] - x_pix, header[f"YCTR_{k}"] - y_pix)
        if d < bestd:
            best, bestd = k - 1, d
    return best


class SpherexL2Adapter:
    archive_id = "irsa-spherex"
    collection = "spherex_qr2"
    #: FITS extensions kept in a slim cutout (PSF reduced to one plane).
    IMAGE_EXTS = ("IMAGE", "FLAGS", "VARIANCE", "ZODI")

    def __init__(self, session: requests.Session | None = None,
                 timeout_s: float = 900.0, source: str = "s3"):
        self._session = session or requests.Session()
        self._timeout = timeout_s
        self._source = source

    # -- discovery ---------------------------------------------------------
    def discovery_query(self, region: ConeRegion, time_range: MjdRange) -> str:
        r = region.radius_deg + HALF_DIAGONAL_DEG
        return (
            f"SELECT {', '.join(OBSCORE_COLS)} FROM spherex.obscore o "
            "JOIN spherex.plane p ON o.obs_publisher_did = p.obs_publisher_did "
            "JOIN spherex.artifact a ON a.planeid = p.planeid "
            f"WHERE CONTAINS(POINT('ICRS',o.s_ra,o.s_dec),"
            f"CIRCLE('ICRS',{region.ra_deg:.6f},{region.dec_deg:.6f},{r:.4f}))=1 "
            f"AND o.t_max >= {time_range.start_mjd_utc:.5f} "
            f"AND o.t_min <= {time_range.stop_mjd_utc:.5f} "
            "AND o.dataproduct_type = 'image' AND o.calib_level = 2 "
            "AND o.data_rights = 'public' AND a.producttype = 'science'"
        )

    def discover(self, region: ConeRegion, time_range: MjdRange,
                 store: SnapshotStore) -> Iterator[Observation]:
        from datetime import datetime, timezone

        query = self.discovery_query(region, time_range)
        request_utc = datetime.now(timezone.utc).isoformat()
        resp = self._session.get(TAP_SYNC, params={"QUERY": query,
                                                   "FORMAT": "CSV"},
                                 timeout=self._timeout)
        resp.raise_for_status()
        if resp.text.lstrip().startswith("<"):
            raise RuntimeError("TAP error: " + resp.text[:500])
        rows = list(csv.DictReader(io.StringIO(resp.text)))
        snapshot = store.store(service_url=TAP_SYNC, query=query,
                               request_utc=request_utc,
                               response_bytes=resp.content,
                               row_count=len(rows),
                               http_status=resp.status_code)
        for row in rows:
            yield self._observation(row, snapshot.snapshot_id, request_utc)

    def _observation(self, row: dict, snapshot_id: str,
                     discovered_utc: str) -> Observation:
        band = row["energy_bandpassname"].replace("SPHEREx-", "")
        detector = int(band[1:])
        t0, t1 = float(row["t_min"]), float(row["t_max"])
        corners = parse_region(row["s_region"])
        uri = row["uri"]
        rel = uri.split("ibe/data/spherex/", 1)[-1]
        em0, em1 = _f(row, "em_min"), _f(row, "em_max")
        return Observation.build(
            archive_id=self.archive_id, collection=row["obs_collection"],
            release=f"l2b-{row.get('provenance_version', '')}",
            native_key={"obs_id": row["obs_id"], "detector": detector},
            band=band,
            wavelength_um=((em0 + em1) / 2.0 * 1e6
                           if em0 is not None and em1 is not None else None),
            t_start_mjd_utc=t0, t_mid_mjd_utc=0.5 * (t0 + t1),
            t_stop_mjd_utc=t1, exptime_s=_f(row, "t_exptime") or (t1 - t0) * 86400,
            corners_icrs_deg=corners,
            wcs={"naxis1": int(_f(row, "s_xel1") or NX),
                 "naxis2": int(_f(row, "s_xel2") or NY),
                 "s_ra": float(row["s_ra"]), "s_dec": float(row["s_dec"]),
                 "ctype1": "RA---TAN-SIP", "ctype2": "DEC--TAN-SIP",
                 "pixel_scale_arcsec": _f(row, "s_pixel_scale")},
            quality_flags={"quality_flag": row.get("quality_flag") or None,
                           "s_resolution_arcsec": _f(row, "s_resolution")},
            products={"l2": {"url": IRSA_ROOT + uri,
                             "s3_url": S3_ROOT + rel,
                             "datalink": row.get("access_url"),
                             "checksum": row.get("contentchecksum"),
                             "content_length": int(_f(row, "contentlength") or 0)}},
            extra={"em_min_um": em0 * 1e6 if em0 else None,
                   "em_max_um": em1 * 1e6 if em1 else None,
                   "planeid": row.get("planeid"),
                   "obs_release_date": row.get("obs_release_date"),
                   "obs_publisher_did": row.get("obs_publisher_did")},
            snapshot_id=snapshot_id, discovered_utc=discovered_utc,
        )

    # -- footprints --------------------------------------------------------
    def nominal_footprint(self, obs: Observation,
                          pad_arcsec: float = 0.0) -> SpherexNominalFootprint:
        return SpherexNominalFootprint(obs, pad_arcsec)

    def exact_footprint(self, obs: Observation,
                        products: ProductSet) -> SpherexExactFootprint:
        p = next(p for p in products.products if p.kind == "l2cut")
        return SpherexExactFootprint(p.path)

    # -- fetch -------------------------------------------------------------
    @staticmethod
    def cutout_name(obs: Observation, cutout: CutoutSpec) -> str:
        return (f"cut{cutout.size_pix}-{obs.native_key['obs_id']}"
                f"D{obs.native_key['detector']}.fits")

    def fetch(self, obs: Observation, kinds: Sequence[str], dest: Path,
              *, cutout: CutoutSpec | None = None) -> ProductSet:
        """Assemble a slim cutout (kind ``l2cut``) by byte-range reads
        from S3 (or IBE). Full-file downloads are deliberately not
        offered (71 MB each). The recorded checksum is the sha256 of the
        assembled file; the archive-published MD5 of the parent file and
        the S3 ETag are written into the primary header."""
        if cutout is None:
            raise ValueError("SPHEREx adapter serves cutouts only")
        if list(kinds) != ["l2cut"]:
            raise ValueError(f"unknown product kinds {kinds}")
        dest.mkdir(parents=True, exist_ok=True)
        path = dest / self.cutout_name(obs, cutout)
        if not path.exists():
            self._assemble_cutout(obs, cutout, path)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        return ProductSet(observation_id=obs.observation_id,
                          products=(LocalProduct(kind="l2cut", path=path,
                                                 checksum=f"sha256:{digest}"),),
                          cutout=cutout)

    def _assemble_cutout(self, obs: Observation, cutout: CutoutSpec,
                         path: Path) -> None:
        from astropy.io import fits
        from astropy.utils.exceptions import AstropyWarning
        from astropy.wcs import WCS

        prod = obs.products["l2"]
        url = prod["s3_url"] if self._source == "s3" else prod["url"]
        size = prod.get("content_length") or None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AstropyWarning)
            fobj = HttpRangeFile(url, self._session, timeout_s=self._timeout,
                                 size=size)
            with fits.open(fobj, lazy_load_hdus=True, memmap=False) as hdul:
                img_hdr = hdul["IMAGE"].header
                wcs = WCS(img_hdr)
                xc, yc = wcs.all_world2pix([[cutout.ra_deg, cutout.dec_deg]], 0)[0]
                if not (np.isfinite(xc) and np.isfinite(yc)):
                    raise FileNotFoundError("cutout centre off detector")
                half = cutout.size_pix // 2
                x0 = int(np.clip(int(round(xc)) - half, 0, NX))
                y0 = int(np.clip(int(round(yc)) - half, 0, NY))
                x1 = int(np.clip(x0 + cutout.size_pix, 0, NX))
                y1 = int(np.clip(y0 + cutout.size_pix, 0, NY))
                if x1 - x0 < 2 or y1 - y0 < 2:
                    raise FileNotFoundError("cutout does not intersect image")
                out = [fits.PrimaryHDU()]
                ph = out[0].header
                for card in hdul[0].header.cards:
                    if card.keyword not in ("SIMPLE", "BITPIX", "NAXIS",
                                            "EXTEND"):
                        ph.append(card)
                ph["SRCURL"] = (url, "parent file")
                ph["SRCMD5"] = (prod.get("checksum") or "", "published checksum of parent")
                ph["SRCETAG"] = (fobj.etag or "", "S3 ETag of parent")
                ph["CUTX0"] = (x0, "0-based parent column of cutout origin")
                ph["CUTY0"] = (y0, "0-based parent row of cutout origin")
                ph["CUTRA"] = (cutout.ra_deg, "requested centre RA")
                ph["CUTDEC"] = (cutout.dec_deg, "requested centre Dec")
                for name in self.IMAGE_EXTS:
                    h = hdul[name]
                    data = h.section[y0:y1, x0:x1]
                    hdr = h.header.copy()
                    for suffix in ("", "A", "W"):
                        if f"CRPIX1{suffix}" in hdr:
                            hdr[f"CRPIX1{suffix}"] -= x0
                            hdr[f"CRPIX2{suffix}"] -= y0
                    out.append(fits.ImageHDU(data=np.ascontiguousarray(data),
                                             header=hdr, name=name))
                psf = hdul["PSF"]
                k = psf_zone_index(psf.header, xc, yc)
                plane = np.ascontiguousarray(psf.section[k])
                phdr = psf.header.copy()
                phdr["PSFZONE"] = (k, "0-based plane index kept from parent")
                phdr["PSFXC"] = (float(psf.header[f"XCTR_{k + 1}"]), "zone centre x")
                phdr["PSFYC"] = (float(psf.header[f"YCTR_{k + 1}"]), "zone centre y")
                for key in list(phdr):
                    if key.startswith(("XCTR_", "YCTR_", "XWID_", "YWID_")):
                        del phdr[key]
                phdr["NAXIS"] = 2
                for key in ("NAXIS3",):
                    if key in phdr:
                        del phdr[key]
                out.append(fits.ImageHDU(data=plane, header=phdr, name="PSF"))
                wave = hdul["WCS-WAVE"]
                out.append(fits.BinTableHDU(data=wave.data, header=wave.header,
                                            name="WCS-WAVE"))
                ph["NREQ"] = (fobj.requests, "range requests used")
                ph["NBYTES"] = (fobj.bytes_fetched, "bytes fetched")
                tmp = path.with_suffix(".tmp")
                fits.HDUList(out).writeto(tmp, overwrite=True)
                tmp.rename(path)


def wavelength_at(wave_table, img_header, x_pix: float, y_pix: float,
                  ) -> tuple[float, float]:
    """(wavelength_um, bandwidth_um) at 0-based cutout pixel (x, y) by
    bilinear interpolation of the 9 x 9 WCS-WAVE table. The table is
    indexed by the image's ``W`` coordinate system (raw 1-based FITS
    pixel coordinates of the parent frame): coord = CRVALiW +
    (p_1based - CRPIXiW) * CDELTiW, which stays valid on cutouts whose
    CRPIXiW were shifted with the data."""
    row = wave_table[0]
    xs = np.asarray(row["X"], dtype=float)
    ys = np.asarray(row["Y"], dtype=float)
    vals = np.asarray(row["VALUES"], dtype=float)  # (ny, nx, 2)
    h = img_header
    xq = h.get("CRVAL1W", 1.0) + (x_pix + 1.0 - h.get("CRPIX1W", 0.0)) * h.get("CDELT1W", 1.0)
    yq = h.get("CRVAL2W", 1.0) + (y_pix + 1.0 - h.get("CRPIX2W", 0.0)) * h.get("CDELT2W", 1.0)
    xq = float(np.clip(xq, xs[0], xs[-1]))
    yq = float(np.clip(yq, ys[0], ys[-1]))
    i = int(np.clip(np.searchsorted(xs, xq) - 1, 0, len(xs) - 2))
    j = int(np.clip(np.searchsorted(ys, yq) - 1, 0, len(ys) - 2))
    fx = (xq - xs[i]) / (xs[i + 1] - xs[i])
    fy = (yq - ys[j]) / (ys[j + 1] - ys[j])
    out = []
    for c in (0, 1):
        v = (vals[j, i, c] * (1 - fx) * (1 - fy) + vals[j, i + 1, c] * fx * (1 - fy)
             + vals[j + 1, i, c] * (1 - fx) * fy + vals[j + 1, i + 1, c] * fx * fy)
        out.append(float(v))
    return out[0], out[1]
