"""Matched-filter forced photometry on image cutouts (plan §4.6): WISE
L1b int/unc/msk triplets and ZTF sci/diff/msk cutouts.

Per cutout we build a PSF-matched flux map and its variance map once;
any trajectory position is then a cheap bilinear lookup. The estimator
is the standard PSF-weighted least-squares point-source amplitude with
a Gaussian PSF approximation:

    f_hat(x0) = sum_i p_i (d_i - b) / sum_i p_i^2
    var(x0)   = sum_i p_i^2 s_i^2 / (sum_i p_i^2)^2

evaluated in closed form via FFT cross-correlation, with fatally masked
pixels excluded (weight renormalized per position) and a sigma-clipped
median background. Gaussian-PSF and constant-background are declared
pilot approximations; they cost sensitivity, never validity, and the
injection stage measures the end-to-end throughput.
"""

from __future__ import annotations

import gzip
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PSF_FWHM_ARCSEC = {"W1": 6.1, "W2": 6.4, "W3": 6.5, "W4": 12.0}
PIX_ARCSEC = 2.75


def _read_fits(path: Path):
    from astropy.io import fits
    from astropy.utils.exceptions import AstropyWarning
    from astropy.wcs import WCS

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rb") as fh:
            with fits.open(fh) as hdul:
                data = np.asarray(hdul[0].data, dtype=np.float64)
                header = hdul[0].header
                wcs = WCS(header)
    return data, header, wcs


def _gaussian_kernel(fwhm_pix: float, radius_pix: int) -> np.ndarray:
    """Unit-sum Gaussian PSF model, so the least-squares amplitude
    ``sum p d / sum p^2`` is the source's TOTAL flux in image units.

    History: before 2026-08-20 the kernel had unit PEAK, so fluxes (and
    every magnitude derived from them) were the peak amplitude, fainter
    than the total by 2.5 log10(2 pi sigma_pix^2) — caught by the ZTF
    asteroid positive control. S/N, thresholds and recovery fractions
    are invariant to this scale; only magnitude labels changed.
    """
    sigma = fwhm_pix / 2.3548
    y, x = np.mgrid[-radius_pix:radius_pix + 1,
                    -radius_pix:radius_pix + 1]
    k = np.exp(-(x * x + y * y) / (2 * sigma * sigma))
    return k / k.sum()


def _xcorr(img: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Same-size cross-correlation via FFT (kernel is symmetric)."""
    from numpy.fft import irfft2, rfft2

    ny, nx = img.shape
    ky, kx = kernel.shape
    py, px = ny + ky - 1, nx + kx - 1
    F = rfft2(img, (py, px)) * rfft2(kernel, (py, px))
    full = irfft2(F, (py, px))
    oy, ox = (ky - 1) // 2, (kx - 1) // 2
    return full[oy:oy + ny, ox:ox + nx]


@dataclass
class FluxMap:
    """PSF-matched flux/variance maps for one cutout, with WCS lookup."""

    flux: np.ndarray        # matched-filter amplitude (DN)
    var: np.ndarray         # variance of the amplitude (DN^2)
    good_frac: np.ndarray   # fraction of PSF weight on unmasked pixels
    wcs: object
    band: str
    magzp: float | None
    mjd: float
    aux: np.ndarray | None = None  # optional per-pixel auxiliary map

    def sample_aux(self, ra_deg, dec_deg):
        """Bilinear sample of ``aux`` (NaN outside / if absent)."""
        ra = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        if self.aux is None:
            return np.full(len(ra), np.nan)
        dec = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        pix = self.wcs.wcs_world2pix(np.stack([ra, dec], axis=1), 0)
        x, y = pix[:, 0], pix[:, 1]
        ny, nx = self.aux.shape
        out = np.full(len(ra), np.nan)
        ok = (x >= 0) & (x <= nx - 1.001) & (y >= 0) & (y <= ny - 1.001) \
            & np.isfinite(x) & np.isfinite(y)
        if ok.any():
            x0 = np.floor(x[ok]).astype(int)
            y0 = np.floor(y[ok]).astype(int)
            fx, fy = x[ok] - x0, y[ok] - y0
            a = self.aux
            out[ok] = (a[y0, x0] * (1 - fx) * (1 - fy)
                       + a[y0, x0 + 1] * fx * (1 - fy)
                       + a[y0 + 1, x0] * (1 - fx) * fy
                       + a[y0 + 1, x0 + 1] * fx * fy)
        return out

    def sample(self, ra_deg, dec_deg):
        """Bilinear sample (flux, var, good_frac) at sky positions.

        Returns NaN outside the cutout. Inputs are arrays or scalars.
        """
        ra = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        dec = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        pix = self.wcs.wcs_world2pix(np.stack([ra, dec], axis=1), 0)
        x, y = pix[:, 0], pix[:, 1]
        ny, nx = self.flux.shape
        out_f = np.full(len(ra), np.nan)
        out_v = np.full(len(ra), np.nan)
        out_g = np.full(len(ra), np.nan)
        ok = (x >= 0) & (x <= nx - 1.001) & (y >= 0) & (y <= ny - 1.001) \
            & np.isfinite(x) & np.isfinite(y)
        if ok.any():
            x0 = np.floor(x[ok]).astype(int)
            y0 = np.floor(y[ok]).astype(int)
            fx, fy = x[ok] - x0, y[ok] - y0

            def bil(a):
                return (a[y0, x0] * (1 - fx) * (1 - fy)
                        + a[y0, x0 + 1] * fx * (1 - fy)
                        + a[y0 + 1, x0] * (1 - fx) * fy
                        + a[y0 + 1, x0 + 1] * fx * fy)

            out_f[ok] = bil(self.flux)
            out_v[ok] = bil(self.var)
            out_g[ok] = bil(self.good_frac)
        return out_f, out_v, out_g


def build_flux_map(int_path: Path, unc_path: Path, msk_full_path: Path,
                   band: str, mjd: float, fatal_mask: int,
                   magzp: float | None = None) -> FluxMap:
    """Build the matched-filter maps for one cutout.

    The full-frame -msk is aligned to the cutout by the integer CRPIX
    offset (IBE cutouts preserve the native pixel grid).
    """
    img, hdr, wcs = _read_fits(int_path)
    unc, _, _ = _read_fits(unc_path)
    msk, mhdr, _ = _read_fits(msk_full_path)

    dx = int(round(mhdr["CRPIX1"] - hdr["CRPIX1"]))
    dy = int(round(mhdr["CRPIX2"] - hdr["CRPIX2"]))
    ny, nx = img.shape
    mcut = np.zeros((ny, nx), dtype=np.int64)
    y0, y1 = max(0, -dy), min(ny, msk.shape[0] - dy)
    x0, x1 = max(0, -dx), min(nx, msk.shape[1] - dx)
    mcut[y0:y1, x0:x1] = msk[y0 + dy:y1 + dy, x0 + dx:x1 + dx]

    good = (np.isfinite(img) & np.isfinite(unc) & (unc > 0)
            & ((mcut & fatal_mask) == 0))
    fwhm_pix = PSF_FWHM_ARCSEC[band] / PIX_ARCSEC
    flux, var, good_frac = matched_filter(img, unc * unc, good, fwhm_pix)
    return FluxMap(flux=flux, var=var, good_frac=good_frac, wcs=wcs,
                   band=band, magzp=magzp, mjd=mjd)


def matched_filter(img: np.ndarray, var_pix: np.ndarray, good: np.ndarray,
                   fwhm_pix: float, subtract_background: bool = True):
    """Core estimator: returns (flux, var, good_frac) maps."""
    if subtract_background:
        bgpix = img[good]
        if bgpix.size:
            med = np.median(bgpix)
            mad = 1.4826 * np.median(np.abs(bgpix - med)) or 1.0
            clip = bgpix[np.abs(bgpix - med) < 4 * mad]
            bg = float(np.median(clip)) if clip.size else float(med)
        else:
            bg = 0.0
    else:
        bg = 0.0
    kernel = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
    d = np.where(good, img - bg, 0.0)
    g = good.astype(float)
    s2 = np.where(good, var_pix, 0.0)
    num = _xcorr(d, kernel)
    denom = _xcorr(g, kernel * kernel)
    varnum = _xcorr(s2, kernel * kernel)
    wsum = _xcorr(g, kernel)
    with np.errstate(divide="ignore", invalid="ignore"):
        flux = num / denom
        var = varnum / (denom * denom)
        good_frac = wsum / kernel.sum()
    bad = (denom <= 0) | ~np.isfinite(flux)
    flux[bad] = np.nan
    var[bad] = np.nan
    return flux, var, good_frac


def _read_fits_any(path: Path):
    """Read the first image HDU (handles fpacked 2-HDU ZTF diff cutouts)."""
    from astropy.io import fits
    from astropy.utils.exceptions import AstropyWarning
    from astropy.wcs import WCS

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        with fits.open(path) as hdul:
            hdu = next(h for h in hdul if h.data is not None)
            data = np.asarray(hdu.data, dtype=np.float64)
            header = hdu.header
            wcs = WCS(header)
    return data, header, wcs


ZTF_PIX_ARCSEC = 1.012


def build_flux_map_ztf(sci_path: Path, diff_path: Path | None,
                       msk_path: Path, fatal_mask: int, band: str,
                       mjd: float) -> FluxMap:
    """ZTF variant. The search image is a pixelwise hybrid: the
    PSF-matched difference image wherever the reference image exists
    (static field removed), and the sky-subtracted science image in the
    strip outside the reference footprint (constant fill value in the
    difference product, e.g. -7440). ``FluxMap.aux`` holds the fraction
    of matched-filter weight drawn from difference-image pixels at each
    position, so downstream stages can separate the two regimes.

    Per-pixel variance = robust background variance of the search image
    (computed separately for the two regimes) + Poisson term from the
    science image (|sci - sky| / GAIN, DN^2). PSF FWHM from the science
    header SEEING; MAGZP from the header. All cutouts share one pixel
    grid (same IBE cutout request), checked via CRPIX.
    """
    sci, hdr, wcs = _read_fits_any(sci_path)
    msk, mhdr, _ = _read_fits_any(msk_path)
    mcut = np.asarray(msk, dtype=np.int64)
    unmasked = np.isfinite(sci) & ((mcut & fatal_mask) == 0)
    gain = float(hdr.get("GAIN", 6.2))
    scipix = sci[unmasked]
    sky = float(np.median(scipix)) if scipix.size else 0.0

    is_diff = np.zeros(sci.shape, dtype=bool)
    img = sci - sky
    if diff_path is not None:
        diff, dhdr, _ = _read_fits_any(diff_path)
        for other in (mhdr, dhdr):
            if (abs(other["CRPIX1"] - hdr["CRPIX1"]) > 0.01
                    or abs(other["CRPIX2"] - hdr["CRPIX2"]) > 0.01
                    or other["NAXIS1"] != hdr["NAXIS1"]):
                raise ValueError("cutout grids differ")
        is_diff = np.isfinite(diff)
        vals, counts = np.unique(np.round(diff[is_diff]), return_counts=True)
        if counts.size and counts.max() > 0.02 * diff.size:
            fill = vals[np.argmax(counts)]
            if abs(fill) > 100:
                is_diff &= np.abs(diff - fill) > 2.0
        img = np.where(is_diff, diff, img)
    good = unmasked & np.isfinite(img)

    def robust_var(pix):
        if pix.size < 50:
            return None
        med = np.median(pix)
        mad = 1.4826 * np.median(np.abs(pix - med))
        return float(mad * mad) if mad > 0 else float(np.var(pix))

    v_diff = robust_var(img[good & is_diff])
    v_sci = robust_var(img[good & ~is_diff])
    v_any = v_diff if v_diff is not None else (v_sci or 1.0)
    bgvar = np.where(is_diff, v_diff if v_diff is not None else v_any,
                     v_sci if v_sci is not None else v_any)
    poisson = np.clip(sci - sky, 0.0, None) / gain
    var_pix = bgvar + poisson
    fwhm_arcsec = float(hdr.get("SEEING", 2.0)) or 2.0
    fwhm_pix = fwhm_arcsec / ZTF_PIX_ARCSEC
    # background already removed in both regimes
    flux, var, good_frac = matched_filter(img, var_pix, good, fwhm_pix,
                                          subtract_background=False)
    kernel = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
    with np.errstate(divide="ignore", invalid="ignore"):
        dfrac = _xcorr((good & is_diff).astype(float), kernel) / np.maximum(
            _xcorr(good.astype(float), kernel), 1e-9)
    return FluxMap(flux=flux, var=var, good_frac=good_frac, wcs=wcs,
                   band=band, magzp=(float(hdr["MAGZP"])
                                     if "MAGZP" in hdr else None),
                   mjd=mjd, aux=np.clip(dfrac, 0.0, 1.0))
