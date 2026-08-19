"""Matched-filter forced photometry on WISE L1b cutouts (plan §4.6).

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
    sigma = fwhm_pix / 2.3548
    y, x = np.mgrid[-radius_pix:radius_pix + 1,
                    -radius_pix:radius_pix + 1]
    k = np.exp(-(x * x + y * y) / (2 * sigma * sigma))
    return k


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
    bgpix = img[good]
    if bgpix.size:
        med = np.median(bgpix)
        mad = 1.4826 * np.median(np.abs(bgpix - med)) or 1.0
        clip = bgpix[np.abs(bgpix - med) < 4 * mad]
        bg = float(np.median(clip)) if clip.size else float(med)
    else:
        bg = 0.0

    fwhm_pix = PSF_FWHM_ARCSEC[band] / PIX_ARCSEC
    kernel = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))

    d = np.where(good, img - bg, 0.0)
    g = good.astype(float)
    s2 = np.where(good, unc * unc, 0.0)

    num = _xcorr(d, kernel)                 # sum p_i d_i
    denom = _xcorr(g, kernel * kernel)      # sum p_i^2 over good pixels
    varnum = _xcorr(s2, kernel * kernel)    # sum p_i^2 s_i^2
    wsum = _xcorr(g, kernel)                # sum p_i over good pixels

    with np.errstate(divide="ignore", invalid="ignore"):
        flux = num / denom
        var = varnum / (denom * denom)
        good_frac = wsum / kernel.sum()
    bad = (denom <= 0) | ~np.isfinite(flux)
    flux[bad] = np.nan
    var[bad] = np.nan

    return FluxMap(flux=flux, var=var, good_frac=good_frac, wcs=wcs,
                   band=band, magzp=magzp, mjd=mjd)
