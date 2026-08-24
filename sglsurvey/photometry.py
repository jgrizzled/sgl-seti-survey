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
    sip: bool = False  # use the iterative SIP inverse (SPHEREx)
    upsample: int = 1  # maps evaluated at 1/upsample-pixel phases
    # v2 (wise v2_plan §4.5 stamp-response trick): the estimator's
    # denominator map sum_i p_i^2 g_i, the (flipped, unit-sum) kernel,
    # the usable-pixel mask and the background level, so that the
    # response to an added source can be computed locally and exactly.
    denom: np.ndarray | None = None
    kernel: np.ndarray | None = None
    good: np.ndarray | None = None
    bg: float | None = None

    def world2pix(self, ra, dec):
        """Sky -> map pixel coordinates (image pixels x ``upsample``)."""
        xy = np.stack([ra, dec], axis=1)
        if self.sip:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pix = self.wcs.all_world2pix(xy, 0, quiet=True)
        else:
            pix = self.wcs.wcs_world2pix(xy, 0)
        return pix * self.upsample

    def sample_aux(self, ra_deg, dec_deg):
        """Bilinear sample of ``aux`` (NaN outside / if absent)."""
        ra = np.atleast_1d(np.asarray(ra_deg, dtype=float))
        if self.aux is None:
            return np.full(len(ra), np.nan)
        dec = np.atleast_1d(np.asarray(dec_deg, dtype=float))
        pix = self.world2pix(ra, dec)
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
        pix = self.world2pix(ra, dec)
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


@dataclass
class WiseCutout:
    """Calibrated WISE L1b cutout arrays, before any estimator runs —
    the seam at which image-level injection happens (v2 plan §4.1)."""

    img: np.ndarray       # DN
    var_pix: np.ndarray   # DN^2
    good: np.ndarray      # usable pixels (finite, unc > 0, no fatal bit)
    header: object
    wcs: object
    mask_bits: np.ndarray  # aligned -msk bitplanes
    frame_origin: tuple[int, int]  # (x, y) full-frame pixel of cutout (0, 0)


def read_wise_cutout(int_path: Path, unc_path: Path, msk_full_path: Path,
                     fatal_mask: int) -> WiseCutout:
    """Read an int/unc cutout pair and align the full-frame -msk to it
    by the integer CRPIX offset (IBE cutouts preserve the native grid).
    ``frame_origin`` is the full-frame pixel position of cutout pixel
    (0, 0), so injection code can pick the focal-plane PRF element."""
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
    return WiseCutout(img=img, var_pix=unc * unc, good=good, header=hdr,
                      wcs=wcs, mask_bits=mcut, frame_origin=(dx, dy))


def flux_map_from_arrays(img: np.ndarray, var_pix: np.ndarray,
                         good: np.ndarray, wcs, band: str, mjd: float,
                         magzp: float | None = None,
                         kernel: np.ndarray | None = None,
                         keep_inputs: bool = False,
                         pix_arcsec: float = PIX_ARCSEC) -> FluxMap:
    """Matched-filter maps from already-read arrays (second half of the
    seam). ``kernel`` overrides the per-band Gaussian. With
    ``keep_inputs`` the map also carries ``denom``, ``kernel``, ``good``
    and ``bg`` for local re-evaluation after an injection.

    ``pix_arcsec`` sets the Gaussian kernel width in pixels. v1 used the
    W1-W3 scale (2.75") for every band; W4 L1b frames are 5.52"/pixel
    (508 x 508 native), so the v1 W4 kernel was twice too wide in
    pixels — a sensitivity loss, recorded in the v3.1 erratum. v2
    callers pass the header scale (``pix_scale_from_header``)."""
    fwhm_pix = PSF_FWHM_ARCSEC[band] / pix_arcsec
    flux, var, good_frac, denom, kern, bg = _matched_filter_full(
        img, var_pix, good, fwhm_pix, True, kernel)
    fm = FluxMap(flux=flux, var=var, good_frac=good_frac, wcs=wcs,
                 band=band, magzp=magzp, mjd=mjd)
    if keep_inputs:
        fm.denom, fm.kernel, fm.good, fm.bg = denom, kern, good, bg
    return fm


def build_flux_map(int_path: Path, unc_path: Path, msk_full_path: Path,
                   band: str, mjd: float, fatal_mask: int,
                   magzp: float | None = None, inject=None,
                   keep_inputs: bool = False,
                   pix_scale_from_header: bool = False) -> FluxMap:
    """Build the matched-filter maps for one cutout.

    ``inject``, if given, is called as ``inject(cutout)`` with the
    :class:`WiseCutout` BEFORE the estimator runs and must return the
    image array to search (typically ``cutout.img + source``); masks and
    variances are untouched, so an injected source on a fatal pixel is
    lost as it would be in reality (v2 plan §4.1).
    """
    cut = read_wise_cutout(int_path, unc_path, msk_full_path, fatal_mask)
    img = cut.img if inject is None else np.asarray(inject(cut), dtype=float)
    pix = PIX_ARCSEC
    if pix_scale_from_header and "PXSCAL2" in cut.header:
        pix = abs(float(cut.header["PXSCAL2"]))
    fm = flux_map_from_arrays(img, cut.var_pix, cut.good, cut.wcs, band,
                              mjd, magzp, keep_inputs=keep_inputs,
                              pix_arcsec=pix)
    fm.pix_arcsec = pix
    if keep_inputs:
        fm.frame_origin = cut.frame_origin
        fm.header = cut.header
    return fm


def matched_filter(img: np.ndarray, var_pix: np.ndarray, good: np.ndarray,
                   fwhm_pix: float, subtract_background: bool = True,
                   kernel: np.ndarray | None = None):
    """Core estimator: returns (flux, var, good_frac) maps.

    ``kernel`` overrides the Gaussian model with an arbitrary unit-sum
    PSF (odd-sized, centred); it is flipped internally so that the FFT
    product computes a cross-correlation for asymmetric profiles.
    """
    return _matched_filter_full(img, var_pix, good, fwhm_pix,
                                subtract_background, kernel)[:3]


def _matched_filter_full(img, var_pix, good, fwhm_pix,
                         subtract_background=True, kernel=None):
    """As :func:`matched_filter` plus (denom, flipped kernel, bg)."""
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
    if kernel is None:
        kernel = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
    else:
        kernel = np.asarray(kernel, dtype=float)
        kernel = kernel[::-1, ::-1] / kernel.sum()
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
    return flux, var, good_frac, denom, kernel, bg


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
                       mjd: float, keep_inputs: bool = False,
                       inject=None) -> FluxMap:
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
    if inject is not None:
        # v2 image-level injection: the source is added to the search
        # image (difference image where it exists, sky-subtracted science
        # image elsewhere) BEFORE the estimator; masks/variances untouched
        img = np.asarray(inject(img, hdr, wcs), dtype=float)

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
    flux, var, good_frac, denom, kern, bg = _matched_filter_full(
        img, var_pix, good, fwhm_pix, subtract_background=False)
    kernel = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
    with np.errstate(divide="ignore", invalid="ignore"):
        dfrac = _xcorr((good & is_diff).astype(float), kernel) / np.maximum(
            _xcorr(good.astype(float), kernel), 1e-9)
    fm = FluxMap(flux=flux, var=var, good_frac=good_frac, wcs=wcs,
                 band=band, magzp=(float(hdr["MAGZP"])
                                   if "MAGZP" in hdr else None),
                 mjd=mjd, aux=np.clip(dfrac, 0.0, 1.0))
    fm.fwhm_arcsec = fwhm_arcsec
    fm.pix_arcsec = ZTF_PIX_ARCSEC
    fm.header = hdr
    if keep_inputs:
        fm.denom, fm.kernel, fm.good, fm.bg = denom, kern, good, bg
    return fm


SPHEREX_PIX_ARCSEC = 6.15
ARCSEC2_TO_SR = (np.pi / 648000.0) ** 2


def spherex_kernel(psf_plane: np.ndarray, oversamp: int = 10,
                   shift_pix: tuple[float, float] = (0.0, 0.0)) -> np.ndarray:
    """Bin a 10x-oversampled SPHEREx PSF plane (101 x 101 at 0.615",
    centre at index 50) to detector sampling, keeping the PSF centre at
    the centre of the middle detector pixel plus ``shift_pix`` (x, y)
    detector pixels (sub-pixel phases). Returns a unit-sum odd kernel
    (11 x 11 for the Level-2 products)."""
    p = np.asarray(psf_plane, dtype=float)
    sx = int(round(shift_pix[0] * oversamp))
    sy = int(round(shift_pix[1] * oversamp))
    if sx or sy:
        q = np.zeros_like(p)
        ys, xs = slice(max(sy, 0), p.shape[0] + min(sy, 0)), slice(max(sx, 0), p.shape[1] + min(sx, 0))
        yd, xd = slice(max(-sy, 0), p.shape[0] + min(-sy, 0)), slice(max(-sx, 0), p.shape[1] + min(-sx, 0))
        q[ys, xs] = p[yd, xd]
        p = q
    n = p.shape[0]
    c = n // 2
    half = oversamp // 2
    # pad so that the central bin spans [c-half, c+half)
    lo = (c - half) % oversamp
    pad_lo = (oversamp - lo) % oversamp
    p = np.pad(p, ((pad_lo, 0), (pad_lo, 0)))
    pad_hi = (oversamp - p.shape[0] % oversamp) % oversamp
    p = np.pad(p, ((0, pad_hi), (0, pad_hi)))
    m = p.shape[0]
    k = p.reshape(m // oversamp, oversamp, m // oversamp, oversamp).sum(axis=(1, 3))
    if k.shape[0] % 2 == 0:  # make odd, keeping the peak central
        iy, ix = np.unravel_index(k.argmax(), k.shape)
        k = k[1:, 1:] if iy > (k.shape[0] - 1) // 2 else k[:-1, :-1]
    k = np.clip(k, 0.0, None)
    return k / k.sum()


def build_flux_map_spherex(cut_path: Path, fatal_mask: int, band: str,
                           mjd: float, use_psf: bool = True,
                           upsample: int = 2, keep_inputs: bool = False,
                           inject=None) -> FluxMap:
    """SPHEREx variant on a slim cutout (adapter ``irsa_spherex``).

    Search image = IMAGE - ZODI (model) with a robust constant residual
    background removed by the estimator; per-pixel variance from the
    VARIANCE extension; usable pixels = finite and no fatal FLAGS bit;
    kernel = the exposure's PSF-zone plane binned to detector sampling
    (Gaussian at the header PSF_FWHM if ``use_psf`` is False).

    Units: the matched-filter amplitude with a unit-sum kernel is the
    source's total surface-brightness sum (MJy/sr summed over pixels);
    it is converted to micro-Jansky with the per-exposure median pixel
    solid angle (header OMEGA_MEDIAN, arcsec^2), so ``flux``/``var`` are
    in uJy / uJy^2 and AB = 23.9 - 2.5 log10(flux). ``magzp`` = 23.9.
    ``aux`` holds the wavelength (um) map is NOT stored here — it is a
    smooth function of detector position handled by the caller.
    """
    from astropy.io import fits
    from astropy.utils.exceptions import AstropyWarning
    from astropy.wcs import WCS

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        with fits.open(cut_path) as hdul:
            hdr = hdul["IMAGE"].header
            img = np.asarray(hdul["IMAGE"].data, dtype=np.float64)
            flags = np.asarray(hdul["FLAGS"].data, dtype=np.int64)
            var = np.asarray(hdul["VARIANCE"].data, dtype=np.float64)
            zodi = np.asarray(hdul["ZODI"].data, dtype=np.float64)
            psf = np.asarray(hdul["PSF"].data, dtype=np.float64)
            wcs = WCS(hdr)
    good = (np.isfinite(img) & np.isfinite(var) & (var > 0)
            & ((flags & fatal_mask) == 0))
    omega_sr = float(hdr.get("OMEGA_MEDIAN", PIX_ARCSEC ** 2)) * ARCSEC2_TO_SR
    to_ujy = omega_sr * 1e12  # MJy/sr * sr -> MJy -> uJy
    fwhm_pix = float(hdr.get("PSF_FWHM", 5.3)) / SPHEREX_PIX_ARCSEC
    resid = img - zodi
    if inject is not None:
        # v2 image-level injection (MJy/sr surface brightness added to
        # the zodi-subtracted image; masks/variance untouched)
        resid = np.asarray(inject(resid, hdr, wcs, psf, omega_sr), dtype=float)
    # Empirical variance: the pipeline VARIANCE plane omits confusion /
    # zodi-model residuals in sparse fields and is conservative by up to
    # ~3x in the deep fields (template reduced chi2 0.12 before this
    # step); rescale to the robust (MAD) scatter of the residual image.
    gp = resid[good]
    var_scale = 1.0
    if gp.size > 100:
        med = np.median(gp)
        mad = 1.4826 * np.median(np.abs(gp - med))
        model = float(np.sqrt(np.median(var[good])))
        if mad > 0 and model > 0:
            var_scale = float((mad / model) ** 2)
    # Under-sampled PSF: evaluate the estimator at UP x UP sub-pixel
    # phases (kernel shifted before binning) and interleave into a map
    # with UP x finer sampling, so that bilinear lookup between nodes
    # loses <~ 0.1 mag instead of up to 0.5 mag at half-pixel offsets
    # (measured with the star control, 2026-08-20).
    UP = upsample if (use_psf and psf.ndim == 2) else 1
    ny, nx = resid.shape
    flux = np.empty((ny * UP, nx * UP))
    v = np.empty_like(flux)
    good_frac = np.empty_like(flux)
    kernel0 = None
    phase_parts = []
    for dy in range(UP):
        for dx in range(UP):
            if UP == 1 and not (use_psf and psf.ndim == 2):
                kern = None
            else:
                kern = spherex_kernel(psf, shift_pix=(dx / UP, dy / UP))
                kernel0 = kern if (dx == 0 and dy == 0) else kernel0
            f_, v_, g_, den_, kf_, bg_ = _matched_filter_full(
                resid, var * var_scale, good, fwhm_pix, True, kern)
            flux[dy::UP, dx::UP] = f_
            v[dy::UP, dx::UP] = v_
            good_frac[dy::UP, dx::UP] = g_
            if keep_inputs:
                phase_parts.append((dy, dx, kf_, den_))
    fm = FluxMap(flux=flux * to_ujy, var=v * to_ujy * to_ujy,
                 good_frac=good_frac, wcs=wcs, band=band, magzp=23.9,
                 mjd=mjd, sip=True, upsample=UP)
    fm.var_scale = var_scale
    fm.kernel = kernel0
    fm.omega_sr = omega_sr
    fm.to_ujy = to_ujy
    fm.psf_plane = psf
    fm.pix_arcsec = SPHEREX_PIX_ARCSEC
    fm.fwhm_arcsec = float(hdr.get("PSF_FWHM", 5.3))
    fm.header = hdr
    if keep_inputs:
        fm.phase_parts = phase_parts   # (dy, dx, flipped kernel, denom) per sub-pixel phase
        fm.good = good
        fm.denom = phase_parts[0][3] if phase_parts else None
    return fm


PS1_PIX_ARCSEC = 0.25


def build_flux_map_ps1(img_path: Path, wt_path: Path | None,
                       msk_full_path: Path, fatal_mask: int, band: str,
                       mjd: float, magzp: float | None = None,
                       keep_inputs: bool = False, inject=None) -> FluxMap:
    """Pan-STARRS1 warp variant (adapter ``mast_ps1``).

    Inputs are a fitscut image cutout (full skycell WCS with shifted
    CRPIX), the matching weight cutout (optional) and the FULL skycell
    mask (fpack), aligned to the cutout by the integer CRPIX offset as in
    the WISE builder. Usable pixels = finite image, no fatal mask bit,
    positive variance.

    Variance: the ``.wt`` plane's convention is checked empirically per
    cutout against the robust background scatter of the image — it is
    used as variance if the two agree within a factor 3 (``var_source =
    "wt"``), as inverse variance if 1/wt agrees (``"1/wt"``), otherwise
    the ZTF-style robust-background + Poisson model (``"robust"``). The
    choice is recorded on the returned map.

    Kernel: Gaussian at the header ``CHIP.SEEING`` FWHM (pixels, clipped
    to 2-16). ``magzp`` defaults to header ``FPA.ZP``; the sampling stage
    replaces it with the in-frame star calibration, which also absorbs
    the Gaussian-vs-true-PSF throughput of this estimator.
    """
    img, hdr, wcs = _read_fits_any(img_path)
    msk, mhdr, _ = _read_fits_any(msk_full_path)
    msk = np.where(np.isfinite(msk), msk, 0).astype(np.int64)
    dx = int(round(mhdr["CRPIX1"] - hdr["CRPIX1"]))
    dy = int(round(mhdr["CRPIX2"] - hdr["CRPIX2"]))
    ny, nx = img.shape
    mcut = np.full((ny, nx), int(fatal_mask), dtype=np.int64)
    y0, y1 = max(0, -dy), min(ny, msk.shape[0] - dy)
    x0, x1 = max(0, -dx), min(nx, msk.shape[1] - dx)
    if y1 > y0 and x1 > x0:
        mcut[y0:y1, x0:x1] = msk[y0 + dy:y1 + dy, x0 + dx:x1 + dx]
    unmasked = np.isfinite(img) & ((mcut & fatal_mask) == 0)

    pix = img[unmasked]
    sky = float(np.median(pix)) if pix.size else 0.0
    mad = 1.4826 * float(np.median(np.abs(pix - sky))) if pix.size else 1.0
    bgvar = mad * mad if mad > 0 else 1.0
    gain = float(hdr.get("CELL.GAIN", hdr.get("HIERARCH CELL.GAIN", 1.0))
                 or 1.0)
    var_source = "robust"
    var_pix = bgvar + np.clip(img - sky, 0.0, None) / gain
    if wt_path is not None:
        wt, whdr, _ = _read_fits_any(wt_path)
        if (abs(whdr["CRPIX1"] - hdr["CRPIX1"]) > 0.01
                or whdr["NAXIS1"] != hdr["NAXIS1"]):
            raise ValueError("weight cutout grid differs from image")
        wok = unmasked & np.isfinite(wt) & (wt > 0)
        if wok.sum() > 1000:
            w_med = float(np.median(wt[wok]))
            if 1 / 3 < w_med / bgvar < 3:
                var_pix = np.where(wok, wt, np.nan)
                var_source = "wt"
            elif 1 / 3 < (1.0 / w_med) / bgvar < 3:
                var_pix = np.where(wok, 1.0 / wt, np.nan)
                var_source = "1/wt"
        unmasked &= np.isfinite(var_pix) & (var_pix > 0)
    good = unmasked
    see = float(hdr.get("CHIP.SEEING", hdr.get("HIERARCH CHIP.SEEING", 5.0))
                or 5.0)
    fwhm_pix = float(np.clip(see, 2.0, 16.0))
    if inject is not None:
        img = np.asarray(inject(img, hdr, wcs), dtype=float)
    flux, var, good_frac, denom, kern, bg = _matched_filter_full(
        img, np.nan_to_num(var_pix, nan=1e30), good, fwhm_pix, True)
    if magzp is None:
        zp = hdr.get("FPA.ZP", hdr.get("HIERARCH FPA.ZP"))
        magzp = float(zp) if zp is not None else None
    fm = FluxMap(flux=flux, var=var, good_frac=good_frac, wcs=wcs,
                 band=band, magzp=magzp, mjd=mjd)
    fm.var_source = var_source
    fm.fwhm_pix = fwhm_pix
    fm.sky = sky
    fm.bg_sigma = mad
    fm.exptime = float(hdr.get("EXPTIME", 0.0) or 0.0)
    fm.mjd_obs = float(hdr.get("MJD-OBS", mjd))
    fm.pix_arcsec = PS1_PIX_ARCSEC
    fm.fwhm_arcsec = fwhm_pix * PS1_PIX_ARCSEC
    fm.header = hdr
    if keep_inputs:
        fm.denom, fm.kernel, fm.good, fm.bg = denom, kern, good, bg
    return fm


DECAM_PIX_ARCSEC = 0.2637


def _coarse_inverse_map(canvas_wcs, x0, x1, y0, y1, ccd_wcs, step=64):
    """Map every canvas pixel in [x0:x1, y0:y1] to CCD pixel coords via
    a coarse grid of exact transforms (TAN -> sky -> TPV inverse) and
    bilinear upsampling — the distortion is smooth, so node spacing of
    ``step`` px keeps the interpolation error << 0.05 px."""
    gx = np.arange(x0, x1 + step, step, dtype=float)
    gy = np.arange(y0, y1 + step, step, dtype=float)
    GX, GY = np.meshgrid(gx, gy)
    sky = canvas_wcs.wcs_pix2world(
        np.stack([GX.ravel(), GY.ravel()], axis=1), 0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pix = ccd_wcs.all_world2pix(sky, 0, quiet=True)
    NX = pix[:, 0].reshape(GX.shape)
    NY = pix[:, 1].reshape(GX.shape)
    xs = np.arange(x0, x1 + 1, dtype=float)
    ys = np.arange(y0, y1 + 1, dtype=float)
    fx = (xs - gx[0]) / step
    fy = (ys - gy[0]) / step
    ix = np.clip(fx.astype(int), 0, len(gx) - 2)
    iy = np.clip(fy.astype(int), 0, len(gy) - 2)
    tx = (fx - ix)[None, :]
    ty = (fy - iy)[:, None]

    def up(a):
        a00 = a[np.ix_(iy, ix)]
        a01 = a[np.ix_(iy, ix + 1)]
        a10 = a[np.ix_(iy + 1, ix)]
        a11 = a[np.ix_(iy + 1, ix + 1)]
        return (a00 * (1 - tx) * (1 - ty) + a01 * tx * (1 - ty)
                + a10 * (1 - tx) * ty + a11 * tx * ty)

    return up(NX), up(NY)


def build_flux_map_decam(ccd_files, dqmask_path: Path, fatal_mask: int,
                         band: str, mjd: float, centre_radec,
                         size_pix: int, magzp: float | None = None,
                         keep_inputs: bool = False, inject=None
                         ) -> FluxMap | None:
    """DECam instcal variant (adapter ``noirlab_decam``).

    ``ccd_files`` is a list of (image_path, wtmap_path|None, extname)
    single-CCD HDU files of ONE exposure; ``dqmask_path`` is the full
    dqmask. Each CCD is matched-filtered in its own TPV frame (an
    ``inject`` callback also runs per CCD, in the image domain, before
    the filter), then flux/var/good/denom are pasted onto a common TAN
    canvas at the native pixel scale by exact-inverse nearest-neighbour
    mapping (hypotheses v1.0 §8.7; placement error <~ 0.15").

    Variance follows the PS1 empirical convention check (wt as variance
    vs inverse variance vs robust background+Poisson). The kernel is a
    Gaussian at the per-CCD header FWHM (pixels); ``magzp`` should be
    the per-frame star calibration — header MAGZERO is unreliable
    (flux_scale_check.json) and is only the fallback.
    """
    from astropy.io import fits
    from astropy.wcs import WCS

    n = int(size_pix)
    cw = WCS(naxis=2)
    cw.wcs.crpix = [n / 2.0 + 0.5, n / 2.0 + 0.5]
    cw.wcs.crval = [float(centre_radec[0]), float(centre_radec[1])]
    cw.wcs.cdelt = [-DECAM_PIX_ARCSEC / 3600.0, DECAM_PIX_ARCSEC / 3600.0]
    cw.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    C_flux = np.full((n, n), np.nan)
    C_var = np.full((n, n), np.nan)
    C_good = np.zeros((n, n))
    C_denom = np.full((n, n), np.nan) if keep_inputs else None
    C_goodpix = np.zeros((n, n), dtype=bool) if keep_inputs else None

    dq = fits.open(dqmask_path)
    dq_index = {h.header.get("EXTNAME"): i
                for i, h in enumerate(dq[1:], start=1)}
    prim = dict(dq[0].header)
    kern = None
    fwhm_pix_used, sky_used, mad_used, vsrc_used = [], [], [], []
    n_pasted = 0
    try:
        for img_path, wt_path, extname in ccd_files:
            if extname not in dq_index:
                continue
            with fits.open(img_path) as ih:
                img = np.asarray(ih[1].data, dtype=float)
                hdr = ih[1].header
                phdr = ih[0].header
                wcs = WCS(hdr)
            msk = np.asarray(dq[dq_index[extname]].data)
            if msk.shape != img.shape:
                # CP version trims can differ by a few rows; crop both
                ny = min(msk.shape[0], img.shape[0])
                nx = min(msk.shape[1], img.shape[1])
                msk, img = msk[:ny, :nx], img[:ny, :nx]
            unmasked = np.isfinite(img) & ((msk.astype(np.int64)
                                            & fatal_mask) == 0)
            pix = img[unmasked]
            if pix.size < 1000:
                continue
            sky = float(np.median(pix))
            mad = 1.4826 * float(np.median(np.abs(pix - sky))) or 1.0
            bgvar = mad * mad
            var_pix = bgvar + np.clip(img - sky, 0.0, None)
            var_source = "robust"
            if wt_path is not None and Path(wt_path).exists():
                with fits.open(wt_path) as wh:
                    wt = np.asarray(wh[1].data, dtype=float)
                wt = wt[:img.shape[0], :img.shape[1]]
                wok = unmasked & np.isfinite(wt) & (wt > 0)
                if wok.sum() > 1000:
                    w_med = float(np.median(wt[wok]))
                    with np.errstate(divide="ignore"):
                        if 1 / 3 < w_med / bgvar < 3:
                            var_pix = np.where(wok, wt, np.nan)
                            var_source = "wt"
                        elif 1 / 3 < (1.0 / w_med) / bgvar < 3:
                            var_pix = np.where(wok, 1.0 / wt, np.nan)
                            var_source = "1/wt"
            good = unmasked & np.isfinite(var_pix) & (var_pix > 0)
            fwhm = (hdr.get("FWHM") or phdr.get("FWHM") or 4.0)
            fwhm_pix = float(np.clip(float(fwhm), 2.0, 16.0))
            if inject is not None:
                img = np.asarray(inject(img, phdr, wcs), dtype=float)
            flux, var, good_frac, denom, k, bg = _matched_filter_full(
                img, np.nan_to_num(var_pix, nan=1e30), good, fwhm_pix,
                True)
            if kern is None:
                kern = k
            # paste: canvas bbox of this CCD -> nearest CCD pixel
            ny, nx = img.shape
            corners = wcs.all_pix2world(
                np.array([[0.5, 0.5], [nx - 0.5, 0.5],
                          [nx - 0.5, ny - 0.5], [0.5, ny - 0.5]]), 0)
            cpix = cw.wcs_world2pix(corners, 0)
            x0 = int(np.clip(np.floor(cpix[:, 0].min()) - 2, 0, n - 1))
            x1 = int(np.clip(np.ceil(cpix[:, 0].max()) + 2, 0, n - 1))
            y0 = int(np.clip(np.floor(cpix[:, 1].min()) - 2, 0, n - 1))
            y1 = int(np.clip(np.ceil(cpix[:, 1].max()) + 2, 0, n - 1))
            if x1 <= x0 or y1 <= y0:
                continue
            CX, CY = _coarse_inverse_map(cw, x0, x1, y0, y1, wcs)
            sx = np.rint(CX).astype(int)
            sy = np.rint(CY).astype(int)
            inb = (sx >= 0) & (sx < nx) & (sy >= 0) & (sy < ny)
            tgt_empty = ~np.isfinite(C_flux[y0:y1 + 1, x0:x1 + 1])
            sel = inb & tgt_empty
            if not sel.any():
                continue
            view = np.s_[y0:y1 + 1, x0:x1 + 1]
            C_flux[view][sel] = flux[sy[sel], sx[sel]]
            C_var[view][sel] = var[sy[sel], sx[sel]]
            C_good[view][sel] = good_frac[sy[sel], sx[sel]]
            if keep_inputs:
                C_denom[view][sel] = denom[sy[sel], sx[sel]]
                C_goodpix[view][sel] = good[sy[sel], sx[sel]]
            n_pasted += 1
            fwhm_pix_used.append(fwhm_pix)
            sky_used.append(sky)
            mad_used.append(mad)
            vsrc_used.append(var_source)
    finally:
        dq.close()
    if n_pasted == 0 or kern is None:
        return None
    if magzp is None:
        mz = prim.get("MAGZERO")
        magzp = float(mz) if mz is not None else None
    fm = FluxMap(flux=C_flux, var=C_var, good_frac=C_good, wcs=cw,
                 band=band, magzp=magzp, mjd=mjd)
    fm.var_source = vsrc_used[0]
    fm.fwhm_pix = float(np.mean(fwhm_pix_used))
    fm.sky = float(np.mean(sky_used))
    fm.bg_sigma = float(np.mean(mad_used))
    fm.exptime = float(prim.get("EXPTIME", 0.0) or 0.0)
    fm.mjd_obs = float(prim.get("MJD-OBS", mjd))
    fm.pix_arcsec = DECAM_PIX_ARCSEC
    fm.fwhm_arcsec = fm.fwhm_pix * DECAM_PIX_ARCSEC
    fm.header = prim
    fm.n_ccds = n_pasted
    if keep_inputs:
        fm.denom, fm.kernel, fm.good, fm.bg = (C_denom, kern, C_goodpix,
                                               float(np.mean(sky_used)))
    return fm
