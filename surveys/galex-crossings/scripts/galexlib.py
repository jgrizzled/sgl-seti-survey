"""Shared machinery for the GALEX crossings survey (freeze v1.0+v1.1,
threshold freeze v1.0): gated aspect/photon fetches, the three frozen
statistics, control constructions, and photon-level injections.

Off-window discipline: every dev-stage caller passes the target-channel
window list to `assert_off_window`; nothing here touches in-window
seconds until the blind confirmatory run.
"""

from __future__ import annotations

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

KM_PER_AU = 1.495978707e8
RSUN_KM = 695_700.0
Z_GRID_AU = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
USABLE_ARCMIN = 33.0
APERTURE_ARCSEC = 8.0
BURST_WIDTHS_S = (0.05, 0.5, 5.0, 50.0)
PERIOD_MIN_S = 0.02
PERIOD_MAX_FRACTION = 1.0 / 3.0
PERIOD_OVERSAMPLE = 5
HTEST_MAX_HARMONICS = 20
#: worst-case LEO line-of-sight drift rate (v_orb / c), freeze s4
ORBIT_DRIFT_RATE = 2.56e-5


# -- geometry -----------------------------------------------------------

def relay_apparent(star_ra_deg, star_dec_deg, z_au, t_mjd):
    t = Time(t_mjd, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", t).xyz.to_value("AU")
    earth = get_body_barycentric("earth", t).xyz.to_value("AU")
    ra = np.radians(star_ra_deg)
    de = np.radians(star_dec_deg)
    u_star = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra),
                       np.sin(de)])
    v = (sun - z_au * u_star) - earth
    v /= np.linalg.norm(v)
    return (float(np.degrees(np.arctan2(v[1], v[0])) % 360.0),
            float(np.degrees(np.arcsin(np.clip(v[2], -1, 1)))))


def ang_arcmin(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = np.radians([ra1, de1, ra2, de2])
    c = (np.sin(d1) * np.sin(d2)
         + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))) * 60.0


def windows_ms(events_tab, target_id, link_direction, r_au,
               mjd_to_ms):
    """All in-era flat-chord windows (ms) for one target-channel."""
    out = []
    m = ((events_tab["target_id"] == target_id)
         & (events_tab["link_direction"] == link_direction))
    for ev in events_tab[m]:
        b = float(ev["b_min_au"])
        if b >= r_au:
            continue
        half_d = (np.sqrt(r_au ** 2 - b ** 2) * KM_PER_AU
                  / float(ev["v_perp_km_s"]) / 86400.0)
        t_ca = float(ev["t_ca_mjd"])
        out.append((mjd_to_ms(t_ca - half_d), mjd_to_ms(t_ca + half_d)))
    return out


def assert_off_window(t0_ms, t1_ms, windows):
    for w0, w1 in windows:
        if t0_ms < w1 and t1_ms > w0:
            raise RuntimeError(
                f"off-window violation: [{t0_ms},{t1_ms}] overlaps "
                f"window [{w0},{w1}]")


def star_track(events_tab, target_id, link_direction="inbound"):
    """Linear (ra, dec)(mjd) interpolator from the per-event star
    positions (the registry's propagated astrometry sampled 2x/yr)."""
    m = ((events_tab["target_id"] == target_id)
         & (events_tab["link_direction"] == link_direction))
    sub = events_tab[m]
    t = np.asarray(sub["t_ca_mjd"], float)
    ra = np.asarray(sub["star_icrs_ra_deg"], float)
    de = np.asarray(sub["star_icrs_dec_deg"], float)
    ix = np.argsort(t)
    t, ra, de = t[ix], ra[ix], de[ix]
    pr = np.polyfit(t - t.mean(), ra, 1)
    pd = np.polyfit(t - t.mean(), de, 1)

    def at(mjd):
        return (float(np.polyval(pr, mjd - t.mean())),
                float(np.polyval(pd, mjd - t.mean())))
    return at


# -- gated fetches ------------------------------------------------------

def sec_id(t_ms):
    return (np.asarray(t_ms, dtype=np.int64) - 995) // 1000


def usable_seconds(client, store, t0_ms, t1_ms, ra, dec):
    """Usable aspect second ids for a position over a visit range:
    boresight <= 33', flag % 2 == 0 (v1.1), deduplicated."""
    rows = client.aspect_range(t0_ms, t1_ms, store)
    out = {}
    for t_ms, ra0, de0, band, flag in rows:
        s = int(sec_id(t_ms))
        if s in out:
            continue
        if flag % 2 != 0:
            continue
        if ang_arcmin(ra0, de0, ra, dec) > USABLE_ARCMIN:
            continue
        out[s] = band
    return out          # {second_id: band_string}


def photons_aperture(client, store, band, ra, dec, t0_ms, t1_ms,
                     usable, r_arcsec=APERTURE_ARCSEC):
    """Photon arrival times (ms, sorted) in an aperture, gated to
    usable seconds whose band includes `band`; photon flag == 0."""
    half = (r_arcsec + 3.0) / 3600.0
    rows = client.photons_box(band, ra, dec, half, t0_ms, t1_ms, store)
    ok_secs = {s for s, b in usable.items() if band in b}
    times = []
    cosd = np.cos(np.radians(dec))
    for t_ms, pra, pde, flag in rows:
        if flag != 0:
            continue
        if int(sec_id(t_ms)) not in ok_secs:
            continue
        d = np.hypot((pra - ra) * cosd, pde - dec) * 3600.0
        if d <= r_arcsec:
            times.append(t_ms)
    return np.array(sorted(times), dtype=np.int64)


# -- frozen statistics --------------------------------------------------

def s_rate(times_ms, t_live_s):
    return len(times_ms) / t_live_s if t_live_s > 0 else np.nan


def s_burst(times_ms, t_live_s):
    """Max boxcar z over the frozen widths (step w/2)."""
    if len(times_ms) == 0 or t_live_s <= 0:
        return 0.0
    t = np.asarray(times_ms, float) / 1000.0
    lam = len(t) / t_live_s
    best = 0.0
    t0, t1 = t.min(), t.max()
    for w in BURST_WIDTHS_S:
        starts = np.arange(t0 - w, t1 + w / 2.0, w / 2.0)
        left = np.searchsorted(t, starts)
        right = np.searchsorted(t, starts + w)
        cmax = int((right - left).max()) if len(starts) else 0
        mu = lam * w
        z = (cmax - mu) / np.sqrt(max(mu, 0.5))
        best = max(best, float(z))
    return best


def period_grid(t_span_s):
    f_min = 1.0 / (t_span_s * PERIOD_MAX_FRACTION)
    f_max = 1.0 / PERIOD_MIN_S
    df = 1.0 / (PERIOD_OVERSAMPLE * t_span_s)
    return np.arange(f_min, f_max + df, df)


def s_period(times_ms, t_span_s, min_photons=10, jitter_rng=None):
    """Max H (de Jager) over the frozen frequency grid; harmonic
    recurrence keeps it to two trig calls per (photon, frequency).
    jitter_rng: de-quantization jitter U(0, 5 ms) added to the 5 ms
    photon ticks before phase folding (amendment v1.2) — without it
    every frequency whose harmonic is commensurate with the 200 Hz
    tick rate accumulates coherent quantization power."""
    t = np.asarray(times_ms, float) / 1000.0
    if len(t) < min_photons:
        return None
    if jitter_rng is not None:
        t = t + jitter_rng.uniform(0.0, 0.005, len(t))
    t = t - t.min()
    freqs = period_grid(t_span_s)
    n = len(t)
    chunk = max(64, int(5e6 // max(n, 1)))
    best = -np.inf
    for i in range(0, len(freqs), chunk):
        f = freqs[i:i + chunk]
        phi = 2.0 * np.pi * np.outer(f, t)          # (F, N)
        c1, s1 = np.cos(phi), np.sin(phi)
        ck, sk = c1.copy(), s1.copy()
        z2 = np.zeros(len(f))
        h = np.full(len(f), -np.inf)
        for k in range(1, HTEST_MAX_HARMONICS + 1):
            z2 += (ck.sum(axis=1) ** 2 + sk.sum(axis=1) ** 2)
            hk = (2.0 / n) * z2 - 4.0 * (k - 1)
            h = np.maximum(h, hk)
            # angle-addition recurrence: cos/sin((k+1)phi)
            ck, sk = ck * c1 - sk * s1, sk * c1 + ck * s1
        best = max(best, float(h.max()))
    return best


# -- injections ---------------------------------------------------------

def inject_persistent(rng, rate_cts_s, secs_ms_sorted):
    """Poisson arrivals at rate over the usable seconds (ms times)."""
    out = []
    for s in secs_ms_sorted:
        k = rng.poisson(rate_cts_s)
        out.extend(s + np.floor(rng.uniform(0, 1000, k) / 5) * 5)
    return np.array(sorted(out), dtype=np.int64)


def inject_pulse(rng, n_photons, width_s, secs_ms_sorted):
    """One boxcar pulse of n photons at a random usable time,
    gated to the usable seconds."""
    base = rng.choice(secs_ms_sorted)
    times = (base + np.floor(
        rng.uniform(0, width_s * 1000.0, n_photons) / 5) * 5
        ).astype(np.int64)
    usable = set((np.asarray(secs_ms_sorted) - 995) // 1000)
    times = times[np.isin((times - 995) // 1000, list(usable))]
    return np.array(sorted(times), dtype=np.int64)


def inject_train(rng, period_s, frac_per_cycle, duty, secs_ms_sorted,
                 drift_rate=None):
    """Periodic boxcar train: mean frac_per_cycle photons per cycle in
    a duty*period window, over the usable span; optional linear
    arrival-time drift (the frozen orbital-smear model)."""
    t0 = secs_ms_sorted[0]
    span_s = (secs_ms_sorted[-1] + 1000 - t0) / 1000.0
    usable = set((np.asarray(secs_ms_sorted) - 995) // 1000)
    out = []
    n_cycles = int(span_s / period_s) + 1
    for i in range(n_cycles):
        k = rng.poisson(frac_per_cycle)
        base = i * period_s
        ts = base + rng.uniform(0, duty * period_s, k)
        out.extend(ts)
    out = np.asarray(sorted(out), float)
    if drift_rate is not None:
        out = out + drift_rate * (out - out[0])
    t_ms = t0 + np.floor(out * 200.0) * 5    # 5 ms ticks
    keep = np.isin((t_ms.astype(np.int64) - 995) // 1000,
                   list(usable))
    return t_ms[keep].astype(np.int64)


def merge_series(*arrs):
    return np.array(sorted(np.concatenate(arrs)), dtype=np.int64)
