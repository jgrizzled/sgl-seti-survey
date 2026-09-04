"""Follow-up to quicklook_v1: the gj-1276 star position showed a 2.7σ
pixel (0.64 mJy, rms 0.24) in the VAST_2257-06 SB52549 cutout
(2023-09-03, the in-window epoch). Is it a persistent faint source or
a one-epoch fluctuation/flare? Same position (per-epoch PM is < 2″,
below the 13″ beam) in every VAST_2257-06 epoch CASDA holds: per-SBID
component catalogue cone (60″) + a SODA cutout with the same
photometry as quicklook_v1. Needs the OPAL login (.env)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import quicklook_v1 as ql  # noqa: E402

from sglsurvey.snapshots import SnapshotStore  # noqa: E402

OUT = ql.RUN / "gj1276_epochs"
TARGET = "gj-1276"
EVENT = "evt-ed23b2322fe3"      # 2023-09-05 A 0.1 AU crossing


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    store = SnapshotStore(OUT)
    session = requests.Session()
    auth = ql.opal_auth()
    assert auth, "OPAL credentials required"
    ev = Table.read(ql.REPO / "crossings" / "universal_v1" / "events.ecsv")
    e = ev[ev["event_id"] == EVENT][0]
    ra, de = float(e["star_icrs_ra_deg"]), float(e["star_icrs_dec_deg"])
    # every restored Stokes-I image of the field
    q = ("SELECT obs_id, filename, access_url, t_min, t_max, quality_level "
         "FROM ivoa.obscore WHERE dataproduct_subtype='cont.restored.t0' "
         "AND filename LIKE 'image.i.VAST_2257-06.%restored.conv.fits' "
         "ORDER BY t_min")
    obs = ql.tap_csv(session, store, ql.CASDA_TAP, q, "field-epochs")
    tables = ql.casda_tables(session, store)
    rows = []
    for r in obs:
        sbid = str(r["obs_id"]).replace("ASKAP-", "")
        pid = str(r["access_url"]).split("ID=")[-1]
        rec = {"sbid": sbid, "product": pid, "t_min": ql.fnum(r["t_min"]),
               "quality": str(r["quality_level"]), "filename": str(r["filename"])}
        cand = [t for t in tables if t.startswith("AS207.") and f"_sb{sbid}_" in t
                and t.endswith("components_v01")]
        if cand:
            cone = ql.cone_generic(session, store, ql.CASDA_TAP, cand[0],
                                   "ra_deg_cont", "dec_deg_cont",
                                   "component_id, ra_deg_cont, dec_deg_cont, flux_peak, "
                                   "flux_peak_err, rms_image", ra, de, 60 / 3600,
                                   label=cand[0])
            rec["catalogue_n_60arcsec"] = None if cone is None else len(cone)
            if cone is not None and len(cone):
                rec["catalogue_peak_mjy"] = [float(x) for x in cone["flux_peak"]]
        else:
            rec["catalogue_n_60arcsec"] = "no_table"
        dl = ql.casda_datalink(session, store, auth, pid)
        if "cutout_token" in dl:
            dest = OUT / f"{TARGET}_{sbid}_{pid}.fits"
            if not dest.exists():
                ql.casda_soda_cutout(session, dl["cutout_token"], ra, de,
                                     ql.CUT_RADIUS_DEG, dest, label=pid)
            if dest.exists():
                rec.update(ql.askap_cutout_photometry(dest, ra, de))
                rec["cutout"] = str(dest.relative_to(ql.REPO))
        else:
            rec["cutout_status"] = dl.get("error", "no_cutout_service")
        print(f"{sbid:6s} {rec.get('date_obs')} q={rec['quality']:9s} "
              f"cat60={rec['catalogue_n_60arcsec']} at={rec.get('value_at_position_mjy')} "
              f"peak={rec.get('peak_in_beam_mjy')} rms={rec.get('rms_annulus_mjy')}",
              flush=True)
        rows.append(rec)
    (ql.RES / "quicklook_gj1276_epochs.json").write_text(json.dumps(rows, indent=1, default=str))
    vals = [r["snr_value_at_position"] for r in rows if r.get("snr_value_at_position") is not None]
    print(f"{len(vals)} epochs with cutouts; SNR at position: "
          f"median {np.median(vals):.2f}, max {max(vals):.2f}, min {min(vals):.2f}")


if __name__ == "__main__":
    main()
