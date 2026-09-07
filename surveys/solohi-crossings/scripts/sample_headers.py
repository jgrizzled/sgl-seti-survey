"""Sample one solohi-1ft L2 frame per ~10 days from the SOAR inventory,
download it (4 MB each) into runs/solohi-crossings/recon/sample/ and
tabulate header facts (XPOSURE, NSUMEXP, OBS_MODE, GAINMODE, DSUN_OBS,
SC_ROLL, CRVAL, DATASAT, VERSION) -> runs/solohi-crossings/coverage/
header_sample_1ft.json."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, requests
from astropy.io import fits
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
INV = REPO / "runs" / "solohi-crossings" / "coverage" / "soar_l2_inventory.json"
DEST = REPO / "runs" / "solohi-crossings" / "recon" / "sample"
OUT = REPO / "runs" / "solohi-crossings" / "coverage" / "header_sample_1ft.json"
KEYS = ["DATE-AVG", "XPOSURE", "NSUMEXP", "OBS_MODE", "OBS_ID", "GAINMODE", "DSUN_OBS", "SC_ROLL",
        "HGLT_OBS", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "DATASAT", "DATAMAX", "VERSION", "VERS_CAL",
        "NBIN", "NAXIS1", "NAXIS2", "TARGET", "STUDY_ID", "DETECTOR"]


def main(step_days=10.0, desc="solohi-1ft"):
    inv = json.load(open(INV)); ci = {c: i for i, c in enumerate(inv["columns"])}
    rows = [r for r in inv["rows"] if r[ci["descriptor"]] == desc]
    mjd = np.array([Time(r[ci["begin_time"]]).mjd for r in rows])
    picks, last = [], -1e9
    for k in np.argsort(mjd):
        if mjd[k] - last >= step_days:
            picks.append(rows[k]); last = mjd[k]
    DEST.mkdir(parents=True, exist_ok=True)
    out = []
    for r in picks:
        fn = DEST / r[ci["filename"]]
        if not fn.exists():
            resp = requests.get("https://soar.esac.esa.int/soar-sl-tap/data",
                                params={"retrieval_type": "LAST_PRODUCT", "data_item_id": r[ci["data_item_id"]],
                                        "product_type": "SCIENCE"}, timeout=600)
            fn.write_bytes(resp.content)
        try:
            h = fits.getheader(fn)
        except Exception as e:
            print("bad", fn, e, file=sys.stderr); continue
        out.append({"file": fn.name, **{k: (float(h[k]) if isinstance(h[k], (int, float)) else h[k]) for k in KEYS if k in h}})
        print(fn.name, h.get("OBS_MODE"), h.get("XPOSURE"), round(h["DSUN_OBS"] / 1.495978707e11, 3), file=sys.stderr, flush=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(len(out), "sampled ->", OUT)


if __name__ == "__main__":
    main()
