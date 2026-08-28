"""Rubin DP2 recon: image access end-to-end at the ross-128 antipode.

SIA (dp2) query -> DataLink resolve -> SODA circle cutout -> verify the
returned FITS opens and carries WCS. Saves the cutout + a probe log.

Output: surveys/rubin/results/recon_image_probe_v0.json
        runs/rubin/recon/ross128_antipode_coadd_cutout.fits
"""
import io
import json
import pathlib
import time

import requests
from astropy.io import fits
from astropy.wcs import WCS

REPO = pathlib.Path(__file__).resolve().parents[3]
OUT = REPO / "surveys" / "rubin" / "results" / "recon_image_probe_v0.json"
CUT = REPO / "runs" / "rubin" / "recon" / "ross128_antipode_coadd_cutout.fits"
RA, DEC = 356.943, -0.788

token = None
for line in (REPO / ".env").read_text().splitlines():
    if line.startswith("RSP_API_TOKEN="):
        token = line.split("=", 1)[1].strip().strip('"')
sess = requests.Session()
sess.headers["Authorization"] = f"Bearer {token}"

out = {"pos": [RA, DEC]}

# 1. SIA query for deep coadds at the position
t0 = time.time()
r = sess.get("https://data.lsst.cloud/api/sia/dp2/query",
             params={"POS": f"CIRCLE {RA} {DEC} 0.01",
                     "DPSUBTYPE": "lsst.deep_coadd", "MAXREC": "10"},
             timeout=60)
r.raise_for_status()
from astropy.io.votable import parse_single_table  # noqa: E402
tab = parse_single_table(io.BytesIO(r.content)).to_table()
out["sia"] = {"elapsed_s": round(time.time() - t0, 1), "n_rows": len(tab),
              "bands": sorted(set(str(b) for b in tab["lsst_band"]))
              if "lsst_band" in tab.colnames else None,
              "colnames_sample": tab.colnames[:25]}
print("SIA rows:", len(tab))

# 2. DataLink resolve on the first r-band row (or first row)
row = None
for cand in tab:
    if str(cand.get("lsst_band", "")) == "r":
        row = cand
        break
row = row if row is not None else tab[0]
dl_url = str(row["access_url"])
out["datalink_url"] = dl_url
t0 = time.time()
r = sess.get(dl_url, timeout=60)
r.raise_for_status()
dl = parse_single_table(io.BytesIO(r.content)).to_table()
out["datalink"] = {
    "elapsed_s": round(time.time() - t0, 1),
    "rows": [{"id_or_service": str(x.get("service_def", "") or ""),
              "semantics": str(x.get("semantics", "")),
              "content_type": str(x.get("content_type", ""))}
             for x in dl],
}
print("datalink rows:", len(dl))

# find the SODA cutout service descriptor in the datalink resources
from astropy.io.votable import parse as vot_parse  # noqa: E402
r2 = sess.get(dl_url, timeout=60)
vot = vot_parse(io.BytesIO(r2.content))
soda_url, dl_id = None, None
for res in vot.resources:
    for p in getattr(res, "params", []):
        if p.name == "accessURL":
            soda_url = p.value
for x in dl:
    if str(x.get("semantics", "")) == "#this":
        dl_id = str(x["ID"]) if "ID" in dl.colnames else None
        out["direct_fits_url"] = str(x["access_url"])
out["soda_url"] = soda_url

# 3. SODA circle cutout (2 arcmin)
if soda_url:
    t0 = time.time()
    r = sess.get(soda_url,
                 params={"ID": dl_id, "CIRCLE": f"{RA} {DEC} 0.017"},
                 timeout=300)
    ct = r.headers.get("content-type", "")
    out["soda"] = {"status": r.status_code, "content_type": ct,
                   "bytes": len(r.content),
                   "elapsed_s": round(time.time() - t0, 1)}
    if r.ok and "fits" in ct:
        CUT.write_bytes(r.content)
        with fits.open(CUT) as hdul:
            sci = next(h for h in hdul if h.data is not None)
            w = WCS(sci.header)
            out["cutout"] = {"n_hdus": len(hdul),
                             "shape": list(sci.data.shape),
                             "wcs_type": list(w.wcs.ctype),
                             "hdu_names": [h.name for h in hdul]}
        print("cutout ok:", out["cutout"])
    else:
        out["soda"]["body_head"] = r.text[:400]

OUT.write_text(json.dumps(out, indent=1, default=str))
print("wrote", OUT.relative_to(REPO))
