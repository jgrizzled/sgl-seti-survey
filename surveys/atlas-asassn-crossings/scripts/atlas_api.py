"""Minimal ATLAS forced-photometry API client for this survey.

Facts baked in from the recon (notes/atlas_asassn_recon_2026-08-25.md):
token auth from ATLAS_API_TOKEN in the repo .env; POST /queue/ returns
201 (queued) or 429 with a machine-parseable wait time; execution is
serial per account (~25-70 min per full-history position, faster for
windowed mjd ranges); result_url is short-lived so results are fetched
to disk immediately and the task DELETEd; result text carries ###
comment prefixes; the Obs column is the per-exposure provenance key.

Snapshot discipline: every fetch writes the raw text plus a sidecar
JSON (query params, task metadata, sha256, timestamps) next to it.
"""
import hashlib
import io
import json
import os
import re
import time
from datetime import datetime, timezone

import pandas as pd
import requests

BASEURL = "https://fallingstar-data.com/forcedphot"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))


def load_token(env_path=None):
    path = env_path or os.path.join(REPO, '.env')
    with open(path) as f:
        for line in f:
            if line.startswith('ATLAS_API_TOKEN'):
                return line.strip().split('=', 1)[1]
    raise RuntimeError(f'ATLAS_API_TOKEN not found in {path}')


def headers(token=None):
    return {'Authorization': f'Token {token or load_token()}',
            'Accept': 'application/json'}


def submit(params, hdrs, max_wait_s=3600):
    """Queue one task; params e.g. {'ra':..,'dec':..,'mjd_min':..,
    'mjd_max':.., 'use_reduced':..}. Honors 429 wait times."""
    waited = 0.0
    while True:
        resp = requests.post(f"{BASEURL}/queue/", headers=hdrs, data=params)
        if resp.status_code == 201:
            return resp.json()
        if resp.status_code == 429 and waited < max_wait_s:
            msg = resp.json().get('detail', '')
            sec = re.findall(r'available in (\d+) second', msg)
            mins = re.findall(r'available in (\d+) minute', msg)
            wait = int(sec[0]) if sec else int(mins[0]) * 60 if mins else 15
            time.sleep(wait + 1)
            waited += wait + 1
            continue
        raise RuntimeError(f'submit {resp.status_code}: {resp.text[:300]}')


def wait_and_fetch(task_url, hdrs, out_path, query_params=None,
                   poll_s=30, timeout_s=6 * 3600, delete=True):
    """Poll a task, save the result text + sidecar snapshot JSON, DELETE
    the task. Returns the sidecar dict."""
    t0 = time.time()
    while True:
        j = requests.get(task_url, headers=hdrs).json()
        if j.get('finishtimestamp'):
            break
        if time.time() - t0 > timeout_s:
            raise TimeoutError(f'{task_url} not finished after {timeout_s}s')
        time.sleep(poll_s)
    return fetch_result(j, hdrs, out_path, query_params, t0=t0, delete=delete)


def fetch_result(j, hdrs, out_path, query_params=None, t0=None, delete=True):
    """Given a finished task's JSON, save result text + sidecar, DELETE
    the task. Returns the sidecar dict."""
    task_url = j['url']
    t0 = t0 or time.time()
    txt = requests.get(j['result_url'], headers=hdrs).text
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(txt)
    sidecar = {
        'query_params': query_params, 'task_meta': j,
        'fetched_utc': datetime.now(timezone.utc).isoformat(),
        'wall_s': round(time.time() - t0, 1),
        'result_sha256': hashlib.sha256(txt.encode()).hexdigest(),
        'result_bytes': len(txt), 'result_path': os.path.relpath(out_path, REPO),
    }
    with open(out_path + '.json', 'w') as f:
        json.dump(sidecar, f, indent=2)
    if delete:
        requests.delete(task_url, headers=hdrs)
    return sidecar


def read_result(path):
    """Parse a saved result text file into a DataFrame."""
    with open(path) as f:
        raw = f.read()
    return pd.read_csv(io.StringIO(raw.replace('###', '')), sep=r'\s+')


# The FAQ cleaning recipe (Rest et al. style), verbatim from the recon.
# Survey-specific masks are frozen in hypotheses.md, not here.
def faq_quality_mask(df):
    return ((df.duJy < 10000) & (df.err == 0)
            & (df.x > 100) & (df.x < 10460) & (df.y > 100) & (df.y < 10460)
            & (df['maj'] < 5) & (df['maj'] > 1.6)
            & (df['min'] < 5) & (df['min'] > 1.6)
            & (df.apfit > -1) & (df.apfit < -0.1)
            & (df.mag5sig > 17) & (df.Sky > 17))
