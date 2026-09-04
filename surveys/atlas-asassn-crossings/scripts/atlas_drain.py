"""Resumable ATLAS forced-photometry queue drain for the coverage stage.

Reads results/task_list_v1.ecsv, submits windowed difference-flux
tasks (comment = task_key, the resume handle), polls, fetches each
result to runs/atlas-asassn-crossings/lc/<task_key>.txt with the
atlas_api sidecar snapshot (query params, task metadata, sha256,
timestamps), DELETEs the server task. The sidecars double as the
done-list; in-flight task URLs persist in inflight.json; on restart,
server tasks whose comment matches an unfinished key are adopted (a
crash between submit and persist cannot duplicate work). Failed
tasks get an error sidecar (never retried without --retry-failed).

Order: channel B before A, t_ca positions before the mini-track ends,
grazing targets first, then by t_ca.

  python atlas_drain.py --limit 1          # time one task
  nohup python atlas_drain.py > runs/.../drain.out &
  python atlas_drain.py --status
"""
import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone

import numpy as np
import requests
from astropy.table import Table

sys.path.insert(0, os.path.dirname(__file__))
import atlas_api  # noqa: E402

REPO = atlas_api.REPO
RUN = os.path.join(REPO, 'runs', 'atlas-asassn-crossings')
LC = os.path.join(RUN, 'lc')
INFLIGHT = os.path.join(RUN, 'inflight.json')
LOG = os.path.join(RUN, 'drain_log.jsonl')
TASKS = os.path.join(os.path.dirname(__file__), '..', 'results', 'task_list_v1.ecsv')
TARGET_ORDER = ['van-maanen', 'wolf-359', 'gj-1276', 'teegarden', 'ross-128',
                'ross-154', 'gj-908']


def now():
    return datetime.now(timezone.utc).isoformat()


def out_path(key):
    return os.path.join(LC, f'{key}.txt')


def sidecar_path(key):
    return out_path(key) + '.json'


def load_inflight():
    if os.path.exists(INFLIGHT):
        with open(INFLIGHT) as f:
            return json.load(f)
    return {}


def save_inflight(d):
    tmp = INFLIGHT + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(d, f, indent=1)
    os.replace(tmp, INFLIGHT)


def log(rec):
    rec['utc'] = now()
    with open(LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')
    print(json.dumps(rec), flush=True)


def priority(row):
    return (0 if row['channel'] == 'B' else 1,
            0 if row['pos_role'] == 'tca' else 1,
            TARGET_ORDER.index(row['target_id']), float(row['t_ca_mjd']))


def params_for(row):
    return {'ra': float(row['ra_deg']), 'dec': float(row['dec_deg']),
            'mjd_min': float(row['mjd_min']), 'mjd_max': float(row['mjd_max']),
            'use_reduced': 'true' if bool(row['use_reduced']) else 'false',
            'comment': str(row['task_key'])}


def failed_keys():
    out = []
    for fn in os.listdir(LC) if os.path.isdir(LC) else []:
        if fn.endswith('.txt.json'):
            with open(os.path.join(LC, fn)) as f:
                if 'error' in json.load(f):
                    out.append(fn[:-len('.txt.json')])
    return out


def list_server_tasks(hdrs):
    url = f'{atlas_api.BASEURL}/queue/?pagesize=100'
    tasks = []
    while url:
        j = requests.get(url, headers=hdrs).json()
        tasks.extend(j.get('results', []))
        url = j.get('next')
    return tasks


def adopt_orphans(hdrs, inflight, keys):
    """Adopt server tasks whose comment is one of our unfinished keys
    but which are not in inflight.json; delete duplicate submissions."""
    seen = set(inflight.values())
    for t in list_server_tasks(hdrs):
        key = t.get('comment')
        if key not in keys or t['url'] in seen:
            continue
        if key in inflight or os.path.exists(sidecar_path(key)):
            log({'event': 'delete_duplicate', 'task_key': key, 'url': t['url']})
            requests.delete(t['url'], headers=hdrs)
            continue
        inflight[key] = t['url']
        seen.add(t['url'])
        log({'event': 'adopt', 'task_key': key, 'url': t['url']})
    save_inflight(inflight)


def poll_once(hdrs, inflight, rows_by_key):
    n_done = 0
    for key, url in list(inflight.items()):
        r = requests.get(url, headers=hdrs)
        if r.status_code == 404:
            log({'event': 'vanished', 'task_key': key, 'url': url})
            del inflight[key]
            save_inflight(inflight)
            continue
        j = r.json()
        if not j.get('finishtimestamp'):
            continue
        t0 = time.time()
        if j.get('error_msg') or not j.get('result_url'):
            os.makedirs(LC, exist_ok=True)
            with open(sidecar_path(key), 'w') as f:
                json.dump({'error': j.get('error_msg') or 'no result_url',
                           'task_meta': j, 'query_params': params_for(rows_by_key[key]),
                           'fetched_utc': now()}, f, indent=2)
            requests.delete(url, headers=hdrs)
            log({'event': 'failed', 'task_key': key, 'error': j.get('error_msg')})
        else:
            sc = atlas_api.fetch_result(j, hdrs, out_path(key),
                                        query_params=params_for(rows_by_key[key]), t0=t0)
            try:
                nrows = int(len(atlas_api.read_result(out_path(key))))
            except Exception as e:  # noqa: BLE001
                nrows = -1
                log({'event': 'parse_warning', 'task_key': key, 'msg': str(e)[:200]})
            log({'event': 'fetched', 'task_key': key, 'rows': nrows,
                 'bytes': sc['result_bytes'], 'sha256': sc['result_sha256'][:12],
                 'start': j.get('starttimestamp'), 'finish': j.get('finishtimestamp'),
                 'submitted': j.get('timestamp')})
        del inflight[key]
        save_inflight(inflight)
        n_done += 1
    return n_done


def status(hdrs, tasks, inflight):
    keys = [str(k) for k in tasks['task_key']]
    done = [k for k in keys if os.path.exists(sidecar_path(k))]
    failed = failed_keys()
    pending = [k for k in keys if k not in done and k not in inflight]
    print(f'tasks {len(keys)}  done {len(done)} (failed {len(failed)})  '
          f'inflight {len(inflight)}  pending {len(pending)}')
    try:
        qp = requests.get(f'{atlas_api.BASEURL}/queuepositions.json', headers=hdrs).json()
        print('server queue positions:', qp)
    except Exception as e:  # noqa: BLE001
        print('queuepositions failed:', e)
    if os.path.exists(LOG):
        walls = []
        for line in open(LOG):
            r = json.loads(line)
            if r.get('event') == 'fetched' and r.get('start') and r.get('finish'):
                a = datetime.fromisoformat(r['start'])
                b = datetime.fromisoformat(r['finish'])
                walls.append((b - a).total_seconds())
        if walls:
            print(f'server exec time per task: median {np.median(walls):.0f}s, '
                  f'max {np.max(walls):.0f}s over {len(walls)} tasks')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=None, help='submit at most N new tasks')
    ap.add_argument('--max-inflight', type=int, default=8)
    ap.add_argument('--poll', type=float, default=30.0)
    ap.add_argument('--only', nargs='*', help='restrict to these task keys')
    ap.add_argument('--status', action='store_true')
    ap.add_argument('--retry-failed', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    tasks = Table.read(TASKS)
    hdrs = atlas_api.headers()
    os.makedirs(LC, exist_ok=True)
    inflight = load_inflight()
    if args.status:
        status(hdrs, tasks, inflight)
        return
    if args.retry_failed:
        for k in failed_keys():
            os.remove(sidecar_path(k))
            log({'event': 'retry_failed', 'task_key': k})

    rows_by_key = {str(r['task_key']): r for r in tasks}
    keys = list(rows_by_key)
    if args.only:
        keys = [k for k in keys if k in set(args.only)]
    unfinished = {k for k in keys if not os.path.exists(sidecar_path(k))}
    adopt_orphans(hdrs, inflight, unfinished)
    pending = sorted((k for k in unfinished if k not in inflight),
                     key=lambda k: priority(rows_by_key[k]))
    if args.limit is not None:
        pending = pending[:args.limit]
    log({'event': 'start', 'pending': len(pending), 'inflight': len(inflight),
         'max_inflight': args.max_inflight, 'dry_run': args.dry_run})
    if args.dry_run:
        for k in pending[:20]:
            print(k, params_for(rows_by_key[k]))
        return

    while pending or inflight:
        while pending and len(inflight) < args.max_inflight:
            k = pending.pop(0)
            p = params_for(rows_by_key[k])
            j = atlas_api.submit(p, hdrs)
            inflight[k] = j['url']
            save_inflight(inflight)
            log({'event': 'submitted', 'task_key': k, 'url': j['url'], 'id': j.get('id')})
        if inflight:
            poll_once(hdrs, inflight, rows_by_key)
        if inflight:
            time.sleep(args.poll)
    log({'event': 'drain_complete'})


if __name__ == '__main__':
    main()
