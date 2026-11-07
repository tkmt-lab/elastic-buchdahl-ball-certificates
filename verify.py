#!/usr/bin/env python3
"""Offline replay. Standard library only; all certificate computations are rerun."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent

def require(ok, message):
    if not ok:
        raise ValueError(message)

def integrity(root):
    names = set()
    for line in (root / 'MANIFEST.sha256').read_text().splitlines():
        digest, name = line.split('  ', 1)
        p = Path(name)
        require(not p.is_absolute() and '..' not in p.parts and name not in names,
                'unsafe or duplicate manifest path')
        names.add(name)
        target = root / p
        require(target.is_file() and not target.is_symlink(), 'missing or linked file: ' + name)
        require(hashlib.sha256(target.read_bytes()).hexdigest() == digest, 'hash mismatch: ' + name)
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and '.git' not in p.relative_to(root).parts}
    require(actual == names | {'MANIFEST.sha256'}, 'manifest inventory mismatch')
    return len(names)

def certificate_identity(root):
    """Check the numerical/checker snapshot identified by the article."""
    expected = {'spec.json', 'claims.json', 'verify.py', 'negative_controls.py'}
    expected.update(p.relative_to(root).as_posix()
                    for p in (root / 'certificates').rglob('*') if p.is_file())
    names = set()
    for line in (root / 'CERTIFICATE-SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ', 1)
        require(name in expected and name not in names, 'unexpected certificate identity path')
        names.add(name)
        target = root / name
        require(target.is_file() and not target.is_symlink() and
                hashlib.sha256(target.read_bytes()).hexdigest() == digest,
                'certificate identity mismatch: ' + name)
    require(names == expected, 'incomplete certificate identity inventory')
    return hashlib.sha256((root / 'CERTIFICATE-SHA256SUMS').read_bytes()).hexdigest()

def validate_b(data, expected):
    require(data['status'] == 'PASS' and data['tubes'] == 1592 and data['interval_bits'] == 192,
            'incomplete certificate B')
    for key in ('center_C1_verified', 'picard_all', 'endpoint_all', 'continuous_physics_all'):
        require(data[key] is True, 'B failed: ' + key)
    for key in ('root_x', 'R', 'compactness', 'compactness_minus_8_over_9', 'pressure_derivative_dx'):
        require(data[key] == expected[key], 'B exact result changed: ' + key)
    scale = 1 << 192
    bounds = [('R', '0.01258970', '0.01259407'), ('compactness', '0.91089951', '0.91180552')]
    for key, lo, hi in bounds:
        a, b = (Q(int(x), scale) for x in data[key])
        require(Q(lo) < a <= b < Q(hi), 'paper interval not enclosing ' + key)
    require(Q(int(data['compactness'][0]), scale) - Q(8, 9) > Q('0.02201062'),
            'compactness margin')
    require(all(int(hi) < 0 for lo, hi in data['pressure_derivative_dx']), 'nontransverse root')

def validate_a(data):
    require(data['row_count'] == data['requested_count'] == len(data['rows']) == 465,
            'incomplete certificate A')
    for key in ('all_C1', 'all_center_enclosed', 'all_endpoints_enclosed'):
        require(data[key] is True, 'A failed: ' + key)
    for i, row in enumerate(data['rows']):
        require(row['step'] == i, 'A step ordering')
        for key in ('picard', 'variational_tube', 'physics', 'pressure', 'C1',
                    'center_enclosed', 'endpoint_enclosed', 'C1_overlap'):
            require(row[key] is True, 'A row failed: ' + str(i) + ':' + key)
    for key in ('polynomial_match', 'own_center_in_source_initial_box'):
        require(data['center'][key] is True, 'A center failed')
    require(all(Q(x) < 1 for key in ('contraction_rows', 'selfmap_rows')
                for x in data['center'][key]), 'A center contraction failed')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--budget-seconds', type=float, default=3600)
    args = parser.parse_args()
    require(__debug__ and not sys.flags.optimize, 'run without -O or -OO')
    require(0 < args.budget_seconds <= 3600, 'invalid budget')
    start = time.monotonic()
    count = integrity(ROOT)
    identity = certificate_identity(ROOT)
    expected = json.loads((ROOT / 'certificates/own-tubes-checked.json').read_text())
    claims = json.loads((ROOT / 'claims.json').read_text())
    bound = next(c for c in claims['claims'] if c['id'] == 'surface-bounds')
    for key, result_key in [('exact_radius', 'R'), ('exact_compactness', 'compactness'), ('exact_root_x', 'root_x')]:
        require(bound[key] == expected[result_key], 'claims mismatch: ' + key)
    require(bound['radius_open_decimal'] == ['0.01258970', '0.01259407'] and
            bound['compactness_open_decimal'] == ['0.91089951', '0.91180552'] and
            bound['margin_above_8_over_9_strict_lower'] == '0.02201062', 'claims display mismatch')
    with tempfile.TemporaryDirectory(prefix='elastic-certificate-') as name:
        scratch = Path(name)
        shutil.copytree(ROOT / 'certificates', scratch / 'certificates')
        work = scratch / 'certificates'
        rows = [json.loads(s) for s in (work / 'stage2-inputs/pilot-steps.jsonl').read_text().splitlines()]
        require(len(rows) == 465 and Q(rows[0]['x_exact']) == Q(1415, 2097152), 'A initial x')
        require(all(r['step'] == i and Q(r['h_exact']) > 0 for i, r in enumerate(rows)), 'A step index/size')
        brows = [json.loads(s) for s in (work / 'stage1-tubes.jsonl').read_text().splitlines()]
        require(len(brows) == 1592, 'B tube count')
        require(all(int(tlo) < int(lo) <= int(hi) < int(thi)
                    for r in brows for (tlo, thi), (lo, hi) in zip(r['tube'], r['Y0'])),
                'B initial box not strictly inside tube')
        for script, result in [('check_own_tubes.py', 'own-tubes-checked.json'),
                               ('replay_producer.py', 'producer-replay.json')]:
            (work / result).unlink(missing_ok=True)
            code = 'import runpy,sys;sys.path.insert(0,sys.argv[1]);runpy.run_path(sys.argv[2],run_name="__main__")'
            remaining = args.budget_seconds - (time.monotonic() - start)
            require(remaining > 0, 'budget exhausted')
            print('Replaying ' + script, flush=True)
            env = dict(os.environ)
            env.pop('PYTHONOPTIMIZE', None)
            subprocess.run([sys.executable, '-I', '-B', '-c', code, str(work), str(work / script)],
                           check=True, timeout=remaining, env=env)
            data = json.loads((work / result).read_text())
            if script == 'check_own_tubes.py':
                validate_b(data, expected)
            else:
                validate_a(data)
    print(json.dumps({'status': 'PASS', 'manifest_files': count, 'certificate_A_tubes': 465,
                      'certificate_A_C1_entries': 4185, 'certificate_B_tubes': 1592,
                      'certificate_identity_sha256': identity,
                      'full_replay': True, 'seconds': round(time.monotonic() - start, 3)}, sort_keys=True))

if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        print('FAIL: ' + str(error), file=sys.stderr)
        raise SystemExit(1)
