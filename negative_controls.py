#!/usr/bin/env python3
"""Fast fail-closed controls; not a substitute for verify.py."""
import copy
import json
from pathlib import Path
import runpy
import tempfile
import sys
V = runpy.run_path(str(Path(__file__).with_name('verify.py')))
root = Path(__file__).resolve().parent

def rejected(fn, *args):
    try:
        fn(*args)
    except ValueError:
        return
    raise RuntimeError('negative control accepted')

a = json.loads((root / 'certificates/producer-replay.json').read_text())
b = json.loads((root / 'certificates/own-tubes-checked.json').read_text())
V['validate_a'](a); V['validate_b'](b, b)
c = copy.deepcopy(a); c['rows'] = c['rows'][:-1]; rejected(V['validate_a'], c)
c = copy.deepcopy(a); c['rows'][17]['C1'] = False; rejected(V['validate_a'], c)
c = copy.deepcopy(a); c['rows'][17]['endpoint_enclosed'] = False; rejected(V['validate_a'], c)
c = copy.deepcopy(a); c['rows'][17]['step'] = 18; rejected(V['validate_a'], c)
c = copy.deepcopy(b); c['compactness'][0] = '0'; rejected(V['validate_b'], c, b)
with tempfile.TemporaryDirectory() as name:
    p = Path(name); (p/'x').write_text('changed')
    (p/'MANIFEST.sha256').write_text('0'*64+'  x\n'); rejected(V['integrity'], p)
    (p/'MANIFEST.sha256').write_text('0'*64+'  ../x\n'); rejected(V['integrity'], p)
# A changed flow derivative must fail actual equation-based replay, not merely a hash.
sys.path.insert(0, str(root/'certificates'))
R = runpy.run_path(str(root/'certificates/replay_producer.py'))
row = json.loads((root/'certificates/stage2-inputs/pilot-steps.jsonl').read_text().splitlines()[0])
row['flow_derivative'][0][0] = {'lo':'0', 'hi':'0'}
changed = R['replay'](row, True)
if changed['C1'] is not False:
    raise RuntimeError('changed derivative accepted by arithmetic replay')
c = copy.deepcopy(a); c['rows'][0] = changed; rejected(V['validate_a'], c)
# The article's frozen identity must also reject a changed numerical input,
# even when someone has regenerated the complete packet manifest.
with tempfile.TemporaryDirectory() as name:
    import shutil
    p = Path(name)
    for line in (root/'CERTIFICATE-SHA256SUMS').read_text().splitlines():
        _, filename = line.split('  ', 1)
        dest = p/filename; dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root/filename, dest)
    shutil.copyfile(root/'CERTIFICATE-SHA256SUMS', p/'CERTIFICATE-SHA256SUMS')
    V['certificate_identity'](p)
    (p/'spec.json').write_text('{}\n')
    import hashlib
    (p/'MANIFEST.sha256').write_text(''.join(
        hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.relative_to(p).as_posix()+'\n'
        for f in sorted(p.rglob('*')) if f.is_file() and f.name!='MANIFEST.sha256'))
    V['integrity'](p)
    rejected(V['certificate_identity'], p)
print(json.dumps({'status':'PASS', 'negative_controls':9, 'includes_arithmetic_replay_mutation':True,
                  'includes_refreshed_manifest_with_changed_scientific_input':True}))
