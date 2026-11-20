"""Reconstruct the default-rule CPM series and verify its physical witnesses.

Python >=3.10, standard library only. Lower bounds rely on archived full
QEDist receipt audits; this script does not rerun the exhaustive searches.
"""
from pathlib import Path
import hashlib
import json
import general_pp as g

ROOT = Path(__file__).resolve().parent


def read(name):
    return json.loads((ROOT / name).read_text())


def bits(row):
    assert len(set(row)) == len(row)
    return sum(1 << q for q in row)


def basis(rows):
    pivots = {}
    for value in rows:
        while value:
            p = value.bit_length() - 1
            if p not in pivots:
                pivots[p] = value
                break
            value ^= pivots[p]
    return pivots


def remainder(value, pivots):
    while value:
        p = value.bit_length() - 1
        if p not in pivots:
            return value
        value ^= pivots[p]
    return 0


for name, expected in read('SHA256SUMS.json').items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
frozen = read('frozen_rule.json')
assert hashlib.sha256((ROOT / 'general_pp.py').read_bytes()).hexdigest() == frozen['generator_sha256']
seed_b = [int(x) for line in (ROOT / 'seed276.cpm').read_text().splitlines()[2:] for x in line.split()]
assert json.loads(json.dumps(g.extract(4, 12, 23, [1] * 96, seed_b))) == frozen['rule']
audit = read('audit_summary.json')
assert audit['status'] == 'PASS' and audit['completed'] == audit['expected'] == 3
assert [(x['P'], x['n'], x['k'], x['d']) for x in audit['results']] == [(23, 276, 98, 14), (37, 444, 154, 16), (53, 636, 218, 20)]
for result in audit['results']:
    P = result['P']
    params, code = read(f'P{P}_parameters.json'), read(f'P{P}.json')
    a, b = g.construct(frozen['rule'], P)
    assert (a, b) == (params['a'], params['b']) and set(a) == {1}
    n = 12 * P
    matrices = []
    for side in (0, 1):
        rows = [0] * (4 * P)
        for i in range(4):
            for l in range(12):
                q = side * 48 + i * 12 + l
                for u in range(P):
                    rows[i * P + (a[q] * u + b[q]) % P] ^= 1 << (l * P + u)
        matrices.append(rows)
    assert matrices == [[bits(row) for row in code[key]] for key in ('hx', 'hz')]
    X, Z = matrices
    assert all((x & z).bit_count() % 2 == 0 for x in X for z in Z)
    bases = [basis(rows) for rows in matrices]
    assert [len(B) for B in bases] == result['ranks']
    assert n - sum(map(len, bases)) == result['k'] == params['k']
    for rows in matrices:
        assert all(row.bit_count() == 12 for row in rows)
        assert all(sum((row >> q) & 1 for row in rows) == 4 for q in range(n))
    for w in result['physical_witnesses']:
        side = 0 if w['side'] == 'X' else 1
        value = bits(w['support'])
        assert 0 < value < (1 << n) and value.bit_count() == w['weight']
        assert all((value & row).bit_count() % 2 == 0 for row in matrices[1 - side])
        assert remainder(value, bases[side])
    assert min(w['weight'] for w in result['physical_witnesses']) == result['d']
    excluded = {s['cutoff'] for s in result['stages'] if s['all_excluded']}
    assert set(range(1, result['d'])).issubset(excluded)
    print(f"PASS P={P}: [[{n},{result['k']},{result['d']}]], physical checks and archived-audit consistency")
