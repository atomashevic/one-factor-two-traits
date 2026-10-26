#!/usr/bin/env python3
"""Recompute Table 2 (tab:2) from Proposition 2 and the Theorem 1 threshold.

Each entry is the smallest secondary-to-primary slope ratio r = s/p at which the
Gaussian one-factor ML approximation of two uncorrelated traits is a composite
(equal block loadings).  With trait correlation rho = 0, residual variance
tau^2 = 1 (probit latent responses) and primary slope a, Proposition 2 gives

    sigma_S^2 = a^2 (1 + r^2),  lambda_eff^2 = sigma_S^2 / (1 + sigma_S^2),
    phi_eff   = 2 r / (1 + r^2),

and Theorem 1 gives the ML threshold

    phi*(lambda, m) = (1 + (m - 1) lambda^2) / (m (2 - lambda^2)),

with limit lambda^2 / (2 - lambda^2) as m -> infinity ("many items").  The
threshold ratio solves phi_eff(r) = phi*(lambda_eff(r), m), found by bisection (standard library only, so the default build needs no SciPy).
The table reports these exact values rounded UP to two decimals, so that each
printed ratio suffices and the ratio .01 below it does not.

The script prints the exact ratios, the rounded-up values and the printed cells
parsed from manuscript/sections/05c_crossloadings.tex, and exits nonzero on any
mismatch (including the text range "from .08 to .22").
"""
from pathlib import Path
import math
import re
import sys



def linspace(lo, hi, n):
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def bisect(f, lo, hi, args=(), tol=1e-14, maxiter=200):
    flo = f(lo, *args)
    for _ in range(maxiter):
        mid = (lo + hi) / 2
        fm = f(mid, *args)
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'manuscript' / 'sections' / '05c_crossloadings.tex'
SLOPES = (0.6, 0.8, 1.0)
ITEMS = (5, 20, math.inf)  # items per trait; inf = "many items"


def phi_star(lam2, m):
    if math.isinf(m):
        return lam2 / (2.0 - lam2)
    return (1.0 + (m - 1.0) * lam2) / (m * (2.0 - lam2))


def margin(r, a, m):
    """phi_eff - phi*; nonnegative exactly when the ML solution is a composite."""
    sig2 = a * a * (1.0 + r * r)
    lam2 = sig2 / (1.0 + sig2)
    return 2.0 * r / (1.0 + r * r) - phi_star(lam2, m)


def threshold(a, m):
    # margin(0) = -phi* < 0 and margin(1) = 1 - phi* > 0 for 0 < lambda^2 < 1.
    grid = linspace(0.0, 1.0, 20001)
    signs = [(v > 0) - (v < 0) for v in (margin(r, a, m) for r in grid)]
    changes = sum(1 for x, y in zip(signs, signs[1:]) if x != y)
    if changes != 1:
        raise AssertionError(f'expected one sign change on [0, 1], found {changes} (a={a}, m={m})')
    return bisect(margin, 0.0, 1.0, args=(a, m))


def round_up(x):
    return math.ceil(round(100.0 * x, 9)) / 100.0


def fmt(x):
    s = f'{x:.2f}'
    return s[1:] if s.startswith('0.') else s


def printed_table():
    tex = SOURCE.read_text()
    body = tex.split(r'\label{tab:2}', 1)[1]
    header = body.split(r'\midrule', 1)[0]
    rows = body.split(r'\midrule', 1)[1].split(r'\bottomrule', 1)[0]
    cells = {}
    for line in rows.splitlines():
        if '&' not in line:
            continue
        parts = [c.strip() for c in line.replace(r'\\', '').split('&')]
        cells[float(parts[0])] = parts[1:]
    return tex, header, cells


def main():
    tex, header, cells = printed_table()
    failures = []
    for token in ('5 items per trait', '20 items per trait', 'Many items'):
        if token not in header:
            failures.append(f'column header "{token}" not found')
    print('Table 2: smallest slope ratio r with phi_eff >= phi*(lambda_eff, m), rho = 0, tau^2 = 1')
    print(f'{"slope":>6} | ' + ' | '.join(f'{("m=" + str(m)) if not math.isinf(m) else "many":>26}' for m in ITEMS))
    rounded_all = []
    for a in SLOPES:
        out = []
        printed = cells.get(a)
        if printed is None or len(printed) != len(ITEMS):
            failures.append(f'row for slope {a} missing or malformed: {printed}')
            printed = [''] * len(ITEMS)
        for m, cell in zip(ITEMS, printed):
            r = threshold(a, m)
            up = round_up(r)
            rounded_all.append(up)
            # The rounded-up ratio suffices and the ratio .01 below does not.
            if margin(up, a, m) < 0 or margin(up - 0.01, a, m) >= 0:
                failures.append(f'rounding check failed at a={a}, m={m}')
            # Once sufficient, larger ratios up to 1 remain sufficient.
            if min(margin(x, a, m) for x in linspace(r + 1e-9, 1.0, 2001)) < 0:
                failures.append(f'margin not nonnegative above threshold at a={a}, m={m}')
            ok = cell == fmt(up)
            if not ok:
                failures.append(f'slope {a}, m={m}: printed {cell!r}, expected {fmt(up)!r} (exact {r:.6f})')
            out.append(f'{r:.6f} -> {fmt(up)} [{cell or "?"}]{"" if ok else " MISMATCH"}')
        print(f'{a:>6} | ' + ' | '.join(f'{c:>26}' for c in out))
    lo, hi = fmt(min(rounded_all)), fmt(max(rounded_all))
    flat = re.sub(r'\s+', ' ', tex)
    if f'from {lo} to {hi}' not in flat:
        failures.append(f'text range "from {lo} to {hi}" not found in {SOURCE.name}')
    else:
        print(f'Text range "from {lo} to {hi}" present.')
    if 'rounded up' not in flat.split(r'\label{tab:2}', 1)[0].rsplit(r'\begin{table}', 1)[-1]:
        failures.append('Table 2 caption does not state that values are rounded up')
    if failures:
        print('FAIL:')
        for f in failures:
            print('  ' + f)
        return 1
    print('PASS: all 9 Table 2 cells equal the exact thresholds rounded up to two decimals')
    return 0


if __name__ == '__main__':
    sys.exit(main())
