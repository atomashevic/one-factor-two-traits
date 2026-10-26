"""The discontinuous switch of Section 5.2: two items loading .90 (block A) and four loading .70 (block B).

Uses the reduced problem of Online Resource 3, Sections 2.2-2.3: every positive root v of the
quartic P is a block-uniform stationary point with psi_A = N/(nu V), psi_B = N/(nu U), and the ML
discrepancy there is G(v) up to a constant that depends on phi only. At 50 digits the script
  1. finds the correlation phi_M at which the two local minima have equal discrepancy;
  2. reports the three positive roots there, the fitted loadings and the coefficient share S_A of
     the trait measured by the shorter block A on each minimum;
  3. confirms with the full 6-item discrepancy, not the reduced profile, that the two minima fit
     equally well at phi_M, and that the full Hessian in (loadings, uniquenesses) is positive
     definite at both minima and has one negative eigenvalue at the saddle;
  4. confirms that the global solution switches from the minimum favouring B to the one favouring
     A as phi increases through phi_M, and finds where the minimum favouring B disappears.
Numerical computation, reported in Online Resource 3, Section 2.6.
Run: python3 checks/unequal_jump_example.py (part of build.py --checks)
"""
from mpmath import mp, mpf, sqrt, log, polyroots, findroot, matrix, det, inverse, eig, diff

mp.dps = 50
failures = 0
mA, mB, lamA, lamB = 2, 4, mpf('0.9'), mpf('0.7')
LA, LB = lamA**2, lamB**2
aA, aB, bA, bB = 1 - LA, 1 - LB, 1 + (mA - 1)*LA, 1 + (mB - 1)*LB
nu, kA, kB = mA + mB - 1, (mA - 1)*(1 - LA), (mB - 1)*(1 - LB)


def report(name, ok):
    global failures
    failures += not ok
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)


def reduced(phi):
    """Positive roots of the quartic and, for each, (v, G(v) - const, alpha, beta, psi_A, psi_B)."""
    c = lamA*lamB*phi; cA, cB = mA*c, mB*c
    e0 = mA*bB + mB*bA - bA*bB - cA*cB
    U = lambda v: bA + cB*v
    V = lambda v: bB + cA/v
    N = lambda v: e0 + cB*kB*v + cA*kA/v
    P = lambda v: v*(v**2*V(v)*(nu*V(v) - N(v)) - U(v)*(nu*U(v) - N(v)))
    # P is a quartic; recover its coefficients from five exact evaluations
    xs = [mpf(k) for k in range(1, 6)]
    A = matrix([[x**j for j in range(4, -1, -1)] for x in xs])
    coef = list(inverse(A)*matrix([P(x) for x in xs]))
    roots = sorted(r.real for r in polyroots(coef, maxsteps=200, extraprec=200)
                   if abs(r.imag) < mpf('1e-30') and r.real > 0)
    out = []
    for v in roots:
        G = nu*log(N(v)) - (mB - 1)*log(U(v)) - (mA - 1)*log(V(v))
        psiA, psiB = N(v)/(nu*V(v)), N(v)/(nu*U(v))
        out.append((v, G, sqrt(1 - psiA), sqrt(1 - psiB), psiA, psiB))
    return out


def share_A(alpha, beta, phi):
    f, g = alpha/lamA, beta/lamB
    wA, wB = (f - phi*g)/(1 - phi**2), (g - phi*f)/(1 - phi**2)
    return wA/(wA + wB)


def population(phi):
    lam = [lamA]*mA + [lamB]*mB
    blk = [0]*mA + [1]*mB
    p = mA + mB
    return matrix([[1 if i == j else lam[i]*lam[j]*(1 if blk[i] == blk[j] else phi) for j in range(p)]
                   for i in range(p)])


def F_full(x, R):
    p = R.rows
    l, psi = x[:p], [mp.e**t for t in x[p:]]
    S = matrix([[l[i]*l[j] + (psi[i] if i == j else 0) for j in range(p)] for i in range(p)])
    M = inverse(S)*R
    return sum(M[i, i] for i in range(p)) - log(det(M)) - p


def full_point(alpha, beta, psiA, psiB):
    return [alpha]*mA + [beta]*mB + [log(psiA)]*mA + [log(psiB)]*mB


def gap(phi):
    sols = reduced(phi)
    assert len(sols) == 3, phi
    return sols[-1][1] - sols[0][1]      # largest root minus smallest root


# 1. the equal-fit correlation
lo, hi = mpf('0.32'), mpf('0.33')
report('the minimum favouring B is global at .32 and the one favouring A at .33', gap(lo) < 0 < gap(hi))
phi_M = findroot(gap, (lo, hi), solver='illinois')
sols = reduced(phi_M)
report(f'three positive roots at phi_M = {mp.nstr(phi_M, 12)}', len(sols) == 3)
report('phi_M rounds to .326, as stated in Section 5.2', round(float(phi_M), 3) == 0.326)

# 2. the two minima and the saddle
(v1, G1, a1, b1, pA1, pB1), (v2, G2, a2, b2, pA2, pB2), (v3, G3, a3, b3, pA3, pB3) = sols
S1, S3 = share_A(a1, b1, phi_M), share_A(a3, b3, phi_M)
for name, (v, G, al, be, pa, pb) in zip(('minimum 1', 'saddle', 'minimum 3'), sols):
    print(f'   {name}: v = {mp.nstr(v, 8)}, alpha = {mp.nstr(al, 6)}, beta = {mp.nstr(be, 6)}, '
          f'S_A = {mp.nstr(share_A(al, be, phi_M), 6)}')
report(f'the two minima have shares S_A = {mp.nstr(S1, 4)} and {mp.nstr(S3, 4)} '
       '(.92 and .10 at two decimals)',
       sorted([round(float(S1), 2), round(float(S3), 2)]) == [0.10, 0.92])

# 3. the full 6-item discrepancy
R = population(phi_M)
Fs = [F_full(full_point(al, be, pa, pb), R) for (_, _, al, be, pa, pb) in sols]
report(f'full discrepancy equal at the two minima ({mp.nstr(Fs[0], 10)}; difference '
       f'{mp.nstr(abs(Fs[0] - Fs[2]), 3)})', abs(Fs[0] - Fs[2]) < mpf('1e-30'))
report('profile G differences match full-discrepancy differences',
       abs((G2 - G1) - (Fs[1] - Fs[0])) < mpf('1e-30') and abs((G3 - G1) - (Fs[2] - Fs[0])) < mpf('1e-30'))


def hessian_signs(x0):
    n = len(x0)
    f = lambda *x: F_full(list(x), R)
    H = matrix(n, n)
    for i in range(n):
        for j in range(i, n):
            order = [0]*n; order[i] += 1; order[j] += 1
            H[i, j] = H[j, i] = diff(f, x0, tuple(order))
    ev = eig(H, left=False, right=False)
    return sum(1 for e in ev if e.real < 0), min(e.real for e in ev)


mp.dps = 30
neg = [hessian_signs(full_point(al, be, pa, pb)) for (_, _, al, be, pa, pb) in sols]
mp.dps = 50
report('full Hessian: minimum, saddle (one negative eigenvalue), minimum',
       [k for k, _ in neg] == [0, 1, 0] and neg[0][1] > 0 and neg[2][1] > 0)

# 4. the switch of the global solution, and the coexistence interval
below, above = reduced(phi_M - mpf('0.01')), reduced(phi_M + mpf('0.01'))
glob_below = min(below, key=lambda s: s[1]); glob_above = min(above, key=lambda s: s[1])
Sb = share_A(glob_below[2], glob_below[3], phi_M - mpf('0.01'))
Sa = share_A(glob_above[2], glob_above[3], phi_M + mpf('0.01'))
report(f'global solution favours B just below phi_M (S_A = {mp.nstr(Sb, 3)}) and A just above '
       f'(S_A = {mp.nstr(Sa, 3)})', Sb < 0.5 < Sa)


def count(phi):
    return len(reduced(phi))


def edge(lo, hi):
    """Bisect between a phi with three positive roots and one with one."""
    want = count(lo)
    for _ in range(60):
        mid = (lo + hi)/2
        lo, hi = (mid, hi) if count(mid) == want else (lo, mid)
    return (lo + hi)/2


grid = [mpf(k)/200 for k in range(1, 200)]
three = [phi for phi in grid if count(phi) == 3]
hi_fold = edge(three[-1], three[-1] + mpf(1)/200)
report(f'two minima coexist for every 0 < phi < {mp.nstr(hi_fold, 6)}, where the one favouring B '
       'meets the saddle', three[0] == grid[0] and three == grid[:len(three)] and phi_M < hi_fold
       and min(reduced(hi_fold - mpf('1e-6')), key=lambda s: s[1])[0] < 1)

print(f'{failures} failures')
raise SystemExit(1 if failures else 0)
