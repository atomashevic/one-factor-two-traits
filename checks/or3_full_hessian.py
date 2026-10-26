"""Full-Hessian check of Online Resource 3, Section 2.4 (contrast sectors and inertia).

At every positive root of the quartic (S2), in designs with one and with three positive roots,
the script builds the analytic Hessian of the ML discrepancy F(l, psi) in all 2p loadings and
uniquenesses, at 30 digits, and checks:
  1. the point is stationary;
  2. in a basis of block means and within-block contrasts, the Hessian has no cross terms
     between the block-mean sector and the contrast sectors, or between the A and B contrast
     sectors, and on each contrast sector it equals I (x) H_A (or H_B), with H_A the 2x2 matrix
     of (S7), for a random orthonormal contrast basis;
  3. the full Hessian is positive definite at a scalar minimum and has exactly one negative
     eigenvalue, and no zero eigenvalue, at a scalar maximum (the saddle between two minima).
Numerical corroboration only; the argument is in Online Resource 3, Section 2.4.
Run: python3 checks/or3_full_hessian.py (part of build.py --checks)
"""
import random
from mpmath import mp, mpf, sqrt, matrix, inverse, polyroots, eigsy, qr

mp.dps = 30
random.seed(20261026)
failures = 0


def report(name, ok):
    global failures
    failures += not ok
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)


def roots(mA, mB, lA, lB, phi):
    """Positive roots v of (S2) with (P'(v), alpha, beta, psi_A, psi_B, tau)."""
    LA, LB = lA**2, lB**2
    aA, aB, bA, bB = 1 - LA, 1 - LB, 1 + (mA - 1)*LA, 1 + (mB - 1)*LB
    nu, kA, kB = mA + mB - 1, (mA - 1)*aA, (mB - 1)*aB
    c = lA*lB*phi; cA, cB = mA*c, mB*c
    e0 = mA*bB + mB*bA - bA*bB - cA*cB
    U = lambda v: bA + cB*v
    V = lambda v: bB + cA/v
    N = lambda v: e0 + cB*kB*v + cA*kA/v
    P = lambda v: v*(v**2*V(v)*(nu*V(v) - N(v)) - U(v)*(nu*U(v) - N(v)))
    xs = [mpf(k) for k in range(1, 6)]
    coef = list(inverse(matrix([[x**j for j in range(4, -1, -1)] for x in xs]))*matrix([P(x) for x in xs]))
    out = []
    for r in polyroots(coef, maxsteps=200, extraprec=200):
        if abs(r.imag) < mpf('1e-20') and r.real > 0:
            v = r.real
            dP = sum((4 - j)*coef[j]*v**(3 - j) for j in range(4))
            psiA, psiB = N(v)/(nu*V(v)), N(v)/(nu*U(v))
            out.append((v, dP, sqrt(1 - psiA), sqrt(1 - psiB), psiA, psiB, nu*U(v)*V(v)/N(v)))
    return sorted(out)


def hessian(R, l, psi):
    p = len(l)
    S = matrix([[l[i]*l[j] + (psi[i] if i == j else 0) for j in range(p)] for i in range(p)])
    K = inverse(S); KRK = K*R*K; G = K - KRK
    grad = max([abs(2*(G*matrix(l))[i]) for i in range(p)] + [abs(G[i, i]) for i in range(p)])
    dirs = []
    for k in range(2*p):
        dl = [mpf(0)]*p
        D = matrix(p, p)
        if k < p:
            dl[k] = mpf(1)
            for j in range(p):
                D[k, j] += l[j]; D[j, k] += l[j]
        else:
            D[k - p, k - p] = mpf(1)
        dirs.append((dl, K*D))
    H = matrix(2*p, 2*p)
    for x in range(2*p):
        for y in range(x, 2*p):
            (dlx, Ax), (dly, Ay) = dirs[x], dirs[y]
            AxAy, AyAx = Ax*Ay, Ay*Ax
            val = sum(-AxAy[i, i] for i in range(p))
            val += sum((AxAy*K*R)[i, i] + (AyAx*K*R)[i, i] for i in range(p))
            val += sum(G[i, j]*(dlx[i]*dly[j] + dly[i]*dlx[j]) for i in range(p) for j in range(p))
            H[x, y] = H[y, x] = val
    return H, grad


def contrast_basis(m, offset, p):
    """1/sqrt(m) on the block and a random orthonormal basis of its contrasts, as columns in R^p."""
    cols = [[mpf(1)/sqrt(m) if offset <= i < offset + m else mpf(0) for i in range(p)]]
    for _ in range(m - 1):
        cols.append([mpf(random.gauss(0, 1)) if offset <= i < offset + m else mpf(0) for i in range(p)])
    Q, _ = qr(matrix(cols).T)
    return [[Q[i, k] for i in range(p)] for k in range(m)]


designs = [(2, 4, '0.9', '0.7', '0.3'), (4, 4, '0.68', '0.72', '0.1'), (2, 3, '0.85', '0.75', '0.2'),
           (3, 3, '0.6', '0.8', '0.01'), (3, 5, '0.6', '0.8', '0.5'), (2, 2, '0.3', '0.8', '0.05')]
for _ in range(4):
    designs.append((random.randint(2, 5), random.randint(2, 5), f'{random.uniform(.3, .9):.3f}',
                    f'{random.uniform(.3, .9):.3f}', f'{random.uniform(.05, .95):.3f}'))

n_points = n_saddles = 0
worst_grad = worst_cross = worst_form = mpf(0)
inertia_ok = True
for mA, mB, lA, lB, phi in designs:
    lA, lB, phi = mpf(lA), mpf(lB), mpf(phi)
    p = mA + mB
    lam = [lA]*mA + [lB]*mB; blk = [0]*mA + [1]*mB
    R = matrix([[1 if i == j else lam[i]*lam[j]*(1 if blk[i] == blk[j] else phi) for j in range(p)]
                for i in range(p)])
    for v, dP, al, be, pA, pB, tau in roots(mA, mB, lA, lB, phi):
        H, grad = hessian(R, [al]*mA + [be]*mB, [pA]*mA + [pB]*mB)
        worst_grad = max(worst_grad, grad)
        A, B = contrast_basis(mA, 0, p), contrast_basis(mB, mA, p)
        basis = [A[0], B[0]] + A[1:] + B[1:]          # block means first, then contrasts
        T = matrix(2*p, 2*p)
        for k, col in enumerate(basis):
            for i in range(p):
                T[i, 2*k] = col[i]; T[p + i, 2*k + 1] = col[i]   # (loading, uniqueness) pair per vector
        Hb = T.T*H*T
        sector = [0, 0] + [1 + k for k in range(p - 2)]  # sector index of each basis vector
        for x in range(2*p):
            for y in range(2*p):
                if sector[x//2] != sector[y//2]:
                    worst_cross = max(worst_cross, abs(Hb[x, y]))
        for k in range(2, p):
            in_A = k < mA + 1
            psi, alpha, a = (pA, al, 1 - lA**2) if in_A else (pB, be, 1 - lB**2)
            S7 = [[2/psi - 2*a/(tau*psi**2), 2*a*alpha/(tau*psi**3)],
                  [2*a*alpha/(tau*psi**3), 1/psi**2 - 2*a*alpha**2/(tau*psi**4)]]
            worst_form = max(worst_form, max(abs(Hb[2*k + i, 2*k + j] - S7[i][j]) for i in range(2) for j in range(2)))
        ev = eigsy(H)[0]
        negatives = sum(1 for e in ev if e < 0)
        expected = 1 if dP > 0 else 0
        inertia_ok &= negatives == expected and min(abs(e) for e in ev) > mpf('1e-12')
        n_points += 1; n_saddles += expected

report(f'{n_points} positive roots in {len(designs)} designs are stationary (max gradient {mp.nstr(worst_grad, 3)})',
       worst_grad < mpf('1e-20'))
report(f'no cross terms between the block-mean and contrast sectors (max {mp.nstr(worst_cross, 3)})',
       worst_cross < mpf('1e-20'))
report(f'each contrast sector equals the 2x2 matrix of (S7) (max difference {mp.nstr(worst_form, 3)})',
       worst_form < mpf('1e-20'))
report(f'full Hessian inertia: positive definite at minima, one negative eigenvalue at the '
       f'{n_saddles} saddles, no zero eigenvalue', inertia_ok and n_saddles > 0)
print(f'{failures} failures')
raise SystemExit(1 if failures else 0)
