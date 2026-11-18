"""Checks for Appendix C (proofs of Theorem 2 and Proposition 1).

Symbolic: the relaxed objective f, its partial derivatives, the boundary
stationarity step, the cone form h, the edge quadratic and the multiplier.
Numerical: the stated minimizers and minimum values against a multistart
minimization of the original bounded problem, for weight ratios below, at and
above (m-1)/m, and the attainment candidate against the escape limit.

These checks verify the identities and cases they compute. They do not
establish the theorem.

Run: python3 checks/least_squares_relaxation.py   (NumPy, SciPy, SymPy)
"""
import itertools
import numpy as np
import sympy as sp
from scipy.optimize import minimize

fails = 0


def report(name, ok):
    global fails
    fails += not ok
    print(('PASS' if ok else 'FAIL') + ': ' + name)


# ---------------------------------------------------------------- symbolic
xi, eta, QA, QB, L, c, w, m = sp.symbols('xi eta Q_A Q_B L c omega m', positive=True)
K = m*(m - 1)
C = K/2
SA = sp.sqrt(QA + K*xi)
SB = sp.sqrt(QB + K*eta)
f = C*((L - xi)**2 + (L - eta)**2) + w*(m**2*c**2 - 2*c*SA*SB + QA*QB)

report('df/dQ_A = omega (Q_B - c S_B/S_A)',
       sp.simplify(sp.diff(f, QA) - w*(QB - c*SB/SA)) == 0)
report('df/dxi = m(m-1)(xi - L - omega c S_B/S_A)',
       sp.simplify(sp.diff(f, xi) - K*(xi - L - w*c*SB/SA)) == 0)

# derivative along the boundary Q_B = m eta
fb = f.subs(QB, m*eta)
dfb = sp.diff(fb, eta)
SAs = sp.Symbol('S_A', positive=True)
target = K*(eta - L) + w*m*QA - w*m*c*SAs/sp.sqrt(eta)
report('boundary stationarity equation',
       sp.simplify(dfb.subs(QA + K*xi, SAs**2) - target.subs(SAs, SA).subs(QA + K*xi, SAs**2)) == 0
       or sp.simplify(dfb - target.subs(SAs, SA)) == 0)
# with S_A = c/sqrt(eta): equation becomes eta = L + omega m xi
eq = (K*(eta - L) + w*m*QA - w*m*(QA + K*xi))
report('boundary equation reduces to eta = L + omega m xi',
       sp.simplify(eq - K*(eta - L - w*m*xi)) == 0)

# cone form for block-uniform loadings
a, b, phi = sp.symbols('alpha beta varphi', positive=True)
unif = f.subs({xi: a**2, eta: b**2, QA: m*a**2, QB: m*b**2, c: L*phi})
u, v = a**2 + b**2, a*b
d = w*m - m + 1
h = C*(u - L)**2 + m*d*v**2 - 2*w*m**2*L*phi*v + C*L**2 + w*m**2*L**2*phi**2
report('f at tight points equals h(u,v)', sp.simplify(sp.expand(unif - h)) == 0)

U, V, t = sp.symbols('u v t', real=True)
hh = C*(U - L)**2 + m*d*V**2 - 2*w*m**2*L*phi*V + C*L**2 + w*m**2*L**2*phi**2
edge = hh.subs({U: 2*t, V: t})
report('edge quadratic', sp.simplify(edge - (K*(L - t)**2 + w*m**2*(L*phi - t)**2)) == 0)
ts = L*(2*C + w*m**2*phi)/(2*C + w*m**2)
report('t_s minimizes the edge quadratic', sp.simplify(sp.diff(edge, t).subs(t, ts)) == 0)
varpi = 2*C*(2*ts - L)
report('multiplier closed form',
       sp.simplify(varpi - K*L*(2*w*m*phi - d)/(m - 1 + w*m)) == 0)
grad = [sp.diff(hh, U).subs({U: 2*ts, V: ts}), sp.diff(hh, V).subs({U: 2*ts, V: ts})]
report('gradient at (2 t_s, t_s) is varpi (1, -2)',
       sp.simplify(grad[0] - varpi) == 0 and sp.simplify(grad[1] + 2*varpi) == 0)
phiw = sp.Rational(1, 2) - (m - 1)/(2*m*w)
report('varpi = 0 and t_s = L/2 at the threshold',
       sp.simplify(varpi.subs(phi, phiw)) == 0 and sp.simplify(ts.subs(phi, phiw) - L/2) == 0)
vstar = w*m*L*phi/d
report('asymmetric minimum value',
       sp.simplify(hh.subs({U: L, V: vstar}) - (C*L**2 - w*m**2*(m - 1)*L**2*phi**2/d)) == 0)
report('symmetric minimum value',
       sp.simplify(edge.subs(t, ts) - w*m**2*(m - 1)*L**2*(1 - phi)**2/(m - 1 + w*m)) == 0)
report('(alpha - beta)^2 on the asymmetric branch',
       sp.simplify(L - 2*vstar - 2*w*m*L*(phiw - phi)/d) == 0)
cand = f.subs({xi: L, eta: L*phi**2, QA: m*L, QB: m*L*phi**2, c: L*phi})
report('attainment candidate value C_m L^2 (1 - phi^2)^2',
       sp.simplify(cand - C*L**2*(1 - phi**2)**2) == 0)
T = sp.Symbol('T', positive=True)
esc = T - 2*c*sp.sqrt(m*T) + m**2*c**2
report('escape limit inf_T = m(m-1)c^2',
       sp.simplify(esc.subs(T, m*c**2) - K*c**2) == 0
       and sp.simplify(sp.diff(esc, T).subs(T, m*c**2)) == 0)


# ---------------------------------------------------------------- numerical
def G(l, mm, Lv, ph, om):
    A, B = l[:mm], l[mm:]
    within = 0.0
    for blk in (A, B):
        P = np.outer(blk, blk)
        iu = np.triu_indices(mm, 1)
        within += np.sum((Lv - P[iu])**2)
    return within + om*np.sum((Lv*ph - np.outer(A, B))**2)


def brute(mm, Lv, ph, om, rng, nstart=60):
    best = None
    for _ in range(nstart):
        x0 = rng.uniform(-1, 1, 2*mm)
        r = minimize(G, x0, args=(mm, Lv, ph, om), method='L-BFGS-B',
                     bounds=[(-1, 1)]*(2*mm), options={'ftol': 1e-15, 'gtol': 1e-11, 'maxiter': 5000})
        if best is None or r.fun < best.fun:
            best = r
    return best


def predicted(mm, Lv, ph, om):
    dd = om*mm - mm + 1
    Cm = mm*(mm - 1)/2
    if dd > 0 and ph < 0.5 - (mm - 1)/(2*mm*om):
        uu, vv = Lv, om*mm*Lv*ph/dd
        al = (np.sqrt(uu + 2*vv) + np.sqrt(uu - 2*vv))/2
        be = (np.sqrt(uu + 2*vv) - np.sqrt(uu - 2*vv))/2
        return Cm*Lv**2 - om*mm**2*(mm - 1)*Lv**2*ph**2/dd, (al, be)
    tsv = Lv*(mm - 1 + om*mm*ph)/(mm - 1 + om*mm)
    return om*mm**2*(mm - 1)*Lv**2*(1 - ph)**2/(mm - 1 + om*mm), (np.sqrt(tsv), np.sqrt(tsv))


rng = np.random.default_rng(20261117)
worst_val = worst_load = 0.0
ncase = 0
for mm, lam, om, ph in itertools.product((2, 3, 4), (.5, .7), (.3, .5, 2/3, .75, 1., 2.), (.02, .08, .2, .45, .8)):
    Lv = lam**2
    val, (al, be) = predicted(mm, Lv, ph, om)
    r = brute(mm, Lv, ph, om, rng)
    x = np.abs(r.x)
    blocks = sorted([x[:mm], x[mm:]], key=lambda z: -z.mean())
    worst_val = max(worst_val, abs(r.fun - val))
    worst_load = max(worst_load, np.max(np.abs(blocks[0] - al)), np.max(np.abs(blocks[1] - be)))
    same_sign = np.all(r.x > 0) or np.all(r.x < 0)
    if not same_sign:
        worst_load = np.inf
    # attainment candidate below the escape limit
    Cm = mm*(mm - 1)/2
    assert Cm*Lv**2*(1 - ph**2)**2 < Cm*Lv**2 + om*mm*(mm - 1)*(Lv*ph)**2
    ncase += 1
report(f'{ncase} cases: minimum value matches multistart search (max diff {worst_val:.1e})', worst_val < 1e-9)
report(f'{ncase} cases: common sign, block-uniform loadings match (max diff {worst_load:.1e})', worst_load < 1e-4)

# zero correlation, d > 0, m = 2: the retained pair (t, L/t) attains C_m L^2
Lv = .49
vals = [G(np.array([tt, Lv/tt, 0, 0]), 2, Lv, 0., 1.) for tt in np.linspace(Lv, 1, 11)]
report('phi = 0, m = 2: retained pairs (t, L/t) all give C_m L^2', np.allclose(vals, Lv**2))
r = brute(2, Lv, 0., 1., rng)
report('phi = 0, m = 2: multistart minimum equals C_m L^2', abs(r.fun - Lv**2) < 1e-10)
r = brute(3, Lv, 0., 1., rng)
x = np.sort(np.abs(r.x))
report('phi = 0, m = 3: one block at lambda, the other at zero',
       abs(r.fun - 3*Lv**2) < 1e-10 and np.allclose(x[3:], .7, atol=1e-4) and np.allclose(x[:3], 0, atol=1e-4))

print('ALL PASS' if not fails else f'{fails} FAILED')
raise SystemExit(1 if fails else 0)
