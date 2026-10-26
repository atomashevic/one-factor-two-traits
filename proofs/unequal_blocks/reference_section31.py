"""Reference numerics for Section 3.1 of 'When one factor represents two traits'.
Reproduces every number the paper states in Section 3.1. Written by the main author's
assistant; the proof agent must NOT reuse this code for its own verification, only
compare against its outputs.  Run: python3 reference_section31.py  (about 1 to 3 minutes)."""
import numpy as np, json
from scipy.optimize import minimize, brentq
rng = np.random.default_rng(2026)

def pop(lA, lB, phi):
    """Two-trait population correlation matrix; lA, lB are lists of loadings (block A on trait 1, block B on trait 2)."""
    lA, lB = np.asarray(lA, float), np.asarray(lB, float); mA, mB = len(lA), len(lB); p = mA + mB
    Lam = np.zeros((p, 2)); Lam[:mA, 0] = lA; Lam[mA:, 1] = lB
    S = Lam @ np.array([[1, phi], [phi, 1]]) @ Lam.T; np.fill_diagonal(S, 1); return S

def F(th, S):  # ML discrepancy, parameters (l, log psi)
    p = S.shape[0]; l = th[:p]; psi = np.exp(th[p:]); Sig = np.outer(l, l) + np.diag(psi)
    A = np.linalg.solve(Sig, S); return np.trace(A) - np.linalg.slogdet(A)[1] - p

def grad(th, S):
    p = S.shape[0]; l = th[:p]; psi = np.exp(th[p:]); Sig = np.outer(l, l) + np.diag(psi)
    Si = np.linalg.inv(Sig); D = Si - Si @ S @ Si; return np.r_[2 * D @ l, np.diag(D) * psi]

def polish(S, x0):
    return minimize(F, x0, args=(S,), jac=grad, method='L-BFGS-B',
                    options={'maxiter': 20000, 'ftol': 1e-16, 'gtol': 1e-13})

def best(S, starts=20):
    p = S.shape[0]; out = None
    for _ in range(starts):
        x0 = np.r_[rng.uniform(-0.1, 1, p), np.log(rng.uniform(.2, 1, p))]
        try: r = polish(S, x0)
        except np.linalg.LinAlgError: continue
        if out is None or r.fun < out.fun: out = r
    return polish(S, out.x)

def share(lbarA, lbarB, lamA, lamB, phi, which=None):
    """Composite weights w solve lbarA = lamA (wA + phi wB), lbarB = lamB (phi wA + wB).
    Returns the share of the followed trait max(w)/sum(w), or of trait `which` (0 = A, 1 = B)."""
    w = np.abs(np.linalg.solve(np.array([[1, phi], [phi, 1]]), [lbarA / lamA, lbarB / lamB]))
    return (w.max() if which is None else w[which]) / w.sum()

res = {}
# 1. Lemma 1: random block-specific designs, full parameter space
worst = 0
for _ in range(60):
    mA, mB = rng.integers(2, 6, 2); la, lb = rng.uniform(.3, .9, 2); phi = rng.uniform(.01, .95)
    b = best(pop([la] * mA, [lb] * mB, phi), 25); l = np.abs(b.x[:mA + mB]); ps = np.exp(b.x[mA + mB:])
    worst = max(worst, l[:mA].std(), l[mA:].std(), ps[:mA].std(), ps[mA:].std())
res['lemma1_max_within_block_spread_60_designs'] = worst

# 2. Mirrored designs: threshold where the antisymmetric Hessian at the swap-symmetric solution is singular
def sym_solution(lam, phi):
    k = len(lam); S = pop(lam, lam, phi); full = lambda t: np.r_[t[:k], t[:k], t[k:], t[k:]]
    f = lambda t: F(full(t), S)
    gr = lambda t: (lambda G: np.r_[G[:k] + G[k:2 * k], G[2 * k:3 * k] + G[3 * k:]])(grad(full(t), S))
    t0 = np.r_[np.array(lam) * .8, np.log(1 - np.array(lam) ** 2 * .64)]
    return full(minimize(f, t0, jac=gr, method='BFGS', options={'gtol': 1e-13}).x), S

def anti_min_eig(lam, phi, h=1e-6):
    th, S = sym_solution(lam, phi); k = len(lam); n = len(th)
    H = np.array([(grad(th + h * e, S) - grad(th - h * e, S)) / (2 * h) for e in np.eye(n)]); H = (H + H.T) / 2
    B = []
    for i in range(k):
        for off in (0, 2 * k):
            v = np.zeros(n); v[off + i] = 1; v[off + k + i] = -1; B.append(v / np.sqrt(2))
    B = np.array(B).T; return np.linalg.eigvalsh(B.T @ H @ B).min()

phistar = lambda L, m: (1 + (m - 1) * L) / (m * (2 - L))
for lam in ([.7] * 4, [.6, .65, .75, .8], [.5, .6, .7, .8]):
    lam = np.array(lam); c = lam ** 2 / (1 - lam ** 2)
    res[f'mirrored_threshold_{list(lam)}'] = brentq(lambda f: anti_min_eig(lam, f), .3, .6, xtol=1e-8)
    res[f'mean_c_prediction_{list(lam)}'] = phistar(c.mean() / (1 + c.mean()), 4)

# 3. Figure 2 values (global minimum), equal block sizes m = 4
for lab, lA, lB in [('mild_.68_.72', .68, .72), ('strong_.60_.80', .6, .8)]:
    for phi in (0.40894, 0.9):
        b = best(pop([lA] * 4, [lB] * 4, phi)); l = np.abs(b.x[:8])
        res[f'share_{lab}_phi{phi}'] = share(l[:4].mean(), l[4:].mean(), lA, lB, phi)

# 4. Weaker-trait local minimum: continuation from the weaker-block solution until it disappears
def fold_point(lA, lB, m=4, step=5e-4):
    x = np.r_[[lA] * m, [0.05] * m, np.log([1 - lA ** 2] * m), np.log([.99] * m)]; last = None
    for phi in np.arange(step, .42, step):
        r = polish(pop([lA] * m, [lB] * m, phi), x); l = np.abs(r.x[:2 * m])
        if l[:m].mean() < l[m:].mean(): return last
        x, last = r.x, phi
    return last
r_loc = polish(pop([.68] * 4, [.72] * 4, .05), np.r_[[.68] * 4, [.05] * 4, np.log([1 - .68 ** 2] * 4), np.log([.99] * 4)])
r_glb = best(pop([.68] * 4, [.72] * 4, .05))
res['mild_phi.05_F_local_vs_global'] = [r_loc.fun, r_glb.fun]
ps = phistar(.49, 4)
for e in (.001, .0025, .005, .01, .02, .04):
    res[f'fold_eps{e}'] = fold_point(.7 - e, .7 + e)

# 5. Tilt at phi* for loadings .7 -/+ eps
for e in (.0005, .001, .002, .004, .008, .016):
    b = best(pop([.7 - e] * 4, [.7 + e] * 4, ps), 30); l = np.abs(b.x[:8])
    res[f'tilt_at_phistar_eps{e}'] = share(l[:4].mean(), l[4:].mean(), .7 - e, .7 + e, ps) - .5

# 6. Unequal block sizes: 2 items at .85 vs 3 items at .75; share of block A (the shorter, more reliable block)
for phi in (.01, .5, .75, .775, .9):
    b = best(pop([.85] * 2, [.75] * 3, phi), 25); l = np.abs(b.x[:5])
    res[f'sizes_2x.85_3x.75_shareA_phi{phi}'] = share(l[:2].mean(), l[2:].mean(), .85, .75, phi, which=0)

for k, v in res.items(): print(k, v)
json.dump({k: (float(v) if np.isscalar(v) else [float(t) for t in v]) for k, v in res.items()},
          open('reference_section31.json', 'w'), indent=1)
