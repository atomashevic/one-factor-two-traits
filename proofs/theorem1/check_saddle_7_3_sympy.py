import sympy as sp
s, phi, m, eps = sp.symbols('s phi m epsilon', positive=True)
n = m - 1; N = 2*m - 1
z = s**2                              # z = l_s^2, the symmetric squared loading
L = z*N/(n + m*phi)                   # invert z = L(n + m*phi)/N, so everything is rational in s
a = 1 - L; b = 1 + n*L; h = m*L*phi
al, be = s + eps, s - eps             # the antisymmetric path (7.1)
x, y = 1 - al**2, 1 - be**2           # uniquenesses keep the unit diagonal
# Exact block decomposition of F for symbolic m:
#  2n within-block contrast directions (R acts as a; Sigma as x on A-contrasts, y on B-contrasts)
#  plus the 2-dim block-mean subspace with R_mean = [[b,h],[h,b]], Sigma_mean = diag(x,y) + m*[al,be][al,be]^T
g = lambda t: t - sp.log(t) - 1
S = sp.Matrix([[x + m*al**2, m*al*be], [m*al*be, y + m*be**2]])
Rm = sp.Matrix([[b, h], [h, b]])
F = n*g(a/x) + n*g(a/y) + (S.inv()*Rm).trace() + sp.log(S.det()) - sp.log(Rm.det()) - 2
d1 = sp.simplify(sp.diff(F, eps).subs(eps, 0))
d3 = sp.simplify(sp.diff(F, eps, 3).subs(eps, 0))
E = sp.simplify(sp.diff(F, eps, 2).subs(eps, 0)/2)
rs = 1 - z; U = 1 + N*z; V = N*rs - 2*n*a
E_claim = 4*n*m*(1 + z)*((2 - L)*z - L)/((1 - z)**3*(1 + N*z))
E_three = 2*n*(a*(1 + 3*z) - rs*(1 + z))/rs**3 + N*(rs - V)/rs**2 + 4*n**2*z*V/(U*rs**2)
print("first derivative at 0 (stationarity):", d1)
print("third derivative at 0 (block-swap symmetry):", d3)
print("E - claimed closed form:", sp.simplify(sp.together(E - E_claim)))
print("E - three-term expression in (7.3):", sp.simplify(sp.together(E - E_three)))
print("checks U = b + h and V = b - h:", sp.simplify(U - (b + h)), sp.simplify(V - (b - h)))
# (7.4): (2-L)z - L = mL(2-L)/N * (phi - phi*), phi* = b/(m(2-L))
print("(7.4) identity:", sp.simplify((2 - L)*z - L - m*L*(2 - L)/N*(phi - b/(m*(2 - L)))))
# sign at phi = 0, in the original parameters (L, m)
Lp = sp.symbols('L', positive=True)
z0 = Lp*n/N
print("(2-L)z - L at phi = 0:", sp.factor(sp.simplify((2 - Lp)*z0 - Lp)))
