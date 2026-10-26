import numpy as np
exec(open('scope_probe.py').read().split('print("Probe 1')[0])
print("Probe 3: mildly unequal blocks, 4 items each, loadings .68 (A) vs .72 (B); exchangeable phi* at .70 is .409")
for phi in np.round(np.arange(0,0.91,0.05),2):
    R=popg([.68]*4,[.72]*4,phi); b=full_best(R,30); l=np.abs(b.x[:8])
    print(f"  phi={phi:.2f}: A {l[:4].mean():.3f}  B {l[4:].mean():.3f}  gap {l[4:].mean()-l[:4].mean():+.3f}")
