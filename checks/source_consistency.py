"""Source links and selected displayed numbers; not a universal scientific check."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import math
import re

root = Path(__file__).resolve().parents[1]
tex = '\n'.join(p.read_text() for p in (root/'manuscript').rglob('*.tex'))
labels = re.findall(r'\\label\{([^}]+)\}', tex)
refs = re.findall(r'\\(?:eqref|ref|autoref)\{([^}]+)\}', tex)
assert len(labels) == len(set(labels)), 'Duplicate labels'
assert not set(refs)-set(labels), 'Unresolved source references'
print(f'PASS: {len(labels)} unique labels and {len(refs)} internal references')
for name in re.findall(r'\\input\{([^}]+)\}', tex):
    assert (root/'manuscript'/f'{name}.tex').is_file(), name
figures = re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', tex)
assert len(figures) == 5
for name in figures:
    assert (root/'manuscript'/name).is_file() or (root/'manuscript/figures'/name).is_file(), name
manifest = json.loads((root/'provenance/preserved-files.json').read_text())
for name, digest in manifest.items():
    assert hashlib.sha256((root/name).read_bytes()).hexdigest() == digest, name
print(f'PASS: {len(manifest)} preserved figure/proof/result files match reconciled inputs')
aux = (root/'build/manuscript/main.aux').read_text()
numbering = {'thm:ml':1,'thm:uls':2,'thm:cusp':3,'lem:uniformity':1,'cor:ml-transition':1,
             'cor:wls':1,'prop:unequal':2,'prop:crossloading':3,'prop:boundary':4,'prop:ceiling':5}
for label, n in numbering.items():
    assert f'\\newlabel{{{label}}}{{{{{n}}}' in aux, label
# Figures are numbered by order of appearance; fig:3 now appears in Appendix G.
for label, i in {'fig:1':1,'fig:2':2,'fig:4':3,'fig:5':4,'fig:3':5}.items():
    assert f'\\newlabel{{{label}}}{{{{{i}}}' in aux, label
print('PASS: theorem/proposition/corollary/lemma and figure numbering')
log = (root/'build/manuscript/main.log').read_text()
assert not re.search(r'LaTeX Warning|Overfull|Underfull|Missing character|undefined', log)
print('PASS: manuscript log has no LaTeX warnings, box warnings, or missing characters')

num = (root/'manuscript/sections/06_numerical.tex').read_text()
table = num.split(r'\label{tab:4}',1)[1].split(r'\bottomrule',1)[0].split(r'\midrule',1)[1]
rows = [line.split('&') for line in table.splitlines() if '&' in line]
with (root/'figures/table4.csv').open() as f:
    data = list(csv.DictReader(f))
assert len(rows) == len(data) == 9
def shown(p):
    # Percentages are rounded to integers, except that values which would round to 0 or 100 keep one decimal.
    v = 100*p; r = round(v)
    return f'{v:.1f}' if (r in (0, 100) and v not in (0.0, 100.0)) else str(r)
keys = ['population_CFI','default_CFI_500','default_CFI_2000','best_found_CFI_500','best_found_CFI_2000']
for row, d in zip(rows, data):
    assert float(row[1]) == float(d['phi'])
    for cell, key in zip(row[2:7], keys):
        assert float(cell) == round(float(d[key]),3), (cell, key)
    percentages = re.findall(r'([\d.]+)\\%',row[7])
    assert percentages == [shown(float(d[f'p_higher_objective_weak_{n}'])) for n in (500,2000)], (row[7], d)
print('PASS: all 45 Table 4 CFI cells and 18 percentages match archived exports')
with (root/'figures/fig3.csv').open() as f:
    f3 = list(csv.DictReader(f))
for est,index,expected in [('WLSMV','cfi_rob',83),('ML','cfi',99)]:
    d = [r for r in f3 if r['est']==est and r['index']==index and float(r['phi'])==.912]
    assert len(d)==1 and round(100*float(d[0]['p_accept']))==expected
assert '83\\% against 99\\%' in num and 'plus .01, .24 and .409' in num
print('PASS: CFI-specific 83%/99% comparison and unequal-loading grid preserved')

def fit(m, L, phi):
    a=1-L; n=m-1; b=1+n*L; q=m*L*phi/b; threshold=b/(m*(2-L)); p=2*m
    FB=-2*n*math.log(a)-math.log(b*b-(m*L*phi)**2)
    if phi < threshold:
        F=-n*math.log(a)-math.log(b)+n*math.log1p(-q*q)
        S=L+(2-L)*q*q
        alpha=(math.sqrt(S+2*q)+math.sqrt(max(0,S-2*q)))/2
        beta=q/alpha
    else:
        k=m*L*(1-phi)/a
        F=(p-1)*math.log1p(k/(p-1))-math.log1p(k)
        alpha=beta=math.sqrt(L*(n+m*phi)/(p-1))
    df=p*(p-3)/2; dfB=p*(p-1)/2
    SRMR=math.sqrt((m*n/2*((L-alpha*alpha)**2+(L-beta*beta)**2)+m*m*(L*phi-alpha*beta)**2)/(m*(2*m+1)))
    return [1-F/FB,1-(F/df)/(FB/dfB),math.sqrt(max(0,F)/df),SRMR,threshold]

def boundary(m,L,i):
    cutoff=[.95,.95,.06,.08][i]
    passes=lambda phi: fit(m,L,phi)[i]>=cutoff if i<2 else fit(m,L,phi)[i]<=cutoff
    if passes(0): return 0
    lo,hi=0.,1-1e-10
    for _ in range(60):
        mid=(lo+hi)/2
        if passes(mid): hi=mid
        else: lo=mid
    return hi

section=(root/'manuscript/sections/04_population_fit.tex').read_text()
table=section.split(r'\midrule',1)[1].split(r'\bottomrule',1)[0]
count=0
for line in table.splitlines():
    if '&' not in line: continue
    cells=[c.strip().replace('\\\\','').strip() for c in line.split('&')]; L=float(cells[0])**2
    if cells[1]=='many':
        cfi=1-2*(1-L)/L*((1-L)**(-.05)-1)
        values=[cfi,cfi,0,max(0,1-.16/L),L/(2-L)]
    else:
        m=int(cells[1])//2
        values=[boundary(m,L,i) for i in range(4)]+[fit(m,L,0)[4]]
    for cell,value in zip(cells[2:],values):
        assert (value==0 if cell.startswith('any') else float(cell)==round(value,2)),(cells,value)
        count+=1
print(f'PASS: all {count} Table 1 cutoff/threshold cells reconstructed from displayed formulas')
