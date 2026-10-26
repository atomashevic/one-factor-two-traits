from pathlib import Path
import csv,gzip,math,collections,json,hashlib,subprocess,shutil,zipfile
pack=Path(__file__).resolve().parents[1]
def read(path):
 with path.open() as f:
  reader=csv.DictReader(f)
  assert len(reader.fieldnames)==len(set(reader.fieldnames)),path
  return list(reader)
f3=read(pack/'figures/fig3.csv');f4=read(pack/'figures/fig4.csv');f5l=read(pack/'figures/fig5_left.csv');f5r=read(pack/'figures/fig5_right.csv');tab=read(pack/'figures/table4.csv')
f3keys=collections.defaultdict(list)
for r in f3:f3keys[(int(r['cell_id']),r['est'])].append(r)
f4keys={(int(r['cell_id']),r['est']):r for r in f4}
tablecells={int(r[key]) for r in tab for key in ['cell_id_500','cell_id_2000']}
f3counts=collections.Counter();f3den=collections.Counter();f4vals=collections.defaultdict(list);f4opposite=collections.Counter()
f5counts=collections.defaultdict(lambda:[0,0]);tabledefaults={};seen=set();negative=0;count=0;cellids=set()
with gzip.open(pack/'evidence/review_fits1f.csv.gz','rt',newline='') as fh:
 reader=csv.DictReader(fh);fitcolumns=reader.fieldnames
 for r in reader:
  count+=1;cell=int(r['cell_id']);cellids.add(cell);est=r['est'];key=(cell,est,int(r['rep']))
  assert key not in seen;seen.add(key)
  assert r['model']=='1F';negative+=float(r['min_psi'])<=0
  if (cell,est) in f3keys:
   for metric in f3keys[(cell,est)]:
    v=float(r[metric['index']]);assert math.isfinite(v)
    cutoff=float(metric['cutoff']);idx=(cell,est,metric['index']);f3den[idx]+=1
    f3counts[idx]+=v>=cutoff if metric['metric']=='cfi' else v<=cutoff
  if (cell,est) in f4keys:
   a=abs(float(r['lA']));b=abs(float(r['lB']));ra=float(r['lamA']);rb=float(r['lamB'])
   f4vals[(cell,est)].append(abs((a/ra-b/rb)/(a/ra+b/rb)))
   f4opposite[(cell,est)]+=float(r['lA'])*float(r['lB'])<0
  if r['design']=='A' and r['N']=='2000' and ((est in ['ML','ULS'] and r['format']=='cont') or(est in ['MLs','WLSMV'] and r['format']=='ord_sym')):
   M=abs(float(r['lA']))-abs(float(r['lB']))
   f5counts[est][0]+=M>.1;f5counts[est][1]+=abs(M)>.1
  if cell in tablecells and est=='ML':tabledefaults[(cell,int(r['rep']))]={k:r[k] for k in ['fx','chisq','df','cfi','mA','mB']}
assert count==217000 and len(cellids)==163 and negative==7
for r in f3:
 idx=(int(r['cell_id']),r['est'],r['index']);assert f3den[idx]==int(r['n']) and f3counts[idx]==int(r['n_pass'])
for r in f4:
 key=(int(r['cell_id']),r['est']);vals=f4vals[key]
 assert len(vals)==int(r['n_C_magnitude_valid'])
 assert abs(math.fsum(vals)/len(vals)-float(r['mean_abs_C_magnitude']))<2e-12
 assert f4opposite[key]==int(r['n_opposite_sign'])
for r in f5r:assert f5counts[r['est']]==[int(r['n_first']),int(r['n_selected'])]
left=collections.defaultdict(lambda:[0,0]);msrows=read(pack/'evidence/figure5_left_replications.csv')
assert len(msrows)==24000
assert len({(r['cell_id'],r['est'],r['rep']) for r in msrows})==24000
for r in msrows:
 flag=r['def_on_weak'] if r['weak']=='A' else r['def_rev_on_weak']
 assert flag in ['TRUE','FALSE'];left[int(r['cell_id'])][0]+=flag=='TRUE';left[int(r['cell_id'])][1]+=1
for r in f5l:assert left[int(r['cell_id'])]==[int(r['n_higher_objective_weak']),int(r['n_flag'])]
pairs=read(pack/'evidence/table4_replications.csv');assert len(pairs)==18000
assert len({(r['cell_id'],r['rep']) for r in pairs})==18000
stats=collections.defaultdict(lambda:[[],[],[]]);min_model=math.inf;min_baseline=math.inf;min_margin=math.inf
for r in pairs:
 cell=int(r['cell_id']);rep=int(r['rep']);N=int(r['N']);d=tabledefaults[(cell,rep)]
 for field in ['chisq','df','cfi']:assert abs(float(d[field])-float(r[field]))<1e-12
 assert abs(float(d['fx'])-float(r['default_fx']))<1e-12
 T=float(r['chisq']);df=float(r['df']);cfi=float(r['cfi']);bestfx=float(r['best_fx'])
 assert 0<cfi<1;model=T-df;baseline=model/(1-cfi)
 assert model>0 and baseline>model
 assert abs(T-2*N*float(r['default_fx']))<1e-9
 bestnc=2*N*bestfx-df;bestcfi=1-max(bestnc,0)/max(bestnc,baseline,0)
 stats[(int(r['table_row']),N)][0].append(cfi);stats[(int(r['table_row']),N)][1].append(bestcfi);stats[(int(r['table_row']),N)][2].append(r['def_on_weak']=='TRUE')
 min_model=min(min_model,model);min_baseline=min(min_baseline,baseline);min_margin=min(min_margin,baseline-model)
for r in tab:
 for N in [500,2000]:
  vals=stats[(int(r['row']),N)];assert len(vals[0])==1000
  for array,key in zip(vals,[f'default_CFI_{N}',f'best_found_CFI_{N}',f'p_higher_objective_weak_{N}']):assert abs(math.fsum(array)/len(array)-float(r[key]))<2e-12

print('PASS: portable keys, dimensions, Figures3–5 values and all18000 Table4 reconstructions.')
