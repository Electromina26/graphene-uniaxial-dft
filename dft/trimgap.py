"""pi-band gap at the three M points (TRIMs) for relaxed structures; uses pz projections."""
import json, numpy as np, os, sys, re, subprocess
import os
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE); import driver as d
from pz import pz_bands
th=sys.argv[1]; eps_list=[float(x) for x in sys.argv[2].split(',')]
R=json.load(open(f'{HERE}/run/th{th}/results.json'))
C=dict(a0=2.4669,Lz=12.0,k=24,ec=50)
res={}
fn_out=f'{HERE}/gap/trim_th{th}.json'
if os.path.exists(fn_out): res=json.load(open(fn_out))
for e in eps_list:
    if f'{e:.3f}' in res: continue
    r=[x for x in R if abs(x['eps']-e)<1e-6]
    if not r: continue
    r=r[0]; wd=f'{HERE}/gap/th{th}_e{e:.3f}'; os.makedirs(wd,exist_ok=True)
    cell=np.array(r['cell']); frac=np.array(r['frac'])
    if not os.path.exists(f'{wd}/scf.out') or 'JOB DONE' not in open(f'{wd}/scf.out').read():
        d.write_in(f'{wd}/scf.in',cell,frac,C,calc='scf'); d.run(f'{wd}/scf.in')
    EF=float(re.findall(r'Fermi energy is\s+(-?\d+\.\d+)',open(f'{wd}/scf.out').read())[-1])
    gaps=[]; sig=[]
    for name,kk in [('M1',(0.5,0)),('M2',(0,0.5)),('M3',(0.5,0.5))]:
        s=open(f'{wd}/scf.in').read().replace("calculation='scf'","calculation='nscf', verbosity='high'").replace("degauss=0.01","degauss=0.01, nbnd=8")
        s=re.sub(r"K_POINTS automatic\n.*\n",lambda m:f"K_POINTS crystal\n1\n{kk[0]} {kk[1]} 0.0 1\n",s)
        open(f'{wd}/{name}.in','w').write(s)
        subprocess.run([d.PW,'-in',f'{wd}/{name}.in'],stdout=open(f'{wd}/{name}.out','w'),stderr=subprocess.STDOUT)
        open(f'{wd}/proj{name}.in','w').write(f"&projwfc\n prefix='g', outdir='{wd}/tmp', filpdos='{wd}/pd{name}'\n/\n")
        subprocess.run([os.environ.get('PROJWFC','projwfc.x'),'-in',f'{wd}/proj{name}.in'],stdout=open(f'{wd}/proj{name}.out','w'),stderr=subprocess.STDOUT)
        b=pz_bands(f'{wd}/proj{name}.out'); pi=[x for x in b if x[2]>0.5 and -15<x[1]-EF<8]
        pi=sorted(pi,key=lambda x:abs(x[1]-EF))[:2]; pi=sorted(pi,key=lambda x:x[1])
        sg=[x for x in b if x[2]<0.5 and abs(x[1]-EF)<3]
        gaps.append((name,round(pi[1][1]-pi[0][1],4))); sig.append((name,[round(x[1]-EF,3) for x in sg]))
    res[f'{e:.3f}']=dict(EF=EF,pigap=gaps,sigma_near_EF=sig)
    print(th,e,gaps,sig,flush=True); json.dump(res,open(fn_out,'w'),indent=1)
