#!/usr/bin/env python3
"""Uniaxial-stress DFT (QE) driver for graphene along direction theta (deg from zigzag).
For each applied stretch eps along n, find transverse stretch et and shear g such that
sigma_tt = sigma_nt = 0 (atoms relaxed at fixed cell each step)."""
import numpy as np, subprocess, os, sys, json, re, time
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PW=os.environ.get('PW','pw.x')
PSEUDO=os.environ.get('PSEUDO_DIR',os.path.join(ROOT,'dft','pseudo'))
Ry_bohr3_to_GPa=14710.507848260711
def cfg():
    return dict(a0=float(os.environ.get('A0','2.4669')), Lz=float(os.environ.get('LZ','12.0')),
                k=int(os.environ.get('KP','24')), ec=float(os.environ.get('EC','50')))
def write_in(fn, cell, frac, C, calc='relax'):
    od=os.path.join(os.path.dirname(os.path.abspath(fn)),'tmp')
    s=f"""&control
 calculation='{calc}', prefix='g', pseudo_dir='{PSEUDO}', outdir='{od}',
 tstress=.true., tprnfor=.true., forc_conv_thr=5e-5, etot_conv_thr=1e-7, nstep=80
/
&system
 ibrav=0, nat=2, ntyp=1, ecutwfc={C['ec']}, ecutrho={8*C['ec']}, occupations='smearing', smearing='mv', degauss=0.01
/
&electrons
 conv_thr=1e-10, mixing_beta=0.4
/
&ions
/
ATOMIC_SPECIES
C 12.011 C.UPF
CELL_PARAMETERS angstrom
"""
    for v in cell: s+=f" {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n"
    s+="ATOMIC_POSITIONS crystal\n"
    for f in frac: s+=f"C {f[0]:.10f} {f[1]:.10f} {f[2]:.10f}\n"
    s+=f"K_POINTS automatic\n{C['k']} {C['k']} 1 0 0 0\n"
    open(fn,'w').write(s)
def run(fn):
    out=fn.replace('.in','.out')
    with open(out,'w') as f: subprocess.run([PW,'-in',fn],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1'))
    txt=open(out).read()
    if 'JOB DONE' not in txt: raise RuntimeError('pw failed '+out)
    E=float(re.findall(r'!\s+total energy\s+=\s+(-?\d+\.\d+)',txt)[-1])*13.605693122994
    blk=txt.split('total   stress')[-1].split('\n')[1:4]
    S=np.array([[float(x) for x in l.split()[:3]] for l in blk])*Ry_bohr3_to_GPa  # QE sign: printed = -sigma? (check)
    fr=txt.split('ATOMIC_POSITIONS (crystal)')
    if len(fr)>1:
        lines=fr[-1].strip().split('\n')[:2]; frac=np.array([[float(x) for x in l.split()[1:4]] for l in lines])
    else: frac=None
    return E,S,frac
def cell_ref(C):
    a=C['a0']; A=np.array([[a,0,0],[a/2,a*np.sqrt(3)/2,0],[0,0,C['Lz']]])
    frac=np.array([[0,0,0.5],[1/3,1/3,0.5]])  # B at (a/2, a/(2 sqrt3)) = (1/3,1/3) frac
    return A,frac
def U_of(th,eps,et,g):
    n=np.array([np.cos(th),np.sin(th)]); t=np.array([-np.sin(th),np.cos(th)])
    U2=np.eye(2)+eps*np.outer(n,n)+et*np.outer(t,t)+g*(np.outer(n,t)+np.outer(t,n))
    U=np.eye(3); U[:2,:2]=U2; return U,n,t
def point(th_deg,eps,guess,C,wd,frac0):
    th=np.radians(th_deg); sym = abs(th_deg%30)<1e-9
    A0,_=cell_ref(C); et,g=guess
    hist=[]
    for it in range(12):
        U,n,t=U_of(th,eps,et,g); cell=(U@A0.T).T
        fn=f"{wd}/e{eps:.4f}_i{it}.in"; write_in(fn,cell,frac0,C)
        E,S,fr=run(fn); sig=-S*C['Lz']*0.1   # N/m, tensile positive
        n3=np.r_[n,0]; t3=np.r_[t,0]
        snn=n3@sig@n3; stt=t3@sig@t3; snt=n3@sig@t3
        hist.append((et,g,stt,snt))
        if fr is not None: frac0=fr
        if abs(stt)<0.02 and (sym or abs(snt)<0.02): break
        # secant/Newton
        if len(hist)>=2 and abs(hist[-1][0]-hist[-2][0])>1e-7:
            dt=(hist[-1][2]-hist[-2][2])/(hist[-1][0]-hist[-2][0]); dt=dt if dt>20 else 350
        else: dt=350.
        et=et-stt/dt
        if not sym:
            if len(hist)>=2 and abs(hist[-1][1]-hist[-2][1])>1e-7:
                dg=(hist[-1][3]-hist[-2][3])/(hist[-1][1]-hist[-2][1]); dg=dg if dg>20 else 290
            else: dg=290.
            g=g-snt/dg
    U,n,t=U_of(th,eps,et,g); cell=(U@A0.T).T
    pos=fr@cell if fr is not None else frac0@cell
    return dict(theta=th_deg,eps=eps,et=et,g=g,E=E,snn=snn,stt=stt,snt=snt,cell=cell.tolist(),frac=frac0.tolist(),pos=pos.tolist(),iters=len(hist)), frac0
if __name__=='__main__':
    th=float(sys.argv[1]); epss=[float(x) for x in sys.argv[2].split(',')]
    C=cfg(); wd=os.path.join(ROOT,"dft","run",f"th{th:g}"); os.makedirs(wd,exist_ok=True)
    res_fn=f"{wd}/results.json"; res=json.load(open(res_fn)) if os.path.exists(res_fn) else []
    done={round(r['eps'],5) for r in res}
    _,frac0=cell_ref(C); et,g=0.0,0.0; prev=None
    for eps in epss:
        if round(eps,5) in done:
            r=[r for r in res if round(r['eps'],5)==round(eps,5)][0]; et,g=r['et'],r['g']; frac0=np.array(r['frac']); prev=(eps,et,g); continue
        # extrapolate guess
        gt=et*(eps/prev[0]) if prev and prev[0]>0 else -0.17*eps
        gg=g*(eps/prev[0]) if prev and prev[0]>0 else 0.0
        t0=time.time(); r,frac0=point(th,eps,(gt,gg),C,wd,np.array(frac0)); r['time']=time.time()-t0
        et,g=r['et'],r['g']; prev=(eps,et,g); res.append(r); json.dump(res,open(res_fn,'w'),indent=1)
        print(f"th={th} eps={eps:.3f} et={et:.5f} g={g:.5f} snn={r['snn']:.3f} stt={r['stt']:.3f} snt={r['snt']:.3f} it={r['iters']} t={r['time']:.0f}s",flush=True)
