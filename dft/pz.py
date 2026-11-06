import re,sys
def pz_bands(fn):
    t=open(fn).read(); out=[]
    for m in re.finditer(r'==== e\(\s*(\d+)\) =\s*(-?\d+\.\d+) eV ====\s*\n(.*?)\|psi\|\^2',t,flags=re.S):
        w=sum(float(a) for a,b in re.findall(r'(\d\.\d+)\*\[#\s*(\d+)\]',m.group(3)) if int(b) in (2,6))
        out.append((int(m.group(1)),float(m.group(2)),round(w,3)))
    return out
if __name__=='__main__':
    for fn in sys.argv[1:]:
        b=pz_bands(fn); pi=[x for x in b if x[2]>0.5]
        print(fn.split('/')[0], [ (x[0],x[1],x[2]) for x in b[2:7]], ' pi-gap at M =', round(pi[1][1]-pi[0][1],4) if len(pi)>=2 else None)
