import sys; sys.path.insert(0, sys.argv[1])
from load import *
cbar=0.679183673469
w,r,l=WT(),RANS(),LL()
def series(x,dl,key,wt=False):
    m=x['Delta']==dl; a=x['AoA'][m]; v=x[key][m]
    d={}
    for aa,vv in zip(a,v): d.setdefault(aa,[]).append(vv)
    d={k:np.nanmean(v) if not np.all(np.isnan(v)) else np.nan for k,v in d.items()}
    if wt:  # WT stored mirrored: plot convention alpha->-alpha, C->-C (CL, CMy); CD unchanged
        s=-1 if key in('CL','CMy','CMy_AC_New','CFz','sigmaCL') else 1
        if key.startswith('sigma'): s=1
        d={-k:s*v for k,v in d.items()}
    return d
WTD=[0,10,20,30,34]; NUMD=[0,10,20,30,35]
def C(name,dl,key):
    if name=='WT': return series(w,dl,{'CM':'CMy_AC_New','sCM':'sigmaCMy_AC_New','sCL':'sigmaCL'}.get(key,key),True)
    x=r if name=='RANS' else l
    return series(x,dl,{'CM':'CMy_AC'}.get(key,key))
print('== WT alpha grid (plot convention):',sorted(C('WT',0,'CL').keys()))
print('== Table 1 (delta=0)')
for n in ['LL','RANS','WT']:
    d=C(n,0,'CL'); a=np.array(sorted(k for k in d if k>=0)); c=np.array([d[k] for k in a])
    i=np.nanargmax(c); print(f'  {n}: CL(12)={d[12]:.4f}  CLmax={c[i]:.4f} at {a[i]:g}')
sw=C('WT',0,'sCL'); print('  WT sigma CL(12)=%.4f'%sw[12])
print('== CLmax per delta (positive alpha branch)')
for k,(dw,dn) in enumerate(zip(WTD,NUMD)):
    out=[]
    for n,dl in [('WT',dw),('RANS',dn),('LL',dn)]:
        d=C(n,dl,'CL'); a=np.array(sorted(k for k in d if k>=0)); c=np.array([d[k] for k in a]); i=np.nanargmax(c)
        out.append(f'{n} d={dl}: {c[i]:.3f}@{a[i]:g}')
    print('  ',' | '.join(out))
# AC positions
print('== AC (fit |a|<=10, delta=0)')
for n,x,orig,cm in [('WT',w,0.8008,'CMy'),('RANS',r,0.228571428571,'CMy'),('LL',l,0.0,'CMy')]:
    m=(x['Delta']==0)&(np.abs(x['AoA'])<=10)
    k=np.polyfit(x['AoA'][m],x[cm][m],1)[0]/np.polyfit(x['AoA'][m],x['CFz'][m],1)[0]
    xa=orig-k*cbar if n=='WT' else orig-k*cbar
    print(f'  {n}: k={k:.4f}  x_AC={xa:.4f} m  (origin {orig})')
