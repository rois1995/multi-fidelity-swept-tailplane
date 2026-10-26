import sys; sys.path.insert(0, sys.argv[1])
exec(open(sys.argv[1]+'/audit.py').read().split("print('== WT alpha grid")[0])
xAC={'WT':0.6426,'RANS':0.6420,'LL':0.7119}
def deriv(n,d0,d1,key,dd=None):
    a0=C(n,d0,key); a1=C(n,d1,key); dd=dd if dd else d1-d0
    return {k:(a1[k]-a0[k])/np.deg2rad(dd) for k in a0 if k in a1}
def mean(dic,rng): return np.mean([dic[a] for a in rng])
rngs={(0,10):range(-10,-3,2),(10,20):range(-10,-3,2)}
print('Method interval  CLd   CMd(own AC)  CMd(WT AC)  -CMd/CLd(own)')
for (d0,d1),rg in rngs.items():
    for n in ['WT','RANS','LL']:
        L=mean(deriv(n,d0,d1,'CL'),rg); M=mean(deriv(n,d0,d1,'CM'),rg)
        # move pole: CM_new = CM_old + (x_new - x_old)/cbar * CFz  (consistent with own-AC transport M_AC=M_O - k CFz, x_AC=x_O - k cbar)
        Fz=mean(deriv(n,d0,d1,'CFz'),rg) if n!='WT' else -mean(deriv(n,d0,d1,'CFz'),rg)*-1
        Mc=M+(xAC['WT']-xAC[n])/cbar*mean(deriv(n,d0,d1,'CFz'),rg)
        print(f'{n:5s} [{d0},{d1}]  {L:.3f}  {M:.3f}   {Mc:.3f}   {-M/L:.3f}')
    # WT error of the mean, as in script vs conservative
    sL=C('WT',d0,'sCL'); sL1=C('WT',d1,'sCL'); e=[np.hypot(sL[a],sL1[a])/np.deg2rad(d1-d0) for a in rg]
    print(f'   WT err CLd: script-style sqrt(sum e^2)/N={np.sqrt(np.sum(np.square(e)))/len(e):.3f}  mean pointwise e={np.mean(e):.3f}')
    sM=C('WT',d0,'sCM'); sM1=C('WT',d1,'sCM'); e=[np.hypot(sM[a],sM1[a])/np.deg2rad(d1-d0) for a in rg]
    print(f'   WT err CMd: script-style={np.sqrt(np.sum(np.square(e)))/len(e):.3f}  mean pointwise={np.mean(e):.3f}')
# sign check of CFz transport for RANS: CM_O vs CM_AC relation
m=(r['Delta']==20)&(r['AoA']==10); print('check RANS: CMy - CMy_AC =',(r['CMy'][m]-r['CMy_AC'][m])[0],' k*CFz=',(-0.6087*r['CFz'][m])[0])
# 30-34 vs 30-35
for n in ['RANS','LL']:
    a=deriv(n,30,35,'CL',4); b=deriv(n,30,35,'CL',5)
    print(n,'[30,35] CLd at a=-10..0 using /4deg:',np.round([a[k] for k in range(-10,1,2)],2),' using /5deg:',np.round([b[k] for k in range(-10,1,2)],2))
