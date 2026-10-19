import sys; sys.path.insert(0, sys.argv[1])
exec(open(sys.argv[1]+'/table3.py').read().split("rngs={")[0])
for d0,d1 in [(0,10),(10,20)]:
    print('WT CLd',d0,d1,{a:round(v,2) for a,v in sorted(deriv('WT',d0,d1,'CL').items()) if -16<=a<=4})
# CDmin
for n,ds in [('WT',WTD),('RANS',NUMD),('LL',NUMD)]:
    out=[]
    for dl in ds:
        d=C(n,dl,'CD'); k=min(d,key=d.get); out.append(f'd{dl}:{d[k]:.4f}@{k:g}')
    print(n,'CDmin',out)
# RANS component lift (flap = UpFlap+DownFlap)
up=load('RANS/Coeffs_DA_0_UpFlap.csv'); dn=load('RANS/Coeffs_DA_0_DownFlap.csv')
for dl in [20,35]:
    a=series(up,dl,'CL'); b=series(dn,dl,'CL')
    print('RANS flap CL d',dl,{k:round(a[k]+b[k],3) for k in sorted(a) if -16<=k<=8})
