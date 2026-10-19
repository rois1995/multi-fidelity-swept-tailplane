import sys; sys.path.insert(0, sys.argv[1])
exec(open(sys.argv[1]+'/table3.py').read().split("rngs={")[0])
targets={('WT',0,10,'CL'):2.100,('WT',10,20,'CL'):1.552,('RANS',0,10,'CL'):2.427,('RANS',10,20,'CL'):2.234,('LL',0,10,'CL'):2.493,('LL',10,20,'CL'):2.416,('WT',0,10,'CM'):-0.521,('WT',10,20,'CM'):-0.384,('LL',0,10,'CM'):-0.912,('LL',10,20,'CM'):-0.977,('RANS',10,20,'CM'):-0.447}
for (n,d0,d1,key),t in targets.items():
    D=deriv(n,d0,d1,key); best=[]
    for lo in range(-24,25,2):
        for hi in range(lo+2,25,2):
            v=np.mean([D[a] for a in range(lo,hi+1,2)]); best.append((abs(v-t),lo,hi,v))
    best.sort(); print(n,d0,d1,key,'target',t,'best:',[(b[1],b[2],round(b[3],3)) for b in best[:3]])
print('pointwise WT CLd [10,20]:',{a:round(v,2) for a,v in sorted(deriv('WT',10,20,'CL').items()) if -16<=a<=2})
print('pointwise RANS CLd [10,20]:',{a:round(v,2) for a,v in sorted(deriv('RANS',10,20,'CL').items()) if -16<=a<=2})
