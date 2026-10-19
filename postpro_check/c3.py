import sys; sys.path.insert(0, sys.argv[1])
exec(open(sys.argv[1]+'/table3.py').read().split("rngs={")[0])
for key in ['CL','CM']:
    D=deriv('WT',20,30,key); print('WT',key,'[20,30]',{a:round(D[a],2) for a in sorted(D) if -12<=a<=8})
    D=deriv('RANS',20,30,key); print('RANS',key,'[20,30]',{a:round(D[a],2) for a in sorted(D) if -12<=a<=8})
