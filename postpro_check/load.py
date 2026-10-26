import csv, numpy as np
from pathlib import Path
ROOT = str(Path(__file__).resolve().parents[2]) + '/'
def load(path):
    with open(ROOT+path) as f:
        rd=csv.DictReader(f); rows=list(rd)
    out={}
    for k in rows[0]:
        v=[]
        for r in rows:
            s=r[k].strip()
            try: v.append(float(s))
            except: v.append(np.nan)
        out[k]=np.array(v)
    return out
WT=lambda: load('WT/WTCoeffs_DA_0_LE_1.csv')
RANS=lambda: load('RANS/Coeffs_DA_0.csv')
LL=lambda: load('DUST/LIFTINGLINE/Coeffs_DA_0.csv')
