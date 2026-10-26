"""Measure upper fixed-region chordwise reversal in the twelve displayed RANS states.

Uses the same projected tangent and exact linear-triangle clipping as the companion
pre-stall paper's revision/phase1b/scripts/{core,analyze}.py. This is a geometric
region ahead of the hinge, not an assertion of native SU2 marker identity.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
sys.dont_write_bytecode = True
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

OUT=Path(__file__).resolve().parents[1]
ROOT=OUT.parent
STATES={0:[12,14,18],10:[8,10,18],20:[6,10,14],35:[8,10,26]}


def negative_linear_fraction(values, epsilon=0):
    """Exact fraction below -epsilon for a scalar linear on a planar triangle."""
    f=np.asarray(values)+epsilon; neg=f<0; count=neg.sum(1)
    out=np.zeros(len(f));out[count==3]=1
    for k in [1,2]:
        ids=np.flatnonzero(count==k)
        if not len(ids):continue
        isolated=np.argmax(neg[ids] if k==1 else ~neg[ids],axis=1)
        j=(isolated+1)%3;l=(isolated+2)%3
        f0=f[ids,isolated];f1=f[ids,j];f2=f[ids,l]
        frac=(-f0/(f1-f0))*(-f0/(f2-f0))
        out[ids]=frac if k==1 else 1-frac
    return out


def local_coordinates(points, delta):
    # Undo the same rigid rotation about the swept 70%-chord hinge as the companion.
    c0=0.9142857142857144; dc=c0*(.4-1)/1.6
    axis=np.array([np.tan(np.deg2rad(30))+.45*dc,1,0]);axis/=np.linalg.norm(axis)
    x,y,z=axis; K=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    a=np.deg2rad(delta); R=np.eye(3)*np.cos(a)+(1-np.cos(a))*np.outer(axis,axis)+np.sin(a)*K
    p=points.copy();c=c0+dc*p[:,1];le=np.tan(np.deg2rad(30))*p[:,1]+.25*(c0-c)
    aft=(p[:,0]-le)/c>.7
    hinge=np.array([.7*c0,0,0]);p[aft]=(p[aft]-hinge)@R+hinge
    c=c0+dc*p[:,1];le=np.tan(np.deg2rad(30))*p[:,1]+.25*(c0-c)
    return (p[:,0]-le)/c,p[:,1]/1.6


def measure(path, delta):
    reader=vtk.vtkXMLUnstructuredGridReader();reader.SetFileName(str(path));reader.Update()
    original=reader.GetOutput()
    assert original.GetNumberOfCells()>0,path
    filt=vtk.vtkDataSetTriangleFilter();filt.SetInputData(original);filt.Update();mesh=filt.GetOutput()
    assert np.all(np.diff(vtk_to_numpy(mesh.GetCells().GetOffsetsArray()))==3)
    tri=vtk_to_numpy(mesh.GetCells().GetConnectivityArray()).reshape(-1,3)
    xyz=vtk_to_numpy(mesh.GetPoints().GetData()).astype(float)
    cf=vtk_to_numpy(mesh.GetPointData().GetArray('Skin_Friction_Coefficient')).astype(float)
    v=xyz[tri];av=np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])/2
    if np.einsum('ij,ij->',v.mean(1),av)<0:av*=-1
    area=np.linalg.norm(av,axis=1);assert np.all(area>0)
    normal=av/area[:,None]
    xi,eta=local_coordinates(v.mean(1),delta)
    tip=(abs(normal[:,1])>.8)&(eta>.98)
    mask=(normal[:,2]>0)&(xi<=.7)&~tip
    tangent=np.array([1.,0,0])-normal[:,0,None]*normal
    norm=np.linalg.norm(tangent,axis=1);assert np.all(norm[mask]>1e-8)
    tangent[mask]/=norm[mask,None]
    projected=np.einsum('ijk,ik->ij',cf[tri[mask]],tangent[mask])
    assert np.isfinite(projected).all()
    assert abs((normal[mask]*tangent[mask]).sum(1)).max()<1e-12
    total=area[mask].sum();assert .6<total<.9,total
    results=[]
    for epsilon in [0,1e-5,1e-4]:
        f=negative_linear_fraction(projected,epsilon)
        assert np.all((f>=0)&(f<=1))
        rev=np.dot(area[mask],f)
        results.append(dict(epsilon=epsilon,area_m2=total,reversed_area_m2=rev,
                            percent=100*rev/total,triangles=int(mask.sum())))
    assert results[0]['percent']>=results[1]['percent']>=results[2]['percent']
    return results


if __name__=='__main__':
    fixture=np.array([[-1,-1,-1],[1,1,1],[-1,1,1],[-1,-1,1],[0,0,0]])
    assert np.allclose(negative_linear_fraction(fixture),[1,0,.25,.75,0])
    assert np.isclose(negative_linear_fraction(np.array([[-2.,2,2]]),1)[0],.0625)
    records=[];sources={}
    for delta,alphas in STATES.items():
        for alpha in alphas:
            rel=f'RANS/DA_0_SA_30_AR_5_TR_0.4_AoA_{alpha}_Delta_{delta}/surface_flow.vtu'
            path=ROOT/rel
            sources[rel]=hashlib.sha256(path.read_bytes()).hexdigest()
            for result in measure(path,delta):
                records.append(dict(delta=delta,alpha=alpha,**result))
            print(f'delta={delta}, alpha={alpha}: '+', '.join('{:.3f}%'.format(r['percent']) for r in records[-3:]),flush=True)
    target=OUT/'postpro_check/fixed_surface_reversal.csv'
    with target.open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]),lineterminator='\n');writer.writeheader();writer.writerows(records)
    (OUT/'postpro_check/fixed_surface_sources.json').write_text(json.dumps(sources,indent=2)+'\n')
    print('PASS: clipping fixtures, tangent projection, area bounds, and threshold monotonicity.')
