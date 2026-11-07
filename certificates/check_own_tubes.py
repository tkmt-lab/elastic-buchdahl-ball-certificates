#!/usr/bin/env python3
from independent import *
ROOT=Path(__file__).resolve().parent

def arctan_inverse(n,N):
    s=sum((Q((-1)**j,(2*j+1)*n**(2*j+1)) for j in range(N)),Q())
    r=Q((-1)**N,(2*N+1)*n**(2*N+1))
    return I(min(s,s+r),max(s,s+r))

def main():
    start=time.perf_counter();x,y,center=center_certificate();rows=[]
    with (ROOT/'stage1-tubes.jsonl').open() as f:
        for line in f:rows.append(json.loads(line))
    mins={}
    for key,value in center['condition_bounds'].items():
        if key not in ('pr','compactness','waves'):mins[key]=int(value[0])
    for j,(lo,hi) in enumerate(center['condition_bounds']['waves']):
        mins[f'wave{j}_lower']=int(lo);mins[f'wave{j}_upper']=SCALE-int(hi)
    last_positive=x;last_positive_index=0
    for idx,r in enumerate(rows):
        rx,h=Q(r['x']),Q(r['h']);assert rx==x and h>0
        yi=[I(int(a),int(b),True) for a,b in r['Y0']]
        tube=[I(int(a),int(b),True) for a,b in r['tube']]
        end=[I(int(a),int(b),True) for a,b in r['Y1']]
        assert all(a.contains(b) for a,b in zip(yi,y))
        xx=I(x,x+h)
        enclosure=[a+I(0,h)*b for a,b in zip(yi,rhs(xx,tube))]
        assert all(a.contains(b) for a,b in zip(tube,enclosure))
        mid=[q.mid() for q in yi]
        rad=[Q(max(q.hi-m.lo,m.lo-q.lo),SCALE) for q,m in zip(yi,mid)]
        nr=propagate_radius(xx,tube,h,rad)
        jets=solution_jets(I(x),mid,6);rem=solution_jets(xx,tube,7)
        calc=[polynomial(c,I(h))+rr[7]*I(h)**7+I(-e,e) for c,rr,e in zip(jets,rem,nr)]
        assert all(a.contains(b) for a,b in zip(end,calc))
        phy=model(xx,tube,True);ep=model(I(x+h),end,True)
        for key,v in phy.items():
            if key not in ('pr','compactness','waves'):
                assert v.lo>0,(idx,key,v);mins[key]=min(mins.get(key,v.lo),v.lo)
        for j,v in enumerate(phy['waves']):
            assert v.lo>0 and v.hi<SCALE
            for s,b in [('lower',v.lo),('upper',SCALE-v.hi)]:
                key=f'wave{j}_{s}';mins[key]=min(mins.get(key,b),b)
        if phy['pr'].lo<=0:assert pressure_derivative(xx,tube).hi<0
        if idx<len(rows)-1:assert ep['pr'].lo>0
        else:assert ep['pr'].hi<0
        if ep['pr'].lo>0:last_positive=x+h;last_positive_index=idx+1
        x+=h;y=end
    rootrows=rows[last_positive_index:];comps=[];dps=[]
    for r in rootrows:
        xx=I(Q(r['x']),Q(r['x'])+Q(r['h']));tube=[I(int(a),int(b),True) for a,b in r['tube']]
        comps.append(model(xx,tube,True)['compactness']);dps.append(pressure_derivative(xx,tube))
    compact=I(min(t.lo for t in comps),max(t.hi for t in comps),True)
    assert compact.lo*9>8*SCALE
    pi=16*arctan_inverse(5,96)-4*arctan_inverse(239,32)
    rr=(I(last_positive,x)/(4*pi*RC)).sqrt()
    result={'status':'PASS','tubes':len(rows),'center_C1_verified':True,'picard_all':True,'endpoint_all':True,'continuous_physics_all':True,
       'root_x':[str(last_positive),str(x)],'R':rr.row(),'R_approx':rr.approx(),'compactness':compact.row(),'compactness_approx':compact.approx(),
       'compactness_minus_8_over_9':(compact-Q(8,9)).row(),'compactness_margin_approx':(compact-Q(8,9)).approx(),
       'pressure_derivative_dx':[I(min(t.lo for t in dps),max(t.hi for t in dps),True).row()],
       'global_minima_approx':{k:float(Q(v,SCALE)) for k,v in mins.items()},'interval_bits':BITS,
       'seconds':time.perf_counter()-start,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (ROOT/'own-tubes-checked.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
