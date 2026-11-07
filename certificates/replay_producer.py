#!/usr/bin/env python3
"""Independent replay of saved next2 tubes; no producer code imports."""
from independent import *
ROOT=Path(__file__).resolve().parent

class D:
    def __init__(self,v=0,g=None):
        if isinstance(v,D):self.v,self.g=v.v,list(v.g)
        else:self.v=I(v);self.g=[I(t) for t in g] if g is not None else [I()]*3
    def __add__(self,b):
        b=D(b);return D(self.v+b.v,[a+c for a,c in zip(self.g,b.g)])
    __radd__=__add__
    def __neg__(self):return D(-self.v,[-v for v in self.g])
    def __sub__(self,b):return self+-D(b)
    def __rsub__(self,b):return D(b)+-self
    def __mul__(self,b):
        b=D(b);return D(self.v*b.v,[a*b.v+self.v*c for a,c in zip(self.g,b.g)])
    __rmul__=__mul__
    def __truediv__(self,b):
        b=D(b);value=self.v/b.v
        return D(value,[(a-value*c)/b.v for a,c in zip(self.g,b.g)])
    def __rtruediv__(self,b):return D(b)/self
    def sqrt(self):
        v=self.v.sqrt();return D(v,[a/(2*v) for a in self.g])
    def __pow__(self,n):
        if n==.5:return self.sqrt()
        assert isinstance(n,int)
        if n==0:return D(1)
        return D(self.v**n,[n*self.v**(n-1)*a for a in self.g])

def physical(v):return [1+v[0],1+v[1]+v[0]/3,1+v[2]+Q(2288,1415)*v[0]]
def local(y):
    r=[q-1 for q in y];return [r[0],r[1]-r[0]/3,r[2]-Q(2288,1415)*r[0]]
def vf(x,v):
    f=rhs(x,physical(v));return [f[0],f[1]-f[0]/3,f[2]-Q(2288,1415)*f[0]]
def vf(x,v):
    # Algebraic cancellation derived from b'(u)=-u^6*a'(u).
    d=8*(1+v[0]);w=v[1]+v[0]/3;u=2*(1+w)
    z=ZC*(1+v[2]+Q(2288,1415)*v[0]);s=x/RC
    a=-Q(9,2)-65/u**4+70/u**2
    ap=260/u**5-140/u**3;b=Q(191,2)-130*u**2+35*u**4
    pr=a*d*d-b;dd=-24*v[1]-24*w*w-8*w**3
    DD=1-2*s*z;sq=DD.sqrt() if isinstance(DD,(I,J)) else DD**.5
    ff0=-ap/(32*a*u*u)*(dd*dd/x+2*z/RC*(d*d+u**6)/(sq*(1+sq)))-d*(z+pr)/(16*RC*DD)
    ff1=dd/(4*u*u*x)+d*z/(2*u*u*RC*sq*(1+sq))
    u2diff=4*w*(2+w)
    adiff=65*u2diff*(u*u+4)/(16*u**4)-70*u2diff/(4*u*u)
    bdiff=u2diff*(-130+35*(u*u+4))
    rhodiff=a*64*v[0]*(2+v[0])+64*adiff+bdiff
    ff2=(rhodiff-RC*(v[2]+Q(2288,1415)*v[0]))/(2*ZC*x)
    return [ff0,ff1-ff0/3,ff2-Q(2288,1415)*ff0]

def parse(d):return I(Q(d['lo']),Q(d['hi']))
def encloses(d,v):return Q(d['lo'])<=Q(v.lo,SCALE) and Q(v.hi,SCALE)<=Q(d['hi'])
def intersects(d,v):return Q(d['lo'])<=Q(v.hi,SCALE) and Q(v.lo,SCALE)<=Q(d['hi'])
def coefficients(x,v,n,dual=False):
    cc=[[D(q,[int(i==j) for j in range(3)]) if dual else q] for i,q in enumerate(v)]
    for order in range(n):
        xx=J([D(x) if dual else x]+([D(1) if dual else I(1)]+[D() if dual else I()]*(order-1) if order else []))
        f=vf(xx,[J(q) for q in cc])
        for j in range(3):cc[j].append(f[j].c[order]/(order+1))
    return cc

def source_center(data):
    x0,own,cert=center_certificate()
    c=[[Q(q) for q in row] for row in data['polynomial']]
    assert c==[[Q(q) for q in row[:5]] for row in cert['coefficients']]
    assert all(encloses(d,v) for d,v in zip(data['initial_box'],own))
    xx=I(0,x0);err=xx**5;P=[polynomial(row,xx) for row in c]
    # Axis-aligned enclosure of P+x^5 L^-1 [-1,1]^3.
    box=[P[0]+err*I(-1,1),P[1]+err*I(-Q(4,3),Q(4,3)),P[2]+err*I(-Q(3703,1415),Q(3703,1415))]
    jets=[J([sum((I(row[m])*comb(m,j)*xx**(m-j) for m in range(j,len(row))),I()) for j in range(6)]) for row in c]
    res=[t.c[5] for t in model(J([xx,I(1)]+[I()]*4),jets)]
    defect=[res[0],res[1]-res[0]/3,res[2]-Q(2288,1415)*res[0]]
    M=[[Q(q) for q in row] for row in cert['M']]
    # Bound |L| |F_Y-M| |L^-1|, without relying on saved derivative numbers.
    Jcols=[]
    for j in range(3):Jcols.append([q.c[1] for q in model(J([xx,I()]),[J([box[i],I(int(i==j))]) for i in range(3)])])
    H=[[ (Jcols[j][i]-M[i][j]).absmax() for j in range(3)] for i in range(3)]
    L=[[Q(1),Q(0),Q(0)],[Q(1,3),Q(1),Q(0)],[Q(2288,1415),Q(0),Q(1)]]
    E=[[sum((L[i][a]*H[a][b]*L[b][j] for a in range(3) for b in range(3)),Q()) for j in range(3)] for i in range(3)]
    bb=[Q(5),Q(13,2),Q(13,2)]
    contractions=[sum(E[i])/bb[i] for i in range(3)]
    selfmap=[(defect[i].absmax()+sum(E[i]))/bb[i] for i in range(3)]
    assert max(contractions+selfmap)<1
    return own,{'polynomial_match':True,'own_center_in_source_initial_box':True,'contraction_rows':[str(v) for v in contractions],'selfmap_rows':[str(v) for v in selfmap]}

def replay(row,full=True):
    x=Q(row['x_exact']);h=Q(row['h_exact']);X=I(x,x+h)
    start=[parse(v) for v in row['v_start']];tube=[parse(v) for v in row['tube']];end=[parse(v) for v in row['v_end']]
    pi=[a+I(0,h)*b for a,b in zip(start,vf(X,tube))]
    assert all(encloses(a,b) for a,b in zip(row['tube'],pi)),('picard',row['step'])
    jd=vf(D(X),[D(tube[i],[int(i==j) for j in range(3)]) for i in range(3)])
    lip=max(sum(v.absmax() for v in q.g) for q in jd)*h;assert lip<1
    W=[[parse(d) for d in r] for r in row['variational_tube']]
    radius=lip/(1-lip)
    assert all(encloses(row['variational_tube'][i][j],I(int(i==j)-radius,int(i==j)+radius)) for i in range(3) for j in range(3)),('variational',row['step'])
    Y=physical(tube);phy=model(X,Y,True);pe=model(I(x+h),physical(end),True)['pr'];ps=model(I(x),physical(start),True)['pr']
    for key,val in phy.items():
        if key not in ('pr','waves','compactness'):assert val.lo>0,(row['step'],key,val)
    assert all(v.lo>0 and v.hi<SCALE for v in phy['waves'])
    dp=pressure_derivative(X,Y)
    if phy['pr'].lo<=0:assert dp.hi<0
    if row['step']<464:assert pe.lo>0
    else:assert ps.lo>0 and pe.hi<0 and dp.hi<0 and phy['compactness'].lo*9>8*SCALE
    result={'step':row['step'],'picard':True,'variational_tube':True,'physics':True,'pressure':True,'C1':None}
    if full:
        # Exact rational midpoint of the stored endpoints; original midpoint is their mean.
        mid=[I((Q(d['lo'])+Q(d['hi']))/2) for d in row['v_start']]
        degree=12
        cc=coefficients(I(x),mid,degree)
        dc=coefficients(I(x),start,degree,True)
        rr=coefficients(X,tube,degree+1,True)
        center=[polynomial(q,I(h))+r[degree+1].v*I(h)**(degree+1) for q,r in zip(cc,rr)]
        E=[[polynomial([q.g[j] for q in dc[i]],I(h))+I(h)**(degree+1)*sum((rr[i][degree+1].g[l]*W[l][j] for l in range(3)),I()) for j in range(3)] for i in range(3)]
        raw=[center[i]+sum((E[i][j]*(start[j]-mid[j]) for j in range(3)),I()) for i in range(3)]
        truncated=[I(max(q.lo,t.lo),min(q.hi,t.hi),True) for q,t in zip(raw,tube)]
        result.update({'C1':all(encloses(row['flow_derivative'][i][j],E[i][j]) for i in range(3) for j in range(3)),
            'center_enclosed':all(encloses(d,v) for d,v in zip(row['flow_center'],center)),
            'endpoint_enclosed':all(encloses(d,v) for d,v in zip(row['v_end'],truncated)),
            'C1_overlap':all(intersects(row['flow_derivative'][i][j],E[i][j]) for i in range(3) for j in range(3))})
    return result

def main(limit=None):
    st=time.perf_counter();files=ROOT/'stage2-inputs'
    rows=[json.loads(s) for s in (files/'pilot-steps.jsonl').read_text().splitlines()]
    own,center=source_center(json.loads((files/'inputs-center.json').read_text()))
    assert all(encloses(d,v) for d,v in zip(rows[0]['v_start'],local(own)))
    for a,b in zip(rows,rows[1:]):
        assert Q(a['x_exact'])+Q(a['h_exact'])==Q(b['x_exact'])
        assert all(Q(v['lo'])<=Q(w['lo'])<=Q(w['hi'])<=Q(v['hi']) for v,w in zip(b['v_start'],a['v_end']))
    ans=[]
    for row in (rows if limit is None else rows[:limit]):
        ans.append(replay(row,True))
        if len(ans)%25==0:print('replayed',len(ans),'seconds',round(time.perf_counter()-st,2),flush=True)
        if time.perf_counter()-st>3300:break
    result={'center':center,'rows':ans,'row_count':len(ans),'requested_count':len(rows),'seconds':time.perf_counter()-st,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      'all_C1':all(a['C1'] for a in ans),'all_center_enclosed':all(a['center_enclosed'] for a in ans),'all_endpoints_enclosed':all(a['endpoint_enclosed'] for a in ans)}
    filename='producer-replay-pilot.json' if limit else 'producer-replay.json'
    (ROOT/filename).write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ('center','rows')})
    if limit:print(ans)
if __name__=='__main__':main()
