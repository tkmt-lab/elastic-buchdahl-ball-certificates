#!/usr/bin/env python3
"""Specification-only checker. Fixed-point integer intervals, no external libs."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json,time,resource,hashlib
BITS=192; SCALE=1<<BITS

class I:
    __slots__=('lo','hi')
    def __init__(self,a=0,b=None,raw=False):
        if raw:self.lo,self.hi=a,b
        elif isinstance(a,I):self.lo,self.hi=a.lo,a.hi
        else:
            a=Q(a);b=a if b is None else Q(b)
            self.lo=(a.numerator*SCALE)//a.denominator
            self.hi=-((-b.numerator*SCALE)//b.denominator)
        assert self.lo<=self.hi
    def __add__(self,b):
        b=I(b);return I(self.lo+b.lo,self.hi+b.hi,True)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo,True)
    def __sub__(self,b):return self+-I(b)
    def __rsub__(self,b):return I(b)+-self
    def __mul__(self,b):
        b=I(b);v=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return I(min(v)//SCALE,-((-max(v))//SCALE),True)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=I(b);assert b.lo*b.hi>0,('zero denominator',b)
        v=[Q(x*SCALE,y) for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return I(min(v).__floor__(),max(v).__ceil__(),True)
    def __rtruediv__(self,b):return I(b)/self
    def __pow__(self,n):
        if n<0:return I(1)/(self**(-n))
        if n==0:return I(1)
        if n%2==0:
            lo=0 if self.lo<=0<=self.hi else min(abs(self.lo),abs(self.hi))
            hi=max(abs(self.lo),abs(self.hi))
            den=SCALE**(n-1)
            return I(lo**n//den,-((-hi**n)//den),True)
        v=[self.lo**n,self.hi**n];den=SCALE**(n-1)
        return I(min(v)//den,-((-max(v))//den),True)
    def sqrt(self):
        assert self.lo>=0
        a=isqrt(self.lo*SCALE);b=isqrt(self.hi*SCALE)
        if b*b<self.hi*SCALE:b+=1
        return I(a,b,True)
    def absmax(self):return Q(max(abs(self.lo),abs(self.hi)),SCALE)
    def mid(self):return I((self.lo+self.hi)//2,(self.lo+self.hi)//2,True)
    def radius(self):return Q(self.hi-self.lo,2*SCALE)
    def contains(self,b):b=I(b);return self.lo<=b.lo and b.hi<=self.hi
    def row(self):return [str(self.lo),str(self.hi)]
    def approx(self):return [float(Q(self.lo,SCALE)),float(Q(self.hi,SCALE))]
    def __repr__(self):return str(self.approx())

class J:
    def __init__(self,c):self.c=list(c);self.n=len(c)-1
    def asj(self,b):return b if isinstance(b,J) else J([b]+[0]*self.n)
    def __add__(self,b):b=self.asj(b);return J([a+c for a,c in zip(self.c,b.c)])
    __radd__=__add__
    def __neg__(self):return J([-a for a in self.c])
    def __sub__(self,b):return self+-self.asj(b)
    def __rsub__(self,b):return self.asj(b)+-self
    def __mul__(self,b):
        b=self.asj(b);return J([sum((self.c[j]*b.c[i-j] for j in range(i+1)),0) for i in range(self.n+1)])
    __rmul__=__mul__
    def inverse(self):
        out=[Q(1)/self.c[0] if isinstance(self.c[0],(int,Q)) else 1/self.c[0]]
        for i in range(1,self.n+1):out.append(-out[0]*sum((self.c[j]*out[i-j] for j in range(1,i+1)),0))
        return J(out)
    def __truediv__(self,b):return self*self.asj(b).inverse()
    def __rtruediv__(self,b):return self.asj(b)*self.inverse()
    def __pow__(self,n):
        if n<0:return self.inverse()**-n
        out=self.asj(1);base=self
        while n:
            if n&1:out=out*base
            n//=2
            if n:base=base*base
        return out
    def sqrt(self):
        a=self.c[0]
        if isinstance(a,I):a=a.sqrt()
        elif isinstance(a,(Q,int)):
            a=Q(a);r,s=isqrt(a.numerator),isqrt(a.denominator);assert r*r==a.numerator and s*s==a.denominator;a=Q(r,s)
        else:a=a**.5
        out=[a]
        for n in range(1,self.n+1):out.append((self.c[n]-sum((out[i]*out[n-i] for i in range(1,n)),0))/(2*a))
        return J(out)

RC=Q(1415,2);ZC=RC/3

def model(x,Y,physical=False):
    d,u,z=8*Y[0],2*Y[1],ZC*Y[2];s=x/RC
    a=-Q(9,2)-65/u**4+70/u**2;b=Q(191,2)-130*u**2+35*u**4
    ap=260/u**5-140/u**3;bp=-260*u+140*u**3
    rho=a*d*d+b;pr=a*d*d-b;q=u*(ap*d*d+bp)/2;pt=pr+q
    pd=2*a*d;pu=ap*d*d-bp;D=1-2*s*z
    sq=D.sqrt() if isinstance(D,(I,J)) else D**.5
    G=d/(u*u*sq)-u
    f=[(q-pu*G/2-s*(rho+pr)*(z+pr)/(2*D))/pd/8,G/4,(rho-3*z)/(2*ZC)]
    if physical:
        DT=-Q(9,2)*d*d-65*u*u+35*(u**4+d*d/u**2)
        DL=a*d*d;K=-65*u*u+35*u**4;KT=-65*u*u+35*d*d/u**2
        waves=[K/DT,(d/u**3)**2*K/DL,KT/DT]
        return {'delta':d,'u':u,'D':D,'pd':pd,'DT':DT,'DL':DL,'pr':pr,'pt':pt,
          'rho-pr':rho-pr,'rho+pr':rho+pr,'rho-pt':rho-pt,'rho+pt':rho+pt,
          'compactness':2*s*z,'waves':waves}
    return f

def solve(A,b):
    n=len(b);m=[list(A[i])+[b[i]] for i in range(n)]
    for j in range(n):
        k=next(i for i in range(j,n) if m[i][j]);m[j],m[k]=m[k],m[j]
        div=m[j][j];m[j]=[q/div for q in m[j]]
        for i in range(n):
            if i!=j:
                div=m[i][j];m[i]=[u-div*v for u,v in zip(m[i],m[j])]
    return [r[-1] for r in m]

def center_coefficients(n):
    zero=J([Q(0),Q(0)]);M=[]
    for j in range(3):
        Y=[J([Q(1),Q(int(i==j))]) for i in range(3)]
        M.append([v.c[1] for v in model(zero,Y)])
    M=[[M[j][i] for j in range(3)] for i in range(3)]
    c=[[Q(1)] for _ in range(3)]
    for j in range(1,n+1):
        xx=J([Q(0),Q(1)]+[Q(0)]*(j-1))
        yy=[J(r+[Q(0)]) for r in c]
        b=[v.c[j] for v in model(xx,yy)]
        cn=solve([[j*int(i==k)-M[i][k] for k in range(3)] for i in range(3)],b)
        for i in range(3):c[i].append(cn[i])
    return M,c

if __name__=='__main__':
    start=time.perf_counter();M,c=center_coefficients(4)
    print('M',[[str(q) for q in r] for r in M]);print('norm',float(max(sum(abs(q) for q in r) for r in M)))
    print('coeff',[[str(q) for q in r] for r in c]);print('seconds',time.perf_counter()-start)

from math import comb

def polynomial(c,x):
    out=x*0
    for v in reversed(c):out=out*x+v
    return out

def center_certificate():
    p=8;M,cs=center_coefficients(p-1);x0=Q(1415,2097152);xx=I(0,x0)
    px=[polynomial(c,xx) for c in cs]
    jets=[J([sum((I(c[m])*comb(m,j)*xx**(m-j) for m in range(j,len(c))),I()) for j in range(p+1)]) for c in cs]
    residual=[q.c[p] for q in model(J([xx,I(1)]+[I()]*(p-1)),jets)]
    R=max(q.absmax() for q in residual)
    def lip(box):
        cols=[]
        for j in range(3):
            vals=[J([box[i],I(int(i==j))]) for i in range(3)]
            cols.append([v.c[1] for v in model(J([xx,I()]),vals)])
        return max(sum(cols[j][i].absmax() for j in range(3)) for i in range(3))
    rough=[a+I(-Q(1,10**12),Q(1,10**12)) for a in px]
    L=lip(rough);assert L<p
    B=2*R/(p-L)
    tube=[a+I(-B*x0**p,B*x0**p) for a in px]
    L=lip(tube);assert L<p and R+B*L<B*p
    endpoint=[polynomial(c,I(x0))+I(-B*x0**p,B*x0**p) for c in cs]
    phy=model(xx,tube,True)
    for k,v in phy.items():
        if k not in ('compactness','waves'):assert v.lo>0,(k,v)
    assert all(q.lo>0 and q.hi<SCALE for q in phy['waves'])
    return x0,endpoint,{'order':p-1,'remainder_power':p,'coefficients':[[str(q) for q in c] for c in cs],
        'M':[[str(q) for q in r] for r in M], 'x0':str(x0),'residual_coefficient_bound':str(R),
        'L_inf':str(L),'B':str(B),'C1_remainder_bound_coefficient':str(L*B+R),
        'contraction_L_over_p':str(L/p),'Y_tube':[v.row() for v in tube],
        'endpoint':[v.row() for v in endpoint], 'condition_bounds':{k:([q.row() for q in v] if k=='waves' else v.row()) for k,v in phy.items()}}

def rhs(x,y):return [q/x for q in model(x,y)]

def solution_jets(x,y,order):
    cc=[[v] for v in y]
    for n in range(order):
        xj=J([x]+([I(1)]+[I()]*(n-1) if n else []))
        vals=rhs(xj,[J(c) for c in cc])
        for j in range(3):cc[j].append(vals[j].c[n]/(n+1))
    return cc

def trial_step(x,y,h,order=7):
    xx=I(x,x+h);hh=I(0,h)
    f0=rhs(I(x),y)
    tube=[yi+hh*fi for yi,fi in zip(y,f0)]
    for _ in range(12):
        tube=[q+I(-q.radius()/8-h*h/100-Q(1,10**42),q.radius()/8+h*h/100+Q(1,10**42)) for q in tube]
        new=[yi+hh*fi for yi,fi in zip(y,rhs(xx,tube))]
        if all(a.contains(b) for a,b in zip(tube,new)):break
        tube=[I(min(a.lo,b.lo),max(a.hi,b.hi),True) for a,b in zip(tube,new)]
    else:return None
    at_start=solution_jets(I(x),y,order-1)
    rem=solution_jets(xx,tube,order)
    end=[polynomial(c,I(h))+rr[order]*I(h)**order for c,rr in zip(at_start,rem)]
    assert all(t.contains(e) for t,e in zip(tube,end))
    return tube,end


def float_pilot(x0,y):
    x=float(x0);Y=[sum(q.approx())/2 for q in y];n=0
    while x<20:
        h=min(x*.03,.002)
        def fun(t,s):return [float(v) for v in rhs(t,s)]
        a=fun(x,Y);b=fun(x+h/2,[v+h*q/2 for v,q in zip(Y,a)])
        c=fun(x+h/2,[v+h*q/2 for v,q in zip(Y,b)]);d=fun(x+h,[v+h*q for v,q in zip(Y,c)])
        Y=[v+h*(i+2*j+2*k+l)/6 for v,i,j,k,l in zip(Y,a,b,c,d)];x+=h;n+=1
        p=model(x,Y,True)
        if p['pr']<0:return {'x':x,'Y':Y,'physical':p,'steps':n}
    return {'x':x,'Y':Y}

def jacobian_rhs(xx,tube):
    cols=[]
    for j in range(3):
        vals=[J([tube[i],I(int(i==j))]) for i in range(3)]
        cols.append([q.c[1] for q in rhs(J([xx,I()]),vals)])
    return [[cols[j][i] for j in range(3)] for i in range(3)]

def propagate_radius(xx,tube,h,radius):
    JJ=jacobian_rhs(xx,tube)
    C=[[Q(JJ[i][j].hi,SCALE) if i==j else JJ[i][j].absmax() for j in range(3)] for i in range(3)]
    sigma=max(Q(0),max(-C[i][i] for i in range(3)))
    B=[[I(h*(C[i][j]+(sigma if i==j else 0))) for j in range(3)] for i in range(3)]
    norm=max(sum(Q(v.hi,SCALE) for v in row) for row in B);assert norm<1
    eye=[[I(int(i==j)) for j in range(3)] for i in range(3)]
    power=eye;E=eye
    for n in range(1,17):
        power=[[sum((power[i][l]*B[l][j] for l in range(3)),I())/n for j in range(3)] for i in range(3)]
        E=[[E[i][j]+power[i][j] for j in range(3)] for i in range(3)]
    from math import factorial
    tail=norm**17/Q(factorial(17))/(1-norm/18)
    t=I(h*sigma);term=I(1);den=I(1)
    for n in range(1,17):term=term*t/n;den=den+term
    factor=1/den  # upper bound e^(-sigma*h)
    out=[factor*sum(((E[i][j]+I(0,tail))*I(radius[j]) for j in range(3)),I()) for i in range(3)]
    return [Q(q.hi,SCALE) for q in out]

def centered_step(x,y,h,order=7):
    xx=I(x,x+h);hh=I(0,h);f0=rhs(I(x),y)
    tube=[yi+hh*fi for yi,fi in zip(y,f0)]
    for _ in range(12):
        tube=[q+I(-q.radius()/8-h*h/100-Q(1,10**42),q.radius()/8+h*h/100+Q(1,10**42)) for q in tube]
        new=[yi+hh*fi for yi,fi in zip(y,rhs(xx,tube))]
        if all(a.contains(b) for a,b in zip(tube,new)):break
        tube=[I(min(a.lo,b.lo),max(a.hi,b.hi),True) for a,b in zip(tube,new)]
    else:return None
    mid=[q.mid() for q in y]
    rad=[Q(max(q.hi-m.lo,m.lo-q.lo),SCALE) for q,m in zip(y,mid)]
    newrad=propagate_radius(xx,tube,h,rad)
    at_mid=solution_jets(I(x),mid,order-1);rem=solution_jets(xx,tube,order)
    end=[polynomial(c,I(h))+r[order]*I(h)**order+I(-e,e) for c,r,e in zip(at_mid,rem,newrad)]
    assert all(t.contains(e) for t,e in zip(tube,end))
    return tube,end

def dyadic_floor(q):
    e=q.numerator.bit_length()-q.denominator.bit_length()
    d=Q(2)**e
    if d>q:d/=2
    return d


def pressure_derivative(xx,tube):
    d,u=8*tube[0],2*tube[1]
    a=-Q(9,2)-65/u**4+70/u**2
    ap=260/u**5-140/u**3;bp=-260*u+140*u**3
    ff=rhs(xx,tube)
    return (2*a*d)*(8*ff[0])+(ap*d*d-bp)*(2*ff[1])
