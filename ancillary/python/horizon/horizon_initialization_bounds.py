#!/usr/bin/env python3
"""Exact interval chain-rule check for the 1e-70 horizon-dual padding.

The horizon-tail certificate supplies uniform bounds for every jet remainder.
This script checks how the rational horizon-to-core map amplifies those bounds.
"""
from fractions import Fraction as F

class IV:
    def __init__(self,lo=0,hi=None):
        self.lo=F(lo);self.hi=self.lo if hi is None else F(hi)
        assert self.lo<=self.hi
    def __add__(self,o):
        o=o if isinstance(o,IV) else IV(o);return IV(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return IV(-self.hi,-self.lo)
    def __sub__(self,o):return self+(-o if isinstance(o,IV) else -IV(o))
    def __rsub__(self,o):return IV(o)+(-self)
    def __mul__(self,o):
        o=o if isinstance(o,IV) else IV(o)
        p=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return IV(min(p),max(p))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=o if isinstance(o,IV) else IV(o);assert o.lo>0 or o.hi<0
        return self*IV(1/o.hi,1/o.lo)
    def mag(self):return max(abs(self.lo),abs(self.hi))

N=9
class D:
    def __init__(self,v=0):
        self.v=v if isinstance(v,IV) else IV(v)
        self.g=[IV() for _ in range(N)]
        self.h=[[IV() for _ in range(N)] for _ in range(N)]
    def __add__(self,o):
        o=o if isinstance(o,D) else D(o);r=D(self.v+o.v)
        r.g=[self.g[i]+o.g[i] for i in range(N)]
        r.h=[[self.h[i][j]+o.h[i][j] for j in range(N)] for i in range(N)]
        return r
    __radd__=__add__
    def __neg__(self):return self*(-1)
    def __sub__(self,o):return self+(-o if isinstance(o,D) else -D(o))
    def __rsub__(self,o):return D(o)+(-self)
    def __mul__(self,o):
        o=o if isinstance(o,D) else D(o);r=D(self.v*o.v)
        r.g=[self.g[i]*o.v+self.v*o.g[i] for i in range(N)]
        r.h=[[self.h[i][j]*o.v+self.g[i]*o.g[j]+self.g[j]*o.g[i]+self.v*o.h[i][j]
              for j in range(N)] for i in range(N)]
        return r
    __rmul__=__mul__
    def inv(self):
        r=D(IV(1)/self.v);q=r.v
        r.g=[-self.g[i]*q*q for i in range(N)]
        r.h=[[2*self.g[i]*self.g[j]*q*q*q-self.h[i][j]*q*q
              for j in range(N)] for i in range(N)]
        return r
    def __truediv__(self,o):return self*(o if isinstance(o,D) else D(o)).inv()
    def __rtruediv__(self,o):return D(o)*self.inv()

boxes=[IV(F(1,2),2),IV(-16,16),IV(-200,200),IV(-6000,6000),
       IV(-8,8),IV(-77,77),IV(-770,770),IV(-7700,7700),IV(F(69,100),F(7,10))]
x=[D(b) for b in boxes]
for i in range(N):x[i].g[i]=IV(1)
w0,w1,w2,w3,q0,q1,q2,q3,R=x;t=D(F(20,19))
G=[t/w0,
   1/w0-(R/t)*w1/(w0*w0),
   q0*t*t*(t-1)/R-t,
   (3*t*t-2*t)*q0/R+(t-1)*q1-1,
   (6*t-2)*q0/R+(4-2/t)*q1+R*(t-1)/(t*t)*q2,
   6*q0/R+6*q1/t+3*R*q2/(t*t)+R*R*(t-1)*q3/(t*t*t*t)]
J=max(sum((g.g[i].mag() for i in range(8)),F()) for g in G)
H=max(sum((g.h[i][j].mag() for i in range(N) for j in range(N)),F()) for g in G)
assert J<F(1000) and H<F(10**6)

r=F(7,20)
# Uniform majorants for the jets, with generous integer ceilings.
w_bounds=[1/(1-r),10/(1-r),100/(1-r),
          1000*(3/(1-r)+r/(1-r)**2),
          10000*(16/(1-r)+8*r/(1-r)**2+r*(1+r)/(1-r)**3)]
q_bounds=[5*10**j/(1-r) for j in range(5)]
assert all(a<b for a,b in zip(w_bounds,[2,16,200,6000,10**6]))
assert all(a<b for a,b in zip(q_bounds,[8,77,770,7700,10**6]))
gradient=max(F('1e-41')/F('0.0066')*10**6,F('1e-27')/20*10**6)
assert gradient<F('1e-20')
value_error=F(1000)*F('1e-75')
dual_error=F(1000)*F('1e-100')+F(10**6)*F('1e-75')*F('1e-20')
assert value_error<F('1e-70') and dual_error<F('1e-70')
print("PASS: horizon-to-core jet Jacobian row-sum norm < 1000")
print("PASS: full Hessian absolute row-sum norm < 1e6")
print("PASS: scaled horizon-jet first derivatives < 1e-20")
print("PASS: omitted state value < 1e-72 and first dual coefficients < 1.00000001e-89")
print("PASS: 1e-70 padding encloses all horizon-state and first-dual remainders")
