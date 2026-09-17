#!/usr/bin/env python3
"""Exact-rational certificate that (A,p,C,d,v) matching forces Y and J matching.

At fixed v=J-a0 Y the conserved constraint is a scalar quadratic in Y.
This script proves its Y derivative is strictly negative on a deliberately
large box containing both the validated core and every certified tail state
at T=40.
"""
from fractions import Fraction as Q
from decimal import Decimal, getcontext
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tail"))
import tail_shape_data as dat
getcontext().prec=60
q=lambda s:Q(str(s))
class I:
    def __init__(self,lo,hi=None):
        self.lo=q(lo) if not isinstance(lo,Q) else lo
        self.hi=self.lo if hi is None else (q(hi) if not isinstance(hi,Q) else hi)
        assert self.lo<=self.hi
    def __add__(self,o):
        o=o if isinstance(o,I) else I(o);return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+(-(o if isinstance(o,I) else I(o)))
    def __rsub__(self,o):return I(o)-self
    def __mul__(self,o):
        o=o if isinstance(o,I) else I(o)
        z=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return I(min(z),max(z))
    __rmul__=__mul__
    def __pow__(self,n):
        assert n>=0
        r=I(1);b=self
        while n:
            if n&1:r=r*b
            n//=2
            if n:b=b*b
        return r

def dec(x):return f"{Decimal(x.numerator)/Decimal(x.denominator):.30E}"
T=I(40);kappa=I(Q(1,2));a0=I(q(dat.a0))
# Much wider than the actual core and tail enclosures.
A=I('1.0324','1.0326')
C=I('-0.5691','-0.5688')
d=I('-1e-14','1e-14')
Y=I('-1e-14','1e-14')
A4=A**4
cJ=I(2)*kappa*T*A4*(T*(d-I(2))-I(3)*C)
# d/dY N(A,p,C,d,J=v+a0Y,Y)
dNY=-I(2)*kappa*T*T*A4*Y+a0*cJ+I(2)*kappa*T*A4*d+I(8)*kappa*T*A4+I(6)*kappa*A4*C
print('Constraint uniqueness in fixed-v coordinates')
print('all arithmetic is exact rational interval arithmetic')
print('coefficient cJ interval [',dec(cJ.lo),',',dec(cJ.hi),']')
print('partial_Y N|v interval [',dec(dNY.lo),',',dec(dNY.hi),']')
assert dNY.hi<0
print('PASS: the constraint is strictly decreasing in Y on the common core/tail box')
print('PASS: matching (A,p,C,d,v) therefore forces Y, and then J=v+a0Y, to match')
