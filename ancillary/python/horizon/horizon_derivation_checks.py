#!/usr/bin/env python3
"""Exact checks accompanying the explicit horizon-majorant derivation.

The all-index inequalities are proved in the manuscript. This script checks
their rational constants, matches the six bounds to the production certificate,
and compares the two recurrences as identities in Q[b] through order 12.
"""
from fractions import Fraction as F
from pathlib import Path
import contextlib
import importlib.util
import io

here = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("majorant", here / "horizon_majorant_certificate.py")
majorant = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(majorant)

assert sum((F(1, k*k) for k in range(1, 6)), F()) + F(1, 5) < F(5, 3)
assert F(1, 8) + F(1, 8) == F(1, 4)
assert sum((F(1, k**4) for k in range(2, 6)), F()) + F(1, 3*5**3) < F(1, 12)
H = sum((F(1, k) for k in range(2, 103)), F())
assert H / 104 < F(1, 10)
assert F(13, 12) > F(104, 103)
for j in range(6):
    # Dropping the final power factor gives an upper bound that decreases
    # with n; its first value therefore controls all n>=181.
    assert F(7,20)*F(182,182-j)<F(9,25)

def six_bounds(n):
    return (F(5,12), F(5,6*(n+4))+F(8,(n+2)**3),
            F(5,2*(n+4))+F(4,(n+2)**2), F(1,30*(n+1)**2),
            F(5,3*(n+2)**2)+F(200,9*(n+2)**3), F(100,9*(n+1)**2))

for n in range(100, 301):
    phi = F((n+2)**2, (n+1)**2)
    bounds = six_bounds(n)
    assert all(x <= phi*y for x, y in zip(majorant.w_parts(n), bounds))
    assert phi*sum(bounds, F()) == majorant.w_asymptotic_bound(n)
    assert majorant.q_ratio(n) <= majorant.q_asymptotic_bound(n)

class P:
    """Small exact polynomial in b, for the recurrence identity regression."""
    def __init__(self, x=0):
        self.c = tuple(x) if isinstance(x, (list, tuple)) else (F(x),)
        while len(self.c)>1 and self.c[-1]==0: self.c=self.c[:-1]
    def __add__(self, other):
        other=other if isinstance(other,P) else P(other)
        n=max(len(self.c),len(other.c))
        return P([(self.c[i] if i<len(self.c) else 0)+(other.c[i] if i<len(other.c) else 0) for i in range(n)])
    __radd__=__add__
    def __neg__(self): return P([-x for x in self.c])
    def __sub__(self,other): return self+(-other if isinstance(other,P) else -P(other))
    def __rsub__(self,other): return P(other)+(-self)
    def __mul__(self,other):
        other=other if isinstance(other,P) else P(other)
        out=[F()]*(len(self.c)+len(other.c)-1)
        for i,x in enumerate(self.c):
            for j,y in enumerate(other.c): out[i+j]+=x*y
        return P(out)
    __rmul__=__mul__
    def __truediv__(self,other): return P([x/F(other) for x in self.c])
    def __eq__(self,other): return self.c==(other if isinstance(other,P) else P(other)).c

N=12; b=P([F(0),F(1)])
w=[P() for _ in range(N+1)]; q=[P() for _ in range(N+2)]
w[0]=P(1);w[1]=1+b;q[0]=P(1);q[1]=-2-3*b;q[2]=1+3*b+3*b*b
for n in range(1,N):
    s=(n+1)*sum(((n+1-k)*q[k]*w[n+1-k] for k in range(1,n+1)),P())
    s+=sum((k*(k+1)*q[k]*w[n+1-k]/6 for k in range(1,n+2)),P())
    w[n+1]=(w[n]/3-s)/(n+1)**2
    j=n-1
    source=sum(((j-i+2)*(j-i+1)*w[i]*w[j-i+2]
                -2*(i+1)*(j-i+1)*w[i+1]*w[j-i+1] for i in range(j+1)),P())
    q[j+3]=-6*source/((j+1)*(j+2)*(j+3)*(j+4))
a=[P() for _ in range(N+1)];g=[P() for _ in range(N+1)]
a[1]=P(1);g[1]=P(1);g[2]=1+b
for ell in range(1,N):
    s=sum(((-1)**i*g[i]*(1+b*a[ell+1-i])
           * (F((ell+1)*(ell+1-i))+F(i*(i+1),6)) for i in range(1,ell+2)),P())
    a[ell+1]=((2*ell*ell+2*ell+1)*a[ell]-ell*ell*a[ell-1]-3*s)/(ell+1)**2
    if ell+2<=N:
        s=sum(((a[i]+a[ell+1-i]*(1+b*a[i]))*(ell+1-i)*(ell-3*i) for i in range(ell+1)),P())
        g[ell+2]=(-1)**(ell+1)*s/F((ell+3)*(ell+2)*(ell+1)*ell,2)
for n in range(1,N+1):
    assert w[n]==1+b*a[n]
    want=(-2-3*b if n==1 else 1+3*b*g[2] if n==2 else 3*b*(-1)**n*g[n])
    assert q[n]==want
# Exact horizon value of the omitted field-equation constraint.
assert -b == q[2]-q[1]*q[1]/3+F(1,3)

print("PASS: rational zeta bounds and all-index harmonic-ratio monotonicity constants")
print("PASS: differentiated geometric-tail ratios are below 9/25 for all n>=181, j<=5")
print("PASS: all six displayed high-order bounds agree with the production majorant")
print("PASS: direct and alpha/gamma recurrences agree in Q[b] through degree in z",N)
print("PASS: horizon constraint vanishes identically in b")
