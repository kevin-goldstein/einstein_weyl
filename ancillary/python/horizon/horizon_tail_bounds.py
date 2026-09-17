from fractions import Fraction as Q
from decimal import Decimal,getcontext
getcontext().prec=80
N=180; z=Q(7,200); L=10

def tail(j,isq=False):
    total=Q(0); last=Q(0)
    pref=5 if isq else 1; power=4 if isq else 2
    for n in range(N+1,N+1001):
        ff=1
        for k in range(j):ff*=n-k
        term=Q(pref*ff*L**n,(n+1)**power)*z**(n-j)
        total+=term;last=term
    # ratio < 0.36 here for j<=5; append factor 2 for all remainder
    return total+2*last

def dec(x): return f'{Decimal(x.numerator)/Decimal(x.denominator):.40E}'
print('N',N,'zmax',float(z))
for kind in (False,True):
 print('Q' if kind else 'W')
 for j in range(6):print(j,dec(tail(j,kind)))
# normalized parameter derivative bounds on tail jets.
delta=Q(66,10000) # 0.0066, below the available Cauchy radius
sb=Q(1,10**41);sr=Q(1,10**27);x=Q(1,20)
for kind in (False,True):
 T=[tail(j,kind) for j in range(6)]
 print('D2 Q' if kind else 'D2 W')
 for j in range(4):
  vals=[T[j], sb*T[j]/delta, sr*x*T[j+1], 2*sb*sb*T[j]/delta**2, sb*sr*x*T[j+1]/delta, sr*sr*x*x*T[j+2]]
  print(j,*[dec(v) for v in vals])
# Simple common enclosures used by horizon_second_derivative.cpp for every
# W/Q jet of order 0,...,3.
common=[Q(1,10**75),Q(1,10**114),Q(1,10**100),Q(1,10**152),Q(1,10**138),Q(1,10**124)]
for kind in (False,True):
 T=[tail(j,kind) for j in range(6)]
 for j in range(4):
  vals=[T[j],sb*T[j]/delta,sr*x*T[j+1],2*sb*sb*T[j]/delta**2,sb*sr*x*T[j+1]/delta,sr*sr*x*x*T[j+2]]
  assert all(v<c for v,c in zip(vals,common))
print('PASS common D2 tail box:',*[dec(c) for c in common])
