#!/usr/bin/env python3
"""Exact-rational Poincare--Miranda certificate for the 5D core/tail match.

The script parses the outward core state, first-variation errors, and uniform
second-variation remainders emitted by verified_core_taylor_model.cpp.  The
tail correction bounds come from centered_tail_uniform_certificate.py.  Every
operation after parsing decimal bounds is exact rational arithmetic. Every
serialized core quantity is widened by one unit in its last printed decimal
place, independently of the producer's rounding mode.
"""
from fractions import Fraction as Q
from decimal import Decimal, getcontext
from pathlib import Path
import re,sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tail"))
import tail_shape_data as dat
getcontext().prec=80
q=lambda s:Q(str(s))
def dec(x,d=26):return f"{Decimal(x.numerator)/Decimal(x.denominator):.{d}E}"

def decimal_unit(token):
    """One unit of the last displayed decimal digit, as an exact rational.

    The C++ producers use decimal scientific or general formatting. Nearest,
    upward, and downward formatting errors are all bounded by this unit.
    Omitted trailing zeros only make this allowance more conservative.
    """
    exponent=Decimal(str(token)).as_tuple().exponent
    return Q(10)**exponent

def parsed_upper(token):
    return q(token)+decimal_unit(token)

def parsed_lower(token):
    return q(token)-decimal_unit(token)

class IV:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        self.lo=q(lo) if not isinstance(lo,Q) else lo
        self.hi=self.lo if hi is None else (q(hi) if not isinstance(hi,Q) else hi)
        assert self.lo<=self.hi
    def __add__(self,o):
        o=o if isinstance(o,IV) else IV(o);return IV(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return IV(-self.hi,-self.lo)
    def __sub__(self,o):return self+(-(o if isinstance(o,IV) else IV(o)))
    def __rsub__(self,o):return IV(o)-self
    def __mul__(self,o):
        o=o if isinstance(o,IV) else IV(o)
        z=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return IV(min(z),max(z))
    __rmul__=__mul__
    def mid(self):return (self.lo+self.hi)/2
    def rad(self):return (self.hi-self.lo)/2
    def mag(self):return max(abs(self.lo),abs(self.hi))

def deriv_shape(P,a,nu):
    out=[Q(0)]*(len(P)+1)
    for j,x in enumerate(P):out[j]-=a*x;out[j+1]+=(nu-j)*x
    while len(out)>1 and out[-1]==0:out.pop()
    return out
def peval(P,r):
    s=Q(0)
    for x in reversed(P):s=s*r+x
    return s

def parse_core(path):
    text=Path(path).read_text()
    m=re.search(r'bootstrap_ratios\s+(\S+)\s+(\S+)',text);assert m
    bootstrap=[parsed_upper(m.group(1)),parsed_upper(m.group(2))];assert max(bootstrap)<1
    m=re.search(r'core_positive_minima A D U\s+(\S+)\s+(\S+)\s+(\S+)',text);assert m
    positive=[parsed_lower(m.group(i)) for i in range(1,4)];assert min(positive)>0
    names=['A','p','C','d','Y','J'];state=[];cb=[];eb=[];cr=[];er=[];rem=[]
    for name in names:
        pat=(rf'^{name} center \[([^,]+),([^\]]+)\] width \S+\n'
             rf' col_bscale center (\S+) error (\S+)\n'
             rf' col_Rscale center (\S+) error (\S+)\n'
             rf' core_quadratic_remainder (\S+)$')
        mm=re.search(pat,text,re.M);assert mm,f'missing core block {name}'
        state.append(IV(parsed_lower(mm.group(1)),parsed_upper(mm.group(2))))
        cb.append(mm.group(3));eb.append(mm.group(4))
        cr.append(mm.group(5));er.append(mm.group(6));rem.append(mm.group(7))
    return state,cb,eb,cr,er,rem,bootstrap,positive

core_path = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[2] / "logs" / "core_reconstructed.txt"
state,cb,eb,cr,er,rem_printed,bootstrap,positive=parse_core(core_path)
A,p,C,d,Y,J=state

T=q(dat.T);rT=1/T
Ainf=q(dat.A0);beta=q(dat.beta0);a0=q(dat.a0);nu0=q(dat.nu0)
eta=q('0.000000000000000018427814530748829568597814537920972986462285544812232204587099609502860722408529466779090574297021034718078465')
PA=[q(x) for x in dat.PA];PC=[q(x) for x in dat.PC]
PY=deriv_shape(deriv_shape(PC,a0,nu0),a0,nu0);PJ=deriv_shape(PY,a0,nu0)
vnorm=peval(PJ,rT)-a0*peval(PY,rT)
PA=[x/vnorm for x in PA];PC=[x/vnorm for x in PC]
Pp=deriv_shape(PA,a0,nu0);Pd=deriv_shape(PC,a0,nu0);PY=deriv_shape(Pd,a0,nu0);PJ=deriv_shape(PY,a0,nu0)
L=[peval(P,rT) for P in (PA,Pp,PC,Pd,PY,PJ)]
LA,Lp,LC,Ld,LY,LJ=L
assert LJ-a0*LY==1

# Match the fifth coordinate v0=J-a0Y.  A separate exact constraint
# monotonicity certificate shows this also forces Y and J individually.
V=J-a0*Y
core=[A,p,C,d,V]
tail0=[IV(Ainf+eta*LA),IV(eta*Lp),IV(beta+eta*LC),IV(eta*Ld),IV(eta)]
phi0=[core[i]-tail0[i] for i in range(5)]

# Convert scaled first-variation columns to physical derivatives.  Printed
# binary128 error radii are widened by a decimal unit and inflated by one
# percent before descaling. A separate unit encloses each printed center.
def col_intervals(c,e,scale):
    out=[]
    for x,r in zip(c,e):
        m=q(x)/q(scale)
        rr=(parsed_upper(r)*q('1.01')+decimal_unit(x))/q(scale)
        out.append(IV(m-rr,m+rr))
    return out
Db6=col_intervals(cb,eb,'1e-41');DR6=col_intervals(cr,er,'1e-27')
Db=[Db6[0],Db6[1],Db6[2],Db6[3],Db6[5]-a0*Db6[4]]
DR=[DR6[0],DR6[1],DR6[2],DR6[3],DR6[5]-a0*DR6[4]]

Z=IV(0)
JI=[
 [Db[0],DR[0],IV(-1),Z,IV(-LA)],
 [Db[1],DR[1],Z,Z,IV(-Lp)],
 [Db[2],DR[2],Z,IV(-1),IV(-LC)],
 [Db[3],DR[3],Z,Z,IV(-Ld)],
 [Db[4],DR[4],Z,Z,IV(-1)],
]
JM=[[x.mid() for x in row] for row in JI]
JR=[[x.rad() for x in row] for row in JI]

def invmat(M):
    n=len(M);A=[list(M[i])+[Q(int(i==j)) for j in range(n)] for i in range(n)]
    for k in range(n):
        piv=next(i for i in range(k,n) if A[i][k]);A[k],A[piv]=A[piv],A[k]
        z=A[k][k];A[k]=[x/z for x in A[k]]
        for i in range(n):
            if i==k:continue
            z=A[i][k]
            if z:A[i]=[A[i][j]-z*A[k][j] for j in range(2*n)]
    return [r[n:] for r in A]
B=invmat(JM)
for i in range(5):
    for j in range(5):assert sum(B[i][k]*JM[k][j] for k in range(5))==Q(int(i==j))

radii=[q('5e-45'),q('5e-31'),q('3e-29'),q('2e-27'),q('2e-29')]
rem6=[parsed_upper(x)*q('1.01') for x in rem_printed]
core_rem=[rem6[0],rem6[1],rem6[2],rem6[3],rem6[5]+a0*rem6[4]]
# Uniform exact tail correction in (A,p,C,d,v0), from the centered-tail ball.
tail_rem=[q('8.9e-32'),q('5.34e-32'),q('1.5225e-29'),q('9.135e-30'),q('5e-30')]
phi_mid=[x.mid() for x in phi0];phi_rad=[x.rad() for x in phi0]
deriv_err=[sum(JR[i][j]*radii[j] for j in range(5)) for i in range(5)]
ferr=[phi_rad[i]+deriv_err[i]+core_rem[i]+tail_rem[i] for i in range(5)]
g0=[sum(B[i][k]*phi_mid[k] for k in range(5)) for i in range(5)]
gerr=[sum(abs(B[i][k])*ferr[k] for k in range(5)) for i in range(5)]
margin=[radii[i]-abs(g0[i])-gerr[i] for i in range(5)]

print('Einstein--Weyl 5D Poincare--Miranda certificate')
print('core source',core_path)
print('all arithmetic after parsing: exact rational arithmetic')
print('each serialized core quantity widened by one unit in its last displayed decimal place')
print('core bootstrap ratios',*[dec(x) for x in bootstrap])
print('core positivity lower bounds A,t+C,U',*[dec(x) for x in positive])
print('stable shape at T=40')
for n,x in zip(['LA','Lp','LC','Ld','LY','LJ'],L):print(' ',n,dec(x))
print('central matching intervals')
for i,x in enumerate(phi0):print(' ',i,'[',dec(x.lo),',',dec(x.hi),']')
print('function error budget before preconditioning')
for i in range(5):
    print(' ',i,'center_rad',dec(phi_rad[i]),'Dcore',dec(deriv_err[i]),'core2',dec(core_rem[i]),'tail',dec(tail_rem[i]),'total',dec(ferr[i]))
print('preconditioned face margins')
for i in range(5):
    print(' ',i,'radius',dec(radii[i]),'|center|',dec(abs(g0[i])),'error',dec(gerr[i]),'margin',dec(margin[i]),'fraction',dec(margin[i]/radii[i]))
assert all(x>0 for x in margin)
print('PASS: every lower face has G_i<0 and every upper face has G_i>0')
print('PASS: Poincare--Miranda gives a zero in the stated 5D box')
