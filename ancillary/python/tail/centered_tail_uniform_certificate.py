#!/usr/bin/env python3
"""Exact-rational uniform a posteriori certificate for the centered Einstein--Weyl tail.

The explicit approximate tail is X0(t)=eta*g(t)*(P_A(1/t),...,P_C'''(1/t)),
where g=exp[-a0(t-T)](t/T)^nu0 and the decimal coefficients in
`tail_shape_data.py` are treated as exact rational numbers.  The script bounds
its full nonlinear ODE residual on t>=T and proves contraction of the integral
operator in a weighted ball around X0.
"""
from fractions import Fraction as Q
from decimal import Decimal,getcontext,localcontext,ROUND_CEILING,ROUND_FLOOR
from math import factorial
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import tail_shape_data as dat
getcontext().prec=80
q=lambda s:Q(str(s))
def dec(x,d=24,lower=False):
 # Both division and final scientific formatting use the required direction.
 # Check the formatted decimal as an exact rational, independently of context.
 with localcontext() as ctx:
  ctx.prec=max(80,d+16)
  ctx.rounding=ROUND_FLOOR if lower else ROUND_CEILING
  out=f"{Decimal(x.numerator)/Decimal(x.denominator):.{d}E}"
 assert q(out)<=x if lower else q(out)>=x
 return out
def trim(P):
 while len(P)>1 and P[-1]==0:P.pop()
 return P
def add(P,Qp):
 n=max(len(P),len(Qp));return trim([(P[i] if i<len(P) else 0)+(Qp[i] if i<len(Qp) else 0) for i in range(n)])
def scale(P,c):return trim([c*x for x in P])
def shift(P,k):return [Q(0)]*k+P
def deriv_shape(P,a,nu):
 # D[g P(r)]/g = (-a+nu r)P-r^2 P_r
 out=[Q(0)]*(len(P)+1)
 for j,x in enumerate(P):
  out[j]-=a*x;out[j+1]+=(nu-j)*x
 return trim(out)
def eval_abs(P,r):
 # sum |p_j| r^j
 s=Q(0);pw=Q(1)
 for x in P:s+=abs(x)*pw;pw*=r
 return s
def peval(P,r):
 s=Q(0)
 for x in reversed(P):s=s*r+x
 return s

T=q(dat.T); rT=1/T; sigma=q('0.6')
A0=q(dat.A0); beta0=q(dat.beta0); a0=q(dat.a0); nu0=q(dat.nu0)
# These finite decimals are the exact rational approximation parameters.  The
# identities a0=1/A0 and nu0=beta0/(2*A0) are NOT assumed by this certificate.
# Their small discrepancies are included when the polynomial residual is formed.
assert abs(a0-1/A0)<q('1e-178')
assert abs(nu0-beta0/(2*A0))<q('1e-178')
assert nu0<=0 and a0>sigma>0
assert a0*T>1 and sigma*T>1
eta=q('0.000000000000000018427814530748829568597814537920972986462285544812232204587099609502860722408529466779090574297021034718078465')
# Five-dimensional Poincare box: only the three tail parameters enter here.
rA=q('3e-29'); rbeta=q('2e-27'); reta=q('2e-29')
Apar_lo=A0-rA; Apar_hi=A0+rA; beta_lo=beta0-rbeta; beta_hi=beta0+rbeta
eta_abs=abs(eta)+reta
PA=[q(x) for x in dat.PA];PC=[q(x) for x in dat.PC]
# Exact re-normalization so that (C'''-a0 C'')(T)=1.
PY=deriv_shape(deriv_shape(PC,a0,nu0),a0,nu0)
PJ=deriv_shape(PY,a0,nu0)
vnorm=peval(PJ,rT)-a0*peval(PY,rT)
PA=scale(PA,1/vnorm);PC=scale(PC,1/vnorm)
Pp=deriv_shape(PA,a0,nu0);Ppp=deriv_shape(Pp,a0,nu0)
Pd=deriv_shape(PC,a0,nu0);PY=deriv_shape(Pd,a0,nu0);PJ=deriv_shape(PY,a0,nu0);PJp=deriv_shape(PJ,a0,nu0)
assert peval(PJ,rT)-a0*peval(PY,rT)==1

onebr=[Q(1),beta0]
# N_F=(1+beta r)Ppp+r(2+beta r)Pp-A r^2 Pd/3-A r PY/6
NF=add(add(add(
    # polynomial multiplication by 1+beta r
    add(Ppp,scale(shift(Ppp,1),beta0)),
    add(scale(shift(Pp,1),2),scale(shift(Pp,2),beta0))),
    scale(shift(Pd,2),-A0/Q(3))),
    scale(shift(PY,1),-A0/Q(6)))
# N_LF=(1+beta r)P_FL
NLF=add(add(scale(shift(Pp,1),-2),scale(shift(Pp,2),-beta0)),add(scale(shift(Pd,2),A0/Q(3)),scale(shift(PY,1),A0/Q(6))))
# N_J=r(1+beta r)PJp-6/A^3 N_LF
NJ=add(add(shift(PJp,1),scale(shift(PJp,2),beta0)),scale(NLF,-Q(6)/(A0**3)))

d0=Q(1)+beta0/T
assert d0>0
K_NF=eval_abs(NF,rT)
assert NJ[0]==0  # otherwise division by r would not be uniformly bounded
K_NJ=sum(abs(c)*T**(1-j) for j,c in enumerate(NJ)) # accounts for 1/r
K_RF_lin=eta_abs*K_NF/d0
K_RJ_lin=eta_abs*K_NJ/d0

# Shape component maxima, valid after multiplication by exp[-a0(t-T)].
KA,Kp,Kc,Kd,KY,KJ=[eval_abs(P,rT) for P in [PA,Pp,PC,Pd,PY,PJ]]
alpha=eta_abs*KA;p=eta_abs*Kp;c=eta_abs*Kc;d=eta_abs*Kd;Y=eta_abs*KY;J=eta_abs*KJ
Alo=Apar_lo-alpha;Ahi=Apar_hi+alpha;D0min=T+beta_lo;Dmin=D0min-c
assert Alo>1 and Dmin>39
# Uniform parameter perturbation of the linearized F operator.  The fixed
# generalized-Yukawa shape uses (A0,beta0), while the algebraic constants range
# over the Poincare box.  The coefficient of p has beta derivative exactly
# +1/(t+beta)^2: the coefficient is -1/t-1/(t+beta).
ratio_par=rA/D0min + A0*rbeta/(D0min*(T+beta0))
K_RF_par=Y*ratio_par/Q(6)+d*ratio_par/(Q(3)*T)+p*rbeta/(D0min**2)
# Quadratic nonlinear correction F-F_linear, uniformly in the parameter box.
common=alpha/Dmin + Apar_hi*c/(Dmin*D0min)
K_F_nl=(Q(2)*p*p/Alo + p*d/Dmin + p*c/(Dmin*D0min)
        + Y*common/Q(6) + d*common/(Q(3)*T))
# Bound the actual linear F uniformly.
beta_abs=max(abs(beta_lo),abs(beta_hi))
K_FL=( (Q(2)*T+beta_abs)*p/(T*D0min)
          + Apar_hi*d/(Q(3)*T*D0min)+Apar_hi*Y/(Q(6)*D0min) )
# Difference of inverse cubes due to the Yukawa profile and to the A_infty box.
inv3dyn=alpha*(Ahi*Ahi+Ahi*Apar_hi+Apar_hi*Apar_hi)/(Alo**3*Apar_lo**3)
inv3par=rA*(Apar_hi*Apar_hi+Apar_hi*A0+A0*A0)/(Apar_lo**3*A0**3)
# The parameter residual retains the SAME exponential rate a0.  In t*RF_par
# and t*F_linear, the remaining coefficient bounds are maximized at T:
# t/(t+beta_lo), t/(t+beta_lo)^2, and
# t/[(t+beta_lo)(t+beta0)] decrease on [T,infinity).
# The last derivative has numerator beta_lo*beta0-t^2.
assert beta_lo<=beta0<=beta_hi<0 and beta_lo*beta0<T*T
# The nonlinear residual contains two profile factors; a0*T>1 gives
# t*exp[-2*a0*(t-T)] <= T*exp[-a0*(t-T)].
K_RJ_par=Q(6)*T*(K_RF_par/(Alo**3)+K_FL*inv3par)
K_RJ_nl=Q(6)*T*(K_FL*inv3dyn + K_F_nl/(Alo**3))
K_RF=K_RF_lin+K_RF_par+K_F_nl
K_RJ=K_RJ_lin+K_RJ_par+K_RJ_nl

# Lipschitz matrix on a weighted ball about X0.
# Perron-near-optimal rational weights.
weights=[q('0.01780'),q('0.01068'),q('3.045'),q('1.827'),q('0.2350'),q('1')]
qA,qp,qC,qd,qu,qv=weights
radius=q('5e-30')
# Center pointwise maxima plus ball maxima.
AloB=Apar_lo-alpha-qA*radius;AhiB=Apar_hi+alpha+qA*radius
Clo=beta_lo-c-qC*radius;Chi=beta_hi+c+qC*radius;Cabs=max(abs(Clo),abs(Chi))
DminB=T+Clo
amin=a0 # central parameter for this certificate
pmax=p+qp*radius;dmax=d+qd*radius
# u=J+aY and v=J-aY for center. Use direct component bounds conservatively.
umax=J+a0*Y+qu*radius;vmax=J+a0*Y+qv*radius
Ymax=(umax+vmax)/(Q(2)*amin)
ratio=T/DminB
FA=Q(2)*pmax*pmax/(AloB*AloB)+Ymax/(Q(6)*DminB)+dmax/(Q(3)*T*DminB)
Fp=Q(4)*pmax/AloB+dmax/DminB+Cabs/(T*DminB)+Q(2)/DminB
FC=AhiB*Ymax/(Q(6)*DminB**2)+AhiB*dmax/(Q(3)*T*DminB**2)+pmax*dmax/DminB**2+pmax/DminB**2
Fd=AhiB/(Q(3)*T*DminB)+pmax/DminB
FY=AhiB/(Q(6)*DminB)
Fabs=Q(2)*pmax*pmax/AloB+AhiB*Ymax/(Q(6)*DminB)+AhiB*dmax/(Q(3)*T*DminB)+pmax*dmax/DminB+Cabs*pmax/(T*DminB)+Q(2)*pmax/DminB
kappa=Q('0.5')
HA=Q(6)*T*pmax*pmax/(kappa*AloB**5)+ratio*Ymax/(Q(2)*kappa*AloB**3)+dmax/(kappa*AloB**3*DminB)+Q(9)*T*Fabs/(kappa*AloB**4)
Hp=Q(12)*T*pmax/(kappa*AloB**4)+Q(3)*ratio*dmax/(kappa*AloB**3)+Q(3)*Cabs/(kappa*AloB**3*DminB)+Q(6)*ratio/(kappa*AloB**3)
HC=T*Ymax/(Q(2)*kappa*AloB**2*DminB**2)+dmax/(kappa*AloB**2*DminB**2)+Q(3)*T*pmax*dmax/(kappa*AloB**3*DminB**2)+Q(3)*T*pmax/(kappa*AloB**3*DminB**2)
Hd=Q(1)/(kappa*AloB**2*DminB)+Q(3)*ratio*pmax/(kappa*AloB**3)
# H_Y=3t/(kappa A^3)*A/[6(t+C)] - a0^2 = t/(2 kappa A^2(t+C))-a0^2
# For C<0, 1 <= t/(t+C) <= T/(T+Clo).  Include the limiting
# endpoint at infinity explicitly; evaluating both endpoints at T is insufficient.
assert Chi<0 and Clo>-T
HY_lo=Q(1)/(Q(2)*kappa*AhiB**2)-a0**2
HY_hi=ratio/(Q(2)*kappa*AloB**2)-a0**2
HY=max(abs(HY_lo),abs(HY_hi))
Yuv=Q(1)/(Q(2)*amin)
JF=[FA,Fp,FC,Fd,FY*Yuv,FY*Yuv];JH=[HA,Hp,HC,Hd,HY*Yuv,HY*Yuv]
M=[[Q(0) for _ in range(6)] for _ in range(6)]
for j in range(6):M[0][j]=JF[j]/sigma**2;M[1][j]=JF[j]/sigma
M[2][4]=M[2][5]=Yuv/sigma**2;M[3][4]=M[3][5]=Yuv/sigma
for j in range(6):M[4][j]=JH[j]/(amin+sigma);M[5][j]=JH[j]/(amin-sigma)
rows=[sum(M[i][j]*weights[j] for j in range(6))/weights[i] for i in range(6)]
contract=max(rows)
# Defect of the integral operator from the differential residuals.
# alpha,p components; c,d are exact kinematic identities; u,v share R_J.
e_lower=q('2.718281828459045')
assert e_lower<sum((Q(1,factorial(k)) for k in range(19)),Q(0))
defects=[K_RF/(a0*a0),K_RF/a0,Q(0),Q(0),K_RJ/(Q(2)*a0),K_RJ/(e_lower*(a0-sigma))]
defect_norm=max(defects[i]/weights[i] for i in range(6))
selfmap=contract+defect_norm/radius

print('Centered Einstein--Weyl tail certificate (uniform Poincare tail box)')
print('all inequalities below use exact rational arithmetic')
print('approximation exponents are exact rational decimal data')
print('residual convention: |R_F(t)| <= K_RF exp[-a0(t-T)], likewise R_J')
print('all printed upper bounds are rounded upward; lower endpoints downward')
print('tail parameter radii',dec(rA),dec(rbeta),dec(reta))
print('exact v(T) normalization error: 0')
print('linear RF bound',dec(K_RF_lin));print('parameter RF bound',dec(K_RF_par));print('nonlinear RF bound',dec(K_F_nl));print('total RF bound',dec(K_RF))
print('linear RJ bound',dec(K_RJ_lin));print('parameter RJ bound',dec(K_RJ_par));print('nonlinear RJ bound',dec(K_RJ_nl));print('total RJ bound',dec(K_RJ))
print('shape maxima at T')
for n,x in zip(['alpha','p','c','d','Y','J'],[alpha,p,c,d,Y,J]):print(' ',n,dec(x))
print('weighted contraction rows')
for i,x in enumerate(rows):print(' ',i,dec(x))
print('contraction',dec(contract))
print('defect norm',dec(defect_norm))
print('radius',dec(radius))
print('defect/radius',dec(defect_norm/radius))
print('self-map',dec(selfmap))
print('H_Y enclosure',dec(HY_lo,lower=True),dec(HY_hi))
assert K_RF<q('4e-33')
assert K_RJ<q('8e-31')
assert contract<q('0.59')
assert defect_norm<q('2e-30')
assert selfmap<1
# Positivity and monotonicity on the complete tail.  The tail norm gives
# |delta p(t)| <= qp*radius*exp[-sigma(t-T)].  Since a0*T>1 and sigma*T>1,
# both t times the approximate p profile and t times its correction are
# maximized at t=T.  Thus U=1-t A'/A is uniformly positive.
assert a0*T>1 and sigma*T>1
Ulo=Q(1)-T*(p+qp*radius)/AloB
assert Clo>-1
Blo=Q(1)  # B=1+(1+C)/(t-1)>1; its global infimum is 1 at infinity.
B_at_T_lo=DminB/(T-Q(1))
print('tail positivity lower bounds A',dec(AloB,lower=True),'t+C',dec(DminB,lower=True),'U',dec(Ulo,lower=True),'B',dec(Blo,lower=True))
print('B(T) lower bound only',dec(B_at_T_lo,lower=True))
assert AloB>0 and DminB>0 and Ulo>0 and Blo>0
print('PASS: unique exact nonlinear tail uniformly over the Poincare tail box')
print('PASS: A>0, t+C>0, 1-t A\'/A>0 throughout t>=T')
