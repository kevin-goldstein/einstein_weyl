#!/usr/bin/env python3
"""Exact symbolic checks of the Einstein--Weyl constraint and its coordinates.

Requires SymPy. No numerical orbit or approximate algebra enters these checks.
The script fails unless all displayed residuals simplify identically to zero.
"""
import sympy as s

t,A,p,C,d,Y,J,kappa,R=s.symbols('t A p C d Y J kappa R',nonzero=True)
a0,v0=s.symbols('a0 v0')
F=2*p**2/A+A*Y/(6*(t+C))+A*d/(3*t*(t+C))-p*d/(t+C)-C*p/(t*(t+C))-2*p/(t+C)
N=(2*kappa*t**2*A**4*d*J-kappa*t**2*A**4*Y**2-4*kappa*t**2*A**4*J-6*kappa*t*A**4*C*J
   +2*kappa*t*A**4*d*Y+8*kappa*t*A**4*Y+6*kappa*A**4*C*Y-4*kappa*A**4*d**2-8*kappa*A**4*d
   +18*t**4*p**2-6*t**3*A*p*d-24*t**3*A*p+18*t**3*C*p**2+6*t**2*A**2*d-18*t**2*A*C*p)

def verify(name,residual):
    result=s.factor(s.cancel(residual))
    assert result==0, f'{name}: nonzero residual {result}'
    print('PASS:',name,'(exact residual 0)')

states=(A,p,C,d,Y,J)
flow=(p,F,d,Y,J,3*t*F/(kappa*A**3))
Ndot=s.diff(N,t)+sum(s.diff(N,x)*fx for x,fx in zip(states,flow))
verify('N derivative equals 4 p N / A',Ndot-4*p*N/A)
E2=-N/(6*A**4)
verify('E2 is constant along the six-state equations',
       s.diff(E2,t)+sum(s.diff(E2,x)*fx for x,fx in zip(states,flow)))

cJ=2*kappa*t*A**4*(t*(d-2)-3*C)
fixed_v_derivative=(-2*kappa*t**2*A**4*Y+a0*cJ+2*kappa*t*A**4*d
                    +8*kappa*t*A**4+6*kappa*A**4*C)
verify('fixed-v derivative used by the scalar matching certificate',
       s.diff(N.subs(J,v0+a0*Y),Y)-fixed_v_derivative)

# Dimensionless Kundt coordinates: Omega=W, H=-z Q, z=R(1-1/t).
# W=t/A and Q=R(t+C)/(t^2(t-1)), so H=-R^2(t+C)/t^3.
Af=s.Function('A')(t);Cf=s.Function('C')(t)
W=t/Af
H=-R**2*(t+Cf)/t**3
Dz=lambda expression:t**2/R*s.diff(expression,t)
Hz=Dz(H);Hzz=Dz(Hz);Hzzz=Dz(Hzz)
kundt_residual=W*Dz(W)*Hz+3*Dz(W)**2*H+W**2-kappa/s.Integer(3)*(Hz*Hzzz-Hzz**2/2+2)
substitution={Af:A,s.diff(Af,t):p,Cf:C,s.diff(Cf,t):d,s.diff(Cf,t,2):Y,s.diff(Cf,t,3):J}
verify('remaining Kundt equation in t coordinates equals -N/(6 A^4)',
       kundt_residual.subs(substitution)-E2)

# Analytic horizon initialization, before dividing by the vanishing t+C.
b,rho=s.symbols('b rho',nonzero=True)
gamma2=(4-rho**2+3*b)/3
H0=0;Hz0=-1;Hzz0=4+6*b;Hzzz0=-6*(1+3*b*gamma2)
W0=1;Wz0=1+b
horizon_residual=(W0*Wz0*Hz0+3*Wz0**2*H0+W0**2
                  -(Hz0*Hzzz0-Hzz0**2/2+2)/(6*rho**2))
verify('general-radius horizon data satisfy the remaining Kundt equation',horizon_residual)
print('All symbolic identities verified with SymPy',s.__version__)
