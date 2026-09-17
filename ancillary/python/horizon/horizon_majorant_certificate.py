#!/usr/bin/env python3
"""Exact-rational certificate for convergence of the Einstein--Weyl horizon series.

For rho_h=1 and complex b with |b| <= 37/100, the direct W,Q coefficient
recurrence is dominated by
    |w_n| <= 10^n/(n+1)^2,
    |q_n| <= 5*10^n/(n+1)^4.
All finite checks and all displayed inequalities use fractions.Fraction.
"""
from fractions import Fraction as Q

L = 10
B = 5
N0 = 100


def q_ratio(n: int) -> Q:
    """Majorant ratio for q_{n+3} / [B L^(n+3)/(n+4)^4]."""
    source = Q(0)
    for i in range(n + 1):
        source += Q((n-i+2)*(n-i+1), (i+1)**2 * (n-i+3)**2)
        source += Q(2*(i+1)*(n-i+1), (i+2)**2 * (n-i+2)**2)
    return Q(6, B*L) * Q((n+4)**3, (n+1)*(n+2)*(n+3)) * source


def w_parts(n: int):
    """Six positive majorant contributions for w_{n+1}/U_{n+1}."""
    A = Q(0)  # z Q W''
    B1 = Q(0) # Q W'
    C = Q(0)  # z Q' W'
    D = Q(1, 3*L*(n+1)**2)
    E = Q(0)  # +(1/3) Q' W
    F = Q(0)  # +(1/6) z Q'' W

    for i in range(1, n):
        j = n - 1 - i
        A += Q(B*(j+2)*(j+1), (i+1)**4 * (j+3)**2)
    for i in range(1, n+1):
        j = n - i
        base = Q(B*(j+1), (i+1)**4 * (j+2)**2)
        B1 += base
        C += i * base
    for i in range(n+1):
        E += Q(B*(i+1), 3*(i+2)**4 * (n-i+1)**2)
    for i in range(n):
        F += Q(B*(i+2)*(i+1), 6*(i+3)**4 * (n-i)**2)

    factor = Q((n+2)**2, (n+1)**2)
    return tuple(factor*x for x in (A, B1, C, D, E, F))


def w_ratio(n: int) -> Q:
    return sum(w_parts(n), Q(0))


# Initial coefficients, uniformly for |b| <= 37/100.
# w0=1, w1=1+b; q0=1, q1=-2-3b, q2=1+3b+3b^2.
assert Q(1) <= Q(1)
assert Q(137,100) <= Q(L,4)
assert Q(1) <= Q(B)
assert Q(311,100) <= Q(B*L,16)
assert Q(25207,10000) <= Q(B*L*L,81)

# Exact finite induction checks.
max_q = (Q(0), -1)
max_w = (Q(0), -1)
for n in range(N0):
    qr = q_ratio(n)
    if qr > max_q[0]:
        max_q = (qr, n)
    assert qr < 1
for n in range(1, N0):
    wr = w_ratio(n)
    if wr > max_w[0]:
        max_w = (wr, n)
    assert wr < 1

# Analytic q-bound for n >= 100.
# H_{n+2}-1 <= log(n+2) <= (n+4)/10, hence
# source <= 5/3 + 2/5 = 31/15.
def q_asymptotic_bound(n: int) -> Q:
    return Q(6, B*L) * Q((n+4)**3, (n+1)*(n+2)*(n+3)) * Q(31,15)

assert q_asymptotic_bound(N0) < 1
# The rational prefactor is decreasing for n >= 0, so n=N0 is the maximum.

# Analytic w-bound for n >= 100.  The six convolution pieces are bounded by
# split sums and zeta estimates zeta(2)<5/3, zeta(3)-1<1/4,
# zeta(4)-1<1/12.  Also 1+log(n+1) <= (n+2)/10 for n>=100.
def w_asymptotic_bound(n: int) -> Q:
    inner = (
        Q(5,12)
        + Q(5, 6*(n+4)) + Q(8, (n+2)**3)
        + Q(5, 2*(n+4)) + Q(4, (n+2)**2)
        + Q(1, 30*(n+1)**2)
        + Q(5, 3*(n+2)**2) + Q(200, 9*(n+2)**3)
        + Q(100, 9*(n+1)**2)
    )
    return Q((n+2)**2, (n+1)**2) * inner

assert w_asymptotic_bound(N0) < 1
# Every term in the inner bound and the prefactor is decreasing, so n=N0
# controls all n>=N0.

# Explicit geometric remainder bounds on |z|<=1/20 (10|z|<=1/2).
# For k=0,1,2,3,4 and n>N, n^(underline k)/(n+1)^2 <= (n+1)^(k-2).
# We bound the polynomial-geometric tail by direct exact summation through
# N+500 plus a final ratio estimate.
def derivative_tail_w(N: int, k: int, z: Q=Q(1,20)) -> Q:
    rho = L*z
    total = Q(0)
    last = Q(0)
    for n in range(N+1, N+501):
        ff = 1
        for j in range(k):
            ff *= n-j
        term = Q(ff * (L**n), (n+1)**2) * (z**(n-k))
        total += term
        last = term
    # for n beyond this point, consecutive ratios are < 3/5 for all k<=4
    assert k <= 4
    total += last * Q(3,2)  # conservative tail after the last included term
    return total

for k in range(5):
    assert derivative_tail_w(80,k) < Q(1,10**12)

print("Einstein--Weyl horizon majorant certificate")
print("all arithmetic: exact rational")
print(f"finite q maximum: {max_q[0]} at n={max_q[1]} = {float(max_q[0]):.12f}")
print(f"finite w maximum: {max_w[0]} at n={max_w[1]} = {float(max_w[0]):.12f}")
print(f"q asymptotic bound at n=100: {float(q_asymptotic_bound(100)):.12f}")
print(f"w asymptotic bound at n=100: {float(w_asymptotic_bound(100)):.12f}")
for k in range(5):
    t = derivative_tail_w(80,k)
    print(f"W derivative order {k} tail, N=80, |z|<=0.05: < {float(t):.6e}")
print("PASS: the horizon series converges absolutely for |z|<0.1, uniformly for |b|<=0.37")
