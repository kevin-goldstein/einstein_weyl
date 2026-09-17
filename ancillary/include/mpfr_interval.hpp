#ifndef EW_RECONSTRUCTED_MPFR_INTERVAL_HPP
#define EW_RECONSTRUCTED_MPFR_INTERVAL_HPP

// Independently written compatibility implementation, 2026-09-14.
// This is not the missing author's header. Only finite, closed intervals are
// supported; undefined operations throw instead of returning uncertified data.
// Requires C++17, MPFR, and GMP. Do not compile with unsafe math optimizations.
#include <mpfr.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstddef>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>

namespace mpi {

inline mpfr_prec_t PREC = 512;

inline mpfr_prec_t checked_precision(mpfr_prec_t p) {
    if (p < MPFR_PREC_MIN || p > MPFR_PREC_MAX)
        throw std::invalid_argument("invalid MPFR precision");
    return p;
}

// A Real is one exact represented dyadic number, not an interval. Decimal
// scalar construction rounds nearest by default; callers that need a bound
// must supply a direction or construct I from the decimal text instead.
struct Real {
    mpfr_t x;
    Real() { mpfr_init2(x, checked_precision(PREC)); mpfr_set_zero(x, 1); }
    explicit Real(mpfr_prec_t p, bool) {
        mpfr_init2(x, checked_precision(p)); mpfr_set_zero(x, 1);
    }
    Real(int n) : Real(static_cast<long>(n)) {}
    Real(long n) : Real(std::max<mpfr_prec_t>(PREC, sizeof(long)*CHAR_BIT), true) {
        mpfr_set_si(x, n, MPFR_RNDN);  // exact at this precision
    }
    explicit Real(double n) : Real(std::max<mpfr_prec_t>(PREC, 53), true) {
        if (!std::isfinite(n)) throw std::domain_error("nonfinite Real");
        mpfr_set_d(x, n, MPFR_RNDN);  // exact binary64 input
    }
    Real(const char* s, mpfr_rnd_t rnd = MPFR_RNDN) : Real() {
        if (!s || mpfr_set_str(x, s, 10, rnd) != 0 || !mpfr_number_p(x)) {
            throw std::invalid_argument("invalid finite decimal Real");
        }
    }
    Real(const std::string& s, mpfr_rnd_t rnd = MPFR_RNDN) : Real(s.c_str(), rnd) {}
    Real(const Real& v) : Real(mpfr_get_prec(v.x), true) { mpfr_set(x, v.x, MPFR_RNDN); }
    Real(Real&& v) noexcept {
        mpfr_init2(x, MPFR_PREC_MIN); mpfr_set_zero(x, 1); mpfr_swap(x, v.x);
    }
    Real& operator=(const Real& v) {
        if (this != &v) {
            mpfr_set_prec(x, mpfr_get_prec(v.x));
            mpfr_set(x, v.x, MPFR_RNDN);
        }
        return *this;
    }
    Real& operator=(Real&& v) noexcept {
        if (this != &v) mpfr_swap(x, v.x);
        return *this;
    }
    ~Real() { mpfr_clear(x); }
};

inline void require_finite(const Real& a) {
    if (!mpfr_number_p(a.x)) throw std::domain_error("nonfinite MPFR result");
}
inline bool lt(const Real& a, const Real& b) {
    require_finite(a); require_finite(b); return mpfr_less_p(a.x, b.x) != 0;
}
inline bool le(const Real& a, const Real& b) {
    require_finite(a); require_finite(b); return mpfr_lessequal_p(a.x, b.x) != 0;
}
inline Real radd(const Real& a, const Real& b, mpfr_rnd_t rnd) {
    Real r; mpfr_add(r.x, a.x, b.x, rnd); require_finite(r); return r;
}
inline Real rsub(const Real& a, const Real& b, mpfr_rnd_t rnd) {
    Real r; mpfr_sub(r.x, a.x, b.x, rnd); require_finite(r); return r;
}
inline Real rmul(const Real& a, const Real& b, mpfr_rnd_t rnd) {
    Real r; mpfr_mul(r.x, a.x, b.x, rnd); require_finite(r); return r;
}
inline Real rdiv(const Real& a, const Real& b, mpfr_rnd_t rnd) {
    if (mpfr_zero_p(b.x)) throw std::domain_error("division by zero Real");
    Real r; mpfr_div(r.x, a.x, b.x, rnd); require_finite(r); return r;
}
inline Real rneg(const Real& a, mpfr_rnd_t rnd) {
    Real r; mpfr_neg(r.x, a.x, rnd); require_finite(r); return r;
}
inline double tod(const Real& a, mpfr_rnd_t rnd = MPFR_RNDN) {
    require_finite(a); return mpfr_get_d(a.x, rnd);
}

struct I {
    Real lo, hi;
    I() : lo(0), hi(0) {}
    I(int n) : lo(n), hi(n) {}
    I(long n) : lo(n), hi(n) {}
    explicit I(double n) : lo(n), hi(n) {}
    I(const char* s) : lo(s, MPFR_RNDD), hi(s, MPFR_RNDU) { validate(); }
    I(const std::string& s) : I(s.c_str()) {}
    I(const Real& a) : lo(a), hi(a) { validate(); }
    I(const Real& a, const Real& b) : lo(a), hi(b) { validate(); }
    void validate() const {
        require_finite(lo); require_finite(hi);
        if (!le(lo, hi)) throw std::invalid_argument("reversed interval endpoints");
    }
};

inline I operator+(const I& a, const I& b) {
    return I(radd(a.lo, b.lo, MPFR_RNDD), radd(a.hi, b.hi, MPFR_RNDU));
}
inline I operator-(const I& a, const I& b) {
    return I(rsub(a.lo, b.hi, MPFR_RNDD), rsub(a.hi, b.lo, MPFR_RNDU));
}
inline I operator-(const I& a) {
    return I(rneg(a.hi, MPFR_RNDD), rneg(a.lo, MPFR_RNDU));
}
inline I operator*(const I& a, const I& b) {
    Real lower[4] = {rmul(a.lo,b.lo,MPFR_RNDD), rmul(a.lo,b.hi,MPFR_RNDD),
                     rmul(a.hi,b.lo,MPFR_RNDD), rmul(a.hi,b.hi,MPFR_RNDD)};
    Real upper[4] = {rmul(a.lo,b.lo,MPFR_RNDU), rmul(a.lo,b.hi,MPFR_RNDU),
                     rmul(a.hi,b.lo,MPFR_RNDU), rmul(a.hi,b.hi,MPFR_RNDU)};
    int il = 0, ih = 0;
    for (int j=1; j<4; ++j) {
        if (lt(lower[j], lower[il])) il = j;
        if (lt(upper[ih], upper[j])) ih = j;
    }
    return I(lower[il], upper[ih]);
}
inline I operator/(const I& a, const I& b) {
    if (mpfr_sgn(b.lo.x) <= 0 && mpfr_sgn(b.hi.x) >= 0)
        throw std::domain_error("interval divisor contains zero");
    // 1/x is decreasing on each component of R\{0}, including x<0.
    I reciprocal(rdiv(Real(1), b.hi, MPFR_RNDD),
                 rdiv(Real(1), b.lo, MPFR_RNDU));
    return a * reciprocal;
}
inline I& operator+=(I& a, const I& b) { a = a+b; return a; }
inline I& operator-=(I& a, const I& b) { a = a-b; return a; }
inline I& operator*=(I& a, const I& b) { a = a*b; return a; }
inline I& operator/=(I& a, const I& b) { a = a/b; return a; }

inline Real width(const I& a) { return rsub(a.hi, a.lo, MPFR_RNDU); }

// The exact midpoint avoids the rounded-midpoint/half-width enclosure bug in
// callers that recenter intervals. Allocate sufficient temporary precision to
// align the two finite binary significands, add exactly, and divide by two.
// Ordinary directed operations remain at PREC, so this does not cause
// permanent precision growth in the ODE recurrence.
inline Real mid(const I& a) {
    a.validate();
    if (mpfr_zero_p(a.lo.x) && mpfr_zero_p(a.hi.x)) return Real(0);
    const mpfr_exp_t el = mpfr_zero_p(a.lo.x) ? mpfr_get_exp(a.hi.x) : mpfr_get_exp(a.lo.x);
    const mpfr_exp_t eh = mpfr_zero_p(a.hi.x) ? el : mpfr_get_exp(a.hi.x);
    const long double needed = static_cast<long double>(std::max(mpfr_get_prec(a.lo.x), mpfr_get_prec(a.hi.x)))
                             + std::fabs(static_cast<long double>(el)-static_cast<long double>(eh)) + 2;
    if (needed > static_cast<long double>(MPFR_PREC_MAX))
        throw std::overflow_error("exact midpoint requires excessive precision");
    Real r(static_cast<mpfr_prec_t>(needed), true);
    if (mpfr_add(r.x, a.lo.x, a.hi.x, MPFR_RNDN) != 0)
        throw std::logic_error("midpoint addition unexpectedly inexact");
    if (mpfr_div_2ui(r.x, r.x, 1, MPFR_RNDN) != 0)
        throw std::underflow_error("midpoint underflow");
    require_finite(r); return r;
}

// Decimal text is generated by MPFR itself with an explicit direction. The
// returned scientific notation is accepted exactly by Python Fraction.
inline std::string str(const Real& a, int digits = 40, mpfr_rnd_t rnd = MPFR_RNDN) {
    require_finite(a);
    if (digits < 1) throw std::invalid_argument("decimal digits must be positive");
    if (mpfr_zero_p(a.x)) return "0";
    mpfr_exp_t exponent;
    char* raw = mpfr_get_str(nullptr, &exponent, 10, static_cast<std::size_t>(digits), a.x, rnd);
    if (!raw) throw std::bad_alloc();
    std::string s(raw); mpfr_free_str(raw);
    const bool negative = !s.empty() && s.front() == '-';
    if (negative) s.erase(s.begin());
    std::string out = negative ? "-" : "";
    out += s.front();
    if (s.size()>1) { out += '.'; out += s.substr(1); }
    out += 'e'; out += std::to_string(exponent-1);
    return out;
}
inline std::string str(const I& a, int digits = 40) {
    a.validate();
    return "[" + str(a.lo, digits, MPFR_RNDD) + "," + str(a.hi, digits, MPFR_RNDU) + "]";
}

} // namespace mpi
#endif
