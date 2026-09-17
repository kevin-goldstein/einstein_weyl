// Focused tests call the actual production helpers. Rename only its entry point
// so including the source cannot run the complete certificate as a side effect.
#define main core_certificate_main_for_tests
#include "../src/verified_core_taylor_model.repaired.cpp"
#undef main
#include <cassert>

static Real exact_q(__float128 x) {
  char text[256];
  quadmath_snprintf(text, sizeof(text), "%.40Qa", x);
  Real out;
  assert(mpfr_set_str(out.x, text, 0, MPFR_RNDN) == 0);
  return out;
}

int main() {
  QuadInterval a(1.0Q, qup(1.0Q), true);
  __float128 center = (a.lo + a.hi) / 2;
  __float128 old_radius = qup(0.5Q * qup(a.hi - a.lo));
  __float128 radius = qradius_about(a, center);
  // An adjacent-endpoint interval has a rounded midpoint at one endpoint.
  assert(center == a.lo);
  assert(old_radius < a.hi - center);
  assert(radius >= a.hi - center);
  assert(qdn(center - radius) <= a.lo && qup(center + radius) >= a.hi);

  for (int denominator = 2; denominator <= 100; denominator++) {
    QuadInterval b(-(__float128)denominator - 1, -(__float128)denominator, true);
    QuadInterval inv = QuadInterval(1) / b;
    assert(inv.lo <= inv.hi);
    Real exact_lo = rdiv(Real(1), Real(-denominator), MPFR_RNDD);
    Real exact_hi = rdiv(Real(1), Real(-denominator - 1), MPFR_RNDU);
    assert(le(exact_q(inv.lo), exact_lo));
    assert(le(exact_hi, exact_q(inv.hi)));
  }

  for (const char *s : {"1.23456789", "-1.23456789", "0.00000000000123456789"}) {
    Real value(s);
    std::string printed = mdownstr(value, 5);
    Real as_decimal(printed, MPFR_RNDU);
    assert(le(as_decimal, value));
  }

  // The actual tube predicate must reject equality at either boundary.
  const QuadInterval tube(-1.0Q, 1.0Q, true);
  assert(strict_quad_subset(QuadInterval(-0.5Q, 0.5Q, true), tube));
  assert(!strict_quad_subset(QuadInterval(-1.0Q, 0.5Q, true), tube));
  assert(!strict_quad_subset(QuadInterval(-0.5Q, 1.0Q, true), tube));
  assert(!strict_quad_subset(tube, tube));
  assert(!strict_quad_subset(QuadInterval(-1.5Q, 0.5Q, true), tube));
  assert(!strict_quad_subset(QuadInterval(-0.5Q, 1.5Q, true), tube));
  assert(strict_quad_subset(QuadInterval(0), tube));
  assert(!strict_quad_subset(QuadInterval(0), QuadInterval(0)));

  std::cout << "PASS: production helpers: rounded-center radius, 99 negative "
               "reciprocal enclosures, directed MPFR minima, 8 strict tube "
               "boundary cases\n";
}
