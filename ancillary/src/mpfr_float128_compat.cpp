// Reconstructed numeric bridge for MPFR builds without binary128 conversion.
// Requires GCC IEEE binary128 and a 128-bit unsigned integer. This implements
// the numeric conversions used by this certificate; MPFR status flags are not
// emulated. All supported directed roundings are performed by MPFR/GMP.
#include <mpfr.h>
#include <quadmath.h>
#include <cstdint>
#include <cstring>
#include <algorithm>

static_assert(sizeof(__float128) == 16, "IEEE binary128 required");
static_assert(__FLT128_MANT_DIG__ == 113 && __FLT128_MAX_EXP__ == 16384,
              "IEEE binary128 exponent and significand required");
using U128 = unsigned __int128;
static constexpr U128 frac_mask = (U128(1) << 112) - 1;
static __float128 from_bits(U128 u) {
    __float128 f;
    std::memcpy(&f, &u, sizeof(f));
    return f;
}
static U128 to_bits(__float128 f) {
    U128 u;
    std::memcpy(&u, &f, sizeof(u));
    return u;
}
static void set_u128(mpz_ptr z, U128 u) {
    std::uint64_t words[2] = {std::uint64_t(u), std::uint64_t(u >> 64)};
    mpz_import(z, 2, -1, sizeof(std::uint64_t), 0, 0, words);
}
static U128 get_u128(mpz_srcptr z) {
    std::uint64_t words[2] = {0, 0};
    size_t count = 0;
    mpz_export(words, &count, -1, sizeof(std::uint64_t), 0, 0, z);
    return U128(words[0]) | (U128(words[1]) << 64);
}
extern "C" int mpfr_set_float128(mpfr_ptr dest, __float128 value,
                                  mpfr_rnd_t rnd) {
    const U128 u = to_bits(value);
    const bool neg = (u >> 127) != 0;
    const unsigned exp = unsigned((u >> 112) & 0x7fff);
    const U128 frac = u & frac_mask;
    if (exp == 0x7fff) {
        if (frac) mpfr_set_nan(dest);
        else mpfr_set_inf(dest, neg ? -1 : 1);
        return 0;
    }
    if (exp == 0 && frac == 0) {
        mpfr_set_zero(dest, neg ? -1 : 1);
        return 0;
    }
    mpz_t z;
    mpz_init(z);
    set_u128(z, frac | (exp ? (U128(1) << 112) : U128(0)));
    if (neg) mpz_neg(z, z);
    const mpfr_exp_t quantum = exp ? mpfr_exp_t(exp) - 16383 - 112 : -16494;
    const int status = mpfr_set_z_2exp(dest, z, quantum, rnd);
    mpz_clear(z);
    return status;
}
static __float128 overflow_value(bool neg, mpfr_rnd_t rnd) {
    const bool infinity = rnd == MPFR_RNDN || rnd == MPFR_RNDA ||
        rnd == MPFR_RNDF || (rnd == MPFR_RNDU && !neg) ||
        (rnd == MPFR_RNDD && neg);
    const U128 sign = U128(neg) << 127;
    return from_bits(sign | (infinity ? U128(0x7fff) << 112 :
                             (U128(0x7ffe) << 112) | frac_mask));
}
extern "C" __float128 mpfr_get_float128(mpfr_srcptr source, mpfr_rnd_t rnd) {
    const bool neg = mpfr_signbit(source) != 0;
    const U128 sign = U128(neg) << 127;
    if (mpfr_nan_p(source)) return from_bits((U128(0x7fff) << 112) | (U128(1) << 111));
    if (mpfr_inf_p(source)) return from_bits(sign | (U128(0x7fff) << 112));
    if (mpfr_zero_p(source)) return from_bits(sign);
    mpfr_exp_t exponent = mpfr_get_exp(source) - 1;
    if (exponent > 16383) return overflow_value(neg, rnd);
    const mpfr_exp_t quantum = std::max<mpfr_exp_t>(exponent - 112, -16494);
    mpfr_t scaled;
    mpfr_init2(scaled, mpfr_get_prec(source));
    // Exact power-of-two scaling within MPFR's much wider default exponent range.
    mpfr_mul_2si(scaled, source, -quantum, MPFR_RNDN);
    mpz_t z;
    mpz_init(z);
    mpfr_get_z(z, scaled, rnd == MPFR_RNDF ? MPFR_RNDN : rnd);
    mpz_abs(z, z);
    U128 sig = get_u128(z);
    mpz_clear(z);
    mpfr_clear(scaled);
    if (sig == 0) return from_bits(sign);
    if (exponent < -16382) {
        // The integer significand includes a possible carry into the minimum normal.
        return from_bits(sign | sig);
    }
    if (sig >= (U128(1) << 113)) {
        sig >>= 1;
        ++exponent;
    }
    if (exponent > 16383) return overflow_value(neg, rnd);
    return from_bits(sign | (U128(exponent + 16383) << 112) | (sig & frac_mask));
}
