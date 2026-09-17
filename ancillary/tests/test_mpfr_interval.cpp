#include "mpfr_interval.hpp"
#include <gmpxx.h>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>
using namespace mpi;

static mpq_class exact(const Real& x) {
    mpq_class q; mpfr_get_q(q.get_mpq_t(), x.x); return q;
}
static I rational_interval(mpq_class a, mpq_class b) {
    if (b<a) std::swap(a,b);
    Real l,h;
    mpfr_set_q(l.x,a.get_mpq_t(),MPFR_RNDD);
    mpfr_set_q(h.x,b.get_mpq_t(),MPFR_RNDU);
    return I(l,h);
}
static void encloses(const I& x, const mpq_class& l, const mpq_class& h) {
    assert(exact(x.lo)<=l); assert(exact(x.hi)>=h);
}
static mpq_class decimal(const std::string& s) {
    const auto ep=s.find_first_of("eE");
    std::string digits=s.substr(0,ep);
    long exp=ep==std::string::npos?0:std::stol(s.substr(ep+1));
    auto dot=digits.find('.');
    if (dot!=std::string::npos) {
        exp-=static_cast<long>(digits.size()-dot-1); digits.erase(dot,1);
    }
    mpz_class n(digits,10),power;
    mpz_ui_pow_ui(power.get_mpz_t(),10,static_cast<unsigned long>(std::abs(exp)));
    mpq_class q = exp>=0 ? mpq_class(n*power) : mpq_class(n,power);
    q.canonicalize(); return q;
}
static void check_serialization(const I& x) {
    for (int digits : {1,2,7,30,50}) {
        std::string s=str(x,digits);
        auto comma=s.find(',');
        auto l=decimal(s.substr(1,comma-1));
        auto h=decimal(s.substr(comma+1,s.size()-comma-2));
        assert(l<=exact(x.lo)); assert(h>=exact(x.hi));
        assert(decimal(str(x.lo,digits,MPFR_RNDD))<=exact(x.lo));
        assert(decimal(str(x.hi,digits,MPFR_RNDU))>=exact(x.hi));
    }
}
int main() {
    std::mt19937_64 rng(0x45575f4d504652ULL);
    unsigned checked=0;
    for (mpfr_prec_t precision : {24,53,113,512}) {
        PREC=precision;
        for (int iteration=0; iteration<350; ++iteration) {
            auto rq=[&](){ mpq_class q(static_cast<long>(rng()%4001)-2000,
                                      static_cast<unsigned long>(rng()%1000)+1);
                          q.canonicalize(); return q; };
            I a=rational_interval(rq(),rq()), b=rational_interval(rq(),rq());
            const mpq_class al=exact(a.lo),ah=exact(a.hi),bl=exact(b.lo),bh=exact(b.hi);
            encloses(a+b,al+bl,ah+bh);
            encloses(a-b,al-bh,ah-bl);
            encloses(-a,-ah,-al);
            std::vector<mpq_class> products={al*bl,al*bh,ah*bl,ah*bh};
            encloses(a*b,*std::min_element(products.begin(),products.end()),
                         *std::max_element(products.begin(),products.end()));
            if (bl>0 || bh<0) {
                std::vector<mpq_class> quotients={al/bl,al/bh,ah/bl,ah/bh};
                encloses(a/b,*std::min_element(quotients.begin(),quotients.end()),
                             *std::max_element(quotients.begin(),quotients.end()));
            } else {
                bool threw=false; try { (void)(a/b); } catch(const std::domain_error&) { threw=true; }
                assert(threw);
            }
            assert(exact(width(a))>=ah-al);
            assert(exact(mid(a))==(al+ah)/2);
            I alias=a; alias+=alias; encloses(alias,2*al,2*ah);
            alias=a; alias*=alias;
            products={al*al,al*ah,ah*al,ah*ah};
            encloses(alias,*std::min_element(products.begin(),products.end()),
                           *std::max_element(products.begin(),products.end()));
            check_serialization(a);
            ++checked;
        }
    }
    PREC=512;
    for (const char* s : {"0.1","-0.1","1e-100","-1e-100","9.9999999e100"}) {
        I x(s); encloses(x,decimal(s),decimal(s)); check_serialization(x);
    }
    // Regression: reciprocal of a negative interval is decreasing too.
    I negative_div=I(1)/I(Real(-4),Real(-2));
    encloses(negative_div,mpq_class(-1,2),mpq_class(-1,4));
    assert(exact(negative_div.lo)==mpq_class(-1,2));
    assert(exact(negative_div.hi)==mpq_class(-1,4));
    // Exact midpoint at a large exponent separation and across zero.
    I gap(Real("1e-200",MPFR_RNDD),Real(1));
    assert(exact(mid(gap))==(exact(gap.lo)+exact(gap.hi))/2);
    I crossing(Real(-1),Real("1e-200",MPFR_RNDU));
    assert(exact(mid(crossing))==(exact(crossing.lo)+exact(crossing.hi))/2);
    Real high("0.123456789012345678901234567890123456789");
    auto q=exact(high); auto p=mpfr_get_prec(high.x);
    PREC=24;
    Real copied=high, assigned; assigned=high;
    assert(exact(copied)==q && exact(assigned)==q);
    assert(mpfr_get_prec(copied.x)==p && mpfr_get_prec(assigned.x)==p);
    high=high; assert(exact(high)==q);
    Real moved=std::move(copied); assert(exact(moved)==q);
    assigned=std::move(moved); assert(exact(assigned)==q);
    for (const char* bad : {"garbage","nan","inf","1.2junk",""}) {
        bool threw=false; try { I x(bad); } catch(const std::invalid_argument&) { threw=true; }
        assert(threw);
    }
    bool reversed=false;
    try { I x(Real(2),Real(1)); } catch(const std::invalid_argument&) { reversed=true; }
    assert(reversed);
    PREC=512;
    std::cout << "PASS: " << checked << " randomized interval pairs; exact rational arithmetic oracle, "
              << "negative division, aliasing, precision changes, exact midpoint, directed decimal parsing/serialization.\n";
    std::cout << "MPFR " << mpfr_get_version() << ", GMP " << gmp_version << "\n";
}
