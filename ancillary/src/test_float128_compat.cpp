#include <mpfr.h>
#include <quadmath.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <random>
extern "C" int mpfr_set_float128(mpfr_ptr, __float128, mpfr_rnd_t);
extern "C" __float128 mpfr_get_float128(mpfr_srcptr, mpfr_rnd_t);
using U128=unsigned __int128;
static U128 bits(__float128 x){U128 u;std::memcpy(&u,&x,16);return u;}
static __float128 value(U128 u){__float128 x;std::memcpy(&x,&u,16);return x;}
int main(){
 mpfr_t x,lo,hi,nxt,tmp;mpfr_inits2(512,x,lo,hi,nxt,tmp,(mpfr_ptr)0);
 std::mt19937_64 rng(729411);
 // Exact round-trip of every class of finite IEEE binary128 encoding.
 for(int i=0;i<10000;i++){
  U128 u=(U128(rng())<<64)|rng();
  if(((u>>112)&0x7fff)==0x7fff)u^=U128(1)<<112;
  __float128 a=value(u);assert(mpfr_set_float128(x,a,MPFR_RNDN)==0);
  assert(bits(mpfr_get_float128(x,MPFR_RNDN))==u);
 }
 const int es[]={-20000,-16495,-16494,-16493,-16383,-16382,-113,-1,0,1,112,16382,16383,16384,20000};
 for(int e:es)for(int j=0;j<100;j++)for(int sign:{-1,1}){
  mpfr_set_ui(x,rng(),MPFR_RNDN);mpfr_div_2ui(x,x,64,MPFR_RNDN);
  mpfr_add_ui(x,x,1,MPFR_RNDN);mpfr_mul_2si(x,x,e,MPFR_RNDN);
  // Add information below the 113-bit significand before directed conversion.
  mpfr_set_ui_2exp(tmp,1,e-160,MPFR_RNDN);mpfr_add(x,x,tmp,MPFR_RNDN);
  if(sign<0)mpfr_neg(x,x,MPFR_RNDN);
  auto a=mpfr_get_float128(x,MPFR_RNDD),b=mpfr_get_float128(x,MPFR_RNDU);
  mpfr_set_float128(lo,a,MPFR_RNDN);mpfr_set_float128(hi,b,MPFR_RNDN);
  assert(mpfr_cmp(lo,x)<=0&&mpfr_cmp(hi,x)>=0);
  if(!isinfq(a)){mpfr_set_float128(nxt,nextafterq(a,HUGE_VALQ),MPFR_RNDN);assert(mpfr_cmp(nxt,x)>=0);}
  if(!isinfq(b)){mpfr_set_float128(nxt,nextafterq(b,-HUGE_VALQ),MPFR_RNDN);assert(mpfr_cmp(nxt,x)<=0);}
 }
 // Known values, signed zeros, and both normal/subnormal tie-to-even cases.
 mpfr_set_float128(x,1.0Q,MPFR_RNDN);assert(mpfr_cmp_ui(x,1)==0);
 mpfr_set_float128(x,-0.0Q,MPFR_RNDN);assert(mpfr_zero_p(x)&&mpfr_signbit(x));
 mpfr_set_ui(x,1,MPFR_RNDN);mpfr_set_ui_2exp(tmp,1,-113,MPFR_RNDN);mpfr_add(x,x,tmp,MPFR_RNDN);
 assert(mpfr_get_float128(x,MPFR_RNDN)==1.0Q);
 mpfr_set_ui_2exp(x,1,-16495,MPFR_RNDN);assert(bits(mpfr_get_float128(x,MPFR_RNDN))==0);
 assert(bits(mpfr_get_float128(x,MPFR_RNDU))==1);
 mpfr_clears(x,lo,hi,nxt,tmp,(mpfr_ptr)0);
 puts("PASS binary128 bridge: 10000 exact round trips; 3000 signed directed cases; boundary/tie checks");
}
