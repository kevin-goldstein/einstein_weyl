#include "mpfr_interval.hpp"
#include <array>
#include <vector>
#include <iostream>
using namespace mpi;
static I inflate_abs(const I&a,const char*e){return I(rsub(a.lo,Real(e),MPFR_RNDD),radd(a.hi,Real(e),MPFR_RNDU));}
struct D2{
 I v; std::array<I,2> g; std::array<std::array<I,2>,2> h;
 D2():v(0){for(auto&x:g)x=I(0);for(auto&r:h)for(auto&x:r)x=I(0);} D2(int x):D2(){v=I(x);} D2(long x):D2(){v=I(x);} D2(const I&x):D2(){v=x;} D2(const char*x):D2(){v=I(x);} };
static D2 operator+(const D2&a,const D2&b){D2 z;z.v=a.v+b.v;for(int i=0;i<2;i++){z.g[i]=a.g[i]+b.g[i];for(int j=0;j<2;j++)z.h[i][j]=a.h[i][j]+b.h[i][j];}return z;}
static D2 operator-(const D2&a,const D2&b){D2 z;z.v=a.v-b.v;for(int i=0;i<2;i++){z.g[i]=a.g[i]-b.g[i];for(int j=0;j<2;j++)z.h[i][j]=a.h[i][j]-b.h[i][j];}return z;}
static D2 operator-(const D2&a){D2 z;z.v=-a.v;for(int i=0;i<2;i++){z.g[i]=-a.g[i];for(int j=0;j<2;j++)z.h[i][j]=-a.h[i][j];}return z;}
static D2 operator*(const D2&a,const D2&b){D2 z;z.v=a.v*b.v;for(int i=0;i<2;i++){z.g[i]=a.g[i]*b.v+a.v*b.g[i];for(int j=0;j<2;j++)z.h[i][j]=a.h[i][j]*b.v+a.g[i]*b.g[j]+a.g[j]*b.g[i]+a.v*b.h[i][j];}return z;}
static D2 inv(const D2&a){D2 z;I iv=I(1)/a.v,iv2=iv*iv,iv3=iv2*iv;z.v=iv;for(int i=0;i<2;i++){z.g[i]=-a.g[i]*iv2;for(int j=0;j<2;j++)z.h[i][j]=I(2)*a.g[i]*a.g[j]*iv3-a.h[i][j]*iv2;}return z;}
static D2 operator/(const D2&a,const D2&b){return a*inv(b);} static D2&operator+=(D2&a,const D2&b){a=a+b;return a;}
static I sym(const char*e){I a(e);return I(rneg(a.hi,MPFR_RNDD),a.hi);}
static D2 common_tail(){D2 z;z.v=sym("1e-75");z.g[0]=sym("1e-114");z.g[1]=sym("1e-100");z.h[0][0]=sym("1e-152");z.h[0][1]=z.h[1][0]=sym("1e-138");z.h[1][1]=sym("1e-124");return z;}
struct EV{D2 v,d1,d2,d3;};
static EV evalp(const std::vector<D2>&c,const D2&z){EV e{D2(0),D2(0),D2(0),D2(0)};for(int n=(int)c.size()-1;n>=0;n--){e.d3=e.d3*z+D2(3)*e.d2;e.d2=e.d2*z+D2(2)*e.d1;e.d1=e.d1*z+e.v;e.v=e.v*z+c[n];}return e;}
static std::array<D2,6> horizon(D2 b,D2 R,int N,const I&t0){std::vector<D2>a(N+1),g(N+1),w(N+1),q(N+1);a[1]=D2(1);g[1]=D2(1);g[2]=D2(1)+b;for(int l=1;l<N;l++){D2 s(0);for(int i=1;i<=l+1;i++)s+=D2(i%2?-1:1)*g[i]*(D2(1)+b*a[l+1-i])*(D2((long)(l+1)*(l+1-i))+D2((long)i*(i+1))/D2(6));a[l+1]=(a[l]*D2((long)(2*l*l+2*l+1))-a[l-1]*D2((long)l*l)-D2(3)*s)/D2((long)(l+1)*(l+1));if(l+2<=N){D2 s2(0);for(int i=0;i<=l;i++)s2+=(a[i]+a[l+1-i]*(D2(1)+b*a[i]))*D2((long)(l+1-i)*(l-3*i));g[l+2]=D2((l+1)%2?-1:1)*s2/(D2("0.5")*D2((long)(l+3)*(l+2)*(l+1)*l));}}
 for(int n=0;n<=N;n++){w[n]=n?D2(1)+b*a[n]:D2(1);if(n==0)q[n]=D2(1);else if(n==1)q[n]=D2(-2)-D2(3)*b;else if(n==2)q[n]=D2(1)+D2(3)*b*g[2];else q[n]=D2(3)*b*D2(n%2?-1:1)*g[n];}
 D2 t(t0),x=D2(1)-D2(1)/t,z=R*x;auto W=evalp(w,z),Q=evalp(q,z);D2 tw=common_tail();W.v+=tw;W.d1+=tw;W.d2+=tw;W.d3+=tw;Q.v+=tw;Q.d1+=tw;Q.d2+=tw;Q.d3+=tw;D2 A=t/W.v,p=D2(1)/W.v-(R/t)*W.d1/(W.v*W.v),C=Q.v*t*t*(t-D2(1))/R-t,d=(D2(3)*t*t-D2(2)*t)*Q.v/R+(t-D2(1))*Q.d1-D2(1),Y=(D2(6)*t-D2(2))*Q.v/R+(D2(4)-D2(2)/t)*Q.d1+R*(t-D2(1))/(t*t)*Q.d2,J=D2(6)*Q.v/R+D2(6)*Q.d1/t+D2(3)*R*Q.d2/(t*t)+R*R*(t-D2(1))*Q.d3/(t*t*t*t);return{A,p,C,d,Y,J};}
int main(){PREC=512;I b0("0.3633018786279300967082111952175988891934060773442491092361561714403219271048605437463631572850729612"),R0("0.6972243957177816041843153424332653277909997011594910385461884399871872297290445630111565633484360785");I rb("5e-45"),rr("5e-31"),sb("1e-41"),sr("1e-27");D2 b(I(rsub(b0.lo,rb.hi,MPFR_RNDD),radd(b0.hi,rb.hi,MPFR_RNDU))),R(I(rsub(R0.lo,rr.hi,MPFR_RNDD),radd(R0.hi,rr.hi,MPFR_RNDU)));b.g[0]=sb;R.g[1]=sr;auto X=horizon(b,R,180,I(20)/I(19));const char*n[6]={"A","p","C","d","Y","J"};for(int i=0;i<6;i++){std::cout<<n[i]<<"\n";for(int a=0;a<2;a++)std::cout<<" g"<<a<<" "<<str(X[i].g[a],45)<<"\n";std::cout<<" h00 "<<str(inflate_abs(X[i].h[0][0],"1e-80"),45)<<"\n";std::cout<<" h01 "<<str(inflate_abs(X[i].h[0][1],"1e-80"),45)<<"\n";std::cout<<" h11 "<<str(inflate_abs(X[i].h[1][1],"1e-80"),45)<<"\n";} }
