#include "mpfr_interval.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace mpi;

static Real powu(const Real&a,unsigned long n,mpfr_rnd_t rnd){
    Real r(1),b=a;
    while(n){if(n&1)r=rmul(r,b,rnd);n>>=1;if(n)b=rmul(b,b,rnd);}return r;
}

// Majorant from the proved horizon coefficient bounds
// |w_n| <= 10^n/(n+1)^2, |q_n| <= 5 10^n/(n+1)^4.
// At |z|<=0.035 the ratio of successive differentiated terms is <0.36;
// the explicit finite sum plus 5/3 times its final term therefore encloses the tail.
static I tailb(int N,int k,const Real& zmax,bool q){
    Real total(0),last(0); long pref=q?5:1; int power=q?4:2;
    for(int n=N+1;n<=N+300;n++){
        long ff=1; for(int j=0;j<k;j++)ff*=n-j;
        Real term=rdiv(
            rmul(rmul(rmul(Real(ff),powu(Real(10),n,MPFR_RNDU),MPFR_RNDU),
                      powu(zmax,n-k,MPFR_RNDU),MPFR_RNDU),Real(pref),MPFR_RNDU),
            powu(Real(n+1),power,MPFR_RNDD),MPFR_RNDU);
        total=radd(total,term,MPFR_RNDU); last=term;
    }
    total=radd(total,rmul(last,Real("1.666666666666666666666666666667"),MPFR_RNDU),MPFR_RNDU);
    return I(rneg(total,MPFR_RNDD),total);
}

struct Eval { I v,d1; };
static Eval eval(const std::vector<I>&c,const I&z){
    Eval e{I(0),I(0)};
    for(int n=(int)c.size()-1;n>=0;n--){e.d1=e.d1*z+e.v;e.v=e.v*z+c[n];}
    return e;
}

int main(){
    PREC=512;
    const I b0("0.3633018786279300967082111952175988891934060773442491092361561714403219271048605437463631572850729612");
    const I R0("0.6972243957177816041843153424332653277909997011594910385461884399871872297290445630111565633484360785");
    const I rb("5e-45"), rR("5e-31");
    const I b(rsub(b0.lo,rb.hi,MPFR_RNDD),radd(b0.hi,rb.hi,MPFR_RNDU));
    const I R(rsub(R0.lo,rR.hi,MPFR_RNDD),radd(R0.hi,rR.hi,MPFR_RNDU));
    if(!lt(b.hi,Real("0.37"))){std::cerr<<"majorant range failure\n";return 2;}
    const int N=180;
    std::vector<I>a(N+1),g(N+1),w(N+1),q(N+1);
    a[1]=I(1);g[1]=I(1);g[2]=I(1)+b;
    for(int l=1;l<N;l++){
        I s(0);
        for(int i=1;i<=l+1;i++)
            s+=(i%2?I(-1):I(1))*g[i]*(I(1)+b*a[l+1-i])*
               (I((long)(l+1)*(l+1-i))+I((long)i*(i+1))/I(6));
        a[l+1]=(a[l]*I((long)(2*l*l+2*l+1))-a[l-1]*I((long)l*l)-I(3)*s)/I((long)(l+1)*(l+1));
        if(l+2<=N){
            I s2(0);
            for(int i=0;i<=l;i++)
                s2+=(a[i]+a[l+1-i]*(I(1)+b*a[i]))*I((long)(l+1-i)*(l-3*i));
            g[l+2]=((l+1)%2?I(-1):I(1))*s2/
                (I("0.5")*I((long)(l+3)*(l+2)*(l+1)*l));
        }
    }
    for(int n=0;n<=N;n++){
        w[n]=n?I(1)+b*a[n]:I(1);
        if(n==0)q[n]=I(1);
        else if(n==1)q[n]=I(-2)-I(3)*b;
        else if(n==2)q[n]=I(1)+I(3)*b*g[2];
        else q[n]=I(3)*b*(n%2?I(-1):I(1))*g[n];
    }
    Real zmax=rdiv(R.hi,Real(20),MPFR_RNDU);
    if(!lt(zmax,Real("0.035"))){std::cerr<<"z range failure\n";return 3;}
    I tw0=tailb(N,0,zmax,false),tw1=tailb(N,1,zmax,false),tq0=tailb(N,0,zmax,true);
    Real minW("1e100"),minWp("1e100"),minQ("1e100");
    const int pieces=256;
    for(int j=0;j<pieces;j++){
        Real zl=rdiv(rmul(zmax,Real(j),MPFR_RNDD),Real(pieces),MPFR_RNDD);
        Real zh=rdiv(rmul(zmax,Real(j+1),MPFR_RNDU),Real(pieces),MPFR_RNDU);
        I z(zl,zh);
        auto W=eval(w,z),Q=eval(q,z);
        W.v+=tw0;W.d1+=tw1;Q.v+=tq0;
        if(lt(W.v.lo,minW))minW=W.v.lo;
        if(lt(W.d1.lo,minWp))minWp=W.d1.lo;
        if(lt(Q.v.lo,minQ))minQ=Q.v.lo;
        if(!lt(Real(0),W.v.lo)||!lt(Real(0),W.d1.lo)||!lt(Real(0),Q.v.lo)){
            std::cerr<<"positivity failure on piece "<<j<<" z="<<str(z,40)<<"\n";
            std::cerr<<"W="<<str(W.v,40)<<" Wz="<<str(W.d1,40)<<" Q="<<str(Q.v,40)<<"\n";
            return 4;
        }
    }
    std::cout<<"Einstein--Weyl horizon positivity certificate\n";
    std::cout<<"parameter box b "<<str(b,55)<<"\n";
    std::cout<<"parameter box R "<<str(R,55)<<"\n";
    std::cout<<"zmax "<<str(I(zmax),45)<<" pieces "<<pieces<<" N "<<N<<"\n";
    std::cout<<"tail W "<<str(tw0,12)<<"\n";
    std::cout<<"tail Wz "<<str(tw1,12)<<"\n";
    std::cout<<"tail Q "<<str(tq0,12)<<"\n";
    std::cout<<"lower bounds W "<<str(I(minW),35)<<" Wz "<<str(I(minWp),35)<<" Q "<<str(I(minQ),35)<<"\n";
    std::cout<<"PASS: W>0, dW/dz>0, Q>0 from the horizon to t=20/19\n";
}
