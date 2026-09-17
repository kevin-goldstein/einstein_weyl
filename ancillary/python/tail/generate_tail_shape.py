# Generate explicit decimal coefficients for corrected one-Yukawa centered tail shape.
import mpmath as mp
mp.mp.dps=220
T=mp.mpf(40);N=40
A=mp.mpf('1.0324864466507795973677035145514604182115606537271511328159699001039011153579008194924211606515855359414467422')
beta=mp.mpf('-0.56897241457440662855219799450966140112676548116077296026491028896680447198418711399863269966073889639436757618')
a=1/A;nu=beta/(2*A)
# correct recurrence

def M(k):
 z=mp.matrix(4,4)
 if k==0:z[1,2]=1;z[2,3]=1;z[3,0]=-12/A**3;z[3,2]=a**2
 else:
  if k==1:z[0,0]=-2
  else:z[0,0]=beta*(-beta)**(k-2)
  if k>=2:z[0,1]=A/3*(-beta)**(k-2)
  z[0,2]=A/6*(-beta)**(k-1);z[3,0]=6*beta/A**3*(-beta)**(k-1);z[3,1]=2/A**2*(-beta)**(k-1);z[3,2]=a**2*(-beta)**k
 return z
Ms=[M(k) for k in range(N+2)];L=-a*mp.eye(4)-Ms[0];y0=mp.matrix([0,1/(2*a*a),-1/(2*a),mp.mpf('.5')])
vals,ER=mp.eig(L.T);ii=min(range(4),key=lambda i:abs(vals[i]));ell=ER[:,ii];ell=ell/mp.fdot(ell,y0)
Aug=mp.matrix(5,5)
for i in range(4):
 for j in range(4):Aug[i,j]=L[i,j]
 Aug[i,4]=y0[i]
for j in range(4):Aug[4,j]=ell[j]
def solve(rhs):
 rr=mp.matrix(5,1)
 for i in range(4):rr[i]=rhs[i]
 return mp.lu_solve(Aug,rr)[:4,:]
yp=[y0,solve(Ms[1]*y0-nu*y0)]
for n in range(2,N+2):
 rb=mp.matrix(4,1)
 for k in range(1,n+1):rb+=Ms[k]*yp[n-k]
 rb-=(nu-n+1)*yp[n-1]
 D=Ms[1]-(nu-n+1)*mp.eye(4);cc=-mp.fdot(ell,rb)/mp.fdot(ell,D*y0);yp[n-1]+=cc*y0
 rr=mp.matrix(4,1)
 for k in range(1,n+1):rr+=Ms[k]*yp[n-k]
 rr-=(nu-n+1)*yp[n-1];yp.append(solve(rr))
ys=yp[:N+1];PA=[];PC=[]
for n in range(N+1):
 if n==0:PA.append(-ys[n][0]/a);PC.append(-ys[n][1]/a)
 else:PA.append(((nu-n+1)*PA[n-1]-ys[n][0])/a);PC.append(((nu-n+1)*PC[n-1]-ys[n][1])/a)
# derivative operator and normalize v(T)=1
def D(P):
 o=[mp.mpf(0)]*(len(P)+1)
 for j,x in enumerate(P):
  o[j]-=a*x;o[j+1]+=(nu-j)*x
 return o
def peval(P,r):
 s=mp.mpf(0)
 for x in reversed(P):s=s*r+x
 return s
PY=D(PC);PJ=D(PY);vT=peval(PJ,1/T)-a*peval(PY,1/T)
PA=[x/vT for x in PA];PC=[x/vT for x in PC]
print('A0 =',mp.nstr(A,180));print('beta0 =',mp.nstr(beta,180));print('a0 =',mp.nstr(a,180));print('nu0 =',mp.nstr(nu,180));print('vT_raw =',mp.nstr(vT,180))
print('PA = [')
for x in PA:print("    '"+mp.nstr(x,180)+"',")
print(']')
print('PC = [')
for x in PC:print("    '"+mp.nstr(x,180)+"',")
print(']')
