import mpmath as mp
mp.mp.dps=100
A=mp.mpf('1.03248644665077959736772032328869905727154373')
beta=mp.mpf('-0.568972414574406628550212757062083159918152678')
a=1/A;nu=beta/(2*A);N=100

def Mk(k):
 M=mp.matrix(4,4)
 if k==0:
  M[1,2]=1;M[2,3]=1;M[3,0]=-12/A**3;M[3,2]=a**2
 else:
  if k==1:M[0,0]=-2
  else:M[0,0]=beta*(-beta)**(k-2)
  if k>=2:M[0,1]=A/3*(-beta)**(k-2)
  M[0,2]=A/6*(-beta)**(k-1)
  M[3,0]=6*beta/A**3*(-beta)**(k-1)
  M[3,1]=2/A**2*(-beta)**(k-1)
  M[3,2]=a**2*(-beta)**k
 return M
M=[Mk(k) for k in range(N+2)]
L=-a*mp.eye(4)-M[0]
y0=mp.matrix([0,1/(2*a*a),-1/(2*a),mp.mpf('.5')])
# left null
vals,ER=mp.eig(L.T);ii=min(range(4),key=lambda i:abs(vals[i]));ell=ER[:,ii];ell=ell/mp.fdot(ell,y0)
# augmented solve L y=rhs, ell.y=0. last multiplier should zero when solvable
Aug=mp.matrix(5,5)
for i in range(4):
 for j in range(4):Aug[i,j]=L[i,j]
 Aug[i,4]=y0[i]
for j in range(4):Aug[4,j]=ell[j]

def solve_gauge(rhs):
 rr=mp.matrix(5,1)
 for i in range(4):rr[i]=rhs[i]
 rr[4]=0
 sol=mp.lu_solve(Aug,rr)
 return sol[:4,:],sol[4]
# order 1 base
rhs=M[1]*y0-nu*y0
print('solv1',mp.nstr(mp.fdot(ell,rhs),30))
yp1,lam=solve_gauge(rhs)
yp=[y0,yp1]
cs=[mp.mpf(0),None]
# for n=2..N+1 determine c_{n-1}, finalize y_{n-1}, then solve y_n particular
for n in range(2,N+2):
 # base rhs using yp[n-1] plus all finalized earlier
 rbase=mp.matrix(4,1)
 for k in range(1,n+1):
  idx=n-k
  if idx==n-1: vec=yp[idx]
  else: vec=yp[idx] # finalized
  rbase += M[k]*vec
 rbase -= (nu-n+1)*yp[n-1]
 D=M[1]-(nu-n+1)*mp.eye(4)
 denom=mp.fdot(ell,D*y0)
 cprev=-mp.fdot(ell,rbase)/denom
 cs[n-1]=cprev
 yp[n-1]=yp[n-1]+cprev*y0
 # recompute rhs with finalized vectors
 rhs=mp.matrix(4,1)
 for k in range(1,n+1):rhs+=M[k]*yp[n-k]
 rhs-=(nu-n+1)*yp[n-1]
 solv=mp.fdot(ell,rhs)
 if n<8:print('n',n,'cprev',mp.nstr(cprev,25),'solv',mp.nstr(solv,15))
 yn,lam=solve_gauge(rhs)
 yp.append(yn);cs.append(None)
# yp 0..N finalized (N finalized at n=N+1)
ys=yp[:N+1]
alpha=[];cc=[]
for n in range(N+1):
 if n==0:alpha.append(-ys[n][0]/a);cc.append(-ys[n][1]/a)
 else:
  alpha.append(((nu-n+1)*alpha[n-1]-ys[n][0])/a)
  cc.append(((nu-n+1)*cc[n-1]-ys[n][1])/a)

def ratios(t,NN):
 P=[mp.mpf(0)]*6;r=1/mp.mpf(t)
 for n in range(NN,-1,-1):
  P[0]=P[0]*r+alpha[n];P[1]=P[1]*r+ys[n][0];P[2]=P[2]*r+cc[n];P[3]=P[3]*r+ys[n][1];P[4]=P[4]*r+ys[n][2];P[5]=P[5]*r+ys[n][3]
 v=P[5]-a*P[4]
 return [x/v for x in P]
print('a',mp.nstr(a,40),'nu',mp.nstr(nu,40))
for t in [500,300,200,100,60,40]:
 for NN in [4,8,12,20,30,40,50,60,70,80,90,100]:
  r=ratios(t,NN)
  if NN in [20,40,60,80,100]:print('t',t,'N',NN,*[mp.nstr(x,18) for x in r])
 print()
# save coefficients as repr maybe
