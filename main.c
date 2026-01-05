#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "basicOperation.h"
#include "QRsolve.h"
#include "newton.h"
#include "rk78.h"
#include "flow.h"


int opham_fdf(  // fills f and df
            int m,double hh,double *xx,double *cg,double *f,double *df,
            double h0rk,double hminrk,double hmaxrk,double tolrk,int nMaxRKStep,
            int (*ham)(double *x,double *h,double *prm),
            int (*field)(int n,double t,double *x,double *f,double *prm),
            double *prm)
{
    int n=2*m;
    double T=xx[0]; // copy xx[0] and store it in a local variable to pass it to flow. This way we do not modify xx
    double *x=malloc((n+n*n)*sizeof(double));
    for(int i=0;i<n;i++)
        x[i]=xx[1+i];
    double Hx=0.;

    if(ham(&xx[1],&Hx,prm)!=0)
    {
        free(x);
        return 1;
    }

    if(df==NULL)
    {
        if(flow(&T,x,&h0rk,T,hminrk,hmaxrk,tolrk,nMaxRKStep,n,field,prm)!=0)  /// does the initial time matter? t=T?
        {
            free(x);
            return 1;
        }
        f[0]=Hx-hh;
        f[1]=dotProduct(n,cg,&xx[1])+cg[n];
        for(int i=0;i<n;i++)
            f[2+i]=x[i]-xx[1+i];
        free(x);
    }

    else // df != NULL
    {
        ICVariational(n,x);
        if(flow(&T,x,&h0rk,T,hminrk,hmaxrk,tolrk,nMaxRKStep,n+n*n,field,prm)!=0) /// initial time = T ?
        {
            printf("opham_fdf: could not finish flow\n");
            free(x);
            return 1;
        }
        colToRowReading(n,n,&x[n]);
        f[0]=Hx-hh;
        f[1]=dotProduct(n,cg,&xx[1])+cg[n];
        for(int i=0;i<n;i++)
            f[2+i]=x[i]-xx[1+i];

        double *fxx=malloc(n*sizeof(double));
        field(n,xx[0],&xx[1],fxx,prm); // fxx = f(IC)

        df[0]=0.;
        df[1*(n+1)+0]=0.;
        for(int j=0;j<m;j++) // put DH into df
        {
            df[0*(n+1)+1 +j]=-fxx[m +j];
            df[1*(n+1)-1 -j]=fxx[m-1 -j];
        }
        for(int j=0;j<n;j++) // put Dg into df
            df[1*(n+1)+1 +j]=cg[j];

        initMat(1,n,fxx); // clear fxx and reuse it
        field(n,T,x,fxx,prm); // since x[i], i=0,1,...,n-1 is phi_T(x) --> fxx = f(phi_T(xx))
        for(int i=0;i<n;i++)
        {
            df[(2+i)*(n+1)+0]=fxx[i];  // put f(phi_T(x)) into df
            for(int j=0;j<n;j++)
                df[(2+i)*(n+1)+1 +j]=x[n +i*n+j]; // put Dphi_T(x) == x[i] for i>=n stored by rows
            df[(2+i)*(n+1)+1+ i]-=1.; // subtract I from Dphi_T(x)
        }
        free(x);
        free(fxx);
        return 0;
    }
    free(x);
    return 1;
}



int pendol (int n, double t, double *x, double *f, double *prm)
{
    f[0]=x[1];
    f[1]=-sin(x[0]);
    if (n==2+2*2)
    {
        f[2]=x[3];
        f[3]=-x[2]*cos(x[0]);
        f[4]=x[5];
        f[5]=-x[4]*cos(x[0]);
    }
    return 0;
}

int hamPendol(double *x,double *h,double *prm)
{
    *h=x[1]*x[1]/2.-cos(x[0]);
    return 0;
}



int opham (  // iterate xx until obtaining a solution
        int m, double hh, double *xx, double tol, int maxit, double *cg,
        double h0rk, double hminrk, double hmaxrk, double tolrk, int nMaxRKStep,
        int (*ham)(double *x, double *h, double *prm),
        int (*field)(int n, double t, double *x, double *f, double *prm),
        double *prm)
{
    int n=2*m;
    double *f=malloc((n+2)*sizeof(double));
    double *df=malloc(((n+2)*(n+1))*sizeof(double));
    double *tau=malloc((n+1)*sizeof(double));
    double *y=malloc((n+1)*sizeof(double));

    for(int k=0;k<maxit;k++)
    {
        if(opham_fdf(m,hh,xx,cg,f,df,h0rk,hminrk,hmaxrk,tolrk,nMaxRKStep,ham,field,prm))
        {
            printf("Error computing the value of f and df\n");
            return 1;
        }
        if (norm2(n+2,f)<tol)  // if the norm of f(x) is small enough, x is a zero of f. Return
        {
            printf("root found! f(x) has norm %.16G\n",norm2(n+2,f));
            free(f);
            free(df);
            free(tau);
            free(y);
            return 0;
        }

        printf("opham: iteration %i, the norm is %.16G\n",k,norm2(n+2,f));
        QRsolve(n+2,n+1,df,f,tau,y,tol);  // if the norm is still large, compute the next iteration
        for(int i=0;i<n+1;i++)
            xx[i]=xx[i]-y[i];

    }
    printf("Convergence error: maximum allowed Newton iterations reached\n");
    free(f);
    free(df);
    free(tau);
    free(y);
    return 1;
}


#define RTBP_M 3
#define RTBP_N 6

#define SQR(x) ((x)*(x))

#define N RTBP_N
#define NV1 42

#define X x[0]
#define Y x[1]
#define Z x[2]
#define PX x[3]
#define PY x[4]
#define PZ x[5]
#define MU (*((double *)prm))

int rtbp (int n, double t, double x[], double f[], double *prm)
{
   double r1x, r2x, r1x2, r2x2, r1ypz, r12, r22, r13, r23, p13, p23,
	  p123, p15, p25, p125, p15x, p25x, p125x,
	  dxpx, dypx, dzpx, dypy,
	  dzpy, dzpz;
   int j;
/* Equations */
   f[0]=PX+Y; f[1]=PY-X; f[2]=PZ;
   r1x=X-MU; r2x=r1x+1; r1x2=r1x*r1x; r2x2=r2x*r2x;
   r1ypz=Y*Y+Z*Z;
   r12=r1x2+r1ypz; r22=r2x2+r1ypz;
   r13=r12*sqrt(r12); r23=r22*sqrt(r22);
   p13=(1-MU)/r13; p23=MU/r23;
   f[3]=PY-(p13*r1x+p23*r2x);
   p123=p13+p23;
   f[4]=-PX-Y*p123;
   f[5]=-Z*p123;
   if (n>=NV1) {
   /* First variational equations */
      p15=p13/r12; p25=p23/r22;
      p125=p15+p25;
      p15x=p15*r1x; p25x=p25*r2x;
      p125x=p15x+p25x;
      dxpx=-p123+3*(p15x*r1x+p25x*r2x);
      dzpx=dypx=3*p125x;
      dypx*=Y; dzpx*=Z;
      dzpy=dypy=3*Y*p125;
      dypy=Y*dypy-p123;
      dzpy*=Z;
      dzpz=-p123+3*Z*Z*p125;
  /* Derivative w.r.t. dependent variables of the field for the flow derivative w.r.t. initial conditions (dependent variable part) */
      for (j=0; j<N; j++) {
#define VARF(i,j) f[(1+(j))*N+(i)]
#define VAR(i,j) x[(1+(j))*N+(i)]
	 VARF(0,j)= VAR(1,j)+VAR(3,j);
	 VARF(1,j)=-VAR(0,j)+VAR(4,j);
	 VARF(2,j)= VAR(5,j);
	 VARF(3,j)=dxpx*VAR(0,j)+dypx*VAR(1,j)+dzpx*VAR(2,j)
	           +VAR(4,j);
	 VARF(4,j)=dypx*VAR(0,j)+dypy*VAR(1,j)+dzpy*VAR(2,j)
	           -VAR(3,j);
	 VARF(5,j)=dzpx*VAR(0,j)+dzpy*VAR(1,j)+dzpz*VAR(2,j);
#undef VAR
#undef VARF
      }
   }
   return 0;
}

/*
 * Hamiltonian of the restricted problem "compatible" with rtbp()
 */

int rtbp_h (double x[], double *h, double *prm) {
   double mu=*((double *)prm), xmmu=X-mu, xmmup1=xmmu+1,
	  r12=SQR(xmmu)+SQR(Y)+SQR(Z),
	  r22=SQR(xmmup1)+SQR(Y)+SQR(Z),
	  r1=sqrt(r12), r2=sqrt(r22), /*r13=r1*r12, r23=r2*r22,*/
	  p1=(1-mu)/r1, p2=mu/r2;
   /*if (dmu!=NULL) *dmu=1/r1-1/r2-xmmu*p1/r12-xmmup1*p2/r22;*/
   *h=.5*(SQR(PX)+SQR(PY)+SQR(PZ))+Y*PX-X*PY-p1-p2;
   return 0;
}

#undef MU



int main()
{

    int n=6;
    int m=n/2;
    double xx[n+1];
    double cg[n+1];
    double h0rk,hminrk,hmaxrk,tolrk=1e-13;
    int nMaxRKStep;

    double hh=-1.500384;//hh=-1.500384
    int maxit=100; // Newton
    double tol=1e-10; // Newton
    double prm=3.040357143*1e-6; // mu
    double T=3.051858;//T=3.051858
    double CI[6]={-0.988950, 0., 0.003235, 0., -0.999225, 0.};//CI[6]={-0.988950, 0., 0.003235, 0., -0.999225, 0.}

    hminrk = 1e-6;
    hmaxrk = 0.5;
    h0rk=0.1;
    nMaxRKStep=1000000; // 1e6 (enough to go with hmin)

    xx[0]=T;
    for(int i=0;i<n;i++)
    {
        xx[1+i]=CI[i];
        cg[i]=0.;
    }
    cg[1]=1.;
    cg[n]=0.;



    if(opham(m,hh,xx,tol,maxit,cg,h0rk,hminrk,hmaxrk,tolrk,nMaxRKStep,rtbp_h,rtbp,&prm)!=0)
    {
        printf("main: opham could not finish!\n");
        return 1;
    }
    printf("\n\n");
    printf("main: finished computing the root and obtained this point as solution\n");
    printMat(1,n+1,xx);

    double solution[7]={
        3.05177804298575,
        -0.9889514941464619,
        0.,
        0.003249420704787417,
        0.,
        -0.9992369544112847,
        0.
    };
    double err[7];
    for(int i=0;i<n+1;i++)
        err[i]=fabs(xx[i]-solution[i]);
    printf("main: finished computing the root and obtained the following\n");
    printf("\n\n");
    printf("pdf solution                 obtained solution           absolute error\n");
    for(int i=1;i<n+1;i++)
    {
        if(solution[i]==0.)
            printf("%.16G                             %.16G       %.16G\n",solution[i],xx[i],err[i]);
        else
            printf("%.16G          %.16G          %.16G\n",solution[i],xx[i],err[i]);
    }
    printf("\n\n");
    printf("T pdf=%.16G     T obtained=%.16G       absolute error=%.16G\n",solution[0],xx[0],err[0]);



    FILE *dades=fopen("ophalo.txt","w");
    double t=xx[0];
    double x[n];
    memcpy(x,&xx[1],n*sizeof(double));
    h0rk=0.1;
    if(flowTXT(&t,x,&h0rk,xx[0],hminrk,hmaxrk,tolrk,nMaxRKStep,n,rtbp,&prm,dades))
    {
        printf("main: error in flowTXT, could not finish");
    }
    return 0;
}


