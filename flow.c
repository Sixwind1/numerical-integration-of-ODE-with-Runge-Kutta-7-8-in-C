#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "basicOperation.h"
#include "rk78.h"


int flow (double *t,double *x,double *h,double T,double hmin,double hmax,double tol,int nMaxStep,int n,int (*field)(int m,double t,double *x,double *f,double *prm),double *prm)
{
    double finalTime=*t+T;
    if(T<0)  // assume the input h is positive
        *h=-(*h); // if we have to integrate backwards, then set h negative
    int numit=0;
    //printf("flow: initial data. ");
    //printf("t= %.8G     x = %.8G        %.8G        %.8G        %.8G        %.8G        %.8G\n",*t,x[0],x[1],x[2],x[3],x[4],x[5]);

    while(*t!=finalTime && numit<nMaxStep)
    {
        if (fabs(*t - finalTime)<=hmax)
            *h=finalTime-(*t);

        if (rk78(t, x, h, hmin, hmax, tol, n, field, prm) != 0)
            {
                printf("flow: RK78 could not finish\n");
                return 1;
            }
        //printf("%.8G     %.8G        %.8G        %.8G        %.8G        %.8G        %.8G\n",*t,x[0],x[1],x[2],x[3],x[4],x[5]);
        numit++;
    }

    printf("flow: the integration has finished and it is at time = %.8G\n",*t);
    if(numit>=nMaxStep)
    {
        printf("flow: the maximum number of allowed RK78 calls has been exhausted, terminating the flow function\n");
        return 1;
    }
    return 0;
}

void ICVariational(int n, double *x)
{
    for(int i=n;i<n+n*n;i++)
        x[i]=0;
    for(int i=0;i<n;i++)
        x[n+i*n+i]=1;
}


int flowTXT (double *t,double *x,double *h,double T,double hmin,double hmax,double tol,int nMaxStep,int n,int (*field)(int m,double t,double *x,double *f,double *prm),double *prm,FILE *data)
{
    double finalTime=*t+T;
    if(T<0)
        *h=-(*h);
    int numit=0;
    //printf("%.8G     %.8G        %.8G        %.8G        %.8G        %.8G        %.8G\n",*t,x[0],x[1],x[2],x[3],x[4],x[5]);
    fprintf(data,"%.8G     %.8G        %.8G        %.8G        %.8G        %.8G        %.8G\n",*t,x[0],x[1],x[2],x[3],x[4],x[5]);

    while(*t!=finalTime && numit<nMaxStep)
    {
        if (fabs(*t - finalTime)<=hmax)
            *h=finalTime-(*t);

        if (rk78(t, x, h, hmin, hmax, tol, n, field, prm) != 0)
            {
                printf("flow: RK78 could not finish\n");
                return 1;
            }
        fprintf(data,"%.8G     %.8G        %.8G        %.8G        %.8G        %.8G        %.8G\n",*t,x[0],x[1],x[2],x[3],x[4],x[5]);
        //printf("%.8G     %.8G        %.8G        %.8G        %.8G        %.8G        %.8G\n",*t,x[0],x[1],x[2],x[3],x[4],x[5]);
        numit++;
    }

    printf("flow: the integration has finished and it is at time = %.8G\n",*t);
    if(numit>=nMaxStep)
    {
        printf("flow: the maximum number of allowed RK78 calls has been exhausted, terminating the flow function\n");
        return 1;
    }
    return 0;
}



