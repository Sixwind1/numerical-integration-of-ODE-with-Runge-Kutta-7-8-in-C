#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QRsolve.h"
#include "basicOperation.h"



int newton(int m, int n, double *x, int (*fdf)(int m,int n,double *x,double *f,double *df),int maxit, double tol) // leaves x as a zero of f
{
    /// input: dimension, initial guess x, function to evaluate f and df at point x, maximum iterations and tolerance
    /// output: returns 1 if it fails, solution in x otherwise
    double *f=malloc(m*sizeof(double));
    double *df=malloc(n*m*sizeof(double));
    double *tau=malloc(n*sizeof(double));
    double *y=malloc(n*sizeof(double));

    for(int k=0;k<maxit;k++)
    {
        if(fdf(m,n,x,f,df)==1) // if the fdf function fails (returns 1), stop the program (return 1)
        {
            printf("Error computing the value of f and df\n");
            return 1;
        }

        //printf("Newton: iteration %i\n",k);
        //printf("x_k = ");
        for(int i=0;i<n;i++)
        {
            //printf("%f ",x[i]);
        }
        //printf("\n");

        if (norm2(m,f)<tol)  // if the norm of f(x) is small enough, we say we found the solution
        {
            //printf("Newton: solution found\n");
            free(f);
            free(df);
            free(tau);
            free(y);
            return 0;
        } // otherwise, perform one iteration
        QRsolve(m,n,df,f,tau,y,tol);
        for(int i=0;i<n;i++)
            x[i]=x[i]-y[i];
    }

    printf("Convergence error: maximum allowed iteration reached\n");
    free(f);
    free(df);
    free(tau);
    free(y);
    return 1;
}

