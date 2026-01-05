
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define sgn(x) ((x)>=0? 1.:-1.)


int rk78(double *t, double x[], double *h,
         double hmin, double hmax, double tol,
         int n, int (*field)(int n, double t, double x[], double f[], double *prm),
         double *prm)
{
    static double c[13] = {
        0., 2. / 27., 1. / 9., 1. / 6., 5. / 12.,
        0.5, 5. / 6., 1. / 6., 2. / 3., 1. / 3.,
        1., 0., 1.
    };

    static double A[78] = {
        2. / 27., 1. / 36., 1. / 12., 1. / 24.,
        0., 1. / 8., 5. / 12., 0., -25. / 16.,
        25. / 16., 0.05, 0., 0., 0.25,
        0.2, -25. / 108., 0., 0., 125. / 108.,
        -65. / 27., 125. / 54., 31. / 300., 0., 0.,
        0., 61. / 225., -2. / 9., 13. / 900., 2.,
        0., 0., -53. / 6., 704. / 45., -107. / 9.,
        67. / 90., 3., -91. / 108., 0., 0.,
        23. / 108., -976. / 135., 311. / 54., -19. / 60., 17. / 6.,
        -1. / 12., 2383. / 4100., 0., 0., -341. / 164.,
        4496. / 1025., -301. / 82., 2133. / 4100., 45. / 82., 45. / 164.,
        18. / 41., 3. / 205., 0., 0., 0.,
        0., -6. / 41., -3. / 205., -3. / 41., 3. / 41.,
        6. / 41., 0., -1777. / 4100., 0., 0.,
        -341. / 164., 4496. / 1025., -289. / 82., 2193. / 4100., 51. / 82.,
        33. / 164., 12. / 41., 0., 1.
    };

    static double b[11] = {
        41. / 840., 0., 0., 0., 0.,
        34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280.,
        41. / 840.
    };

    static double bp[13] = {
        0., 0., 0., 0., 0.,
        34. / 105., 9. / 35., 9. / 35., 9. / 280., 9. / 280.,
        0., 41. / 840., 41. / 840.
    };

    int indexA;
    double  coef, err;
    double k[13*n];
    double rk7[n];
    double rk8[n];

    do
    {
        // computation of the k's
        indexA = 0;
        for(int i=0;i<13;i++) // this routine computes the k_i's
        {
            memcpy(rk7,x,n*sizeof(double)); // initialize rk7 = x (would it be better to initialize with zeros and then add x?)

            for(int j=0;j<i;j++) // rk7 +=  sum(h*aij*K_j)
            {
                if(A[indexA]!=0.) // A[indexA]==0 --> coef==0 --> rk7 is not modified. We avoid this computation to reduce error
                {
                    coef=A[indexA]* (*h); // coef = h*aij
                    for (int l=0;l<n;l++) // rk7 = x+(h*aij*kj)
                        rk7[l]+=coef*k[j*n+l];
                }
                indexA++;
            }

            // k_i = f(t_field,rk7)

            if (i==0 || i==11) // c[i]==0 for these coefficients
            {
                if(field(n,(*t),rk7,&k[i*n],prm))
                {
                    printf("RK78: could not evaluate the field");
                    return 1;
                }
            }
            else
            {
                if(field(n,(*t)+c[i]* (*h),rk7,&k[i*n],prm))
                {
                    printf("RK78: could not evaluate the field");
                    return 1;
                }
            }
        }

        // computation of rk7, rk8 and error in L1 norm (sensitive to nearly zero values)
        err=0;
        for(int l=0;l<n;l++)
        {
            rk7[l]=x[l];  // rk7 = rk8 = x
            rk8[l]=x[l];

            for(int i=5;i<11 && i!=10;i++) // we ALMOST compute rk7=x+sum(h*k_i*b_i) and rk8=x+sum(h*k_i*b'_i); we skip some indices that are zero or convenient to skip
            {
                coef= (*h)*k[i*n+l];  // reuse memory: coef = h*k_i
                rk7[l]+=coef*b[i];  // rk7 = x + h*k_i*b_i
                rk8[l]+=coef*bp[i]; // rk8 = x + h*k_i*b'_i
            }
            rk8[l]+= (*h)*(bp[11]*k[11*n+l] + bp[12]*k[12*n+l]); // complete rk8: missing stages 12 and 13
            rk7[l]+= (*h)*(  b[0]*k[0*n+l]  +   b[10]*k[10*n+l]); // complete rk7: missing stages 1 and 10

            err+=fabs(rk8[l]-rk7[l]); // L1 norm of rk8 - rk7
        }
        // check if we can exit the while loop
        if(err<=tol || fabs(*h)<=hmin) // note it must be <=; we always try with hmin as long as the step formula gives a very small step
            break;

        // otherwise, both are too large, update the step and repeat the routine
        (*h) *= .9*pow(tol/err,.125);

        if(fabs(*h)<hmin) // if the step is too small, take the minimum accepted step
            (*h)=sgn(*h)*hmin;

    }while(1);

    // there are only 2 cases when exiting the loop:
    // err<=tol and |h|>=hmin   or   err>tol and |h|==hmin (by construction)
    // in both cases we take a step (i.e., we accept the error if err>tol, only printing a warning)

    (*t)+=(*h);
    memcpy(x,rk8,n*sizeof(double));

    // update the step
    /// a peculiar phenomenon is that if you do NOT give a sufficiently small minimum error,
    /// the algorithm fails and does not converge.
    if(err<tol) // small error --> large step. We do not want an excessively large step --> avoid very small error (lower bound for error)
        err=((err)<(tol/10) ? tol/10:err);
    else
    {
        printf("\nrk78: Warning! Cannot compute with the requested precision (tol)\n");
        printf("       because it requires the next h to be %.16G and we are at t=%.16G\n",*h,*t);
        printf("Nevertheless, we ignore this error and continue integrating\n");
    }

    (*h) *= .9*pow(tol/err,.125);

    if(fabs(*h)<hmin)
    {
        (*h)=sgn(*h)*hmin;
    }
    else if(fabs(*h)>hmax)
    {
        (*h)=sgn(*h)*hmax;
    }

    return 0;
}

