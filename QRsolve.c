#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int QRsolve(int m,int n, double* a, double* b, double* tau,double* x,double tol) // matrix a is stored as a vector, a[i,j]=a[i*n+j]
{
    /// input: dimensions, matrix A, vector b, allowed tolerance, memory to store tau's and the solution x
    /// output: matrix A becomes upper triangular, with u_k stored in the lower part. b becomes the associated vector of the upper triangular A.
    ///         tau is filled with intermediate data, x is filled with the solution of the system

    double sk=0; double beta=0; // memory to store variables inside Householder iterations

    // Householder: obtain the upper triangular matrix, with u_k in the lower part, and all tau's filled
    for (int k=0;k<n;k++)
    {
        sk=0;
        for (int i=k;i<m;i++)
            sk+=a[i*n+k]*a[i*n+k];
        if(sk<tol)  // if the squared norm of the vector is very small, stop the algorithm due to ill-conditioning
        {
            printf("in QRsolve: ERROR, nearly singular matrix \n");
            return -1;
        }
        sk=sqrt(sk);
        if(a[k*n+k]>0)
            sk=-sk;
        tau[k]=a[k*n+k]-sk;  // not the final value; we are only reusing memory space to store results
        for(int i=k+1;i<m;i++) // normalize u_k; start from k+1 because the first element (k-th position) is always 1, saving memory
            a[i*n+k]=a[i*n+k]/tau[k];
        tau[k]=-tau[k]/sk;
        a[k*n+k]=sk;

        for(int j=k+1;j<n;j++)
        {
            beta=a[k*n+j];
            for(int i=k+1;i<m;i++)
                beta+=a[i*n+k]*a[i*n+j];
            beta=beta*tau[k];
            a[k*n+j]-=beta;
            for(int i=k+1;i<m;i++)
                a[i*n+j]-=beta*a[i*n+k];

        }
    }

    // compute transpose(Q)*b
    for (int k=0;k<n;k++)
    {
        beta=b[k];
        for(int i=k+1;i<m;i++)
            beta+=a[i*n+k]*b[i];
        beta=beta*tau[k];
        b[k]-=beta;
        for(int i=k+1;i<m;i++)
            b[i]-=beta*a[i*n+k];
    }

    // solve the upper triangular system by substitution
    for(int i=n-1;i>=0;i--)
    {
        x[i]=b[i];
        for (int j=i+1;j<n;j++)
            x[i]-=a[i*n+j]*x[j];

        x[i]=x[i]/a[i*n+i];
    }

    return 0;
}
