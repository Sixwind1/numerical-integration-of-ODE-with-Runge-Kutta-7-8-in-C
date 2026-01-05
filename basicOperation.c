#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void printMat(int n,int m,double A[m*n])
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            printf("%.16G ",A[i*n+j]);
        }
        printf("\n");
    }
}

double dotProduct(int dim,double v[dim],double w[dim])
{
    double resultat=0.;
    for(int i=0;i<dim;i++)
        resultat+=v[i]*w[i];
    return resultat;
}

double norm2(int dim,double v[dim])
{
    double resultat=0.;
    for(int i=0;i<dim;i++)
    {
        resultat+=v[i]*v[i];
    }
    return sqrt(resultat);
}

void initMat(int m,int n, double *A)
{
    for(int i=0;i<m;i++)
    {
         for(int j=0;j<n;j++)
             A[i*n+j]=0.;
    }
}

void colToRowReading(int n,int m, double *a)
{
    double b[n*m];
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
            b[i*n+j]=a[j*m+i];
    }
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
            a[i*n+j]=b[i*n+j];
    }

}


