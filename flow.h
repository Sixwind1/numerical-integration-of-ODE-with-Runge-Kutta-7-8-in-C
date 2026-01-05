


// remember that aij=x[n+j*n+i] where i,j are between 0 and n-1

int flow (double *t,double *x,double *h,double T,double hmin,double hmax,double tol,int npasmax,int n,int (*camp)(int m,double t,double *x,double *f,double *prm),double *prm);

void ICVariational(int n, double *x);

int flowTXT (double *t,double *x,double *h,double T,double hmin,double hmax,double tol,int npasmax,int n,int (*camp)(int m,double t,double *x,double *f,double *prm),double *prm,FILE *dades);
