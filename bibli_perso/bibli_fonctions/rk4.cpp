#include "bibli_fonctions.h"
void rk4(void(*sd)(double *,double,double *,int),double *q,double t,double dt,int n)
{
  int i,k,p,PM=4;
  static const double c2=1./2,c3=1./3,c6=1./6;
  /* Allocations et initialisations */
  double **a=(double **)malloc(PM*sizeof(double *)); 
  for(i=0;i<PM;i++) a[i]=(double *)malloc(PM*sizeof(double));
  ini_D_2(a,PM,PM,c2,0.,0.,0.,0.,c2,0.,0.,0.,0.,1.,0.,c6,c3,c3,c6);
  double *b=(double *)malloc(PM*sizeof(double));
  ini_D_1(b,PM,0.,c2,c2,1.);
  double **y=(double **)malloc((PM+1)*sizeof(double *)); 
  for(i=0; i<PM+1; i++) y[i]=(double *)malloc(n*sizeof(double));
  double **z=(double **)malloc(PM*sizeof(double *)); 
  for(i=0; i<PM; i++) z[i]=(double *)malloc(n*sizeof(double));
  /* Calcul */
  for(i=0; i<n; i++) y[0][i]=q[i];
  for(p=1; p<=PM; p++)
    {
      sd(y[p-1],t+b[p-1]*dt,z[p-1],n);
      for(i=0; i<n; i++) y[p][i]=q[i];
      for(k=0; k<p; k++) {for(i=0; i<n; i++) y[p][i]=y[p][i]+dt*a[p-1][k]*z[k][i];}
    }   
  for(i=0; i<n; i++) q[i]=y[PM][i];
  /* Desallocations */
  f_D_2(a,PM,PM); f_D_1(b,PM); f_D_2(y,PM+1,n); f_D_2(z,PM,n); 
}
