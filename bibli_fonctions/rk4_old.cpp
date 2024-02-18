#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
#define PM 4
void rk4_old(void(*sd)(double *,double,double *,int),double q[],double t,double dt,int n)
{
  /* Declarations et initialisations */
  int i,k,p;
  double *y[PM+1],*z[PM];
  static const double a[PM][PM]={{1./2,0.,0.,0.},{0.,1./2,0.,0.},{0.,0.,1.,0.},{1./6,1./3,1./3,1./6}};
  static const double b[PM]={0.,1./2,1./2,1.};
  /* Allocations */
  for(i=0; i<PM+1; i++) y[i]=(double *)malloc(n*sizeof(double));
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
  for(i=0; i<PM+1; i++) {free(y[i]); y[i]=NULL;}
  for(i=0; i<PM; i++) {free(z[i]); z[i]=NULL;}
}
