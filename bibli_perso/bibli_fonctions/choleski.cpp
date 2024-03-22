#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
int choleski(double *a,double *b,double *c,double *x,double *y,int n)
{ 
  double *al,*be,nu; int i;
  al=(double *)malloc(n*sizeof(double)); be=(double *)malloc(n*sizeof(double));
  if(b[0]==0.){printf("Echec 1 de choleski"); return 1;} 
  al[0]=-c[0]/b[0]; be[0]=y[0]/b[0]; 
  for(i=1; i<n; i++)
    {
      nu=a[i]*al[i-1]+b[i]; if(nu==0.){printf("Echec 2 de choleski"); return 2;}
      al[i]=-c[i]/nu; be[i]=(y[i]-a[i]*be[i-1])/nu;
    }
  x[n-1]=be[n-1];
  for(i=n-2; i>=0; i--) x[i]=al[i]*x[i+1]+be[i];
  return 0;
}
/* #define K 4
int main(void)
{
  int i;
  double p[K]={0.,1.,6.,2.};
  double q[K]={2.,1.,-2.,-3.};
  double r[K]={1.,-3.,1.,0.};
  double v[K]={7.,-10.,7.,13.};
  double u[K];
  choleski(p,q,r,u,v,K);
  for(i=0; i<K; i++) printf("%lg ",u[i]);
  printf("\n");
} */
