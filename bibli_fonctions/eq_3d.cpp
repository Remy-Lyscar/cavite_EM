#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
static double rac(double x)
{
  double y,tt=1./3;
  y=pow(fabs(x),tt);
  if(x<0.) y=-y;
  return y;
}
int eq_3d(double a,double b,double c,double d,double *x)
{
  double r,s,t,p,q,u,v,phi,rr,rrr,rs3,dps3;
  dps3=2*M_PI/3;
  if(a==0.){printf("Equation de degre < 3"); return -1;}
  r=b/a; s=c/a; t=d/a; p=s-r*r/3; q=2*r*r*r/27-r*s/3+t; rs3=r/3;
  u=q*q/4+p*p*p/27;
  if(u>0) {x[0]=rac(-q/2+sqrt(u))+rac(-q/2-sqrt(u))-rs3; return 1;}
  if(u==0.) {v=2*rac(-q/2); x[0]=v-rs3; x[1]=x[2]=-v/2-rs3; return 2;}
  rr=sqrt(-p*p*p/27); rrr=2*sqrt(-p/3); phi=acos(-q/2/rr);
  x[0]=rrr*cos(phi/3)-rs3; x[1]=rrr*cos(phi/3+dps3)-rs3; x[2]=rrr*cos(phi/3-dps3)-rs3; 
  return 3;
}

