#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <bibli_fonctions.h>
double gauss(double t)
{
  double c=1./sqrt(2.*acos(-1.));
  return (c*exp(-0.5*t*t));
}
int main(void)
{
  double x,ee=1.e-10; int np=100;
  printf("x=? ");scanf("%lg",&x);
  printf("Integrale de -infini a %lg=%15.12lg\n",x,0.5+integ(0.,x,gauss,np)); /* juste pour comparer */ 
  printf("Integrale de -infini a %lg=%15.12lg\n",x,0.5+integ_preci(0.,x,gauss,np,ee)); 
  return (0);
}
