/* Calcul d'une integrale par la méthode des rectangles au milieu.
   Il faut fournir les bornes inf et sup, la fonction, le nombre de rectangles */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
double integ(double a,double b,double (*f)(double),int n)
{
  double h,s; int i;
  h=(b-a)/n;
  for(s=0.,i=1; i<=n; i++) s=s+(*f)(a+(i-0.5)*h);
  return (h*s);
}




