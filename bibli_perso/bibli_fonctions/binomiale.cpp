/* Fournit la valeur de la loi binomiale pour une proba p, un nombre d'epreuves n et un nombre de realisations k */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
double binomiale(double p,int n,int k)
{
  int i; 
  double b,q=1.-p;
  b=pow(q,n);
  for(i=1;i<=k;i++) b=b*(n-i+1)/i*p/q; 
  return b;
}
