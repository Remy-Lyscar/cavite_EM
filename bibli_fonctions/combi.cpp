/* Calcul de c(n,k) en reels, pas d'avertissement en cas de depassement */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
double combi(int n,int k)
{
  int c,j;
  if(k>n-k) k=n-k;
  for(c=1,j=1;j<=k;j++) c=c*(n-j+1)/j;
  return c;
}
