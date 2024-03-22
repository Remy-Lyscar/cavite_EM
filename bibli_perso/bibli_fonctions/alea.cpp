#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
static const int q=1664713,p=1288;
static int n=1; 
void germe(int g)
{
  if((n=g%q)==0) {printf("Prendre une valeur du germe >=1 et <=%d\n",q-1); exit (0);} /* on met g%q et non g pour ne pas risquer de depasser */ 
}                                                                                 /* la valeur max permise pour un int : 2^31-1 */
double alea()
{
  n=n*p%q;
  return (double)(n-1)/(q-2);
}
