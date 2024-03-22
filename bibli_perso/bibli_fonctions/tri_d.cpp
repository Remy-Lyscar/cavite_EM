#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
void tri_d(double n[],int nmax)
{
  double nn;
  int i,j;
  for(j=1; j<nmax; j++)
    {
      nn=n[j];
      for(i=j-1; i>=0; i--)
	{
	  if(n[i]<=nn) goto vega;
	  n[i+1]=n[i];
	}
      i=-1;
      vega : n[i+1]=nn;
    }
}
  
