/* Calcul d'une integrale par la méthode des rectangles au milieu.
   Il faut fournir les bornes inf et sup, la fonction, le nombre de rectangles 
   et une tolerance eps. 
   Le nombre de rectangles fourni n'est qu'un point de depart.
   Il est augmente jusqu'à ce que la variation relative du resultat soit inferieure 
   en valeur absolue a eps.  
   On fait appel a la fonction integ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
double integ_preci(double p,double q,double (*f)(double),int k,double eps)
{
  double i,ip,d;
  ip=integ(p,q,(*f),k);
  do
    {
      k=2*k; 
      i=integ(p,q,(*f),k);
      d=fabs((i-ip)/ip);
      ip=i;
    }
  while(d>eps);
  printf("Nombre final de rectangles :%d\n",k);
  return (i);
}



