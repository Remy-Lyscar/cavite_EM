#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
#define D_t double
#define Int_t int
/* ***************************************************************************** */
/* Resolution d'une equation du second degre, la valeur de la fonction est : */
/*   3 si a=b=c=0 (tout x est racine) */
/*   0 si a=b=0 et c!=0 ou si a!=0 et delta<0 (pas de racines) */
/*   1 si a=0 et b!=0 ou si a!=0 et delta=0 (une racine unique renvoyee dans *xp et *xs) */
/*   2 si a!=0 et delta>0 (deux racines renvoyees dans xp et xs) */
Int_t r2(D_t a, D_t b, D_t c, D_t *xp, D_t *xs)
{
  D_t delta,rd ;
  /* fprintf(fich,"a=%lg b=%lg c=%lg\n",a,b,c); */
  if(a==0.) 
    {
      if(b==0.) 
	{
	  if(c==0.) return 3 ;
	  else return 0 ;
	} 
      else 
	{
	  *xp=-c/b ; *xs=*xp ; return 1 ;
	}
    } 
  else
    {
      delta=b*b-4.*a*c ;
      /* fprintf(fich,"delta=%lg\n",delta); */
      if(delta<0.) return 0 ;
      if(delta==0.) {*xp=-b/2./a ; *xs=*xp ; return 1 ;}
      else {rd=sqrt(delta) ; *xp=(-b-rd)/2./a ; *xs=(-b+rd)/2./a ; return 2 ;}
    } 
}
/* ***************************************************************************** */

