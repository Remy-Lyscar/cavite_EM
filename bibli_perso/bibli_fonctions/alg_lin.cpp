#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bibli_fonctions.h"
/* Mise a zero d'un vecteur a n elements */
void zero_v(double *a,int n)
{
  int i;
  for(i=0; i<n; i++) {a[i]=0.;}
}
/* ----------------------------------------------------------- */
/* Mise a zero d'une matrice m,n */
void zero(double **a,int m,int n)
{
  int i,j;
  for(i=0; i<m; i++) {for(j=0; j<n; j++) a[i][j]=0.;}
}
/* ----------------------------------------------------------- */
/* Fonction imprimant une matrice quelconque */ 
void mat_imp(double **x,int l,int c)
{
  int i,j;
  for(i=0; i<l; i++) 
    {
      for(j=0; j<c; j++) printf("%e ",x[i][j]); 
      printf("\n");
    } 
  printf("\n");
}
/* ********************************************************************** */
void mat_prod(double **a,double **b,double **ab,int m,int n,int p)
     /* Fait le produit ab de la matrice a(m,n) par la matrice b(n,p) */
     /* Le resultat est dans ab(m,p) */
{
  int i,j,k;
  for(i=0; i<m; i++) 
    {
      for(k=0; k<p; k++) 
	{
	  ab[i][k]=0.;
	  for(j=0; j<n; j++) ab[i][k]+=a[i][j]*b[j][k];
	}
    }
}
/* ********************************************************************** */
int *ivector(int nl,int nh)
     /* Allocate an int vector with subscript range v[nl...nh] */ 
{
  int *v;
  v=(int *)malloc((size_t)((nh-nl+2)*sizeof(int)));
  if(!v) printf("Allocation failure in ivector()");
  return(v-nl+1);
}
/* ********************************************************************** */
void free_ivector(int *v,int nl,int nh)
     /* Free an int vector allocated with ivector() */ 
{
  free((char *)(v+nl-1));
}
/* ********************************************************************** */
void swap(double *x,double *y)
     /* Echange de deux quantites */
{
  double z;
  z=*x; *x=*y; *y=z;
}
/* *********************************************************************** */
void g_j(double **a,int n,double **b,int m)
     /* Par rapport a Numerical Recipes j'ai diminue de 1 tous les indices de a et b pour passer des matrices */
     /* en argument effectif dont l'indice commence a 0 a leurs matrices dont l'indice commence a 1 */
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv; /* temp; */
  indxc=ivector(1,n);
  indxr=ivector(1,n);  
  ipiv=ivector(1,n);
  for(j=1; j<=n; j++) ipiv[j]=0;
  for(i=1; i<=n; i++)
    {
      big=0.0;
      for(j=1; j<=n; j++)
	if(ipiv[j]!=1)
	  for(k=1; k<=n; k++)
	    {
	      if(ipiv[k]==0)
		{
		  if(fabs(a[j-1][k-1])>=big)
		    {
		      big=fabs(a[j-1][k-1]);
		      irow=j;
		      icol=k;
		    }
		}
	      else if(ipiv[k]>1) printf("g_j : matrice singuliere-1");
	    }	   
      ++(ipiv[icol]);
      if(irow!=icol)
	{
	  for(l=1; l<=n; l++) swap(&a[irow-1][l-1],&a[icol-1][l-1]);
	  for(l=1; l<=m; l++) swap(&b[irow-1][l-1],&b[icol-1][l-1]);
	}
      indxr[i]=irow;
      indxc[i]=icol;
      if(a[icol-1][icol-1]==0.0 ) printf("g_j : matrice singuliere-2");
      pivinv=1.0/a[icol-1][icol-1];
      a[icol-1][icol-1]=1.0;
      for(l=1; l<=n; l++) a[icol-1][l-1]*=pivinv;
      for(l=1; l<=m; l++) b[icol-1][l-1]*=pivinv;
      for(ll=1; ll<=n; ll++)
	if(ll!=icol)
	  {
	    dum=a[ll-1][icol-1];
	    a[ll-1][icol-1]=0.0;
	    for(l=1; l<=n; l++) a[ll-1][l-1]-=a[icol-1][l-1]*dum;
	    for(l=1; l<=m; l++) b[ll-1][l-1]-=b[icol-1][l-1]*dum;
	  }
    }
  for(l=n; l>=1; l--)
    {
      if(indxr[l]!=indxc[l]) for(k=1; k<=n; k++) swap(&a[k-1][indxr[l]-1],&a[k-1][indxc[l]-1]);
    }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}
/* ***************************************************************************** */
/* Egalite de deux matrices colonnes */
void egal(double u[],double up[],int n)
{
  int i;
  for(i=0; i<n; i++) u[i]=up[i];
}
/* ----------------------------------------------------------- */
/* Fait le produit de matrices AB, A matrice m,n, B matrice colonne n */
void pmc(double **a,double *b,double *c,int m,int n)
{
  int i,j; double s;
  for(i=0; i<m; i++)
    { 
      for(s=0.,j=0; j<n; j++) s=s+a[i][j]*b[j]; 
      c[i]=s;
    }
}
/* ----------------------------------------------------------- */
void mat_inv(double **a,double **inv_a,int n)
     /* Calcul l'inverse d'une matrice n,n */
     /* Rappel : g_j(x,n,y,n) renvoit l'inverse du x d'entree dans le x de sortie */
{
  int i,j;
  double **b;
  b=D_2(n,n);
  zero(b,n,n); /* seconds membres bidons */
  for(i=0; i<n; i++) for(j=0; j<n; j++) inv_a[i][j]=a[i][j]; /* sauvegarde de a[][] */
  g_j(inv_a,n,b,n);
}
/* ********************************************************************** */
