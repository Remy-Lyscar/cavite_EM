#include"bibli_fonctions.h"
//-------------------------------------------------------------------------------------------------------
double p_s(double *u,double *v)
{
  // Retourne le produit scalaire de u et v
  int i,n=3;
  double s;
  for(s=0.,i=0;i<n;i++) s+=u[i]*v[i];
  return s;
}
//-------------------------------------------------------------------------------------------------------
void p_v(double *w,double *u,double *v)
{
  // Fournit dans w les composantes cartésiennes du produit vectoriel u v
  int i,j,k,n=3;
  for(i=0;i<n;i++)
    {
      j=(i+1)%3; k=(i+2)%3;
      w[i]=u[j]*v[k]-u[k]*v[j];
    }
}
//-------------------------------------------------------------------------------------------------------
int p_v_n(double *w,double *u,double *v)
{
  // Vaut 0 si le produit vectoriel u v est nul, 1 sinon
  // S'il n'est pas nul fournit dans w les composantes cartésiennes du produit vectoriel u v normé à 1
  int i,j,k,n=3; double nor;
  for(nor=0,i=0;i<n;i++)
    {
      j=(i+1)%3; k=(i+2)%3;
      w[i]=u[j]*v[k]-u[k]*v[j];
      nor+=w[i]*w[i];
    }
  nor=sqrt(nor);
  if(nor==0.) return 0; 
  else{for(i=0;i<n;i++){w[i]/=nor;} return 1;}
}
//-------------------------------------------------------------------------------------------------------
int tet_phi_vect(double *u,double *tet,double *phi)
{
  // Calcule les angles sphériques theta et phi d'un vecteur dont on connaît les composantes cartésiennes
  double r,rr;
  r=sqrt(p_s(u,u)); if(r==0.) return 0;
  *tet=acos(u[2]/r);
  rr=r*sin(*tet); if(rr==0.) return 0;
  *phi=acos(u[0]/rr); 
  if(u[1]<0.) *phi=2.*M_PI-*phi;
  return 1;
}
//-------------------------------------------------------------------------------------------------------
// w=au produit d'un nombre par un vecteur
void s_v(double *w,double a,double *u)
{
  int i,n=3;
  for(i=0;i<n;i++) w[i]=a*u[i];
}
//-------------------------------------------------------------------------------------------------------
// w=au+bv combinaison linéaire de deux vecteurs
void c_l(double *w,double a,double *u,double b,double *v)
{
  int i,n=3;
  for(i=0;i<n;i++) w[i]=a*u[i]+b*v[i];
}
//-------------------------------------------------------------------------------------------------------
// w=u-v différence de deux vecteurs
void d_v(double *w,double *u,double *v)
{
  int i,n=3;
  for(i=0;i<n;i++) w[i]=u[i]-v[i];
}
//-------------------------------------------------------------------------------------------------------
