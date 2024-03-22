#include"bibli_fonctions.h"
void cercle(double xc,double yc,double r,double tmin,double tmax,fstream &fich)
{
  int i,np=360; double tet,dtet;
  dtet=(tmax-tmin)/(np-1);
  fich << endl;
  for(i=0;i<np;i++)
    {
      tet=tmin+i*dtet;
      fich << xc+r*cos(tet) << " " << yc+r*sin(tet) << endl;
    }
}
