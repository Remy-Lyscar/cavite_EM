// Cette fonction sert à tracer des histogrammes simples yi=f(xi) avec Python pour np valeurs xi équi-espacées de xmin à xmax : xi=xmin+(xmax-xmin)/(np-1)*i avec 0 <= i < np
// Les rectangles de l'histogramme sont centrés sur les valeurs de xi, les marges sont calculées automatiquement
// Il faut fournir :
// - xmin et xmax
// - np
// - un tableau contenant les yi
#include "bibli_fonctions.h"
void histo_simple_py(double xmin,double xmax,int np,double *y)
{
  int i;
  double x,dx,dxs2,ymin,ymax,dy,mr=20.; // mr détermine la marge au-dessous de ymin et au-dessus de ymax pour que l'histogramme ne touche pas le cadre du graphe
  dx=(xmax-xmin)/(np-1); dxs2=dx/2.;
  // Détermination de ymin et ymax
  ymin=y[0]; ymax=y[0]; 
  for(i=1;i<np;i++){if(y[i]<ymin) ymin=y[i]; if(y[i]>ymax) ymax=y[i];}
  dy=(ymax-ymin)/mr;
  // Écriture des points de contour de l'histogramme dans un fichier
  fstream res("histo_simple_py.res",ios::out);
  for(i=0;i<np;i++)
    {
      x=xmin+i*dx;
      res << x-dxs2 << " " << y[i] << endl; 
      res << x+dxs2 << " " << y[i] << endl; 
    }
  res.close();
  // Tracé par Python à partir du fichier, ajustement automatique des marges
  ostringstream pyth;
  pyth
    << "l=loadtxt('histo_simple_py.res')\n"
    << "x=l[:,0]\n"
    << "y=l[:,1]\n"
    << "plot(x,y)\n"
    << "xlim(" << xmin-dx << "," << xmax+dx << ")\n"
    << "ylim(" << ymin-dy << "," << ymax+dy << ")\n"
    ;
  make_plot_py(pyth);
  system("rm histo_simple_py.res");
}
/*
int main()
{
  int np=8; double xmin=1.,xmax=8.;
  double *y=D_1(np);                 
  ini_D_1(y,np,1.,4.,2.6,-8.,3.,5.,9.,0.);
  histo_simple_py(xmin,xmax,np,y);
  return 0;
}
*/

