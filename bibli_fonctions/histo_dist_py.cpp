// Cette fonction trace un histogramme de distribution statistique avec Python
// La v. a. dont on trace la distribution est comprise entre les valeurs xmin et xmax entre lesquelles il y a nt tranches d'égales largeur dx, donc dx=(xmax-xmin)/nt 
// Les rectangles de l'histogramme sont centrés sur les valeurs de xmin+i*dx+dx/2 0<=i<=nt-1, les marges sont calculées automatiquement
// Il faut fournir :
// - xmin et xmax
// - nt
// - un tableau de doubles contenant le nombre de valeurs dans chaque tranche i 0<=i<=nt-1
// - une fonction de type double et ayant un seul argument double représentant la distribution théorique si on veut tracer aussi cette dernière, la valeur NULL sinon
// - le nombre de points souhaité pour cette distribution théorique (n'importe quelle valeur, non prise en compte, si l'argument précédent est NULL) : np points équi-espacés de xmin à xmax inclus donc ddx=(xmax-xmin)/(np-1)
#include "bibli_fonctions.h"
void histo_dist_py(double xmin,double xmax,int nt,double *y,double (*dist_th)(double),int np)
{
  int i;
  double x,dx,ddx,ymin,ymax,dy,mr=20.,s; // mr détermine la marge au-dessous de ymin et au-dessus de ymax pour que l'histogramme ne touche pas le cadre du graphe
  dx=(xmax-xmin)/nt;
  // Détermination de ymin et ymax
  ymin=y[0]; ymax=y[0]; 
  for(i=1;i<nt;i++){if(y[i]<ymin) ymin=y[i]; if(y[i]>ymax) ymax=y[i];}
  dy=(ymax-ymin)/mr;
  // Écriture des points du contour de l'histogramme dans un fichier
  fstream res1("histo_dist_py_1.res",ios::out);
  for(i=0;i<nt;i++)
    { 
      x=xmin+i*dx;
      res1 << x << " " << y[i] << endl; 
      res1 << x+dx << " " << y[i] << endl; 
    }
  res1.close();
  // Tracé par Python à partir du fichier, ajustement automatique des marges
  ostringstream pyth;
  pyth
    << "l1=loadtxt('histo_dist_py_1.res')\n"
    << "x1=l1[:,0]\n"
    << "y1=l1[:,1]\n"
    << "plot(x1,y1)\n"
    << "xlim(" << xmin-dx << "," << xmax+dx << ")\n"
    << "ylim(" << ymin-dy << "," << ymax+dy << ")\n"
    ;
  if(dist_th!=NULL)
    {
      // Écriture des points de la distribution théorique dans un autre fichier
      for(s=0.,i=0;i<nt;i++){s+=y[i];} s*=dx; // pour la normalisation
      //      cout << "s=" << s << endl;
      fstream res2("histo_dist_py_2.res",ios::out);
      ddx=(xmax-xmin)/(np-1);
      for(i=0;i<np;i++)
	{ 
	  x=xmin+i*ddx;
	  res2 << x << " " << s*dist_th(x) << endl; 
	}
      res2.close();
      pyth
	<< "l2=loadtxt('histo_dist_py_2.res')\n"
	<< "x2=l2[:,0]\n"
	<< "y2=l2[:,1]\n"
	<< "plot(x2,y2)\n"
	;
    }
  make_plot_py(pyth);
  if(dist_th!=NULL) system("rm histo_dist_py_1.res histo_dist_py_2.res"); else system("rm histo_dist_py_1.res"); 
}
/*
int main()
{
  int i,ni=1000000,t,nt=5,np=100; double x,xmin=0.,xmax=10.,dx=(xmax-xmin)/nt,u;
  l_=1.2;
  double *h=D_1(nt);
  for(i=1;i<=ni;i++)
    {
      do {u=drand48();} while(u==0.);
      x=-l_*log(u);
      t=floor((x-xmin)/dx);
      if(t<0) t=0; else if(t>nt-1) t=nt-1;
      h[t]+=1;
    }
  histo_dist_py(xmin,xmax,nt,h,dist,np);
  //  histo_dist_py(xmin,xmax,nt,h,NULL,np);
  return 0;
}
*/


