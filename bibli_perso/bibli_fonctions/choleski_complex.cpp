#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <complex>
using namespace std; // important de mettre cette ligne avant la suivante
#include "bibli_fonctions.h"
int choleski_complex(complex<double> *a,complex<double> *b,complex<double> *c,complex<double> *x,complex<double> *y,int n)
{ 
  complex<double> *al,*be,nu,zero(0.,0.); int i;
  al=new complex<double>[n]; be=new complex<double>[n];
  if(b[0]==zero){printf("Echec 1 de choleski"); return 1;} 
  al[0]=-c[0]/b[0]; be[0]=y[0]/b[0]; 
  for(i=1; i<n; i++)
    {
      nu=a[i]*al[i-1]+b[i]; if(nu==zero){printf("Echec 2 de choleski"); return 2;}
      al[i]=-c[i]/nu; be[i]=(y[i]-a[i]*be[i-1])/nu;
    }
  x[n-1]=be[n-1];
  for(i=n-2; i>=0; i--) x[i]=al[i]*x[i+1]+be[i];
  delete al; delete be;
  return(0);
}
#if 0
#define K 4
int main(void)
{
  int i;
  complex<double> I(0.,1.);
  complex<double> p[K]={0.+I*0.,6.+I*2.,4.+I*8.,3.+I*9.};
  complex<double> q[K]={4.+I*7.,3.+I*7.,6.+I*5.,4.+I*3.};
  complex<double> r[K]={5.+I*8.,4.+I*3.,3.+I*9.,0.+I*0.};
  complex<double> v[K]={1.+I*8.,4.+I*3.,2.+I*6.,3.+I*7.};
  complex<double> u[K];
  choleski_complex(p,q,r,u,v,K);
  //for(i=0; i<K; i++)   cout << "real=" << real(u[i]) << " imag=" << imag(u[i]) << endl; 
  for(i=0; i<K; i++) cout << i << " " << u[i] << endl;
  // Verification :
  cout << q[0]*u[0]+r[0]*u[1] << v[0] << endl;
  for(i=1; i<K-1; i++) cout << p[i]*u[i-1]+q[i]*u[i]+r[i]*u[i+1] << v[i] << endl;
  cout << p[K-1]*u[K-2]+q[K-1]*u[K-1] << v[K-1] << endl; 
  return(0);
}
#endif





