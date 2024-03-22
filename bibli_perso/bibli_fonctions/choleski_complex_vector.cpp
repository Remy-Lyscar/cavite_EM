#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <complex>
#include <vector>
using namespace std; // important de mettre cette ligne avant la suivante
#include "bibli_fonctions.h"
#define vcd vector<complex<double> >
int choleski_complex_vector(vcd a,vcd b,vcd c,vcd &x,vcd y,int n)
{ 
  complex<double> nu,zero(0.,0.); int i;
  vcd al(n); vcd be(n);
  if(b[0]==zero){printf("Echec 1 de choleski"); return 1;} 
  al[0]=-c[0]/b[0]; be[0]=y[0]/b[0]; 
  for(i=1; i<n; i++)
    {
      nu=a[i]*al[i-1]+b[i]; if(nu==zero){printf("Echec 2 de choleski"); return 2;}
      al[i]=-c[i]/nu; be[i]=(y[i]-a[i]*be[i-1])/nu; 
    }
  x[n-1]=be[n-1];
  for(i=n-2; i>=0; i--) x[i]=al[i]*x[i+1]+be[i];
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
  vcd vp(p,p+K); vcd vq(q,q+K); vcd vr(r,r+K); vcd vv(v,v+K); 
  vcd vu(K);
  choleski_complex_vector(vp,vq,vr,vu,vv,K);
  for(i=0; i<K; i++) cout << vu[i] << endl;
  // Verification :
  cout << vq[0]*vu[0]+vr[0]*vu[1] << vv[0] << endl;
  for(i=1; i<K-1; i++) cout << vp[i]*vu[i-1]+vq[i]*vu[i]+vr[i]*vu[i+1] << vv[i] << endl;
  cout << vp[K-1]*vu[K-2]+vq[K-1]*vu[K-1] << vv[K-1] << endl; 
  return(0);
}
#endif





