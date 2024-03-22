#include "bibli_fonctions.h"
void multi_courbes_py(string s,int n)
{
ostringstream pyth;
pyth
<< "A=loadtxt('" << s << "')\n"
<< "for i in range(1," << n+1 << ") :\n"
<< "    plot(A[:,0],A[:,i])\n"
;
make_plot_py(pyth);
}
#if 0
int main()
{
  int n=102;
  multi_courbes_py("exemple_multi_courbes_py.dat",n);
  return 0;
}
#endif
