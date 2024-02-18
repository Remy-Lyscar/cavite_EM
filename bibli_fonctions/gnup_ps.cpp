#include"bibli_fonctions.h"
void gnup_ps(string sss)
{
  fstream gnu("gnup_ps.gnu",ios::out);
  gnu << "set nokey" << endl;
  gnu << sss.c_str() << endl;
  gnu << "pause -1" << endl;
  gnu << "set term postscript" << endl; 
  gnu << "set output 'gnup_ps.eps'" << endl;
  gnu << "replot" << endl;
  gnu.close();
  system("gnuplot gnup_ps.gnu");
  //  system("rm gnup_ps.gnu"); 
}
