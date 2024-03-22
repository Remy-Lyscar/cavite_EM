#include"bibli_fonctions.h"
void gnup(string sss)
{
  fstream gnu("gnup.gnu",ios::out);
  gnu << "set nokey" << endl;
  gnu << sss.c_str() << endl;
  gnu << "pause -1" << endl;
  gnu.close();
  system("gnuplot gnup.gnu");
  //  system("rm gnup.gnu"); 
}
