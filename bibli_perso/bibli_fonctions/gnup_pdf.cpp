#include"bibli_fonctions.h"
void gnup_pdf(string sss)
{
  fstream gnu("gnup_pdf.gnu",ios::out);
  gnu << "set nokey" << endl;
  gnu << sss.c_str() << endl;
  gnu << "pause -1" << endl;
  gnu << "set term pdf" << endl; 
  gnu << "set output 'gnup_pdf.pdf'" << endl;
  gnu << "replot" << endl;
  gnu.close();
  system("gnuplot gnup_pdf.gnu");
  //  system("rm gnup_pdf.gnu"); 
}
