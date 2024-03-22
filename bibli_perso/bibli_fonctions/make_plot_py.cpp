#include "bibli_fonctions.h"
void make_plot_py(ostringstream &ss)
{
  ostringstream run_py;
  run_py
    << "from numpy import *\n"
    << "from matplotlib.pyplot import *\n"
    << "from mpl_toolkits.mplot3d import Axes3D\n"
    << ss.str()
    << "show()"
    ;
  Py_Initialize();
  PyRun_SimpleString(run_py.str().c_str());
  Py_Finalize();
}
