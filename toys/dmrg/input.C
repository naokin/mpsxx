#include <iostream>
#include <iomanip>
#include "input.h"

std::ostream& operator<< (std::ostream& ost, const DmrgInput& input)
{
  using std::endl;
  ost.setf(std::ios::fixed, std::ios::floatfield);
  ost.precision(8);
  ost << "\t##################################################"
      <<   "##################################################"       << endl;
  ost << "\t\tPRINTING INPUT PARAMETERS"                              << endl;
  ost << "\t\t\tL  = " << input.N_sites                               << endl;
  ost << "\t\t\tSz = " << static_cast<double>(input.N_phys_index-1)/2 << endl;
  ost << "\t\t\tM  = " << input.N_max_states                          << endl;
  ost << "\t\t\tJ  = " << input.J                                     << endl;
  ost << "\t\t\tJz = " << input.Jz                                    << endl;
  ost << "\t\t\thz = " << input.hz                                    << endl;
  ost << "\t##################################################"
      <<   "##################################################"       << endl;
  return ost;
}
