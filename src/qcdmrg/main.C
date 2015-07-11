#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include <qcdmrg/qcdmrg.h>
#include <utility/pout.h>

void print_help_msg ()
{
  using std::cout;
  using std::endl;

  cout << " usage qcdmrg [-i file] [-o file] [-s path] " << endl;
  cout << "   -i file : read parameters from input (searching ""conf"" by default) " << endl;
  cout << "   -o file : specify output file name (stdout is used by default) " << endl;
  cout << "   -s path : specify scratch directory (current directory by default) " << endl;

  return;
}

int main(int argc, char* argv[])
{
#ifdef _ENABLE_MPI
  boost::mpi::environment env;
#endif // _ENABLE_MPI

  using std::cout;
  using std::endl;

  std::string f_inp = "conf";
  std::string f_out;
  std::string p_scr = ".";

  if(argc == 1) { print_help_msg(); return 1; }

  for(int iarg = 1; iarg < argc; ++iarg)
  {
    if(strcmp(argv[iarg],"-i") == 0) f_inp = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) f_out = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) p_scr = argv[++iarg];
    if(strcmp(argv[iarg],"-h") == 0) print_help_msg();
  }

  std::ifstream ifs(f_inp.c_str());

  mpsxx::qcdmrg::DMRG_input input(ifs, p_scr);

  //
  // assign cout as alias to f_out
  //
  std::streambuf *backup;
  backup = cout.rdbuf();
  std::ofstream ofs;
  if(f_out.size() > 0) {
    ofs.open(f_out.c_str());
    cout.rdbuf(ofs.rdbuf());
  }

  pout << "\t****************************************************************************************************" << endl;
  pout << "\t\t\t\tMPSXX::PROTOTYPE::DMRG::OPTIMIZATION "                                                          << endl;
  pout << "\t****************************************************************************************************" << endl;

  //
  // dmrg optimization
  //
  mpsxx::qcdmrg::rundmrg(input);

  pout << endl;

  cout.rdbuf(backup);
  fout.close();

  return 0;
}
