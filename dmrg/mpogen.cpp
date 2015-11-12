#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>

#include <legacy/DENSE/TArray.h>

#include "mpidefs.h"

#include "fileio.h"
#include "gen_qc_operators.h"
#include "parsing_integral.h"
#include "compress_qc_mpos.h"

int main (int argc, char* argv[])
{
#ifndef _SERIAL
  boost::mpi::environment env(argc,argv);
#endif
  using std::endl;
  using namespace mpsxx;

  Communicator world;

  int iprint = 0;
  bool enable_swap_sweep = false;
  bool do_compress = false;

  std::string f_dump = "FCIDUMP";
  std::string prefix = ".";
  std::string f_rord;
  std::string f_iout;

  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) f_dump = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) f_iout = argv[++iarg];
    if(strcmp(argv[iarg],"-r") == 0) f_rord = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) prefix = argv[++iarg];
    if(strcmp(argv[iarg],"-x") == 0) enable_swap_sweep = true;
    if(strcmp(argv[iarg],"-c") == 0) do_compress = true;
    if(strcmp(argv[iarg],"-v") == 0) iprint = 1;
  }

  //
  // assign cout as alias to ost_iout
  //
  std::streambuf *backup;
  backup = std::cout.rdbuf();
  std::ofstream ost_iout;
  if(f_iout.size() > 0) {
    std::ostringstream oss;
    oss << f_iout << "." << world.rank();
    ost_iout.open(oss.str().c_str());
    std::cout.rdbuf(ost_iout.rdbuf());
  }

  int Norbs;
  int Nelec;
  double Ecore;
  btas::TArray<double,2> oneint;
  btas::TArray<double,4> twoint;

  pout << "\t****************************************************************************************************" << endl;
  pout << "\t\t\t\tMPSXX::PROTOTYPE::MPO GENERATOR FOR QC "                                                        << endl;
  pout << "\t****************************************************************************************************" << endl;
  pout << endl;
  pout << "\t====================================================================================================" << endl;
  pout << "\t\tLoading integral information from FCIDUMP "                                                         << endl;
  pout << "\t====================================================================================================" << endl;
  pout << endl;

  std::ifstream ist_dump(f_dump.c_str());
  if(!f_rord.empty()) {
    std::ifstream ist_rord(f_rord.c_str());
    std::vector<int> reorder;
    parsing_reorder(ist_rord,reorder);
    parsing_fcidump(ist_dump,Norbs,Nelec,Ecore,oneint,twoint,reorder);
  }
  else {
    parsing_fcidump(ist_dump,Norbs,Nelec,Ecore,oneint,twoint);
  }

  pout << "\t====================================================================================================" << endl;
  pout << "\t\tGenerating QC MPOs from 1- and 2-particle integrals"                                                << endl;
  pout << "\t====================================================================================================" << endl;
  pout << endl;

  gen_qc_operators(Norbs,Ecore,oneint,twoint,"mpo",prefix,enable_swap_sweep);

  if(do_compress) compress_qc_mpos(Norbs,"mpo",prefix);

  std::cout.rdbuf(backup);
  ost_iout.close();

  return 0;
}
