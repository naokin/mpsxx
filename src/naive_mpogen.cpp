#include <iostream>
#include <iomanip>
#include <vector>
#include <cstring>

#include <btas/DENSE/TArray.h>
#include <btas/QSPARSE/QSTArray.h>

#include "mpidefs.h"

#include "parsing_integral.h"
#include "fileio.h"

#include "gen_qc_naive_mpos.h"
#include "compress_qc_mpos.h"

int main(int argc, char* argv[])
{
#ifndef _SERIAL
  boost::mpi::environment env(argc,argv);
#endif
  using std::cout;
  using std::endl;
  using namespace mpsxx;

  Communicator world;

  size_t iprint = 0;

  std::string f_dump = "FCIDUMP";
  std::string prefix = ".";
  std::string f_rord;
  std::string f_iout;

  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) f_dump = argv[++iarg];
    if(strcmp(argv[iarg],"-o") == 0) f_iout = argv[++iarg];
    if(strcmp(argv[iarg],"-r") == 0) f_rord = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) prefix = argv[++iarg];
    if(strcmp(argv[iarg],"-v") == 0) iprint = 1;
  }

  //
  // assign cout as alias to ost_iout
  //
  std::streambuf *backup;
  backup = cout.rdbuf();
  std::ofstream ost_iout;
  if(f_iout.size() > 0) {
    std::ostringstream oss;
    oss << f_iout << "." << world.rank();
    ost_iout.open(oss.str().c_str());
    cout.rdbuf(ost_iout.rdbuf());
  }

  int Norbs;
  int Nelec;
  double Ecore;
  btas::TArray<double,2> oneint;
  btas::TArray<double,4> twoint;

  cout << "\t****************************************************************************************************" << endl;
  cout << "\t\t\t\tMPSXX::PROTOTYPE::MPO GENERATOR FOR QC "                                                        << endl;
  cout << "\t****************************************************************************************************" << endl;
  cout << endl;
  cout << "\t====================================================================================================" << endl;
  cout << "\t\tLoading integral information from FCIDUMP "                                                         << endl;
  cout << "\t====================================================================================================" << endl;
  cout << endl;

  if(world.rank() == 0) {
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
  }
#ifndef _SERIAL
  boost::mpi::broadcast(world,Norbs,0);
  boost::mpi::broadcast(world,Nelec,0);
  boost::mpi::broadcast(world,Ecore,0);
  boost::mpi::broadcast(world,oneint,0);
  boost::mpi::broadcast(world,twoint,0);
#endif

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tGenerating QC MPOs from 1- and 2-particle integrals"                                                << endl;
  cout << "\t====================================================================================================" << endl;
  cout << endl;

  std::vector<int> groups;

  std::vector<std::vector<btas::QSTArray<double,4,fermion>>> mpos;
  groups = gen_qc_naive_mpos(Norbs,Ecore,oneint,twoint,mpos);

  std::cout << "\t\t" << std::setw(6) << mpos.size() << " operators have generated." << std::endl;
//std::fill(groups.begin(),groups.end(),1);

  std::vector<std::vector<btas::QSTArray<double,4,fermion>>> comp;
  groups = compress_qc_mpos(groups,mpos,comp);

  // deallocation
  std::vector<std::vector<btas::QSTArray<double,4,fermion>>>().swap(mpos);

  for(size_t i = 0; i < Norbs; ++i) {
    std::vector<btas::QSTArray<double,4,fermion>> impo(groups.size());
    for(size_t g = 0; g < groups.size(); ++g) {
      impo[g] = comp[g][i];
    }
    save(impo,getfile("mpo",prefix,i));
  }

  std::cout.rdbuf(backup);
  ost_iout.close();

  return 0;
}
