#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "GAInput.h"
#include "GAOptimize.h"
#include "fiedler.h"
using namespace std;

#ifndef SERIAL
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#endif

int main(int argc, char* argv[])
{
#ifndef SERIAL
  mpi::environment env(argc, argv);
  mpi::communicator world;
  if(world.rank() == 0) cout << "Parallel GA simulation" << endl;
#endif

  string confFileName;
  string dumpFileName;

  for(int i = 1; i < argc; ++i)
  {
    if(strcmp(argv[i], "-config")   == 0) confFileName = argv[++i];
    if(strcmp(argv[i], "-integral") == 0) dumpFileName = argv[++i];
  }

  ifstream confFile(confFileName.c_str());
  ifstream dumpFile(dumpFileName.c_str());

  std::string tmp;
  getline(dumpFile,tmp);

  vector<int> num;

  int ind = 0;

  while(tmp[ind + 12] != ','){

     num.push_back(tmp[ind + 12] - 48);

     ind++;

  }

  int n = 0;

  for(int i = 0;i < num.size();++i)
     n += num[i]*pow(10,num.size() - 1 - i);

  std::vector<int> fiedlerv(n);

  for(int i = 0;i < n;++i)
     fiedlerv[i] = i;

  genetic::Cell final = genetic::gaordering(confFile, dumpFile, fiedlerv);

#ifndef SERIAL
  if(world.rank() == 0)
#endif
  {
    cout << "##################### MINIMUM GENE REP. #####################" << endl;
    cout << "Gene with MinValue = " << final << endl;
    cout << "Effective Distance = " << sqrt(final.Fitness()) << endl;

    cout << "#################### DMRG REORDER FORMAT ####################" << endl;

    vector<int> gaorder(final.Gen().Sequence());

    ofstream out("reorder.dat");

    for(int i = 0; i < n - 1; ++i) 
       out << gaorder[i] + 1 << ",";
    out  << gaorder[n - 1]  + 1 << endl;
  }

  return 0;
}
