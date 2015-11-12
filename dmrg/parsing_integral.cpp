//
// parsing integral files
//

#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <boost/algorithm/string.hpp>

#include "mpidefs.h"
#include "parsing_integral.h"

std::vector<std::string> gettoken(std::ifstream& fin)
{
  // read msg from fin
  std::string msg;
  std::getline(fin, msg);
  boost::trim(msg);
  // split msg into token
  std::vector<std::string> tok;
  boost::split(tok, msg, boost::is_any_of("=, \t"), boost::token_compress_on);
  return tok;
}

//
// reading orbital ordering
void parsing_reorder
(std::ifstream& frord, std::vector<int>& reorder)
{
  Communicator world;
  if(world.rank() == 0) {
    std::vector<std::string> tok;
    tok = gettoken(frord);
    std::vector<int> reindex;
    for(int i = 0; i < tok.size(); ++i) reindex.push_back(atoi(tok[i].c_str())-1);
    reorder.resize(reindex.size(), 0);
    for(int i = 0; i < reindex.size(); ++i) reorder[reindex[i]] = i;
  }
#ifndef _SERIAL
  boost::mpi::broadcast(world,reorder,0);
#endif
}

//
// reading 1-particle integrals from general format
//
void parsing_oneint
(std::ifstream& fint1, int& norbs, btas::TArray<double,2>& oneint)
{
  fint1 >> norbs;
  oneint.resize(norbs, norbs);
  oneint.fill(0.0);
  size_t ix, jx;
  double value;
  while(fint1 >> ix >> jx >> value) {
    oneint(ix, jx) = value;
  }
}

void parsing_oneint
(std::ifstream& fint1, int& norbs, btas::TArray<double,2>& oneint, const std::vector<int>& reorder)
{
  fint1 >> norbs;
  oneint.resize(norbs, norbs);
  oneint.fill(0.0);
  size_t ix, jx;
  double value;
  while(fint1 >> ix >> jx >> value) {
    size_t ir = reorder[ix];
    size_t jr = reorder[jx];
    oneint(ir, jr) = value;
  }
}

//
// reading 2-particle integrals from general format
//
void parsing_twoint
(std::ifstream& fint2, int& norbs, btas::TArray<double,4>& twoint)
{
  fint2 >> norbs;
  twoint.resize(norbs, norbs, norbs, norbs);
  twoint.fill(0.0);
  size_t ix, jx, kx, lx;
  double value;
  while(fint2 >> ix >> jx >> kx >> lx >> value) {
    twoint(ix, jx, kx, lx) = value;
  }
}

void parsing_twoint
(std::ifstream& fint2, int& norbs, btas::TArray<double,4>& twoint, const std::vector<int>& reorder)
{
  fint2 >> norbs;
  twoint.resize(norbs, norbs, norbs, norbs);
  twoint.fill(0.0);
  size_t ix, jx, kx, lx;
  double value;
  while(fint2 >> ix >> jx >> kx >> lx >> value) {
    size_t ir = reorder[ix];
    size_t jr = reorder[jx];
    size_t kr = reorder[kx];
    size_t lr = reorder[lx];
    twoint(ir, jr, kr, lr) = value;
  }
}

bool is_local_int1e (const int& norbs, const int& i, const int& j)
{
  Communicator world;
  size_t iproc = world.rank();
  size_t nproc = world.size();
  size_t irank = 0;

  int n = 1+norbs/2;
//float f = sqrt(static_cast<float>(nproc));
//int n = static_cast<int>(norbs*f/(1.0+f));
  if(i > j) {
    if(i >= n && j < n) {
      irank = j %nproc;
    }
    else {
      int ij = i*(i+1)/2+j;
      irank = ij%nproc;
    }
  }
  else {
    if(j >= n && i < n) {
      irank = i %nproc;
    }
    else {
      int ji = j*(j+1)/2+i;
      irank = ji%nproc;
    }
  }

  return (iproc == irank);
}

bool is_local_int2e (const int& norbs, const int& i, const int& j, const int& k, const int& l)
{
  Communicator world;
  size_t iproc = world.rank();
  size_t nproc = world.size();
  size_t irank = 0;

  std::vector<int> index = { i, j, k, l };
  std::sort(index.begin(),index.end());
  const int& ix = index[0];
  const int& jx = index[1];
  const int& kx = index[2];
  const int& lx = index[3];

  int n = 1+norbs/2;
//float f = sqrt(static_cast<float>(nproc));
//int n = static_cast<int>(norbs*f/(1.0+f));
  if     (jx >= n) {
    irank = ix%nproc;
  }
  else if(kx >= n) {
    int ji = jx*(jx+1)/2+ix;
    irank = ji%nproc;
  }
  else {
    irank = lx%nproc;
  }

  return (iproc == irank);
}

//
// reading integrals from molpro FCIDUMP file
//
void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::TArray<double,2>& oneint, btas::TArray<double,4>& twoint)
{
  Communicator world;

  std::vector<std::string> tok;
  norbs = 0;
  nelec = 0;
  // get first line
  tok = gettoken(fdump);
  for(int i = 0; i < tok.size(); ++i) {
    if(tok[i] == "NORB" ) norbs = atoi(tok[++i].c_str());
    if(tok[i] == "NELEC") nelec = atoi(tok[++i].c_str());
  }
  // allocate integral array
  oneint.resize(norbs, norbs);
  oneint.fill(0.0);
  twoint.resize(norbs, norbs, norbs, norbs);
  twoint.fill(0.0);
  // find &END control
  while(!fdump.eof()) {
    tok = gettoken(fdump);
    if(tok[0] == "&END" || tok[0] == "/") break;
  }
  // read integrals
  ecore = 0.0;
  while((tok = gettoken(fdump)).size() > 1) {
    double value = atof(tok[0].c_str());
    int i = atoi(tok[1].c_str()) - 1;
    int j = atoi(tok[2].c_str()) - 1;
    int k = atoi(tok[3].c_str()) - 1;
    int l = atoi(tok[4].c_str()) - 1;

    if(i <  0 && j <  0 && world.rank() == 0) {
      ecore = value;
    }
    else if(k <  0 && l <  0) {
      if(is_local_int1e(norbs,i,j)) {
        oneint(i, j) = value;
        oneint(j, i) = value;
      }
    }
    else {
      if(is_local_int2e(norbs,i,j,k,l)) {
        twoint(i, k, j, l) = value;
        twoint(i, l, j, k) = value;
        twoint(j, k, i, l) = value;
        twoint(j, l, i, k) = value;
        twoint(k, i, l, j) = value;
        twoint(k, j, l, i) = value;
        twoint(l, i, k, j) = value;
        twoint(l, j, k, i) = value;
      }
    }
  }
}

//
// reading integrals from molpro FCIDUMP file with reordering
//
void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::TArray<double,2>& oneint, btas::TArray<double,4>& twoint, const std::vector<int>& reorder)
{
  std::vector<std::string> tok;
  norbs = 0;
  nelec = 0;
  // get first line
  tok = gettoken(fdump);
  for(int i = 0; i < tok.size(); ++i) {
    if(tok[i] == "NORB" ) norbs = atoi(tok[++i].c_str());
    if(tok[i] == "NELEC") nelec = atoi(tok[++i].c_str());
  }
  // allocate integral array
  oneint.resize(norbs, norbs);
  twoint.resize(norbs, norbs, norbs, norbs);
  // find &END control
  while(!fdump.eof()) {
    tok = gettoken(fdump);
    if(tok[0] == "&END" || tok[0] == "/") break;
  }
  // read integrals
  while((tok = gettoken(fdump)).size() > 1) {
    double value = atof(tok[0].c_str());
    int i = atoi(tok[1].c_str()) - 1;
    int j = atoi(tok[2].c_str()) - 1;
    int k = atoi(tok[3].c_str()) - 1;
    int l = atoi(tok[4].c_str()) - 1;
    if(i <  0 && j <  0) {
      ecore = value;
    }
    else if(k <  0 && l <  0) {
      int ix = reorder[i];
      int jx = reorder[j];
      oneint(ix, jx) = value;
      oneint(jx, ix) = value;
    }
    else {
      int ix = reorder[i];
      int jx = reorder[j];
      int kx = reorder[k];
      int lx = reorder[l];
      twoint(ix, kx, jx, lx) = value;
      twoint(ix, lx, jx, kx) = value;
      twoint(jx, kx, ix, lx) = value;
      twoint(jx, lx, ix, kx) = value;
      twoint(kx, ix, lx, jx) = value;
      twoint(kx, jx, lx, ix) = value;
      twoint(lx, ix, kx, jx) = value;
      twoint(lx, jx, kx, ix) = value;
    }
  }
}

//
// writing integrals in molpro FCIDUMP format
//
void writing_fcidump
(std::ofstream& fdump, const int& norbs, const int& nelec, const double& ecore, const btas::TArray<double,2>& oneint, const btas::TArray<double,4>& twoint)
{
  using std::setw;
  using std::endl;
  fdump << " &FCI NORB=" << setw(3) << norbs << ",NELEC=" << setw(2) << nelec << ",MS2= 0," << endl;
  fdump << "  ORBSYM=";
  for(int i = 0; i < norbs; ++i) fdump << "1,";
  fdump << endl;
  fdump << "  ISYM=1" << endl;
  fdump << " &END" << endl;
  fdump.precision(20);
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j <= i; ++j) {
      for(int k = 0; k <= i; ++k) {
        int lmax = (i == k) ? j : k;
        for(int l = 0; l <= lmax; ++l) {
          const double& val = twoint(i, k, j, l);
          if(std::fabs(val) > 1.0e-16) {
            fdump << setw(28) << std::scientific << val
                  << setw(4) << i+1 << setw(4) << j+1 << setw(4) << k+1 << setw(4) << l+1 << endl;
          }
        }
      }
    }
  }
  for(int i = 0; i < norbs; ++i) {
    for(int j = 0; j <= i; ++j) {
      const double& val = oneint(i, j);
      if(std::fabs(val) > 1.0e-16) {
        fdump << setw(28) << std::scientific << val << setw(4) << i+1 << setw(4) << j+1 << "   0   0" << endl;
      }
    }
  }
  fdump << setw(28) << std::scientific << ecore << "   0   0   0   0" << endl;
}
