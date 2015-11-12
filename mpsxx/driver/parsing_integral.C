//
// parsing integral files
//

#include <iomanip>
#include <string>
#include <cmath>
#include <boost/algorithm/string.hpp>

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
  std::vector<std::string> tok;
  tok = gettoken(frord);
  std::vector<int> reindex;
  for(int i = 0; i < tok.size(); ++i) reindex.push_back(atoi(tok[i].c_str())-1);
  reorder.resize(reindex.size(), 0);
  for(int i = 0; i < reindex.size(); ++i) reorder[reindex[i]] = i;
}

//
// reading 1-particle integrals from general format
//
void parsing_oneint
(std::ifstream& fint1, int& norbs, btas::DArray<2>& oneint)
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
(std::ifstream& fint1, int& norbs, btas::DArray<2>& oneint, const std::vector<int>& reorder)
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
(std::ifstream& fint2, int& norbs, btas::DArray<4>& twoint)
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
(std::ifstream& fint2, int& norbs, btas::DArray<4>& twoint, const std::vector<int>& reorder)
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

//
// reading integrals from molpro FCIDUMP file
//
void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::DArray<2>& oneint, btas::DArray<4>& twoint)
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
      oneint(i, j) = value;
      oneint(j, i) = value;
    }
    else {
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

//
// reading integrals from molpro FCIDUMP file with reordering
//
void parsing_fcidump
(std::ifstream& fdump, int& norbs, int& nelec, double& ecore, btas::DArray<2>& oneint, btas::DArray<4>& twoint, const std::vector<int>& reorder)
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
(std::ofstream& fdump, const int& norbs, const int& nelec, const double& ecore, const btas::DArray<2>& oneint, const btas::DArray<4>& twoint)
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
