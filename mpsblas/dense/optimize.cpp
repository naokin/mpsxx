#include <iostream>
#include <iomanip>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include "optimize.h"
#include "initialize.hpp"
#include "davidson.hpp"
#include "driver.hpp"

void check_mpos (const mpsxx::MPOs<double>& mpos)
{
  using namespace mpsxx;

  std::cout.precision(2);
  std::cout.setf(std::ios::fixed,std::ios::floatfield);
  for(size_t k = 0; k < mpos.size(); ++k) {

    std::cout << "\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << "\tSITE [ " << std::setw(2) << k << " ] " << std::endl;
    std::cout << "\t~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    for(size_t i = 0; i < mpos[k].shape(0); ++i) {
      for(size_t p = 0; p < mpos[k].shape(1); ++p) {
        std::cout << "\t";
        for(size_t j = 0; j < mpos[k].shape(3); ++j) {
          for(size_t q = 0; q < mpos[k].shape(2); ++q) {
            std::cout << std::setw(6) << mpos[k](i,p,q,j);
          }
          std::cout << "  ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
}

double mpsxx::optimize_1site (
        bool forward,
  const mpsxx::MPO  <double>& mpo0,
  const mpsxx::BLOCK<double>& lop0,
  const mpsxx::BLOCK<double>& rop0,
        mpsxx::MPS  <double>& mps0,
        mpsxx::MPS  <double>& mps1,
        mpsxx::BLOCK<double>& xop1)
{
  boost::function<void(const MPS<double>&, MPS<double>&)>
  fmult = boost::bind(compute_sigma_vector<double>,mpo0,lop0,rop0,_1,_2);

  MPS<double> diag;
  compute_diagonal_elements(mpo0,lop0,rop0,diag);

  double energy = davidson::diagonalize(fmult,diag,mps0);

  btas::TArray<double,2> g;
  btas::TArray<double,3> temp(mps0);
  canonicalize(forward,temp,mps0,g,0);

  temp = mps1;
  mps1.clear();
  xop1.clear();
  if(forward) {
    btas::Gemm(CblasNoTrans,CblasNoTrans,1.0,g,temp,1.0,mps1);
    renormalize(1,mpo0,lop0,mps0,mps0,xop1);
  }
  else {
    btas::Gemm(CblasNoTrans,CblasNoTrans,1.0,temp,g,1.0,mps1);
    renormalize(0,mpo0,rop0,mps0,mps0,xop1);
  }

  return energy;
}

double mpsxx::sweep (
  const mpsxx::MPOs<double>& mpos,
        mpsxx::MPSs<double>& mpss,
        std::vector<mpsxx::BLOCK<double>>& lops,
        std::vector<mpsxx::BLOCK<double>>& rops,
  const size_t& M)
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::fixed;

  size_t N = mpos.size();

  double emin = 1.0e8;

  cout.precision(16);

  // fowrad sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  for(int i = 0; i < N-1; ++i) {
    // diagonalize
    double eswp = optimize_1site(1,mpos[i],lops[i],rops[i],mpss[i],mpss[i+1],lops[i+1]);
    if(eswp < emin) emin = eswp;
    // print result
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
  }

  // backward sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  for(int i = N-1; i > 0; --i) {
    // diagonalize
    double eswp = optimize_1site(0,mpos[i],lops[i],rops[i],mpss[i],mpss[i-1],rops[i-1]);
    if(eswp < emin) emin = eswp;
    // print result
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
  }

  return emin;
}

double mpsxx::dmrg (
  const mpsxx::MPOs<double>& mpos,
        mpsxx::MPSs<double>& mpss,
  const size_t& M)
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::fixed;
  using std::scientific;

//check_mpos(mpos);

  std::vector<BLOCK<double>> lops;
  std::vector<BLOCK<double>> rops;

  initialize(mpos,mpss,lops,rops,M);

  size_t MX_ITER = 100;

  double esav = 1.0e8;

  cout.precision(16);

  for(int iter = 0; iter < MX_ITER; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = mpsxx::sweep(mpos,mpss,lops,rops,M);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << "\t====================================================================================================" << endl;
    cout << endl;
    esav = eswp;
    if(fabs(edif) < 1.0e-8) break;
  }

  return esav;
}
