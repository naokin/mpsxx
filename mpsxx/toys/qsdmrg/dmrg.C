#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "FermiQuantum.h"
namespace btas { typedef FermiQuantum Quantum; };

#include <legacy/QSDArray.h>
#include <legacy/QSDblas.h>
#include "btas_template_specialize.h"

#include "dmrg.h"
#include "driver.h"
#include "davidson.h"
using namespace btas;

//
// random number generator
//
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX) - 1.0; }

//
// Heisenberg model
//
void prototype::Heisenberg::construct_mpo(MpStorages& sites, const HeisenbergModel& params)
{
  int    Nz = params.Nz;
  double J  = params.J;
  double Jz = params.Jz;
  double Hz = params.Hz;

  int    L  = sites.size();
  int    d  = Nz + 1;
  double sz = static_cast<double>(Nz)/2;

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tCONSTRUCT MATRIX PRODUCT OPERATORS (MPOs) "                                                         << endl;
  cout.precision(4);
  cout << "\t\t\t+ coupling coefficient J  : " << setw(8) << fixed << J  << endl; 
  cout << "\t\t\t+ coupling coefficient Jz : " << setw(8) << fixed << Jz << endl;
  cout << "\t\t\t+ coupling coefficient Hz : " << setw(8) << fixed << Hz << endl;
  cout << "\t\t\t+ coupling coefficient Sz : " << setw(8) << fixed << sz << endl;
  cout << "\t====================================================================================================" << endl;

  Qshapes qp; // physical index
  for(int i = 0; i < d; ++i) {
    int iz = Nz - 2 * i;
    qp.push_back(Quantum(0, iz));
  }

  Qshapes qz; // 0 quantum number
  qz.push_back(Quantum(0,  0));

  Qshapes qi; // quantum index comes in
  qi.push_back(Quantum(0,  0)); // I
  qi.push_back(Quantum(0, +2)); // S- (from S+)
  qi.push_back(Quantum(0, -2)); // S+ (from S-)
  qi.push_back(Quantum(0,  0)); // Sz
  qi.push_back(Quantum(0,  0)); // I

  Qshapes qo; // quantum index comes out
  qo.push_back(Quantum(0,  0)); // I
  qo.push_back(Quantum(0, -2)); // S+ (to S-)
  qo.push_back(Quantum(0, +2)); // S- (to S+)
  qo.push_back(Quantum(0,  0)); // Sz
  qo.push_back(Quantum(0,  0)); // I

  // resize & set to 0
  sites[ 0 ].mpo.resize(Quantum::zero(), TinyVector<Qshapes, 4>(qz, qp,-qp, qo));
  for(int i = 1; i < L-1; ++i) {
    sites[i].mpo.resize(Quantum::zero(), TinyVector<Qshapes, 4>(qi, qp,-qp, qo));
  }
  sites[L-1].mpo.resize(Quantum::zero(), TinyVector<Qshapes, 4>(qi, qp,-qp, qz));

  // construct mpos for spin-hamiltonian
  double mz = sz;
  for(int m = 0; m < d; ++m) {
    // set block elements
    DArray<4> data_Id(1, 1, 1, 1); data_Id = 1.0;
    DArray<4> data_Mz(1, 1, 1, 1); data_Mz = mz;
    DArray<4> data_Hz(1, 1, 1, 1); data_Hz = Hz * mz;
    DArray<4> data_Jz(1, 1, 1, 1); data_Jz = Jz * mz;
    // insert blocks
    sites[ 0 ].mpo.insert(shape(0, m, m, 0), data_Hz); // Hz  Sz
    sites[ 0 ].mpo.insert(shape(0, m, m, 3), data_Jz); // Jz  Sz
    sites[ 0 ].mpo.insert(shape(0, m, m, 4), data_Id); //     I
    for(int i = 1; i < L-1; ++i) {
      sites[i].mpo.insert(shape(0, m, m, 0), data_Id); //     I
      sites[i].mpo.insert(shape(3, m, m, 0), data_Mz); //     Sz
      sites[i].mpo.insert(shape(4, m, m, 0), data_Hz); // Hz  Sz
      sites[i].mpo.insert(shape(4, m, m, 3), data_Jz); // Jz  Sz
      sites[i].mpo.insert(shape(4, m, m, 4), data_Id); //     I
    }
    sites[L-1].mpo.insert(shape(0, m, m, 0), data_Id); //     I
    sites[L-1].mpo.insert(shape(3, m, m, 0), data_Mz); //     Sz
    sites[L-1].mpo.insert(shape(4, m, m, 0), data_Hz); // Hz  Sz
    mz -= 1.0;
  }
  double mz_plus  = sz - 1.0;
  double mz_minus = sz;
  for(int m = 0; m < d - 1; ++m) {
    double c_plus  = sqrt(sz*(sz+1.0)-mz_plus *(mz_plus +1.0));
    double c_minus = sqrt(sz*(sz+1.0)-mz_minus*(mz_minus-1.0));
    // set block elements
    DArray<4> data_Sp(1, 1, 1, 1); data_Sp = c_plus;
    DArray<4> data_Sm(1, 1, 1, 1); data_Sm = c_minus;
    DArray<4> data_Jp(1, 1, 1, 1); data_Jp = c_plus  * J / 2;
    DArray<4> data_Jm(1, 1, 1, 1); data_Jm = c_minus * J / 2;
    // insert blocks
    sites[ 0 ].mpo.insert(shape(0, m,   m+1, 1), data_Jp); // J/2 S+
    sites[ 0 ].mpo.insert(shape(0, m+1, m,   2), data_Jm); // J/2 S-
    for(int i = 1; i < L-1; ++i) {
      sites[i].mpo.insert(shape(1, m+1, m,   0), data_Sm); //     S-
      sites[i].mpo.insert(shape(2, m,   m+1, 0), data_Sp); //     S+
      sites[i].mpo.insert(shape(4, m,   m+1, 1), data_Jp); // J/2 S+
      sites[i].mpo.insert(shape(4, m+1, m,   2), data_Jm); // J/2 S-
    }
    sites[L-1].mpo.insert(shape(1, m+1, m,   0), data_Sm); //     S-
    sites[L-1].mpo.insert(shape(2, m,   m+1, 0), data_Sp); //     S+
    mz_plus  -= 1.0;
    mz_minus -= 1.0;
  }
}

//
// Hubbard model
//
void prototype::Hubbard::construct_mpo(MpStorages& sites, const HubbardModel& params)
{
  double t  = params.t;
  double U  = params.U;

  int    L  = sites.size();

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tCONSTRUCT MATRIX PRODUCT OPERATORS (MPOs) "                                                         << endl;
  cout.precision(4);
  cout << "\t\t\t+ coupling coefficient t  : " << setw(8) << fixed << t  << endl; 
  cout << "\t\t\t+ coupling coefficient U  : " << setw(8) << fixed << U  << endl;
  cout << "\t====================================================================================================" << endl;

  Qshapes qp; // physical index { |0>, |a>, |b>, |ab> }
  qp.push_back(Quantum( 0,  0));
  qp.push_back(Quantum( 1,  1));
  qp.push_back(Quantum( 1, -1));
  qp.push_back(Quantum( 2,  0));

  Qshapes qz; // 0 quantum number
  qz.push_back(Quantum( 0,  0));

  Qshapes qi; // quantum index comes in
  qi.push_back(Quantum( 0,  0)); // 0
  qi.push_back(Quantum( 1,  1)); // a-
  qi.push_back(Quantum(-1, -1)); // a+
  qi.push_back(Quantum( 1, -1)); // b-
  qi.push_back(Quantum(-1,  1)); // b+
  qi.push_back(Quantum( 0,  0)); // I

  Qshapes qo; // quantum index comes out
  qo.push_back(Quantum( 0,  0)); // I
  qo.push_back(Quantum(-1, -1)); // a+
  qo.push_back(Quantum( 1,  1)); // a-
  qo.push_back(Quantum(-1,  1)); // b+
  qo.push_back(Quantum( 1, -1)); // b-
  qo.push_back(Quantum( 0,  0)); // 0

  // resize & set to 0
  sites[ 0 ].mpo.resize(Quantum::zero(), TinyVector<Qshapes, 4>(qz, qp,-qp, qo));
  for(int i = 1; i < L-1; ++i) {
    sites[i].mpo.resize(Quantum::zero(), TinyVector<Qshapes, 4>(qi, qp,-qp, qo));
  }
  sites[L-1].mpo.resize(Quantum::zero(), TinyVector<Qshapes, 4>(qi, qp,-qp, qz));

  // set block elements
  DArray<4> data_Ip(1, 1, 1, 1); data_Ip = 1.0;
  DArray<4> data_Im(1, 1, 1, 1); data_Im =-1.0;
  DArray<4> data_tp(1, 1, 1, 1); data_tp = t;
  DArray<4> data_tm(1, 1, 1, 1); data_tm =-t;
  DArray<4> data_Un(1, 1, 1, 1); data_Un = U;
  // insert blocks
  sites[ 0 ].mpo.insert(shape(0, 3, 3, 0), data_Un); //  U
  sites[ 0 ].mpo.insert(shape(0, 1, 0, 1), data_tp); //  t a+
  sites[ 0 ].mpo.insert(shape(0, 3, 2, 1), data_tm); //  t a+ [x(-1) due to P(a,b)]
  sites[ 0 ].mpo.insert(shape(0, 0, 1, 2), data_tm); // -t a-
  sites[ 0 ].mpo.insert(shape(0, 2, 3, 2), data_tp); // -t a- [x(-1) due to P(a,b)]
  sites[ 0 ].mpo.insert(shape(0, 2, 0, 3), data_tp); //  t b+
  sites[ 0 ].mpo.insert(shape(0, 3, 1, 3), data_tp); //  t b+
  sites[ 0 ].mpo.insert(shape(0, 0, 2, 4), data_tm); // -t b-
  sites[ 0 ].mpo.insert(shape(0, 1, 3, 4), data_tm); // -t b-
  sites[ 0 ].mpo.insert(shape(0, 0, 0, 5), data_Ip); //  I
  sites[ 0 ].mpo.insert(shape(0, 1, 1, 5), data_Ip); //  I
  sites[ 0 ].mpo.insert(shape(0, 2, 2, 5), data_Ip); //  I
  sites[ 0 ].mpo.insert(shape(0, 3, 3, 5), data_Ip); //  I
  for(int i = 1; i < L-1; ++i) {
    sites[i].mpo.insert(shape(0, 0, 0, 0), data_Ip); //  I
    sites[i].mpo.insert(shape(0, 1, 1, 0), data_Ip); //  I
    sites[i].mpo.insert(shape(0, 2, 2, 0), data_Ip); //  I
    sites[i].mpo.insert(shape(0, 3, 3, 0), data_Ip); //  I
    sites[i].mpo.insert(shape(1, 0, 1, 0), data_Ip); //  a-
    sites[i].mpo.insert(shape(1, 2, 3, 0), data_Im); //  a- [x(-1) due to P(a,b)]
    sites[i].mpo.insert(shape(2, 1, 0, 0), data_Ip); //  a+
    sites[i].mpo.insert(shape(2, 3, 2, 0), data_Im); //  a+ [x(-1) due to P(a,b)]
    sites[i].mpo.insert(shape(3, 0, 2, 0), data_Ip); //  b-
    sites[i].mpo.insert(shape(3, 1, 3, 0), data_Ip); //  b-
    sites[i].mpo.insert(shape(4, 2, 0, 0), data_Ip); //  b+
    sites[i].mpo.insert(shape(4, 3, 1, 0), data_Ip); //  b+
    sites[i].mpo.insert(shape(5, 3, 3, 0), data_Un); //  U
    sites[i].mpo.insert(shape(5, 1, 0, 1), data_tp); //  t a+
    sites[i].mpo.insert(shape(5, 3, 2, 1), data_tm); //  t a+ [x(-1) due to P(a,b)]
    sites[i].mpo.insert(shape(5, 0, 1, 2), data_tm); // -t a-
    sites[i].mpo.insert(shape(5, 2, 3, 2), data_tp); // -t a- [x(-1) due to P(a,b)]
    sites[i].mpo.insert(shape(5, 2, 0, 3), data_tp); //  t b+
    sites[i].mpo.insert(shape(5, 3, 1, 3), data_tp); //  t b+
    sites[i].mpo.insert(shape(5, 0, 2, 4), data_tm); // -t b-
    sites[i].mpo.insert(shape(5, 1, 3, 4), data_tm); // -t b-
    sites[i].mpo.insert(shape(5, 0, 0, 5), data_Ip); //  I
    sites[i].mpo.insert(shape(5, 1, 1, 5), data_Ip); //  I
    sites[i].mpo.insert(shape(5, 2, 2, 5), data_Ip); //  I
    sites[i].mpo.insert(shape(5, 3, 3, 5), data_Ip); //  I
  }
  sites[L-1].mpo.insert(shape(0, 0, 0, 0), data_Ip); //  I
  sites[L-1].mpo.insert(shape(0, 1, 1, 0), data_Ip); //  I
  sites[L-1].mpo.insert(shape(0, 2, 2, 0), data_Ip); //  I
  sites[L-1].mpo.insert(shape(0, 3, 3, 0), data_Ip); //  I
  sites[L-1].mpo.insert(shape(1, 0, 1, 0), data_Ip); //  a-
  sites[L-1].mpo.insert(shape(1, 2, 3, 0), data_Im); //  a- [x(-1) due to P(a,b)]
  sites[L-1].mpo.insert(shape(2, 1, 0, 0), data_Ip); //  a+
  sites[L-1].mpo.insert(shape(2, 3, 2, 0), data_Im); //  a+ [x(-1) due to P(a,b)]
  sites[L-1].mpo.insert(shape(3, 0, 2, 0), data_Ip); //  b-
  sites[L-1].mpo.insert(shape(3, 1, 3, 0), data_Ip); //  b-
  sites[L-1].mpo.insert(shape(4, 2, 0, 0), data_Ip); //  b+
  sites[L-1].mpo.insert(shape(4, 3, 1, 0), data_Ip); //  b+
  sites[L-1].mpo.insert(shape(5, 3, 3, 0), data_Un); //  U

  // taking operator parity P(O(l),<n|)
  std::vector<int> indx1(1, 0); // left mpo index [0]
  std::vector<int> indx2(1, 1); // phys bra index [1]
  for(int i = 0; i < L; ++i) {
    sites[i].mpo.parity(indx1, indx2);
  }
}

void prototype::set_quantum_blocks(const MpStorages& sites, const Quantum& qt, std::vector<Qshapes>& qb, int QMAX_SIZE)
{
  int L = sites.size();

  // physical index
  Qshapes qp;
  // 0 quantum number
  Qshapes qz(1, Quantum::zero());

  qb.resize(L);

  // quantum blocks from the entire Fock space
  qb[0] =-sites[0].mpo.qshape(2);
  for(int i = 1; i < L-1; ++i) {
    qp    =-sites[i].mpo.qshape(2);
    qb[i] = qb[i-1] & qp; // get unique elements of { q(left) x q(phys) }
  }
  qp      =-sites[L-1].mpo.qshape(2);
  qb[L-1] = Qshapes(1, qt);

  // reduce zero quantum blocks
  for(int i = L-1; i > 0; --i) {
    qp =-sites[i].mpo.qshape(2);
    Qshapes& ql = qb[i-1];
    Qshapes& qr = qb[i];

    // check non-zero for each ql index
    Qshapes::iterator lt = ql.begin();
    while(lt != ql.end()) {
      bool non_zero = false;
      for(int p = 0; p < qp.size(); ++p) {
        for(int r = 0; r < qr.size(); ++r) {
          non_zero |= (qr[r] == (qp[p] * (*lt)));
        }
      }
      if(non_zero)
        ++lt;
      else
        ql.erase(lt);
    }
    // further reduction
    if(QMAX_SIZE > 0 && ql.size() > QMAX_SIZE) {
      int offs = (ql.size() - QMAX_SIZE) / 2;
      ql = Qshapes(ql.begin()+offs, ql.begin()+offs+QMAX_SIZE);
    }
  }
}

void prototype::initialize(MpStorages& sites, const Quantum& qt, int M)
{
  int L = sites.size();

  // physical index
  Qshapes qp;
  Dshapes dp;
  // left state index
  Qshapes ql;
  Dshapes dl;
  // right state index
  Qshapes qr;
  Dshapes dr;
  // 0 quantum number
  Qshapes qz(1, Quantum::zero());
  Dshapes dz(qz.size(), 1);

  // non-zero quantum numbers for each site
  std::vector<Qshapes> qb;

  int max_size = 20;
  set_quantum_blocks(sites, qt, qb, max_size);

  //
  // create random wavefunction
  //

  int M0 = 1;
  int Mx = M;

  TinyVector<Qshapes, 3> qshape;
  TinyVector<Dshapes, 3> dshape;

  qr = qz;
  dr = Dshapes(qr.size(), 1);

  for(int i = 0; i < L-1; ++i) {
    // physical index is taken from mpo's ket index
    qp =-sites[i].mpo.qshape(2);
    dp = Dshapes(qp.size(), 1);
    // left index equals to previous right index
    ql = qr;
    dl = dr;
    // non-zero quantum numbers for site i
    qr = qb[i];
    dr = Dshapes(qr.size(), M0);

    qshape = TinyVector<Qshapes, 3>( ql, qp,-qr);
    dshape = TinyVector<Dshapes, 3>( dl, dp, dr);
    sites[i].wfnc[0].resize(Quantum::zero(), qshape, dshape);
    sites[i].wfnc[0] = rgen;
    sites[i].save(i);
  }

  qp =-sites[L-1].mpo.qshape(2);
  dp = Dshapes(qp.size(), 1);
  ql = qr;
  dl = dr;
  qr = qz;
  dr = dz;
  qshape = TinyVector<Qshapes, 3>( ql, qp,-qr);
  dshape = TinyVector<Dshapes, 3>( dl, dp, dr);
  sites[L-1].wfnc[0].resize(qt, qshape, dshape);
  sites[L-1].wfnc[0] = rgen;

  //
  // canonicalize & renormalize
  //

  qshape = TinyVector<Qshapes, 3>( qz, qz, qz);
  dshape = TinyVector<Dshapes, 3>( dz, dz, dz);

  sites[L-1].ropr[0].resize(Quantum::zero(), qshape, dshape);
  sites[L-1].ropr[0] = 1.0;

  for(int i = L-1; i > 0; --i) {
    sites[i-1].load(i-1);
    util::Normalize(sites[i].wfnc[0]);
    Canonicalize(0, sites[i].wfnc[0], sites[i].rmps[0], Mx);
    QSDcopy(sites[i-1].wfnc[0], sites[i-1].lmps[0]);
    sites[i-1].wfnc[0].clear();
    ComputeGuess(0, sites[i].rmps[0], sites[i].wfnc[0], sites[i-1].lmps[0], sites[i-1].wfnc[0]);
    sites[i-1].ropr[0].clear();
    Renormalize (0, sites[i].mpo, sites[i].ropr[0], sites[i].rmps[0], sites[i].rmps[0], sites[i-1].ropr[0]);
    sites[i].save(i);
  }

  util::Normalize(sites[0].wfnc[0]);
  sites[0].lopr[0].resize(Quantum::zero(), qshape, dshape);
  sites[0].lopr[0] = 1.0;
  sites[0].save(0);
}

double prototype::optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M)
{
  boost::function<void(const QSDArray<3>&, QSDArray<3>&)>
  f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, sysdot.lopr[0], sysdot.ropr[0], _1, _2);
  QSDArray<3> diag(sysdot.wfnc[0].q(), sysdot.wfnc[0].qshape());
  ComputeDiagonal(sysdot.mpo, sysdot.lopr[0], sysdot.ropr[0], diag);
  double energy = davidson::diagonalize(f_contract, diag, sysdot.wfnc[0]);

  if(forward) {
    Canonicalize(1, sysdot.wfnc[0], sysdot.lmps[0], M);
    envdot.wfnc[0].clear();
    ComputeGuess(1, sysdot.lmps[0], sysdot.wfnc[0], envdot.rmps[0], envdot.wfnc[0]);
    envdot.lopr[0].clear();
    Renormalize (1, sysdot.mpo,  sysdot.lopr[0], sysdot.lmps[0], sysdot.lmps[0], envdot.lopr[0]);
  }
  else {
    Canonicalize(0, sysdot.wfnc[0], sysdot.rmps[0], M);
    envdot.wfnc[0].clear();
    ComputeGuess(0, sysdot.rmps[0], sysdot.wfnc[0], envdot.lmps[0], envdot.wfnc[0]);
    envdot.ropr[0].clear();
    Renormalize (0, sysdot.mpo,  sysdot.ropr[0], sysdot.rmps[0], sysdot.rmps[0], envdot.ropr[0]);
  }

  return energy;
}

double prototype::optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M)
{
  QSDArray<4> wfnc;
  QSDArray<4> diag;
  boost::function<void(const QSDArray<4>&, QSDArray<4>&)> f_contract;
  if(forward) {
    QSDgemm(NoTrans, NoTrans, 1.0, sysdot.wfnc[0], envdot.rmps[0], 1.0, wfnc);
    f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, envdot.mpo, sysdot.lopr[0], envdot.ropr[0], _1, _2);
    diag.resize(wfnc.q(), wfnc.qshape());
    ComputeDiagonal(sysdot.mpo, envdot.mpo, sysdot.lopr[0], envdot.ropr[0], diag);
  }
  else {
    QSDgemm(NoTrans, NoTrans, 1.0, envdot.lmps[0], sysdot.wfnc[0], 1.0, wfnc);
    f_contract = boost::bind(ComputeSigmaVector, envdot.mpo, sysdot.mpo, envdot.lopr[0], sysdot.ropr[0], _1, _2);
    diag.resize(wfnc.q(), wfnc.qshape());
    ComputeDiagonal(envdot.mpo, sysdot.mpo, envdot.lopr[0], sysdot.ropr[0], diag);
  }

//cout << "====================================================================================================" << endl;
//cout << "debug: optimize_twosite $ wfnc: " << wfnc << endl;
//cout << "====================================================================================================" << endl;
//cout << "debug: optimize_twosite $ diag: " << diag << endl;
//cout << "====================================================================================================" << endl;

  double energy = davidson::diagonalize(f_contract, diag, wfnc);

  if(forward) {
    Canonicalize(1,        wfnc, sysdot.lmps[0], envdot.wfnc[0], M);
    envdot.lopr[0].clear();
    Renormalize (1, sysdot.mpo,  sysdot.lopr[0], sysdot.lmps[0], sysdot.lmps[0], envdot.lopr[0]);
  }
  else {
    Canonicalize(0,        wfnc, sysdot.rmps[0], envdot.wfnc[0], M);
    envdot.ropr[0].clear();
    Renormalize (0, sysdot.mpo,  sysdot.ropr[0], sysdot.rmps[0], sysdot.rmps[0], envdot.ropr[0]);
  }

  return energy;
}

double prototype::dmrg_sweep(MpStorages& sites, DMRG_ALGORITHM algo, int M)
{
  int    L    = sites.size();
  double emin = 1.0e8;
  // fowrad sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  sites[0].load(0);
  for(int i = 0; i < L-1; ++i) {
    sites[i+1].load(i+1);
    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(1, sites[i], sites[i+1], M);
    else                eswp = optimize_twosite(1, sites[i], sites[i+1], M);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    sites[i].save(i);
  }
  // backward sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  for(int i = L-1; i > 0; --i) {
    sites[i-1].load(i-1);
    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(0, sites[i], sites[i-1], M);
    else                eswp = optimize_twosite(0, sites[i], sites[i-1], M);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    sites[i].save(i);
  }
  sites[0].save(0);
  return emin;
}

double prototype::dmrg(const DmrgInput& input, MpStorages& sites)
{
  const int L = input.N_sites;
  const int M = input.N_max_states;

  sites.resize(L, MpSite("state-0", 1, input.prefix));

  switch(input.model) {
    case HEISENBERG:
      Heisenberg::construct_mpo(sites, input.heisenberg);
      initialize(sites, FermiQuantum(0, input.N_spins), M);
      break;
    case HUBBARD:
      Hubbard::construct_mpo(sites, input.hubbard);
      initialize(sites, FermiQuantum(input.N_elecs, input.N_spins), M);
      break;
    default:
      cout << "\tUnknown model type specified" << endl;
      return 0.0;
  }

  double esav = 1.0e8;

  //
  // two-site optimization
  //
  for(int iter = 0; iter < 20; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(sites, TWOSITE, M);
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

  //
  // one-site optimization
  //
  for(int iter = 0; iter < 20; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(sites, ONESITE, M);
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

void prototype::analysis(const DmrgInput& input, MpStorages& sites)
{
  const int L = input.N_sites;

  // fowrad sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  for(int i = 0; i < L-1; ++i) {
    sites[i].load(i);
    cout << "\tanalyzing transfer operator at site [ " << setw(2) << i << " ] " << endl;
//  AnalyzeTransferOperator(sites[i].wfnc[0], sites[i].wfnc[0]);
    AnalyzeTransferOperator(sites[i].lmps[0], sites[i].lmps[0]);
    sites[i].free();
  }
}

