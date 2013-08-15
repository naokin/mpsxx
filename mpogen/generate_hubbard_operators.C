#include <iostream>
#include <iomanip>

#include <vector>

#include <btas/QSPARSE/QSDArray.h>

#include <MpSite.h>
#include <driver/fileio.h>

#include "generate_hubbard_operators.h"
#include "generate_site_operator.h"
#include "prime_operators.h"

void mpsxx::fermionic::generate_hubbard_operators
(mpsxx::MpOperators<mpsxx::fermionic::Quantum>& mpos, const double& t, const double& u, const std::string& prefix)
{
  using std::cout;
  using std::endl;
  using std::setw;
  using std::fixed;

  size_t N = mpos.size();

  cout << "\t====================================================================================================" << endl;
  cout << "\t\tCONSTRUCT FERMIONIC-HUBBARD MATRIX PRODUCT OPERATORS (MPOs) "                                       << endl;
  cout.precision(4);
  cout << "\t\t\t+ coupling coefficient t  : " << setw(8) << fixed << t  << endl; 
  cout << "\t\t\t+ coupling coefficient U  : " << setw(8) << fixed << u  << endl;
  cout << "\t====================================================================================================" << endl;

  btas::Qshapes<Quantum> qz; // 0 quantum number
  qz.push_back(Quantum( 0,  0));

  btas::Qshapes<Quantum> qi; // quantum index comes in
  qi.push_back(Quantum( 0,  0)); // 0
  qi.push_back(Quantum( 1,  1)); // a-
  qi.push_back(Quantum(-1, -1)); // a+
  qi.push_back(Quantum( 1, -1)); // b-
  qi.push_back(Quantum(-1,  1)); // b+
  qi.push_back(Quantum( 0,  0)); // I

  btas::Qshapes<Quantum> qo; // quantum index comes out
  qo.push_back(Quantum( 0,  0)); // I
  qo.push_back(Quantum(-1, -1)); // a+
  qo.push_back(Quantum( 1,  1)); // a-
  qo.push_back(Quantum(-1,  1)); // b+
  qo.push_back(Quantum( 1, -1)); // b-
  qo.push_back(Quantum( 0,  0)); // 0

  // resize & set to 0
  mpos[ 0 ].resize(Quantum::zero(), make_array(qz, MpSite<Quantum>::quanta(),-MpSite<Quantum>::quanta(), qo));
  for(int i = 1; i < N-1; ++i)
    mpos[i].resize(Quantum::zero(), make_array(qi, MpSite<Quantum>::quanta(),-MpSite<Quantum>::quanta(), qo));
  mpos[N-1].resize(Quantum::zero(), make_array(qi, MpSite<Quantum>::quanta(),-MpSite<Quantum>::quanta(), qz));

  // set block elements
  btas::DArray<4> data_Ip(1, 1, 1, 1); data_Ip = 1.0;
  btas::DArray<4> data_Im(1, 1, 1, 1); data_Im =-1.0;
  btas::DArray<4> data_tp(1, 1, 1, 1); data_tp = t;
  btas::DArray<4> data_tm(1, 1, 1, 1); data_tm =-t;
  btas::DArray<4> data_Un(1, 1, 1, 1); data_Un = u;
  // insert blocks
  mpos[ 0 ].insert(btas::shape(0, 3, 3, 0), data_Un); //  U
  mpos[ 0 ].insert(btas::shape(0, 1, 0, 1), data_tp); //  t a+
  mpos[ 0 ].insert(btas::shape(0, 3, 2, 1), data_tm); //  t a+ [x(-1) due to P(a,b)]
  mpos[ 0 ].insert(btas::shape(0, 0, 1, 2), data_tm); // -t a-
  mpos[ 0 ].insert(btas::shape(0, 2, 3, 2), data_tp); // -t a- [x(-1) due to P(a,b)]
  mpos[ 0 ].insert(btas::shape(0, 2, 0, 3), data_tp); //  t b+
  mpos[ 0 ].insert(btas::shape(0, 3, 1, 3), data_tp); //  t b+
  mpos[ 0 ].insert(btas::shape(0, 0, 2, 4), data_tm); // -t b-
  mpos[ 0 ].insert(btas::shape(0, 1, 3, 4), data_tm); // -t b-
  mpos[ 0 ].insert(btas::shape(0, 0, 0, 5), data_Ip); //  I
  mpos[ 0 ].insert(btas::shape(0, 1, 1, 5), data_Ip); //  I
  mpos[ 0 ].insert(btas::shape(0, 2, 2, 5), data_Ip); //  I
  mpos[ 0 ].insert(btas::shape(0, 3, 3, 5), data_Ip); //  I
  for(int i = 1; i < N-1; ++i) {
    mpos[i].insert(btas::shape(0, 0, 0, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(0, 1, 1, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(0, 2, 2, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(0, 3, 3, 0), data_Ip); //  I
    mpos[i].insert(btas::shape(1, 0, 1, 0), data_Ip); //  a-
    mpos[i].insert(btas::shape(1, 2, 3, 0), data_Im); //  a- [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(2, 1, 0, 0), data_Ip); //  a+
    mpos[i].insert(btas::shape(2, 3, 2, 0), data_Im); //  a+ [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(3, 0, 2, 0), data_Ip); //  b-
    mpos[i].insert(btas::shape(3, 1, 3, 0), data_Ip); //  b-
    mpos[i].insert(btas::shape(4, 2, 0, 0), data_Ip); //  b+
    mpos[i].insert(btas::shape(4, 3, 1, 0), data_Ip); //  b+
    mpos[i].insert(btas::shape(5, 3, 3, 0), data_Un); //  U
    mpos[i].insert(btas::shape(5, 1, 0, 1), data_tp); //  t a+
    mpos[i].insert(btas::shape(5, 3, 2, 1), data_tm); //  t a+ [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(5, 0, 1, 2), data_tm); // -t a-
    mpos[i].insert(btas::shape(5, 2, 3, 2), data_tp); // -t a- [x(-1) due to P(a,b)]
    mpos[i].insert(btas::shape(5, 2, 0, 3), data_tp); //  t b+
    mpos[i].insert(btas::shape(5, 3, 1, 3), data_tp); //  t b+
    mpos[i].insert(btas::shape(5, 0, 2, 4), data_tm); // -t b-
    mpos[i].insert(btas::shape(5, 1, 3, 4), data_tm); // -t b-
    mpos[i].insert(btas::shape(5, 0, 0, 5), data_Ip); //  I
    mpos[i].insert(btas::shape(5, 1, 1, 5), data_Ip); //  I
    mpos[i].insert(btas::shape(5, 2, 2, 5), data_Ip); //  I
    mpos[i].insert(btas::shape(5, 3, 3, 5), data_Ip); //  I
  }
  mpos[N-1].insert(btas::shape(0, 0, 0, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(0, 1, 1, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(0, 2, 2, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(0, 3, 3, 0), data_Ip); //  I
  mpos[N-1].insert(btas::shape(1, 0, 1, 0), data_Ip); //  a-
  mpos[N-1].insert(btas::shape(1, 2, 3, 0), data_Im); //  a- [x(-1) due to P(a,b)]
  mpos[N-1].insert(btas::shape(2, 1, 0, 0), data_Ip); //  a+
  mpos[N-1].insert(btas::shape(2, 3, 2, 0), data_Im); //  a+ [x(-1) due to P(a,b)]
  mpos[N-1].insert(btas::shape(3, 0, 2, 0), data_Ip); //  b-
  mpos[N-1].insert(btas::shape(3, 1, 3, 0), data_Ip); //  b-
  mpos[N-1].insert(btas::shape(4, 2, 0, 0), data_Ip); //  b+
  mpos[N-1].insert(btas::shape(4, 3, 1, 0), data_Ip); //  b+
  mpos[N-1].insert(btas::shape(5, 3, 3, 0), data_Un); //  U

  // taking operator parity P(O(l),<n|)
  std::vector<int> indx1(1, 0); // left mpo index [0]
  std::vector<int> indx2(1, 1); // phys bra index [1]
  for(int i = 0; i < N; ++i) {
    mpos[i].parity(indx1, indx2);
//  save(mpos[i], get_mpofile(prefix, FERMION_HUBBARD, i));
    save(mpos[i], get_mpofile(prefix, i));
    mpos[i].clear();
  }
}

/*

{
  // left  operator indices
  boundary_opinfo l_ops(0, N);
  // right operator indices
  boundary_opinfo s_ops;
  // right operator indices
  boundary_opinfo r_ops;

  std::vector<size_t> l_indxs;
  std::vector<size_t> r_indxs;
  for(size_t i = N-1; i > 0; --i) r_indxs.push_back(i);

  size_t nnz = 0;
  for(size_t i = 0; i < N; ++i) {
    std::cout << "\t\t====================================================================================================" << std::endl;
    std::cout << "\t\t\tSITE [ " << std::setw(3) << i << " ] " << std::endl;
    std::cout << "\t\t----------------------------------------------------------------------------------------------------" << std::endl;
    // get boundary operators
    s_ops.reset(i); r_ops.reset(i+1, N, enable_swap_sweep_dir);
    std::cout << "\t\t\tL-BLOCK: " << std::setw(3) << l_indxs.size() << " sites ( " << std::setw(6) << l_ops.size() << " ops. ) " << std::endl;
    std::cout << "\t\t\tR-BLOCK: " << std::setw(3) << r_indxs.size() << " sites ( " << std::setw(6) << r_ops.size() << " ops. ) " << std::endl;
    std::cout << "\t\t----------------------------------------------------------------------------------------------------" << std::endl;
    // resize mpo array
    mpos[i].clear();
    bool swap_sweep_dir = (l_ops.direction() != r_ops.direction());
    std::cout << "\t\t\tresizing site operator array..." << std::flush;
    if(enable_swap_sweep_dir) {
      if     (swap_sweep_dir)
        mpos[i].resize(Quantum::zero(), make_array( l_ops.get_qshape(), MpSite<Quantum>::quanta(),-MpSite<Quantum>::quanta(), r_ops.get_qshape()));
      else if(l_ops.direction() == boundary_opinfo::FORWARD)
        mpos[i].resize(Quantum::zero(), make_array( l_ops.get_qshape(), MpSite<Quantum>::quanta(),-MpSite<Quantum>::quanta(),-r_ops.get_qshape()));
      else
        mpos[i].resize(Quantum::zero(), make_array(-l_ops.get_qshape(), MpSite<Quantum>::quanta(),-MpSite<Quantum>::quanta(), r_ops.get_qshape()));
    }
    else {
        mpos[i].resize(Quantum::zero(), make_array( l_ops.get_qshape(), MpSite<Quantum>::quanta(),-MpSite<Quantum>::quanta(),-r_ops.get_qshape()));
    }
    std::cout << "done" << std::endl;
    // 'dot with sys' in Block code
    size_t nnz_local = 0;
    if(l_ops.direction() == boundary_opinfo::FORWARD) {
      std::cout << "\t\t\tgenerating site operators (fwd)..." << std::flush;
      for(auto l = l_ops.begin(); l != l_ops.end(); ++l) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_l*O_s
          std::vector<BIT_OPERATOR_TYPE> ls_ops = get_product_ops(l->first, s->first, r_indxs, swap_sweep_dir);
          // create operator
          for(size_t j = 0; j < ls_ops.size(); ++j) {
            auto r = r_ops.find(ls_ops[j]);
            if(r != r_ops.end()) {
              ++nnz_local;
//            if(nnz_local %  100 == 0) std::cout << nnz_local / 100 << "..." << std::flush;
//            if(nnz_local % 1000 == 0) std::cout << std::endl << "\t\t\t                               ..." << std::flush;
//            std::cout << "DEBUG[generate_qc_operators]: " << translate(l->first) << " x " << translate(s->first) << " -> " << translate(r->first) << std::endl;
              generate_site_operator(mpos[i], l->first, l->second, s->first, r->first, r->second, oneint, twoint);
            }
          }
        }
      }
    }
    // 'dot with env'
    else {
      std::cout << "\t\t\tgenerating site operators (bwd)..." << std::flush;
      for(auto r = r_ops.begin(); r != r_ops.end(); ++r) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_s*O_r
          std::vector<BIT_OPERATOR_TYPE> rs_ops = get_product_ops(r->first, s->first, l_indxs, swap_sweep_dir);
          // create operator
          for(size_t j = 0; j < rs_ops.size(); ++j) {
            auto l = l_ops.find(rs_ops[j]);
            if(l != l_ops.end()) {
              ++nnz_local;
//            if(nnz_local %  100 == 0) std::cout << nnz_local / 100 << "..." << std::flush;
//            if(nnz_local % 1000 == 0) std::cout << std::endl << "\t\t\t                               ..." << std::flush;
//            std::cout << "DEBUG[generate_qc_operators]: " << translate(l->first) << " <- " << translate(s->first) << " x " << translate(r->first) << std::endl;
              generate_site_operator(mpos[i], l->first, l->second, s->first, r->first, r->second, oneint, twoint);
            }
          }
        }
      }
    }
    nnz += nnz_local;
    std::cout << "done ( " << nnz_local << " ops. are generated ) " << std::endl;
    l_ops = r_ops;
    l_indxs.push_back(i); r_indxs.pop_back();
    std::cout << "\t\t\tsaving site operator array..." << std::flush;
    save(mpos[i], get_mpofile(prefix, MOLECULAR, i));
    mpos[i].clear();
    std::cout << "done" << std::endl;
  }
    std::cout << "\t\t====================================================================================================" << std::endl;
    std::cout << "\t\t\tTotal number of operator elements: " << nnz << std::endl;
}

*/
