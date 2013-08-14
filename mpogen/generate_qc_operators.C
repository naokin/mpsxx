#include <iostream>
#include <iomanip>

#include <vector>

#include <btas/QSPARSE/QSDArray.h>

#include <MpSite.h>
#include <driver/fileio.h>

#include "generate_qc_operators.h"
#include "generate_site_operator.h"
#include "boundary_opinfo.h"
#include "get_product_ops.h"
#include "prime_operators.h"

void mpsxx::fermionic::generate_qc_operators
(mpsxx::MpOperators<mpsxx::fermionic::Quantum>& mpos, const btas::DArray<2>& oneint, const btas::DArray<4>& twoint, bool enable_swap_sweep_dir, const std::string& prefix)
{
  size_t N = mpos.size();

  // check integral size
  assert(oneint.shape() == btas::shape(N,N));
  assert(twoint.shape() == btas::shape(N,N,N,N));

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

