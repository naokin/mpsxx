#include <iostream>
#include <iomanip>

#include <vector>
#include <numeric>

#include <btas/QSPARSE/QSTArray.h>

#include "mpsite.h"
#include "fileio.h"

#include "gen_qc_operators.h"
#include "gen_site_operator.h"
#include "boundary_opinfo.h"
#include "get_product_ops.h"
#include "prime_operators.h"

void mpsxx::gen_qc_operators (
  const size_t& N,
  const double& Ecore,
  const btas::TArray<double,2>& oneint,
  const btas::TArray<double,4>& twoint,
  const std::string& opname,
  const std::string& prefix,
  const bool& ENABLE_SWAP_SWEEP)
{
  // check integral size
  assert(oneint.shape() == btas::shape(N,N));
  assert(twoint.shape() == btas::shape(N,N,N,N));

  // MPO tensor
  btas::QSTArray<double,4,fermion> mpo;

  // left  operator indices
  boundary_opinfo l_ops(0,N);
  // right operator indices
  boundary_opinfo s_ops;
  // right operator indices
  boundary_opinfo r_ops;

  std::vector<size_t> l_indxs;
  std::vector<size_t> r_indxs;
  for(size_t i = N-1; i > 0; --i) r_indxs.push_back(i);

  size_t nnz = 0;
  for(int i = 0; i < N; ++i) {
    std::cout << "\t====================================================================================================" << std::endl;
    std::cout << "\t\tSITE [ " << std::setw(3) << i << " ] " << std::endl;
    std::cout << "\t----------------------------------------------------------------------------------------------------" << std::endl;
    // get boundary operators
    s_ops.reset(i); r_ops.reset(i+1,N,ENABLE_SWAP_SWEEP);
    std::cout << "\t\tL-BLOCK: " << std::setw(3) << l_indxs.size() << " sites ( " << std::setw(6) << l_ops.size() << " ops. ) " << std::endl;
    std::cout << "\t\tR-BLOCK: " << std::setw(3) << r_indxs.size() << " sites ( " << std::setw(6) << r_ops.size() << " ops. ) " << std::endl;
    std::cout << "\t----------------------------------------------------------------------------------------------------" << std::endl;

    // resize mpo array
    std::cout << "\t\tresizing site operator array..." << std::flush;

    mpo.clear();

    bool swap_sweep_dir = ENABLE_SWAP_SWEEP && (l_ops.direction() != r_ops.direction());

    btas::Dshapes l_dshape(l_ops.size(),0);
    btas::Dshapes s_dshape(MpSite<fermion>::quanta().size(),1);
    btas::Dshapes r_dshape(r_ops.size(),0);

    if(ENABLE_SWAP_SWEEP) {
      if(swap_sweep_dir)
        mpo.resize(fermion::zero(),
                   btas::make_array(l_ops.get_qshape(),MpSite<fermion>::quanta(),-MpSite<fermion>::quanta(),r_ops.get_qshape()),
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);
      else if(l_ops.direction() == boundary_opinfo::FORWARD)
        mpo.resize(fermion::zero(),
                   btas::make_array(l_ops.get_qshape(),MpSite<fermion>::quanta(),-MpSite<fermion>::quanta(),-r_ops.get_qshape()),
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);
      else
        mpo.resize(fermion::zero(),
                   btas::make_array(-l_ops.get_qshape(),MpSite<fermion>::quanta(),-MpSite<fermion>::quanta(),r_ops.get_qshape()),
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);
    }
    else {
        mpo.resize(fermion::zero(),
                   btas::make_array(l_ops.get_qshape(),MpSite<fermion>::quanta(),-MpSite<fermion>::quanta(),-r_ops.get_qshape()),
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);
    }

    std::cout << "done" << std::endl;
//  std::cout << "DEBUG :: l_ops " << std::endl << l_ops << std::endl;
//  std::cout << "DEBUG :: s_ops " << std::endl << s_ops << std::endl;
//  std::cout << "DEBUG :: r_ops " << std::endl << r_ops << std::endl;
    // 'dot with sys' in Block code
    size_t nnz_local = 0;
    if(l_ops.direction() == boundary_opinfo::FORWARD) {
      std::cout << "\t\tgenerating site operators (fwd)..." << std::flush;
      for(auto l = l_ops.begin(); l != l_ops.end(); ++l) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_l*O_s
          std::vector<mpogen::BIT_OPERATOR_TYPE> ls_ops = mpogen::get_product_ops(l->first,s->first,r_indxs,swap_sweep_dir);
          // create operator
          for(size_t j = 0; j < ls_ops.size(); ++j) {
            auto r = r_ops.find(ls_ops[j]);
            if(r != r_ops.end()) {
              ++nnz_local;
              gen_site_operator(mpo,l->first,l->second,s->first,r->first,r->second,Ecore,oneint,twoint);
            }
          }
        }
      }
    }
    // 'dot with env'
    else {
      std::cout << "\t\tgenerating site operators (bwd)..." << std::flush;
      for(auto r = r_ops.begin(); r != r_ops.end(); ++r) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_s*O_r
          std::vector<mpogen::BIT_OPERATOR_TYPE> rs_ops = mpogen::get_product_ops(r->first,s->first,l_indxs,swap_sweep_dir);
          // create operator
          for(size_t j = 0; j < rs_ops.size(); ++j) {
            auto l = l_ops.find(rs_ops[j]);
            if(l != l_ops.end()) {
              ++nnz_local;
              gen_site_operator(mpo,l->first,l->second,s->first,r->first,r->second,Ecore,oneint,twoint);
            }
          }
        }
      }
    }
    nnz += nnz_local;
    std::cout << "done ( " << nnz_local << " ops. are generated ) " << std::endl;

    std::cout << "\t\tremoving zero contributed terms (fwd)..." << std::flush;

    btas::TVector<btas::Dshapes,4> mpo_dshape = mpo.check_net_dshape();
    btas::TVector<btas::Dshapes,4> nz_indices;
    for(size_t x = 0; x < mpo_dshape[0].size(); ++x) nz_indices[0].push_back(x);
    for(size_t x = 0; x < mpo_dshape[1].size(); ++x) nz_indices[1].push_back(x);
    for(size_t x = 0; x < mpo_dshape[2].size(); ++x) nz_indices[2].push_back(x);
    for(size_t x = 0; x < mpo_dshape[3].size(); ++x) if(mpo_dshape[3][x] > 0) nz_indices[3].push_back(x);

    mpo = mpo.subarray(nz_indices);
    r_ops.clean(mpo_dshape[3]);

    std::cout << "done ( non-zero elements: " << mpo.nnz() << " ) " << std::endl;

    l_ops = r_ops;
    l_indxs.push_back(i); r_indxs.pop_back();

    std::cout << "\t\tsaving site operator array..." << std::flush;
    save(mpo,getfile(opname,prefix,i));
    mpo.clear();
    std::cout << "done" << std::endl;
  }

  btas::Dshapes r_nz_index(1,0);
  for(int i = N-1; i >= 0; --i) {
    std::cout << "\t====================================================================================================" << std::endl;
    std::cout << "\t\tSITE [ " << std::setw(3) << i << " ] " << std::endl;
    std::cout << "\t----------------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "\t\tloading site operator array..." << std::flush;
    load(mpo,getfile(opname,prefix,i));
    std::cout << "done" << std::endl;

    std::cout << "\t\tremoving zero contributed terms (bwd)..." << std::flush;

    btas::TVector<btas::Dshapes,4> mpo_dshape = mpo.check_net_dshape();
    btas::TVector<btas::Dshapes,4> nz_indices;

    for(size_t x = 0; x < mpo_dshape[0].size(); ++x) if(mpo_dshape[0][x] > 0) nz_indices[0].push_back(x);
    for(size_t x = 0; x < mpo_dshape[1].size(); ++x) nz_indices[1].push_back(x);
    for(size_t x = 0; x < mpo_dshape[2].size(); ++x) nz_indices[2].push_back(x);
                                                     nz_indices[3] = r_nz_index;

    mpo = mpo.subarray(nz_indices);
    btas::Dshapes nz_index_save = nz_indices[0];

    mpo_dshape = mpo.check_net_dshape();
    nz_indices[0].clear(); r_nz_index.clear();
    for(size_t x = 0; x < mpo_dshape[0].size(); ++x)
      if(mpo_dshape[0][x] > 0) {
        nz_indices[0].push_back(x);
        r_nz_index.push_back(nz_index_save[x]);
      }

    nz_indices[3].clear();
    for(size_t x = 0; x < mpo_dshape[3].size(); ++x) nz_indices[3].push_back(x);

    mpo = mpo.subarray(nz_indices);

    size_t ldim = std::accumulate(mpo_dshape[0].begin(),mpo_dshape[0].end(),0);
    size_t rdim = std::accumulate(mpo_dshape[3].begin(),mpo_dshape[3].end(),0);
    std::cout << "done" << std::endl;
    std::cout << "\t\tL::" << std::setw(6) << ldim << ", R::" << std::setw(6) << rdim << " ( elements: " << mpo.nnz() << " ) " << std::endl;

    std::cout << "\t\tsaving site operator array..." << std::flush;
    save(mpo,getfile(opname,prefix,i));
    mpo.clear();
    std::cout << "done" << std::endl;
  }

    std::cout << "\t====================================================================================================" << std::endl;
    std::cout << "\t\tTotal number of operator elements: " << nnz << std::endl;
}

