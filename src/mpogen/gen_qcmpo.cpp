#include <iostream>
#include <vector>
#include <cstring>

#include <legacy/QSPARSE/QSTArray.h>

#include <driver/fileio.h>
#include "gen_qc_operators.h"
#include "driver/parsing_integral.h"
#include "compress_qc_mpos.h"

/// Function which generates MPO tensors of quantum chemistry Hamiltonian.
/// ENABLE_SWAP_SWEEP_DIR :
/// This produces MPOs which have essentially the same concept of complementary operators.
int mpsxx::gen_qcmpo (
        size_t Norbs,
        size_t Nelec,
        double Ecore,
  const mpsxx::OneIntArray& oneint,
  const mpsxx::TwoIntArray& twoint,
  const std::string& prefix = ".",
        bool ENABLE_SWAP_SWEEP_DIR = false)
{
  assert(oneint.extent() == btas::shape(Norbs,Norbs));
  assert(twoint.extent() == btas::shape(Norbs,Norbs,Norbs,Norbs));

  // left  block operators info
  BoundaryOpInfo l_ops(0,N,ENABLE_SWAP_SWEEP_DIR);
  // dot   block operators info 
  BoundaryOpInfo s_ops;
  // right block operators info
  BoundaryOpInfo r_ops;

  // l_indxs contains sites on the left  of a certain boundary
  std::vector<size_t> l_indxs;
  // l_indxs contains sites on the right of a certain boundary
  std::vector<size_t> r_indxs;
  for(size_t i = N-1; i > 0; --i) r_indxs.push_back(i);

  size_t nnz = 0;

  // MPO storage
//MPO_vector<double,fermion> mpos(Norbs);
  btas::QSTArray<double,4,fermion> mpo; // i-th MPO

  btas::Qshapes qs = { fermion(0, 0), fermion(1, 1), fermion(1,-1), fermion(2, 0) };

  for(int i = 0; i < N; ++i) {

//  mpo.clear();

    // set boundary operators
    s_ops.reset(i);
    r_ops.reset(i+1,N,ENABLE_SWAP_SWEEP_DIR);

    bool swap_comp_block = (l_ops.comp() != r_ops.comp());

//  std::vector<size_t> l_dshape(l_ops.size(),0);
//  std::vector<size_t> r_dshape(r_ops.size(),0);
    btas::Dshapes l_dshape(l_ops.size(),0);
    btas::Dshapes s_dshape(4,1);
    btas::Dshapes r_dshape(r_ops.size(),0);

    if(ENABLE_SWAP_SWEEP_DIR) {
      if (swap_comp_block)
//      mpo.resize(i, fermion::zero(), l_ops.qshape(), r_ops.qshape());
        mpo.resize(fermion::zero(),
                   btas::make_array(l_ops.qshape(),qs,-qs,r_ops.qshape())
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);

      else if(l_ops.comp() == BoundaryOpInfo::RIGHT)
//      mpo.resize(i, fermion::zero(), l_ops.qshape(),-r_ops.qshape());
        mpo.resize(fermion::zero(),
                   btas::make_array(l_ops.qshape(),qs,-qs,-r_ops.qshape())
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);

      else
//      mpo.resize(i, fermion::zero(),-l_ops.qshape(), r_ops.qshape());
        mpo.resize(fermion::zero(),
                   btas::make_array(-l_ops.qshape(),qs,-qs,r_ops.qshape())
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);
    }
    else {
//      mpo.resize(i, fermion::zero(), l_ops.qshape(),-r_ops.qshape());
        mpo.resize(fermion::zero(),
                   btas::make_array(l_ops.qshape(),qs,-qs,-r_ops.qshape())
                   btas::make_array(l_dshape,s_dshape,s_dshape,r_dshape),false);
    }

    size_t nnz_local = 0;

    // 'dot with sys' in Block code
    if(l_ops.comp() == BoundaryOpInfo::RIGHT) {
      for(auto l = l_ops.begin(); l != l_ops.end(); ++l) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // Take a direct product O_l x O_s
          std::vector<DMRG_OpType> ls_ops = gen_product_operators(l->first,s->first,r_indxs,swap_comp_block);
          // create operator
          for(size_t j = 0; j < ls_ops.size(); ++j) {
            auto r = r_ops.find(ls_ops[j]);
            if(r != r_ops.end()) {
              ++nnz_local;
              gen_site_operator(mpo,l->first,l->second,s->first,r->first,r->second,oneint,twoint);
            }
          }
        }
      }
    }
    // 'dot with env'
    else {
      for(auto r = r_ops.begin(); r != r_ops.end(); ++r) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // Take a direct product operator O_s*O_r
          std::vector<DMRG_OpType> rs_ops = gen_product_operators(r->first,s->first,l_indxs,swap_comp_block);
          // create operator
          for(size_t j = 0; j < rs_ops.size(); ++j) {
            auto l = l_ops.find(rs_ops[j]);
            if(l != l_ops.end()) {
              ++nnz_local;
              gen_site_operator(mpo,l->first,l->second,s->first,r->first,r->second,oneint,twoint);
            }
          }
        }
      }
    }

    nnz += nnz_local;

    // Remove zero contributions

//  std::array<std::vector<size_t>,2> net_shape = check_net_shape(mpo);
//  std::array<std::vector<size_t>,2> net_index;
//  for(size_t k = 0; k < net_shape[0].size(); ++k)                         net_index[0].push_back(k);
//  for(size_t k = 0; k < net_shape[1].size(); ++k) if(net_shape[1][k] > 0) net_index[1].push_back(k);
    std::array<btas::Dshapes,4> net_shape = mpo.check_net_dshape();
    std::array<btas::Dshapes,4> net_index;
    for(size_t k = 0; k < net_shape[0].size(); ++k)                         net_index[0].push_back(k);
    for(size_t k = 0; k < net_shape[1].size(); ++k)                         net_index[1].push_back(k);
    for(size_t k = 0; k < net_shape[2].size(); ++k)                         net_index[2].push_back(k);
    for(size_t k = 0; k < net_shape[3].size(); ++k) if(net_shape[3][k] > 0) net_index[3].push_back(k);

//  mpo = mpo.submatrix(net_index);
//  r_ops.clean(net_shape[1]);
    mpo = mpo.subarray(net_index);
    r_ops.clean(net_shape[3]);

    l_ops = r_ops;
    l_indxs.push_back(i); r_indxs.pop_back();

    // saving i-th MPO
    // ...
    // clear i-th MPO to free memory
    save(mpo,get_mpofile(prefix,i));
    mpo.clear();
  }

//std::vector<size_t> r_net_index(1,0);
  btas::Dshapes r_net_index(1,0);

  for(int i = N-1; i >= 0; --i) {
    // loading i-th MPO
    // ...
    load(mpo,get_mpofile(prefix,i));

//  std::array<std::vector<size_t>,2> net_shape = check_net_shape(mpo);
//  std::array<std::vector<size_t>,2> net_index;
//  net_index[1] = r_net_index;
//  for(size_t k = 0; k < net_shape[0].size(); ++k) if(net_shape[0][k] > 0) net_index[0].push_back(k);
    std::array<Dshapes,4> net_shape = mpo.check_net_dshape();
    std::array<Dshapes,4> net_index;
    for(size_t k = 0; k < net_shape[0].size(); ++k) if(net_shape[0][k] > 0) net_index[0].push_back(k);
    for(size_t k = 0; k < net_shape[0].size(); ++k)                         net_index[1].push_back(k);
    for(size_t k = 0; k < net_shape[0].size(); ++k)                         net_index[2].push_back(k);
                                                                            net_index[3] = r_net_index;

//  mpo = mpo.submatrix(net_index);
    mpo = mpo.subarray(net_index);

    net_shape = check_net_shape(mpo);

    std::vector<size_t> index_save = net_index[0];
    net_index[0].clear();

    r_net_index.clear();

    for(size_t k = 0; k < net_shape[0].size(); ++k)
      if(net_shape[0][k] > 0) {
        net_index[0].push_back(k);
        r_net_index.push_back(index_save[k]);
      }

//  net_index[1].clear();
//  for(size_t k = 0; k < net_shape[1].size(); ++k) net_index[1].push_back(k);
    net_index[3].clear();
    for(size_t k = 0; k < net_shape[3].size(); ++k) net_index[3].push_back(k);

//  mpo = mpo.submatrix(net_index);
    mpo = mpo.subarray(net_index);

    // saving i-th MPO
    // ...
    // clear i-th MPO to free memory
    save(mpo,get_mpofile(prefix,i));
    mpo.claer();
  }

  return 0;
}

int display_help ()
{
  using std::cout;
  using std::endl;
  cout << "Program gen_qcmpo : generate MPO formalism of quantum chemistry Hamiltonian." << endl;
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "\t-h  print this help messages." << endl;
  cout << "\t-f  specify integral file as an input (default FCIDUMP)." << endl;
  cout << "\t-s  directory prefix where the MPOs will be stored." << endl;
  cout << "\t-x  enable swap sweep direction, which is essentially the same as " << endl;
  cout << "\t    complementary operators formalism." <<endl;
  return 0;
}

int main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  using namespace mpsxx;

  bool ENABLE_SWAP_SWEEP_DIR = false;

  std::string fname_fcidump = "FCIDUMP";
  std::string prefix = ".";
  std::string fname_reorder;

  for(int iarg = 1; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-f") == 0) fname_fcidump = argv[++iarg];
    if(strcmp(argv[iarg],"-r") == 0) fname_reorder = argv[++iarg];
    if(strcmp(argv[iarg],"-s") == 0) prefix = argv[++iarg];
    if(strcmp(argv[iarg],"-x") == 0) ENABLE_SWAP_SWEEP_DIR = true;
    if(strcmp(argv[iarg],"-h") == 0) { display_help(); return 0 }
  }

  op::QC::setTables();

  size_t Norbs;
  size_t Nelec;
  double Ecore;

  OneIntArray oneint; // btas::Tensor<double,2,CblasRowMajor>
  TwoIntArray twoint; // btas::Tensor<double,4,CblasRowMajor>

  std::ifstream ist_fcidump(fname_fcidump.c_str());

  if(!fname_reorder.empty()) {
    std::ifstream ist_reorder(fname_reorder.c_str());
    std::vector<size_t> reorder;
    parsing_reorder(ist_reorder,reorder);
    parsing_fcidump(ist_fcidump,Norbs,Nelec,Ecore,oneint,twoint,reorder);
  }
  else {
    parsing_fcidump(ist_fcidump,Norbs,Nelec,Ecore,oneint,twoint);
  }

  gen_qcmpo(Norbs,Nelec,Ecore,oneint,twoint,prefix,ENABLE_SWAP_SWEEP_DIR);

  return 0;
}
