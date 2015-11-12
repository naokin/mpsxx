#include <iostream>
#include <vector>
#include <cstring>

#include <legacy/QSPARSE/QSTArray.h>

#include <driver/fileio.h>
#include <driver/parsing_integral.h>

#include "BoundaryOpInfo.h"
#include "gen_qcmpo_utils.h"

namespace mpsxx {

/// Function which generates MPO tensors of quantum chemistry Hamiltonian.
/// ENABLE_SWAP_SWEEP_DIR :
/// This produces MPOs which have essentially the same concept of complementary operators.
int gen_qcmpo (
        size_t Norbs,
        size_t Nelec,
        double Ecore,
  const OneIntArray& oneint,
  const TwoIntArray& twoint,
  const std::string& prefix = ".",
        bool ENABLE_SWAP_SWEEP_DIR = false)
{
//std::cout << "DEBUG :: gen_qcmpo :: 01" << std::endl;
  assert(oneint.shape() == btas::shape(Norbs,Norbs));
  assert(twoint.shape() == btas::shape(Norbs,Norbs,Norbs,Norbs));

  // left  block operators info
  BoundaryOpInfo l_ops(0,Norbs,ENABLE_SWAP_SWEEP_DIR);
  // dot   block operators info 
  BoundaryOpInfo s_ops;
  // right block operators info
  BoundaryOpInfo r_ops;

  // l_indxs contains sites on the left  of a certain boundary
  std::vector<int> l_indxs;
  // l_indxs contains sites on the right of a certain boundary
  std::vector<int> r_indxs;
  for(size_t i = Norbs-1; i > 0; --i) r_indxs.push_back(i);

  size_t nnz = 0;

  for(int i = 0; i < Norbs; ++i) {

std::cout << "DEBUG :: gen_qcmpo :: 02-" << i << std::endl;
//  mpo.clear();

    // set boundary operators
    s_ops.reset(i);
    r_ops.reset(i+1,Norbs,ENABLE_SWAP_SWEEP_DIR);

    bool swap_comp_block = (l_ops.comp() != r_ops.comp());

    size_t nnz_local = 0;

    // 'dot with sys' in Block code
    if(l_ops.comp() == BoundaryOpInfo::RIGHT) {
std::cout << "DEBUG :: gen_qcmpo :: 03-A-" << i << std::endl;
      for(auto l = l_ops.begin(); l != l_ops.end(); ++l) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // Take a direct product O_l x O_s
//std::cout << "DEBUG :: gen_qcmpo :: 04-A-" << i << std::endl;
          std::vector<DMRG_OpType> ls_ops = gen_product_operators(l->first,s->first,r_indxs,swap_comp_block);
          // create operator
          for(size_t j = 0; j < ls_ops.size(); ++j) {
//std::cout << "DEBUG :: gen_qcmpo :: 05-A-" << i << std::endl;
            auto r = r_ops.find(ls_ops[j]);
            if(r != r_ops.end()) {
              ++nnz_local;
//            gen_site_operator(mpo,l->first,l->second,s->first,r->first,r->second,oneint,twoint);
            }
          }
//std::cout << "DEBUG :: gen_qcmpo :: 06-A-" << i << std::endl;
        }
      }
    }
    // 'dot with env'
    else {
//std::cout << "DEBUG :: gen_qcmpo :: 03-B-" << i << std::endl;
      for(auto r = r_ops.begin(); r != r_ops.end(); ++r) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // Take a direct product operator O_s*O_r
          std::vector<DMRG_OpType> rs_ops = gen_product_operators(r->first,s->first,l_indxs,swap_comp_block);
          // create operator
          for(size_t j = 0; j < rs_ops.size(); ++j) {
            auto l = l_ops.find(rs_ops[j]);
            if(l != l_ops.end()) {
              ++nnz_local;
//            gen_site_operator(mpo,l->first,l->second,s->first,r->first,r->second,oneint,twoint);
            }
          }
        }
      }
    }

    nnz += nnz_local;

  }

  return 0;
}

} // namespace mpsxx

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
    if(strcmp(argv[iarg],"-h") == 0) { display_help(); return 0; }
  }

  op::QC::setTables();

  int Norbs;
  int Nelec;
  double Ecore;

  OneIntArray oneint; // btas::Tensor<double,2,CblasRowMajor>
  TwoIntArray twoint; // btas::Tensor<double,4,CblasRowMajor>

  std::ifstream ist_fcidump(fname_fcidump.c_str());

  if(!fname_reorder.empty()) {
std::cout << "DEBUG :: 01-A" << std::endl;
    std::ifstream ist_reorder(fname_reorder.c_str());
    std::vector<int> reorder;
    parsing_reorder(ist_reorder,reorder);
    parsing_fcidump(ist_fcidump,Norbs,Nelec,Ecore,oneint,twoint,reorder);
  }
  else {
std::cout << "DEBUG :: 01-B" << std::endl;
    parsing_fcidump(ist_fcidump,Norbs,Nelec,Ecore,oneint,twoint);
  }

std::cout << "DEBUG :: 02" << std::endl;
  mpsxx::gen_qcmpo(Norbs,Nelec,Ecore,oneint,twoint,prefix,ENABLE_SWAP_SWEEP_DIR);

  return 0;
}
