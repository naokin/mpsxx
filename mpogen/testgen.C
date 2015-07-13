#include <iostream>
#include <vector>
#include <cstring>

#include "bit_operator_type.h"
#include "boundary_opinfo.h"
#include "get_product_ops.h"

//std::string get_mpofile(const std::string& prefix, const int& index)
//{
//  std::ostringstream filename;
//  filename << prefix << "/qc_mpo_site_" << index << /* mpigetrank() << */ ".tmp";
//  return filename.str();
//}

int main(int argc, char* argv[])
{
  using namespace mpsxx;

  size_t N = 10;
  bool enable_swap_sweep_dir = false;
  int iprint = 0;
  for(int iarg = 0; iarg < argc; ++iarg) {
    if(strcmp(argv[iarg],"-n") == 0) N = atoi(argv[++iarg]);
    if(strcmp(argv[iarg],"-s") == 0) enable_swap_sweep_dir = true;
    if(strcmp(argv[iarg],"-v") == 0) iprint = 1;
  }

  // left  operator indices
  boundary_opinfo l_ops(0, N);
  // right operator indices
  boundary_opinfo s_ops;
  // right operator indices
  boundary_opinfo r_ops;

  std::vector<size_t> l_indxs;
  std::vector<size_t> r_indxs;
  for(size_t i = N-1; i > 0; --i) r_indxs.push_back(i);

  size_t total_params = 0;
  for(size_t i = 0; i < N; ++i) {
    std::cout << "====================================================================================================" << std::endl;
    // get boundary operators
    s_ops.reset(i); r_ops.reset(i+1, N, enable_swap_sweep_dir);
    bool swap_sweep_dir = (l_ops.direction() != r_ops.direction());
    // 'dot with sys' in Block code
    size_t nnz = 0;
    if(l_ops.direction() == boundary_opinfo::FORWARD) {
      for(auto l = l_ops.begin(); l != l_ops.end(); ++l) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_l*O_s
          std::vector<mpogen::BIT_OPERATOR_TYPE> ls_ops;
          ls_ops = mpogen::get_product_ops(l->first, s->first, r_indxs, swap_sweep_dir);
          // create operator
          for(size_t j = 0; j < ls_ops.size(); ++j) {
            auto r = r_ops.find(ls_ops[j]);
            if(r != r_ops.end()) {
              if(iprint > 0)
                std::cout << mpogen::translate(l->first) << " x " << mpogen::translate(s->first) << " -> " << mpogen::translate(r->first) << std::endl;
              ++nnz;
            }
          }
        }
      }
    }
    // 'dot with env'
    else {
      for(auto r = r_ops.begin(); r != r_ops.end(); ++r) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_s*O_r
          std::vector<mpogen::BIT_OPERATOR_TYPE> rs_ops = mpogen::get_product_ops(r->first, s->first, l_indxs, swap_sweep_dir);
          // create operator
          for(size_t j = 0; j < rs_ops.size(); ++j) {
            auto l = l_ops.find(rs_ops[j]);
            if(l != l_ops.end()) {
              if(iprint > 0)
                std::cout << mpogen::translate(l->first) << " <- " << mpogen::translate(s->first) << " x " << mpogen::translate(r->first) << std::endl;
              ++nnz;
            }
          }
        }
      }
    }
    std::cout << "\tL: " << l_indxs.size() << " ( " << l_ops.size() << " ) " << std::endl;
    std::cout << "\tR: " << r_indxs.size() << " ( " << r_ops.size() << " ) " << std::endl;
    std::cout << "\tH: " << nnz << std::endl;
    l_ops = r_ops;
    l_indxs.push_back(i); r_indxs.pop_back();
    total_params += nnz;
  }
  std::cout << "====================================================================================================" << std::endl;
  std::cout << "\tTotal: " << total_params << std::endl;
}
