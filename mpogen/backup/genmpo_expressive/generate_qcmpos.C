#include <vector>
#include <cstring>

#include <btas/QSPARSE/QSDArray.h>

#include "boundary_opinfo.h"
#include "get_product_ops.h"

using mpsxx::fermionic::Quantum;

std::string get_mpofile(const std::string& prefix, const int& index)
{
  std::stringstream filename;
  filename << prefix << "/qc_mpo_site_" << index << /* mpigetrank() << */ ".tmp";
  return filename.str();
}

void mpsxx::fermionic::generate_mpoperators
(mpsxx::MpOperators& mpos, const btas::DArray<2>& oneint, const btas::DArray<4>& twoint)
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

  for(size_t i = 0; i < N; ++i) {
    // get boundary operators
    s_ops.reset(i); r_ops.reset(i+1,N);
    // resize mpo array
    mpos[i].clear();
    mpos[i].resize(Fermion::zero(), make_array(l_ops.get_qshape(), fock(),-fock(),-r_ops.get_qshape()));
    // 'dot with sys' in Block code
    if(l_ops.direction == boundary_opinfo::FORWARD) {
      for(auto l = l_ops.begin(); l != l_ops.end(); ++l) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_l*O_s
          if(r_ops.direction == boundary_opinfo::FORWARD)
            std::vector<BIT_OPERATOR_TYPE> ls_ops = get_product_ops(l->first, s->first, r_indxs, false);
          else // swapping sweep direction
            std::vector<BIT_OPERATOR_TYPE> ls_ops = get_product_ops(l->first, s->first, r_indxs, true);
          // create operator
          for(size_t j = 0; j < ls_ops.size(); ++j) {
            auto r = r_ops.find(ls_ops[j]);
            if(r != r_ops.end())
              set_site_operator(mpos[i], l->first, l->second, s->first, r->first, r->second, oneint, twoint);
          }
        }
      }
    }
    // 'dot with env'
    else {
      for(auto r = r_ops.begin(); r != r_ops.end(); ++r) {
        for(auto s = s_ops.begin(); s != s_ops.end(); ++s) {
          // prod. operator O_s*O_r
          std::vector<BIT_OPERATOR_TYPE> rs_ops = get_product_ops(r->first, s->first, l_indxs, false);
          // create operator
          for(size_t j = 0; j < rs_ops.size(); ++j) {
            auto l = l_ops.find(rs_ops[j]);
            if(l != l_ops.end())
              set_site_operator(mpos[i], l->first, l->second, s->first, r->first, r->second, oneint, twoint);
          }
        }
      }
    }
    l_ops = r_ops;
    l_indxs.push_back(i); r_indxs.pop_back();
    btas::save(mpos[i], get_mpofile(prefix, i));
    mpos[i].clear();
  }
}

void mpsxx:fermionic::generate_site_operator
(btas::QSDArray<4, mpsxx::fermionic::Quantum>& op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& l_op, const size_t& l_index,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& s_op,
 const mpsxx::fermionic::BIT_OPERATOR_TYPE& r_op, const size_t& r_index,
 const btas::DArray<2>& oneint, const btas::DArray<4>& twoint)
{
  int sx;    // site  index
  int ix,jx; // left  block indices
  int kx,lx; // right block indices

  switch(s_op & TYPE) {
    case IDEN:
      if(l_op == r_op || l_op == (r_op ^ COMP) || l_op == (r_op ^ CONJ_S_1 ^ COMP)) {
// NOTE: more precise checking
//    if(l_op == r_op || ((l_op & IDEN) && l_op == (r_op ^ COMP))
//                    || ((l_op & NORMAL & TYPE) == SINGLE_1 && l_op == (r_op ^ CONJ_S_1 ^ COMP))) {
        prime_op_generator<_identity>(op, l_index, r_index, 1.0);
      }
      else if((l_op & TYPE) == DOUBLE && (r_op & TYPE) == DOUBLE) {
        ix = (l_op & INDEX & FIRST) >> INDEX_SHIFT;
        jx = (l_op & INDEX & SECOND);
        kx = (r_op & INDEX & FIRST) >> INDEX_SHIFT;
        lx = (r_op & INDEX & SECOND);
        // FIXME: compute parity somehow
        double parity = 1.0;
        prime_op_generator<_identity>(op, l_index, r_index, parity*twoint(ix,jx,kx,lx));
      }
      break;

    case HAM:
      sx = (s_op & INDEX & SECOND);
      prime_op_generator<_cre_a_des_a>(op, l_index, r_index, oneint(sx,sx));
      prime_op_generator<_cre_b_des_b>(op, l_index, r_index, oneint(sx,sx));
      prime_op_generator<_cre_a_des_a_cre_b_des_b(op, l_index, r_index, twoint(sx,sx,sx,sx)); 
      break;

    case SINGLE_1:
      switch(s_op & MASK) {
        // FIXME: compute parity somehow
        double parity = 1.0;
        case CRE_A_1:
          prime_op_generator<_cre_a>(op, l_index, r_index, parity);
          break;
        case CRE_B_1:
          prime_op_generator<_cre_b>(op, l_index, r_index, parity);
          break;
        case DES_A_1:
          prime_op_generator<_des_a>(op, l_index, r_index, parity);
          break;
        case DES_B_1:
          prime_op_generator<_des_b>(op, l_index, r_index, parity);
          break;
        default:
      }
      break;

    case SINGLE_1 | COMP:
      switch(s_op & MASK) {
        // site  index
        sx = (s_op & INDEX & FIRST);
        // left  index
        if((l_op & TYPE) == SINGLE_1)
          ix = (l_op & INDEX & FIRST);
        else
          ix = -1;
        // right index
        if((r_op & TYPE) == SINGLE_1)
          kx = (r_op & INDEX & FIRST);
        else
          kx = -1;
        // FIXME: compute parity somehow
        double parity = 1.0;
        case CRE_A_1:
          if     (ix > 0) {
            prime_op_generator<_cre_a>(op, l_index, r_index, parity*oneint(ix,sx));
            prime_op_generator<_cre_a_cre_b_des_b>(op, l_index, r_index, parity*twoint(ix,sx,sx,sx));
          }
          else if(kx > 0) {
            prime_op_generator<_cre_a>(op, l_index, r_index, parity*oneint(sx,kx));
            prime_op_generator<_cre_a_cre_b_des_b>(op, l_index, r_index, parity*twoint(sx,sx,sx,kx));
          }
          break;
        case CRE_B_1:
          if     (ix > 0) {
            prime_op_generator<_cre_b>(op, l_index, r_index, parity*oneint(ix,sx));
            prime_op_generator<_cre_b_cre_a_des_a>(op, l_index, r_index, parity*twoint(ix,sx,sx,sx));
          }
          else if(kx > 0) {
            prime_op_generator<_cre_b>(op, l_index, r_index, parity*oneint(sx,kx));
            prime_op_generator<_cre_b_cre_a_des_a>(op, l_index, r_index, parity*twoint(sx,sx,sx,kx));
          }
          break;
        case DES_A_1:
          if     (ix > 0) {
            prime_op_generator<_des_a>(op, l_index, r_index, parity*oneint(ix,sx));
            prime_op_generator<_cre_b_des_b_des_a>(op, l_index, r_index, parity*twoint(ix,sx,sx,sx));
          }
          else if(kx > 0) {
            prime_op_generator<_des_a>(op, l_index, r_index, parity*oneint(sx,kx));
            prime_op_generator<_cre_b_des_b_des_a>(op, l_index, r_index, parity*twoint(sx,sx,sx,kx));
          }
          break;
        case DES_B_1:
          if     (ix > 0) {
            prime_op_generator<_des_b>(op, l_index, r_index, parity*oneint(ix,sx));
            prime_op_generator<_cre_a_des_a_des_b>(op, l_index, r_index, parity*twoint(ix,sx,sx,sx));
          }
          else if(kx > 0) {
            prime_op_generator<_des_b>(op, l_index, r_index, parity*oneint(sx,kx));
            prime_op_generator<_cre_a_des_a_des_b>(op, l_index, r_index, parity*twoint(sx,sx,sx,kx));
          }
          break;
        default:
      }
      break;

    case DOUBLE:
      break;

    default:
  }
}


