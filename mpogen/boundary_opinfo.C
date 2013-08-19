
#include <iostream>
#include <vector>

#include "boundary_opinfo.h"

void mpsxx::fermionic::boundary_opinfo::reset(size_t index)
{
  using namespace bit_operator;
  m_bn_dir = NOBOUND;
  m_bn_ops.clear();

  size_t nop = 0;
  // I
  m_bn_ops.insert(std::make_pair(IDEN, nop++));
  // H
  m_bn_ops.insert(std::make_pair((HAM | (index << INDEX_SHIFT)), nop++));

  // Cre
  m_bn_ops.insert(std::make_pair((CRE_A | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((CRE_B | (index << INDEX_SHIFT)), nop++));
  // Des
  m_bn_ops.insert(std::make_pair((DES_A | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((DES_B | (index << INDEX_SHIFT)), nop++));

  // CreComp <= CreCreDes
  m_bn_ops.insert(std::make_pair((COMP | CRE_A | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((COMP | CRE_B | (index << INDEX_SHIFT)), nop++));
  // DesComp <= CreDesDes
  m_bn_ops.insert(std::make_pair((COMP | DES_A | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((COMP | DES_B | (index << INDEX_SHIFT)), nop++));

///* DEBUG */ return;

  // CreCre, CreDes, DesCre, and DesDes
  m_bn_ops.insert(std::make_pair((CRE_A_DES_A | (index << INDEX_SHIFT) | index), nop++));
  m_bn_ops.insert(std::make_pair((CRE_B_DES_B | (index << INDEX_SHIFT) | index), nop++));

  m_bn_ops.insert(std::make_pair((CRE_A_CRE_B | (index << INDEX_SHIFT) | index), nop++));
  m_bn_ops.insert(std::make_pair((CRE_A_DES_B | (index << INDEX_SHIFT) | index), nop++));
  m_bn_ops.insert(std::make_pair((DES_A_CRE_B | (index << INDEX_SHIFT) | index), nop++));
  m_bn_ops.insert(std::make_pair((DES_A_DES_B | (index << INDEX_SHIFT) | index), nop++));

  assert(m_bn_ops.size() == nop);
  return;
}

void mpsxx::fermionic::boundary_opinfo::reset(size_t L, size_t N, bool _enable_swap_sweep_dir)
{
  using namespace bit_operator;
  size_t R = N-L;
  size_t N_estm;
  std::vector<size_t> loop_indices;
  if(!_enable_swap_sweep_dir || L < R) { /* loop indices = 0...L-1 */
    m_bn_dir = FORWARD;
    for(size_t i = 0; i < L; ++i) loop_indices.push_back(i);
    // Estimated size of operators
    N_estm = 2+4*N+2*L+4*L*(2*L-1);
  }
  else      { /* loop indices = L...N-1 */
    m_bn_dir = BACKWARD;
    for(size_t i = L; i < N; ++i) loop_indices.push_back(i);
    // Estimated size of operators
    N_estm = 2+4*N+2*R+4*R*(2*R-1);
  }
  m_bn_ops.clear();

  size_t nop = 0;
  // I
  if(_enable_swap_sweep_dir || R > 0)
    m_bn_ops.insert(std::make_pair(IDEN, nop++));
  // H
  if((!_enable_swap_sweep_dir || R > 0) && L > 0)
    m_bn_ops.insert(std::make_pair(HAM , nop++));
  if(L == 0 || R == 0) return;

  BIT_OPERATOR_TYPE comp_type = (!_enable_swap_sweep_dir || L < R) ? ZERO : COMP;
  for(size_t i = 0; i < L; ++i) {
    // Cre
    m_bn_ops.insert(std::make_pair((comp_type | CRE_A | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | CRE_B | (i << INDEX_SHIFT)), nop++));
    // Des
    m_bn_ops.insert(std::make_pair((comp_type | DES_A | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | DES_B | (i << INDEX_SHIFT)), nop++));
  }

  comp_type ^= COMP;
  for(size_t i = L; i < N; ++i) {
    // CreComp
    m_bn_ops.insert(std::make_pair((comp_type | CRE_A | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | CRE_B | (i << INDEX_SHIFT)), nop++));
    // DesComp
    m_bn_ops.insert(std::make_pair((comp_type | DES_A | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | DES_B | (i << INDEX_SHIFT)), nop++));
  }

///* DEBUG */ return;

  // CreCre, CreDes, DesCre, and DesDes
  for(size_t i = 0; i < loop_indices.size(); ++i) {
    size_t ix = loop_indices[i];
    m_bn_ops.insert(std::make_pair((CRE_A_DES_A | (ix << INDEX_SHIFT) | ix), nop++));
    m_bn_ops.insert(std::make_pair((CRE_B_DES_B | (ix << INDEX_SHIFT) | ix), nop++));

    m_bn_ops.insert(std::make_pair((CRE_A_CRE_B | (ix << INDEX_SHIFT) | ix), nop++));
    m_bn_ops.insert(std::make_pair((CRE_A_DES_B | (ix << INDEX_SHIFT) | ix), nop++));
    m_bn_ops.insert(std::make_pair((DES_A_CRE_B | (ix << INDEX_SHIFT) | ix), nop++));
    m_bn_ops.insert(std::make_pair((DES_A_DES_B | (ix << INDEX_SHIFT) | ix), nop++));

    for(size_t j = 0; j < i; ++j) {
      size_t iz = ix;
      size_t jx = loop_indices[j];
      if(iz > jx) std::swap(iz, jx);
      m_bn_ops.insert(std::make_pair((CRE_A_CRE_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((CRE_A_CRE_B | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((CRE_B_CRE_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((CRE_B_CRE_B | (iz << INDEX_SHIFT) | jx), nop++));

      m_bn_ops.insert(std::make_pair((CRE_A_DES_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((CRE_A_DES_B | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((CRE_B_DES_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((CRE_B_DES_B | (iz << INDEX_SHIFT) | jx), nop++));

      m_bn_ops.insert(std::make_pair((DES_A_CRE_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((DES_A_CRE_B | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((DES_B_CRE_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((DES_B_CRE_B | (iz << INDEX_SHIFT) | jx), nop++));

      m_bn_ops.insert(std::make_pair((DES_A_DES_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((DES_A_DES_B | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((DES_B_DES_A | (iz << INDEX_SHIFT) | jx), nop++));
      m_bn_ops.insert(std::make_pair((DES_B_DES_B | (iz << INDEX_SHIFT) | jx), nop++));
    }
  }

  assert(m_bn_ops.size() == nop);
  return;
}

btas::Qshapes<mpsxx::fermionic::Quantum> mpsxx::fermionic::boundary_opinfo::get_qshape() const
{
  size_t nop = m_bn_ops.size();
  btas::Qshapes<Quantum> q(nop);
  for(const_iterator it = m_bn_ops.begin(); it != m_bn_ops.end(); ++it)
    q.at(it->second) = get_quantum(it->first);
  return std::move(q);
}

void mpsxx::fermionic::boundary_opinfo::clean(const btas::Dshapes& _dn_shape)
{
  std::map<size_t, BIT_OPERATOR_TYPE> _bn_map;
  for(auto it = m_bn_ops.begin(); it != m_bn_ops.end();) {
    if(_dn_shape[it->second] > 0) {
      _bn_map.insert(std::make_pair(it->second, it->first));
      ++it;
    }
    else {
      it = m_bn_ops.erase(it);
    }
  }
  size_t nnz = 0;
  for(auto it = _bn_map.begin(); it != _bn_map.end(); ++it, ++nnz)
    m_bn_ops.find(it->second)->second = nnz;
}

std::ostream& operator<< (std::ostream& ost, const mpsxx::fermionic::boundary_opinfo& info)
{
  for(auto it = info.begin(); it != info.end(); ++it)
    ost << "\t[ " << it->second << " ] : " << mpsxx::fermionic::translate(it->first) << std::endl;
  return ost;
}

