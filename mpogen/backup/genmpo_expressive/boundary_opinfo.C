
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
  m_bn_ops.insert(std::make_pair((HAM | index), nop++));

  // Cre
  m_bn_ops.insert(std::make_pair((CRE_A_1 | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((CRE_B_1 | (index << INDEX_SHIFT)), nop++));
  // Des
  m_bn_ops.insert(std::make_pair((DES_A_1 | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((DES_B_1 | (index << INDEX_SHIFT)), nop++));

  // CreComp <= CreCreDes
  m_bn_ops.insert(std::make_pair((COMP | CRE_A_1 | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((COMP | CRE_B_1 | (index << INDEX_SHIFT)), nop++));
  // DesComp <= CreDesDes
  m_bn_ops.insert(std::make_pair((COMP | DES_A_1 | (index << INDEX_SHIFT)), nop++));
  m_bn_ops.insert(std::make_pair((COMP | DES_B_1 | (index << INDEX_SHIFT)), nop++));

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

void mpsxx::fermionic::boundary_opinfo::reset(size_t L, size_t N)
{
  using namespace bit_operator;
  size_t R = N-L;
  size_t N_estm;
  std::vector<size_t> loop_indices;
  if(L < R) { /* loop indices = 0...L-1 */
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
  m_bn_ops.insert(std::make_pair(IDEN, nop++));
  if(L == 0 || R == 0) return;
  // H
  m_bn_ops.insert(std::make_pair(HAM , nop++));

  BIT_OPERATOR_TYPE comp_type = (L < R) ? ZERO : COMP;
  for(size_t i = 0; i < L; ++i) {
    // Cre
    m_bn_ops.insert(std::make_pair((comp_type | CRE_A_1 | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | CRE_B_1 | (i << INDEX_SHIFT)), nop++));
    // Des
    m_bn_ops.insert(std::make_pair((comp_type | DES_A_1 | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | DES_B_1 | (i << INDEX_SHIFT)), nop++));
  }

  comp_type ^= COMP;
  for(size_t i = L; i < N; ++i) {
    // CreComp
    m_bn_ops.insert(std::make_pair((comp_type | CRE_A_1 | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | CRE_B_1 | (i << INDEX_SHIFT)), nop++));
    // DesComp
    m_bn_ops.insert(std::make_pair((comp_type | DES_A_1 | (i << INDEX_SHIFT)), nop++));
    m_bn_ops.insert(std::make_pair((comp_type | DES_B_1 | (i << INDEX_SHIFT)), nop++));
  }

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

std::ostream& operator<< (std::ostream& ost, const mpsxx::fermionic::boundary_opinfo& info)
{
  for(auto it = info.begin(); it != info.end(); ++it)
    ost << mpsxx::fermionic::translate(it->first) << std::endl;
  return ost;
}

