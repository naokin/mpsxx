#ifndef __MPSXX_MP_BLAS_HPP
#define __MPSXX_MP_BLAS_HPP

namespace mpsxx {

/// MPO = MPO x MPO :: {O} = {A} x {B}
template<class Q, class MatrixA, class MatrixB, class MatrixO>
void mult (const MPO<Q,MatrixA>& a, const MPO<Q,MatrixB>& b, MPO<Q,MatrixO>& o)
{
  // check physical quantum numbers are the same
  assert(std::equal(a.quanta().begin(),a.quanta().end(),b.quanta().begin()));
  //
  size_t n = a.quanta().size();
  // check/resize v
  if(!o.empty()) {
    assert(std::equal(a.quanta().begin(),a.quanta().end(),o.quanta().begin()));
  }
  else {
    o.resize(a.q()+b.q(),a.quanta());
  }
  // contract physical indices
  for(size_t i = 0; i < n; ++i)
    for(size_t j = 0; j < n; ++j)
      for(size_t k = 0; k < n; ++k)
        dmult(a(i,k),b(k,j),o(i,j)); // direct product
}

/// MPS = MPO x MPS :: |l'> = {O} x |l>
template<class Q, class MatrixO, class MatrixC, class MatrixV>
void mult (const MPO<Q,MatrixO>& o, const MPS<Q,MatrixC>& c, MPS<Q,MatrixV>& v)
{
  // check physical quantum numbers are the same
  assert(std::equal(o.quanta().begin(),o.quanta().end(),c.quanta().begin()));
  //
  size_t n = c.quanta().size();
  // check/resize v
  if(!v.empty()) {
    assert(std::equal(o.quanta().begin(),o.quanta().end(),v.quanta().begin()));
  }
  else {
    v.resize(o.q()+c.q(),c.quanta());
  }
  // contract physical indices
  for(size_t i = 0; i < n; ++i)
    for(size_t j = 0; j < n; ++j)
      dmult(o(i,j),c(j),v(i)); // direct product
}

/// MPS = MPS x MPO :: <l'| = <l| x {O}
template<class Q, class MatrixC, class MatrixO, class MatrixV>
void mult (const MPS<Q,MatrixC>& c, const MPO<Q,MatrixO>& o, MPS<Q,MatrixV>& v)
{
  // check physical quantum numbers are the same
  assert(std::equal(o.quanta().begin(),o.quanta().end(),c.quanta().begin()));
  //
  size_t n = c.quanta().size();
  // check/resize v
  if(!v.empty()) {
    assert(std::equal(o.quanta().begin(),o.quanta().end(),v.quanta().begin()));
  }
  else {
    v.resize(o.q()+c.q(),c.quanta());
  }
  // contract physical indices
  for(size_t i = 0; i < n; ++i)
    for(size_t j = 0; j < n; ++j)
      dmult(c(j),o(j,i),v(i)); // direct product
}

} // namespace mpsxx

#endif // __MPSXX_MP_BLAS_HPP
