#ifndef __MPSXX_MP_BASE_HPP
#define __MPSXX_MP_BASE_HPP

#include <vector>

namespace mpsxx {

template<class Q, class Matrix>
class MP_base
{
protected:

  std::vector<Q> quanta_; ///< array of physical quantum numbers

  std::vector<Matrix> store_; ///< matrix elements for each physical quantum number

}; // class MP_base

} // namespace mpsxx

#endif // __MPSXX_MP_BASE_HPP
