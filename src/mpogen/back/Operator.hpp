#ifndef __MPSXX_MPO_OPERATOR_HPP
#define __MPSXX_MPO_OPERATOR_HPP

#include <BaseOperator.hpp>

namespace mpsxx {

struct DerivedOperator : public BaseOperator
{
   DerivedOperator (const BaseOperator& b, op::CATEGORY c = op::I, op::SPINCASE s = op::N, long i = -1, long j = -1)
   :  BaseOperator(c, s, i, j)
   {
   }

private:

   BaseOperator op_rep_;

};

} // namespace mpsxx

mpsxx::DerivedOperator operator* (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
}

#endif // __MPSXX_MPO_OPERATOR_HPP
