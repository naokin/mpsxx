#ifndef __MPSXX_MPO_ADD_OP_COMPONENT_HPP
#define __MPSXX_MPO_ADD_OP_COMPONENT_HPP

namespace mpsxx {

template<typename Tp, class Q>
void AddOpComponent (
            btas::QSTArray<Tp, 4, Q>& mpo,
      const BaseOperator& lOp,
      const BaseOperator& sOp,
      const BaseOperator& rOp,
      const btas::TArray<Tp, 2>& oneint,
      const btas::TArray<Tp, 4>& twoint)
{
}

template<typename Tp>
void AddOpComponent (
            btas::TArray<Tp, 4>& mpo,
      const BaseOperator& lOp,
      const BaseOperator& sOp,
      const BaseOperator& rOp,
      const btas::TArray<Tp, 2>& oneint,
      const btas::TArray<Tp, 4>& twoint)
{
}

} // namespace mpsxx

#endif // __MPSXX_MPO_ADD_OP_COMPONENT_HPP
