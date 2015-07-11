#ifndef __MPSXX_GENERIC_RENORMALIZE_HPP
#define __MPSXX_GENERIC_RENORMALIZE_HPP

namespace mpsxx {

template<class _MPS, class _MPO, class _Tensor>
void renormalize (
      const MPS_SWEEP_DIR& dir,
      const _MPO& mpo,
      const _Ops& xop,
      const _MPS& bra,
      const _MPS& ket,
            _Ops& yop)
{
   if(dir == ForwardSweep)
   {
      yop("bj,wj,kj") += bra("bi,pi,bj")*mpo("wi,pi,qi,wj")*xop("bi,wi,ki")*ket("ki,qi,kj");
   }
   else
   {
      yop("bi,wi,ki") += bra("bi,pi,bj")*mpo("wi,pi,qi,wj")*xop("bj,wj,kj")*ket("ki,qi,kj");
   }
}

template<class _MPS, class _MPO, class _Tensor>
void renormalize (
      const MPS_SWEEP_DIR& dir,
      const _Ops& xop,
      const _MPS& bra,
      const _MPS& ket,
            _Ops& yop)
{
   if(dir == ForwardSweep)
   {
      yop("bj,kj") += bra("bi,pi,bj")*xop("bi,ki")*ket("ki,qi,kj");
   }
   else
   {
      yop("bi,ki") += bra("bi,pi,bj")*xop("bj,kj")*ket("ki,qi,kj");
   }
}

} // namespace mpsxx

#endif // __MPSXX_GENERIC_RENORMALIZE_HPP
