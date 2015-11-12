#ifndef __MPSXX_MPOGEN_ADD_OPERATORS_H
#define __MPSXX_MPOGEN_ADD_OPERATORS_H

#include <symmetry/fermion.h>
#include "Category.h"
//#include "MPO.h"
#include <legacy/QSPARSE/QSTArray.h>

namespace mpsxx {

template<op::QC::CATEGORY OpCat> struct OpComponent;

template<>
struct OpComponent<op::QC::I>
{
  static constexpr unsigned int size = 4u;
  static constexpr unsigned int bra[size] = { 0, 1, 2, 3 };
  static constexpr unsigned int ket[size] = { 0, 1, 2, 3 };
  static const double value[size]; //  1, 1, 1, 1
};

template<>
struct OpComponent<op::QC::CREA>
{
  static constexpr unsigned int size = 2u;
  static constexpr unsigned int bra[size] = { 1, 3 };
  static constexpr unsigned int ket[size] = { 0, 2 };
  static const double value[size]; //  1,-1
};

template<>
struct OpComponent<op::QC::CREB>
{
  static constexpr unsigned int size = 2u;
  static constexpr unsigned int bra[size] = { 2, 3 };
  static constexpr unsigned int ket[size] = { 0, 1 };
  static const double value[size]; //  1, 1
};

template<>
struct OpComponent<op::QC::DESA>
{
  static constexpr unsigned int size = 2u;
  static constexpr unsigned int bra[size] = { 0, 2 };
  static constexpr unsigned int ket[size] = { 1, 3 };
  static const double value[size]; //  1,-1
};

template<>
struct OpComponent<op::QC::DESB>
{
  static constexpr unsigned int size = 2u;
  static constexpr unsigned int bra[size] = { 0, 1 };
  static constexpr unsigned int ket[size] = { 2, 3 };
  static const double value[size]; //  1, 1
};

template<>
struct OpComponent<op::QC::CREA_CREB>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 3 };
  static constexpr unsigned int ket[size] = { 0 };
  static const double value[size]; // -1
};

template<>
struct OpComponent<op::QC::CREB_CREA>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 3 };
  static constexpr unsigned int ket[size] = { 0 };
  static const double value[size]; //  1
};

template<>
struct OpComponent<op::QC::DESA_DESB>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 0 };
  static constexpr unsigned int ket[size] = { 3 };
  static const double value[size]; //  1
};

template<>
struct OpComponent<op::QC::DESB_DESA>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 0 };
  static constexpr unsigned int ket[size] = { 3 };
  static const double value[size]; // -1
};

template<>
struct OpComponent<op::QC::CREA_DESA>
{
  static constexpr unsigned int size = 2u;
  static constexpr unsigned int bra[size] = { 1, 3 };
  static constexpr unsigned int ket[size] = { 1, 3 };
  static const double value[size]; //  1, 1
};

template<>
struct OpComponent<op::QC::CREB_DESB>
{
  static constexpr unsigned int size = 2u;
  static constexpr unsigned int bra[size] = { 2, 3 };
  static constexpr unsigned int ket[size] = { 2, 3 };
  static const double value[size]; //  1, 1
};

template<>
struct OpComponent<op::QC::CREA_DESB>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 1 };
  static constexpr unsigned int ket[size] = { 2 };
  static const double value[size]; //  1
};

template<>
struct OpComponent<op::QC::CREB_DESA>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 2 };
  static constexpr unsigned int ket[size] = { 1 };
  static const double value[size]; //  1
};

template<>
struct OpComponent<op::QC::CREA_CREB_DESB>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 3 };
  static constexpr unsigned int ket[size] = { 2 };
  static const double value[size]; // -1
};

template<>
struct OpComponent<op::QC::CREB_CREA_DESA>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 3 };
  static constexpr unsigned int ket[size] = { 1 };
  static const double value[size]; //  1
};

template<>
struct OpComponent<op::QC::CREB_DESB_DESA>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 2 };
  static constexpr unsigned int ket[size] = { 3 };
  static const double value[size]; // -1
};

template<>
struct OpComponent<op::QC::CREA_DESA_DESB>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 1 };
  static constexpr unsigned int ket[size] = { 3 };
  static const double value[size]; //  1
};

template<>
struct OpComponent<op::QC::CREA_DESA_CREB_DESB>
{
  static constexpr unsigned int size = 1u;
  static constexpr unsigned int bra[size] = { 3 };
  static constexpr unsigned int ket[size] = { 3 };
  static const double value[size]; //  1
};

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<op::QC::CATEGORY OpCat>
//void add_operators (MPO<double,fermion>& mpo, const size_t& l, const size_t& r, const double& scale = 1.0)
void add_operators (btas::QSTArray<double,4,fermion>& mpo, const size_t& l, const size_t& r, const double& scale = 1.0)
{
  if(std::fabs(scale) < 1.0e-16) return;

//MPO<double,fermion>::block_type block(1,1,1,1);
  btas::TArray<double,4> block(1,1,1,1);
  for(size_t i = 0; i < OpComponent<OpCat>::size; ++i) {
    size_t p = OpComponent<OpCat>::bra[i];
    size_t q = OpComponent<OpCat>::ket[i];
//  mpo.is_allowed(l,p,q,r);
//  if(mpo.is_local(l,p,q,r)) {
//    if((mpo(p,q).qnum_array(0)[l].p() & 1)
//    && (mpo(p,q).qnum_array(1)[p].p() & 1))
//      block(0,0) =-OpComponent<OpCat>::value[i]*scale;
//    else
//      block(0,0) = OpComponent<OpCat>::value[i]*scale;
//    mpo(p,q)(l,r) = block;
      if(mpo.qshape(0)[l].parity() && mpo.qshape(1)[p].parity())
        block =-OpComponent<OpCat>::value[i]*scale;
      else
        block = OpComponent<OpCat>::value[i]*scale;
      mpo.insert(btas::shape(l,p,q,r),block);
//  }
  }
}

} // namespace mpsxx

#endif // __MPSXX_MPOGEN_ADD_OPERATORS_H
