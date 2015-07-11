#ifndef __MPSXX_MPO_INTEGRAL_COMPONENT_HPP
#define __MPSXX_MPO_INTEGRAL_COMPONENT_HPP

namespace mpsxx {

/// Return 1-particle integral component to construct complementary operator
template<typename Tp>
Tp Int1eComponent (const BaseOperator& lsrOp, const btas::TArray<Tp, 2>& oneint)
{
   assert(lsrOp.index().size() == 2);

   size_t i = lsrOp.index()[0];
   size_t j = lsrOp.index()[1];

   Tp value = static_cast<Tp>(0);

   if     (lsrOp.type() == BaseOperator::CREA_DESA ||
           lsrOp.type() == BaseOperator::CREB_DESB)
   {
      value += oneint(i, j);
   }
   else if(lsrOp.type() == BaseOperator::DESA_CREA ||
           lsrOp.type() == BaseOperator::DESB_CREB)
   {
      value -= oneint(i, j);
   }

   return value;
}

/// Return 2-particle integral component to construct complementary operator
///
/// Independent integral contributions
///
/// F: C[i,a] x C[j,a] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * D[k,a] x D[l,a] : 5
/// E: C[i,a] x C[j,b] -- ( - V[i,j,k,l] )              * D[k,a] x D[l,b] : 4
///                    -- ( + V[i,j,l,k] )              * D[k,b] x D[l,a] : 1
/// B: C[i,b] x C[j,a] -- ( - V[i,j,k,l] )              * D[k,b] x D[l,a] : 1
///                    -- ( + V[i,j,l,k] )              * D[k,a] x D[l,b] : 4
/// A: C[i,b] x C[j,b] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * D[k,b] x D[l,b] : 0
///
/// D: C[i,a] x D[j,a] -- ( - V[i,j,l,k] + V[i,k,j,l] ) * C[k,a] x D[l,a] : D
///                    -- ( + V[i,j,k,l] - V[i,k,j,l] ) * D[k,a] x C[l,a] : 7
///                    -- ( + V[i,k,j,l] )              * C[k,b] x D[l,b] : 8
///                    -- ( - V[i,k,j,l] )              * D[k,b] x C[l,b] : 2
/// C: C[i,a] x D[j,b] -- ( - V[i,j,l,k] )              * C[k,b] x D[l,a] : 9
///                    -- ( + V[i,j,k,l] )              * D[k,a] x C[l,b] : 6
/// 9: C[i,b] x D[j,a] -- ( - V[i,j,l,k] )              * C[k,a] x D[l,b] : C
///                    -- ( + V[i,j,k,l] )              * D[k,b] x C[l,a] : 3
/// 8: C[i,b] x D[j,b] -- ( - V[i,j,l,k] + V[i,k,j,l] ) * C[k,b] x D[l,b] : 8
///                    -- ( + V[i,j,k,l] - V[i,k,j,l] ) * D[k,b] x C[l,b] : 2
///                    -- ( + V[i,k,j,l] )              * C[k,a] x D[l,a] : D
///                    -- ( - V[i,k,j,l] )              * D[k,a] x C[l,a] : 7
///
/// 7: D[i,a] x C[j,a] -- ( + V[i,j,k,l] - V[i,k,j,l] ) * C[k,a] x D[l,a] : D
///                    -- ( - V[i,j,l,k] + V[i,k,j,l] ) * D[k,a] x C[l,a] : 7
///                    -- ( - V[i,k,j,l] )              * C[k,b] x D[l,b] : 8
///                    -- ( + V[i,k,j,l] )              * D[k,b] x C[l,b] : 2
/// 6: D[i,a] x C[j,b] -- ( + V[i,j,k,l] )              * C[k,a] x D[l,b] : C
///                    -- ( - V[i,j,l,k] )              * D[k,b] x C[l,a] : 3
/// 3: D[i,b] x C[j,a] -- ( + V[i,j,k,l] )              * C[k,b] x D[l,a] : 9
///                    -- ( - V[i,j,l,k] )              * D[k,a] x C[l,b] : 6
/// 2: D[i,b] x C[j,b] -- ( + V[i,j,k,l] - V[i,k,j,l] ) * C[k,b] x D[l,b] : 8
///                    -- ( - V[i,j,l,k] + V[i,k,j,l] ) * D[k,b] x C[l,b] : 2
///                    -- ( - V[i,k,j,l] )              * C[k,a] x D[l,a] : D
///                    -- ( + V[i,k,j,l] )              * D[k,a] x C[l,a] : 7
///
/// 5: D[i,a] x D[j,a] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * C[k,a] x C[l,a] : F
/// 4: D[i,a] x D[j,b] -- ( - V[i,j,k,l] )              * C[k,a] x C[l,b] : E
///                    -- ( + V[i,j,l,k] )              * C[k,b] x C[l,a] : B
/// 1: D[i,b] x D[j,a] -- ( - V[i,j,k,l] )              * C[k,b] x C[l,a] : B
///                    -- ( + V[i,j,l,k] )              * C[k,a] x C[l,b] : E
/// 0: D[i,b] x D[j,b] -- ( + V[i,j,l,k] - V[i,j,k,l] ) * C[k,b] x C[l,b] : A
///
template<typename Tp>
Tp Int2eComponent (const BaseOperator& lsrOp, const btas::TArray<Tp, 4>& twoint)
{
   assert(lsrOp.index().size() == 4);

   size_t i = lsrOp.index()[0];
   size_t j = lsrOp.index()[1];
   size_t k = lsrOp.index()[2];
   size_t l = lsrOp.index()[3];

   Tp value = static_cast<Tp>(0);

   switch(lsrOp.type())
   {
      case BaseOperator::CREA_CREA_DESA_DESA:
      case BaseOperator::CREA_CREB_DESB_DESA:
      case BaseOperator::CREB_CREA_DESA_DESB:
      case BaseOperator::CREB_CREB_DESB_DESB:

      case BaseOperator::DESA_DESA_CREA_CREA:
      case BaseOperator::DESA_DESB_CREB_CREA:
      case BaseOperator::DESB_DESA_CREA_CREB:
      case BaseOperator::DESB_DESB_CREB_CREB:

         value += twoint(i, j, l, k); break;

      case BaseOperator::CREA_DESA_CREA_DESA:
      case BaseOperator::CREA_DESB_CREB_DESA:
      case BaseOperator::CREB_DESA_CREA_DESB:
      case BaseOperator::CREB_DESB_CREB_DESB:

      case BaseOperator::DESA_CREA_DESA_CREA:
      case BaseOperator::DESA_CREB_DESB_CREA:
      case BaseOperator::DESB_CREA_DESA_CREB:
      case BaseOperator::DESB_CREB_DESB_CREB:

         value -= twoint(i, j, l, k); break;

      default: break;
   }

   switch(lsrOp.type())
   {
      case BaseOperator::CREA_DESA_DESA_CREA:
      case BaseOperator::CREA_DESB_DESA_CREB:
      case BaseOperator::CREB_DESA_DESB_CREA:
      case BaseOperator::CREB_DESB_DESB_CREB:

      case BaseOperator::DESA_CREA_CREA_DESA:
      case BaseOperator::DESA_CREB_CREA_DESB:
      case BaseOperator::DESB_CREA_CREB_DESA:
      case BaseOperator::DESB_CREB_CREB_DESB:

         value += twoint(i, j, k, l); break;

      case BaseOperator::CREA_CREA_DESA_DESA:
      case BaseOperator::CREA_CREB_DESA_DESB:
      case BaseOperator::CREB_CREA_DESB_DESA:
      case BaseOperator::CREB_CREB_DESB_DESB:

      case BaseOperator::DESA_DESA_CREA_CREA:
      case BaseOperator::DESA_DESB_CREA_CREB:
      case BaseOperator::DESB_DESA_CREB_CREA:
      case BaseOperator::DESB_DESB_CREB_CREB:

         value -= twoint(i, j, k, l); break;

      default: break;
   }

   switch(lsrOp.type())
   {
      case BaseOperator::CREA_DESA_CREA_DESA:
      case BaseOperator::CREA_DESA_CREB_DESB:
      case BaseOperator::CREB_DESB_CREA_DESA:
      case BaseOperator::CREB_DESB_CREB_DESB:

      case BaseOperator::DESA_CREA_DESA_CREA:
      case BaseOperator::DESA_CREA_DESB_CREB:
      case BaseOperator::DESB_CREB_DESA_CREA:
      case BaseOperator::DESB_CREB_DESB_CREB:

         value += twoint(i, k, j, l); break;

      case BaseOperator::CREA_DESA_DESA_CREA:
      case BaseOperator::CREA_DESA_DESB_CREB:
      case BaseOperator::CREB_DESB_DESA_CREA:
      case BaseOperator::CREB_DESB_DESB_CREB:

      case BaseOperator::DESA_CREA_CREA_DESA:
      case BaseOperator::DESA_CREA_CREB_DESB:
      case BaseOperator::DESB_CREB_CREA_DESA:
      case BaseOperator::DESB_CREB_CREB_DESB:

         value -= twoint(i, k, j, l); break;

      default: break;
   }

   return value;
}

} // namespace mpsxx

#endif // __MPSXX_MPO_INTEGRAL_COMPONENT_HPP
