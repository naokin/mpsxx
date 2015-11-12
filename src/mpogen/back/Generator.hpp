#ifndef __MPSXX_MPO_GENERATOR_HPP
#define __MPSXX_MPO_GENERATOR_HPP

namespace mpsxx {

/// Generate MPO of the QC Hamiltonian
template<typename Tp, class Q>
MPO<Tp, Q> GenerateQCMPO (
      const std::vector<btas::Qshapes<Q>>& qN, ///< Physical Q# for each site
      const btas::TArray<Tp, 2>& oneint,
      const btas::TArray<Tp, 4>& twoint)
{
   // Get # of sites
   const size_t N = qN.size();

   // Left block operators
   BoundaryOperators lOps(0, N); // {0,...,N} in right block
   // Dot block operators
   BoundaryOperators sOps;
   // Right block operators
   BoundaryOperators rOps;

   size_t nnz = 0;

   for(size_t i = 0; i < N; ++i)
   {
      sOps.reset(i);
      rOps.reset(i+1, N); // {i+1,...,N} in right block

      bool do_swap_sweep = (lOps.loop_block() != rOps.loop_block());

      size_t nnz_local = 0;

      if(lOps.loop_block() == BoundaryOperators::Left)
      {
         // Dot with sys
         for(auto lOp = lOps.begin(); lOp != lOps.end(); ++lOp)
         {
            for(auto sOp = sOps.begin(); sOp != sOps.end(); ++sOp)
            {
               if(!do_swap_sweep)
               {
                  BaseOperators lsOp = (*lOp)*(*sOp);
                  auto rOp = rOps.find(lsOp);
                  if(rOp != rOps.end())
                  {
                     AddOpComponent(mpo[i], lOp, sOp, rOp, oneint, twoint);
                  }
               }
               else
               {
                  std::vector<BaseOperators> lsOps = Complement((*lOp)*(*sOp), rIndex);
                  for(auto lsOp = lsOps.begin(); lsOp != lsOps.end(); ++lsOp)
                  {
                     auto rOp = rOps.find(lsOp);
                     if(rOp != rOps.end())
                     {
                        AddOpComponent(mpo[i], lOp, sOp, rOp, oneint, twoint);
                     }
                  }
               }
            }
         }
      }
      else
      {
         // Dot with env
      }
   }
}

} // namespace mpsxx

#endif // __MPSXX_MPO_GENERATOR_HPP
