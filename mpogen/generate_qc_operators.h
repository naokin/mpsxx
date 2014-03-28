#ifndef _MPSXX_CXX11_GENERATE_QC_OPERATORS_H
#define _MPSXX_CXX11_GENERATE_QC_OPERATORS_H 1

#include <string>

#include <btas/DENSE/DArray.h>

#include <MPSblas.h>
#include <symmetry/Fermion/Quantum.h>

namespace mpsxx     {

   namespace fermionic {

      void generate_qc_operators
         (MPO<double, Quantum>& mpos, const btas::DArray<2>& oneint, const btas::DArray<4>& twoint, bool enable_swap_sweep_dir = false, const std::string& prefix = "./");

   }; // namespace fermionic

}; // namespace mpsxx

#endif // _MPSXX_CXX11_GENERATE_QC_OPERATORS_H
