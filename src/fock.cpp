
#include "fock.h"

void mpsxx::twopdm_to_onepdm
(const size_t& Norbs,
 const size_t& Nelec,
 const btas::DArray<4>& twopdm,
       btas::DArray<2>& onepdm)
{
  onepdm.resize(Norbs, Norbs);
  onepdm.fill(0.0);
  double scale = 2.0 / (Nelec - 1);
  for(size_t i = 0; i < Norbs; ++i) {
    for(size_t j = 0; j < Norbs; ++j) {
      for(size_t k = 0; k < Norbs; ++k) {
        onepdm(i, j) += twopdm(i, k, k, j);
      }
      onepdm(i, j) *= scale;
    }
  }
}

void mpsxx::compute_fock_matrix
(const size_t& Norbs,
 const btas::DArray<2>& oneint,
 const btas::DArray<4>& twoint,
 const btas::DArray<2>& onepdm,
 const btas::DArray<4>& twopdm,
       btas::DArray<2>& fock)
{
  fock.resize(Norbs, Norbs);
  fock.fill(0.0);
  for(size_t i = 0; i < Norbs; ++i) {
    for(size_t j = 0; j < Norbs; ++j) {
      double& fij = fock(i, j);
      for(size_t p = 0; p < Norbs; ++p) {
        fij += oneint(i, p) * onepdm(j, p);
        for(size_t q = 0; q < Norbs; ++q) {
          for(size_t r = 0; r < Norbs; ++r) {
            fij += 2.0 * twoint(i, p, r, q) * twopdm(j, p, q, r);
          }
        }
      }
    }
  }
}

void mpsxx::compute_gradient
(const size_t& Norbs,
 const btas::DArray<2>& fock,
       btas::DArray<2>& grad)
{
  grad.resize(Norbs, Norbs);
  grad.fill(0.0);
  for(size_t i = 0; i < Norbs; ++i) {
    for(size_t j = 0; j < Norbs; ++j) {
      grad(i, j) = 2.0 * (fock(i, j) - fock(j, i));
    }
  }
}

void mpsxx::compute_hessian
(const size_t& Norbs,
 const btas::DArray<2>& fock,
 const btas::DArray<2>& oneint,
 const btas::DArray<4>& twoint,
 const btas::DArray<2>& onepdm,
 const btas::DArray<4>& twopdm,
       btas::DArray<4>& hess)
{
  hess.resize(Norbs, Norbs, Norbs, Norbs);
  hess.fill(0.0);

  btas::DArray<4> hw(Norbs, Norbs, Norbs, Norbs);
  for(size_t i = 0; i < Norbs; ++i) {
    for(size_t j = 0; j < Norbs; ++j) {
      for(size_t k = 0; k < Norbs; ++k) {
        hw(i, k, k, j) = fock(i, j) + fock(j, i);
      }
    }
  }
  for(size_t i = 0; i < Norbs; ++i) {
    for(size_t j = 0; j < Norbs; ++j) {
      for(size_t k = 0; k < Norbs; ++k) {
        for(size_t l = 0; l < Norbs; ++l) {
          double wijkl = 0.0;
          for(size_t p = 0; p < Norbs; ++p) {
            for(size_t q = 0; q < Norbs; ++q) {
              wijkl += twopdm(k, p, q, i) * twoint(j, p, l, q);
              wijkl += twopdm(k, p, i, q) * twoint(j, l, p, q);
              wijkl += twopdm(k, i, p, q) * twoint(j, l, p, q);
            }
          }
          hw(i, j, k, l) += 2.0 * onepdm(i, k) * oneint(j, l) + 4.0 * wijkl;
        }
      }
    }
  }
  for(size_t i = 0; i < Norbs; ++i) {
    for(size_t j = 0; j < Norbs; ++j) {
      for(size_t k = 0; k < Norbs; ++k) {
        for(size_t l = 0; l < Norbs; ++l) {
          hess(i, j, k, l) = hw(i, j, k, l) - hw(j, i, k, l) - hw(i, j, l, k) + hw(j, i, l, k);
        }
      }
    }
  }
}

double mpsxx::compute_energy
(const size_t& Norbs,
 const btas::DArray<2>& oneint,
 const btas::DArray<4>& twoint,
 const btas::DArray<2>& onepdm,
 const btas::DArray<4>& twopdm)
{
  return btas::Ddot(oneint, onepdm) + btas::Ddot(twoint, twopdm);
}

