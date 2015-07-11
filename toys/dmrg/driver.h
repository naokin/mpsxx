#ifndef PROTOTYPE_DMRG_DRIVER_H
#define PROTOTYPE_DMRG_DRIVER_H

#include <cstdlib>
#include <legacy/DArray.h>

#define DEBUG(msg) { cout << "DEBUG:" << msg << endl; }

// canonicalize / 1-site algorithm
double canonicalize(bool forward, const btas::DArray<3>& wfn0,
                                        btas::DArray<3>& mps0,
                                        btas::DArray<3>& wfn1, int M = 0);

// canonicalize / 1-site algorithm / from system density matrix
double canonicalize(bool forward, const btas::DArray<4>& sdm0,
                                        btas::DArray<1>& val0,
                                        btas::DArray<3>& mps0,
                                        btas::DArray<3>& nul0, int M = 0);

// renormalize block-operator
void   renormalize (bool forward, const btas::DArray<4>& mpo0,
                                  const btas::DArray<3>& str0,
                                  const btas::DArray<3>& bra0,
                                  const btas::DArray<3>& ket0,
                                        btas::DArray<3>& str1);

// renormalize block-overlap
void   renormalize (bool forward, const btas::DArray<2>& str0,
                                  const btas::DArray<3>& bra0,
                                  const btas::DArray<3>& ket0,
                                        btas::DArray<2>& str1);

// compute diagonal elements
void   compute_h_diagonal        (const btas::DArray<4>& mpo0,
                                  const btas::DArray<3>& lstr,
                                  const btas::DArray<3>& rstr,
                                        btas::DArray<3>& diag);

// compute sigma vector
void   compute_sigma_vector      (const btas::DArray<4>& mpo0,
                                  const btas::DArray<3>& lstr,
                                  const btas::DArray<3>& rstr,
                                  const btas::DArray<3>& mps0,
                                        btas::DArray<3>& sgv0);

#endif // PROTOTYPE_DMRG_DRIVER_H
