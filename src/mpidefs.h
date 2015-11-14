#ifndef __MPSXX_MPI_DEFS_H
#define __MPSXX_MPI_DEFS_H

#include <iostream>

#ifndef _SERIAL

#include <boost/mpi.hpp>
typedef boost::mpi::communicator Communicator;

#define pout if (Communicator().rank() == 0) std::cout

#define xout std::cout << "DEBUG [ " << Communicator().rank() << " ] :: "

#else

/// Dummy communicator for serial run.
struct Communicator {
  size_t size () const { return 1; }
  size_t rank () const { return 0; }
};

#define pout std::cout

#define xout std::cout << "DEBUG :: "

#endif

#endif // __MPSXX_MPI_DEFS_H
