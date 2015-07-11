
#include <qcdmrg/build_random_system.h>

void mpsxx::qcdmrg::Build_random_system (const mpsxx::qcdmrg::DMRG_input& input)
{
#ifdef _ENABLE_MPI
  boost::mpi::communicator world;
#endif // _ENABLE_MPI

  Build_QC_operator(input);
  Build_random_mps(input);

#ifdef _ENABLE_MPI
  world.barrier();
#endif // _ENABLE_MPI
}

void mpsxx::qcdmrg::Build_random_mps (const mpsxx::qcdmrg::DMRG_input& input)
{
#ifdef _ENABLE_MPI
  boost::mpi::communicator world;
  if(world.rank() > 0) return;
#endif // _ENABLE_MPI

  const size_t& N = input.N_sites();
  const size_t& M = input.M_start();

  MPS_skeleton<Quanta> Dq(N);

  // NOTE:
  // void MPO_component<double>::load(std::string, std::string, int, int);
  // [0] name suffix
  // [1] scratch file path
  // [2] process #
  // [3] site #
  // return: working rep. of i-th MPO

  for(size_t i = 0; i < N; ++i)
  {
    MPO_component<double> mpo;
    mpo.load("qc", input.prefix(), 0, i);

    Dq[i].phys() = mpo.index();
  }

  Dq.set_quanta(M);

  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  // NOTE:
  // void MPS_component<double>::load(std::string, std::string, int, int, int);
  // void MPS_component<double>::save(std::string, std::string, int, int, int);
  // [0] name suffix
  // [1] scratch file path
  // [2] process #
  // [3] site #
  // [4] state #
  // return: working rep. of i-th MPO

  MPS_component<double> wfnc;

  for(size_t i = 0; i < N-1; ++i)
  {
    wfnc.resize(Dq[i].left(), Dq[i].phys(), Dq[i].right());
    wfnc.fill(std::bind(dist, mt19937()));
    wfnc.normalize();
    // NOTE:
    // Tensor<double> MPS_component<double>::canonicalize(int, int);
    // [0] 0:right, 1:left
    // [1] # states to be kept, 0 means no truncation
    // return: residual tensor, ie. gauge matrix
    wfnc.canonicalize(1, 0);
    wfnc.save("lmps", input.prefix(), 0, i, 0);
  }

  wfnc.resize(Dq[N-1].left(), Dq[N-1].phys(), Dq[N-1].right());
  wfnc.fill(std::bind(dist, mt19937()));

  for(size_t i = N-1; i > 0; --i)
  {
    wfnc.normalize();
    wfnc.save("wfnc", input.prefix(), 0, i, 0);

    Gauge_matrix<double> g = wfnc.canonicalize(0, 0);
    wfnc.save("rmps", input.prefix(), 0, i, 0);

    MPS_component<double> lmps;
    lmps.load("lmps", input.prefix(), 0, i-1, 0);

    wfnc = lmps * g;
  }

  wfnc.normalize();
  wfnc.save("wfnc", input.prefix(), 0, 0, 0);
}

