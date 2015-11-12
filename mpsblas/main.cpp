#include <iostream>
#include <iomanip>

int main (int argc, char* argv[])
{
  using namespace mpsxx;

  size_t N = 100;
  size_t D =  20;

  MPS_vector<double> wav; // wavefunction
  MPO_vector<double> ham; // hamiltonian
  MPO_vector<double> hpm; // hamiltonian (preconditioned)

  make_random_mps(N,wav,D);
  make_heisenberg(N,ham,0.5); // non-local spin model (to be solved)
  make_heisenberg(N,hpm); // local spin model (preconditioned?)

  left_normalize(wav);

  size_t MAX_ITER = 10;
  size_t MAX_RITX = 20;

  size_t conv = 0;
  size_t root = 1;
  size_t iter = 0;

  std::vector<MPS_vector<double>> trial(MAX_RITZ);

  while(conv < root && iter < MAX_ITER) {
    for(size_t m = 1; m <= MAX_RITZ; ++m) {
      // compute small Hamiltonian matrix
      btas::Tensor<double,2,CblasRowMajor> H(m,m);
      btas::Tensor<double,2,CblasRowMajor> S(m,m);
      for(size_t i = 0; i < m; ++i) {
        H(i,i) = expectation(trial[i],ham,trial[i]);
        S(i,i) = overlap(trial[i],trial[i]);
        for(size_t j = 0; j < i; ++j) {
          double Hij = expectation(trial[i],ham,trial[j]);
          H(i,j) = Hij;
          H(j,i) = Hij;
          double Sij = overlap(trial[i],trial[j]);
          S(i,j) = Sij;
          S(j,i) = Sij;
        }
      }
      // solve eigenvalue problem to obtain Ritz value & vector
      btas::Tensor<double,2,CblasRowMajor> rvec;
      btas::Tensor<double,1,CblasRowMajor> rval;
      btas::syev('V','U',H,rval,rvec);
      e = rval(0);

      std::vector<MPS_vector<double>> trial_save(m);
      for(size_t i = 0; i < m; ++i) {
        trial_save[i] = trial[i];
        scale(rvec(i,i),trial[i]);
      }
      for(size_t i = 0; i < m; ++i)
        for(size_t j = 0; j < m; ++j)
          if(i != j) scaled_add(rvec(j,i),trial_save[j],trial[i]);
    }
  }
  return 0;
}
