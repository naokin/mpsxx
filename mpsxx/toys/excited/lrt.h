

{
  for(int i = L-1; i > 0; --i) {
    wfnc1[i] = random_gen;
    Dnormalize(wfnc1[i]);
    DArray<4> dens1;
    compute_perturbed_density(0, wfnc0[i], wfnc1[i], dens1);
    compute_perturbed_mpstate(0, dens1, lval0[i], lmps0[i], lnul0[i], lmps1[i]);
    compute_perturbed_opblock(0, mpo[i], lstr0[i], lstr1[i], lmps0[i], lmps1[i], lstr1[i+1]);
  }
}


double MpsDot(const MpStates& wfnc0, const MpStates& wfnc1)
{
  double value;
  for(int i = 0; i < wfnc0.size(); ++i) value += Ddot(wfnc0[i], wfnc1[i]);
  return value;
}

void MpsAxpy(double alpha, const MpStates& wfnc0, MpStates& wfnc1)
{
  for(int i = 0; i < wfnc0.size(); ++i) Daxpy(alpha, wfnc0[i], wfnc1[i]);
}

