
//! Random number generator
double rgen() { return 2.0*(static_cast<double>(rand())/RAND_MAX)-1.0; }

void set_quantum_blocks(const MpStorages& sites, const Quantum& qt, std::vector<Qshapes<>>& qb, int QMAX_SIZE)
{
  int L = sites.size();

  // physical index
  Qshapes<> qp;
  // 0 quantum number
  Qshapes<> qz(1, Quantum::zero());

  qb.resize(L);

  // quantum blocks from the entire Fock space
  qb[0] =-sites[0].mpo.qshape(2);
  for(int i = 1; i < L-1; ++i) {
    qp    =-sites[i].mpo.qshape(2);
    qb[i] = qb[i-1] & qp; // get unique elements of { q(left) x q(phys) }
  }
  qp      =-sites[L-1].mpo.qshape(2);
  qb[L-1] = Qshapes<>(1, qt);

  // reduce zero quantum blocks
  for(int i = L-1; i > 0; --i) {
    qp =-sites[i].mpo.qshape(2);
    Qshapes<>& ql = qb[i-1];
    Qshapes<>& qr = qb[i];

    // check non-zero for each ql index
    Qshapes<>::iterator lt = ql.begin();
    while(lt != ql.end()) {
      bool non_zero = false;
      for(int p = 0; p < qp.size(); ++p) {
        for(int r = 0; r < qr.size(); ++r) {
          non_zero |= (qr[r] == (qp[p] * (*lt)));
        }
      }
      if(non_zero)
        ++lt;
      else
        ql.erase(lt);
    }
    // further reduction
    if(QMAX_SIZE > 0 && ql.size() > QMAX_SIZE) {
      int offs = (ql.size() - QMAX_SIZE) / 2;
      ql = Qshapes<>(ql.begin()+offs, ql.begin()+offs+QMAX_SIZE);
    }
  }
}

void prototype::initialize(MpStorages& sites, const Quantum& qt, int M)
{
  int L = sites.size();

  // physical index
  Qshapes<> qp;
  Dshapes   dp;
  // left state index
  Qshapes<> ql;
  Dshapes   dl;
  // right state index
  Qshapes<> qr;
  Dshapes   dr;
  // 0 quantum number
  Qshapes<> qz(1, Quantum::zero());
  Dshapes   dz(qz.size(), 1);

  // non-zero quantum numbers for each site
  std::vector<Qshapes<>> qb;

  int max_size = 20;
  set_quantum_blocks(sites, qt, qb, max_size);

  //
  // create random wavefunction
  //

  int M0 = 1;
  int Mx = M;

  TVector<Qshapes<>, 3> qshape;
  TVector<Dshapes,   3> dshape;

  qr = qz;
  dr = Dshapes(qr.size(), 1);

  for(int i = 0; i < L-1; ++i) {
    // physical index is taken from mpo's ket index
    qp =-sites[i].mpo.qshape(2);
    dp = Dshapes(qp.size(), 1);
    // left index equals to previous right index
    ql = qr;
    dl = dr;
    // non-zero quantum numbers for site i
    qr = qb[i];
    dr = Dshapes(qr.size(), M0);

    qshape = make_array( ql, qp,-qr);
    dshape = make_array( dl, dp, dr);
    sites[i].wfnc[0].resize(Quantum::zero(), qshape, dshape);
    sites[i].wfnc[0].generate(rgen);
    sites[i].save(i);
  }

  qp =-sites[L-1].mpo.qshape(2);
  dp = Dshapes(qp.size(), 1);
  ql = qr;
  dl = dr;
  qr = qz;
  dr = dz;
  qshape = make_array( ql, qp,-qr);
  dshape = make_array( dl, dp, dr);
  sites[L-1].wfnc[0].resize(qt, qshape, dshape);
  sites[L-1].wfnc[0].generate(rgen);

  //
  // canonicalize & renormalize
  //

  qshape = make_array( qz, qz, qz);
  dshape = make_array( dz, dz, dz);

  sites[L-1].ropr[0].resize(Quantum::zero(), qshape, dshape);
  sites[L-1].ropr[0] = 1.0;

  for(int i = L-1; i > 0; --i) {
    sites[i-1].load(i-1);
    QSDnormalize(sites[i].wfnc[0]);
    Canonicalize(0, sites[i].wfnc[0], sites[i].rmps[0], Mx);
    QSDcopy(sites[i-1].wfnc[0], sites[i-1].lmps[0]);
    sites[i-1].wfnc[0].clear();
    ComputeGuess(0, sites[i].rmps[0], sites[i].wfnc[0], sites[i-1].lmps[0], sites[i-1].wfnc[0]);
    sites[i-1].ropr[0].clear();
    Renormalize (0, sites[i].mpo, sites[i].ropr[0], sites[i].rmps[0], sites[i].rmps[0], sites[i-1].ropr[0]);
    sites[i].save(i);
  }

  QSDnormalize(sites[0].wfnc[0]);
  sites[0].lopr[0].resize(Quantum::zero(), qshape, dshape);
  sites[0].lopr[0] = 1.0;
  sites[0].save(0);
}

double prototype::optimize_onesite(bool forward, MpSite& sysdot, MpSite& envdot, int M)
{
  boost::function<void(const QSDArray<3>&, QSDArray<3>&)>
  f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, sysdot.lopr[0], sysdot.ropr[0], _1, _2);
  QSDArray<3> diag(sysdot.wfnc[0].q(), sysdot.wfnc[0].qshape());
  ComputeDiagonal(sysdot.mpo, sysdot.lopr[0], sysdot.ropr[0], diag);
  double energy = davidson::diagonalize(f_contract, diag, sysdot.wfnc[0]);

  if(forward) {
    Canonicalize(1, sysdot.wfnc[0], sysdot.lmps[0], M);
    envdot.wfnc[0].clear();
    ComputeGuess(1, sysdot.lmps[0], sysdot.wfnc[0], envdot.rmps[0], envdot.wfnc[0]);
    envdot.lopr[0].clear();
    Renormalize (1, sysdot.mpo,  sysdot.lopr[0], sysdot.lmps[0], sysdot.lmps[0], envdot.lopr[0]);
  }
  else {
    Canonicalize(0, sysdot.wfnc[0], sysdot.rmps[0], M);
    envdot.wfnc[0].clear();
    ComputeGuess(0, sysdot.rmps[0], sysdot.wfnc[0], envdot.lmps[0], envdot.wfnc[0]);
    envdot.ropr[0].clear();
    Renormalize (0, sysdot.mpo,  sysdot.ropr[0], sysdot.rmps[0], sysdot.rmps[0], envdot.ropr[0]);
  }

  return energy;
}

double prototype::optimize_twosite(bool forward, MpSite& sysdot, MpSite& envdot, int M)
{
  QSDArray<4> wfnc;
  QSDArray<4> diag;
  boost::function<void(const QSDArray<4>&, QSDArray<4>&)> f_contract;
  if(forward) {
    QSDgemm(NoTrans, NoTrans, 1.0, sysdot.wfnc[0], envdot.rmps[0], 1.0, wfnc);
    f_contract = boost::bind(ComputeSigmaVector, sysdot.mpo, envdot.mpo, sysdot.lopr[0], envdot.ropr[0], _1, _2);
    diag.resize(wfnc.q(), wfnc.qshape());
    ComputeDiagonal(sysdot.mpo, envdot.mpo, sysdot.lopr[0], envdot.ropr[0], diag);
  }
  else {
    QSDgemm(NoTrans, NoTrans, 1.0, envdot.lmps[0], sysdot.wfnc[0], 1.0, wfnc);
    f_contract = boost::bind(ComputeSigmaVector, envdot.mpo, sysdot.mpo, envdot.lopr[0], sysdot.ropr[0], _1, _2);
    diag.resize(wfnc.q(), wfnc.qshape());
    ComputeDiagonal(envdot.mpo, sysdot.mpo, envdot.lopr[0], sysdot.ropr[0], diag);
  }

  double energy = davidson::diagonalize(f_contract, diag, wfnc);

  if(forward) {
    Canonicalize(1,        wfnc, sysdot.lmps[0], envdot.wfnc[0], M);
    envdot.lopr[0].clear();
    Renormalize (1, sysdot.mpo,  sysdot.lopr[0], sysdot.lmps[0], sysdot.lmps[0], envdot.lopr[0]);
  }
  else {
    Canonicalize(0,        wfnc, sysdot.rmps[0], envdot.wfnc[0], M);
    envdot.ropr[0].clear();
    Renormalize (0, sysdot.mpo,  sysdot.ropr[0], sysdot.rmps[0], sysdot.rmps[0], envdot.ropr[0]);
  }

  return energy;
}

double prototype::dmrg_sweep(MpStorages& sites, DMRG_ALGORITHM algo, int M)
{
  int    L    = sites.size();
  double emin = 1.0e8;
  // fowrad sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  sites[0].load(0);
  for(int i = 0; i < L-1; ++i) {
    sites[i+1].load(i+1);
    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(1, sites[i], sites[i+1], M);
    else                eswp = optimize_twosite(1, sites[i], sites[i+1], M);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    sites[i].save(i);
  
  // backward sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tBACKWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  for(int i = L-1; i > 0; --i) {
    sites[i-1].load(i-1);
    // diagonalize
    double eswp;
    if(algo == ONESITE) eswp = optimize_onesite(0, sites[i], sites[i-1], M);
    else                eswp = optimize_twosite(0, sites[i], sites[i-1], M);
    if(eswp < emin) emin = eswp;
    // print result
    cout.precision(16);
    cout << "\t\t\tEnergy = " << setw(24) << fixed << eswp << endl;
    sites[i].save(i);
  }
  sites[0].save(0);
  return emin;
}

double prototype::dmrg(const DmrgInput& input, MpStorages& sites)
{
  const int L = input.N_sites;
  const int M = input.N_max_states;

  sites.resize(L, MpSite("state-0", 1, input.prefix));

  switch(input.model) {
    case HEISENBERG:
      Heisenberg::construct_mpo(sites, input.heisenberg);
      initialize(sites, FermiQuantum(0, input.N_spins), M);
      break;
    case HUBBARD:
      Hubbard::construct_mpo(sites, input.hubbard);
      initialize(sites, FermiQuantum(input.N_elecs, input.N_spins), M);
      break;
    default:
      cout << "\tUnknown model type specified" << endl;
      return 0.0;
  }

  double esav = 1.0e8;

  //
  // two-site optimization
  //
  for(int iter = 0; iter < 20; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(sites, TWOSITE, M);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << "\t====================================================================================================" << endl;
    cout << endl;
    esav = eswp;
    if(fabs(edif) < 1.0e-8) break;
  }

  //
  // one-site optimization
  //
  for(int iter = 0; iter < 20; ++iter) {
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] "   << endl;
    cout << "\t====================================================================================================" << endl;
    double eswp = dmrg_sweep(sites, ONESITE, M);
    double edif = eswp - esav;
    cout << "\t====================================================================================================" << endl;
    cout << "\t\tSWEEP ITERATION [ " << setw(4) << iter << " ] FINISHED" << endl;
    cout.precision(16);
    cout << "\t\t\tSweep Energy = " << setw(24) << fixed << eswp << " ( delta E = ";
    cout.precision(2);
    cout << setw(8) << scientific << edif << " ) " << endl;
    cout << "\t====================================================================================================" << endl;
    cout << endl;
    esav = eswp;
    if(fabs(edif) < 1.0e-8) break;
  }

  return esav;
}

void prototype::analysis(const DmrgInput& input, MpStorages& sites)
{
  const int L = input.N_sites;

  // fowrad sweep
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "\t\t\tFORWARD SWEEP" << endl;
  cout << "\t++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  for(int i = 0; i < L-1; ++i) {
    sites[i].load(i);
    cout << "\tanalyzing transfer operator at site [ " << setw(2) << i << " ] " << endl;
//  AnalyzeTransferOperator(sites[i].wfnc[0], sites[i].wfnc[0]);
    AnalyzeTransferOperator(sites[i].lmps[0], sites[i].lmps[0]);
    sites[i].free();
  }
}

