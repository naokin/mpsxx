


// Terms from \hat{H}\ket{\Psi_0} are stored as a vector

// walkers
// distributed
std::vector<Walker<MPS>> walkers;

// H[0] contains I
// replicated
std::vector<DpOperator<Q>> H;

// H_Psi0[0] contains Psi0
// replicated
std::vector<MpState<Q>> H_Psi0;

// t: time step. "current time" = t*dt, "total time" = tmax*dt
for (size_t t = 0; t < tmax; ++t) {

   // ==============================================================================================
   // propagation steps
   // ==============================================================================================

   double sum_of_walker_coeff = 0.0;

   // iwk: i-th walker, nwk: # of walkers
   for (size_t iwk = 0; iwk < nwk; ++iwk) {

      double& Siwk = walkerOverlap[iwk]

      // m: m-th term of "H" or "H_Psi0"
      // ipdf: integrated prob. dist. func.
      size_t n_term = H_Psi0.size();
      pdf[0] = 1.0;
      for (size_t m = 1; m < n_term; ++m) {
         // dt: duration, Sn: overlap b/w Psi0 & walkers
         pdf[m] = pdf[m-1]+std::max(-dt*H_Psi0[m]*walkers[iwk]/Siwk, 0.0);
      }
      double norm = pdf[n_term-1];

      // draw random number [0,1) and multiplies to pdf[n_term-1] (i.e. norm of pdfs)
      double pdfx = rgen()*pdf[n_term-1];
      // binary search to get term to be evolved
      size_t j_term = std::lower_bound(pdf.begin(), pdf.end(), pdfx) - pdf.begin();

      // if (j_term == 0) do nothing!

      if (j_term > 0) {
         // apply H[j_term]
         walkers[iwk] = walkers[iwk]*H[j_term];
         // normalize
         Normalize(1, walkers[iwk]); // left-normalization
         // compt. overlap
         Siwk = walkers[iwk]*H_Psi0[0];
         if (Siwk < 0.0) {
            scal(-1.0, walkers[iwk]);
            Siwk = -Siwk;
         }
         walkers[iwk].Coeff() *= norm;
         sum_of_walker_coeff += walkers[iwk].Coeff();
      }
   }

   // ==============================================================================================
   // energy calc.
   // ==============================================================================================

   // ==============================================================================================
   // population control
   // ==============================================================================================
}
