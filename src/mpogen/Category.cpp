#include "Category.h"

namespace mpsxx {
namespace op {
namespace QC {

fermion QNumTable[SIZE];

std::string LabTable[SIZE];

CATEGORY ConjTable[SIZE];

CATEGORY MultTable[SIZE][SIZE];

} // namespace QC

RDM::CATEGORY QC2RDM[QC::SIZE];

} // namespace op
} // namespace mpsxx

// Right-Operator Table (swap_comp_block == false)
// -----+------------------------------+---------------------+
// L \ S| I    H    Ck   Dk   Rl   Sl  | CkCk CkDk DkCk DkDk |
// -----+------------------------------+---------------------+
// I    | I    H    Ck   Dk   Rl   Sl  | CkCk CkDk DkCk DkDk |
// H    | H    0    0    0    0    0   | 0    0    0    0    |
// Ci   | Ci   0    CiCk CiDk 0    H   | 0    Rl   Rl   Sl   |
// Di   | Di   0    DiCk DiDk H    0   | Rl   Sl   Sl   0    |
// Rl   | Rl   0    0    H    0    0   | 0    0    0    0    |
// Sl   | Sl   0    H    0    0    0   | 0    0    0    0    |
// -----+------------------------------+---------------------+
// CiCj | CiCj 0    0    Rl   0    0   | 0    0    0    H    |
// CiDj | CiDj 0    Rl   Sl   0    0   | 0    H    H    0    |
// DiCj | DiCj 0    Rl   Sl   0    0   | 0    H    H    0    |
// DiDj | DiDj 0    Sl   0    0    0   | H    0    0    0    |
// -----+------------------------------+---------------------+

void mpsxx::op::QC::setTables ()
{
  // set QNumTable

  QNumTable[I] = fermion( 0, 0);

  QNumTable[CREA] = fermion( 1, 1);
  QNumTable[CREB] = fermion( 1,-1);
  QNumTable[DESA] = fermion(-1,-1);
  QNumTable[DESB] = fermion(-1, 1);

  QNumTable[CREA_CREA] = fermion( 2, 2);
  QNumTable[CREA_CREB] = fermion( 2, 0);
  QNumTable[CREB_CREA] = fermion( 2, 0);
  QNumTable[CREB_CREB] = fermion( 2,-2);
  QNumTable[CREA_DESA] = fermion( 0, 0);
  QNumTable[CREA_DESB] = fermion( 0, 2);
  QNumTable[CREB_DESA] = fermion( 0,-2);
  QNumTable[CREB_DESB] = fermion( 0, 0);
  QNumTable[DESA_CREA] = fermion( 0, 0);
  QNumTable[DESA_CREB] = fermion( 0,-2);
  QNumTable[DESB_CREA] = fermion( 0, 2);
  QNumTable[DESB_CREB] = fermion( 0, 0);
  QNumTable[DESA_DESA] = fermion(-2,-2);
  QNumTable[DESA_DESB] = fermion(-2, 0);
  QNumTable[DESB_DESA] = fermion(-2, 0);
  QNumTable[DESB_DESB] = fermion(-2, 2);

  QNumTable[H] = fermion( 0, 0);

  QNumTable[CREA_COMP] = fermion( 1, 1);
  QNumTable[CREB_COMP] = fermion( 1,-1);
  QNumTable[DESA_COMP] = fermion(-1,-1);
  QNumTable[DESB_COMP] = fermion(-1, 1);

  // set LabTable

  LabTable[I] = "I";

  LabTable[CREA] = "CREA";
  LabTable[CREB] = "CREB";
  LabTable[DESA] = "DESA";
  LabTable[DESB] = "DESB";

  LabTable[CREA_CREA] = "CREA_CREA";
  LabTable[CREA_CREB] = "CREA_CREB";
  LabTable[CREB_CREA] = "CREB_CREA";
  LabTable[CREB_CREB] = "CREB_CREB";
  LabTable[CREA_DESA] = "CREA_DESA";
  LabTable[CREA_DESB] = "CREA_DESB";
  LabTable[CREB_DESA] = "CREB_DESA";
  LabTable[CREB_DESB] = "CREB_DESB";
  LabTable[DESA_CREA] = "DESA_CREA";
  LabTable[DESA_CREB] = "DESA_CREB";
  LabTable[DESB_CREA] = "DESB_CREA";
  LabTable[DESB_CREB] = "DESB_CREB";
  LabTable[DESA_DESA] = "DESA_DESA";
  LabTable[DESA_DESB] = "DESA_DESB";
  LabTable[DESB_DESA] = "DESB_DESA";
  LabTable[DESB_DESB] = "DESB_DESB";

  LabTable[H] = "H";

  LabTable[CREA_COMP] = "CREA_COMP";
  LabTable[CREB_COMP] = "CREB_COMP";
  LabTable[DESA_COMP] = "DESA_COMP";
  LabTable[DESB_COMP] = "DESB_COMP";

  // set ConjTable

  ConjTable[I] = H;

  ConjTable[CREA] = DESA_COMP;
  ConjTable[CREB] = DESB_COMP;
  ConjTable[DESA] = CREA_COMP;
  ConjTable[DESB] = CREB_COMP;

  ConjTable[CREA_CREA] = DESA_DESA;
  ConjTable[CREA_CREB] = DESB_DESA;
  ConjTable[CREB_CREA] = DESA_DESB;
  ConjTable[CREB_CREB] = DESB_DESB;
  ConjTable[CREA_DESA] = CREA_DESA;
  ConjTable[CREA_DESB] = CREB_DESA;
  ConjTable[CREB_DESA] = CREA_DESB;
  ConjTable[CREB_DESB] = CREB_DESB;
  ConjTable[DESA_CREA] = DESA_CREA;
  ConjTable[DESA_CREB] = DESB_CREA;
  ConjTable[DESB_CREA] = DESA_CREB;
  ConjTable[DESB_CREB] = DESB_CREB;
  ConjTable[DESA_DESA] = CREA_CREA;
  ConjTable[DESA_DESB] = CREB_CREA;
  ConjTable[DESB_DESA] = CREA_CREB;
  ConjTable[DESB_DESB] = CREB_CREB;

  ConjTable[H] = I;

  ConjTable[CREA_COMP] = DESA;
  ConjTable[CREB_COMP] = DESB;
  ConjTable[DESA_COMP] = CREA;
  ConjTable[DESB_COMP] = CREB;

  // set MultTable

  for(int l = 0; l < SIZE; ++l)
    for(int s = 0; s < SIZE; ++s)
      MultTable[l][s] = NONE;

  // I(L) x O(S) -> O(LS)
  // O(L) x I(S) -> O(LS)
  for(int o = 0; o < SIZE; ++o) {
    MultTable[I][o] = static_cast<CATEGORY>(o);
    MultTable[o][I] = static_cast<CATEGORY>(o);
  }

  // Ci(L) x Cj(S) -> Cij(LS)
  MultTable[CREA][CREA] = CREA_CREA;
  MultTable[CREA][CREB] = CREA_CREB;
  MultTable[CREB][CREA] = CREB_CREA;
  MultTable[CREB][CREB] = CREB_CREB;
  MultTable[CREA][DESA] = CREA_DESA;
  MultTable[CREA][DESB] = CREA_DESB;
  MultTable[CREB][DESA] = CREB_DESA;
  MultTable[CREB][DESB] = CREB_DESB;
  MultTable[DESA][CREA] = DESA_CREA;
  MultTable[DESA][CREB] = DESA_CREB;
  MultTable[DESB][CREA] = DESB_CREA;
  MultTable[DESB][CREB] = DESB_CREB;
  MultTable[DESA][DESA] = DESA_DESA;
  MultTable[DESA][DESB] = DESA_DESB;
  MultTable[DESB][DESA] = DESB_DESA;
  MultTable[DESB][DESB] = DESB_DESB;

  // Ci(L) x Si(S) -> H(LS)
  MultTable[CREA][DESA_COMP] = H;
  MultTable[CREB][DESB_COMP] = H;
  MultTable[DESA][CREA_COMP] = H;
  MultTable[DESB][CREB_COMP] = H;

  // Si(L) x Ci(S) -> H(LS)
  MultTable[CREA_COMP][DESA] = H;
  MultTable[CREB_COMP][DESB] = H;
  MultTable[DESA_COMP][CREA] = H;
  MultTable[DESB_COMP][CREB] = H;

  // Ci(L) x Cjk(S) -> Sl(LS)
  MultTable[CREA][CREA_DESA] = CREA_COMP;
  MultTable[CREA][CREB_DESA] = CREB_COMP;
  MultTable[CREA][CREB_DESB] = CREA_COMP;
  MultTable[CREA][DESA_CREA] = CREA_COMP;
  MultTable[CREA][DESA_CREB] = CREB_COMP;
  MultTable[CREA][DESB_CREB] = CREA_COMP;
  MultTable[CREA][DESA_DESA] = DESA_COMP;
  MultTable[CREA][DESA_DESB] = DESB_COMP;
  MultTable[CREA][DESB_DESA] = DESB_COMP;

  MultTable[CREB][CREA_DESA] = CREB_COMP;
  MultTable[CREB][CREA_DESB] = CREA_COMP;
  MultTable[CREB][CREB_DESB] = CREB_COMP;
  MultTable[CREB][DESA_CREA] = CREB_COMP;
  MultTable[CREB][DESB_CREA] = CREA_COMP;
  MultTable[CREB][DESB_CREB] = CREB_COMP;
  MultTable[CREB][DESA_DESB] = DESA_COMP;
  MultTable[CREB][DESB_DESA] = DESA_COMP;
  MultTable[CREB][DESB_DESB] = DESB_COMP;

  MultTable[DESA][CREA_CREA] = CREA_COMP;
  MultTable[DESA][CREA_CREB] = CREB_COMP;
  MultTable[DESA][CREB_CREA] = CREB_COMP;
  MultTable[DESA][CREA_DESA] = DESA_COMP;
  MultTable[DESA][CREA_DESB] = DESB_COMP;
  MultTable[DESA][CREB_DESB] = DESA_COMP;
  MultTable[DESA][DESA_CREA] = DESA_COMP;
  MultTable[DESA][DESB_CREA] = DESB_COMP;
  MultTable[DESA][DESB_CREB] = DESA_COMP;

  MultTable[DESB][CREA_CREB] = CREA_COMP;
  MultTable[DESB][CREB_CREA] = CREA_COMP;
  MultTable[DESB][CREB_CREB] = CREB_COMP;
  MultTable[DESB][CREA_DESA] = DESB_COMP;
  MultTable[DESB][CREB_DESA] = DESA_COMP;
  MultTable[DESB][CREB_DESB] = DESB_COMP;
  MultTable[DESB][DESA_CREA] = DESB_COMP;
  MultTable[DESB][DESA_CREB] = DESA_COMP;
  MultTable[DESB][DESB_CREB] = DESB_COMP;

  // Cij(L) x Ck(S) -> Sl(LS)
  MultTable[CREA_DESA][CREA] = CREA_COMP;
  MultTable[CREB_DESA][CREA] = CREB_COMP;
  MultTable[CREB_DESB][CREA] = CREA_COMP;
  MultTable[DESA_CREA][CREA] = CREA_COMP;
  MultTable[DESA_CREB][CREA] = CREB_COMP;
  MultTable[DESB_CREB][CREA] = CREA_COMP;
  MultTable[DESA_DESA][CREA] = DESA_COMP;
  MultTable[DESA_DESB][CREA] = DESB_COMP;
  MultTable[DESB_DESA][CREA] = DESB_COMP;

  MultTable[CREA_DESA][CREB] = CREB_COMP;
  MultTable[CREA_DESB][CREB] = CREA_COMP;
  MultTable[CREB_DESB][CREB] = CREB_COMP;
  MultTable[DESA_CREA][CREB] = CREB_COMP;
  MultTable[DESB_CREA][CREB] = CREA_COMP;
  MultTable[DESB_CREB][CREB] = CREB_COMP;
  MultTable[DESA_DESB][CREB] = DESA_COMP;
  MultTable[DESB_DESA][CREB] = DESA_COMP;
  MultTable[DESB_DESB][CREB] = DESB_COMP;

  MultTable[CREA_CREA][DESA] = CREA_COMP;
  MultTable[CREA_CREB][DESA] = CREB_COMP;
  MultTable[CREB_CREA][DESA] = CREB_COMP;
  MultTable[CREA_DESA][DESA] = DESA_COMP;
  MultTable[CREA_DESB][DESA] = DESB_COMP;
  MultTable[CREB_DESB][DESA] = DESA_COMP;
  MultTable[DESA_CREA][DESA] = DESA_COMP;
  MultTable[DESB_CREA][DESA] = DESB_COMP;
  MultTable[DESB_CREB][DESA] = DESA_COMP;

  MultTable[CREA_CREB][DESB] = CREA_COMP;
  MultTable[CREB_CREA][DESB] = CREA_COMP;
  MultTable[CREB_CREB][DESB] = CREB_COMP;
  MultTable[CREA_DESA][DESB] = DESB_COMP;
  MultTable[CREB_DESA][DESB] = DESA_COMP;
  MultTable[CREB_DESB][DESB] = DESB_COMP;
  MultTable[DESA_CREA][DESB] = DESB_COMP;
  MultTable[DESA_CREB][DESB] = DESA_COMP;
  MultTable[DESB_CREB][DESB] = DESB_COMP;

  // Cij(L) x Ckl(S) -> H
  MultTable[CREA_CREA][DESA_DESA] = H;
  MultTable[CREA_CREB][DESA_DESB] = H;
  MultTable[CREA_CREB][DESB_DESA] = H;
  MultTable[CREB_CREA][DESA_DESB] = H;
  MultTable[CREB_CREA][DESB_DESA] = H;
  MultTable[CREB_CREB][DESB_DESB] = H;

  MultTable[CREA_DESA][CREA_DESA] = H;
  MultTable[CREA_DESB][CREB_DESA] = H;
  MultTable[CREB_DESA][CREA_DESB] = H;
  MultTable[CREB_DESB][CREB_DESB] = H;
  MultTable[CREA_DESA][DESA_CREA] = H;
  MultTable[CREA_DESB][DESA_CREB] = H;
  MultTable[CREB_DESA][DESB_CREA] = H;
  MultTable[CREB_DESB][DESB_CREB] = H;

  MultTable[DESA_CREA][DESA_CREA] = H;
  MultTable[DESA_CREB][DESB_CREA] = H;
  MultTable[DESB_CREA][DESA_CREB] = H;
  MultTable[DESB_CREB][DESB_CREB] = H;
  MultTable[DESA_CREA][CREA_DESA] = H;
  MultTable[DESA_CREB][CREA_DESB] = H;
  MultTable[DESB_CREA][CREB_DESA] = H;
  MultTable[DESB_CREB][CREB_DESB] = H;

  MultTable[DESA_DESA][CREA_CREA] = H;
  MultTable[DESA_DESB][CREA_CREB] = H;
  MultTable[DESA_DESB][CREB_CREA] = H;
  MultTable[DESB_DESA][CREA_CREB] = H;
  MultTable[DESB_DESA][CREB_CREA] = H;
  MultTable[DESB_DESB][CREB_CREB] = H;
}
