#ifndef __MPSXX_MPOGEN_CATEGORY_H
#define __MPSXX_MPOGEN_CATEGORY_H

#include <string>
#include <symmetry/fermion.h>

namespace mpsxx {
namespace op {

/// For quantum chemistry Hamiltonian
namespace QC {

enum CATEGORY {

  I = 0u,

  CREA,
  CREB,
  DESA,
  DESB,

  CREA_CREA,
  CREA_CREB,
  CREB_CREA,
  CREB_CREB,
  CREA_DESA,
  CREA_DESB,
  CREB_DESA,
  CREB_DESB,
  DESA_CREA,
  DESA_CREB,
  DESB_CREA,
  DESB_CREB,
  DESA_DESA,
  DESA_DESB,
  DESB_DESA,
  DESB_DESB,

  H,

  CREA_COMP,
  CREB_COMP,
  DESA_COMP,
  DESB_COMP,

  SIZE,

  CREA_CREB_DESB,
  CREB_CREA_DESA,
  CREA_DESA_DESB,
  CREB_DESB_DESA,

  CREA_DESA_CREB_DESB,

  NONE = 0xffffffff
};

extern fermion QNumTable[SIZE];

extern std::string LabTable[SIZE];

extern CATEGORY ConjTable[SIZE];

extern CATEGORY MultTable[SIZE][SIZE];

void setTables ();

} // namespace QC

/// For RDM calculation (has not yet implemented).
namespace RDM {

typedef unsigned int CATEGORY;

static const CATEGORY CREA = 0x00000000; // 00
static const CATEGORY CREB = 0x00000001; // 01
static const CATEGORY DESA = 0x00000002; // 10
static const CATEGORY DESB = 0x00000003; // 11

static const CATEGORY MAX_ORDER = 8;

static const CATEGORY NONE = 0xffffffff;

static const CATEGORY SIZE = (1u << (2u*MAX_ORDER)); // = 4^MAX_ORDER

} // namespace RDM

extern RDM::CATEGORY QC2RDM[QC::SIZE];

} // namespace op
} // namespace mpsxx

#endif // __MPSXX_MPOGEN_CATEGORY_H
