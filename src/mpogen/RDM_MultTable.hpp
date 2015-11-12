#include <iostream>
#include <iomanip>
#include <string>

#include "Category.h"

namespace mpsxx {
namespace op {
namespace RDM {

/// Table of multiplied operator category
static unsigned int MultTable[SIZE][SIZE];

template<unsigned int ORDER1st, unsigned int ORDERMAX> struct setMultTable_recursive_helper_1;
template<unsigned int ORDER1st, unsigned int ORDER2nd> struct setMultTable_recursive_helper_2;

template<unsigned int ORDER1st, unsigned int ORDERMAX>
struct setMultTable_recursive_helper_1
{
  static void make ()
  {
    setMultTable_recursive_helper_2<ORDER1st,ORDERMAX-ORDER1st>::make();
    setMultTable_recursive_helper_1<ORDER1st+1u,ORDERMAX>::make();
  }
};

template<unsigned int ORDERMAX>
struct setMultTable_recursive_helper_1<ORDERMAX,ORDERMAX>
{
  static void make () { }
};

template<unsigned int ORDER1st, unsigned int ORDER2nd>
struct setMultTable_recursive_helper_2
{
  static void make ()
  {
    for(unsigned int op0 = CREA; op0 <= DESB; ++op0)
      setMultTable_recursive_helper_2<ORDER1st-1u,ORDER2nd>::make(op0);
  }

  static void make (const unsigned int& op0)
  {
    for(unsigned int op1 = CREA; op1 <= DESB; ++op1)
      setMultTable_recursive_helper_2<ORDER1st-1u,ORDER2nd>::make((op0 << 2u) | op1);
  }
};

template<unsigned int ORDER2nd>
struct setMultTable_recursive_helper_2<0u,ORDER2nd>
{
  static void make (const unsigned int& op0)
  {
    for(unsigned int op1 = CREA; op1 <= DESB; ++op1)
      setMultTable_recursive_helper_2<0u,ORDER2nd-1u>::make(op0, op1);
  }

  static void make (const unsigned int& op0, const unsigned int& op1)
  {
    for(unsigned int op2 = CREA; op2 <= DESB; ++op2)
      setMultTable_recursive_helper_2<0u,ORDER2nd-1u>::make(op0, (op1 << 2u) | op2);
  }
};

template<>
struct setMultTable_recursive_helper_2<0u,0u>
{
  static void make (const unsigned int& op0, const unsigned int& op1)
  {
    unsigned int opmod = op1;
    unsigned int shift;
    for(shift = 1u; shift < MAX_ORDER; ++shift) {
      opmod /= 4u; if(opmod == 0u) break;
    }
    MultTable[op0][op1] = ((op0 << (2u*shift)) | op1);
  }
};

template<unsigned int ORDER>
void setMultTable_recursive ()
{
  setMultTable_recursive<ORDER-1>();
  setMultTable_recursive_helper_1<1u,ORDER>::make();
}

template<>
void setMultTable_recursive<2u> ()
{
  setMultTable_recursive_helper_1<1u,2u>::make();
}

void setMultTable ()
{
  for(unsigned int i = 0u; i < SIZE; ++i)
    for(unsigned int j = 0u; j < SIZE; ++j)
      MultTable[i][j] = NONE;

  setMultTable_recursive<MAX_ORDER>();
}

std::string get_label (unsigned int op)
{
  static const char label_[][8u] = { "CREA", "CREB", "DESA", "DESB" };

  if(op == NONE) return std::string("NONE");

  std::string str(label_[op%4u]);
  op /= 4u;
  for(; op > 0u; op /= 4u) {
    std::string tmp(label_[op%4u]);
    str = tmp + "_" + str;
  }

  return str;
}

} // namespace RDM
} // namespace op
} // namespace mpsxx

int main ()
{
  op::setMultTable();

  for(unsigned int i = 1u; i < op::MAX_ORDER; ++i)
    for(unsigned int j = 1u; j <= op::MAX_ORDER-i; ++j)
      for(unsigned int iop = 0u; iop < (1u << (2u*i)); ++iop)
        for(unsigned int jop = 0u; jop < (1u << (2u*j)); ++jop) {
          std::cout << op::get_label(iop) << " x " << op::get_label(jop) << " = " << op::get_label(op::MultTable[iop][jop]) << std::endl;
        }

  return 0;
}
