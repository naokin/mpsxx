#ifndef __MPSXX_MPO_BASE_OPERATOR_HPP
#define __MPSXX_MPO_BASE_OPERATOR_HPP

namespace mpsxx {

namespace op {

enum CATEGORY
{
      H = 0x100000,
      I = 0x000000,

      CRE = 0x010001,
      DES = 0x010000,

      CRE_COMP = 0x110001,
      DES_COMP = 0x110000,

      CRE_CRE = 0x020011,
      CRE_DES = 0x020010,
      DES_CRE = 0x020001,
      DES_DES = 0x020000,

      CRE_CRE_COMP = 0x120011,
      CRE_DES_COMP = 0x120010,
      DES_CRE_COMP = 0x120001,
      DES_DES_COMP = 0x120000,

      CRE_CRE_DES = 0x030110,
      CRE_DES_CRE = 0x030101,
      DES_CRE_CRE = 0x030011,

      CRE_DES_DES = 0x030100,
      DES_CRE_DES = 0x030010,
      DES_DES_CRE = 0x030001,

      CRE_CRE_DES_DES = 0x041100,
      CRE_DES_CRE_DES = 0x041010,
      CRE_DES_DES_CRE = 0x041001,

      DES_DES_CRE_CRE = 0x040011,
      DES_CRE_DES_CRE = 0x040101,
      DES_CRE_CRE_DES = 0x040110
};

CATEGORY multiply (CATEGORY c1, CATEGORY c2)
{
   // if one is H, other must be I
   if(c1 == H) assert(c2 == I);
   if(c2 == H) assert(c1 == I);

   // comp. op. not to be multiplied by comp. op.
   assert(!(c1 & c2 & 0xf00000));

   n1 = (c1 & 0x0f0000) >> 16;
   n2 = (c2 & 0x0f0000) >> 16;
   n3 = n1 + n2;

   // only allow up to 4-index operator
   assert(n3 <= 4);

   CATEGORY c3 = ((c1 | c2) & 0xf00000) | (n3 << 16) | (((c1 & 0x00ffff) << (n2*4)) & 0x00ffff) | (c2 & 0x00ffff);

   return c3;
}

enum SPINCASE
{
      N = 0x000000,

      A = 0x010001,
      B = 0x010000,

      AA = 0x020011,
      AB = 0x020010,
      BA = 0x020001,
      BB = 0x020000,

      AAA = 0x030111,
      AAB = 0x030110,
      ABA = 0x030101,
      BAA = 0x030011,
      ABB = 0x030100,
      BAB = 0x030010,
      BBA = 0x030001,
      BBB = 0x030000,

      AAAA = 0x041111,
      AABB = 0x041100,
      ABAB = 0x041010,
      ABBA = 0x041001,
      BBAA = 0x040011,
      BABA = 0x040101,
      BAAB = 0x040110,
      BBBB = 0x040000
};

SPINCASE multiply (SPINCASE s1, SPINCASE s2)
{
   n1 = (s1 & 0x0f0000) >> 16;
   n2 = (s2 & 0x0f0000) >> 16;
   n3 = n1 + n2;

   // only allow up to 4-index operator
   assert(n3 <= 4);

   SPINCASE s3 = (n3 << 16) | (((s1 & 0x00ffff) << (n2*4)) & 0x00ffff) | (s2 & 0x00ffff);

   return s3;
}

} // namespace op

struct BaseOperator
{
   op::CATEGORY category;

   op::SPINCASE spincase;

   long index[4];

   explicit
   BaseOperator (op::CATEGORY c = op::I, op::SPINCASE s = op::N, long i = -1, long j = -1, long k = -1, long l = -1)
   {
      long ct = (c & 0x0f0000) >> 16;
      long st = (s & 0x0f0000) >> 16;

      assert(ct == st);

      index[0] = i;
      index[1] = j;
   }
};

} // namespace mpsxx

inline bool operator== (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
   if(x.category != y.category)
      return false;
   if(x.spincase != y.spincase)
      return false;

   return (x.index[0] == y.index[0] && x.index[1] == y.index[1]);
}

inline bool operator!= (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
   if(x.category != y.category)
      return true;
   if(x.spincase != y.spincase)
      return true;

   return (x.index[0] != y.index[0] || x.index[1] != y.index[1]);
}

inline bool operator< (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
   if(x.category != y.category)
      return x.category < y.category;
   if(x.spincase != y.spincase)
      return x.spincase < y.spincase;

   return (x.index[0] != y.index[0]) ? (x.index[0] < y.index[0]) : (x.index[1] < y.index[1]);
}

inline bool operator<= (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
   if(x.category != y.category)
      return x.category < y.category;
   if(x.spincase != y.spincase)
      return x.spincase < y.spincase;

   return (x.index[0] != y.index[0]) ? (x.index[0] < y.index[0]) : (x.index[1] <= y.index[1]);
}

inline bool operator> (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
   if(x.category != y.category)
      return x.category > y.category;
   if(x.spincase != y.spincase)
      return x.spincase > y.spincase;

   return (x.index[0] != y.index[0]) ? (x.index[0] > y.index[0]) : (x.index[1] > y.index[1]);
}

inline bool operator>= (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
   if(x.category != y.category)
      return x.category > y.category;
   if(x.spincase != y.spincase)
      return x.spincase > y.spincase;

   return (x.index[0] != y.index[0]) ? (x.index[0] > y.index[0]) : (x.index[1] >= y.index[1]);
}

inline mpsxx::BaseOperator operator* (const mpsxx::BaseOperator& x, const mpsxx::BaseOperator& y)
{
}

#endif // __MPSXX_MPO_BASE_OPERATOR_HPP
