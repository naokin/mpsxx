#ifndef __MPSXX_MPSITE_HPP
#define __MPSXX_MPSITE_HPP 1

// TODO: Test for spin-1/2 system

namespace mpsxx {

class MpSite {

   unsigned long size_;

   std::vector<unsigned long> spin_up_;

   std::vector<unsigned long> spin_down_;

};

} // namespace mpsxx

#endif // __MPSXX_MPSITE_HPP
