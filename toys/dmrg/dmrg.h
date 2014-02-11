#ifndef PROTOTYPE_DMRG_H
#define PROTOTYPE_DMRG_H

#include <iostream>
#include <vector>
#include <btas/DArray.h>
#include "input.h"

//
// MpStates: Mps ( l, n, r )
//
typedef std::vector< btas::DArray<3> > MpStates;
//
// MpOprtrs: Mpo ( l, n', n, r )
//
typedef std::vector< btas::DArray<4> > MpOprtrs;
//
// Storages: Str ( l', l ) ( or, Str ( r', r ) )
//
typedef std::vector< btas::DArray<3> > Storages;
//
// DnWeight: Dnw ( l ) ( or, Dnw ( r ) )
//
typedef std::vector< btas::DArray<1> > DnWeight;
//
// Overlaps: Ovs ( l', l ) ( or, Ovs ( r', r ) )
//
typedef std::vector< btas::DArray<2> > Overlaps;


void mpo_init(const DmrgInput& input, MpOprtrs& mpos);

void dmrg   (std::ostream& fout, const DmrgInput& input, const MpOprtrs& mpos, MpStates& wfns,
             DnWeight& lval, MpStates& lmps, MpStates& lnul, Storages& lstr,
             DnWeight& rval, MpStates& rmps, MpStates& rnul, Storages& rstr);

void dmrg   (std::ostream& fout, const DmrgInput& input, const MpOprtrs& mpos, MpStates& wfns,
             DnWeight& lval, MpStates& lmps, MpStates& lnul, Storages& lstr,
             DnWeight& rval, MpStates& rmps, MpStates& rnul, Storages& rstr,
             const std::vector<MpStates>& wfns0, const std::vector<MpStates>& lmps0, const std::vector<MpStates>& rmps0);

void dmrglrt(std::ostream& fout, const DmrgInput& input, const MpOprtrs& mpos, MpStates& wfns,
             DnWeight& lval, MpStates& lmps, MpStates& lnul, Storages& lstr,
             DnWeight& rval, MpStates& rmps, MpStates& rnul, Storages& rstr, int nroot = 3);

#endif // PROTOTYPE_DMRG_H
