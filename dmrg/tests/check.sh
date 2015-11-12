#!/bin/bash

g++ intgen.cpp -o intgen.x
./intgen.x <<< 6

/home100/naokin/Block-smp/block.spin_adapted block.conf.0 >& block.out.0
/home100/naokin/Block-smp/block.spin_adapted block.conf.2 >& block.out.2
#/home100/naokin/Block-smp/block.spin_adapted block.conf.4 >& block.out.4

mpogen -f FCIDUMP -x >& log
rundmrg -i dmrg.conf -o dmrg.out

grep "Sweep Energy" block.out.0 block.out.2 dmrg.out

#
