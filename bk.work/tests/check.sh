#!/bin/bash

g++ intgen.cpp -o intgen.x
./intgen.x <<< 4

/home100/naokin/Block-smp/block.spin_adapted block.conf >& block.out
mpogen -f FCIDUMP >& log
rundmrg -i dmrg.conf -o dmrg.out

#
