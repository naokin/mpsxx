#!/bin/bash

#generate an ordering
#../genetic/gaopt -config ../genetic/ga.conf -integral FCIDUMP

# generate QC-MPO and store on disk
../../mpogen/mpogen.x -f FCIDUMP -r reorder.dat

# DMRG opt. with 2-site algo'
../dmrg.x -i twosite.conf
# DMRG opt. with 1-site algo' restarting from 2-site calc.
../dmrg.x -i onesite.conf

#remove temporary files
rm *.tmp
