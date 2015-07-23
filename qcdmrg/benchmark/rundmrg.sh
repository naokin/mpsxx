#$ -cwd
#$ -pe smp 6
#$ -S /bin/bash
#$ -N DMRG

if [ ! -e /work/naokin/qcdmrg ]; then
mkdir -p /work/naokin/qcdmrg
fi

#generate an ordering
#../genetic/gaopt -config ../genetic/ga.conf -integral FCIDUMP

OUT=dmrg.log

# generate QC-MPO and store on disk
../../mpogen/mpogen.x -f FCIDUMP -r reorder.dat -s /work/naokin/qcdmrg/ >& $OUT
#../../mpogen/naive_mpogen.x -f FCIDUMP -r reorder.dat -s /work/naokin/qcdmrg/ >& $OUT

# DMRG opt. with 2-site algo'
../dmrg.x -i twosite.conf -s /work/naokin/qcdmrg/ 1>> $OUT 2>> $OUT
# DMRG opt. with 1-site algo' restarting from 2-site calc.
../dmrg.x -i onesite.conf -s /work/naokin/qcdmrg/ 1>> $OUT 2>> $OUT

#remove temporary files
rm -f /work/naokin/qcdmrg/*.tmp
mv /work/naokin/qcdmrg/* .
rmdir /work/naokin/qcdmrg
