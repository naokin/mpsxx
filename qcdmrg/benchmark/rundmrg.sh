#$ -cwd
#$ -pe smp 6
#$ -S /bin/bash
#$ -N DMRG

if [ ! -e /work/naokin/qcdmrg ]; then
mkdir -p /work/naokin/qcdmrg
fi

#generate an ordering
#../genetic/gaopt -config ../genetic/ga.conf -integral FCIDUMP

if [ ! -d state ]; then
mkdir state
fi

if [ ! -d hamiltonian ]; then
mkdir hamiltonian
fi

# generate QC-MPO and store on disk
../../mpogen/mpogen.x -f FCIDUMP -r reorder.dat -s /work/naokin/qcdmrg/ >& dmrg.log

# DMRG opt. with 2-site algo'
../dmrg.x -i twosite.conf -s /work/naokin/qcdmrg/ >& dmrg.log.2site
# DMRG opt. with 1-site algo' restarting from 2-site calc.
../dmrg.x -i onesite.conf -s /work/naokin/qcdmrg/ >& dmrg.log.1site

#remove temporary files
#rm -f /work/naokin/qcdmrg/*.tmp

#mv /work/naokin/qcdmrg/* .
#rmdir /work/naokin/qcdmrg
