CXX=/homec/naokin/gnu/gcc/4.8.2/bin/g++ -std=c++11
CXXFLAGS=-g -O3 -openmp -D_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=
BLASINC=-I/home100/opt/intel/mkl/include
BLASLIB=-L/home100/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential

BOOSTINC=-I/homec/naokin/boost/1.54.0/include
BOOSTLIB=-L/homec/naokin/boost/1.54.0/lib -lboost_serialization

BTASDIR= $(HOME)/devel/btas
BTASINC=-I. -I.. -I$(BTASDIR)/include
BTASLIB= $(BTASDIR)/lib/libbtas.a

INCLUDEFLAGS= $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS= $(BLASLIB) $(BOOSTLIB)

SRC_DMRG = dmrg.C main.C ../driver/fileio.C

OBJ_DMRG = $(SRC_DMRG:.C=.o)

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: dmrg

dmrg	: $(OBJ_DMRG)
	$(CXX) $(CXXFLAGS) -o dmrg.x $(OBJ_DMRG) $(BTASLIB) $(LIBRARYFLAGS)

clean	:
	rm *.o; rm *.x; rm *.a; rm ../driver/*.o;
