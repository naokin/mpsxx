CXX=g++ -std=c++11
CXXFLAGS=-g -O3 -openmp -D_DEFAULT_QUANTUM -D_HAS_CBLAS -D_HAS_INTEL_MKL

BLASDIR=
BLASINC=
BLASLIB=-lmkl_intel_lp64 -lmkl_core -lmkl_sequential
#BLASLIB=-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lguide -lpthread

BOOSTLIB= -lboost_serialization

BTASDIR=/home/bright/bestanden/programmas/btas
BTASINC=-I. -I.. -I$(BTASDIR)/include
BTASLIB= $(BTASDIR)/lib/libbtas.a

INCLUDEFLAGS= $(BLASINC) $(BOOSTINC) $(BTASINC)
LIBRARYFLAGS= $(BLASLIB) $(BOOSTLIB)

SRC_DMRG = dmrg.C main.C btas_template_specialize.C ../driver/fileio.C

OBJ_DMRG = $(SRC_DMRG:.C=.o)

.C.o	:
	$(CXX) $(CXXFLAGS) $(INCLUDEFLAGS) -c $< -o $@

all	: dmrg

dmrg	: $(OBJ_DMRG)
	$(CXX) $(CXXFLAGS) -o dmrg.x $(OBJ_DMRG) $(BTASLIB) $(LIBRARYFLAGS)

clean	:
	rm *.o; rm *.x; rm *.a;
