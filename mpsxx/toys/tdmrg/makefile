BOOSTDIR = /home/naokin/boost/1.52.0
BOOSTINCLUDE = -I$(BOOSTDIR)/include
BOOSTLIB = -L$(BOOSTDIR)/lib -lboost_serialization

BLITZDIR = /home/naokin/blitz/0.9
BLITZINCLUDE = -I$(BLITZDIR)/include
BLITZLIB = -L$(BLITZDIR)/lib -lblitz

BTASDIR = /home/naokin/btas-0.2.2
BTASINCLUDE = -I$(BTASDIR)


MKLDIR = /opt/intel/mkl/10.2.4.032
MKLINCLUDE = -I$(MKLDIR)/include
MKLLIB = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#MKLLIB = -L$(MKLDIR)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread

#CXX = /home/naokin/gnu/gcc/4.7.1/bin/g++ -O2 -std=c++0x
CXX = g++44

FLAGS = $(BOOSTINCLUDE) $(BLITZINCLUDE) $(BTASINCLUDE) $(MKLINCLUDE)
LIBS  = $(BOOSTLIB) $(BLITZLIB) $(MKLLIB)

.C.o :
	$(CXX) $(FLAGS) -c $< -o $@

all	: tdmrg.x

tdmrg.x : tdmrg.o dmrg.o real_tdmrg.o imag_tdmrg.o
	$(CXX) $(FLAGS) -o tdmrg.x tdmrg.o dmrg.o real_tdmrg.o imag_tdmrg.o $(LIBS) $(BTASDIR)/libbtas.a

clean:
	rm *.o *.x

#
