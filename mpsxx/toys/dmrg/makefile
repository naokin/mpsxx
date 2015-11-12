BOOSTDIR = /home/naokin/boost/1.52.0
BOOSTINCLUDE = -I$(BOOSTDIR)/include
BOOSTLIB = -L$(BOOSTDIR)/lib -lboost_serialization

BLITZDIR = /home/naokin/blitz/0.10
BLITZINCLUDE = -I$(BLITZDIR)/include
BLITZLIB = -L$(BLITZDIR)/lib -lblitz

BTASDIR = /home/naokin/btas-0.2.1
BTASINCLUDE = -I$(BTASDIR) -I$(BTASDIR)/btas -DMKL_CBLAS


MKLDIR = /opt/intel/mkl/10.2.4.032
MKLINCLUDE = -I$(MKLDIR)/include
#MKLLIB = -L$(MKLDIR)/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
MKLLIB = -L$(MKLDIR)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide -lpthread

CXX = g++44 -O3

FLAGS = $(BOOSTINCLUDE) $(BLITZINCLUDE) $(BTASINCLUDE) $(MKLINCLUDE)
LIBS  = $(BOOSTLIB) $(BLITZLIB) $(MKLLIB)

.C.o	:
	$(CXX) $(FLAGS) -c $< -o $@

all	: dmrg.x

dmrg.x	: main.o dmrg.o driver.o driverlrt.o input.o
	$(CXX) $(FLAGS) -o dmrg.x main.o dmrg.o driver.o driverlrt.o input.o $(LIBS) $(BTASDIR)/src/libbtas.a

clean:
	rm *.o *.x

#
