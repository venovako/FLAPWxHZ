SHELL=/bin/bash
ARCH=$(shell uname)
RM=rm -rfv
AR=ar
ARFLAGS=rsv
FC=gfortran
ifeq ($(ARCH),Darwin)
CC=clang
else # Linux
CC=gcc
endif # ?Darwin
CPUFLAGS=-DUSE_GNU -DUSE_X64
FORFLAGS=-cpp $(CPUFLAGS) -fdefault-integer-8 -ffree-line-length-none -fopenmp -fstack-arrays #-DHAVE_IMAGINARY
C11FLAGS=$(CPUFLAGS) -DFORTRAN_INTEGER_KIND=8
ifeq ($(ARCH),Darwin)
C11FLAGS += -std=gnu17 -pthread
else # Linux
C11FLAGS += -std=gnu11 -fopenmp
endif # ?Darwin
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native
DBGFLAGS=-DNDEBUG
ifeq ($(ARCH),Darwin)
OPTFFLAGS=$(OPTFLAGS) -Wa,-q -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller
OPTCFLAGS=$(OPTFLAGS) -integrated-as
DBGFFLAGS=$(DBGFLAGS) -fopt-info-optimized-vec
DBGCFLAGS=$(DBGFLAGS)
else # Linux
OPTFLAGS += -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS += -fopt-info-optimized-vec
DBGFFLAGS=$(DBGFLAGS) -pedantic -Wall -Wextra -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
DBGCFLAGS=$(DBGFLAGS)
endif # ?Darwin
FPUFLAGS=-ffp-contract=fast
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS) -fno-math-errno
OPTFFLAGS += -DMKL_DIRECT_CALL
ifdef NEST
OPTFFLAGS += -DMKL_DIRECT_CALL$(NEST)
endif # ?NEST
else # DEBUG
OPTFLAGS=-Og -march=native
ifeq ($(ARCH),Darwin)
OPTFFLAGS=$(OPTFLAGS) -Wa,-q
OPTCFLAGS=$(OPTFLAGS) -integrated-as
else # Linux
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
endif # ?Darwin
DBGFLAGS=-g
DBGFFLAGS=$(DBGFLAGS) -fcheck=all -finit-local-zero -finit-real=snan -finit-derived -pedantic -Wall -Wextra -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
DBGCFLAGS=$(DBGFLAGS) -ftrapv
FPUFLAGS=-ffp-contract=fast
FPUFFLAGS=$(FPUFLAGS) -ffpe-trap=invalid,zero,overflow
FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-DUSE_MKL -DMKL_ILP64 -DMKL_NEST$(NEST) -I. -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
ifeq ($(ARCH),Darwin)
ifeq ($(NEST),_SEQ)
LDFLAGS=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
else # parallel
LDFLAGS=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5
endif # ?NEST
else # Linux
ifeq ($(NEST),_SEQ)
LDFLAGS=-L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64 -Wl,-rpath=${MKLROOT}/../compiler/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core
else # parallel
LDFLAGS=-L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -L${MKLROOT}/../compiler/lib/intel64 -Wl,-rpath=${MKLROOT}/../compiler/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core
endif # ?NEST
endif # ?Darwin
LDFLAGS += -lpthread -lm -ldl
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUCFLAGS)
