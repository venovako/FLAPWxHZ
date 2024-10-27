SHELL=/bin/bash
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
AR=ar
ARFLAGS=rsv
ifndef MARCH
ifeq ($(shell uname -m),ppc64le)
MARCH=mcpu=native -mpower8-fusion -mtraceback=full
else # !ppc64le
MARCH=march=native
endif # ?ppc64le
endif # !MARCH
CPUFLAGS=-DUSE_GNU -DUSE_X64 -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -fvect-cost-model=unlimited -$(MARCH) -fopenmp
FORFLAGS=-cpp $(CPUFLAGS) -fdefault-integer-8 -ffree-line-length-none -fstack-arrays
C18FLAGS=$(CPUFLAGS) -std=gnu18
ifeq ($(ARCH),Darwin)
ifndef GNU
GNU=-8
endif # !GNU
endif # Darwin
CC=gcc$(GNU)
FC=gfortran$(GNU)
DBGFLAGS=-Wall -Wextra
FPUFLAGS=-ffp-contract=fast
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS += -DNDEBUG
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # Darwin
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFFLAGS=$(DBGFLAGS) -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
DBGCFLAGS=$(DBGFLAGS)
FPUFLAGS += -fno-math-errno
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O$(DEBUG)
DBGFLAGS += -$(DEBUG)
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # ?Darwin
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFFLAGS=$(DBGFLAGS) -fcheck=array-temps -finit-local-zero -finit-real=snan -finit-derived -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all #-fcheck=all
DBGCFLAGS=$(DBGFLAGS) #-ftrapv
FPUFFLAGS=$(FPUFLAGS) #-ffpe-trap=invalid,zero,overflow
FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-DMKL_ILP64 -I. -I../../../JACSD/vn
ifeq ($(ARCH),Linux)
LIBFLAGS += -D_GNU_SOURCE
endif # Linux
LDFLAGS=-rdynamic -static-libgcc -static-libgfortran -static-libquadmath -L../../../JACSD -lvn$(DEBUG)
ifdef MKLROOT
LIBFLAGS += -DUSE_MKL -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
ifdef NDEBUG
LIBFLAGS += -DMKL_DIRECT_CALL
endif # NDEBUG
ifeq ($(ARCH),Darwin)
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../../compiler/latest/mac/compiler/lib -Wl,-rpath,${MKLROOT}/../../compiler/latest/mac/compiler/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lgomp
else # Linux
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp
endif # ?Darwin
else # !MKLROOT
ifndef LAPACK
LAPACK=$(HOME)/lapack-ilp64
endif # !LAPACK
LDFLAGS += -L$(LAPACK) -ltmglib -llapack -lrefblas
endif # ?MKLROOT
LDFLAGS += -lpthread -lm -ldl $(shell if [ -L /usr/lib64/libmemkind.so ]; then echo '-lmemkind'; fi)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
