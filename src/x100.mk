SHELL=/bin/bash
ARCH=$(shell uname)
RM=rm -rfv
AR=xiar
ifdef NDEBUG
ARFLAGS=-lib rsv
else # DEBUG
ARFLAGS=-qnoipo -lib rsv
endif # ?NDEBUG
FC=ifort
CC=icc
CPUFLAGS=-DUSE_INTEL -DUSE_X100 -mmic -fexceptions
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -threads #-DHAVE_IMAGINARY
C11FLAGS=$(CPUFLAGS) -DFORTRAN_INTEGER_KIND=8 -std=c11
ifdef NDEBUG
OPTFLAGS=-fast
OPTFFLAGS=$(OPTFLAGS) -DMKL_DIRECT_CALL
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-DNDEBUG -qopt-report=5 -traceback -diag-disable=10397
DBGFFLAGS=$(DBGFLAGS)
DBGCFLAGS=$(DBGFLAGS)-diag-disable=1572,2547
FPUFLAGS=-fma -fimf-precision=high -fimf-use-svml=true
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O0
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-g -debug emit_column -debug extended -debug inline-debug-info -debug parallel -debug pubnames -traceback -diag-disable=10397
DBGFFLAGS=$(DBGFLAGS) -debug-parameters all -check all -warn all
DBGCFLAGS=$(DBGFLAGS) -check=stack,uninit -w3 -diag-disable=1572,2547
FPUFLAGS=-fp-model strict -fp-stack-check -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-precision=high
FPUFFLAGS=$(FPUFLAGS) -assume ieee_fpe_flags
FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-DUSE_MKL -DMKL_ILP64 -I. -I${MKLROOT}/include/mic/ilp64 -I${MKLROOT}/include -qopenmp
LDFLAGS=-L${MKLROOT}/lib/mic -Wl,-rpath=${MKLROOT}/lib/mic -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
LDFLAGS += -lpthread -lm -ldl
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUCFLAGS)
