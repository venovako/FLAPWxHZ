SHELL=/bin/bash
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ifndef FP
ifdef NDEBUG
FP=precise
else # DEBUG
FP=strict
endif # ?NDEBUG
endif # !FP
RM=rm -rfv
AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=icc
FC=ifort
ifndef MARCH
MARCH=Host
# COMMON-AVX512 for KNLs
endif # !MARCH
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -qopt-multi-version-aggressive -qopt-zmm-usage=high -vec-threshold0 -qopenmp -x$(MARCH)
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -threads
C18FLAGS=$(CPUFLAGS) -std=gnu18
FPUFLAGS=-fp-model $(FP) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fno-math-errno
ifeq ($(FP),strict)
FPUFLAGS += -fp-stack-check
else # !strict
FPUFLAGS += -fimf-use-svml=true
endif # ?strict
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
ifeq ($(FP),strict)
FPUFFLAGS += -assume ieee_fpe_flags
endif # strict
DBGFLAGS=-traceback -diag-disable=10397,10441
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -inline-level=2
OPTFFLAGS=$(OPTFLAGS) -DMKL_DIRECT_CALL
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS += -DNDEBUG -qopt-report=5
DBGFFLAGS=$(DBGFLAGS)
DBGCFLAGS=$(DBGFLAGS) -w3 -diag-disable=1572,2547
else # DEBUG
OPTFLAGS=-O0
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames
ifneq ($(ARCH),Darwin)
DBGFLAGS += -debug parallel
endif # Linux
DBGFFLAGS=$(DBGFLAGS) -debug-parameters all -check all -warn all
DBGCFLAGS=$(DBGFLAGS) -check=stack,uninit -w3 -diag-disable=1572,2547,10441
endif # ?NDEBUG
LIBFLAGS=-DUSE_MKL -DMKL_ILP64 -I. -I../../../JACSD/vn -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
LDFLAGS=-rdynamic -L../../../JACSD -lvn$(DEBUG)
ifeq ($(ARCH),Darwin)
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
else # Linux
LIBFLAGS += -D_GNU_SOURCE
LDFLAGS += -static-libgcc -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
endif # ?Darwin
LDFLAGS += -lpthread -lm -ldl
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
