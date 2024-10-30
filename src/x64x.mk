SHELL=/bin/bash
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ifndef FP
FP=precise
endif # !FP
RM=rm -rfv
AR=ar
ARFLAGS=rsv
CC=icx
FC=ifx
ifndef MARCH
MARCH=Host
# common-avx512 for KNLs
endif # !MARCH
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -mprefer-vector-width=512 -vec-threshold0 -qopenmp -x$(MARCH)
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -threads
C18FLAGS=$(CPUFLAGS) -std=c18
FPUFLAGS=-fp-model=$(FP) -fp-speculation=safe -fma -fprotect-parens -no-ftz -fimf-precision=high
ifneq ($(FP),strict)
FPUFLAGS += -fimf-use-svml=true
endif # !strict
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
ifeq ($(FP),strict)
FPUFFLAGS += -assume ieee_fpe_flags
endif # strict
DBGFLAGS=-traceback
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -fno-math-errno -inline-level=2
OPTFFLAGS=$(OPTFLAGS) -DMKL_DIRECT_CALL
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS += -DNDEBUG -qopt-report=3
DBGFFLAGS=$(DBGFLAGS)
DBGCFLAGS=$(DBGFLAGS)
else # DEBUG
OPTFLAGS=-O0
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug parallel -debug pubnames
DBGFFLAGS=$(DBGFLAGS) -debug-parameters all -check all -warn all
DBGCFLAGS=$(DBGFLAGS)
endif # ?NDEBUG
LIBFLAGS=-D_GNU_SOURCE -DUSE_MKL -DMKL_ILP64 -I. -I../../../JACSD/vn -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
LDFLAGS=-rdynamic -static-libgcc -L../../../JACSD -lvn$(DEBUG) -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -ldl
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
