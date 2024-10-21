SHELL=/bin/bash
ARCH=$(shell uname)
ifneq ($(ARCH),Linux)
$(error NVIDIA HPC SDK is supported on Linux only)
endif # !Linux
ifndef WP
WP=8
endif # !WP
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
AR=ar
ARFLAGS=rsv
ifndef MARCH
MARCH=native
endif # !MARCH
CPUFLAGS=-DUSE_NVIDIA -DUSE_X64 -DKIND_QUAD=$(WP) -m64 -mp -KPIC -Mframe -Meh_frame -Minfo -tp=$(MARCH) -nvmalloc -traceback
FORFLAGS=$(CPUFLAGS) -i8 -Mdclchk -Mlarge_arrays -Mrecursive -Mstack_arrays
C18FLAGS=$(CPUFLAGS) -c18
CC=nvc
FC=nvfortran
FPUFLAGS=-Kieee -Mfma -Mnodaz -Mnoflushz -Mnofpapprox -Mnofprelaxed -Mno-recip-div
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-DNDEBUG
DBGFFLAGS=$(DBGFLAGS)
DBGCFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O0
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-g -Mbounds -Mchkstk
DBGFFLAGS=$(DBGFLAGS)
DBGCFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-D_GNU_SOURCE -DMKL_ILP64 -I. -I../../../JACSD/vn
LDFLAGS=-Wl,-E -static-nvidia -L../../../JACSD -lvn$(DEBUG)
ifdef MKLROOT
LIBFLAGS += -DUSE_MKL -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
# define MKL=pgi_thread for a threaded MKL
ifndef MKL
MKL=sequential
endif # !MKL
LDFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_$(MKL) -lmkl_core
else # !MKLROOT
LDFLAGS += -L$(HOME)/lapack_ilp64 -ltmglib -llapack -lrefblas
endif # ?MKLROOT
LDFLAGS += -pgf90libs -lpthread -lm -ldl $(shell if [ -L /usr/lib64/libmemkind.so ]; then echo '-lmemkind'; fi)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
