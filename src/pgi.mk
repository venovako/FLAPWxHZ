SHELL=/bin/bash
ARCH=$(shell uname)
ifneq ($(ARCH),Linux)
$(error PGI is supported on Linux only)
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
CPUFLAGS=-DUSE_PGI -DUSE_X64 -DKIND_QUAD=$(WP) -DOLD_OMP -m64 -mp -KPIC -Mframe -Meh_frame -Minfo -Mnollvm
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -Mnoinline -Minstrument=functions
endif # PROFILE
CPUFLAGS += -DTSC_FREQ_HZ=$(shell if [ `if [ -r /etc/redhat-release ]; then grep -c 'release 7' /etc/redhat-release; else echo 0; fi` = 1 ]; then echo `dmesg | grep 'TSC clocksource calibration' | cut -d':' -f3 | cut -d' ' -f2 | sed 's/\.//g'`000; else echo 0; fi)ull
FORFLAGS=$(CPUFLAGS) -i8 -Mdclchk -Mlarge_arrays -Mrecursive -Mstack_arrays
C11FLAGS=$(CPUFLAGS) -c11
CC=pgcc
FC=pgfortran
FPUFLAGS=-Kieee -Mfma -Mnodaz -Mnoflushz -Mnofpapprox -Mnofprelaxed
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
DBGFLAGS=-g -Mbounds -Mchkstk -traceback
DBGFFLAGS=$(DBGFLAGS)
DBGCFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
# define MKL=pgi_thread for a threaded MKL
ifndef MKL
MKL=sequential
endif # !MKL
LIBFLAGS=-D_GNU_SOURCE -DUSE_MKL -DMKL_ILP64 -I. -I../../../JACSD/jstrat -I../../../JACSD/vn -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
LDFLAGS=-L../../../JACSD -lvn$(PROFILE)$(DEBUG) -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_$(MKL) -lmkl_core -pgf90libs -lpthread -lm -ldl $(shell if [ -L /usr/lib64/libmemkind.so ]; then echo '-lmemkind'; fi)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUCFLAGS)
