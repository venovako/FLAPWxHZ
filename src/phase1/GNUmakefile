ifeq ($(CPU),x64) # Xeon / Intel Fortran
include ../x64.mk
MKFS=GNUmakefile ../x64.mk
else ifeq ($(CPU),x200) # Knights Landing / Intel Fortran
include ../x200.mk
MKFS=GNUmakefile ../x200.mk
else ifeq ($(CPU),llvm) # Flang
include ../llvm.mk
MKFS=GNUmakefile ../llvm.mk
else # GNU Fortran
ifndef CPU
CPU=gnu
endif # !CPU
include ../gnu.mk
MKFS=GNUmakefile ../gnu.mk
endif # ?CPU

.PHONY: all help clean

all: phase1.exe

help:
	@echo "gmake [CPU=x64|x100|x200|llvm|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

phase1.exe: phase1.o ../shared/libFLAPWxHZ.a ../../../JACSD/libjstrat.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) phase1.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -ljstrat -lvn $(LDFLAGS)

phase1.o: phase1.F90 par_work.F90 bio.F90 readcl.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/dtypes.mod ../shared/hif_hz.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/jstrat -I../../../JACSD/vn -c phase1.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
