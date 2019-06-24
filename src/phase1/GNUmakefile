ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: phase1.exe

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

phase1.exe: phase1.o ../shared/libFLAPWxHZ.a ../../../JACSD/libjstrat.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) phase1.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -ljstrat -lvn $(LDFLAGS)

phase1.o: phase1.F90 par_work.F90 bio.F90 readcl.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/dtypes.mod ../shared/hif_hz.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod $(MKFS)
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
