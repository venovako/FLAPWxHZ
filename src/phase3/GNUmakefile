ifndef COMPILER
COMPILER=gnu
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: phase3.exe

help:
	@echo "gmake [COMPILER=gnu|x64|x64x|x200|nvidia] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

phase3.exe: phase3.o ../shared/libFLAPWxHZ$(PROFILE)$(DEBUG).a ../../../JACSD/libjstrat$(PROFILE)$(DEBUG).a ../../../JACSD/libvn$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) phase3.o -o$@ -L../shared -lFLAPWxHZ$(PROFILE)$(DEBUG) -L../../../JACSD -ljstrat$(PROFILE)$(DEBUG) $(LDFLAGS)

phase3.o: phase3.F90 par_work.F90 bio.F90 readcl.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/dtypes.mod ../shared/hif_hz.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared -I../../../JACSD/jstrat $(FFLAGS) -c phase3.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.opt.yaml
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
