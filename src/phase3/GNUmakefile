ifndef COMPILER
COMPILER=x64x
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: phase3.exe

help:
	@echo "gmake [COMPILER=x64x|x200|x64|gnu|nvidia] [NDEBUG=optimization_level] [all|clean|help]"

phase3.exe: phase3.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libjstrat$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) phase3.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) -L../../../JACSD -ljstrat$(DEBUG) $(LDFLAGS)

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
