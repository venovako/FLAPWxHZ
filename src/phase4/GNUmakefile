ifndef COMPILER
COMPILER=x64x
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: phase4.exe phase4d.exe

help:
	@echo "gmake [COMPILER=x64x|x200|x64|gnu|nvidia] [NDEBUG=optimization_level] [all|clean|help]"

phase4.exe: phase4.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) phase4.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

phase4.o: phase4.F90 readcl.F90 bio.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c phase4.F90

phase4d.exe: phase4d.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) phase4d.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

phase4d.o: phase4d.F90 readcl.F90 biod.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c phase4d.F90

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
