ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: phase4.exe phase4d.exe

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

phase4.exe: phase4.o ../shared/libFLAPWxHZ$(PROFILE)$(DEBUG).a ../../../JACSD/libvn$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) phase4.o -o$@ -L../shared -lFLAPWxHZ$(PROFILE)$(DEBUG) $(LDFLAGS)

phase4.o: phase4.F90 readcl.F90 bio.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c phase4.F90

phase4d.exe: phase4d.o ../shared/libFLAPWxHZ$(PROFILE)$(DEBUG).a ../../../JACSD/libvn$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) phase4d.o -o$@ -L../shared -lFLAPWxHZ$(PROFILE)$(DEBUG) $(LDFLAGS)

phase4d.o: phase4d.F90 readcl.F90 biod.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c phase4d.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
