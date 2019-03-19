ifeq ($(CPU),x64) # Xeon / Intel Fortran
include ../x64.mk
MKFS=GNUmakefile ../x64.mk
else ifeq ($(CPU),x200) # Knights Landing / Intel Fortran
include ../x200.mk
MKFS=GNUmakefile ../x200.mk
else # !CPU
CPU=x86_64 # GNU Fortran
include ../gnu.mk
MKFS=GNUmakefile ../gnu.mk
endif # ?CPU

.PHONY: all help clean

all: libFLAPWxHZ.a

help:
	@echo "make [CPU=x64|x100|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libFLAPWxHZ.a: binio.o blas_utils.o dtypes.o hif_hz.o jqr.o params.o timer.o $(MKFS) #binio1.o ximaginary.o
	$(AR) $(ARFLAGS) $@ binio.o blas_utils.o dtypes.o hif_hz.o jqr.o params.o timer.o #binio1.o ximaginary.o

binio.o binio.mod: binio.F90 ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c binio.F90

binio1.o binio1.mod: binio1.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio1.F90

blas_utils.o blas_utils.mod: blas_utils.F90 blas_zrotm_xx.F90 params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS) #ximaginary.mod
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c blas_utils.F90

dtypes.o dtypes.mod: dtypes.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtypes.F90

hif_hz.o hif_hz.mod: hif_hz.F90 hif.F90 hz.F90 janal.F90 zsetzi.F90 zjscl.F90 zimul.F90 zjmul.F90 zcopyh.F90 zprep_blks.F90 zupdate.F90 zxfer2.F90 zupxfer.F90 zhzl1.F90 zhzl1sa.F90 zhzl2.F90 blas_utils.mod dtypes.mod timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/jstrat -I../../../JACSD/vn -c hif_hz.F90

jqr.o jqr.mod: jqr.F90 blas_utils.mod params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c jqr.F90

params.o params.mod: params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 params.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c timer.F90

ximaginary.o ximaginary.mod: ximaginary.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c ximaginary.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
