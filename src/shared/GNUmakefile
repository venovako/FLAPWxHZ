ifndef COMPILER
COMPILER=x64x
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: libFLAPWxHZ$(DEBUG).a

help:
	@echo "gmake [COMPILER=x64x|x200|x64|gnu|nvidia] [NDEBUG=optimization_level] [all|clean|help]"

libFLAPWxHZ$(DEBUG).a: binio.o blas_utils.o dtypes.o hif_hz.o jqr.o params.o timer.o $(MKFS)
	$(AR) $(ARFLAGS) $@ binio.o blas_utils.o dtypes.o hif_hz.o jqr.o params.o timer.o

binio.o binio.mod: binio.F90 ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio.F90

blas_utils.o blas_utils.mod: blas_utils.F90 blas_zrotm_xx.F90 params.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c blas_utils.F90

dtypes.o dtypes.mod: dtypes.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtypes.F90

hif_hz.o hif_hz.mod: hif_hz.F90 hif.F90 hz.F90 janal.F90 zsetzi.F90 zjscl.F90 zimul.F90 zjmul.F90 zcopyh.F90 zprep_blks.F90 zupdate.F90 zxfer2.F90 zupxfer.F90 zhzl1.F90 zhzl1sa.F90 zhzl2.F90 blas_utils.mod dtypes.mod timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) -I../../../JACSD/jstrat $(FFLAGS) -c hif_hz.F90

jqr.o jqr.mod: jqr.F90 blas_utils.mod params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c jqr.F90

params.o params.mod: params.F90 ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 params.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c timer.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
