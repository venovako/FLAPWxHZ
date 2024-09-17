ifndef COMPILER
COMPILER=x64x
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: hif_test.exe jstrat_test.exe err_test.exe derr_test.exe dggsvd3t.exe rot_test.exe zhegvt.exe zggsvd3t.exe zgent.exe zhzl1sat.exe zksqrt2t.exe zjht.exe zjh2.exe binvec.exe

help:
	@echo "gmake [COMPILER=x64x|x200|x64|gnu|nvidia] [NDEBUG=optimization_level] [all|clean|help]"

hif_test.exe: hif_test.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libjstrat$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) hif_test.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) -L../../../JACSD -ljstrat$(DEBUG) $(LDFLAGS)

hif_test.o: hif_test.F90 ../shared/blas_utils.mod ../shared/dtypes.mod ../shared/hif_hz.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared -I../../../JACSD/jstrat $(FFLAGS) -c hif_test.F90

jstrat_test.exe: jstrat_test.o ../../../JACSD/libjstrat$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) jstrat_test.o -o$@ -L../../../JACSD -ljstrat$(DEBUG) $(LDFLAGS)

jstrat_test.o: jstrat_test.F90 ../../../JACSD/jstrat/jstrat_f.mod $(MKFS)
	$(FC) -I../../../JACSD/jstrat $(FFLAGS) -c jstrat_test.F90

err_test.exe: err_test.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) err_test.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

err_test.o: err_test.F90 ../shared/binio.mod ../shared/params.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c err_test.F90

derr_test.exe: derr_test.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) derr_test.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

derr_test.o: derr_test.F90 ../shared/binio.mod ../shared/params.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c derr_test.F90

dggsvd3t.exe: dggsvd3t.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) dggsvd3t.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

dggsvd3t.o: dggsvd3t.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c dggsvd3t.F90

rot_test.exe: rot_test.o ../shared/libFLAPWxHZ$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) rot_test.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

rot_test.o: rot_test.F90 ../shared/params.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c rot_test.F90

zhegvt.exe: zhegvt.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zhegvt.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

zhegvt.o: zhegvt.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c zhegvt.F90

zggsvd3t.exe: zggsvd3t.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zggsvd3t.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

zggsvd3t.o: zggsvd3t.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c zggsvd3t.F90

zgent.exe: zgent.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zgent.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

zgent.o: zgent.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c zgent.F90

zhzl1sat.exe: zhzl1sat.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libjstrat$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zhzl1sat.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) -L../../../JACSD -ljstrat$(DEBUG) $(LDFLAGS)

zhzl1sat.o: zhzl1sat.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/dtypes.mod ../shared/hif_hz.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_timer_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared -I../../../JACSD/jstrat $(FFLAGS) -c zhzl1sat.F90

zksqrt2t.exe: zksqrt2t.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zksqrt2t.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

zksqrt2t.o: zksqrt2t.F90 ../shared/jqr.mod ../shared/blas_utils.mod ../shared/params.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c zksqrt2t.F90

zjht.exe: zjht.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a  $(MKFS)
	$(FC) $(FFLAGS) zjht.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

zjht.o: zjht.F90 ../shared/jqr.mod ../shared/blas_utils.mod ../shared/params.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c zjht.F90

zjh2.exe: zjh2.o ../shared/libFLAPWxHZ$(DEBUG).a ../../../JACSD/libvn$(DEBUG).a  $(MKFS)
	$(FC) $(FFLAGS) zjh2.o -o$@ -L../shared -lFLAPWxHZ$(DEBUG) $(LDFLAGS)

zjh2.o: zjh2.F90 ../shared/jqr.mod ../shared/blas_utils.mod ../shared/params.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_imaginary_f.mod ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) -I../shared $(FFLAGS) -c zjh2.F90

binvec.exe: binvec.o $(MKFS)
	$(CC) $(CFLAGS) binvec.o -o$@ $(LDFLAGS)

binvec.o: binvec.c $(MKFS)
	$(CC) $(CFLAGS) -c binvec.c

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
