ifeq ($(CPU),x64) # Xeon / Intel Fortran
include ../x64.mk
MKFS=GNUmakefile ../x64.mk
else ifeq ($(CPU),x100) # Knights Corner / Intel Fortran
include ../x100.mk
MKFS=GNUmakefile ../x100.mk
else ifeq ($(CPU),x200) # Knights Landing / Intel Fortran
include ../x200.mk
MKFS=GNUmakefile ../x200.mk
else # !CPU
CPU=x86_64 # GNU Fortran
include ../gnu.mk
MKFS=GNUmakefile ../gnu.mk
endif # ?CPU

.PHONY: all help clean

all: hif_test.exe jstrat_test.exe err_test.exe derr_test.exe rot_test.exe zhegvt.exe zggsvd3t.exe zgent.exe zhzl1sat.exe zksqrt2t.exe zjht.exe binvec.exe

help:
	@echo "make [CPU=x64|x100|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

hif_test.exe: hif_test.o ../shared/libFLAPWxHZ.a ../../../JACSD/libjstrat.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) hif_test.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -ljstrat -lvn $(LDFLAGS)

hif_test.o: hif_test.F90 ../shared/blas_utils.mod ../shared/dtypes.mod ../shared/hif_hz.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/jstrat -I../../../JACSD/vn -c hif_test.F90

jstrat_test.exe: jstrat_test.o ../../../JACSD/libjstrat.a $(MKFS)
	$(FC) $(FFLAGS) jstrat_test.o -o$@ -L../../../JACSD -ljstrat $(LDFLAGS)

jstrat_test.o: jstrat_test.F90 ../../../JACSD/jstrat/jstrat_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/jstrat -c jstrat_test.F90

err_test.exe: err_test.o ../shared/libFLAPWxHZ.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) err_test.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

err_test.o: err_test.F90 ../shared/binio.mod ../shared/params.mod ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/vn -c err_test.F90

derr_test.exe: derr_test.o ../shared/libFLAPWxHZ.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) derr_test.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

derr_test.o: derr_test.F90 ../shared/binio.mod ../shared/params.mod ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/vn -c derr_test.F90

rot_test.exe: rot_test.o $(MKFS)
	$(FC) $(FFLAGS) rot_test.o -o$@

rot_test.o: rot_test.F90 $(MKFS)
	$(FC) $(FFLAGS) -c rot_test.F90

zhegvt.exe: zhegvt.o ../shared/libFLAPWxHZ.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) zhegvt.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zhegvt.o: zhegvt.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/vn -c zhegvt.F90

zggsvd3t.exe: zggsvd3t.o ../shared/libFLAPWxHZ.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) zggsvd3t.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zggsvd3t.o: zggsvd3t.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/vn -c zggsvd3t.F90

zgent.exe: zgent.o ../shared/libFLAPWxHZ.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) zgent.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zgent.o: zgent.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/params.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/vn -c zgent.F90

zhzl1sat.exe: zhzl1sat.o ../shared/libFLAPWxHZ.a ../../../JACSD/libjstrat.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) zhzl1sat.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -ljstrat -lvn $(LDFLAGS)

zhzl1sat.o: zhzl1sat.F90 ../shared/binio.mod ../shared/blas_utils.mod ../shared/dtypes.mod ../shared/hif_hz.mod ../shared/params.mod ../shared/timer.mod ../../../JACSD/jstrat/jstrat_f.mod ../../../JACSD/vn/vn_binio_f.mod ../../../JACSD/vn/vn_blas_f.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/jstrat -I../../../JACSD/vn -c zhzl1sat.F90

zksqrt2t.exe: zksqrt2t.o ../shared/libFLAPWxHZ.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) zksqrt2t.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zksqrt2t.o: zksqrt2t.F90 ../shared/jqr.mod ../shared/blas_utils.mod ../shared/params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/vn -c zksqrt2t.F90

zjht.exe: zjht.o ../shared/libFLAPWxHZ.a ../../../JACSD/libvn.a  $(MKFS)
	$(FC) $(FFLAGS) zjht.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zjht.o: zjht.F90 ../shared/jqr.mod ../shared/blas_utils.mod ../shared/params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../shared -I../../../JACSD/vn -c zjht.F90

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
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
