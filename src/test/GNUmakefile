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

all: hif_test.exe jacstr_test.exe err_test.exe rot_test.exe zhegvt.exe zggsvd3t.exe zgent.exe

help:
	@echo "make [CPU=x64|x100|x200] [NDEBUG=0|1|2|3|4|5] [NEST=_SEQ] [all|clean|help]"

# TODO: specify the module dependencies for .o targets.

hif_test.exe: hif_test.o ../shared/libFLAPWxHZ.a $(MKFS)
	$(FC) $(FFLAGS) hif_test.o -o$@ -L../shared -lFLAPWxHZ $(LDFLAGS)

hif_test.o: hif_test.F90 $(MKFS)
	$(FC) $(FFLAGS) -I../shared -c hif_test.F90

jacstr_test.exe: jacstr_test.o ../shared/libFLAPWxHZ.a ../../../JACSD/libjstrat.a $(MKFS)
	$(FC) $(FFLAGS) jacstr_test.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -ljstrat $(LDFLAGS)

jacstr_test.o: jacstr_test.F90 $(MKFS)
	$(FC) $(FFLAGS) -I../shared -c jacstr_test.F90

err_test.exe: err_test.o ../shared/libFLAPWxHZ.a $(MKFS)
	$(FC) $(FFLAGS) err_test.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

err_test.o: err_test.F90 $(MKFS)
	$(FC) $(FFLAGS) -I../shared -c err_test.F90

rot_test.exe: rot_test.o $(MKFS)
	$(FC) $(FFLAGS) rot_test.o -o$@

rot_test.o: rot_test.F90 $(MKFS)
	$(FC) $(FFLAGS) -c rot_test.F90

zhegvt.exe: zhegvt.o ../shared/libFLAPWxHZ.a $(MKFS)
	$(FC) $(FFLAGS) zhegvt.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zhegvt.o: zhegvt.F90 $(MKFS)
	$(FC) $(FFLAGS) -I../shared -c zhegvt.F90

zggsvd3t.exe: zggsvd3t.o ../shared/libFLAPWxHZ.a $(MKFS)
	$(FC) $(FFLAGS) zggsvd3t.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zggsvd3t.o: zggsvd3t.F90 $(MKFS)
	$(FC) $(FFLAGS) -I../shared -c zggsvd3t.F90

zgent.exe: zgent.o ../shared/libFLAPWxHZ.a $(MKFS)
	$(FC) $(FFLAGS) zgent.o -o$@ -L../shared -lFLAPWxHZ -L../../../JACSD -lvn $(LDFLAGS)

zgent.o: zgent.F90 $(MKFS)
	$(FC) $(FFLAGS) -I../shared -c zgent.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
