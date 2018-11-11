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

all: libFLAPWxHZ.a

help:
	@echo "make [CPU=x64|x100|x200] [NDEBUG=0|1|2|3|4|5] [NEST=_SEQ] [all|clean|help]"

libFLAPWxHZ.a: binio.o binio1.o blas_utils.o dtypes.o hif_hz.o jacstr.o jqr.o params.o timer.o timer_.o vn_binio.o $(MKFS) #ximaginary.o
	$(AR) $(ARFLAGS) $@ binio.o binio1.o blas_utils.o dtypes.o hif_hz.o jacstr.o jqr.o params.o timer.o timer_.o vn_binio.o #ximaginary.o

binio.o binio.mod: binio.F90 params.mod vn_binio.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio.F90

binio1.o binio1.mod: binio1.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio1.F90

blas_utils.o blas_utils.mod: blas_utils.F90 blas_zrotm_xx.F90 params.mod $(MKFS) #ximaginary.mod
	$(FC) $(FFLAGS) -c blas_utils.F90

dtypes.o dtypes.mod: dtypes.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtypes.F90

hif_hz.o hif_hz.mod: hif_hz.F90 hif.F90 hz.F90 janal.F90 zsetzi.F90 zjscl.F90 zimul.F90 zjmul.F90 zcopyh.F90 zprep_blks.F90 zupdate.F90 zxfer2.F90 zupxfer.F90 zhzl1.F90 zhzl1s.F90 zhzl1p.F90 zhzl2.F90 blas_utils.mod dtypes.mod jacstr.mod timer.mod $(MKFS)
	$(FC) $(FFLAGS) -c hif_hz.F90

jacstr.o jacstr.mod: jacstr.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c jacstr.F90

jqr.o jqr.mod: jqr.F90 blas_utils.mod $(MKFS)
	$(FC) $(FFLAGS) -c jqr.F90

params.o params.mod: params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 timer_.h $(MKFS)
	$(FC) $(FFLAGS) -c timer.F90

timer_.o: timer_.c timer_.h $(MKFS)
	$(CC) -D_GNU_SOURCE $(CFLAGS) -c timer_.c

vn_binio.o vn_binio.mod: vn_binio.F90 $(MKFS)
	$(FC) $(FFLAGS) -c vn_binio.F90

ximaginary.o ximaginary.mod: ximaginary.F90 $(MKFS)
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
