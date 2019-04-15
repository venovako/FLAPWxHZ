ifeq ($(CPU),x64) # Xeon / Intel C
include ../x64.mk
MKFS=GNUmakefile ../x64.mk
else ifeq ($(CPU),x200) # Knights Landing / Intel C
include ../x200.mk
MKFS=GNUmakefile ../x200.mk
else ifeq ($(CPU),llvm) # Clang
include ../llvm.mk
MKFS=GNUmakefile ../llvm.mk
else # GNU C or Clang
ifndef CPU
CPU=gnu
endif # !CPU
include ../gnu.mk
MKFS=GNUmakefile ../gnu.mk
endif # ?CPU

.PHONY: all help clean

all: phase0.exe

help:
	@echo "gmake [CPU=x64|x100|x200|llvm|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

phase0.exe: phase0.o $(MKFS)
	$(CC) $(CFLAGS) phase0.o -o$@ $(LDFLAGS)

phase0.o: phase0.c main.c split_data.c load_data.c $(MKFS)
	$(CC) $(CFLAGS) -c phase0.c

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
