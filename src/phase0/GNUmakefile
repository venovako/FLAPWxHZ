ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: phase0.exe

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

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
