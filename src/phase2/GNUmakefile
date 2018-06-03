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

all:

help:
	@echo "make [CPU=x64|x100|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

clean:
