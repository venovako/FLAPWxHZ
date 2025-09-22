ifndef COMPILER
COMPILER=x64x
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean PHASE0 PHASE1 PHASE2 PHASE3 PHASE4 SHARED TEST

all: PHASE0 PHASE1 PHASE2 PHASE3 PHASE4 SHARED TEST

help:
	@echo "gmake [COMPILER=x64x|x200|gnu|nvidia] [MARCH=...] [NDEBUG=0|1|2|3|...|g] [all|clean|help]"

PHASE0: SHARED $(MKFS)
ifdef NDEBUG
	pushd phase0 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd phase0 && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

PHASE1: SHARED $(MKFS)
ifdef NDEBUG
	pushd phase1 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd phase1 && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

PHASE2: SHARED $(MKFS)
ifdef NDEBUG
	pushd phase2 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd phase2 && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

PHASE3: SHARED $(MKFS)
ifdef NDEBUG
	pushd phase3 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd phase3 && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

PHASE4: SHARED $(MKFS)
ifdef NDEBUG
	pushd phase4 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd phase4 && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

SHARED: $(MKFS)
ifdef NDEBUG
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd shared && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

TEST: SHARED $(MKFS)
ifdef NDEBUG
	pushd test && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd test && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

clean:
ifdef NDEBUG
	pushd test && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd phase4 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd phase3 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd phase2 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd phase1 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd phase0 && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
else # DEBUG
	pushd test && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd phase4 && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd phase3 && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd phase2 && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd phase1 && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd phase0 && $(MAKE) COMPILER=$(COMPILER) clean && popd
endif # ?NDEBUG
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
