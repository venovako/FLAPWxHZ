ifndef CPU
CPU=gnu
endif # !CPU
include $(CPU).mk
MKFS=GNUmakefile $(CPU).mk

.PHONY: all help clean PHASE0 PHASE1 PHASE2 PHASE3 PHASE4 SHARED TEST

all: PHASE0 PHASE1 PHASE2 PHASE3 PHASE4 SHARED TEST

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

PHASE0: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd phase0 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase0 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd phase0 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase0 && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

PHASE1: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd phase1 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase1 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd phase1 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase1 && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

PHASE2: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd phase2 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase2 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd phase2 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase2 && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

PHASE3: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd phase3 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase3 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd phase3 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase3 && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

PHASE4: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd phase4 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase4 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd phase4 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd phase4 && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

SHARED: $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

TEST: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd test && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd test && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd test && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd test && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

clean:
ifdef NDEBUG
ifdef PROFILE
	pushd test && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd phase4 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd phase3 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd phase2 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd phase1 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd phase0 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
else # !PROFILE
	pushd test && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd phase4 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd phase3 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd phase2 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd phase1 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd phase0 && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd test && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd phase4 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd phase3 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd phase2 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd phase1 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd phase0 && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
else # !PROFILE
	pushd test && $(MAKE) CPU=$(CPU) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) clean && popd
	pushd phase4 && $(MAKE) CPU=$(CPU) clean && popd
	pushd phase3 && $(MAKE) CPU=$(CPU) clean && popd
	pushd phase2 && $(MAKE) CPU=$(CPU) clean && popd
	pushd phase1 && $(MAKE) CPU=$(CPU) clean && popd
	pushd phase0 && $(MAKE) CPU=$(CPU) clean && popd
endif # ?PROFILE
endif # ?NDEBUG
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
