#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =
FC      = ifort
OPTSC   =  -c -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -module mod
OPTSL   =  -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)DIV_TEST: $(MKDIRS) $(DOBJ)div_test.o
	@rm -f $(filter-out $(DOBJ)div_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DIV_TEST
$(DEXE)FIELD_TEST: $(MKDIRS) $(DOBJ)field_test.o
	@rm -f $(filter-out $(DOBJ)field_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) FIELD_TEST
$(DEXE)GRAD_TEST: $(MKDIRS) $(DOBJ)grad_test.o
	@rm -f $(filter-out $(DOBJ)grad_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) GRAD_TEST
$(DEXE)CURL_TEST: $(MKDIRS) $(DOBJ)curl_test.o
	@rm -f $(filter-out $(DOBJ)curl_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) CURL_TEST
$(DEXE)DOMAIN_TEST: $(MKDIRS) $(DOBJ)domain_test.o
	@rm -f $(filter-out $(DOBJ)domain_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DOMAIN_TEST
$(DEXE)SBP_OPERATORS_TEST: $(MKDIRS) $(DOBJ)sbp_operators_test.o
	@rm -f $(filter-out $(DOBJ)sbp_operators_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) SBP_OPERATORS_TEST
$(DEXE)CENTRAL_OPERATORS_TEST: $(MKDIRS) $(DOBJ)central_operators_test.o
	@rm -f $(filter-out $(DOBJ)central_operators_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) CENTRAL_OPERATORS_TEST

#compiling rules
$(DOBJ)domain_mod.o: src/domain_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_mod.o: src/grad_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)differential_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_differential_operator_mod.o: src/sbp_differential_operator_mod.f90 \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_mod.o: src/div_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)differential_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_mod.o: src/curl_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)differential_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)central_differential_operator_mod.o: src/central_differential_operator_mod.f90 \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)field_mod.o: src/field_mod.f90 \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)differential_operator_mod.o: src/differential_operator_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)const_mod.o: src/const_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_test.o: src/tests/div_test.f90 \
	$(DOBJ)div_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)field_test.o: src/tests/field_test.f90 \
	$(DOBJ)field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_test.o: src/tests/grad_test.f90 \
	$(DOBJ)grad_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_test.o: src/tests/curl_test.f90 \
	$(DOBJ)curl_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_test.o: src/tests/domain_test.f90 \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_operators_test.o: src/tests/sbp_operators_test.f90 \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)central_operators_test.o: src/tests/central_operators_test.f90 \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
all: $(addprefix $(DEXE),$(EXES))
