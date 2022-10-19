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
$(DEXE)STVEC_SWE_TEST: $(MKDIRS) $(DOBJ)stvec_swe_test.o
	@rm -f $(filter-out $(DOBJ)stvec_swe_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) STVEC_SWE_TEST
$(DEXE)FIELD_TEST: $(MKDIRS) $(DOBJ)field_test.o
	@rm -f $(filter-out $(DOBJ)field_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) FIELD_TEST
$(DEXE)MESH_TEST: $(MKDIRS) $(DOBJ)mesh_test.o
	@rm -f $(filter-out $(DOBJ)mesh_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MESH_TEST
$(DEXE)ADVECTIVE_CALCULATE_TEST: $(MKDIRS) $(DOBJ)advective_calculate_test.o
	@rm -f $(filter-out $(DOBJ)advective_calculate_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) ADVECTIVE_CALCULATE_TEST
$(DEXE)SWE_ADVECTION_OPERATOR_TEST: $(MKDIRS) $(DOBJ)swe_advection_operator_test.o
	@rm -f $(filter-out $(DOBJ)swe_advection_operator_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) SWE_ADVECTION_OPERATOR_TEST
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
$(DOBJ)timescheme_mod.o: src/timescheme_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)domain_mod.o: src/domain_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_mod.o: src/operator_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_mod.o: src/grad_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)differential_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swe_advective_operator_mod.o: src/swe_advective_operator_mod.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grad_mod.o \
	$(DOBJ)div_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_differential_operator_mod.o: src/sbp_differential_operator_mod.f90 \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_mod.o: src/stvec_mod.f90 \
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

$(DOBJ)vec_math_mod.o: src/vec_math_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
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

$(DOBJ)rk4_mod.o: src/rk4_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_swe_mod.o: src/stvec_swe_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)explicit_euler_mod.o: src/explicit_Euler_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
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

$(DOBJ)stvec_swe_test.o: src/tests/stvec_swe_test.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)field_test.o: src/tests/field_test.f90 \
	$(DOBJ)field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)mesh_test.o: src/tests/mesh_test.f90 \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)advective_calculate_test.o: src/tests/advective_calculate_test.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)swe_advective_operator_mod.o \
	$(DOBJ)rk4_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swe_advection_operator_test.o: src/tests/swe_advection_operator_test.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)swe_advective_operator_mod.o
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
