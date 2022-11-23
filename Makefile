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
$(DEXE)DIV_SAT_TEST: $(MKDIRS) $(DOBJ)div_sat_test.o
	@rm -f $(filter-out $(DOBJ)div_sat_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DIV_SAT_TEST
$(DEXE)TEST_0_HORIZONTAL_ADVECTION: $(MKDIRS) $(DOBJ)test_0_horizontal_advection.o
	@rm -f $(filter-out $(DOBJ)test_0_horizontal_advection.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_0_HORIZONTAL_ADVECTION
$(DEXE)GRAD_SAT_TEST: $(MKDIRS) $(DOBJ)grad_sat_test.o
	@rm -f $(filter-out $(DOBJ)grad_sat_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) GRAD_SAT_TEST
$(DEXE)DIV_TEST: $(MKDIRS) $(DOBJ)div_test.o
	@rm -f $(filter-out $(DOBJ)div_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) DIV_TEST
$(DEXE)MULTI_GRID_FIELD_TEST: $(MKDIRS) $(DOBJ)multi_grid_field_test.o
	@rm -f $(filter-out $(DOBJ)multi_grid_field_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MULTI_GRID_FIELD_TEST
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
$(DEXE)INTERPOLATION_TEST: $(MKDIRS) $(DOBJ)interpolation_test.o
	@rm -f $(filter-out $(DOBJ)interpolation_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) INTERPOLATION_TEST
$(DEXE)TEST_3_GEOSTROPHIC_CYCLONE: $(MKDIRS) $(DOBJ)test_3_geostrophic_cyclone.o
	@rm -f $(filter-out $(DOBJ)test_3_geostrophic_cyclone.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_3_GEOSTROPHIC_CYCLONE
$(DEXE)ADVECTIVE_CALCULATE_TEST: $(MKDIRS) $(DOBJ)advective_calculate_test.o
	@rm -f $(filter-out $(DOBJ)advective_calculate_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) ADVECTIVE_CALCULATE_TEST
$(DEXE)MULTI_DOMAIN_TEST: $(MKDIRS) $(DOBJ)multi_domain_test.o
	@rm -f $(filter-out $(DOBJ)multi_domain_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) MULTI_DOMAIN_TEST
$(DEXE)SWE_ADVECTION_OPERATOR_TEST: $(MKDIRS) $(DOBJ)swe_advection_operator_test.o
	@rm -f $(filter-out $(DOBJ)swe_advection_operator_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) SWE_ADVECTION_OPERATOR_TEST
$(DEXE)TEST_1_GAUSSIAN_HILL: $(MKDIRS) $(DOBJ)test_1_gaussian_hill.o
	@rm -f $(filter-out $(DOBJ)test_1_gaussian_hill.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_1_GAUSSIAN_HILL
$(DEXE)VEC_MATH_TEST: $(MKDIRS) $(DOBJ)vec_math_test.o
	@rm -f $(filter-out $(DOBJ)vec_math_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) VEC_MATH_TEST
$(DEXE)TEST_2_GEOSTROPHIC_BALANCE: $(MKDIRS) $(DOBJ)test_2_geostrophic_balance.o
	@rm -f $(filter-out $(DOBJ)test_2_geostrophic_balance.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_2_GEOSTROPHIC_BALANCE
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
$(DEXE)TIMESCHEME_TEST: $(MKDIRS) $(DOBJ)timescheme_test.o
	@rm -f $(filter-out $(DOBJ)timescheme_test.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TIMESCHEME_TEST
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

$(DOBJ)read_write_mod.o: src/read_write_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)initial_condition_mod.o: src/initial_condition_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_mod.o: src/stvec_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)multi_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vec_math_mod.o: src/vec_math_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)multi_grid_field_mod.o: src/multi_grid_field_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)multi_domain_mod.o: src/multi_domain_mod.f90 \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)field_mod.o: src/field_mod.f90 \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)const_mod.o: src/const_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_swe_mod.o: src/stvec_swe_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)stvec_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)operator_mod.o: src/operators/operator_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)multi_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swe_advective_operator_mod.o: src/operators/swe_advective_operator_mod.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)grad_mod.o \
	$(DOBJ)div_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swe_vect_inv_operator_mod.o: src/operators/swe_vect_inv_operator_mod.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)grad_mod.o \
	$(DOBJ)div_mod.o \
	$(DOBJ)curl_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)horizontal_advection_operator_mod.o: src/operators/horizontal_advection_operator_mod.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)grad_mod.o \
	$(DOBJ)div_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_sat_test.o: src/tests/div_SAT_test.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)div_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_0_horizontal_advection.o: src/tests/test_0_horizontal_advection.f90 \
	$(DOBJ)initial_condition_mod.o \
	$(DOBJ)swe_advective_operator_mod.o \
	$(DOBJ)swe_vect_inv_operator_mod.o \
	$(DOBJ)horizontal_advection_operator_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timesheme_factory_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)explicit_euler_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)read_write_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_sat_test.o: src/tests/grad_SAT_test.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)grad_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o
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

$(DOBJ)multi_grid_field_test.o: src/tests/multi_grid_field_test.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)stvec_swe_test.o: src/tests/stvec_swe_test.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)field_test.o: src/tests/field_test.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolation_test.o: src/tests/interpolation_test.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)read_write_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)interpolation_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_3_geostrophic_cyclone.o: src/tests/test_3_geostrophic_cyclone.f90 \
	$(DOBJ)initial_condition_mod.o \
	$(DOBJ)swe_vect_inv_operator_mod.o \
	$(DOBJ)swe_advective_operator_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)curl_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timesheme_factory_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)read_write_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)advective_calculate_test.o: src/tests/advective_calculate_test.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)swe_advective_operator_mod.o \
	$(DOBJ)rk4_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)multi_domain_test.o: src/tests/multi_domain_test.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)swe_advection_operator_test.o: src/tests/swe_advection_operator_test.f90 \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)swe_advective_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_1_gaussian_hill.o: src/tests/test_1_gaussian_hill.f90 \
	$(DOBJ)initial_condition_mod.o \
	$(DOBJ)swe_advective_operator_mod.o \
	$(DOBJ)swe_vect_inv_operator_mod.o \
	$(DOBJ)horizontal_advection_operator_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timesheme_factory_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)explicit_euler_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)read_write_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)vec_math_test.o: src/tests/vec_math_test.f90 \
	$(DOBJ)div_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)test_2_geostrophic_balance.o: src/tests/test_2_geostrophic_balance.f90 \
	$(DOBJ)initial_condition_mod.o \
	$(DOBJ)swe_vect_inv_operator_mod.o \
	$(DOBJ)swe_advective_operator_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)central_differential_operator_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timesheme_factory_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)read_write_mod.o
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

$(DOBJ)timescheme_test.o: src/tests/timescheme_test.f90 \
	$(DOBJ)initial_condition_mod.o \
	$(DOBJ)swe_advective_operator_mod.o \
	$(DOBJ)swe_vect_inv_operator_mod.o \
	$(DOBJ)horizontal_advection_operator_mod.o \
	$(DOBJ)sbp_differential_operator_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timesheme_factory_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)explicit_euler_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)stvec_swe_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)read_write_mod.o
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

$(DOBJ)timescheme_mod.o: src/timeschemes/timescheme_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)multi_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)timesheme_factory_mod.o: src/timeschemes/timesheme_factory_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)explicit_euler_mod.o \
	$(DOBJ)rk4_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)rk4_mod.o: src/timeschemes/rk4_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)stvec_swe_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)explicit_euler_mod.o: src/timeschemes/explicit_Euler_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)multi_domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sat_mod.o: src/differential_operators/SAT_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)interpolation_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)grad_mod.o: src/differential_operators/grad_mod.f90 \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)sat_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)sbp_differential_operator_mod.o: src/differential_operators/sbp_differential_operator_mod.f90 \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)div_mod.o: src/differential_operators/div_mod.f90 \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)sat_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)curl_mod.o: src/differential_operators/curl_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)multi_domain_mod.o \
	$(DOBJ)multi_grid_field_mod.o \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)sat_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)central_differential_operator_mod.o: src/differential_operators/central_differential_operator_mod.f90 \
	$(DOBJ)differential_operator_mod.o \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)differential_operator_mod.o: src/differential_operators/differential_operator_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC)  $< -o $@

$(DOBJ)interpolation_mod.o: src/differential_operators/interpolation_mod.f90 \
	$(DOBJ)field_mod.o \
	$(DOBJ)multi_grid_field_mod.o
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
