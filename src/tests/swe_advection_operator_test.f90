program swe_advection_operator_test
use stvec_swe_mod, only: stvec_swe_t
use swe_advective_operator_mod, only: swe_advective_operator_t
implicit none

  type(swe_advective_operator_t) :: swe
  type(stvec_swe_t)              :: state

end program swe_advection_operator_test
