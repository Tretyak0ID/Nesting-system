module swe_state_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
use state_mod, only: state_t
implicit none

  type, extends(state_t) :: swe_state_t

    type(field_t) u, v, h

  contains

  end type swe_state_t

contains

end module swe_state_mod
