module state_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
implicit none

  type, abstract :: state_t

    type(field_t) field

  contains

  end type state_t

contains

end module state_mod
