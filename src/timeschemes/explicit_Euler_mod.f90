module explicit_Euler_mod

use stvec_mod,        only : stvec_t
use timescheme_mod,   only : timescheme_t
use operator_mod,     only : operator_t
use multi_domain_mod, only : multi_domain_t

implicit none

type, public, extends(timescheme_t) :: explicit_Euler_t
  class(stvec_t), allocatable :: tendency
contains
  procedure, public :: step => step_explicit_Euler
end type explicit_Euler_t

contains

subroutine step_explicit_Euler(this, v0, operator, multi_domain, dt)

  class(explicit_Euler_t), intent(inout) :: this
  class(stvec_t),          intent(inout) :: v0
  class(operator_t),       intent(inout) :: operator
  type(multi_domain_t),    intent(in)    :: multi_domain
  real(kind=8),            intent(in)    :: dt

  call operator%apply(this%tendency, v0, multi_domain)
  call v0%update(dt, this%tendency, multi_domain)

end subroutine step_explicit_Euler

end module explicit_Euler_mod
