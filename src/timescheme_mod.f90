module timescheme_mod

use operator_mod,  only : operator_t
use stvec_mod,     only : stvec_t
use domain_mod,      only : domain_t
implicit none

type, abstract, public :: timescheme_t
contains
  procedure(step),  deferred :: step
end type timescheme_t

abstract interface
  subroutine step(this, v0, operator, domain, dt)
    import operator_t, timescheme_t, stvec_t, domain_t

    class(timescheme_t), intent(inout) :: this
    class(stvec_t),      intent(inout) :: v0
    class(operator_t),   intent(inout) :: operator
    type(domain_t),        intent(in)    :: domain
    real(kind=8),        intent(in)    :: dt
  end subroutine step
end interface

contains

end module timescheme_mod
