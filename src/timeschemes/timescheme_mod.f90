module timescheme_mod

use operator_mod,     only : operator_t
use stvec_mod,        only : stvec_t
use multi_domain_mod, only : multi_domain_t
implicit none

type, abstract, public :: timescheme_t
contains
  procedure(step),  deferred :: step
end type timescheme_t

abstract interface
  subroutine step(this, v0, operator, multi_domain, dt)
    import operator_t, timescheme_t, stvec_t, multi_domain_t

    class(timescheme_t), intent(inout) :: this
    class(stvec_t),      intent(inout) :: v0
    class(operator_t),   intent(inout) :: operator
    type(multi_domain_t),intent(in)    :: multi_domain
    real(kind=8),        intent(in)    :: dt
  end subroutine step
end interface

contains

end module timescheme_mod
