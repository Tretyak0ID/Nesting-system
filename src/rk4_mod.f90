module rk4_mod

use stvec_mod,      only : stvec_t
use timescheme_mod, only : timescheme_t
use operator_mod,   only : operator_t
use domain_mod,       only : domain_t

implicit none

private

type, extends(timescheme_t), public :: rk4_t
    class(stvec_t), allocatable :: k1, k2, k3, k4, y
contains
    procedure, public :: step => step_rk4
end type rk4_t

contains

subroutine step_rk4(this, v0, operator, domain, dt)

    class(rk4_t),      intent(inout) :: this
    class(stvec_t),    intent(inout) :: v0
    class(operator_t), intent(inout) :: operator
    type(domain_t),      intent(in)    :: domain
    real(kind=8),      intent(in)    :: dt

    call operator%apply(this%k1, v0, domain)

    call this%y%assign(1.0_8, v0, 0.5_8*dt, this%k1, domain)
    call operator%apply(this%k2, this%y, domain)

    call this%y%assign(1.0_8, v0, 0.5_8*dt, this%k2, domain)
    call operator%apply(this%k3, this%y, domain)

    call this%y%assign(1.0_8, v0, dt, this%k3, domain)
    call operator%apply(this%k4, this%y, domain)

    call v0%update(dt/6.0_8, this%k1, dt/3.0_8, this%k2, domain)
    call v0%update(dt/3.0_8, this%k3, dt/6.0_8, this%k4, domain)

end subroutine step_rk4

end module rk4_mod
