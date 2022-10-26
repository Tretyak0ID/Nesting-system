module timescheme_factory_mod

use stvec_mod,         only : stvec_t
use timescheme_mod,    only : timescheme_t
use explicit_Euler_mod,only : explicit_Euler_t
use rk4_mod,           only : rk4_t

implicit none

contains

subroutine create_timescheme(timescheme, v, timescheme_name)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),      intent(in) :: v
    character(len=*),    intent(in) :: timescheme_name

    select case(timescheme_name)
    case("explicit_Euler")
        call create_explicit_Euler_timescheme(timescheme, v)
    case("rk4")
        call create_rk4_timescheme(timescheme, v)
    case default
    end select

end subroutine create_timescheme

subroutine create_explicit_Euler_timescheme(timescheme, v)
    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v

    type(explicit_Euler_t), allocatable :: Euler_timescheme

    allocate(Euler_timescheme)
    call v%create_similar(Euler_timescheme%tendency)
    call move_alloc(Euler_timescheme, timescheme)

end subroutine create_explicit_Euler_timescheme

subroutine create_rk4_timescheme(timescheme, v)

    class(timescheme_t), allocatable, intent(out) :: timescheme
    class(stvec_t),                   intent(in)  :: v

    type(rk4_t), allocatable :: rk4

    allocate(rk4)

    call v%create_similar(rk4%k1)
    call v%create_similar(rk4%k2)
    call v%create_similar(rk4%k3)
    call v%create_similar(rk4%k4)
    call v%create_similar(rk4%y)

    call move_alloc(rk4, timescheme)
end subroutine create_rk4_timescheme

end module timescheme_factory_mod
