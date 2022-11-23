module differential_operator_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
implicit none

  type, abstract :: differential_operator_t

      character(len=5), public :: name

    contains

      procedure(apply_i), deferred :: apply

  end type differential_operator_t

abstract interface
  subroutine apply_i(this, out, in, domain, direction)

    import differential_operator_t, domain_t, field_t

    class    (differential_operator_t), intent(in)  :: this
    type     (field_t),                 intent(inout) :: out
    type     (field_t),                 intent(in)  :: in
    type     (domain_t),                intent(in)  :: domain
    character(len=1),                   intent(in)  :: direction

  end subroutine apply_i
end interface


contains

end module differential_operator_mod
