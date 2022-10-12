module operator_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
implicit none

  type, abstract :: operator_t

    contains

      procedure(apply_i), deferred :: apply

  end type operator_t

abstract interface
  subroutine apply_i(this, out, in, domain)
    import operator_t, domain_t, field_t

    class    (operator_t), intent(in) :: this
    type     (field_t),    intent(inout) :: out
    type     (field_t),    intent(in)  :: in
    type     (domain_t),   intent(in)  :: domain

  end subroutine apply_i
end interface


contains

end module operator_mod
