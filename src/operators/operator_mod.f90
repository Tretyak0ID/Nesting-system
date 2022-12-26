module operator_mod
use stvec_mod,        only: stvec_t
use multi_domain_mod, only: multi_domain_t
implicit none

  type, abstract :: operator_t
  !Abstract type for RHS operators
    contains

      procedure(apply_i), deferred :: apply

  end type operator_t

abstract interface
  subroutine apply_i(this, out, in, multi_domain)
    import operator_t, multi_domain_t, stvec_t

    class    (operator_t),     intent(inout) :: this
    class    (stvec_t),        intent(inout) :: out
    class    (stvec_t),        intent(inout) :: in
    type     (multi_domain_t), intent(in)    :: multi_domain

  end subroutine apply_i
end interface


contains

end module operator_mod
