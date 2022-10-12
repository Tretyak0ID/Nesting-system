module grad_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
use differential_operator_mod, only: differential_operator_t
implicit none

contains

  subroutine calc_grad(gx, gy, in, domain, diff_opx, diff_opy)
    type  (field_t),                 intent(inout) :: gx, gy
    type  (field_t),                 intent(in)    :: in
    type  (domain_t),                intent(in)    :: domain
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy

    call diff_opx%apply(gx, in, domain, 'x')
    call diff_opy%apply(gy, in, domain, 'y')
  end subroutine

end module grad_mod
