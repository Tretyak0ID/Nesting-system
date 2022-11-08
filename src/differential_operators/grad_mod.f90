module grad_mod
use multi_domain_mod,          only : multi_domain_t
use multi_grid_field_mod,      only : multi_grid_field_t
use differential_operator_mod, only : differential_operator_t
use SAT_mod,                   only : sbp_SAT_penalty_two_block
implicit none

contains

  subroutine calc_grad(gx, gy, in, domain, diff_opx, diff_opy)
    type  (multi_grid_field_t),      intent(inout) :: gx, gy
    type  (multi_grid_field_t),      intent(in)    :: in
    type  (multi_domain_t),          intent(in)    :: domain
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy

    integer(kind=8) :: n, m
    do n = 1, domain%num_sub_x
      do m = 1, domain%num_sub_y
        call diff_opx%apply(gx%subfields(n, m), in%subfields(n, m), domain%subdomains(n, m), 'x')
        call diff_opy%apply(gy%subfields(n, m), in%subfields(n, m), domain%subdomains(n, m), 'y')
      end do
    end do

    call sbp_SAT_penalty_two_block(gx, in, 'x', domain, 'sbp21')
    call sbp_SAT_penalty_two_block(gy, in, 'y', domain, 'sbp21')

  end subroutine

end module grad_mod
