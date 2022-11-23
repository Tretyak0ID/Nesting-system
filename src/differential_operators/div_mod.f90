module div_mod
use multi_domain_mod,          only : multi_domain_t
use multi_grid_field_mod,      only : multi_grid_field_t
use differential_operator_mod, only : differential_operator_t
use SAT_mod,                   only : sbp_SAT_penalty_two_block
implicit none

contains

  subroutine calc_div(div, inx, iny, domain, diff_opx, diff_opy)
    type  (multi_grid_field_t),      intent(inout) :: div
    type  (multi_grid_field_t),      intent(in)    :: inx, iny
    type  (multi_domain_t),          intent(in)    :: domain
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy

    type(multi_grid_field_t)   :: gx_buff
    type(multi_grid_field_t)   :: gy_buff
    integer(kind=8) :: i, j, n, m

    call gx_buff%init(domain)
    call gy_buff%init(domain)

    do n = 1, domain%num_sub_x
      do m = 1, domain%num_sub_y
        call diff_opx%apply(gx_buff%subfields(n, m), inx%subfields(n, m), domain%subdomains(n, m), 'x')
        call diff_opy%apply(gy_buff%subfields(n, m), iny%subfields(n, m), domain%subdomains(n, m), 'y')
      end do
    end do

    call sbp_SAT_penalty_two_block(gx_buff, inx, 'x', domain, diff_opx%name)
    call sbp_SAT_penalty_two_block(gy_buff, iny, 'y', domain, diff_opy%name)
    call div%assign(1.0_8, gx_buff, 1.0_8, gy_buff, domain)

  end subroutine

end module div_mod
