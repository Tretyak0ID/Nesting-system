module curl_mod
use field_mod,                 only: field_t
use domain_mod,                only: domain_t
use multi_domain_mod,          only : multi_domain_t
use multi_grid_field_mod,      only : multi_grid_field_t
use differential_operator_mod, only : differential_operator_t
use SAT_mod,                   only : sbp_SAT_penalty_two_block
implicit none

contains

  subroutine calc_curl(curl, inx, iny, domain, diff_opx, diff_opy)
    type  (multi_grid_field_t),      intent(inout) :: curl
    type  (multi_grid_field_t),      intent(in)    :: inx, iny
    type  (multi_domain_t),          intent(in)    :: domain
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy
    type(multi_grid_field_t)   :: gx_buff
    type(multi_grid_field_t)   :: gy_buff
    integer(kind=8)            :: i, j, n, m

    call gx_buff%init(domain)
    call gy_buff%init(domain)

    do n = 1, domain%num_sub_x
      do m = 1, domain%num_sub_y
        call diff_opx%apply(gy_buff%subfields(n, m), inx%subfields(n, m), domain%subdomains(n, m), 'y')
        call diff_opy%apply(gx_buff%subfields(n, m), iny%subfields(n, m), domain%subdomains(n, m), 'x')
      end do
    end do

    call sbp_SAT_penalty_two_block(gx_buff, iny, 'x', domain, 'sbp21')
    call sbp_SAT_penalty_two_block(gy_buff, inx, 'y', domain, 'sbp21')

    do n = 1, domain%num_sub_x
      do m = 1, domain%num_sub_y
        do i = domain%subdomains(n, m)%is, domain%subdomains(n, m)%ie
          do j = domain%subdomains(n, m)%js, domain%subdomains(n, m)%je
            curl%subfields(n, m)%f(i, j) = gx_buff%subfields(n, m)%f(i, j) - gy_buff%subfields(n, m)%f(i, j)
          end do
        end do
      end do
    end do

  end subroutine

end module curl_mod
