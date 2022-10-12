module div_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
use differential_operator_mod, only: differential_operator_t
implicit none

contains

  subroutine calc_div(div, in, domain, diff_opx, diff_opy)
    type  (field_t),                 intent(inout) :: div
    type  (field_t),                 intent(in)    :: in
    type  (domain_t),                intent(in)    :: domain
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy
    type(field_t)   :: gx_buff
    type(field_t)   :: gy_buff
    integer(kind=8) :: i, j

    call gx_buff%init_on_domain(domain)
    call gy_buff%init_on_domain(domain)

    call diff_opx%apply(gx_buff, in, domain, 'x')
    call diff_opy%apply(gy_buff, in, domain, 'y')

    do i = 0, domain%nx
      do j = 0, domain%ny
        div%f(i, j) = gx_buff%f(i, j) + gy_buff%f(i, j)
      end do
    end do
  end subroutine

end module div_mod
