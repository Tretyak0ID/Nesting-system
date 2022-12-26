module laplacian_mod
use multi_domain_mod,          only : multi_domain_t
use multi_grid_field_mod,      only : multi_grid_field_t
use differential_operator_mod, only : differential_operator_t
use SAT_mod,                   only : sbp_SAT_penalty_two_block_diffusion
implicit none

contains

  subroutine calc_laplacian(out, in, multi_domain, coefs, diff2_op)
    type (multi_grid_field_t),      intent(inout) :: out
    type (multi_grid_field_t),      intent(in)    :: in
    type (multi_domain_t),          intent(in)    :: multi_domain
    real (kind=8), allocatable,     intent(in)    :: coefs(:, :)
    class(differential_operator_t), intent(in)    :: diff2_op

    type(multi_grid_field_t)  :: dinx, diny
    integer(kind=4) :: n, m

    call dinx%init(multi_domain)
    call diny%init(multi_domain)

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        call diff2_op%apply(dinx%subfields(n, m), in%subfields(n, m), multi_domain%subdomains(n, m), 'x')
        call diff2_op%apply(diny%subfields(n, m), in%subfields(n, m), multi_domain%subdomains(n, m), 'y')
        call out%subfields(n, m)%assign(coefs(n, m), dinx%subfields(n, m), coefs(n, m), diny%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

    !call out%assign(coefs(1, 1), dinx, coefs(1, 1), diny, multi_domain)

    call sbp_SAT_penalty_two_block_diffusion(out, in, multi_domain, coefs, diff2_op%name)

  end subroutine calc_laplacian

end module laplacian_mod
