module SAT_mod
use domain_mod,           only : domain_t
use field_mod,            only : field_t
use multi_domain_mod,     only : multi_domain_t
use multi_grid_field_mod, only : multi_grid_field_t
use interpolation_mod,    only : interp_1d_sbp21_2to1_ratio, interp_1d_sbp42_2to1_ratio, identity
implicit none

contains

subroutine sbp_SAT_penalty_two_block(tend, in, direction, domains, diff_method)

  type(multi_grid_field_t), intent(inout) :: tend
  type(multi_grid_field_t), intent(in)    :: in
  type(multi_domain_t),     intent(in)    :: domains
  character(len=1),         intent(in)    :: direction
  character(len=5),         intent(in)    :: diff_method

  type(field_t)   :: layer_ls, layer_le, layer_rs, layer_re
  integer(kind=4) :: n, m, k
  real(kind=8)    :: df, h0

  if (direction == 'x') then
    do n = 2, domains%num_sub_x
      do m = 1, domains%num_sub_y
        !n-m-th subdomain
        if(domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
          !blocks n and n-1 size 1:1
          h0 = 1.0_8 / 2.0_8
          call layer_ls%init(0, 0, domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je)
          call layer_le%init(0, 0, domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je)
          call layer_rs%init(0, 0, domains%subdomains(n, m)%js, domains%subdomains(n, m)%je)
          call layer_re%init(0, 0, domains%subdomains(n, m)%js, domains%subdomains(n, m)%je)

          do k = layer_ls%js, layer_ls%je
            layer_ls%f(0, k) = in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%is, k)
            layer_le%f(0, k) = in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie, k)
            layer_rs%f(0, k) = in%subfields(n, m)%f(in%subfields(n, m)%is, k)
            layer_re%f(0, k) = in%subfields(n, m)%f(in%subfields(n, m)%ie, k)
          end do

          do k = layer_ls%js, layer_ls%je
            df = (layer_re%f(0, k) - layer_ls%f(0, k)) / (domains%subdomains(n - 1, m)%dx * h0)
            tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) - (df) / 2.0_8

            df = (layer_le%f(0, k) - layer_rs%f(0, k)) / (domains%subdomains(n - 1, m)%dx * h0)
            tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) - (df) / 2.0_8

            df = (layer_le%f(0, k) - layer_rs%f(0, k)) / (domains%subdomains(n, m)%dx * h0)
            tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) - (df) / 2.0_8

            df = (layer_re%f(0, k) - layer_ls%f(0, k)) / (domains%subdomains(n, m)%dx * h0)
            tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) - (df) / 2.0_8
          end do

        end if
        !end n-m-sobdmain
      end do
    end do
  end if

end subroutine sbp_SAT_penalty_two_block

end module SAT_mod
