module SAT_mod
use domain_mod,           only : domain_t
use field_mod,            only : field_t
use multi_domain_mod,     only : multi_domain_t
use multi_grid_field_mod, only : multi_grid_field_t
use interpolation_mod,    only : interp_1d_sbp21_2to1_ratio, interp_1d_sbp42_2to1_ratio, interp_identity
implicit none

contains

subroutine sbp_SAT_penalty_two_block(tend, in, direction, domains, diff_method)
  !penalty topics for block docking interfaces and periodicals for the first-order differentiation operator
  type(multi_grid_field_t), intent(inout) :: tend
  type(multi_grid_field_t), intent(in)    :: in
  type(multi_domain_t),     intent(in)    :: domains
  character(len=1),         intent(in)    :: direction
  character(len=7),         intent(in)    :: diff_method

  type(field_t)   :: layer_ls, layer_le, layer_rs, layer_re, interp_layer_ls, interp_layer_le, interp_layer_rs, interp_layer_re
  integer(kind=4) :: n, m, k
  real(kind=8)    :: df, h0

  if (diff_method == 'sbp21_1') then
    h0 = 1.0_8 / 2.0_8
  else if (diff_method == 'sbp42_1') then
    h0 = 17.0_8 / 48.0_8
  end if

  if (direction == 'x') then
    do n = 2, domains%num_sub_x
      do m = 1, domains%num_sub_y
        !n-m-th subdomain
        call layer_ls%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 0, 0)
        call layer_le%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 0, 0)
        call layer_rs%init(domains%subdomains(n, m)%js, domains%subdomains(n, m)%je, 0, 0)
        call layer_re%init(domains%subdomains(n, m)%js, domains%subdomains(n, m)%je, 0, 0)

        do k = layer_ls%is, layer_ls%ie
          layer_ls%f(k, 0) = in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%is, k)
          layer_le%f(k, 0) = in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie, k)
        end do

        do k = layer_rs%is, layer_rs%ie
          layer_rs%f(k, 0) = in%subfields(n, m)%f(in%subfields(n, m)%is, k)
          layer_re%f(k, 0) = in%subfields(n, m)%f(in%subfields(n, m)%ie, k)
        end do

        if(domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
          !blocks n and n-1 size 1:1

          do k = layer_ls%is, layer_ls%ie
            df = (layer_re%f(k, 0) - layer_ls%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
            tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) - (df) / 2.0_8

            df = (layer_le%f(k, 0) - layer_rs%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
            tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) - (df) / 2.0_8

            df = (layer_le%f(k, 0) - layer_rs%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
            tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) - (df) / 2.0_8

            df = (layer_re%f(k, 0) - layer_ls%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
            tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) - (df) / 2.0_8
          end do

        else

          call layer_ls%create_similar(interp_layer_re)
          call layer_le%create_similar(interp_layer_rs)
          call layer_rs%create_similar(interp_layer_le)
          call layer_re%create_similar(interp_layer_ls)

          if (diff_method == 'sbp21_1') then
            if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
              call interp_1d_sbp21_2to1_ratio(layer_ls, interp_layer_ls, 'coarse2fine')
              call interp_1d_sbp21_2to1_ratio(layer_le, interp_layer_le, 'coarse2fine')
              call interp_1d_sbp21_2to1_ratio(layer_rs, interp_layer_rs, 'fine2coarse')
              call interp_1d_sbp21_2to1_ratio(layer_re, interp_layer_re, 'fine2coarse')
            else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
              call interp_1d_sbp21_2to1_ratio(layer_ls, interp_layer_ls, 'fine2coarse')
              call interp_1d_sbp21_2to1_ratio(layer_le, interp_layer_le, 'fine2coarse')
              call interp_1d_sbp21_2to1_ratio(layer_rs, interp_layer_rs, 'coarse2fine')
              call interp_1d_sbp21_2to1_ratio(layer_re, interp_layer_re, 'coarse2fine')
            end if
          else if (diff_method == 'sbp42_1') then
            if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
              call interp_1d_sbp42_2to1_ratio(layer_ls, interp_layer_ls, 'coarse2fine')
              call interp_1d_sbp42_2to1_ratio(layer_le, interp_layer_le, 'coarse2fine')
              call interp_1d_sbp42_2to1_ratio(layer_rs, interp_layer_rs, 'fine2coarse')
              call interp_1d_sbp42_2to1_ratio(layer_re, interp_layer_re, 'fine2coarse')
            else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
              call interp_1d_sbp42_2to1_ratio(layer_ls, interp_layer_ls, 'fine2coarse')
              call interp_1d_sbp42_2to1_ratio(layer_le, interp_layer_le, 'fine2coarse')
              call interp_1d_sbp42_2to1_ratio(layer_rs, interp_layer_rs, 'coarse2fine')
              call interp_1d_sbp42_2to1_ratio(layer_re, interp_layer_re, 'coarse2fine')
            end if
          else
            exit
          end if

          do k = layer_ls%is, layer_ls%ie
            df = (interp_layer_re%f(k, 0) - layer_ls%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
            tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) - (df) / 2.0_8

            df = (layer_le%f(k, 0) - interp_layer_rs%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
            tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) - (df) / 2.0_8
          end do

          do k = layer_rs%is, layer_rs%ie
            df = (interp_layer_le%f(k, 0) - layer_rs%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
            tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) - (df) / 2.0_8

            df = (layer_re%f(k, 0) - interp_layer_ls%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
            tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) - (df) / 2.0_8
          end do

        end if
        !end n-m-sobdmain
      end do
    end do
  end if

end subroutine sbp_SAT_penalty_two_block



subroutine sbp_SAT_penalty_two_block_diffusion(tend, in, domains, coefs, diff_method)
  !penalty topics for block docking interfaces and periodicals
  !for a second-order differentiation operator taking into account coefficients for diffusion
  type(multi_grid_field_t),  intent(inout) :: tend
  type(multi_grid_field_t),  intent(in)    :: in
  type(multi_domain_t),      intent(in)    :: domains
  real(kind=8), allocatable, intent(in)    :: coefs(:, :)
  character(len=5),          intent(in)    :: diff_method

  type(field_t)   :: layer_ls, layer_le, layer_rs, layer_re, interp_layer_ls, interp_layer_le, interp_layer_rs, interp_layer_re
  integer(kind=4) :: n, m, k
  real(kind=8)    :: df, h0

  h0 = 1.0_8 / 2.0_8
  if (diff_method == 'sbp21_2') then
    h0 = 1.0_8 / 2.0_8
  else if (diff_method == 'sbp42_2') then
    h0 = 17.0_8 / 48.0_8
  end if

  do n = 2, domains%num_sub_x
    do m = 1, domains%num_sub_y
      !n-m-th subdomain
      call layer_ls%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 0, 0)
      call layer_le%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 0, 0)
      call layer_rs%init(domains%subdomains(n, m)%js, domains%subdomains(n, m)%je, 0, 0)
      call layer_re%init(domains%subdomains(n, m)%js, domains%subdomains(n, m)%je, 0, 0)

      do k = layer_ls%is, layer_ls%ie
        layer_ls%f(k, 0) = (3.0_8 * in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%is, k) - 4.0_8 * in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%is + 1, k) + in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%is + 2, k)) / 2.0_8 / domains%subdomains(n - 1, m)%dx
        layer_le%f(k, 0) = (3.0_8 * in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie, k) - 4.0_8 * in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie - 1, k) + in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie - 2, k)) / 2.0_8 / domains%subdomains(n - 1, m)%dx
      end do

      do k = layer_rs%is, layer_rs%ie
        layer_rs%f(k, 0) = (3.0_8 * in%subfields(n, m)%f(in%subfields(n, m)%is, k) - 4.0_8 * in%subfields(n, m)%f(in%subfields(n, m)%is + 1, k) + in%subfields(n, m)%f(in%subfields(n, m)%is + 2, k)) / 2.0_8 / domains%subdomains(n, m)%dx
        layer_re%f(k, 0) = (3.0_8 * in%subfields(n, m)%f(in%subfields(n, m)%ie, k) - 4.0_8 * in%subfields(n, m)%f(in%subfields(n, m)%ie - 1, k) + in%subfields(n, m)%f(in%subfields(n, m)%ie - 2, k)) / 2.0_8 / domains%subdomains(n, m)%dx
      end do

      if(domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
        !blocks n and n-1 size 1:1

        do k = layer_ls%is, layer_ls%ie
          df = (coefs(n, m) * layer_re%f(k, 0) + coefs(n - 1, m) * layer_ls%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) - (df) / 2.0_8

          df = (coefs(n - 1, m) * layer_rs%f(k, 0) + coefs(n, m) * layer_le%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) - (df) / 2.0_8

          df = (coefs(n - 1, m) * layer_rs%f(k, 0) + coefs(n, m) * layer_le%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
          tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) - (df) / 2.0_8

          df = (coefs(n, m) * layer_re%f(k, 0) + coefs(n - 1, m) * layer_ls%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
          tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) - (df) / 2.0_8
        end do

      else

        call layer_ls%create_similar(interp_layer_re)
        call layer_le%create_similar(interp_layer_rs)
        call layer_rs%create_similar(interp_layer_le)
        call layer_re%create_similar(interp_layer_ls)

          if (diff_method == 'sbp21_2') then
            if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
              call interp_1d_sbp21_2to1_ratio(layer_ls, interp_layer_ls, 'coarse2fine')
              call interp_1d_sbp21_2to1_ratio(layer_le, interp_layer_le, 'coarse2fine')
              call interp_1d_sbp21_2to1_ratio(layer_rs, interp_layer_rs, 'fine2coarse')
              call interp_1d_sbp21_2to1_ratio(layer_re, interp_layer_re, 'fine2coarse')
            else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
              call interp_1d_sbp21_2to1_ratio(layer_ls, interp_layer_ls, 'fine2coarse')
              call interp_1d_sbp21_2to1_ratio(layer_le, interp_layer_le, 'fine2coarse')
              call interp_1d_sbp21_2to1_ratio(layer_rs, interp_layer_rs, 'coarse2fine')
              call interp_1d_sbp21_2to1_ratio(layer_re, interp_layer_re, 'coarse2fine')
            end if
          else if (diff_method == 'sbp42_2') then
            if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
              call interp_1d_sbp42_2to1_ratio(layer_ls, interp_layer_ls, 'coarse2fine')
              call interp_1d_sbp42_2to1_ratio(layer_le, interp_layer_le, 'coarse2fine')
              call interp_1d_sbp42_2to1_ratio(layer_rs, interp_layer_rs, 'fine2coarse')
              call interp_1d_sbp42_2to1_ratio(layer_re, interp_layer_re, 'fine2coarse')
            else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
              call interp_1d_sbp42_2to1_ratio(layer_ls, interp_layer_ls, 'fine2coarse')
              call interp_1d_sbp42_2to1_ratio(layer_le, interp_layer_le, 'fine2coarse')
              call interp_1d_sbp42_2to1_ratio(layer_rs, interp_layer_rs, 'coarse2fine')
              call interp_1d_sbp42_2to1_ratio(layer_re, interp_layer_re, 'coarse2fine')
            end if
          else
            exit
          end if

        do k = layer_ls%is, layer_ls%ie
          df = (coefs(n, m) * interp_layer_re%f(k, 0) + coefs(n - 1, m) * layer_ls%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%is, k) - (df) / 2.0_8

          df = (coefs(n - 1, m) * layer_rs%f(k, 0) + coefs(n, m) * interp_layer_le%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) - (df) / 2.0_8
        end do

        do k = layer_rs%is, layer_rs%ie
          df = (coefs(n - 1, m) * interp_layer_rs%f(k, 0) + coefs(n, m) * layer_le%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
          tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) - (df) / 2.0_8

          df = (coefs(n, m) * layer_re%f(k, 0) + coefs(n - 1, m) * interp_layer_ls%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
          tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%ie, k) - (df) / 2.0_8
        end do

      end if
      !end n-m-sobdmain
    end do
  end do

end subroutine sbp_SAT_penalty_two_block_diffusion

end module SAT_mod
