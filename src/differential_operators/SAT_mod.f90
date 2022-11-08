module SAT_mod
use domain_mod,           only : domain_t
use field_mod,            only : field_t
use multi_domain_mod,     only : multi_domain_t
use multi_grid_field_mod, only : multi_grid_field_t
use interpolation_mod,    only : interp_1d_sbp21_2to1_ratio, interp_1d_sbp42_2to1_ratio, identity
implicit none

contains

subroutine sbp_SAT_penalty_two_block(tend, in, direction, multi_domain, diff_method)

  type(multi_grid_field_t), intent(inout) :: tend
  type(multi_grid_field_t), intent(in)    :: in
  type(multi_domain_t),     intent(in)    :: multi_domain
  character(len=1),         intent(in)    :: direction
  character(len=5),         intent(in)    :: diff_method

  type(field_t)   :: ff, df, layer_ls, layer_le, layer_rs, layer_re
  integer(kind=4) :: i, j, k
  real(kind=8)    :: h0

  if (diff_method == 'sbp21' .and. direction == 'x') then
    h0 = 1.0_8 / 2.0_8
    do i = 1, multi_domain%num_sub_x !iteration on subdomains
      do j = 1, multi_domain%num_sub_y

        if (i > 1) then
          if(multi_domain%subdomains(i - 1, j)%ny == multi_domain%subdomains(i, j)%ny) then !left : right <=> 1:1

            call layer_ls%init(in%subfields(i - 1, j)%js, in%subfields(i - 1, j)%je, 0, 0)
            call layer_re%init(in%subfields(i, j)%js, in%subfields(i, j)%je, 0, 0)
            call layer_rs%init(in%subfields(i, j)%js, in%subfields(i, j)%je, 0, 0)
            call layer_le%init(in%subfields(i - 1, j)%js, in%subfields(i - 1, j)%je, 0, 0)

            do k = layer_ls%js, layer_ls%je
              layer_ls%f(k, 0) = in%subfields(i - 1, j)%f(in%subfields(i - 1, j)%ie, k)
              layer_le%f(k, 0) = in%subfields(i - 1, j)%f(in%subfields(i - 1, j)%is, k)
            end do

            do k = layer_rs%js, layer_rs%je
              layer_rs%f(k, 0) = in%subfields(i, j)%f(in%subfields(i, j)%is, k)
              layer_re%f(k, 0) = in%subfields(i, j)%f(in%subfields(i, j)%ie, k)
            end do

            !right y interface
            do k = in%subfields(i - 1, j)%js, in%subfields(i - 1, j)%je
              tend%subfields(i - 1, j)%f(tend%subfields(i - 1, j)%ie, k) = tend%subfields(i - 1, j)%f(tend%subfields(i - 1, j)%ie, k) \

            end do

            !left y interface
            do k = in%subfields(i, j)%js, in%subfields(i, j)%je
              tend%subfields(i, j)%f(tend%subfields(i, j)%is, k) = tend%subfields(i, j)%f(tend%subfields(i, j)%is, k)
            end do

            !right y interface
            do k = in%subfields(i - 1, j)%js, in%subfields(i - 1, j)%je
              tend%subfields(i - 1, j)%f(tend%subfields(i, j)%is, k) = tend%subfields(i - 1, j)%f(tend%subfields(i, j)%is, k)
            end do

            !left y interface
            do k = in%subfields(i, j)%js, in%subfields(i, j)%je
              tend%subfields(i, j)%f(tend%subfields(i, j)%ie, k) = tend%subfields(i, j)%f(tend%subfields(i, j)%ie, k)
            end do

          end if
        end if

      end do
    end do

  else if (diff_method == 'sbp42' .and. direction == 'x') then
    continue
  end if

end subroutine sbp_SAT_penalty_two_block

end module SAT_mod
