module SAT_mod
use domain_mod,           only : domain_t
use field_mod,            only : field_t
use multi_domain_mod,     only : multi_domain_t
use multi_grid_field_mod, only : multi_grid_field_t
use interpolation_mod,    only : interp_identity, interp_MC2order_2to1ratio, interp_MC2order_2to1ratio_periodic, interp_MC4order_2to1ratio, interp_MC4order_2to1ratio_periodic
use boundary_methods_mod, only : apply_sbp21_2_boundary_method, apply_sbp42_2_boundary_method, apply_sbp42_2_boundary_method_second_terms, apply_sbp21_2_boundary_method_second_terms
implicit none

contains

subroutine sbp_SAT_penalty_two_block(tend, in, direction, domains, diff_method)
  !penalty topics for block docking interfaces and periodicals for the first-order differentiation operator
  type(multi_grid_field_t), intent(inout) :: tend
  type(multi_grid_field_t), intent(in)    :: in
  type(multi_domain_t),     intent(in)    :: domains
  character(len=1),         intent(in)    :: direction
  character(len=7),         intent(in)    :: diff_method

  type(field_t) :: s1, eNx !boundary layers
  type(field_t) :: interp_s1, interp_eNx !interpolation buffers for boundary layers
  type(field_t) :: le, rs !intermediate layers
  type(field_t) :: interp_le, interp_rs !interpolation buffers intermediate layers

  integer(kind=4) :: n, m, k, flag = 1
  real(kind=8)    :: df, h0

  if (diff_method == 'sbp21_1') then
    h0 = 1.0_8 / 2.0_8
  else if (diff_method == 'sbp42_1') then
    h0 = 17.0_8 / 48.0_8
  else
    h0 = 0.0_8
  end if

  if(direction == 'x' .and. (h0 /= 0.0_8)) then
    do m = 1, domains%num_sub_y
      !interpolation into intermediate layers
      do n = 2, domains%num_sub_x
        call le%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 0, 0)
        call rs%init(domains%subdomains(n, m)%js, domains%subdomains(n, m)%je, 0, 0)

        call le%create_similar(interp_rs)
        call rs%create_similar(interp_le)

        do k = le%is, le%ie
          le%f(k, 0) = in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie, k)
        end do

        do k = rs%is, rs%ie
          rs%f(k, 0) = in%subfields(n, m)%f(in%subfields(n, m)%is, k)
        end do

        if (diff_method == 'sbp21_1') then
          if (domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_identity(le, interp_le, 'coarse2fine')
            call interp_identity(rs, interp_rs, 'fine2coarse')
          else if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_MC2order_2to1ratio(le, interp_le, 'coarse2fine')
            call interp_MC2order_2to1ratio(rs, interp_rs, 'fine2coarse')
          else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
            call interp_MC2order_2to1ratio(le, interp_le, 'fine2coarse')
            call interp_MC2order_2to1ratio(rs, interp_rs, 'coarse2fine')
          end if
        else if (diff_method == 'sbp42_1') then
          if (domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_identity(le, interp_le, 'coarse2fine')
            call interp_identity(rs, interp_rs, 'fine2coarse')
          else if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_MC4order_2to1ratio(le, interp_le, 'coarse2fine')
            call interp_MC4order_2to1ratio(rs, interp_rs, 'fine2coarse')
          else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
            call interp_MC4order_2to1ratio(le, interp_le, 'fine2coarse')
            call interp_MC4order_2to1ratio(rs, interp_rs, 'coarse2fine')
          end if
        else
          exit
        end if

        do k = le%is, le%ie
          df = (le%f(k, 0) - interp_rs%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) - (df) / 2.0_8
        end do

        do k = rs%is, rs%ie
          df = (interp_le%f(k, 0) - rs%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
          tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) - (df) / 2.0_8
        end do

      end do

      !boundary layers interpolation
      call s1%init(domains%subdomains(1, m)%js, domains%subdomains(1, m)%je, 0, 0)
      call eNx%init(domains%subdomains(domains%num_sub_x, m)%js, domains%subdomains(domains%num_sub_x, m)%je, 0, 0)
      call s1%create_similar(interp_eNx)
      call eNx%create_similar(interp_s1)

      do k = domains%subdomains(1, m)%js, domains%subdomains(1, m)%je
        s1%f(k, 0) = in%subfields(1, m)%f(in%subfields(1, m)%is, k)
      end do
      do k = domains%subdomains(domains%num_sub_x, m)%js, domains%subdomains(domains%num_sub_x, m)%je
        eNx%f(k, 0) = in%subfields(domains%num_sub_x, m)%f(in%subfields(domains%num_sub_x, m)%ie, k)
      end do

      if (diff_method == 'sbp21_1') then
        if (domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
          call interp_identity(s1,   interp_s1,   'coarse2fine')
          call interp_identity(eNx,  interp_eNx,  'fine2coarse')
        else if (2 * domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
          call interp_MC2order_2to1ratio(s1,  interp_s1,  'coarse2fine')
          call interp_MC2order_2to1ratio(eNx, interp_eNx, 'fine2coarse')
        else if (domains%subdomains(1, m)%ny == 2 * domains%subdomains(domains%num_sub_x, m)%ny) then
          call interp_MC2order_2to1ratio(s1,  interp_s1,  'fine2coarse')
          call interp_MC2order_2to1ratio(eNx, interp_eNx, 'coarse2fine')
        end if
      else if (diff_method == 'sbp42_1') then
        if (domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
          call interp_identity(s1,   interp_s1,   'coarse2fine')
          call interp_identity(eNx,  interp_eNx,  'fine2coarse')
        else if (2 * domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
          call interp_MC4order_2to1ratio(s1,  interp_s1,  'coarse2fine')
          call interp_MC4order_2to1ratio(eNx, interp_eNx, 'fine2coarse')
        else if (domains%subdomains(1, m)%ny == 2 * domains%subdomains(domains%num_sub_x, m)%ny) then
          call interp_MC4order_2to1ratio(s1,  interp_s1,  'fine2coarse')
          call interp_MC4order_2to1ratio(eNx, interp_eNx, 'coarse2fine')
        end if
      else
        exit
      end if

      do k = s1%is, s1%ie
        df = (interp_eNx%f(k, 0) - s1%f(k, 0)) / (domains%subdomains(1, m)%dx * h0)
        tend%subfields(1, m)%f(tend%subfields(1, m)%is, k) = tend%subfields(1, m)%f(tend%subfields(1, m)%is, k) - (df) / 2.0_8
      end do

      do k = eNx%is, eNx%ie
        df = (eNx%f(k, 0) - interp_s1%f(k, 0)) / (domains%subdomains(domains%num_sub_x, m)%dx * h0)
        tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie, k) = tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie, k) - (df) / 2.0_8
      end do
    end do

  else if (direction == 'y' .and. (h0 /= 0.0_8)) then

    do n = 1, domains%num_sub_x
      !interpolation into intermediate layers
      do m = 2, domains%num_sub_y
        call le%init(domains%subdomains(n, m - 1)%is, domains%subdomains(n, m - 1)%ie, 0, 0)
        call rs%init(domains%subdomains(n, m)%is, domains%subdomains(n, m)%ie, 0, 0)

        call le%create_similar(interp_rs)
        call rs%create_similar(interp_le)

        do k = le%is, le%ie
          le%f(k, 0) = in%subfields(n, m - 1)%f(k, in%subfields(n, m - 1)%je)
        end do

        do k = rs%is, rs%ie
          rs%f(k, 0) = in%subfields(n, m)%f(k, in%subfields(n, m)%js)
        end do

        if (diff_method == 'sbp21_1') then
          if (domains%subdomains(n, m - 1)%ny == domains%subdomains(n, m)%ny) then
            call interp_identity(le, interp_le, 'coarse2fine')
            call interp_identity(rs, interp_rs, 'fine2coarse')
          else if (2 * domains%subdomains(n, m - 1)%ny == domains%subdomains(n, m)%ny) then
            call interp_MC2order_2to1ratio(le, interp_le, 'coarse2fine')
            call interp_MC2order_2to1ratio(rs, interp_rs, 'fine2coarse')
          else if (domains%subdomains(n, m - 1)%ny == 2 * domains%subdomains(n, m)%ny) then
            call interp_MC2order_2to1ratio(le, interp_le, 'fine2coarse')
            call interp_MC2order_2to1ratio(rs, interp_rs, 'coarse2fine')
          end if
        else if (diff_method == 'sbp42_1') then
          if (domains%subdomains(n, m - 1)%ny == domains%subdomains(n, m)%ny) then
            call interp_identity(le, interp_le, 'coarse2fine')
            call interp_identity(rs, interp_rs, 'fine2coarse')
          else if (2 * domains%subdomains(n, m - 1)%ny == domains%subdomains(n, m)%ny) then
            call interp_MC4order_2to1ratio(le, interp_le, 'coarse2fine')
            call interp_MC4order_2to1ratio(rs, interp_rs, 'fine2coarse')
          else if (domains%subdomains(n, m - 1)%ny == 2 * domains%subdomains(n, m)%ny) then
            call interp_MC4order_2to1ratio(le, interp_le, 'fine2coarse')
            call interp_MC4order_2to1ratio(rs, interp_rs, 'coarse2fine')
          end if
        else
          exit
        end if

        do k = le%is, le%ie
          df = (le%f(k, 0) - interp_rs%f(k, 0)) / (domains%subdomains(n, m - 1)%dy * h0)
          tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je) = tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je) - (df) / 2.0_8
        end do

        do k = rs%is, rs%ie
          df = (interp_le%f(k, 0) - rs%f(k, 0)) / (domains%subdomains(n, m)%dy * h0)
          tend%subfields(n, m)%f(k, tend%subfields(n, m)%js) = tend%subfields(n, m)%f(k, tend%subfields(n, m)%js) - (df) / 2.0_8
        end do

      end do

      !boundary layers interpolation
      call s1%init(domains%subdomains(n, 1)%is, domains%subdomains(n, 1)%ie, 0, 0)
      call eNx%init(domains%subdomains(n, domains%num_sub_y)%is, domains%subdomains(n, domains%num_sub_y)%ie, 0, 0)
      call s1%create_similar(interp_eNx)
      call eNx%create_similar(interp_s1)

      do k = domains%subdomains(n, 1)%is, domains%subdomains(n, 1)%ie
        s1%f(k, 0) = in%subfields(n, 1)%f(k, in%subfields(n, 1)%js)
      end do
      do k = domains%subdomains(n, domains%num_sub_y)%is, domains%subdomains(n, domains%num_sub_y)%ie
        eNx%f(k, 0) = in%subfields(n, domains%num_sub_y)%f(k, in%subfields(n, domains%num_sub_y)%je)
      end do

      if (diff_method == 'sbp21_1') then
        if (domains%subdomains(n, 1)%ny == domains%subdomains(n, domains%num_sub_y)%ny) then
          call interp_identity(s1,   interp_s1,   'coarse2fine')
          call interp_identity(eNx,  interp_eNx,  'fine2coarse')
        else if (2 * domains%subdomains(n, 1)%ny == domains%subdomains(n, domains%num_sub_y)%ny) then
          call interp_MC2order_2to1ratio(s1,  interp_s1,  'coarse2fine')
          call interp_MC2order_2to1ratio(eNx, interp_eNx, 'fine2coarse')
        else if (domains%subdomains(n, 1)%ny == 2 * domains%subdomains(n, domains%num_sub_y)%ny) then
          call interp_MC2order_2to1ratio(s1,  interp_s1,  'fine2coarse')
          call interp_MC2order_2to1ratio(eNx, interp_eNx, 'coarse2fine')
        end if
      else if (diff_method == 'sbp42_1') then
        if (domains%subdomains(n, 1)%ny == domains%subdomains(n, domains%num_sub_y)%ny) then
          call interp_identity(s1,   interp_s1,   'coarse2fine')
          call interp_identity(eNx,  interp_eNx,  'fine2coarse')
        else if (2 * domains%subdomains(n, 1)%ny == domains%subdomains(n, domains%num_sub_y)%ny) then
          call interp_MC4order_2to1ratio(s1,  interp_s1,  'coarse2fine')
          call interp_MC4order_2to1ratio(eNx, interp_eNx, 'fine2coarse')
        else if (domains%subdomains(n, 1)%ny == 2 * domains%subdomains(n, domains%num_sub_y)%ny) then
          call interp_MC4order_2to1ratio(s1,  interp_s1,  'fine2coarse')
          call interp_MC4order_2to1ratio(eNx, interp_eNx, 'coarse2fine')
        end if
      else
        exit
      end if

      do k = s1%is, s1%ie
        df = (interp_eNx%f(k, 0) - s1%f(k, 0)) / (domains%subdomains(n, 1)%dy * h0)
        tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js) = tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js) - (df) / 2.0_8
      end do

      do k = eNx%is, eNx%ie
        df = (eNx%f(k, 0) - interp_s1%f(k, 0)) / (domains%subdomains(n, domains%num_sub_y)%dy * h0)
        tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je) = tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je) - (df) / 2.0_8
      end do
    end do
  end if

end subroutine sbp_SAT_penalty_two_block


!-----------------------------------------
!------------D-I-F-F-U-S-I-O-N------------
!-----------------------------------------
subroutine sbp_SAT_penalty_two_block_diffusion(tend, in, domains, coefs, direction, diff_method)
  !input parameters
  type(multi_grid_field_t),  intent(inout) :: tend
  type(multi_grid_field_t),  intent(in)    :: in
  type(multi_domain_t),      intent(in)    :: domains
  real(kind=8), allocatable, intent(in)    :: coefs(:, :)
  character(len=1),          intent(in)    :: direction
  character(len=7),          intent(in)    :: diff_method

  !routines parameter
  type(field_t) :: ft_ls, ft_le, ft_rs, ft_re, st_ls, st_le, st_rs, st_re !ft/st = first/second term, l/r = left/right, s/e = start/end
  type(field_t) :: ft_ils, ft_ile, ft_irs, ft_ire, st_ils, st_ile, st_irs, st_ire, stl, str, stlb, strb !i = interp field
  integer(kind=4) :: n, m, k, flag = 1
  real(kind=8)    :: ft_df, h0, q

  if      (diff_method == 'sbp21_2') then
    h0 =  1.0_8 / 2.0_8
    q  =  2.5_8
  else if (diff_method == 'sbp42_2') then
    h0 = 17.0_8 / 48.0_8
    q  = 4.0_8
  else
    h0 = 0.0_8
    q  = 0.0_8
  end if

  if (direction == 'x') then
    do m = 1, domains%num_sub_y
    do n = 2, domains%num_sub_x

        !create
        call ft_le%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 0, 0)
        call ft_rs%init(domains%subdomains(n, m)%js,     domains%subdomains(n, m)%je, 0, 0)
        call st_le%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 0, 0)
        call st_rs%init(domains%subdomains(n, m)%js,     domains%subdomains(n, m)%je, 0, 0)
        call ft_le%create_similar(ft_irs)
        call ft_rs%create_similar(ft_ile)
        call st_le%create_similar(st_irs)
        call st_rs%create_similar(st_ile)
        call stl%init(domains%subdomains(n - 1, m)%js, domains%subdomains(n - 1, m)%je, 1, 4)
        call str%init(domains%subdomains(n, m)%js, domains%subdomains(n, m)%je, 1, 4)

        !feeling
        do k = ft_le%is, ft_le%ie
          ft_le%f(k, 0) = in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie, k)
          st_le%f(k, 0) = in%subfields(n - 1, m)%f(in%subfields(n - 1, m)%ie, k)
        end do
        do k = ft_rs%is, ft_rs%ie
          ft_rs%f(k, 0) = in%subfields(n, m)%f(in%subfields(n, m)%is, k)
          st_rs%f(k, 0) = in%subfields(n, m)%f(in%subfields(n, m)%is, k)
        end do

        !boundary method
        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method(ft_le, in%subfields(n - 1, m), domains%subdomains(n - 1, m), direction, 'e')
          call apply_sbp21_2_boundary_method(ft_rs, in%subfields(n, m), domains%subdomains(n, m), direction, 's')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method(ft_le, in%subfields(n - 1, m), domains%subdomains(n - 1, m), direction, 'e')
          call apply_sbp42_2_boundary_method(ft_rs, in%subfields(n, m), domains%subdomains(n, m), direction, 's')
        end if

        !interpolation
        if (diff_method == 'sbp21_2') then
          if (domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_identity(st_le, st_ile, 'coarse2fine')
            call interp_identity(st_rs, st_irs, 'fine2coarse')
            call interp_identity(ft_le, ft_ile, 'coarse2fine')
            call interp_identity(ft_rs, ft_irs, 'fine2coarse')
          else if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_MC2order_2to1ratio(st_le, st_ile, 'coarse2fine')
            call interp_MC2order_2to1ratio(st_rs, st_irs, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_le, ft_ile, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_rs, ft_irs, 'fine2coarse')
          else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
            call interp_MC2order_2to1ratio(st_le, st_ile, 'fine2coarse')
            call interp_MC2order_2to1ratio(st_rs, st_irs, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_le, ft_ile, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_rs, ft_irs, 'coarse2fine')
          end if
        else if (diff_method == 'sbp42_2') then
          if (domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_identity(st_le, st_ile, 'coarse2fine')
            call interp_identity(st_rs, st_irs, 'fine2coarse')
            call interp_identity(ft_le, ft_ile, 'coarse2fine')
            call interp_identity(ft_rs, ft_irs, 'fine2coarse')
          else if (2 * domains%subdomains(n - 1, m)%ny == domains%subdomains(n, m)%ny) then
            call interp_MC4order_2to1ratio(st_le, st_ile, 'coarse2fine')
            call interp_MC4order_2to1ratio(st_rs, st_irs, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_le, ft_ile, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_rs, ft_irs, 'fine2coarse')
          else if (domains%subdomains(n - 1, m)%ny == 2 * domains%subdomains(n, m)%ny) then
            call interp_MC4order_2to1ratio(st_le, st_ile, 'fine2coarse')
            call interp_MC4order_2to1ratio(st_rs, st_irs, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_le, ft_ile, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_rs, ft_irs, 'coarse2fine')
          end if
        else
          exit
        end if

        do k = ft_le%is, ft_le%ie
          st_le%f(k, 0) = st_le%f(k, 0) - st_irs%f(k, 0)
        end do
        do k = ft_rs%is, ft_rs%ie
          st_rs%f(k, 0) = st_rs%f(k, 0) - st_ile%f(k, 0)
        end do

        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method_second_terms(stl, st_le, domains%subdomains(n - 1, m), 'x', 'e')
          call apply_sbp21_2_boundary_method_second_terms(str, st_rs, domains%subdomains(n, m), 'x', 's')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method_second_terms(stl, st_le, domains%subdomains(n - 1, m), 'x', 'e')
          call apply_sbp42_2_boundary_method_second_terms(str, st_rs, domains%subdomains(n, m), 'x', 's')
        end if

!apply xi
        do k = ft_le%is, ft_le%ie
          ft_df = (coefs(n - 1, m) * ft_le%f(k, 0) + coefs(n, m) * ft_irs%f(k, 0) + (coefs(n - 1, m)+coefs(n, m)) * q / domains%subdomains(n - 1, m)%dx / 2.0_8 * st_le%f(k, 0)) / (domains%subdomains(n - 1, m)%dx * h0)
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie, k) - ft_df / 2.0_8 + coefs(n - 1, m) * stl%f(k, 1) / 2.0_8
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie - 1, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie - 1, k) + coefs(n - 1, m) * stl%f(k, 2) / 2.0_8
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie - 2, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie - 2, k) + coefs(n - 1, m) * stl%f(k, 3) / 2.0_8
          tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie - 3, k) = tend%subfields(n - 1, m)%f(tend%subfields(n - 1, m)%ie - 3, k) + coefs(n - 1, m) * stl%f(k, 4) / 2.0_8
        end do

        do k = ft_rs%is, ft_rs%ie
          ft_df = (coefs(n, m) * ft_rs%f(k, 0) + coefs(n - 1, m) * ft_ile%f(k, 0) + (coefs(n - 1, m)+coefs(n, m)) * q / domains%subdomains(n, m)%dx / 2.0_8 * st_rs%f(k, 0)) / (domains%subdomains(n, m)%dx * h0)
          tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is, k) - ft_df / 2.0_8 + coefs(n, m) * str%f(k, 1) / 2.0_8
          tend%subfields(n, m)%f(tend%subfields(n, m)%is + 1, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is + 1, k) + coefs(n, m) * str%f(k, 2) / 2.0_8
          tend%subfields(n, m)%f(tend%subfields(n, m)%is + 2, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is + 2, k) + coefs(n, m) * str%f(k, 3) / 2.0_8
          tend%subfields(n, m)%f(tend%subfields(n, m)%is + 3, k) = tend%subfields(n, m)%f(tend%subfields(n, m)%is + 3, k) + coefs(n, m) * str%f(k, 4) / 2.0_8
        end do
    end do
    !bound
        call ft_ls%init(domains%subdomains(1, m)%js, domains%subdomains(1, m)%je, 0, 0)
        call ft_re%init(domains%subdomains(domains%num_sub_x, m)%js,     domains%subdomains(domains%num_sub_x, m)%je, 0, 0)
        call st_ls%init(domains%subdomains(1, m)%js, domains%subdomains(1, m)%je, 0, 0)
        call st_re%init(domains%subdomains(domains%num_sub_x, m)%js,     domains%subdomains(domains%num_sub_x, m)%je, 0, 0)
        call ft_ls%create_similar(ft_ire)
        call ft_re%create_similar(ft_ils)
        call st_ls%create_similar(st_ire)
        call st_re%create_similar(st_ils)
        call stlb%init(domains%subdomains(1, m)%js, domains%subdomains(1, m)%je, 1, 4)
        call strb%init(domains%subdomains(domains%num_sub_x, m)%js, domains%subdomains(domains%num_sub_x, m)%je, 1, 4)
        !feeling
        do k = ft_ls%is, ft_ls%ie
          ft_ls%f(k, 0) = in%subfields(1, m)%f(in%subfields(1, m)%is, k)
          st_ls%f(k, 0) = in%subfields(1, m)%f(in%subfields(1, m)%is, k)
        end do
        do k = ft_re%is, ft_re%ie
          ft_re%f(k, 0) = in%subfields(domains%num_sub_x, m)%f(in%subfields(domains%num_sub_x, m)%ie, k)
          st_re%f(k, 0) = in%subfields(domains%num_sub_x, m)%f(in%subfields(domains%num_sub_x, m)%ie, k)
        end do
        !boundary method
        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method(ft_ls, in%subfields(1, m), domains%subdomains(1, m), direction, 's')
          call apply_sbp21_2_boundary_method(ft_re, in%subfields(domains%num_sub_x, m), domains%subdomains(domains%num_sub_x, m), direction, 'e')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method(ft_ls, in%subfields(1, m), domains%subdomains(1, m), direction, 's')
          call apply_sbp42_2_boundary_method(ft_re, in%subfields(domains%num_sub_x, m), domains%subdomains(domains%num_sub_x, m), direction, 'e')
        end if
        !interpolation
        if (diff_method == 'sbp21_2') then
          if (domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
            call interp_identity(st_ls, st_ils, 'coarse2fine')
            call interp_identity(st_re, st_ire, 'fine2coarse')
            call interp_identity(ft_ls, ft_ils, 'coarse2fine')
            call interp_identity(ft_re, ft_ire, 'fine2coarse')
          else if (2 * domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
            call interp_MC2order_2to1ratio(st_ls, st_ils, 'coarse2fine')
            call interp_MC2order_2to1ratio(st_re, st_ire, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_ls, ft_ils, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_re, ft_ire, 'fine2coarse')
          else if (domains%subdomains(1, m)%ny == 2 * domains%subdomains(domains%num_sub_x, m)%ny) then
            call interp_MC2order_2to1ratio(st_ls, st_ils, 'fine2coarse')
            call interp_MC2order_2to1ratio(st_re, st_ire, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_ls, ft_ils, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_re, ft_ire, 'coarse2fine')
          end if
        else if (diff_method == 'sbp42_2') then
          if (domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
            call interp_identity(st_ls, st_ils, 'coarse2fine')
            call interp_identity(st_re, st_ire, 'fine2coarse')
            call interp_identity(ft_ls, ft_ils, 'coarse2fine')
            call interp_identity(ft_re, ft_ire, 'fine2coarse')
          else if (2 * domains%subdomains(1, m)%ny == domains%subdomains(domains%num_sub_x, m)%ny) then
            call interp_MC4order_2to1ratio(st_ls, st_ils, 'coarse2fine')
            call interp_MC4order_2to1ratio(st_re, st_ire, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_ls, ft_ils, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_re, ft_ire, 'fine2coarse')
          else if (domains%subdomains(1, m)%ny == 2 * domains%subdomains(domains%num_sub_x, m)%ny) then
            call interp_MC4order_2to1ratio(st_ls, st_ils, 'fine2coarse')
            call interp_MC4order_2to1ratio(st_re, st_ire, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_ls, ft_ils, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_re, ft_ire, 'coarse2fine')
          end if
        else
          exit
        end if
        do k = ft_ls%is, ft_ls%ie
          st_ls%f(k, 0) = st_ls%f(k, 0) - st_ire%f(k, 0)
        end do
        do k = ft_re%is, ft_re%ie
          st_re%f(k, 0) = st_re%f(k, 0) - st_ils%f(k, 0)
        end do

        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method_second_terms(stlb, st_ls, domains%subdomains(1, m), 'x', 's')
          call apply_sbp21_2_boundary_method_second_terms(strb, st_re, domains%subdomains(domains%num_sub_x, m), 'x', 'e')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method_second_terms(stlb, st_ls, domains%subdomains(1, m), 'x', 's')
          call apply_sbp42_2_boundary_method_second_terms(strb, st_re, domains%subdomains(domains%num_sub_x, m), 'x', 'e')
        end if
!apply xb
        do k = ft_ls%is, ft_ls%ie
          ft_df = (coefs(1, m) * ft_ls%f(k, 0) + coefs(domains%num_sub_x, m) * ft_ire%f(k, 0) + (coefs(1, m)+coefs(domains%num_sub_x, m)) * q / domains%subdomains(1, m)%dx / 2.0_8 * st_ls%f(k, 0)) / (domains%subdomains(1, m)%dx * h0)
          tend%subfields(1, m)%f(tend%subfields(1, m)%is, k) = tend%subfields(1, m)%f(tend%subfields(1, m)%is, k) - ft_df / 2.0_8 + coefs(1, m) * stlb%f(k ,1) / 2.0_8
          tend%subfields(1, m)%f(tend%subfields(1, m)%is + 1, k) = tend%subfields(1, m)%f(tend%subfields(1, m)%is + 1, k) + coefs(1, m) * stlb%f(k, 2) / 2.0_8
          tend%subfields(1, m)%f(tend%subfields(1, m)%is + 2, k) = tend%subfields(1, m)%f(tend%subfields(1, m)%is + 2, k) + coefs(1, m) * stlb%f(k, 3) / 2.0_8
          tend%subfields(1, m)%f(tend%subfields(1, m)%is + 3, k) = tend%subfields(1, m)%f(tend%subfields(1, m)%is + 3, k) + coefs(1, m) * stlb%f(k, 4) / 2.0_8
        end do
        do k = ft_re%is, ft_re%ie
          ft_df = (coefs(domains%num_sub_x, m) * ft_re%f(k, 0) + coefs(1, m) * ft_ils%f(k, 0) + (coefs(1, m)+coefs(domains%num_sub_x, m)) * q / domains%subdomains(domains%num_sub_x, m)%dx / 2.0_8 * st_re%f(k, 0)) / (domains%subdomains(domains%num_sub_x, m)%dx * h0)
          tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie, k) = tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie, k) - ft_df / 2.0_8 + coefs(domains%num_sub_x, m) * strb%f(k, 1) / 2.0_8
          tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie - 1, k) = tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie - 1, k) + coefs(domains%num_sub_x, m) * strb%f(k, 2) / 2.0_8
          tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie - 2, k) = tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie - 2, k) + coefs(domains%num_sub_x, m) * strb%f(k, 3) / 2.0_8
          tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie - 3, k) = tend%subfields(domains%num_sub_x, m)%f(tend%subfields(domains%num_sub_x, m)%ie - 3, k) + coefs(domains%num_sub_x, m) * strb%f(k, 4) / 2.0_8
        end do
    end do
  else if (direction == 'y') then
    do n = 1, domains%num_sub_x
    do m = 2, domains%num_sub_y
        !create
        call ft_le%init(domains%subdomains(n, m - 1)%is, domains%subdomains(n, m - 1)%ie, 0, 0)
        call ft_rs%init(domains%subdomains(n, m)%is,     domains%subdomains(n, m)%ie, 0, 0)
        call st_le%init(domains%subdomains(n, m - 1)%is, domains%subdomains(n, m - 1)%ie, 0, 0)
        call st_rs%init(domains%subdomains(n, m)%is,     domains%subdomains(n, m)%ie, 0, 0)
        call ft_le%create_similar(ft_irs)
        call ft_rs%create_similar(ft_ile)
        call st_le%create_similar(st_irs)
        call st_rs%create_similar(st_ile)
        call stl%init(domains%subdomains(n, m - 1)%is, domains%subdomains(n, m - 1)%ie, 1, 4)
        call str%init(domains%subdomains(n, m)%is, domains%subdomains(n, m)%ie, 1, 4)

        !feeling
        do k = ft_le%is, ft_le%ie
          ft_le%f(k, 0) = in%subfields(n, m - 1)%f(k, in%subfields(n, m - 1)%je)
          st_le%f(k, 0) = in%subfields(n, m - 1)%f(k, in%subfields(n, m - 1)%je)
        end do
        do k = ft_rs%is, ft_rs%ie
          ft_rs%f(k, 0) = in%subfields(n, m)%f(k, in%subfields(n, m)%js)
          st_rs%f(k, 0) = in%subfields(n, m)%f(k, in%subfields(n, m)%js)
        end do
        !boundary method
        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method(ft_le, in%subfields(n, m - 1), domains%subdomains(n, m - 1), direction, 'e')
          call apply_sbp21_2_boundary_method(ft_rs, in%subfields(n, m), domains%subdomains(n, m), direction, 's')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method(ft_le, in%subfields(n, m - 1), domains%subdomains(n, m - 1), direction, 'e')
          call apply_sbp42_2_boundary_method(ft_rs, in%subfields(n, m), domains%subdomains(n, m), direction, 's')
        end if
        !interpolation
        if (diff_method == 'sbp21_2') then
          if (domains%subdomains(n, m - 1)%nx == domains%subdomains(n, m)%nx) then
            call interp_identity(st_le, st_ile, 'coarse2fine')
            call interp_identity(st_rs, st_irs, 'fine2coarse')
            call interp_identity(ft_le, ft_ile, 'coarse2fine')
            call interp_identity(ft_rs, ft_irs, 'fine2coarse')
          else if (2 * domains%subdomains(n, m - 1)%nx == domains%subdomains(n, m)%nx) then
            call interp_MC2order_2to1ratio(st_le, st_ile, 'coarse2fine')
            call interp_MC2order_2to1ratio(st_rs, st_irs, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_le, ft_ile, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_rs, ft_irs, 'fine2coarse')
          else if (domains%subdomains(n, m - 1)%nx == 2 * domains%subdomains(n, m)%nx) then
            call interp_MC2order_2to1ratio(st_le, st_ile, 'fine2coarse')
            call interp_MC2order_2to1ratio(st_rs, st_irs, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_le, ft_ile, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_rs, ft_irs, 'coarse2fine')
          end if
        else if (diff_method == 'sbp42_2') then
          if (domains%subdomains(n, m - 1)%nx == domains%subdomains(n, m)%nx) then
            call interp_identity(st_le, st_ile, 'coarse2fine')
            call interp_identity(st_rs, st_irs, 'fine2coarse')
            call interp_identity(ft_le, ft_ile, 'coarse2fine')
            call interp_identity(ft_rs, ft_irs, 'fine2coarse')
          else if (2 * domains%subdomains(n, m - 1)%ny == domains%subdomains(n, m)%nx) then
            call interp_MC4order_2to1ratio(st_le, st_ile, 'coarse2fine')
            call interp_MC4order_2to1ratio(st_rs, st_irs, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_le, ft_ile, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_rs, ft_irs, 'fine2coarse')
          else if (domains%subdomains(n, m - 1)%nx == 2 * domains%subdomains(n, m)%nx) then
            call interp_MC4order_2to1ratio(st_le, st_ile, 'fine2coarse')
            call interp_MC4order_2to1ratio(st_rs, st_irs, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_le, ft_ile, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_rs, ft_irs, 'coarse2fine')
          end if
        else
          exit
        end if

        do k = ft_le%is, ft_le%ie
          st_le%f(k, 0) = st_le%f(k, 0) - st_irs%f(k, 0)
        end do
        do k = ft_rs%is, ft_rs%ie
          st_rs%f(k, 0) = st_rs%f(k, 0) - st_ile%f(k, 0)
        end do

        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method_second_terms(stl, st_le, domains%subdomains(n, m - 1), 'y', 'e')
          call apply_sbp21_2_boundary_method_second_terms(str, st_rs, domains%subdomains(n, m), 'y', 's')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method_second_terms(stl, st_le, domains%subdomains(n, m - 1), 'y', 'e')
          call apply_sbp42_2_boundary_method_second_terms(str, st_rs, domains%subdomains(n, m), 'y', 's')
        end if

!apply yi
        do k = ft_le%is, ft_le%ie
          ft_df = (coefs(n, m - 1) * ft_le%f(k, 0) + coefs(n, m) * ft_irs%f(k, 0) + (coefs(n, m)+coefs(n, m - 1)) * q / domains%subdomains(n, m - 1)%dy / 2.0_8 * st_le%f(k, 0)) / (domains%subdomains(n, m - 1)%dy * h0)
          tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je) = tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je) - ft_df / 2.0_8 + coefs(n, m - 1) * stl%f(k, 1) / 2.0_8
          tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je - 1) = tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je - 1) + coefs(n, m - 1) * stl%f(k, 2) / 2.0_8
          tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je - 2) = tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je - 2) + coefs(n, m - 1) * stl%f(k, 3) / 2.0_8
          tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je - 3) = tend%subfields(n, m - 1)%f(k, tend%subfields(n, m - 1)%je - 3) + coefs(n, m - 1) * stl%f(k, 4) / 2.0_8
        end do
        do k = ft_rs%is, ft_rs%ie
          ft_df = (coefs(n, m) * ft_rs%f(k, 0) + coefs(n, m - 1) * ft_ile%f(k, 0) + (coefs(n, m)+coefs(n, m - 1)) * q / domains%subdomains(n, m)%dy / 2.0_8 * st_rs%f(k, 0)) / (domains%subdomains(n, m)%dy * h0)
          tend%subfields(n, m)%f(k, tend%subfields(n, m)%js) = tend%subfields(n, m)%f(k, tend%subfields(n, m)%js) - ft_df / 2.0_8 + coefs(n, m) * str%f(k, 1) / 2.0_8
          tend%subfields(n, m)%f(k, tend%subfields(n, m)%js + 1) = tend%subfields(n, m)%f(k, tend%subfields(n, m)%js + 1) + coefs(n, m) * str%f(k, 2) / 2.0_8
          tend%subfields(n, m)%f(k, tend%subfields(n, m)%js + 2) = tend%subfields(n, m)%f(k, tend%subfields(n, m)%js + 2) + coefs(n, m) * str%f(k, 3) / 2.0_8
          tend%subfields(n, m)%f(k, tend%subfields(n, m)%js + 3) = tend%subfields(n, m)%f(k, tend%subfields(n, m)%js + 3) + coefs(n, m) * str%f(k, 4) / 2.0_8
        end do
    end do
        !bound
        call ft_ls%init(domains%subdomains(n, 1)%is, domains%subdomains(n, 1)%ie, 0, 0)
        call ft_re%init(domains%subdomains(n, domains%num_sub_y)%is,     domains%subdomains(n, domains%num_sub_y)%ie, 0, 0)
        call st_ls%init(domains%subdomains(n, 1)%is, domains%subdomains(n, 1)%ie, 0, 0)
        call st_re%init(domains%subdomains(n, domains%num_sub_y)%is,     domains%subdomains(n, domains%num_sub_y)%ie, 0, 0)
        call ft_ls%create_similar(ft_ire)
        call ft_re%create_similar(ft_ils)
        call st_ls%create_similar(st_ire)
        call st_re%create_similar(st_ils)
        call stlb%init(domains%subdomains(n, 1)%is, domains%subdomains(n, 1)%ie, 1, 4)
        call strb%init(domains%subdomains(n, domains%num_sub_y)%is, domains%subdomains(n, domains%num_sub_y)%ie, 1, 4)
        !feeling
        do k = ft_ls%is, ft_ls%ie
          ft_ls%f(k, 0) = in%subfields(n, 1)%f(k, in%subfields(n, 1)%js)
          st_ls%f(k, 0) = in%subfields(n, 1)%f(k, in%subfields(n, 1)%js)
        end do
        do k = ft_re%is, ft_re%ie
          ft_re%f(k, 0) = in%subfields(n, domains%num_sub_y)%f(k, in%subfields(n, domains%num_sub_y)%je)
          st_re%f(k, 0) = in%subfields(n, domains%num_sub_y)%f(k, in%subfields(n, domains%num_sub_y)%je)
        end do
        !boundary method
        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method(ft_ls, in%subfields(n, 1), domains%subdomains(n, 1), direction, 's')
          call apply_sbp21_2_boundary_method(ft_re, in%subfields(n, domains%num_sub_y), domains%subdomains(n, domains%num_sub_y), direction, 'e')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method(ft_ls, in%subfields(n, 1), domains%subdomains(n, 1), direction, 's')
          call apply_sbp42_2_boundary_method(ft_re, in%subfields(n, domains%num_sub_y), domains%subdomains(n, domains%num_sub_y), direction, 'e')
        end if
        !interpolation
        if (diff_method == 'sbp21_2') then
          if (domains%subdomains(n, 1)%nx == domains%subdomains(n, domains%num_sub_y)%nx) then
            call interp_identity(st_ls, st_ils, 'coarse2fine')
            call interp_identity(st_re, st_ire, 'fine2coarse')
            call interp_identity(ft_ls, ft_ils, 'coarse2fine')
            call interp_identity(ft_re, ft_ire, 'fine2coarse')
          else if (2 * domains%subdomains(n, 1)%nx == domains%subdomains(n, domains%num_sub_y)%nx) then
            call interp_MC2order_2to1ratio(st_ls, st_ils, 'coarse2fine')
            call interp_MC2order_2to1ratio(st_re, st_ire, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_ls, ft_ils, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_re, ft_ire, 'fine2coarse')
          else if (domains%subdomains(n, 1)%nx == 2 * domains%subdomains(n, domains%num_sub_y)%nx) then
            call interp_MC2order_2to1ratio(st_ls, st_ils, 'fine2coarse')
            call interp_MC2order_2to1ratio(st_re, st_ire, 'coarse2fine')
            call interp_MC2order_2to1ratio(ft_ls, ft_ils, 'fine2coarse')
            call interp_MC2order_2to1ratio(ft_re, ft_ire, 'coarse2fine')
          end if
        else if (diff_method == 'sbp42_2') then
          if (domains%subdomains(n, 1)%nx == domains%subdomains(n, domains%num_sub_y)%nx) then
            call interp_identity(st_ls, st_ils, 'coarse2fine')
            call interp_identity(st_re, st_ire, 'fine2coarse')
            call interp_identity(ft_ls, ft_ils, 'coarse2fine')
            call interp_identity(ft_re, ft_ire, 'fine2coarse')
          else if (2 * domains%subdomains(n, 1)%nx == domains%subdomains(n, domains%num_sub_y)%nx) then
            call interp_MC4order_2to1ratio(st_ls, st_ils, 'coarse2fine')
            call interp_MC4order_2to1ratio(st_re, st_ire, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_ls, ft_ils, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_re, ft_ire, 'fine2coarse')
          else if (domains%subdomains(n, 1)%nx == 2 * domains%subdomains(n, domains%num_sub_y)%nx) then
            call interp_MC4order_2to1ratio(st_ls, st_ils, 'fine2coarse')
            call interp_MC4order_2to1ratio(st_re, st_ire, 'coarse2fine')
            call interp_MC4order_2to1ratio(ft_ls, ft_ils, 'fine2coarse')
            call interp_MC4order_2to1ratio(ft_re, ft_ire, 'coarse2fine')
          end if
        else
          exit
        end if
        do k = ft_ls%is, ft_ls%ie
          st_ls%f(k, 0) = st_ls%f(k, 0) - st_ire%f(k, 0)
        end do
        do k = ft_re%is, ft_re%ie
          st_re%f(k, 0) = st_re%f(k, 0) - st_ils%f(k, 0)
        end do

        if (diff_method == 'sbp21_2') then
          call apply_sbp21_2_boundary_method_second_terms(stlb, st_ls, domains%subdomains(n, 1), 'y', 's')
          call apply_sbp21_2_boundary_method_second_terms(strb, st_re, domains%subdomains(n, domains%num_sub_y), 'y', 'e')
        else if (diff_method == 'sbp42_2') then
          call apply_sbp42_2_boundary_method_second_terms(stlb, st_ls, domains%subdomains(n, 1), 'y', 's')
          call apply_sbp42_2_boundary_method_second_terms(strb, st_re, domains%subdomains(n, domains%num_sub_y), 'y', 'e')
        end if
!apply yb
        do k = ft_ls%is, ft_ls%ie
          ft_df = (coefs(n, 1) * ft_ls%f(k, 0) + coefs(n, domains%num_sub_y) * ft_ire%f(k, 0) + (coefs(n, 1) + coefs(n, domains%num_sub_y)) * q / domains%subdomains(n, 1)%dy / 2.0_8 * st_ls%f(k, 0)) / (domains%subdomains(n, 1)%dy * h0)
          tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js) = tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js) - ft_df / 2.0_8 + coefs(n, 1) * stlb%f(k, 1) / 2.0_8
          tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js + 1) = tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js + 1) + coefs(n, 1) * stlb%f(k, 2) / 2.0_8
          tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js + 2) = tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js + 2) + coefs(n, 1) * stlb%f(k, 3) / 2.0_8
          tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js + 3) = tend%subfields(n, 1)%f(k, tend%subfields(n, 1)%js + 3) + coefs(n, 1) * stlb%f(k, 4) / 2.0_8
        end do
        do k = ft_re%is, ft_re%ie
          ft_df = (coefs(n, domains%num_sub_y) * ft_re%f(k, 0) + coefs(n, 1) * ft_ils%f(k, 0) + (coefs(n, 1)+coefs(n, domains%num_sub_y)) * q / domains%subdomains(n, domains%num_sub_y)%dy / 2.0_8 * st_re%f(k, 0)) / (domains%subdomains(n, domains%num_sub_y)%dy * h0)
          tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je) = tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je) - ft_df / 2.0_8 + coefs(n, domains%num_sub_y) * strb%f(k, 1) / 2.0_8
          tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je - 1) = tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je - 1) + coefs(n, domains%num_sub_y) * strb%f(k, 2) / 2.0_8
          tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je - 2) = tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je - 2) + coefs(n, domains%num_sub_y) * strb%f(k, 3) / 2.0_8
          tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je - 3) = tend%subfields(n, domains%num_sub_y)%f(k, tend%subfields(n, domains%num_sub_y)%je - 3) + coefs(n, domains%num_sub_y) * strb%f(k, 4) / 2.0_8
        end do
    end do
  end if

end subroutine sbp_SAT_penalty_two_block_diffusion

end module SAT_mod