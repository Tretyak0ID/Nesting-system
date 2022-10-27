module initial_conditions_mod
  use field_mod,            only : field_t
  use domain_mod,           only : domain_t
  use multi_domain_mod,     only : multi_domain_t
  use multi_grid_field_mod, only : multi_grid_field_t
  use stvec_mod,            only : stvec_t
  use stvec_swe_mod,        only : stvec_swe_t
  use const_mod,            only : Earth_grav, pcori, pi
implicit none

contains

  subroutine swm_gaussian_hill(vec, multi_domain, h_mean, kx, ky)

    type(stvec_swe_t),    intent(inout) :: vec
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean, kx, ky

    real(kind=8)    :: dx, dy
    integer(kind=8) :: i, j, n, m

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
            dx = (multi_domain%subdomains(n, m)%x(i) - 0.8_8 * (maxval(multi_domain%global_domain%x) - minval(multi_domain%global_domain%x)) / 2.0_8) / multi_domain%global_domain%x(multi_domain%global_domain%ie)
            dy = (multi_domain%subdomains(n, m)%y(j) - 0.8_8 * (maxval(multi_domain%global_domain%y) - minval(multi_domain%global_domain%y)) / 2.0_8) / multi_domain%global_domain%y(multi_domain%global_domain%je)
            vec%h%subfields(n, m)%f(i, j) = h_mean + 0.1_8 * h_mean * exp(- kx * (dx ** 2.0_8) - ky * (dy ** 2.0_8))
          end do
        end do
      end do
    end do

  end subroutine swm_gaussian_hill



  subroutine swm_geostrophic_balance(vec, multi_domain, h_mean, scale_h, scale_sigma)

    type(stvec_swe_t),    intent(inout) :: vec
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean, scale_h, scale_sigma

    type(multi_grid_field_t) :: dhdr, vtan
    real(kind=8)             :: dx, dy
    integer(kind=8)          :: i, j, n, m

    call dhdr%init(multi_domain)
    call vtan%init(multi_domain)

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
            dx = (multi_domain%subdomains(n, m)%x(i) - (maxval(multi_domain%global_domain%x) - minval(multi_domain%global_domain%x)) / 2.0_8)
            dy = (multi_domain%subdomains(n, m)%y(j) - (maxval(multi_domain%global_domain%y) - minval(multi_domain%global_domain%y)) / 2.0_8)

            dhdr%subfields(n, m)%f(i, j) = h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8) * scale_h * 2.0_8 * sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma ** 2.0_8
            vtan%subfields(n, m)%f(i, j) = (- sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori + sqrt((sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori) ** 2.0_8 + 4.0_8 * Earth_grav * sqrt(dx ** 2.0_8 + dy ** 2.0_8) * dhdr%subfields(n, m)%f(i, j))) / 2.0_8

            vec%h%subfields(n, m)%f(i, j) = h_mean - (h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8)) * scale_h
            vec%u%subfields(n, m)%f(i, j) = - vtan%subfields(n, m)%f(i, j) * cos(pi / 2.0_8 - atan2(dy, dx))
            vec%v%subfields(n, m)%f(i, j) =   vtan%subfields(n, m)%f(i, j) * sin(pi / 2.0_8 - atan2(dy, dx))
          end do
        end do
      end do
    end do

  end subroutine swm_geostrophic_balance



  subroutine swm_geostrophic_cyclone(vec, multi_domain, h_mean, scale_h, scale_sigma)

    type(stvec_swe_t),    intent(inout) :: vec
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean, scale_h, scale_sigma

    type(multi_grid_field_t) :: dhdr, vtan, h0, u0, v0
    real(kind=8)             :: dx, dy
    integer(kind=8)          :: i, j, n, m

    call dhdr%init(multi_domain)
    call vtan%init(multi_domain)
    call h0%init(multi_domain)
    call u0%init(multi_domain)
    call v0%init(multi_domain)

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
            dx = (multi_domain%subdomains(n, m)%x(i) - (maxval(multi_domain%global_domain%x) - minval(multi_domain%global_domain%x)) / 2.0_8)
            dy = (multi_domain%subdomains(n, m)%y(j) - (maxval(multi_domain%global_domain%y) - minval(multi_domain%global_domain%y)) / 2.0_8)

            dhdr%subfields(n, m)%f(i, j) = h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8) * scale_h * 2.0_8 * sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma ** 2.0_8
            vtan%subfields(n, m)%f(i, j) = (- sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori + sqrt((sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori) ** 2.0_8 + 4.0_8 * Earth_grav * sqrt(dx ** 2.0_8 + dy ** 2.0_8) * dhdr%subfields(n, m)%f(i, j))) / 2.0_8
            h0%subfields(n, m)%f(i, j) = pcori  / Earth_grav * cos((multi_domain%subdomains(n, m)%y(j) - 10.0_8 ** 7.0_8) / (6371.22_8 * 1000.0_8)) * 10.0_8 * 6371.22_8 * 1000.0_8
            u0%subfields(n, m)%f(i, j) = sin((multi_domain%subdomains(n, m)%y(j) - 10.0_8 ** 7.0_8) / (6371.22_8 * 1000.0_8)) * 10.0_8

            vec%h%subfields(n, m)%f(i, j) = h_mean - (h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8)) * scale_h + h0%subfields(n, m)%f(i, j)
            vec%u%subfields(n, m)%f(i, j) = - vtan%subfields(n, m)%f(i, j) * cos(pi / 2.0_8 - atan2(dy, dx)) + u0%subfields(n, m)%f(i, j)
            vec%v%subfields(n, m)%f(i, j) =   vtan%subfields(n, m)%f(i, j) * sin(pi / 2.0_8 - atan2(dy, dx))
          end do
        end do
      end do
    end do

  end subroutine swm_geostrophic_cyclone

end module initial_conditions_mod
