module initial_conditions_mod
  use field_mod,     only : field_t
  use domain_mod,    only : domain_t
  use stvec_mod,     only : stvec_t
  use stvec_swe_mod, only : stvec_swe_t
  use const_mod,     only : Earth_grav, pcori, pi
implicit none

contains

  subroutine swm_gaussian_hill(vec, domain, h_mean, kx, ky)

    type(stvec_swe_t), intent(inout) :: vec
    type(domain_t),    intent(in)    :: domain
    real(kind=8),      intent(in)    :: h_mean, kx, ky

    real(kind=8)    :: dx, dy
    integer(kind=8) :: i, j

    do i = domain%is, domain%ie
      do j = domain%js, domain%je
        dx = (domain%x(i) - 0.8_8 * (maxval(domain%x) - minval(domain%x)) / 2.0_8) / domain%x(domain%ie)
        dy = (domain%y(j) - 0.8_8 * (maxval(domain%y) - minval(domain%y)) / 2.0_8) / domain%y(domain%je)
        vec%h%f(i, j) = h_mean + 0.1_8 * h_mean * exp(- kx * (dx ** 2.0_8) - ky * (dy ** 2.0_8))
      end do
    end do

  end subroutine swm_gaussian_hill



  subroutine swm_geostrophic_balance(vec, domain, h_mean, scale_h, scale_sigma)

    type(stvec_swe_t), intent(inout) :: vec
    type(domain_t),    intent(in)    :: domain
    real(kind=8),      intent(in)    :: h_mean, scale_h, scale_sigma

    type(field_t)   :: dhdr, vtan
    real(kind=8)    :: dx, dy
    integer(kind=8) :: i, j

    call dhdr%init_on_domain(domain)
    call vtan%init_on_domain(domain)

    do i = domain%is, domain%ie
      do j = domain%js, domain%je
        dx = (domain%x(i) - (maxval(domain%x) - minval(domain%x)) / 2.0_8)
        dy = (domain%y(j) - (maxval(domain%y) - minval(domain%y)) / 2.0_8)

        dhdr%f(i, j) = h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8) * scale_h * 2.0_8 * sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma ** 2.0_8
        vtan%f(i, j) = (- sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori + sqrt((sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori) ** 2.0_8 + 4.0_8 * Earth_grav * sqrt(dx ** 2.0_8 + dy ** 2.0_8) * dhdr%f(i, j))) / 2.0_8

        vec%h%f(i, j) = h_mean - (h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8)) * scale_h
        vec%u%f(i, j) = - vtan%f(i, j) * cos(pi / 2.0_8 - atan2(dy, dx))
        vec%v%f(i, j) =   vtan%f(i, j) * sin(pi / 2.0_8 - atan2(dy, dx))
      end do
    end do

  end subroutine swm_geostrophic_balance



  subroutine swm_geostrophic_cyclone(vec, domain, h_mean, scale_h, scale_sigma)

    type(stvec_swe_t), intent(inout) :: vec
    type(domain_t),    intent(in)    :: domain
    real(kind=8),      intent(in)    :: h_mean, scale_h, scale_sigma

    type(field_t)   :: dhdr, vtan, h0, u0, v0
    real(kind=8)    :: dx, dy
    integer(kind=8) :: i, j

    call dhdr%init_on_domain(domain)
    call vtan%init_on_domain(domain)
    call h0%init_on_domain(domain)
    call u0%init_on_domain(domain)
    call v0%init_on_domain(domain)

    do i = domain%is, domain%ie
      do j = domain%js, domain%je
        dx = (domain%x(i) - (maxval(domain%x) - minval(domain%x)) / 2.0_8)
        dy = (domain%y(j) - (maxval(domain%y) - minval(domain%y)) / 2.0_8)

        dhdr%f(i, j) = h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8) * scale_h * 2.0_8 * sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma ** 2.0_8
        vtan%f(i, j) = (- sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori + sqrt((sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori) ** 2.0_8 + 4.0_8 * Earth_grav * sqrt(dx ** 2.0_8 + dy ** 2.0_8) * dhdr%f(i, j))) / 2.0_8

        h0%f(i, j) = pcori  / Earth_grav * cos((domain%y(j) - 10.0_8 ** 7.0_8) / (6371.22_8 * 1000.0_8)) * 10.0_8 * 6371.22_8 * 1000.0_8
        u0%f(i, j) = sin((domain%y(j) - 10.0_8 ** 7.0_8) / (6371.22_8 * 1000.0_8)) * 10.0_8

        vec%h%f(i, j) = h_mean - (h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8)) * scale_h + h0%f(i, j)
        vec%u%f(i, j) = - vtan%f(i, j) * cos(pi / 2.0_8 - atan2(dy, dx)) + u0%f(i, j)
        vec%v%f(i, j) =   vtan%f(i, j) * sin(pi / 2.0_8 - atan2(dy, dx)) + v0%f(i, j)
      end do
    end do

  end subroutine swm_geostrophic_cyclone

end module initial_conditions_mod
