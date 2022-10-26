module initial_conditions_mod
  use field_mod,     only : field_t
  use domain_mod,    only : domain_t
  use stvec_mod,     only : stvec_t
  use stvec_swe_mod, only : stvec_swe_t
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

end module initial_conditions_mod
