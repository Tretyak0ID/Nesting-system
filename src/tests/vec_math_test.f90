program vec_math_test
use div_mod,                           only: calc_div
use sbp_differential_operator_mod,     only: sbp21_t, sbp42_t
use central_differential_operator_mod, only: central2_t, central4_t
use field_mod,                         only: field_t
use domain_mod,                        only: domain_t
use const_mod,                         only: pi
use vec_math_mod,                      only: calc_mass_field, calc_c_norm_field, calc_sqrt_l2_norm_field
implicit none

  type(field_t)    :: in_field, div, div_field
  type(domain_t)   :: domain
  type(sbp21_t)    :: sbp21
  type(sbp42_t)    :: sbp42
  type(central2_t) :: central2
  type(central4_t) :: central4

  real(kind=8)     :: c, l2, mass
  type(field_t)    :: diff_m
  integer          :: i, j

  call domain%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call in_field%init(0, 64, 0, 64)
  call div%init(0, 64, 0, 64)
  call div_field%init(0, 64, 0, 64)
  call diff_m%init_on_domain(domain)

  do j = domain%is, domain%ie
    do i = domain%js, domain%je
      in_field%f(i, j)  = sin(domain%y(j)) * sin(domain%x(i))
      div_field%f(i, j) = cos(domain%x(i)) * sin(domain%y(j)) + cos(domain%y(j)) * sin(domain%x(i))
    end do
  end do

  call calc_div(div, in_field, in_field, domain, sbp42, sbp42)

  do j = domain%is, domain%ie
    do i = domain%js, domain%je
      diff_m%f(i, j) = div_field%f(i, j) - div%f(i, j)
    end do
  end do

  mass = calc_mass_field(diff_m, domain)
  c = calc_c_norm_field(diff_m, domain)
  l2 = calc_sqrt_l2_norm_field(diff_m, domain)

  print *, 'divergention error mass-norm:', mass
  print *, 'divergention error c-norm:', c
  print *, 'divergention error l2-norm:', l2

end program vec_math_test
