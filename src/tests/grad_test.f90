program grad_test
use grad_mod,                          only: calc_grad
use sbp_differential_operator_mod,     only: sbp21_t, sbp42_t
use central_differential_operator_mod, only: central2_t, central4_t
use field_mod,                         only: field_t
use domain_mod,                        only: domain_t
use const_mod,                         only: pi
implicit none

  type(field_t)    :: in_field, gx, gy, gx_field, gy_field
  type(domain_t)   :: domain
  type(sbp21_t)    :: sbp21
  type(sbp42_t)    :: sbp42
  type(central2_t) :: central2
  type(central4_t) :: central4
  integer          :: i, j

  call domain%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call in_field%init(0, 64, 0, 64)
  call gx%init(0, 64, 0, 64)
  call gy%init(0, 64, 0, 64)
  call gx_field%init(0, 64, 0, 64)
  call gy_field%init(0, 64, 0, 64)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%domain_y(j)) * sin(domain%domain_x(i))
      gx_field%f(i, j) = cos(domain%domain_x(i)) * sin(domain%domain_y(j))
      gy_field%f(i, j) = cos(domain%domain_y(j)) * sin(domain%domain_x(i))
    end do
  end do

  call calc_grad(gx, gy, in_field, domain, sbp42, sbp42)

  print *, 'grad gx:', abs(gx_field%f(10, :) - gx%f(10, :))
  print *, 'grad gy:', abs(gy_field%f(:, 10) - gy%f(:, 10))

end program grad_test
