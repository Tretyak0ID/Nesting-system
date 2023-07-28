program central_operators_test
use central_differential_operator_mod, only: central2_t, central4_t
use field_mod,         only: field_t
use domain_mod,        only: domain_t
use const_mod,         only: pi

  type(field_t)  :: in_field, out_field, cos_field
  type(domain_t) :: domain
  type(central2_t)  :: central2
  integer        :: i, j

  call domain%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call in_field%init(0, 64, 0, 64)
  call out_field%init(0, 64, 0, 64)
  call cos_field%init(0, 64, 0, 64)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%y(j))
      cos_field%f(i, j) = cos(domain%y(j))
    end do
  end do

  call central2%apply(out_field, in_field, domain, 'y')

  print *, cos_field%f(0, :) - out_field%f(0, :)

end program central_operators_test
