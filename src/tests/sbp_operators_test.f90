program sbp_operators_test
use sbp_differential_operator_mod, only: sbp21_t, sbp42_t
use field_mod,         only: field_t
use domain_mod,        only: domain_t
use const_mod,         only: pi
use read_write_mod,    only: write_field

  type(field_t)  :: in_field, out_field, cos_field
  type(domain_t) :: domain
  type(sbp42_t)  :: sbp42
  integer        :: i, j

  call domain%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call in_field%init(0, 64, 0, 64)
  call out_field%init(0, 64, 0, 64)
  call cos_field%init(0, 64, 0, 64)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = domain%x(i) ** 2.0_8
      cos_field%f(i, j) = 2.0_8 * domain%x(i)
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'y')

  call write_field(out_field, domain, './data/sbp_test_y.dat', 1)
  call write_field(out_field, domain, './data/cos_field_y.dat', 1)

  print *, out_field%f(:, 1) - cos_field%f(:, 1)

  !-------------------------------------------------

  call domain%init(0.0_8, 2.0_8 * pi, 0, 128, 0.0_8, 2.0_8 * pi, 0, 128)
  call in_field%init(0, 128, 0, 128)
  call out_field%init(0, 128, 0, 128)
  call cos_field%init(0, 128, 0, 128)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, out_field%f(:, 1) - cos_field%f(:, 1)

  call write_field(out_field, domain, './data/sbp_test_x.dat', 1)
  call write_field(out_field, domain, './data/cos_field_x.dat', 1)

end program sbp_operators_test
