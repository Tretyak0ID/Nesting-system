program sbp21_2_test
use sbp_differential_operator_mod, only: sbp21_t, sbp42_t, sbp21_2_t
use field_mod,         only: field_t
use domain_mod,        only: domain_t
use const_mod,         only: pi

  type(field_t)  :: in_field, out_field, cos_field
  type(domain_t) :: domain
  type(sbp21_2_t)  :: sbp21_2
  integer        :: i, j

  !-------------------------------------------------

  call domain%init(0.0_8, 2.0_8 * pi, 0, 128, 0.0_8, 2.0_8 * pi, 0, 128)
  call in_field%init(0, 128, 0, 128)
  call out_field%init(0, 128, 0, 128)
  call cos_field%init(0, 128, 0, 128)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j)  =  sin(domain%y(j))
      cos_field%f(i, j) =  -sin(domain%y(j))
    end do
  end do

  call sbp21_2%apply(out_field, in_field, domain, 'y')

  print *, cos_field%f(:, 0) - out_field%f(:, 0)
  print *, cos_field%f(:, 10) - out_field%f(:, 10)
  print *, cos_field%f(:, 128) - out_field%f(:, 128)

end program sbp21_2_test
