program interpolation_test
use field_mod,         only : field_t
use read_write_mod,    only : write_field
use domain_mod,        only : domain_t
use interpolation_mod, only : interp_1d_sbp21_2to1_ratio
use const_mod,         only : pi
implicit none

type(field_t) :: f1, f2, f
type(domain_t) :: domain1, domain2

integer(kind=8) :: i, j

call domain1%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 1.0_8, 0, 0)
call domain2%init(0.0_8, 2.0_8 * pi, 0, 32, 0.0_8, 1.0_8, 0, 0)
call f1%init_on_domain(domain1)
call f2%init_on_domain(domain2)
call f%init_on_domain(domain1)

do i = f2%is, f2%ie
  f2%f(i, 0) = sin(domain2%x(i))
end do

call interp_1d_sbp21_2to1_ratio(f2, f1, 'coarse2fine')

do i = f%is, f%ie
  f%f(i, 0) = sin(domain1%x(i)) - f1%f(i, 0)
end do

call write_field(f, domain1, './data/test_interp.dat')

print *, f%f

end program interpolation_test
