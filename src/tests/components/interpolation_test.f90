program interpolation_test
use field_mod,         only : field_t
use read_write_mod,    only : write_field
use domain_mod,        only : domain_t
use interpolation_mod, only : interp_MC2order_2to1ratio, interp_MC4order_2to1ratio
use const_mod,         only : pi
implicit none

type(field_t) :: f1, f2, f
type(domain_t) :: domain1, domain2

integer(kind=8) :: i, j

call domain1%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 1.0_8, 0, 0)
call domain2%init(0.0_8, 2.0_8 * pi, 0, 128, 0.0_8, 1.0_8, 0, 0)
call f1%init_on_domain(domain1)
call f2%init_on_domain(domain2)
call f%init_on_domain(domain2)

do i = 16, 20 
  f1%f(i, 0) = 1.0_8
end do

call interp_MC4order_2to1ratio(f1, f2, 'coarse2fine')

call write_field(f1, domain1, './data/test_interp3.dat')
call write_field(f2, domain2, './data/test_interp4.dat')

end program interpolation_test
