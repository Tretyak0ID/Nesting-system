module vec_math_mod

use field_mod,  only : field_t
use domain_mod, only : domain_t

implicit none

contains

function calc_mass_field(field, domain) result(out)

  type(field_t),  intent(in) :: field
  type(domain_t), intent(in) :: domain
  real(kind=8)               :: out

  integer(kind=4) :: i, j

  out = 0.0_8;
  do i = domain%is, domain%ie
    do j = domain%js, domain%je
      out = out + abs(field%f(i, j)) * domain%dx * domain%dy
    end do
  end do

end function calc_mass_field

function calc_c_norm_field(field, domain) result(out)

  type(field_t),  intent(in) :: field
  type(domain_t), intent(in) :: domain
  real(kind=8)               :: out

  out = maxval(abs(field%f))

end function calc_c_norm_field

function calc_sqrt_l2_norm_field(field, domain) result(out)

  type(field_t),  intent(in) :: field
  type(domain_t), intent(in) :: domain
  real(kind=8)               :: out

  integer(kind=4) :: i, j

  out = 0.0_8;
  do i = domain%is, domain%ie
    do j = domain%js, domain%je
      out = out + field%f(i ,j) ** 2 * domain%dx * domain%dy
    end do
  end do

end function calc_sqrt_l2_norm_field

end module vec_math_mod
