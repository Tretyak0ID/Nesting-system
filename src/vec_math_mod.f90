module vec_math_mod

use field_mod,  only : field_t
use domain_mod, only : domain_t
use multi_domain_mod, only : multi_domain_t
use multi_grid_field_mod, only : multi_grid_field_t

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

function calc_l1_norm_domain(field, domain) result(out)

  type(multi_grid_field_t), intent(in) :: field
  type(multi_domain_t),     intent(in) :: domain
  real(kind=8)                         :: out

  integer(kind=4) :: i, j, n, m

  out = 0.0_8;
  do m = 1, domain%num_sub_y
  do n = 1, domain%num_sub_x

  do i = domain%subdomains(n,m)%is, domain%subdomains(n,m)%ie
    do j = domain%subdomains(n,m)%js, domain%subdomains(n,m)%je
      out = out + abs(field%subfields(n,m)%f(i, j)) * domain%subdomains(n,m)%dx * domain%subdomains(n,m)%dy
    end do
  end do

  end do
  end do

end function calc_l1_norm_domain

function calc_c_norm_domain(field, domain) result(out)

  type(multi_grid_field_t),  intent(in) :: field
  type(multi_domain_t), intent(in) :: domain
  real(kind=8)               :: out

  integer(kind=8) :: n, m

  out = 0.0_8
  do m = 1, domain%num_sub_y
  do n = 1, domain%num_sub_x
    if (maxval(abs(field%subfields(n,m)%f)) > out) then
      out = maxval(abs(field%subfields(n,m)%f))
    end if
  end do
  end do

end function calc_c_norm_domain

function calc_sqrt_l2_norm_domain(field, domain) result(out)

  type(multi_grid_field_t),  intent(in) :: field
  type(multi_domain_t), intent(in)      :: domain
  real(kind=8)                          :: out

  integer(kind=4) :: i, j, n, m

  out = 0.0_8;
  do m = 1, domain%num_sub_y
  do n = 1, domain%num_sub_x
  do i = domain%subdomains(n,m)%is, domain%subdomains(n,m)%ie
    do j = domain%subdomains(n,m)%js, domain%subdomains(n,m)%je
      out = out + field%subfields(n,m)%f(i ,j) ** 2
    end do
  end do
  end do
  end do
   out = sqrt(out)

end function calc_sqrt_l2_norm_domain

end module vec_math_mod
