program laplacian_test
  use domain_mod,           only : domain_t
  use field_mod,            only : field_t
  use multi_domain_mod,     only : multi_domain_t
  use multi_grid_field_mod, only : multi_grid_field_t
  use laplacian_mod,        only : calc_laplacian
  use sbp_differential_operator_mod, only : sbp21_2_t
  use const_mod,                     only : pi
  use vec_math_mod,                  only : calc_c_norm_field
  implicit none

  type(multi_domain_t)         :: multi_domain
  type(multi_grid_field_t)     :: in_field, lap_field, lap
  type(domain_t)               :: global_domain
  type(sbp21_2_t)              :: sbp21_2
  integer(kind=4), allocatable :: deg(:, :)
  integer(kind=4)              :: n, m, i, j
  real(kind=8)                 :: error_left, error_right
  real(kind=8),    allocatable :: coefs(:, :)

  allocate(deg(1:2, 1:1))

  deg(1, 1) = 1
  deg(2, 1) = 1

  allocate(coefs(1:2, 1:1))
  coefs(1, 1) = 1.0_8
  coefs(2, 1) = 1.0_8

  sbp21_2%name = 'sbp21_2'

  call global_domain%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call multi_domain%init(global_domain, 2, 1, deg)
  call in_field%init(multi_domain)
  call lap_field%init(multi_domain)
  call lap%init(multi_domain)

  do n = 1, multi_domain%num_sub_x
    do m = 1, multi_domain%num_sub_y
      do j = 0, multi_domain%subdomains(n, m)%ny
        do i = 0, multi_domain%subdomains(n, m)%nx
          in_field%subfields(n, m)%f(i, j)  = exp(-2.0_8 * pi * multi_domain%subdomains(n, m)%y(j) - 2.0_8 * pi * multi_domain%subdomains(n, m)%x(i))
          lap_field%subfields(n, m)%f(i, j) = 8.0_8 * pi ** 2  * exp(-2.0_8 * pi * multi_domain%subdomains(n, m)%y(j) - 2.0_8 * pi * multi_domain%subdomains(n, m)%x(i))
        end do
      end do
    end do
  end do

  call calc_laplacian(lap, in_field, multi_domain, coefs, sbp21_2)

  call lap%update_s1v1(-1.0_8, lap_field, multi_domain)

  error_left  = calc_c_norm_field(lap%subfields(1, 1), multi_domain%subdomains(1, 1))
  error_right = calc_c_norm_field(lap%subfields(2, 1), multi_domain%subdomains(2, 1))

  print *, 'left grid field : ',  error_left
  print *, 'right grid field : ', error_right
  print *, lap%subfields(1, 1)%f(:, :)

end program laplacian_test
