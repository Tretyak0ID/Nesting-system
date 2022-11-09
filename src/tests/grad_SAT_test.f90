program grad_SAT_test
  use domain_mod,           only : domain_t
  use field_mod,            only : field_t
  use multi_domain_mod,     only : multi_domain_t
  use multi_grid_field_mod, only : multi_grid_field_t
  use grad_mod,             only : calc_grad
  use sbp_differential_operator_mod, only : sbp21_t
  use const_mod,                     only : pi
  use vec_math_mod,                  only : calc_c_norm_field
  implicit none

  type(multi_domain_t)         :: multi_domain
  type(multi_grid_field_t)     :: in_field, gx_field, gy_field, gx, gy
  type(domain_t)               :: global_domain
  type(sbp21_t)                :: sbp21
  integer(kind=4), allocatable :: deg(:, :)
  integer(kind=4)              :: n, m, i, j
  real(kind=8)                 :: error_x_left, error_x_right, error_y_left, error_y_right

  allocate(deg(1:2, 1:1))

  deg(1, 1) = 1
  deg(2, 1) = 1

  call global_domain%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call multi_domain%init(global_domain, 2, 1, deg)
  call in_field%init(multi_domain)
  call gx_field%init(multi_domain)
  call gy_field%init(multi_domain)
  call gx%init(multi_domain)
  call gy%init(multi_domain)

  do n = 1, multi_domain%num_sub_x
    do m = 1, multi_domain%num_sub_y
      do j = 0, multi_domain%subdomains(n, m)%ny
        do i = 0, multi_domain%subdomains(n, m)%nx
          in_field%subfields(n, m)%f(i, j) = sin(multi_domain%subdomains(n, m)%y(j)) * sin(multi_domain%subdomains(n, m)%x(i))
          gx_field%subfields(n, m)%f(i, j) = cos(multi_domain%subdomains(n, m)%x(i)) * sin(multi_domain%subdomains(n, m)%y(j))
          gy_field%subfields(n, m)%f(i, j) = cos(multi_domain%subdomains(n, m)%y(j)) * sin(multi_domain%subdomains(n, m)%x(i))
        end do
      end do
    end do
  end do

  call calc_grad(gx, gy, in_field, multi_domain, sbp21, sbp21)

  call gx%update_s1v1(-1.0_8, gx_field, multi_domain)
  call gy%update_s1v1(-1.0_8, gy_field, multi_domain)

  error_x_left  = calc_c_norm_field(gx%subfields(1, 1), multi_domain%subdomains(1, 1))
  error_x_right = calc_c_norm_field(gx%subfields(2, 1), multi_domain%subdomains(2, 1))
  error_y_left  = calc_c_norm_field(gy%subfields(1, 1), multi_domain%subdomains(1, 1))
  error_y_right = calc_c_norm_field(gy%subfields(2, 1), multi_domain%subdomains(2, 1))

  print *, 'x left grid field : ',  error_x_left
  print *, 'x right grid field : ', error_x_right
  print *, 'y left grid field : ',  error_y_left
  print *, 'y right grid field : ', error_y_right
end program grad_SAT_test
