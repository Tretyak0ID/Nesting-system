program multi_grid_field
  use domain_mod,           only : domain_t
  use field_mod,            only : field_t
  use multi_domain_mod,     only : multi_domain_t
  use multi_grid_field_mod, only : multi_grid_field_t
  implicit none

  type(multi_domain_t)         :: multi_domain
  type(multi_grid_field_t)     :: milti_grid_field
  type(domain_t)               :: global_domain
  integer(kind=4), allocatable :: deg(:, :)
  integer(kind=4)              :: n, m

  allocate(deg(1:2, 1:1))

  deg(1, 1) = 1
  deg(2, 1) = 2

  call global_domain%init(-100.0_8, 100.0_8, -10, 10, -100.0_8, 100.0_8, -10, 10)
  call multi_domain%init(global_domain, 2, 1, deg)
  call milti_grid_field%init(multi_domain)

  print *, 'left grid field : ',  milti_grid_field%subfields(1, 1)%f
  print *, 'right grid field : ', milti_grid_field%subfields(2, 1)%f

end program multi_grid_field
