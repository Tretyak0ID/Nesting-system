program multi_domain_test
  use domain_mod,       only : domain_t
  use multi_domain_mod, only : multi_domain_t
  implicit none

  type(multi_domain_t)         :: multi_domain
  type(domain_t)               :: global_domain
  integer(kind=4), allocatable :: deg(:, :)
  integer(kind=4)              :: n, m

  allocate(deg(1:2, 1:1))

  deg(1, 1) = 1
  deg(2, 1) = 2

  call global_domain%init(-100.0_8, 100.0_8, -10, 10, -100.0_8, 100.0_8, -10, 10)
  call multi_domain%init(global_domain, 2, 1, deg)

  print *, 'first x: ',  multi_domain%subdomains(1, 1)%x
  print *, 'second x: ', multi_domain%subdomains(2, 1)%x

  print *, 'first y: \n',  multi_domain%subdomains(1, 1)%y
  print *, 'second y: ', multi_domain%subdomains(2, 1)%y


end program multi_domain_test
