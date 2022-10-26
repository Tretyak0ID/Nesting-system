program field_test
use field_mod, only: field_t
use domain_mod, only : domain_t

  type(field_t) :: f1, f2
  type(domain_t) :: domain

  call domain%init(0.0_8, 1.0_8, -2, 2, 0.0_8, 1.0_8, -2, 2)
  call f1%init_on_domain(domain)

  f1%f(1, -1) = 10.0_8
  call f2%assign(2.0_8, f1, domain)

  print *, f2%f

end program field_test
