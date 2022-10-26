program stvec_swe_test
use stvec_swe_mod, only: stvec_swe_t
use field_mod,     only: field_t
use domain_mod,    only: domain_t
implicit none

  type(field_t) :: u, v, h
  type(stvec_swe_t) :: state1, state2
  type(domain_t) :: domain

  call domain%init(0.0_8, 1.0_8, 0, 10, 0.0_8, 1.0_8, 0, 10)

  call u%init_on_domain(domain)
  call v%init_on_domain(domain)
  call h%init_on_domain(domain)

  u%f(0, 0) = 1
  v%f(0, 0) = 2
  h%f(0, 0) = 3

  state1%u = u
  state1%v = v
  state1%h = h

  call state2%assign(1.0_8, state1, domain)

  print *, "u1 field: ", state1%u%f
  print *, "v1 field: ", state1%v%f
  print *, "h1 field: ", state1%h%f

  print *, "u2 field: ", state2%u%f
  print *, "v2 field: ", state2%v%f
  print *, "h2 field: ", state2%h%f

end program stvec_swe_test
