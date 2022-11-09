program stvec_swe_test
use stvec_swe_mod,        only: stvec_swe_t
use field_mod,            only: field_t
use multi_grid_field_mod, only: multi_grid_field_t
use domain_mod,           only: domain_t
use multi_domain_mod,     only: multi_domain_t
implicit none

  type(multi_grid_field_t) :: u, v, h
  type(stvec_swe_t)        :: state1, state2
  type(domain_t)           :: global_domain
  type(multi_domain_t)     :: multi_domain
  integer(kind=4), allocatable :: deg(:, :)

  allocate(deg(1:2, 1:1))

  deg(1, 1) = 1
  deg(2, 1) = 2

  call global_domain%init(0.0_8, 1.0_8, 0, 10, 0.0_8, 1.0_8, 0, 10)
  call multi_domain%init(global_domain, 2, 1, deg)

  call u%init(multi_domain)
  call v%init(multi_domain)
  call h%init(multi_domain)

  u%subfields(1, 1)%f(0, 0) = 1
  v%subfields(1, 1)%f(0, 0) = 1
  h%subfields(1, 1)%f(0, 0) = 1

  state1%u = u
  state1%v = v
  state1%h = h

  call state1%create_similar(state2)

  print *, "u2 field: ", state2%u%subfields(2, 1)%f
  print *, "v2 field: ", state2%v%subfields(2, 1)%f
  print *, "h2 field: ", state2%h%subfields(2, 1)%f

end program stvec_swe_test
