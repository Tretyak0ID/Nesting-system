program swe_state_test
use swe_state_mod, only: swe_state_t
use field_mod, only: field_t
implicit none

  type(field_t) u, v, h
  type(swe_state_t) state

  call u%init(0, 10, 0, 10)
  call v%init(0, 10, 0, 10)
  call h%init(0, 10, 0, 10)

  u%f(0, 0) = 1
  v%f(0, 0) = 2
  h%f(0, 0) = 3

  state%u = u
  state%v = v
  state%h = h

  print *, "u field: ", state%u%f
  print *, "v field: ", state%v%f
  print *, "h field: ", state%h%f

end program swe_state_test
