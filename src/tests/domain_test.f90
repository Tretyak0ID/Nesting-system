program domain_test
use domain_mod, only: domain_t

  type(domain_t) :: domain
  call domain%init(0.0_8, 1.0_8, 0, 10, 0.0_8, 1.0_8, 0, 20)

  print *, "xs = ", domain%xs, "xe = ", domain%xe
  print *, "ys = ", domain%ys, "ye = ", domain%ye
  print *, "nx = ", domain%nx, "ny = ", domain%ny
  print *, "dx = ", domain%dx, " dy = ", domain%dy
  print *, "x = ", domain%x
  print *, "y = ", domain%y

end program domain_test
