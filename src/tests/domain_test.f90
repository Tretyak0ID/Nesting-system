program domain_test
use domain_mod, only: domain_t

  type(domain_t) :: domain
  call domain%init(0.0_8, 1.0_8, 10, 0.0_8, 1.0_8, 20)

  print *, "xs = ", domain%xs, "xe = ", domain%xe
  print *, "ys = ", domain%ys, "ye = ", domain%ye
  print *, "nx = ", domain%nx, "ny = ", domain%ny
  print *, "dx = ", domain%dx, " dy = ", domain%dy
  print *, "mesh_x = ", domain%mesh_x
  print *, "mesh_y = ", domain%mesh_y

end program domain_test
