program mesh_test
use mesh_mod, only: mesh_t

  type(mesh_t) :: mesh
  call mesh%init(0.0_8, 1.0_8, 0, 10, 0.0_8, 1.0_8, 0, 20)

  print *, "xs = ", mesh%xs, "xe = ", mesh%xe
  print *, "ys = ", mesh%ys, "ye = ", mesh%ye
  print *, "nx = ", mesh%nx, "ny = ", mesh%ny
  print *, "dx = ", mesh%dx, " dy = ", mesh%dy
  print *, "mesh_x = ", mesh%mesh_x
  print *, "mesh_y = ", mesh%mesh_y

end program mesh_test
