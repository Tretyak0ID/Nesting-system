program sbp_operators_test
use sbp_differential_operator_mod, only: sbp21_t, sbp42_t
use field_mod,         only: field_t
use mesh_mod,        only: mesh_t
use const_mod,         only: pi

  type(field_t)  :: in_field, out_field, cos_field
  type(mesh_t) :: mesh
  type(sbp42_t)  :: sbp42
  integer        :: i, j

  call mesh%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call in_field%init(0, 64, 0, 64)
  call out_field%init(0, 64, 0, 64)
  call cos_field%init(0, 64, 0, 64)

  do j = 0, mesh%ny
    do i = 0, mesh%nx
      in_field%f(i, j) = sin(mesh%mesh_y(j))
      cos_field%f(i, j) = cos(mesh%mesh_y(j))
    end do
  end do

  call sbp42%apply(out_field, in_field, mesh, 'y')

  print *, cos_field%f(10, :) - out_field%f(10, :)

  !-------------------------------------------------

  call mesh%init(0.0_8, 2.0_8 * pi, 0, 128, 0.0_8, 2.0_8 * pi, 0, 128)
  call in_field%init(0, 128, 0, 128)
  call out_field%init(0, 128, 0, 128)
  call cos_field%init(0, 128, 0, 128)

  do j = 0, mesh%ny
    do i = 0, mesh%nx
      in_field%f(i, j) = sin(mesh%mesh_y(j))
      cos_field%f(i, j) = cos(mesh%mesh_y(j))
    end do
  end do

  call sbp42%apply(out_field, in_field, mesh, 'y')

  print *, cos_field%f(10, :) - out_field%f(10, :)

end program sbp_operators_test
