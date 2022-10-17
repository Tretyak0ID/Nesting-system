program div_test
use div_mod,                           only: calc_div
use sbp_differential_operator_mod,     only: sbp21_t, sbp42_t
use central_differential_operator_mod, only: central2_t, central4_t
use field_mod,                         only: field_t
use mesh_mod,                        only: mesh_t
use const_mod,                         only: pi
implicit none

  type(field_t)    :: in_field, div, div_field
  type(mesh_t)   :: mesh
  type(sbp21_t)    :: sbp21
  type(sbp42_t)    :: sbp42
  type(central2_t) :: central2
  type(central4_t) :: central4
  integer          :: i, j

  call mesh%init(0.0_8, 2.0_8 * pi, 0, 64, 0.0_8, 2.0_8 * pi, 0, 64)
  call in_field%init(0, 64, 0, 64)
  call div%init(0, 64, 0, 64)
  call div_field%init(0, 64, 0, 64)

  do j = 0, mesh%ny
    do i = 0, mesh%nx
      in_field%f(i, j)  = sin(mesh%mesh_y(j)) * sin(mesh%mesh_x(i))
      div_field%f(i, j) = cos(mesh%mesh_x(i)) * sin(mesh%mesh_y(j)) + cos(mesh%mesh_y(j)) * sin(mesh%mesh_x(i))
    end do
  end do

  call calc_div(div, in_field, in_field, mesh, sbp42, sbp42)

  print *, 'divergention error:', maxval(abs(div_field%f - div%f))

end program div_test
