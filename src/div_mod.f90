module div_mod
use field_mod,  only: field_t
use mesh_mod, only: mesh_t
use differential_operator_mod, only: differential_operator_t
implicit none

contains

  subroutine calc_div(div, inx, iny, mesh, diff_opx, diff_opy)
    type  (field_t),                 intent(inout) :: div
    type  (field_t),                 intent(in)    :: inx, iny
    type  (mesh_t),                  intent(in)    :: mesh
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy
    type(field_t)   :: gx_buff
    type(field_t)   :: gy_buff
    integer(kind=8) :: i, j

    call gx_buff%init_on_mesh(mesh)
    call gy_buff%init_on_mesh(mesh)

    call diff_opx%apply(gx_buff, inx, mesh, 'x')
    call diff_opy%apply(gy_buff, iny, mesh, 'y')

    call div%assign(1.0_8, gx_buff, 1.0_8, gy_buff, mesh)

  end subroutine

end module div_mod
