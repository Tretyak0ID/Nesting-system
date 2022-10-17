module div_mod
use field_mod,  only: field_t
use mesh_mod, only: mesh_t
use differential_operator_mod, only: differential_operator_t
implicit none

contains

  subroutine calc_div(div, inx, iny, mesh, diff_opx, diff_opy)
    type  (field_t),                 intent(inout) :: div
    type  (field_t),                 intent(in)    :: inx, iny
    type  (mesh_t),                intent(in)    :: mesh
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy
    type(field_t)   :: gx_buff
    type(field_t)   :: gy_buff
    integer(kind=8) :: i, j

    call gx_buff%init_on_mesh(mesh)
    call gy_buff%init_on_mesh(mesh)

    call diff_opx%apply(gx_buff, inx, mesh, 'x')
    call diff_opy%apply(gy_buff, iny, mesh, 'y')

    do i = mesh%sindx, mesh%eindx
      do j = mesh%sindy, mesh%eindy
        div%f(i, j) = gx_buff%f(i, j) + gy_buff%f(i, j)
      end do
    end do
  end subroutine

end module div_mod
