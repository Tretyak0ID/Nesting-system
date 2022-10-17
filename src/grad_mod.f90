module grad_mod
use field_mod,  only: field_t
use mesh_mod, only: mesh_t
use differential_operator_mod, only: differential_operator_t
implicit none

contains

  subroutine calc_grad(gx, gy, in, mesh, diff_opx, diff_opy)
    type  (field_t),                 intent(inout) :: gx, gy
    type  (field_t),                 intent(in)    :: in
    type  (mesh_t),                intent(in)    :: mesh
    class (differential_operator_t), intent(in)    :: diff_opx, diff_opy

    call diff_opx%apply(gx, in, mesh, 'x')
    call diff_opy%apply(gy, in, mesh, 'y')
  end subroutine

end module grad_mod
