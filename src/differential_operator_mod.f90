module differential_operator_mod
use field_mod,  only: field_t
use mesh_mod, only: mesh_t
implicit none

  type, abstract :: differential_operator_t

    contains

      procedure(apply_i), deferred :: apply

  end type differential_operator_t

abstract interface
  subroutine apply_i(this, out, in, mesh, direction)

    import differential_operator_t, mesh_t, field_t

    class    (differential_operator_t), intent(in)  :: this
    type     (field_t),                 intent(inout) :: out
    type     (field_t),                 intent(in)  :: in
    type     (mesh_t),                intent(in)  :: mesh
    character(len=1),                   intent(in)  :: direction

  end subroutine apply_i
end interface


contains

end module differential_operator_mod
