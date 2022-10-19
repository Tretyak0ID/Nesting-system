module operator_mod
use stvec_mod,           only: stvec_t
use mesh_mod,          only: mesh_t
implicit none

  type, abstract :: operator_t

    contains

      procedure(apply_i), deferred :: apply

  end type operator_t

abstract interface
  subroutine apply_i(this, out, in, mesh)
    import operator_t, mesh_t, stvec_t

    class    (operator_t), intent(inout) :: this
    class    (stvec_t),    intent(inout) :: out
    class    (stvec_t),    intent(inout) :: in
    type     (mesh_t),   intent(in)      :: mesh

  end subroutine apply_i
end interface


contains

end module operator_mod
