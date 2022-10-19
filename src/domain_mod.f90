module domain_mod

use mesh_mod, only : mesh_t
implicit none

type, public :: domain_t

  type(mesh_t), public    :: mesh
  real(kind=8), public    :: et, dt
  integer(kind=8), public :: nt

contains

  procedure :: init

end type domain_t

contains

  subroutine init(this, et, nt, mesh)
    class(domain_t), intent(inout) :: this
    real(kind=8),    intent(in)    :: et
    integer(kind=8), intent(in)    :: nt
    type(mesh_t),    intent(in)    :: mesh

    this%et   = et
    this%nt   = nt
    this%dt   = et / nt
    this%mesh = mesh

  end subroutine init

end module domain_mod
