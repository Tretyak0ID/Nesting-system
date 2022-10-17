module stvec_swe_mod
use field_mod,  only: field_t
use mesh_mod, only: mesh_t
use stvec_mod,  only: stvec_t
implicit none

type, extends(stvec_t), public :: stvec_swe_t

  type(field_t) :: h, u, v

contains
  !update
  procedure, public :: update_s1v1 !st = st + s1*v1
  procedure, public :: update_s1v1s2v2 !st = st + s1*v1 + s2*v2
  !assign
  procedure, public :: assign_s1v1 !st = s1*v1
  procedure, public :: assign_s1v1s2v2 !st = s1*v1 + s2*v2

end type stvec_swe_t

contains

subroutine update_s1v1(this, scalar1, v1, mesh)

  class(stvec_swe_t), intent(inout) :: this
  real(kind=8),       intent(in)    :: scalar1
  class(stvec_t), intent(in)    :: v1
  type(mesh_t),     intent(in)    :: mesh

  select type (v1)
  class is (stvec_swe_t)
    call this%h%update(scalar1, v1%h, mesh)
    call this%u%update(scalar1, v1%u, mesh)
    call this%v%update(scalar1, v1%v, mesh)
  class default
  end select

end subroutine update_s1v1

subroutine update_s1v1s2v2(this, scalar1, v1, scalar2, v2, mesh)

  class(stvec_swe_t),  intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),  intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  class(stvec_t),  intent(in)    :: v2
  type(mesh_t),      intent(in)    :: mesh

  select type (v1)
  class is (stvec_swe_t)
    select type (v2)
    class is (stvec_swe_t)
      call this%h%update(scalar1, v1%h, scalar2, v2%h, mesh)
      call this%u%update(scalar1, v1%u, scalar2, v2%u, mesh)
      call this%v%update(scalar1, v1%v, scalar2, v2%v, mesh)
    class default
    end select
  class default
  end select

end subroutine update_s1v1s2v2

subroutine assign_s1v1(this, scalar1, v1, mesh)

  class(stvec_swe_t),  intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),  intent(in)    :: v1
  type(mesh_t),      intent(in)    :: mesh

  select type (v1)
  class is (stvec_swe_t)
    call this%h%assign(scalar1, v1%h, mesh)
    call this%u%assign(scalar1, v1%u, mesh)
    call this%v%assign(scalar1, v1%v, mesh)
  class default
  end select

end subroutine assign_s1v1

subroutine assign_s1v1s2v2(this, scalar1, v1, scalar2, v2, mesh)

  class(stvec_swe_t),  intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),  intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  class(stvec_t),  intent(in)    :: v2
  type(mesh_t),      intent(in)    :: mesh

  select type (v1)
  class is (stvec_swe_t)
    select type (v2)
    class is (stvec_swe_t)
      call this%h%assign(scalar1, v1%h, scalar2, v2%h, mesh)
      call this%u%assign(scalar1, v1%u, scalar2, v2%u, mesh)
      call this%v%assign(scalar1, v1%v, scalar2, v2%v, mesh)
    class default
    end select
  class default
  end select

end subroutine assign_s1v1s2v2

end module stvec_swe_mod
