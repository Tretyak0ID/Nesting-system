module field_mod
use mesh_mod, only: mesh_t
implicit none

type, public :: field_t

  real(kind=8), allocatable :: f(:, :)

contains

  procedure, public :: init
  procedure, public :: init_on_mesh
  !update
  procedure, public :: update_s1       => update_field_s1 !v = v + s1*unity
  procedure, public :: update_s1v1     => update_field_s1v1 !v = v + s1*v1
  procedure, public :: update_s1v1s2v2 => update_field_s1v1s2v2 !v = v + s1*v1+s2*v2
  generic           :: update          => update_s1, update_s1v1, update_s1v1s2v2
  !assign
  procedure, public :: assign_s1       => assign_field_s1 !v = s1
  procedure, public :: assign_v1       => assign_field_v1 !v = v1
  procedure, public :: assign_s1v1     => assign_field_s1v1 !v = s1*v1
  procedure, public :: assign_s1v1s2v2 => assign_field_s1v1s2v2 !v = s1*v1+s2*v2
  generic           :: assign          => assign_s1v1, assign_s1, assign_v1, assign_s1v1s2v2

end type field_t

contains

subroutine init(this, sindx, eindx, sindy, eindy)

  integer(kind=4),  intent(in)  :: sindx, eindx, sindy, eindy
  class  (field_t), intent(out) :: this

  allocate(this%f(sindx : eindx, sindy : eindy))

end subroutine init

subroutine init_on_mesh(this, mesh)

  class  (field_t), intent(out) :: this
  type(mesh_t),   intent(in)  :: mesh

  allocate(this%f(mesh%sindx : mesh%eindx, mesh%sindy : mesh%eindy))

end subroutine init_on_mesh

!update
subroutine update_field_s1(this, scalar1, mesh)

  class(field_t), intent(inout) :: this
  real(kind=8),   intent(in)    :: scalar1
  type(mesh_t), intent(in)    :: mesh
  integer(kind=8) :: i, j

  do i = mesh%sindx, mesh%eindx
    do j = mesh%sindy, mesh%eindy
      this%f(i, j) = this%f(i, j) + scalar1
    end do
  end do

end subroutine update_field_s1

subroutine update_field_s1v1(this, scalar1, v1, mesh)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  type(mesh_t),      intent(in)    :: mesh
  integer(kind=8) :: i, j

  do i = mesh%sindx, mesh%eindx
    do j = mesh%sindy, mesh%eindy
      this%f(i, j) = this%f(i, j) + scalar1 * v1%f(i, j)
    end do
  end do

end subroutine update_field_s1v1

subroutine update_field_s1v1s2v2(this, scalar1, v1, scalar2, v2, mesh)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  type(field_t),       intent(in)    :: v2
  type(mesh_t),      intent(in)    :: mesh
  integer(kind=8) :: i, j

  do i = mesh%sindx, mesh%eindx
    do j = mesh%sindy, mesh%eindy
      this%f(i, j) = this%f(i, j) + scalar1 * v1%f(i, j) + scalar2 * v2%f(i, j)
    end do
  end do

end subroutine update_field_s1v1s2v2

!assign
subroutine assign_field_s1(this, scalar1, mesh)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(mesh_t),      intent(in)    :: mesh
  integer(kind=8) :: i, j

  do i = mesh%sindx, mesh%eindx
    do j = mesh%sindy, mesh%eindy
      this%f(i, j) = scalar1
    end do
  end do

end subroutine assign_field_s1

subroutine assign_field_v1(this, v1, mesh)

  class(field_t),      intent(inout) :: this
  type(field_t),       intent(in)    :: v1
  type(mesh_t),      intent(in)    :: mesh
  integer(kind=8) :: i, j

  do i = mesh%sindx, mesh%eindx
    do j = mesh%sindy, mesh%eindy
      this%f(i, j) = v1%f(i, j)
    end do
  end do

end subroutine assign_field_v1

subroutine assign_field_s1v1(this, scalar1, v1, mesh)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  type(mesh_t),      intent(in)    :: mesh
  integer(kind=8) :: i, j

  do i = mesh%sindx, mesh%eindx
    do j = mesh%sindy, mesh%eindy
      this%f(i, j) = scalar1 * v1%f(i, j)
    end do
  end do

end subroutine assign_field_s1v1

subroutine assign_field_s1v1s2v2(this, scalar1, v1, scalar2, v2, mesh)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  type(field_t),       intent(in)    :: v2
  type(mesh_t),      intent(in)    :: mesh
  integer(kind=8) :: i, j

  do i = mesh%sindx, mesh%eindx
    do j = mesh%sindy, mesh%eindy
      this%f(i, j) = scalar1 * v1%f(i, j) + scalar2 * v2%f(i, j)
    end do
  end do

end subroutine assign_field_s1v1s2v2

end module field_mod
