module central_differential_operator_mod
use differential_operator_mod, only : differential_operator_t
use field_mod,                 only: field_t
use mesh_mod,                only: mesh_t
implicit none

  type, extends(differential_operator_t) :: central2_t

  contains

    procedure :: apply => apply_central2

  end type central2_t

  type, extends(differential_operator_t) :: central4_t

  contains

    procedure :: apply => apply_central4

  end type central4_t

contains

  subroutine apply_central2(this, out, in, mesh, direction)
    class    (central2_t), intent(in)    :: this
    type     (field_t),    intent(inout) :: out
    type     (field_t),    intent(in)    :: in
    type     (mesh_t),   intent(in)    :: mesh
    character(len=1),      intent(in)    :: direction
    integer  (kind=8)                    :: i, j

    if (direction == 'x') then
      do i = mesh%sindx + 1, mesh%eindx - 1
        do j = mesh%sindy, mesh%eindy
          out%f(i, j) = (in%f(i + 1, j) - in%f(i - 1, j)) / (2.0_8 * mesh%dx)
        end do
      end do
      do j = mesh%sindy, mesh%eindy
        out%f(0, j)         = (in%f(1, j) - in%f(mesh%nx - 1, j)) / (2.0_8 * mesh%dx)
        out%f(mesh%nx, j) = out%f(0, j)
      end do

    else if(direction == 'y') then
      do i = mesh%sindx, mesh%eindx
        do j = mesh%sindy + 1, mesh%eindy - 1
          out%f(i, j) = (in%f(i, j + 1) - in%f(i, j - 1)) / (2.0_8 * mesh%dy)
        end do
      end do
      do i = mesh%sindx, mesh%eindx
        out%f(i, 0)             = (in%f(i, 1) - in%f(i, mesh%ny - 1)) / (2.0_8 * mesh%dy)
        out%f(i, mesh%ny)      = out%f(i, 0)
      end do

    else
      print *, 'Error in diff_central4. Wrong direction value'
    end if

  end subroutine apply_central2

  subroutine apply_central4(this, out, in, mesh, direction)
    class    (central4_t), intent(in)    :: this
    type     (field_t),    intent(inout) :: out
    type     (field_t),    intent(in)    :: in
    type     (mesh_t),   intent(in)    :: mesh
    character(len=1),      intent(in)    :: direction
    integer  (kind=8)                    :: i, j

    if (direction == 'x') then
      do i = mesh%sindx + 2, mesh%eindx - 2
        do j = mesh%sindy, mesh%eindy
          out%f(i, j) = (in%f(i - 2, j) - 8.0_8 * in%f(i - 1, j) + 8.0_8 * in%f(i + 1, j) - in%f(i + 2, j)) / (12.0_8 * mesh%dx)
        end do
      end do
      do j = mesh%sindy, mesh%eindy
        out%f(0, j)             = (in%f(mesh%nx - 2, j) - 8.0_8 * in%f(mesh%nx - 1, j) + 8.0_8 * in%f(1, j) - in%f(2, j)) / (12.0_8 * mesh%dx)
        out%f(1, j)             = (in%f(mesh%nx - 1, j) - 8.0_8 * in%f(0, j) + 8.0_8 * in%f(2, j) - in%f(3, j)) / (12.0_8 * mesh%dx)
        out%f(mesh%nx - 1, j) = (in%f(mesh%nx - 3, j) - 8.0_8 * in%f(mesh%nx - 2, j) + 8.0_8 * in%f(0, j) - in%f(1, j)) / (12.0_8 * mesh%dx)
        out%f(mesh%nx, j)      = out%f(0, j)
      end do

    else if(direction == 'y') then
      do i = mesh%sindx, mesh%eindx
        do j = mesh%sindy + 2, mesh%eindy - 2
          out%f(i, j) = (in%f(i, j - 2) - 8.0_8 * in%f(i, j - 1) + 8.0_8 * in%f(i, j + 1) - in%f(i, j + 2)) / (12.0_8 * mesh%dy)
        end do
      end do
      do i = mesh%sindx, mesh%eindx
        out%f(i, 0)             = (in%f(i, mesh%ny - 2) - 8.0_8 * in%f(i, mesh%ny - 1) + 8.0_8 * in%f(i, 1) - in%f(i, 2)) / (12.0_8 * mesh%dy)
        out%f(i, 1)             = (in%f(i, mesh%ny - 1) - 8.0_8 * in%f(i, 0) + 8.0_8 * in%f(i, 2) - in%f(i, 3)) / (12.0_8 * mesh%dy)
        out%f(i, mesh%ny - 1) = (in%f(i, mesh%ny - 3) - 8.0_8 * in%f(i, mesh%ny - 2) + 8.0_8 * in%f(i, 0) - in%f(i, 1)) / (12.0_8 * mesh%dy)
        out%f(i, mesh%ny)      = out%f(i, 0)
      end do

    else
      print *, 'Error in diff_central4. Wrong direction value'
    end if

  end subroutine apply_central4

end module central_differential_operator_mod
