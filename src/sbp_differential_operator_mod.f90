module sbp_differential_operator_mod
use differential_operator_mod, only : differential_operator_t
use field_mod,  only: field_t
use mesh_mod, only: mesh_t
implicit none

  type, extends(differential_operator_t) :: sbp21_t

  contains

    procedure :: apply => apply_sbp21

  end type sbp21_t

  type, extends(differential_operator_t) :: sbp42_t

  contains

    procedure :: apply => apply_sbp42

  end type sbp42_t

contains

  subroutine apply_sbp21(this, out, in, mesh, direction)

    class    (sbp21_t),  intent(in)    :: this
    type     (field_t),  intent(inout) :: out
    type     (field_t),  intent(in)    :: in
    type     (mesh_t), intent(in)    :: mesh
    character(len=1),    intent(in)    :: direction
    integer  (kind=8)                  :: i, j

    if (direction == 'x') then
      do i = mesh%sindx + 1, mesh%eindx - 1
        do j = mesh%sindy, mesh%eindy
          out%f(i, j) = (in%f(i + 1, j) - in%f(i - 1, j)) / 2.0_8 / mesh%dx
        end do
      end do
      do j = mesh%sindy, mesh%eindy
        out%f(mesh%nx, j) = (in%f(mesh%eindx, j) - in%f(mesh%eindx - 1, j)) / mesh%dx
        out%f(0, j)         = (in%f(1, j) - in%f(0, j)) / mesh%dx
      end do

    else if (direction == 'y') then
      do i = mesh%sindx, mesh%eindx
        do j = mesh%sindy + 1, mesh%eindy - 1
          out%f(i, j) = (in%f(i, j + 1) - in%f(i, j - 1)) / 2.0_8 / mesh%dy
        end do
      end do
      do i = mesh%sindx, mesh%eindx
        out%f(i, mesh%ny) = (in%f(i, mesh%eindy) - in%f(i, mesh%eindy - 1)) / mesh%dy
        out%f(i, 0)         = (in%f(i, 1) - in%f(i, 0)) / mesh%dy
      end do

    else
      print *, 'Error in diff_sbp21. Wrong direction value'
    end if

  end subroutine apply_sbp21

  subroutine apply_sbp42(this, out, in, mesh, direction)
    class    (sbp42_t),  intent(in)    :: this
    type     (field_t),  intent(inout) :: out
    type     (field_t),  intent(in)    :: in
    type     (mesh_t), intent(in)    :: mesh
    character(len=1),    intent(in)    :: direction
    integer  (kind=8)                  :: i, j

    if (direction == 'x') then
      do i = mesh%sindx + 4, mesh%eindx - 3
        do j = mesh%sindy, mesh%eindy
          out%f(i, j) = (in%f(i - 2, j) - 8.0_8 * in%f(i - 1, j) + 8.0_8 * in%f(i + 1, j) - in%f(i + 2, j)) / (12.0_8 * mesh%dx)
        end do
      end do
      do j = mesh%sindy, mesh%eindy
        out%f(0, j) = (-24.0_8 / 17.0_8 * in%f(0, j) + 59.0_8 / 34.0_8 * in%f(1, j) - 4.0_8 / 17.0_8 * in%f(2, j) - 3.0_8 / 34.0_8 * in%f(3, j)) / (mesh%dx)
        out%f(1, j) = (-in%f(0, j) + in%f(2, j)) / (2.0_8 * mesh%dx)
        out%f(2, j) = (4.0_8 / 43.0_8 * in%f(0, j) - 59.0_8 / 86.0_8 * in%f(1, j) + 59.0_8 / 86.0_8 * in%f(3, j) - 4.0_8 / 43.0_8 * in%f(4, j)) / (mesh%dx)
        out%f(3, j) = (3.0_8 / 98.0_8 * in%f(0, j) - 59.0_8 / 98.0_8 * in%f(2, j) + 32.0_8 / 49.0_8 * in%f(4, j) - 4.0_8 / 49.0_8 * in%f(5, j)) / (mesh%dx)

        out%f(mesh%eindx, j) = (24.0_8 / 17.0_8 * in%f(mesh%eindx, j) - 59.0_8 / 34.0_8 * in%f(mesh%eindx - 1, j) + 4.0_8 / 17.0_8 * in%f(mesh%eindx - 2, j) + 3.0_8 / 34.0_8 * in%f(mesh%eindx - 3, j)) / (mesh.dx)
        out%f(mesh%eindx - 1, j) = (in%f(mesh%eindx, j) - in%f(mesh%eindx - 2, j)) / (2.0_8 * mesh%dx)
        out%f(mesh%eindx - 2, j) = (-4.0_8 / 43.0_8 * in%f(mesh%eindx, j) + 59.0_8 / 86.0_8 * in%f(mesh%eindx - 1, j) - 59.0_8 / 86.0_8 * in%f(mesh%eindx - 3, j) + 4.0_8 / 43.0_8 * in%f(mesh%eindx - 4 , j)) / (mesh%dx)
        out%f(mesh%eindx - 3, j) = (-3.0_8 / 98.0_8 * in%f(mesh%eindx, j) + 59.0_8 / 98.0_8 * in%f(mesh%eindx - 2, j) - 32.0_8 / 49.0_8 * in%f(mesh%eindx - 4, j) + 4.0_8 / 49.0_8 * in%f(mesh%eindx - 5, j)) / (mesh%dx)
      end do

    else if(direction == 'y') then
      do i = mesh%sindx, mesh%eindx
        do j = mesh%sindy + 4, mesh%eindy - 3
          out%f(i, j) = (in%f(i, j - 2) - 8.0_8 * in%f(i, j - 1) + 8.0_8 * in%f(i, j + 1) - in%f(i, j + 2)) / (12.0_8 * mesh%dy)
        end do
      end do
      do i = mesh%sindx, mesh%eindx
        out%f(i, 0) = (-24.0_8 / 17.0_8 * in%f(i, 0) + 59.0_8 / 34.0_8 * in%f(i, 1) - 4.0_8 / 17.0_8 * in%f(i, 2) - 3.0_8 / 34.0_8 * in%f(i, 3)) / (mesh%dy)
        out%f(i, 1) = (-in%f(i, 0) + in%f(i, 2)) / (2.0_8 * mesh%dy)
        out%f(i, 2) = (4.0_8 / 43.0_8 * in%f(i, 0) - 59.0_8 / 86.0_8 * in%f(i, 1) + 59.0_8 / 86.0_8 * in%f(i, 3) - 4.0_8 / 43.0_8 * in%f(i, 4)) / (mesh%dy)
        out%f(i, 3) = (3.0_8 / 98.0_8 * in%f(i, 0) - 59.0_8 / 98.0_8 * in%f(i, 2) + 32.0_8 / 49.0_8 * in%f(i, 4) - 4.0_8 / 49.0_8 * in%f(i, 5)) / (mesh%dy)

        out%f(i, mesh%eindy) = (24.0_8 / 17.0_8 * in%f(i, mesh%eindy) - 59.0_8 / 34.0_8 * in%f(i, mesh%eindy - 1) + 4.0_8 / 17.0_8 * in%f(i, mesh%eindy - 2) + 3.0_8 / 34.0_8 * in%f(i, mesh%eindy - 3)) / (mesh.dy)
        out%f(i, mesh%eindy - 1) = (in%f(i, mesh%eindy) - in%f(i, mesh%eindy - 2)) / (2.0_8 * mesh%dy)
        out%f(i, mesh%eindy - 2) = (-4.0_8 / 43.0_8 * in%f(i, mesh%eindy) + 59.0_8 / 86.0_8 * in%f(i, mesh%eindy - 1) - 59.0_8 / 86.0_8 * in%f(i, mesh%eindy - 3) + 4.0_8 / 43.0_8 * in%f(i, mesh%eindy - 4 )) / (mesh%dy)
        out%f(i, mesh%eindy - 3) = (-3.0_8 / 98.0_8 * in%f(i, mesh%eindy) + 59.0_8 / 98.0_8 * in%f(i, mesh%eindy - 2) - 32.0_8 / 49.0_8 * in%f(i, mesh%eindy - 4) + 4.0_8 / 49.0_8 * in%f(i, mesh%eindy - 5)) / (mesh%dy)
      end do

    else
      print *, 'Error in diff_sbp42. Wrong direction value'
    end if
  end subroutine apply_sbp42

end module sbp_differential_operator_mod
