module mesh_mod
implicit none

type, public :: mesh_t

  real   (kind=8) :: xs, xe, ys, ye, dx, dy
  integer(kind=4) :: sindx, eindx, sindy, eindy, nx, ny
  real   (kind=8), allocatable :: mesh_x(:), mesh_y(:)

contains

  procedure :: init

end type mesh_t

contains

  subroutine init(this, xs, xe, sindx, eindx, ys, ye, sindy, eindy)
    class  (mesh_t), intent(out) :: this
    real   (kind=8),   intent(in)  :: xs, xe, ys, ye
    integer(kind=4),   intent(in)  :: sindx, eindx, sindy, eindy
    integer(kind=4) :: i, j

    this%xs = xs
    this%xe = xe
    this%sindx = sindx
    this%eindx = eindx
    this%nx = eindx - sindx
    this%ys = ys
    this%ye = ye
    this%sindy = sindy
    this%eindy = eindy
    this%ny = eindy - sindy
    this%dx = (xe - xs) / this%nx
    this%dy = (ye - ys) / this%ny

    allocate(this%mesh_x(0 : this%nx))
    allocate(this%mesh_y(0 : this%ny))

    do i = sindx, eindx
      this%mesh_x(i) = xs + i * this%dx
    end do

    do j = sindy, eindy
      this%mesh_y(j) = ys + j * this%dy
    end do

  end subroutine init

end module mesh_mod
