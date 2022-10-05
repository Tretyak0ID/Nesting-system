module domain_mod
implicit none

type, public :: domain_t

  real   (kind=8) :: xs, xe, ys, ye, dx, dy
  integer(kind=4) :: nx, ny
  real   (kind=8), allocatable :: mesh_x(:), mesh_y(:)

contains

  procedure :: init

end type domain_t

contains

  subroutine init(this, xs, xe, nx, ys, ye, ny)
    class  (domain_t), intent(out) :: this
    real   (kind=8),   intent(in)  :: xs, xe, ys, ye
    integer(kind=4),   intent(in)  :: nx, ny
    integer(kind=4) :: i, j

    this%xs = xs
    this%xe = xe
    this%nx = nx
    this%ys = ys
    this%ye = ye
    this%ny = ny
    this%dx = (xe - xs) / nx
    this%dy = (ye - ys) / ny

    allocate(this%mesh_x(0 : nx))
    allocate(this%mesh_y(0 : ny))

    do i = 0, nx
      this%mesh_x(i) = xs + i * this%dx
    end do

    do j = 0, ny
      this%mesh_y(j) = ys + j * this%dy
    end do

  end subroutine init

end module domain_mod
