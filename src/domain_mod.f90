module domain_mod
implicit none

type, public :: domain_t

  real   (kind=8) :: xs, xe, ys, ye, dx, dy
  integer(kind=4) :: is, ie, js, je, nx, ny
  real   (kind=8), allocatable :: domain_x(:), domain_y(:)

contains

  procedure :: init

end type domain_t

contains

  subroutine init(this, xs, xe, is, ie, ys, ye, js, je)
    class  (domain_t), intent(out) :: this
    real   (kind=8),   intent(in)  :: xs, xe, ys, ye
    integer(kind=4),   intent(in)  :: is, ie, js, je
    integer(kind=4) :: i, j

    this%xs = xs
    this%xe = xe
    this%is = is
    this%ie = ie
    this%nx = ie - is
    this%ys = ys
    this%ye = ye
    this%js = js
    this%je = je
    this%ny = je - js
    this%dx = (xe - xs) / this%nx
    this%dy = (ye - ys) / this%ny

    allocate(this%domain_x(0 : this%nx))
    allocate(this%domain_y(0 : this%ny))

    do i = is, ie
      this%domain_x(i) = xs + i * this%dx
    end do

    do j = js, je
      this%domain_y(j) = ys + j * this%dy
    end do

  end subroutine init

end module domain_mod
