module field_mod
implicit none

type, public :: field_t

  real(kind=8), allocatable :: f(:, :)

contains

  procedure :: init

end type field_t

contains

  subroutine init(this, sindx, eindx, sindy, eindy)
    integer(kind=4),  intent(in)  :: sindx, eindx, sindy, eindy
    class  (field_t), intent(out) :: this

    allocate(this%f(sindx : eindx, sindy : eindy))
  end subroutine init

end module field_mod
