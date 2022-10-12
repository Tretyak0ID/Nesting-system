module field_mod
use domain_mod, only: domain_t
implicit none

type, public :: field_t

  real(kind=8), allocatable :: f(:, :)

contains

  procedure :: init
  procedure :: init_on_domain

end type field_t

contains

  subroutine init(this, sindx, eindx, sindy, eindy)
    integer(kind=4),  intent(in)  :: sindx, eindx, sindy, eindy
    class  (field_t), intent(out) :: this

    allocate(this%f(sindx : eindx, sindy : eindy))
  end subroutine init

  subroutine init_on_domain(this, domain)
    class  (field_t), intent(out) :: this
    type(domain_t),   intent(in)  :: domain

    allocate(this%f(0 : domain%nx, 0 : domain%ny))
  end subroutine init_on_domain

end module field_mod
