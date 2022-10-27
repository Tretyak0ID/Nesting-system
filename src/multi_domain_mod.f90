module multi_domain_mod

  use domain_mod, only : domain_t
  implicit none

  type, public :: multi_domain_t

    type(domain_t)               :: global_domain
    type(domain_t), allocatable  :: subdomains(:, :)
    integer(kind=4)              :: num_sub_x, num_sub_y
    integer(kind=4), allocatable :: degree_condensation(:, :)

  contains

    procedure, public :: init => init_multi_domain

  end type multi_domain_t

contains

  subroutine init_multi_domain(this, global_domain, num_sub_x, num_sub_y, degree_condensation)

    class(multi_domain_t),        intent(inout) :: this
    class(domain_t),              intent(in)    :: global_domain
    integer(kind=4),              intent(in)    :: num_sub_x, num_sub_y
    integer(kind=4), allocatable, intent(in)    :: degree_condensation(:, :)

    integer(kind=4) :: i, j, n, m, ie, je
    real(kind=8)    :: xs, xe, ys, ye

    this%global_domain = global_domain
    this%num_sub_x     = num_sub_x
    this%num_sub_y     = num_sub_y

    allocate(this%subdomains(1 : num_sub_x, 1 : num_sub_y))
    allocate(this%degree_condensation(1 : num_sub_x, 1 : num_sub_y))

    this%degree_condensation = degree_condensation

    do n = 1, num_sub_x
      do m = 1, num_sub_y
        xs = global_domain%xs + (n - 1) * abs(global_domain%xe - global_domain%xs) / num_sub_x
        xe = global_domain%xs + n * abs(global_domain%xe - global_domain%xs) / num_sub_x
        ys = global_domain%ys + (m - 1) * abs(global_domain%ye - global_domain%ys) / num_sub_y
        ye = global_domain%ys + m * abs(global_domain%ye - global_domain%ys) / num_sub_y
        ie = this%degree_condensation(n, m) * global_domain%nx
        je = this%degree_condensation(n, m) * global_domain%ny
        call this%subdomains(n, m)%init(xs, xe, 0, ie, ys, ye, 0, je)
      end do
    end do

  end subroutine init_multi_domain

end module multi_domain_mod
