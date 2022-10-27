module multi_grid_field_mod
  use field_mod,        only : field_t
  use domain_mod,       only : domain_t
  use multi_domain_mod, only : multi_domain_t
  implicit none

  type, public :: multi_grid_field_t

    type(field_t), allocatable   :: subfields(:, :)

  contains

    procedure, public :: init => init_multi_grid_field

  end type multi_grid_field_t

contains

  subroutine init_multi_grid_field(this, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    type(multi_domain_t),      intent(in)    :: multi_domain

    integer(kind=4) ::  n, m
    allocate(this%subfields(1 : multi_domain%num_sub_x, 1 : multi_domain%num_sub_y))

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        call this%subfields(n, m)%init_on_domain(multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine init_multi_grid_field

end module multi_grid_field_mod
