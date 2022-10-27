module multi_grid_field_mod
  use field_mod,        only : field_t
  use domain_mod,       only : domain_t
  use multi_domain_mod, only : multi_domain_t
  implicit none

  type, public :: multi_grid_field_t

    type(field_t), allocatable   :: subfields(:, :)
    integer(kind=4)              :: num_sub_x, num_sub_y

  contains

    procedure, public :: init => init_multi_grid_field
    procedure, public :: copy => copy_multi_grid_field
    !update
    procedure, public :: update_s1       => update_multi_grid_field_s1 !v = v + s1*unity
    procedure, public :: update_s1v1     => update_multi_grid_field_s1v1 !v = v + s1*v1
    procedure, public :: update_s1v1s2v2 => update_multi_grid_field_s1v1s2v2 !v = v + s1*v1+s2*v2
    procedure, public :: update_s1v1v2   => update_multi_grid_field_s1v1v2 !V = v + v1*v2
    generic           :: update          => update_s1, update_s1v1v2, update_s1v1, update_s1v1s2v2
    !assign
    procedure, public :: assign_s1       => assign_multi_grid_field_s1 !v = s1
    procedure, public :: assign_v1       => assign_multi_grid_field_v1 !v = v1
    procedure, public :: assign_s1v1     => assign_multi_grid_field_s1v1 !v = s1*v1
    procedure, public :: assign_s1v1v2   => assign_multi_grid_field_s1v1v2 !v = v1*v2
    procedure, public :: assign_s1v1s2v2 => assign_multi_grid_field_s1v1s2v2 !v = s1*v1+s2*v2
    generic           :: assign          => assign_s1v1, assign_s1v1v2, assign_s1, assign_v1, assign_s1v1s2v2

  end type multi_grid_field_t

contains

  subroutine init_multi_grid_field(this, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) ::  n, m

    allocate(this%subfields(1 : multi_domain%num_sub_x, 1 : multi_domain%num_sub_y))
    this%num_sub_x = multi_domain%num_sub_x
    this%num_sub_y = multi_domain%num_sub_y

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        call this%subfields(n, m)%init_on_domain(multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine init_multi_grid_field

  subroutine copy_multi_grid_field(this, fin)

    class(multi_grid_field_t), intent(inout) :: this
    class(multi_grid_field_t), intent(in)    :: fin
    integer(kind=4) ::  n, m

    allocate(this%subfields(1 : fin%num_sub_x, 1 : fin%num_sub_y))
    this%num_sub_x = fin%num_sub_x
    this%num_sub_y = fin%num_sub_y

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%copy(fin%subfields(n, m))
      end do
    end do

  end subroutine copy_multi_grid_field

  !update
  subroutine update_multi_grid_field_s1(this, scalar1, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    real(kind=8),              intent(in)    :: scalar1
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%update(scalar1, multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine update_multi_grid_field_s1

  subroutine update_multi_grid_field_s1v1(this, scalar1, v1, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    real(kind=8),              intent(in)    :: scalar1
    type(multi_grid_field_t),  intent(in)    :: v1
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%update(scalar1, v1%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine update_multi_grid_field_s1v1

  subroutine update_multi_grid_field_s1v1v2(this, scalar1, f1, f2, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    real(kind=8),              intent(in)    :: scalar1
    type(multi_grid_field_t),  intent(in)    :: f1, f2
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%update(scalar1, f1%subfields(n, m), f2%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine update_multi_grid_field_s1v1v2

  subroutine update_multi_grid_field_s1v1s2v2(this, scalar1, v1, scalar2, v2, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    real(kind=8),              intent(in)    :: scalar1
    type(multi_grid_field_t),  intent(in)    :: v1
    real(kind=8),              intent(in)    :: scalar2
    type(multi_grid_field_t),  intent(in)    :: v2
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%update(scalar1, v1%subfields(n, m), scalar2, v2%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine update_multi_grid_field_s1v1s2v2

  !assign
  subroutine assign_multi_grid_field_s1(this, scalar1, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    real(kind=8),              intent(in)    :: scalar1
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%assign(scalar1, multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine assign_multi_grid_field_s1

  subroutine assign_multi_grid_field_v1(this, v1, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    type(multi_grid_field_t),  intent(in)    :: v1
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%assign(v1%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine assign_multi_grid_field_v1

  subroutine assign_multi_grid_field_s1v1(this, scalar1, v1, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    real(kind=8),              intent(in)    :: scalar1
    type(multi_grid_field_t),  intent(in)    :: v1
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%assign(scalar1, v1%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine assign_multi_grid_field_s1v1

  subroutine assign_multi_grid_field_s1v1v2(this, scalar1, f1, f2, multi_domain)

    class(multi_grid_field_t), intent(inout) :: this
    real(kind=8),              intent(in)    :: scalar1
    type(multi_grid_field_t),  intent(in)    :: f1, f2
    type(multi_domain_t),      intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%assign(scalar1, f1%subfields(n, m), f2%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine assign_multi_grid_field_s1v1v2

  subroutine assign_multi_grid_field_s1v1s2v2(this, scalar1, v1, scalar2, v2, multi_domain)

    class(multi_grid_field_t),intent(inout) :: this
    real(kind=8),             intent(in)    :: scalar1
    type(multi_grid_field_t), intent(in)    :: v1
    real(kind=8),             intent(in)    :: scalar2
    type(multi_grid_field_t), intent(in)    :: v2
    type(multi_domain_t),     intent(in)    :: multi_domain
    integer(kind=4) :: n, m

    do n = 1, this%num_sub_x
      do m = 1, this%num_sub_y
        call this%subfields(n, m)%assign(scalar1, v1%subfields(n, m), scalar2, v2%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

  end subroutine assign_multi_grid_field_s1v1s2v2

end module multi_grid_field_mod
