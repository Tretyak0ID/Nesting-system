module field_mod
use domain_mod, only: domain_t
implicit none

type, public :: field_t

  real(kind=8), allocatable :: f(:, :)
  integer(kind=4)           :: is, ie, js, je

contains

  procedure, public :: init
  procedure, public :: init_on_domain
  procedure, public :: init_real
  procedure, public :: copy
  procedure, public :: create_similar
  !update
  procedure, public :: update_s1       => update_field_s1 !v = v + s1*unity
  procedure, public :: update_s1v1     => update_field_s1v1 !v = v + s1*v1
  procedure, public :: update_s1v1s2v2 => update_field_s1v1s2v2 !v = v + s1*v1+s2*v2
  procedure, public :: update_s1v1v2   => update_field_s1v1v2 !V = v + v1*v2
  generic           :: update          => update_s1, update_s1v1v2, update_s1v1, update_s1v1s2v2
  !assign
  procedure, public :: assign_s1       => assign_field_s1 !v = s1
  procedure, public :: assign_v1       => assign_field_v1 !v = v1
  procedure, public :: assign_s1v1     => assign_field_s1v1 !v = s1*v1
  procedure, public :: assign_s1v1v2   => assign_field_s1v1v2 !v = v1*v2
  procedure, public :: assign_s1v1s2v2 => assign_field_s1v1s2v2 !v = s1*v1+s2*v2
  generic           :: assign          => assign_s1v1, assign_s1v1v2, assign_s1, assign_v1, assign_s1v1s2v2

end type field_t

contains

subroutine init(this, is, ie, js, je)

  integer(kind=4),  intent(in)  :: is, ie, js, je
  class  (field_t), intent(out) :: this

  integer(kind=4) :: i, j

  allocate(this%f(is : ie, js : je))

  this%is = is
  this%ie = ie
  this%js = js
  this%je = je

  do i = is, ie
    do j = js, je
      this%f(i, j) = 0.0_8
    end do
  end do

end subroutine init

subroutine init_on_domain(this, domain)

  class  (field_t), intent(out) :: this
  type(domain_t),   intent(in)  :: domain

  integer(kind=8) :: i, j

  allocate(this%f(domain%is : domain%ie, domain%js : domain%je))

  this%is = domain%is
  this%ie = domain%ie
  this%js = domain%js
  this%je = domain%je

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = 0.0_8
    end do
  end do

end subroutine init_on_domain

subroutine init_real(this, r, domain)

  class  (field_t), intent(inout) :: this
  real(kind=8),     intent(in)    :: r
  type(domain_t),   intent(in)    :: domain

  integer(kind=8) :: i, j

  allocate(this%f(domain%is : domain%ie, domain%js : domain%je))

  this%is = domain%is
  this%ie = domain%ie
  this%js = domain%js
  this%je = domain%je

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = r
    end do
  end do

end subroutine init_real

subroutine copy(this, fin)

  class(field_t), intent(inout) :: this
  class(field_t), intent(in)    :: fin

  integer(kind=4) :: i, j

  allocate(this%f(fin%is : fin%ie, fin%js : fin%je))

  this%is = fin%is
  this%ie = fin%ie
  this%js = fin%js
  this%je = fin%je

  do j = fin%js, fin%je
    do i = fin%is, fin%ie
      this%f(i, j) = fin%f(i, j)
    end do
  end do

end subroutine copy

subroutine create_similar(this, destination)

  class(field_t), intent(in)    :: this
  class(field_t), intent(inout) :: destination

  integer(kind=4) :: i, j

  call destination%init(this%is, this%ie, this%js, this%je)

end subroutine create_similar

!update
subroutine update_field_s1(this, scalar1, domain)

  class(field_t), intent(inout) :: this
  real(kind=8),   intent(in)    :: scalar1
  type(domain_t),   intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = this%f(i, j) + scalar1
    end do
  end do

end subroutine update_field_s1

subroutine update_field_s1v1(this, scalar1, v1, domain)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  type(domain_t),      intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = this%f(i, j) + scalar1 * v1%f(i, j)
    end do
  end do

end subroutine update_field_s1v1

subroutine update_field_s1v1v2(this, scalar1, f1, f2, domain)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: f1, f2
  type(domain_t),        intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = this%f(i, j) + scalar1 * f2%f(i, j) * f1%f(i, j)
    end do
  end do

end subroutine update_field_s1v1v2

subroutine update_field_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  type(field_t),       intent(in)    :: v2
  type(domain_t),      intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = this%f(i, j) + scalar1 * v1%f(i, j) + scalar2 * v2%f(i, j)
    end do
  end do

end subroutine update_field_s1v1s2v2

!assign
subroutine assign_field_s1(this, scalar1, domain)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(domain_t),      intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = scalar1
    end do
  end do

end subroutine assign_field_s1

subroutine assign_field_v1(this, v1, domain)

  class(field_t),      intent(inout) :: this
  type(field_t),       intent(in)    :: v1
  type(domain_t),      intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = v1%f(i, j)
    end do
  end do

end subroutine assign_field_v1

subroutine assign_field_s1v1(this, scalar1, v1, domain)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  type(domain_t),      intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = scalar1 * v1%f(i, j)
    end do
  end do

end subroutine assign_field_s1v1

subroutine assign_field_s1v1v2(this, scalar1, f1, f2, domain)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: f1, f2
  type(domain_t),        intent(in)  :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = scalar1 * f2%f(i, j) * f1%f(i, j)
    end do
  end do

end subroutine assign_field_s1v1v2

subroutine assign_field_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

  class(field_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(field_t),       intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  type(field_t),       intent(in)    :: v2
  type(domain_t),      intent(in)    :: domain
  integer(kind=8) :: i, j

  do j = domain%js, domain%je
    do i = domain%is, domain%ie
      this%f(i, j) = scalar1 * v1%f(i, j) + scalar2 * v2%f(i, j)
    end do
  end do

end subroutine assign_field_s1v1s2v2

end module field_mod
