module stvec_swe_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
use stvec_mod,  only: stvec_t
implicit none

type, extends(stvec_t), public :: stvec_swe_t

  type(field_t) :: h, u, v

contains

  procedure, public :: create_similar
  procedure, public :: copy
  !update
  procedure, public :: update_s1v1 !st = st + s1*v1
  procedure, public :: update_s1v1v2 !st = st + s1*v1*v2
  procedure, public :: update_s1v1s2v2 !st = st + s1*v1 + s2*v2
  !assign
  procedure, public :: assign_s1v1 !st = s1*v1
  procedure, public :: assign_s1v1v2 !st = s1*v1*v2
  procedure, public :: assign_s1v1s2v2 !st = s1*v1 + s2*v2

end type stvec_swe_t

contains

subroutine create_similar(this, destination)
  class(stvec_swe_t),          intent(in)    :: this
  class(stvec_t), allocatable, intent(inout) :: destination

  allocate(stvec_swe_t :: destination)
  select type (destination)
  class is (stvec_swe_t)
    call destination%h%init(this%h%is, this%h%ie, this%h%js, this%h%je)
    call destination%u%init(this%u%is, this%u%ie, this%u%js, this%u%je)
    call destination%v%init(this%v%is, this%v%ie, this%v%js, this%v%je)
  end select

end subroutine create_similar

subroutine copy(this, fin, domain)

  class(stvec_swe_t), intent(inout) :: this
  class(stvec_t),     intent(in)    :: fin
  class(domain_t),     intent(in)    :: domain

  select type (fin)
  class is (stvec_swe_t)
    call this%u%copy(fin%u)
    call this%v%copy(fin%v)
    call this%h%copy(fin%h)
  class default
  end select

end subroutine copy

subroutine update_s1v1(this, scalar1, v1, domain)

  class(stvec_swe_t), intent(inout) :: this
  real(kind=8),       intent(in)    :: scalar1
  class(stvec_t),     intent(in)    :: v1
  class(domain_t),     intent(in)    :: domain

  select type (v1)
  class is (stvec_swe_t)
    call this%h%update(scalar1, v1%h, domain)
    call this%u%update(scalar1, v1%u, domain)
    call this%v%update(scalar1, v1%v, domain)
  class default
  end select

end subroutine update_s1v1

subroutine update_s1v1v2(this, scalar1, v1, v2, domain)

  class(stvec_swe_t), intent(inout) :: this
  real(kind=8),       intent(in)    :: scalar1
  class(stvec_t),     intent(in)    :: v1, v2
  class(domain_t),    intent(in)   :: domain

  select type (v1)
  class is (stvec_swe_t)
    select type (v2)
    class is (stvec_swe_t)
      call this%h%update(scalar1, v1%h, v2%h, domain)
      call this%u%update(scalar1, v1%u, v2%u, domain)
      call this%v%update(scalar1, v1%v, v2%v, domain)
    class default
    end select
  class default
  end select

end subroutine update_s1v1v2

subroutine update_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

  class(stvec_swe_t),  intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),      intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  class(stvec_t),      intent(in)    :: v2
  class(domain_t),     intent(in)  :: domain

  select type (v1)
  class is (stvec_swe_t)
    select type (v2)
    class is (stvec_swe_t)
      call this%h%update(scalar1, v1%h, scalar2, v2%h, domain)
      call this%u%update(scalar1, v1%u, scalar2, v2%u, domain)
      call this%v%update(scalar1, v1%v, scalar2, v2%v, domain)
    class default
    end select
  class default
  end select

end subroutine update_s1v1s2v2

subroutine assign_s1v1(this, scalar1, v1, domain)

  class(stvec_swe_t),  intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),      intent(in)    :: v1
  class(domain_t),     intent(in)    :: domain

  select type (v1)
  class is (stvec_swe_t)
    call this%h%assign(scalar1, v1%h, domain)
    call this%u%assign(scalar1, v1%u, domain)
    call this%v%assign(scalar1, v1%v, domain)
  class default
  end select

end subroutine assign_s1v1

subroutine assign_s1v1v2(this, scalar1, v1, v2, domain)

  class(stvec_swe_t), intent(inout) :: this
  real(kind=8),       intent(in)    :: scalar1
  class(stvec_t),     intent(in)    :: v1, v2
  class(domain_t),    intent(in)  :: domain

  select type (v1)
  class is (stvec_swe_t)
    select type (v2)
    class is (stvec_swe_t)
      call this%h%assign(scalar1, v1%h, v2%h, domain)
      call this%u%assign(scalar1, v1%u, v2%u, domain)
      call this%v%assign(scalar1, v1%v, v2%v, domain)
    class default
    end select
  class default
  end select

end subroutine assign_s1v1v2

subroutine assign_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

  class(stvec_swe_t),  intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),      intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  class(stvec_t),      intent(in)    :: v2
  class(domain_t),    intent(in)  :: domain

  select type (v1)
  class is (stvec_swe_t)
    select type (v2)
    class is (stvec_swe_t)
      call this%h%assign(scalar1, v1%h, scalar2, v2%h, domain)
      call this%u%assign(scalar1, v1%u, scalar2, v2%u, domain)
      call this%v%assign(scalar1, v1%v, scalar2, v2%v, domain)
    class default
    end select
  class default
  end select

end subroutine assign_s1v1s2v2

end module stvec_swe_mod
