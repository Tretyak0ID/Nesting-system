module stvec_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
implicit none

type, abstract, public :: stvec_t

contains
    !updates
  procedure, public :: update_s1v1     => update_stvec_s1v1 !v = v + s1*v1
  procedure, public :: update_s1v1s2v2 => update_stvec_s1v1s2v2 !v = v + s1*v1+s2*v2
  generic           :: update          => update_s1v1, update_s1v1s2v2
    !assign
  procedure, public :: assign_s1         => assign_stvec_s1 !v = s1
  procedure, public :: assign_s1v1       => assign_stvec_s1v1 !v = s1*v1
  procedure, public :: assign_s1v1s2v2   => assign_stvec_s1v1s2v2 !v = s1*v1+s2*v2
  generic           :: assign            => assign_s1v1, assign_s1, assign_s1v1s2v2

end type stvec_t

contains

subroutine update_stvec_s1v1(this, scalar1, v1, domain)

  class(stvec_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),      intent(in)    :: v1
  type(domain_t),      intent(in)    :: domain

end subroutine update_stvec_s1v1

subroutine update_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

  class(stvec_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),      intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  class(stvec_t),      intent(in)    :: v2
  type(domain_t),      intent(in)    :: domain

end subroutine update_stvec_s1v1s2v2

subroutine assign_stvec_s1(this, scalar1, domain)

  class(stvec_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  type(domain_t),      intent(in)    :: domain

end subroutine assign_stvec_s1

subroutine assign_stvec_s1v1(this, scalar1, v1, domain)

  class(stvec_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),      intent(in)    :: v1
  type(domain_t),      intent(in)    :: domain

end subroutine assign_stvec_s1v1

subroutine assign_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, domain)

  class(stvec_t),      intent(inout) :: this
  real(kind=8),        intent(in)    :: scalar1
  class(stvec_t),      intent(in)    :: v1
  real(kind=8),        intent(in)    :: scalar2
  class(stvec_t),      intent(in)    :: v2
  type(domain_t),      intent(in)    :: domain

end subroutine assign_stvec_s1v1s2v2

end module stvec_mod