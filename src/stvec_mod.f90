module stvec_mod
use field_mod,            only: field_t
use domain_mod,           only: domain_t
use multi_grid_field_mod, only : multi_grid_field_t
use multi_domain_mod,     only : multi_domain_t
implicit none

type, abstract, public :: stvec_t

contains

  procedure, public :: copy           => copy_stvec
  procedure, public :: create_similar => create_similar_stvec
    !updates
  procedure, public :: update_s1v1     => update_stvec_s1v1 !v = v + s1*v1
  procedure, public :: update_s1v1s2v2 => update_stvec_s1v1s2v2 !v = v + s1*v1+s2*v2
  procedure, public :: update_s1v1v2   => update_stvec_s1v1v2 !v = v + v1*v2
  generic           :: update          => update_s1v1, update_s1v1v2, update_s1v1s2v2
    !assign
  procedure, public :: assign_s1         => assign_stvec_s1 !v = s1
  procedure, public :: assign_s1v1       => assign_stvec_s1v1 !v = s1*v1
  procedure, public :: assign_s1v1s2v2   => assign_stvec_s1v1s2v2 !v = s1*v1+s2*v2
  procedure, public :: assign_s1v1v2     => assign_stvec_s1v1v2     !v = v1*v2
  generic           :: assign            => assign_s1v1, assign_s1, assign_s1v1s2v2, assign_s1v1v2

end type stvec_t

contains

subroutine copy_stvec(this, fin, multi_domain)

  class(stvec_t), intent(inout)     :: this
  class(stvec_t), intent(in)        :: fin
  class(multi_domain_t), intent(in) :: multi_domain

end subroutine copy_stvec

subroutine create_similar_stvec(this, destination)

  class(stvec_t), intent(in)    :: this
  class(stvec_t), intent(inout) :: destination

end subroutine create_similar_stvec

subroutine update_stvec_s1v1(this, scalar1, v1, multi_domain)

  class(stvec_t),        intent(inout) :: this
  real(kind=8),          intent(in)    :: scalar1
  class(stvec_t),        intent(in)    :: v1
  class(multi_domain_t), intent(in)    :: multi_domain

end subroutine update_stvec_s1v1

subroutine update_stvec_s1v1v2(this, scalar1, v1, v2, multi_domain)

  class(stvec_t),        intent(inout) :: this
  real(kind=8),          intent(in)    :: scalar1
  class(stvec_t),        intent(in)    :: v1
  class(stvec_t),        intent(in)    :: v2
  class(multi_domain_t), intent(in)    :: multi_domain

end subroutine update_stvec_s1v1v2

subroutine update_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, multi_domain)

  class(stvec_t),        intent(inout) :: this
  real(kind=8),          intent(in)    :: scalar1
  class(stvec_t),        intent(in)    :: v1
  real(kind=8),          intent(in)    :: scalar2
  class(stvec_t),        intent(in)    :: v2
  class(multi_domain_t), intent(in)    :: multi_domain

end subroutine update_stvec_s1v1s2v2

subroutine assign_stvec_s1(this, scalar1, multi_domain)

  class(stvec_t),        intent(inout) :: this
  real(kind=8),          intent(in)    :: scalar1
  class(multi_domain_t), intent(in)    :: multi_domain

end subroutine assign_stvec_s1

subroutine assign_stvec_s1v1(this, scalar1, v1, multi_domain)

  class(stvec_t),        intent(inout) :: this
  real(kind=8),          intent(in)    :: scalar1
  class(stvec_t),        intent(in)    :: v1
  class(multi_domain_t), intent(in)    :: multi_domain

end subroutine assign_stvec_s1v1

subroutine assign_stvec_s1v1v2(this, scalar1, v1, v2, multi_domain)

  class(stvec_t),        intent(inout) :: this
  real(kind=8),          intent(in)    :: scalar1
  class(stvec_t),        intent(in)    :: v1, v2
  class(multi_domain_t), intent(in)    :: multi_domain

end subroutine assign_stvec_s1v1v2

subroutine assign_stvec_s1v1s2v2(this, scalar1, v1, scalar2, v2, multi_domain)

  class(stvec_t),        intent(inout) :: this
  real(kind=8),          intent(in)    :: scalar1
  class(stvec_t),        intent(in)    :: v1
  real(kind=8),          intent(in)    :: scalar2
  class(stvec_t),        intent(in)    :: v2
  class(multi_domain_t), intent(in)    :: multi_domain

end subroutine assign_stvec_s1v1s2v2

end module stvec_mod
