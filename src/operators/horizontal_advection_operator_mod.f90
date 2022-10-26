module horizontal_advection_operator_mod
  use stvec_swe_mod,             only: stvec_swe_t
  use stvec_mod,                 only: stvec_t
  use operator_mod,              only: operator_t
  use differential_operator_mod, only: differential_operator_t
  use field_mod,                 only: field_t
  use domain_mod,                only: domain_t
  use grad_mod,                  only: calc_grad
  use div_mod,                   only: calc_div
  use const_mod,                 only: Earth_grav, pcori
implicit none

type, public, extends(operator_t) :: horizontal_advection_operator_t

  class(differential_operator_t), allocatable :: diff_opx
  class(differential_operator_t), allocatable :: diff_opy
  !work fields for operator
  type(field_t)                 :: div, gx, gy

contains

  procedure :: init => init_ha
  procedure :: apply => apply_ha

end type horizontal_advection_operator_t

contains

  subroutine init_ha(this, diff_opx, diff_opy)

    class(horizontal_advection_operator_t), intent(inout) :: this
    class(differential_operator_t),         intent(in)    :: diff_opx, diff_opy
    this%diff_opx = diff_opx
    this%diff_opy = diff_opy

  end subroutine init_ha

  subroutine apply_ha(this, out, in, domain)

    class (horizontal_advection_operator_t), intent(inout) :: this
    class (stvec_t),                         intent(inout) :: out
    class (stvec_t),                         intent(inout) :: in
    type  (domain_t),                        intent(in)    :: domain

    select type (out)
    class is (stvec_swe_t)
      select type(in)
      class is (stvec_swe_t)

        call this%gx%init_on_domain(domain)
        call this%gy%init_on_domain(domain)

        call calc_grad(this%gx, this%gy, in%h, domain, this%diff_opx, this%diff_opy)

        call out%h%assign(-1.0_8, in%u, this%gx, domain)
        call out%h%update(-1.0_8, in%v, this%gy, domain)

      class default
      end select
    class default
    end select

  end subroutine apply_ha

end module horizontal_advection_operator_mod
