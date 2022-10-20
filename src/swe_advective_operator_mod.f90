module swe_advective_operator_mod
use stvec_swe_mod,             only: stvec_swe_t
use stvec_mod,                 only: stvec_t
use operator_mod,              only: operator_t
use differential_operator_mod, only: differential_operator_t
use field_mod,                 only: field_t
use domain_mod,                  only: domain_t
use grad_mod,                  only: calc_grad
use div_mod,                   only: calc_div
use const_mod,                 only: Earth_grav, pcori
implicit none

  type, public, extends(operator_t) :: swe_advective_operator_t

    class(differential_operator_t), allocatable :: diff_opx
    class(differential_operator_t), allocatable :: diff_opy
    !work fields for operator
    type(field_t)                 :: div, gx, gy, curl !for h field
    type(field_t)                 :: gux, guy, gvx, gvy !velocity field grad
    type(field_t)                 :: hu, hv, divmf !mass fluxes in continuty eq

  contains

    procedure :: init => init_swe
    procedure :: apply => apply_swe_advective

  end type swe_advective_operator_t

contains

  subroutine init_swe(this, diff_opx, diff_opy)

    class(swe_advective_operator_t), intent(inout) :: this
    class(differential_operator_t),  intent(in)    :: diff_opx, diff_opy
    this%diff_opx = diff_opx
    this%diff_opy = diff_opy

  end subroutine init_swe

  subroutine apply_swe_advective(this, out, in, domain)

    class (swe_advective_operator_t), intent(inout) :: this
    class (stvec_t),                  intent(inout) :: out
    class (stvec_t),                  intent(inout) :: in
    type  (domain_t),                   intent(in)    :: domain

    integer       :: i, j

    select type (out)
    class is (stvec_swe_t)
      select type(in)
      class is (stvec_swe_t)

        call this%divmf%init_on_domain(domain)
        call this%gx%init_on_domain(domain)
        call this%gy%init_on_domain(domain)
        call this%gux%init_on_domain(domain)
        call this%guy%init_on_domain(domain)
        call this%gvx%init_on_domain(domain)
        call this%gvy%init_on_domain(domain)

        call this%hu%assign(1.0_8, in%h, in%u, domain)
        call this%hv%assign(1.0_8, in%h, in%v, domain)

        call calc_div(this%divmf, this%hu, this%hv, domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gx, this%gy, in%h, domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gux, this%guy, in%u, domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gvx, this%gvy, in%v, domain, this%diff_opx, this%diff_opy)

        !dudt calculate
        call out%u%assign(-1.0_8, in%u, this%gux, domain)
        call out%u%update(-1.0_8, in%v, this%guy, domain)
        call out%u%update(-1.0_8 * Earth_grav, this%gx, domain)
        call out%u%update(pcori, in%v, domain)
        !dvdt calculate
        call out%v%assign(-1.0_8, in%u, this%gvx, domain)
        call out%v%update(-1.0_8, in%v, this%gvy, domain)
        call out%v%update(-1.0_8 * Earth_grav, this%gy, domain)
        call out%v%update(-1._8 * pcori, in%u, domain)
        !dhdt calculate
        call out%h%assign(-1.0_8, this%divmf, domain)

      class default
      end select
    class default
    end select

  end subroutine apply_swe_advective

end module swe_advective_operator_mod
