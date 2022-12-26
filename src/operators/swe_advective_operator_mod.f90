module swe_advective_operator_mod
use stvec_swe_mod,             only: stvec_swe_t
use stvec_mod,                 only: stvec_t
use operator_mod,              only: operator_t
use differential_operator_mod, only: differential_operator_t
use multi_grid_field_mod,      only: multi_grid_field_t
use multi_domain_mod,          only: multi_domain_t
use grad_mod,                  only: calc_grad
use div_mod,                   only: calc_div
use const_mod,                 only: Earth_grav, pcori
implicit none

  type, public, extends(operator_t) :: swe_advective_operator_t
  !RHS operator for shallow water equation advective form
    class(differential_operator_t), allocatable :: diff_opx
    class(differential_operator_t), allocatable :: diff_opy
    !work fields for operator
    type(multi_grid_field_t) :: div, gx, gy, curl !for h field
    type(multi_grid_field_t) :: gux, guy, gvx, gvy !velocity field grad
    type(multi_grid_field_t) :: hu, hv, divmf !mass fluxes in continuty eq

  contains

    procedure :: init => init_swe
    procedure :: apply => apply_swe_advective

  end type swe_advective_operator_t

contains

  subroutine init_swe(this, diff_opx, diff_opy, multi_domain)

    class(swe_advective_operator_t), intent(inout) :: this
    class(differential_operator_t),  intent(in)    :: diff_opx, diff_opy
    class(multi_domain_t),           intent(in)    :: multi_domain

    call this%div%init(multi_domain)
    call this%hu%init(multi_domain)
    call this%hv%init(multi_domain)
    call this%gx%init(multi_domain)
    call this%gy%init(multi_domain)
    call this%gux%init(multi_domain)
    call this%guy%init(multi_domain)
    call this%gvx%init(multi_domain)
    call this%gvy%init(multi_domain)

    this%diff_opx = diff_opx
    this%diff_opy = diff_opy

  end subroutine init_swe

  subroutine apply_swe_advective(this, out, in, multi_domain)

    class (swe_advective_operator_t), intent(inout) :: this
    class (stvec_t),                  intent(inout) :: out
    class (stvec_t),                  intent(inout) :: in
    type  (multi_domain_t),           intent(in)    :: multi_domain

    select type (out)
    class is (stvec_swe_t)
      select type(in)
      class is (stvec_swe_t)

        call this%hu%assign(-1.0_8, in%h, in%u, multi_domain)
        call this%hv%assign(-1.0_8, in%h, in%v, multi_domain)

        call calc_div(this%div, this%hu, this%hv, multi_domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gx, this%gy, in%h, multi_domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gux, this%guy, in%u, multi_domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gvx, this%gvy, in%v, multi_domain, this%diff_opx, this%diff_opy)

        !dudt calculate
        call out%u%assign(-1.0_8, in%u, this%gux, multi_domain)
        call out%u%update(-1.0_8, in%v, this%guy, multi_domain)
        call out%u%update(-1.0_8 * Earth_grav, this%gx, multi_domain)
        call out%u%update(pcori, in%v, multi_domain)
        !dvdt calculate
        call out%v%assign(-1.0_8, in%u, this%gvx, multi_domain)
        call out%v%update(-1.0_8, in%v, this%gvy, multi_domain)
        call out%v%update(-1.0_8 * Earth_grav, this%gy, multi_domain)
        call out%v%update(-1.0_8 * pcori, in%u, multi_domain)
        !dhdt calculate
        call out%h%assign(1.0_8, this%div, multi_domain)

      class default
      end select
    class default
    end select

  end subroutine apply_swe_advective

end module swe_advective_operator_mod
