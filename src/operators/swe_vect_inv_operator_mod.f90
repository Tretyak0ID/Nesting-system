module swe_vect_inv_operator_mod
use stvec_swe_mod,             only: stvec_swe_t
use stvec_mod,                 only: stvec_t
use operator_mod,              only: operator_t
use differential_operator_mod, only: differential_operator_t
use multi_grid_field_mod,      only: multi_grid_field_t
use multi_domain_mod,          only: multi_domain_t
use grad_mod,                  only: calc_grad
use div_mod,                   only: calc_div
use curl_mod,                  only: calc_curl
use const_mod,                 only: Earth_grav, pcori
implicit none

  type, public, extends(operator_t) :: swe_vect_inv_operator_t
  !RHS operator for shallow water equation in vector invariant form
    class(differential_operator_t), allocatable :: diff_opx
    class(differential_operator_t), allocatable :: diff_opy
    !work fields for operator
    type(multi_grid_field_t) :: div, gx, gy, curl !for h field
    type(multi_grid_field_t) :: kin_energy, gh_kin_energy
    type(multi_grid_field_t) :: hu, hv !mass fluxes in continuty eq

  contains

    procedure :: init => init_vecinv
    procedure :: apply => apply_swe_vect_inv

  end type swe_vect_inv_operator_t

contains

  subroutine init_vecinv(this, diff_opx, diff_opy, multi_domain)

    class(swe_vect_inv_operator_t), intent(inout) :: this
    class(differential_operator_t), intent(in)    :: diff_opx, diff_opy
    class(multi_domain_t),          intent(in)    :: multi_domain

    call this%div%init(multi_domain)
    call this%curl%init(multi_domain)
    call this%gx%init(multi_domain)
    call this%gy%init(multi_domain)
    call this%kin_energy%init(multi_domain)
    call this%gh_kin_energy%init(multi_domain)
    call this%hu%init(multi_domain)
    call this%hv%init(multi_domain)

    this%diff_opx = diff_opx
    this%diff_opy = diff_opy

  end subroutine init_vecinv

  subroutine apply_swe_vect_inv(this, out, in, multi_domain)

    class (swe_vect_inv_operator_t),  intent(inout) :: this
    class (stvec_t),                  intent(inout) :: out
    class (stvec_t),                  intent(inout) :: in
    type  (multi_domain_t),           intent(in)    :: multi_domain

    select type (out)
    class is (stvec_swe_t)
      select type(in)
      class is (stvec_swe_t)

        call this%hu%assign(-1.0_8, in%h, in%u, multi_domain)
        call this%hv%assign(-1.0_8, in%h, in%v, multi_domain)

        call this%kin_energy%assign(0.50_8, in%u, in%u, multi_domain)
        call this%kin_energy%update(0.50_8, in%v, in%v, multi_domain)
        call this%gh_kin_energy%assign(Earth_grav, in%h, 1.0_8, this%kin_energy, multi_domain)

        call calc_div(this%div, this%hu, this%hv, multi_domain, this%diff_opx, this%diff_opy)
        call calc_curl(this%curl, in%u, in%v, multi_domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gx, this%gy, this%gh_kin_energy, multi_domain, this%diff_opx, this%diff_opy)

        !u field
        call out%u%assign(-1.0_8, this%gx, multi_domain)
        call out%u%update(pcori, in%v, multi_domain)
        call out%u%update(1.0_8, this%curl, in%v, multi_domain)
        !v field
        call out%v%assign(-1.0_8, this%gy, multi_domain)
        call out%v%update(-pcori, in%u, multi_domain)
        call out%v%update(-1.0_8, this%curl, in%u, multi_domain)
        !h field
        call out%h%assign(1.0_8, this%div, multi_domain)

        !forcing
        !call out%u%update(-0.0_8 * pcori, multi_domain)
        !call out%v%update( 10.0_8 * pcori, multi_domain)

      class default
      end select
    class default
    end select

  end subroutine

end module swe_vect_inv_operator_mod
