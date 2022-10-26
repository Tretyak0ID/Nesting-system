module swe_vect_inv_operator_mod
use stvec_swe_mod,             only: stvec_swe_t
use stvec_mod,                 only: stvec_t
use operator_mod,              only: operator_t
use differential_operator_mod, only: differential_operator_t
use field_mod,                 only: field_t
use domain_mod,                only: domain_t
use grad_mod,                  only: calc_grad
use div_mod,                   only: calc_div
use curl_mod,                  only: calc_curl
use const_mod,                 only: Earth_grav, pcori
implicit none

  type, public, extends(operator_t) :: swe_vect_inv_operator_t

    class(differential_operator_t), allocatable :: diff_opx
    class(differential_operator_t), allocatable :: diff_opy
    !work fields for operator
    type(field_t)                 :: div, gx, gy, curl !for h field
    type(field_t)                 :: kin_energy, gh_kin_energy
    type(field_t)                 :: hu, hv !mass fluxes in continuty eq

  contains

    procedure :: init => init_vecinv
    procedure :: apply => apply_swe_vect_inv

  end type swe_vect_inv_operator_t

contains

  subroutine init_vecinv(this, diff_opx, diff_opy)

    class(swe_vect_inv_operator_t), intent(inout) :: this
    class(differential_operator_t), intent(in)    :: diff_opx, diff_opy
    this%diff_opx = diff_opx
    this%diff_opy = diff_opy

  end subroutine init_vecinv

  subroutine apply_swe_vect_inv(this, out, in, domain)

    class (swe_vect_inv_operator_t),  intent(inout) :: this
    class (stvec_t),                  intent(inout) :: out
    class (stvec_t),                  intent(inout) :: in
    type  (domain_t),                 intent(in)    :: domain

    select type (out)
    class is (stvec_swe_t)
      select type(in)
      class is (stvec_swe_t)

        call this%div%init_on_domain(domain)
        call this%curl%init_on_domain(domain)
        call this%gx%init_on_domain(domain)
        call this%gy%init_on_domain(domain)
        call this%kin_energy%init_on_domain(domain)
        call this%gh_kin_energy%init_on_domain(domain)
        call this%hu%init_on_domain(domain)
        call this%hv%init_on_domain(domain)

        call this%hu%assign(1.0_8, in%h, in%u, domain)
        call this%hv%assign(1.0_8, in%h, in%v, domain)

        call this%kin_energy%assign(0.5_8, in%u, in%u, domain)
        call this%kin_energy%update(0.5_8, in%v, in%v, domain)
        call this%gh_kin_energy%assign(Earth_grav, in%h, 1.0_8, this%kin_energy, domain)

        call calc_div(this%div, this%hu, this%hv, domain, this%diff_opx, this%diff_opy)
        call calc_curl(this%curl, in%u, in%v, domain, this%diff_opx, this%diff_opy)
        call calc_grad(this%gx, this%gy, this%gh_kin_energy, domain, this%diff_opx, this%diff_opy)

        !u field
        call out%u%assign(-1.0_8, this%gx, domain)
        call out%u%update(pcori, in%v, domain)
        call out%u%update(1.0_8, this%curl, in%v, domain)
        !v field
        call out%v%assign(-1.0_8, this%gy, domain)
        call out%v%update(-pcori, in%u, domain)
        call out%v%update(-1.0_8, this%curl, in%u, domain)
        !h field
        call out%h%assign(-1.0_8, this%div, domain)

      class default
      end select
    class default
    end select

  end subroutine

end module swe_vect_inv_operator_mod
