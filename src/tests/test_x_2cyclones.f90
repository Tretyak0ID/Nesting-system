program test_5_barotropic_instability
use initial_conditions_mod,            only : set_two_cyclones_initial_conditions
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use diffusion_operator_mod,            only : diffusion_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t, sbp21_2_t
use central_differential_operator_mod, only : central2_t, central4_t
use curl_mod,                          only : calc_curl
use timescheme_mod,                    only : timescheme_t
use timescheme_factory_mod,            only : create_timescheme
use rk4_mod,                           only : rk4_t
use field_mod,                         only : field_t
use multi_grid_field_mod,              only : multi_grid_field_t
use domain_mod,                        only : domain_t
use multi_domain_mod,                  only : multi_domain_t
use stvec_swe_mod,                     only : stvec_swe_t
use const_mod,                         only : Earth_radii, Earth_grav, pcori, pi
use read_write_mod,                    only : write_field
implicit none

  type(multi_domain_t)             :: multi_domain
  type(domain_t)                   :: domain
  type(multi_grid_field_t)         :: curl
  type(stvec_swe_t)                :: state
  type(swe_vect_inv_operator_t)    :: op
  type(diffusion_operator_t)       :: diffusion
  type(sbp21_t)                    :: sbp21
  type(sbp42_t)                    :: sbp42
  type(sbp21_2_t)                  :: sbp21_2
  type(central2_t)                 :: central2
  type(central4_t)                 :: central4
  class(timescheme_t), allocatable :: timescheme, explicit_Euler
  integer(kind=4),     allocatable :: deg(:, :)
    real(kind=8),      allocatable :: coefs(:, :)

  real(kind=8), parameter :: LX     = 4000e3_8, LY = 4000e3_8
  integer(kind=4), parameter :: Nx = 512, Ny = 512
  real(kind=8), parameter    :: H_MEAN = 1e4_8
  real(kind=8), parameter    :: T_max  = 128.0_8*3600.0_8, dt = 25.0_8, tau_wr = 3600.0_8
  integer(kind=4), parameter :: Nt     = int(T_max/dt)
  integer(kind=4) :: t, nzap, irec

  nzap = int(tau_wr/dt)

  allocate(deg(1:1, 1:1))
  deg(1, 1) = 1
  !deg(2, 1) = 1

  sbp21%name    = 'sbp21_1'
  sbp42%name    = 'sbp42_1'
  central2%name = 'cent2_1'
  central4%name = 'cent4_1'
  sbp21_2%name  = 'sbp21_2'

  call domain%init(0.0_8, LX, 0, Nx, 0.0_8, LY, 0, Ny)
  call multi_domain%init(domain, 1, 1, deg)
  call curl%init(multi_domain)
  call state%h%init(multi_domain)
  call state%u%init(multi_domain)
  call state%v%init(multi_domain)
  call op%init(sbp42, central4, multi_domain)

  allocate(coefs(1:1, 1:1))
  coefs(1, 1) = multi_domain%subdomains(1, 1)%dx ** 2.0_8 / 2.0_8 / dt
  !coefs(2, 1) = multi_domain%subdomains(2, 1)%dx ** 2.0_8 / 2.0_8 / dt
  call diffusion%init(sbp21_2, coefs, multi_domain)

  call create_timescheme(timescheme, state, 'rk4')
  call create_timescheme(explicit_Euler, state, 'explicit_Euler')
  call set_two_cyclones_initial_conditions(state, multi_domain, H_MEAN)

  do t = 0, Nt
     print *, "step", t
    if(mod(t,nzap)==0) then
        irec = t/nzap+1
        print *, "write t=", t,irec
        call write_field(state%u%subfields(1, 1), multi_domain%subdomains(1, 1), './data/two_cyclones_u.dat', irec)
        call write_field(state%v%subfields(1, 1), multi_domain%subdomains(1, 1), './data/two_cyclones_v.dat', irec)
        call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/two_cyclones_h.dat', irec)
        call calc_curl(curl, state%u, state%v, multi_domain, sbp42, central4)
        call write_field(curl%subfields(1, 1), multi_domain%subdomains(1, 1), './data/two_cyclones_curl.dat', irec)
    end if
    call timescheme%step(state, op, multi_domain, dt)
    call explicit_Euler%step(state, diffusion, multi_domain, dt)
  end do

end program test_5_barotropic_instability
