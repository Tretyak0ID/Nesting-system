program test_3_geostriphic_cyclone_discontinuum
use initial_conditions_mod,            only : swm_geostrophic_cyclone
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use diffusion_operator_mod,            only : diffusion_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t, sbp21_2_t
use central_differential_operator_mod, only : central2_t, central4_t
use curl_mod,                          only : calc_curl
use grad_mod,                          only : calc_grad
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
  type(multi_grid_field_t)         :: curl, gux, guy, gvx, gvy
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
  real(kind=8),        allocatable :: coefs(:, :)

  real(kind=8)    :: LX     = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii
  real(kind=8)    :: H_MEAN = 10.0_8 ** 4.0_8
  integer(kind=4) :: Nt     = 180 * 64, t
  real(kind=8)    :: T_max  = 40.0_8 * 3600.0_8 * 24.0_8, dt, u0 = 10.0_8, v0 = 0.0_8

  allocate(deg(1:2, 1:1))
  deg(1, 1) = 2
  deg(2, 1) = 1
  dt = T_max / Nt

  sbp21%name    = 'sbp21_1'
  sbp42%name    = 'sbp42_1'
  central2%name = 'cent2_1'
  central4%name = 'cent4_1'
  sbp21_2%name  = 'sbp21_2'

  call domain%init(0.0_8, LX, 0, 128, 0.0_8, LY, 0, 128)
  call multi_domain%init(domain, 2, 1, deg)
  call curl%init(multi_domain)
  call state%h%init(multi_domain)
  call state%u%init(multi_domain)
  call state%v%init(multi_domain)
  call state%u%create_similar(gux)
  call state%u%create_similar(guy)
  call state%v%create_similar(gvx)
  call state%v%create_similar(gvy)
  call op%init(sbp42, central4, multi_domain)

  allocate(coefs(1:2, 1:1))
  coefs(1, 1) = multi_domain%subdomains(1, 1)%dx ** 2.0_8 / dt / 2.0_8
  coefs(2, 1) = multi_domain%subdomains(2, 1)%dx ** 2.0_8 / dt / 2.0_8
  call diffusion%init(sbp21_2, coefs, multi_domain)

  call create_timescheme(timescheme, state, 'rk4')
  call create_timescheme(explicit_Euler, state, 'explicit_Euler')
  call swm_geostrophic_cyclone(state, multi_domain, H_MEAN, 22.0_8 * 10.0_8 ** (-3.0_8), 1000000.0_8, 0.0_8, u0, v0, 'const_discontinuum')

  do t = 0, Nt
    if (mod(t, 100) == 0) print *, 'step: ',  t
    if (mod(t, 50) == 0) call calc_curl(curl, state%u, state%v, multi_domain, sbp42, central4)
    !call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test3h_left.dat', t + 1)
    !call write_field(state%h%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test3h_right.dat', t + 1)
    if (mod(t, 50) == 0) call write_field(curl%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test4curl_left.dat', t / 50 + 1)
    if (mod(t, 50) == 0) call write_field(curl%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test4curl_right.dat', t / 50 + 1)
    call timescheme%step(state, op, multi_domain, dt)
    call explicit_Euler%step(state, diffusion, multi_domain, dt)
  end do

end program test_3_geostriphic_cyclone_discontinuum
