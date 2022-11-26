program test_5_barotropic_instability
use initial_conditions_mod,            only : swm_baratropic_instability
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use swe_advective_operator_mod,        only : swe_advective_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t
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
  type(sbp21_t)                    :: sbp21
  type(sbp42_t)                    :: sbp42
  type(central2_t)                 :: central2
  type(central4_t)                 :: central4
  class(timescheme_t), allocatable :: timescheme
  integer(kind=4),     allocatable :: deg(:, :)

  real(kind=8)    :: LX     = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii
  real(kind=8)    :: H_MEAN = 10.0_8 ** 4.0_8
  integer(kind=4) :: Nt     = 180 * 128, t
  real(kind=8)    :: T_max  = 40.0_8 * 3600.0_8 * 24.0_8, dt

  allocate(deg(1:2, 1:1))
  deg(1, 1) = 2
  deg(2, 1) = 1
  dt = T_max / Nt

  sbp21%name    = 'sbp21'
  sbp42%name    = 'sbp42'
  central2%name = 'cent2'
  central4%name = 'cent4'

  call domain%init(0.0_8, LX, 0, 256, 0.0_8, LY, 0, 256)
  call multi_domain%init(domain, 2, 1, deg)
  call curl%init(multi_domain)
  call state%h%init(multi_domain)
  call state%u%init(multi_domain)
  call state%v%init(multi_domain)
  call op%init(sbp42, central4, multi_domain)

  call create_timescheme(timescheme, state, 'rk4')
  call swm_baratropic_instability(state, multi_domain, H_MEAN, 50.0_8)

  do t = 0, Nt
    if (mod(t, 100) == 0) print *, 'step: ',  t
    if (mod(t, 50) == 0) call calc_curl(curl, state%u, state%v, multi_domain, sbp42, central4)
    !call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test3h_left.dat', t + 1)
    !call write_field(state%h%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test3h_right.dat', t + 1)
    if (mod(t, 50) == 0) call write_field(curl%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test5curl_left.dat', t / 50 + 1)
    if (mod(t, 50) == 0) call write_field(curl%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test5curl_right.dat', t / 50 + 1)
    call timescheme%step(state, op, multi_domain, dt)
  end do

end program test_5_barotropic_instability