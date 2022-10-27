program test_2_geostriphic_balance
use initial_conditions_mod,            only : swm_geostrophic_balance
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use central_differential_operator_mod, only : central2_t, central4_t
use timescheme_mod,                    only : timescheme_t
use timescheme_factory_mod,            only : create_timescheme
use rk4_mod,                           only : rk4_t
use domain_mod,                        only : domain_t
use stvec_swe_mod,                     only : stvec_swe_t
use const_mod,                         only : Earth_radii, Earth_grav, pcori, pi
use read_write_mod,                    only : write_field
implicit none

  type(domain_t)                   :: domain
  type(stvec_swe_t)                :: state
  type(swe_vect_inv_operator_t)    :: op
  type(central2_t)                 :: central2
  type(central4_t)                 :: central4
  class(timescheme_t), allocatable :: timescheme

  real(kind=8)    :: LX     = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii
  real(kind=8)    :: H_MEAN = 10.0_8 ** 4.0_8
  integer(kind=4) :: Nt     = 180 * 2, t
  real(kind=8)    :: T_max  = 5.0_8 * 3600.0_8 * 24.0_8, dt
  dt = T_max / Nt

  call domain%init(0.0_8, LX, 0, 128, 0.0_8, LY, 0, 128)
  call state%h%init_on_domain(domain)
  call state%u%init_on_domain(domain)
  call state%v%init_on_domain(domain)
  call op%init(central4, central4)

  call create_timescheme(timescheme, state, 'rk4')
  call swm_geostrophic_balance(state, domain, H_MEAN, 22.0_8 * 10.0_8 ** (-3.0_8), 1000000.0_8)

  do t = 0, Nt
    if (mod(t, 100) == 0) print *, 'step: ',  t
    call write_field(state%h, domain, './data/test2h.dat', t + 1)
    call timescheme%step(state, op, domain, dt)
  end do

end program test_2_geostriphic_balance
