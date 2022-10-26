program test_1_gaussian_hill
use initial_conditions_mod,            only : swm_gaussian_hill
use swe_advective_operator_mod,        only : swe_advective_operator_t
use horizontal_advection_operator_mod, only : horizontal_advection_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t
use central_differential_operator_mod, only : central2_t, central4_t
use timescheme_mod,                    only : timescheme_t
use timescheme_factory_mod,            only : create_timescheme
use rk4_mod,                           only : rk4_t
use explicit_Euler_mod,                only : explicit_Euler_t
use domain_mod,                        only : domain_t
use stvec_swe_mod,                     only : stvec_swe_t
use const_mod,                         only : Earth_radii, Earth_grav, pcori, pi
use read_write_mod,                    only : write_field
implicit none

type(domain_t)                 :: domain
type(stvec_swe_t)              :: state
type(horizontal_advection_operator_t) :: op
type(sbp21_t)                  :: sbp21
type(sbp42_t)                  :: sbp42
type(central2_t)               :: central2
type(central4_t)               :: central4
class(timescheme_t), allocatable :: timescheme

real(kind=8) :: LX = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii
real(kind=8) :: H_MEAN = 10.0_8 ** 4.0_8
integer :: Nt = 180, t
real(kind=8) :: T_max = 2.0_8 * 3600.0_8 * 24.0_8, dt
dt = T_max / Nt

call domain%init(0.0_8, LX, 0, 128, 0.0_8, LY, 0, 128)
call state%h%init_on_domain(domain)
call state%u%init_real(100.0_8, domain)
call state%v%init_on_domain(domain)
call op%init(sbp21, sbp21)

call create_timescheme(timescheme, state, 'rk4')
call swm_gaussian_hill(state, domain, H_MEAN, 50.0_8, 50.0_8)

do t = 0, Nt
  print *, t
  call write_field(state%h, domain, './data/test1.b', t + 1)
  call timescheme%step(state, op, domain, dt)
end do

end program test_1_gaussian_hill
