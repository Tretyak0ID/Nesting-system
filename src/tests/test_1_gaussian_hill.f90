program test_1_gaussian_hill
use initial_conditions_mod,            only : swm_gaussian_hill
use swe_advective_operator_mod,        only : swe_advective_operator_t
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use horizontal_advection_operator_mod, only : horizontal_advection_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t
use timescheme_mod,                    only : timescheme_t
use timescheme_factory_mod,            only : create_timescheme
use rk4_mod,                           only : rk4_t
use explicit_Euler_mod,                only : explicit_Euler_t
use domain_mod,                        only : domain_t
use multi_domain_mod,                  only : multi_domain_t
use stvec_swe_mod,                     only : stvec_swe_t
use const_mod,                         only : Earth_radii, Earth_grav, pcori, pi
use read_write_mod,                    only : write_field
implicit none

type(domain_t)                   :: domain
type(multi_domain_t)             :: multi_domain
type(stvec_swe_t)                :: state
type(horizontal_advection_operator_t)   :: op
type(sbp21_t)                    :: sbp21
type(sbp42_t)                    :: sbp42
class(timescheme_t), allocatable :: timescheme
integer(kind=4), allocatable :: deg(:, :)

real(kind=8)    :: LX     = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii
real(kind=8)    :: H_MEAN = 10.0_8 ** 4.0_8
integer(kind=4) :: Nt     = 180 * 2, t
real(kind=8)    :: T_max  = 5.0_8 * 3600.0_8 * 24.0_8, dt
allocate(deg(1:2, 1:1))
deg(1, 1) = 1
deg(2, 1) = 1
dt = T_max / Nt

call domain%init(0.0_8, LX, 0, 128, 0.0_8, LY, 0, 128)
call multi_domain%init(domain, 2, 1, deg)
call state%h%init(multi_domain)
call state%u%init(multi_domain)
call state%v%init(multi_domain)
call op%init(sbp21, sbp21, multi_domain)

call create_timescheme(timescheme, state, 'rk4')
call swm_gaussian_hill(state, multi_domain, H_MEAN, 100.0_8, 100.0_8)

do t = 0, Nt
  print *, 'step: ',  t
  call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test1h.dat', t + 1)
  call timescheme%step(state, op, multi_domain, dt)
end do

end program test_1_gaussian_hill
