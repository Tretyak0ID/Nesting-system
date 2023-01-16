program diffusion_test
use initial_conditions_mod,            only : set_swm_gaussian_hill, set_swm_rotor_velocity
use diffusion_operator_mod,            only : diffusion_operator_t
use sbp_differential_operator_mod,     only : sbp21_2_t
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

type(domain_t)                        :: domain
type(multi_domain_t)                  :: multi_domain
type(stvec_swe_t)                     :: state
type(diffusion_operator_t)            :: op
type(sbp21_2_t)                       :: sbp21
class(timescheme_t), allocatable      :: timescheme
integer(kind=4),     allocatable      :: deg(:, :)
real(kind=8),        allocatable      :: coefs(:, :)

real(kind=8)    :: LX     = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii
real(kind=8)    :: H_MEAN = 10.0_8 ** 4.0_8
integer(kind=4) :: Nt     = 180 * 16, t, i, j, n, m
real(kind=8)    :: T_max  = 20.0_8 * 3600.0_8 * 24.0_8, dt

allocate(deg(1:2, 1:1))
allocate(coefs(1:2, 1:1))
deg(1, 1) = 2
deg(2, 1) = 1
dt = T_max / Nt

sbp21%name = 'sbp21_2'

call domain%init(0.0_8, LX, 0, 128, 0.0_8, LY, 0, 128)
call multi_domain%init(domain, 2, 1, deg)
call state%h%init(multi_domain)
call state%u%init(multi_domain)
call state%v%init(multi_domain)

coefs(1, 1) = 100000000.0_8
coefs(2, 1) = 100000000.0_8

call op%init(sbp21, coefs, multi_domain)

call create_timescheme(timescheme, state, 'rk4')
call set_swm_gaussian_hill(state, multi_domain, H_MEAN, 1000.0_8, 1000.0_8, 0)

do t = 0, Nt
  if (mod(t, 100) == 0) print *, 'step: ',  t
  if (mod(t, 20) == 0) call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/diffusion_left.dat', t/20 + 1)
  if (mod(t, 20) == 0) call write_field(state%h%subfields(2, 1), multi_domain%subdomains(2, 1), './data/diffusion_right.dat', t/20 + 1)
  call timescheme%step(state, op, multi_domain, dt)
end do

end program diffusion_test
