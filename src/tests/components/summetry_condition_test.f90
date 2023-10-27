program summetry_condition_test
use initial_conditions_mod,            only : set_swm_gaussian_hill, set_swm_rotor_velocity
use diffusion_operator_mod,            only : diffusion_operator_t
use sbp_differential_operator_mod,     only : sbp21_2_t, sbp42_2_t
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
type(sbp42_2_t)                       :: sbp42
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

sbp42%name = 'sbp42_2'

call domain%init(0.0_8, LX, 0, 128, 0.0_8, LY, 0, 128)
call multi_domain%init(domain, 2, 1, deg)
call state%h%init(multi_domain)
call state%u%init(multi_domain)
call state%v%init(multi_domain)

do n = 1, 2
    do m = 1, 1
      coefs(n, m) = multi_domain%subdomains(n, m)%dx ** 2.0_8 / sqrt(dt) / 50.0_8
    end do
 end do



end program summetry_condition_test