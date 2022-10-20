program test_1_gaussian_hill
use initial_conditions_mod,            only : swm_gaussian_hill
use swe_advective_operator_mod,        only : swe_advective_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t
use central_differential_operator_mod, only : central2_t, central4_t
use domain_mod,                        only : domain_t
use stvec_swe_mod,                     only : stvec_swe_t
use const_mod,                         only : Earth_radii, Earth_grav, pcori, pi
implicit none

type(domain_t)                 :: domain
type(stvec_swe_t)              :: in, out
type(swe_advective_operator_t) :: op
type(sbp21_t)                  :: sbp21
type(sbp42_t)                  :: sbp42
type(central2_t)               :: central2
type(central4_t)               :: central4

real(kind=8) :: LX = 2 * pi * Earth_radii, LY = 2 * pi * Earth_radii
real(kind=8) :: H_MEAN = 10 ** 4
integer :: Nt = 180, t
real(kind=8) :: T_max = 5_8 * 3600_8 * 24.0_8, dt
dt = T_max / Nt

call domain%init(0.0_8, LX, 0, 128, 0.0_8, LY, 0, 128)
call in%h%init_on_domain(domain)
call in%u%init_on_domain(domain)
call in%v%init_on_domain(domain)
call out%h%init_on_domain(domain)
call out%u%init_on_domain(domain)
call out%v%init_on_domain(domain)
call op%init(sbp21, sbp21)

call swm_gaussian_hill(in, domain, H_MEAN, 50.0_8, 50.0_8)

do t = 0, Nt
  call op%apply(out, in, domain)
end do

print *, maxval(in%h%f)

end program test_1_gaussian_hill
