program zero_vorticity_gradient_test
use initial_conditions_mod,            only : set_swm_gaussian_hill, set_swm_rotor_velocity
use swe_advective_operator_mod,        only : swe_advective_operator_t
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use horizontal_advection_operator_mod, only : horizontal_advection_operator_t
use central_differential_operator_mod, only : central2_t, central4_t
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
use vec_math_mod,                      only : calc_sqrt_l2_norm_domain, calc_c_norm_domain, calc_l1_norm_domain
use grad_mod,                          only : calc_grad
use curl_mod,                          only : calc_curl
use multi_grid_field_mod,              only : multi_grid_field_t
implicit none

  type(domain_t)                        :: domain
  type(multi_domain_t)                  :: multi_domain
  type(stvec_swe_t)                     :: state
  type(sbp21_t)                         :: sbp21
  type(sbp42_t)                         :: sbp42
  type(central2_t)                      :: central2
  type(central4_t)                      :: central4
  integer(kind=4),     allocatable      :: deg(:, :)
  real(kind=8)                          :: l2, c, l1
  type(multi_grid_field_t)              :: gx, gy, curl

  !test constants
  real(kind=8)    :: LX = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii, H_MEAN = 10.0_8 ** 4.0_8
  real(kind=8)    :: T_max  = 20.0_8 * 3600.0_8 * 24.0_8, dt, max_v = 100.0_8, Kx = 1e3_8, Ky = 5e1_8
  integer(kind=4) :: Nt = 180 * 16 + 2, Nx = 128, Ny = 128, num_sub_x = 2, num_sub_y = 1
  integer(kind=4) :: t, t_step_disp = 500, t_step_rec = 1
  dt = T_max / Nt

  allocate(deg(1:num_sub_x, 1:num_sub_y))
  deg(1, 1) = 1
  if (num_sub_x > 1) then
    deg(2, 1) = 2
  end if

  !domain and dynamic operator init
  sbp21%name = 'sbp21_1'
  sbp42%name = 'sbp42_1'
  central2%name = 'cent2_1'
  central4%name = 'cent4_1'
  call domain%init(0.0_8, LX, 0, Nx, 0.0_8, LY, 0, Ny)
  call multi_domain%init(domain, num_sub_x, num_sub_y, deg)
  call state%h%init(multi_domain)
  call state%u%init(multi_domain)
  call state%v%init(multi_domain)
  call curl%init(multi_domain)
  call gx%init(multi_domain)
  call gy%init(multi_domain)

  !set initial conditions
  call set_swm_gaussian_hill(state, multi_domain, H_MEAN, Kx, Ky, 1)
  call set_swm_rotor_velocity(state, multi_domain, max_v)

  do t = 0, Nt
    !step display
    if (mod(t, t_step_disp) == 0) print *, 'step: ',  t

    !recording
    call calc_grad(gx, gy, state%h, multi_domain, sbp42, sbp42)
    call calc_curl(curl, gx, gy, multi_domain, sbp42, sbp42)
    l2 = calc_sqrt_l2_norm_domain(curl, multi_domain)
    c  = calc_c_norm_domain(curl, multi_domain)

    print *, 'test_0_horizontal_advection successfully completed'
    print *, 'l2-norm: ', l2
    print *, 'c-norm: ', c

    !calculate
    
  end do

end program zero_vorticity_gradient_test