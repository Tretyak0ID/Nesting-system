program test_3_geostriphic_cyclone_discontinuum
use initial_conditions_mod,            only : set_swm_geostrophic_cyclone
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

  type(multi_domain_t)             :: multi_domain, multi_domain_left, multi_domain_right
  type(domain_t)                   :: domain, domain_left, domain_right
  type(multi_grid_field_t)         :: curl, curl_left, curl_right
  type(stvec_swe_t)                :: state, state_left, state_right
  type(swe_vect_inv_operator_t)    :: op, op_left, op_right
  type(diffusion_operator_t)       :: diffusion, diffusion_left, diffusion_right
  type(sbp21_2_t)                  :: sbp21_2
  type(sbp21_t)                    :: sbp21
  type(sbp42_t)                    :: sbp42
  type(central2_t)                 :: central2
  type(central4_t)                 :: central4
  class(timescheme_t), allocatable :: timescheme, explicit_Euler, timescheme_left, timescheme_right, explicit_Euler_left, explicit_Euler_right
  integer(kind=4),     allocatable :: deg(:, :)
  real(kind=8),        allocatable :: coefs(:, :)

  !---------------------multi-block-scale---------------------

  !test constants
  real(kind=8)    :: LX = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii, H_MEAN = 10.0_8 ** 4.0_8
  real(kind=8)    :: T_max  = 30.0_8 * 3600.0_8 * 24.0_8, dt, scale_h = 22.0e-3_8, scale_sigma = 1e6_8, h0 = 0.0_8, u0 = 10.0_8, v0 = 0.0_8 
  integer(kind=4) :: Nt = 180 * 32, Nx = 192, Ny = 192, num_sub_x = 3, num_sub_y = 3
  integer(kind=4) :: t, n, m, t_step_disp = 500, t_step_rec = 50
  dt = T_max / Nt

  allocate(deg(1:num_sub_x, 1:num_sub_y))
  deg(1, 1) = 1
  if (num_sub_x > 1 .or. num_sub_y > 1) then
    deg(1, 2) = 1
    deg(1, 3) = 1
    deg(2, 1) = 1
    deg(2, 2) = 2
    deg(2, 3) = 1
    deg(3, 1) = 1
    deg(3, 2) = 1
    deg(3, 3) = 1
  end if

  !domain and dynamic operator init
  sbp21%name    = 'sbp21_1'
  sbp42%name    = 'sbp42_1'
  central2%name = 'cent2_1'
  central4%name = 'cent4_1'
  sbp21_2%name  = 'sbp21_2'
  call domain%init(0.0_8, LX, 0, Nx, 0.0_8, LY, 0, Ny)
  call multi_domain%init(domain, num_sub_x, num_sub_y, deg)
  call curl%init(multi_domain)
  call state%h%init(multi_domain)
  call state%u%init(multi_domain)
  call state%v%init(multi_domain)
  call op%init(sbp42, sbp42, pcori, multi_domain)

  !diffusion operator init
  allocate(coefs(1:num_sub_x, 1:num_sub_y))
  do n = 1, num_sub_x
    do m = 1, num_sub_y
      coefs(n, m) = multi_domain%subdomains(n, m)%dx ** 2.0_8 / dt / 2.0_8
    end do
  end do
  call diffusion%init(sbp21_2, coefs, multi_domain)

  !time scheme init
  call create_timescheme(timescheme, state, 'rk4')
  call create_timescheme(explicit_Euler, state, 'explicit_Euler')

  !initial conditions init
  call set_swm_geostrophic_cyclone(state, multi_domain, H_MEAN, scale_h, scale_sigma, h0, u0, v0, 'const_discontinuum')

  do t = 0, Nt
    !step display
    if (mod(t, t_step_disp) == 0) print *, 'step: ',  t

    !recording
    if (mod(t, t_step_rec) == 0) call calc_curl(curl, state%u, state%v, multi_domain, sbp42, sbp42)
    if (num_sub_x > 1 .or. num_sub_y > 1) then 
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test42_11.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(1, 2), multi_domain%subdomains(1, 2), './data/test42_12.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(1, 3), multi_domain%subdomains(1, 3), './data/test42_13.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test42_21.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(2, 2), multi_domain%subdomains(2, 2), './data/test42_22.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(2, 3), multi_domain%subdomains(2, 3), './data/test42_23.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(3, 1), multi_domain%subdomains(3, 1), './data/test42_31.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(3, 2), multi_domain%subdomains(3, 2), './data/test42_32.dat', t / t_step_rec + 1)
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(3, 3), multi_domain%subdomains(3, 3), './data/test42_33.dat', t / t_step_rec + 1)
    else
      if (mod(t, t_step_rec) == 0) call write_field(curl%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test4curl.dat', t / t_step_rec + 1)
    end if

    !calculate
    call timescheme%step(state, op, multi_domain, dt)
    call explicit_Euler%step(state, diffusion, multi_domain, dt)
  end do

  print *, 'test_4_geostrophic_cyclone_discontinuum successfully completed'

end program test_3_geostriphic_cyclone_discontinuum
