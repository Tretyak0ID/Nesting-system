program test_2_geostriphic_balance
use initial_conditions_mod,            only : set_swm_geostrophic_balance
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use swe_advective_operator_mod,        only : swe_advective_operator_t
use horizontal_advection_operator_mod, only : horizontal_advection_operator_t
use diffusion_operator_mod,            only : diffusion_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t, sbp21_2_t, sbp42_2_t
use central_differential_operator_mod, only : central2_t, central4_t
use timescheme_mod,                    only : timescheme_t
use timescheme_factory_mod,            only : create_timescheme
use rk4_mod,                           only : rk4_t
use domain_mod,                        only : domain_t
use multi_domain_mod,                  only : multi_domain_t
use multi_grid_field_mod,              only : multi_grid_field_t
use stvec_swe_mod,                     only : stvec_swe_t
use const_mod,                         only : Earth_radii, Earth_grav, pcori, pi
use read_write_mod,                    only : write_field
use vec_math_mod,                      only : calc_sqrt_l2_norm_domain, calc_c_norm_domain, calc_l1_norm_domain
implicit none

  type(multi_domain_t)             :: multi_domain
  type(multi_grid_field_t)         :: buff1, buff2, buff3
  type(domain_t)                   :: domain
  type(stvec_swe_t)                :: state, diff_state
  type(swe_vect_inv_operator_t)    :: op
  type(sbp21_t)                    :: sbp21
  type(sbp42_t)                    :: sbp42
  type(diffusion_operator_t)       :: diffusion
  type(sbp21_2_t)                  :: sbp21_2
  type(sbp42_2_t)                  :: sbp42_2
  type(central2_t)                 :: central2
  type(central4_t)                 :: central4
  class(timescheme_t), allocatable :: timescheme, explicit_Euler
  integer(kind=4),     allocatable :: deg(:, :)
  real(kind=8),        allocatable :: coefs(:, :)
  real(kind=8) :: l2_norm, c_norm, l1_norm


  !test constants
  real(kind=8)    :: LX = 2.0_8 * pi * Earth_radii, LY = 2.0_8 * pi * Earth_radii, H_MEAN = 10.0_8 ** 4.0_8
  real(kind=8)    :: T_max  = 10.0_8 * 3600.0_8 * 24.0_8, dt, scale_h = 22.0e-3_8, scale_sigma = 1.0e6_8
  integer(kind=4) :: Nt = 180 * 16 * 2, Nx = 96 * 2, Ny = 96 * 2, num_sub_x = 3, num_sub_y = 3
  integer(kind=4) :: t, t_step_disp = 500, t_step_rec = 10, n, m
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

  !domain and dynamic operatior init
  sbp21%name = 'sbp21_1'
  sbp42%name = 'sbp42_1'
  central2%name = 'bp21_1'
  central4%name = 'bp42_1'
  sbp21_2%name  = 'sbp21_2'
  call domain%init(0.0_8, LX, 0, Nx, 0.0_8, LY, 0, Ny)
  call multi_domain%init(domain, num_sub_x, num_sub_y, deg)
  call state%h%init(multi_domain)
  call state%u%init(multi_domain)
  call state%v%init(multi_domain)
  call diff_state%h%init(multi_domain)
  call diff_state%u%init(multi_domain)
  call diff_state%v%init(multi_domain)
  call op%init(sbp42, sbp42, pcori, multi_domain)

  !diffusion operator init
  allocate(coefs(1:num_sub_x, 1:num_sub_y))
  do n = 1, num_sub_x
    do m = 1, num_sub_y
      coefs(n, m) = multi_domain%subdomains(n, m)%dx ** 2.0_8 / dt / 2.0_8
    end do
  end do
  call diffusion%init(sbp42_2, coefs, multi_domain)

  !time scheme init
  call create_timescheme(timescheme, state, 'rk4')
  call create_timescheme(explicit_Euler, state, 'explicit_Euler')

  !initial conditions
  call set_swm_geostrophic_balance(state, multi_domain, H_MEAN, scale_h, scale_sigma)

  do t = 0, Nt
    !step display
    ! if (mod(t, t_step_disp) == 0) print *, 'step: ',  t

    !recording
    ! if (num_sub_x > 1 .or. num_sub_y > 1) then 
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test23_11.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(1, 2), multi_domain%subdomains(1, 2), './data/test23_12.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(1, 3), multi_domain%subdomains(1, 3), './data/test23_13.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test23_21.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(2, 2), multi_domain%subdomains(2, 2), './data/test23_22.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(2, 3), multi_domain%subdomains(2, 3), './data/test23_23.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(3, 1), multi_domain%subdomains(3, 1), './data/test23_31.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(3, 2), multi_domain%subdomains(3, 2), './data/test23_32.dat', t / t_step_rec + 1)
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(3, 3), multi_domain%subdomains(3, 3), './data/test23_33.dat', t / t_step_rec + 1)
    ! else
    !   if (mod(t, t_step_rec) == 0) call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test2h.dat', t / t_step_rec + 1)
    ! end if

    !calculate
    call timescheme%step(state, op, multi_domain, dt)
    call explicit_Euler%step(state, diffusion, multi_domain, dt)
  end do

  call set_swm_geostrophic_balance(diff_state, multi_domain, H_MEAN, scale_h, scale_sigma)
  call diff_state%update(-1.0_8, state, multi_domain)
  l2_norm = calc_sqrt_l2_norm_domain(diff_state%h, multi_domain)
  c_norm = calc_c_norm_domain(diff_state%h, multi_domain)
  l1_norm = calc_l1_norm_domain(diff_state%h, multi_domain)
  print *, 'test_2_geostrophic_balance successfully completed'
  print *, 'l_2 norm: ', l2_norm
  print *, 'c norm: '  ,  c_norm
  print *, 'l_1 norm: ', l1_norm

end program test_2_geostriphic_balance
