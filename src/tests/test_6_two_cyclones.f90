program test_6_two_cyclones
use initial_conditions_mod,            only : set_two_cyclones_initial_conditions
use swe_vect_inv_operator_mod,         only : swe_vect_inv_operator_t
use diffusion_operator_mod,            only : diffusion_operator_t
use sbp_differential_operator_mod,     only : sbp21_t, sbp42_t, sbp21_2_t, sbp42_2_t
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
  type(diffusion_operator_t)       :: diffusion
  type(sbp21_t)                    :: sbp21
  type(sbp42_t)                    :: sbp42
  type(sbp21_2_t)                  :: sbp21_2
  type(sbp42_2_t)                  :: sbp42_2
  type(central2_t)                 :: central2
  type(central4_t)                 :: central4
  class(timescheme_t), allocatable :: timescheme, explicit_Euler
  integer(kind=4),     allocatable :: deg(:, :)
  real(kind=8),        allocatable :: coefs(:, :)

  real(kind=8),    parameter :: LX = 4000e3_8, LY = 4000e3_8
  integer(kind=4), parameter :: Nx = 96, Ny = 96
  integer(kind=8), parameter :: num_sub_x = 3, num_sub_y = 3
  real(kind=8),    parameter :: H_MEAN = 1e4_8
  real(kind=8),    parameter :: T_max  = 256.0_8*3600.0_8, dt = 25.0_8, tau_wr = 3600.0_8
  integer(kind=4), parameter :: Nt     = int(T_max/dt)
  integer(kind=4)            :: t, nzap, irec, n, m

  nzap = int(tau_wr/dt)

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

  sbp21%name    = 'sbp21_1'
  sbp42%name    = 'sbp42_1'
  central2%name = 'cent2_1'
  central4%name = 'cent4_1'
  sbp21_2%name  = 'sbp21_2'
  sbp42_2%name  = 'sbp42_2'

  call domain%init(0.0_8, LX, 0, Nx, pcori, LY, 0, Ny)
  call multi_domain%init(domain, num_sub_x, num_sub_y, deg)
  call curl%init(multi_domain)
  call state%h%init(multi_domain)
  call state%u%init(multi_domain)
  call state%v%init(multi_domain)
  call op%init(sbp42, sbp42, 0.62_8 * 10e-4, multi_domain)

  allocate(coefs(1:num_sub_x, 1:num_sub_y))
  do n = 1, num_sub_x
    do m = 1, num_sub_y
      coefs(n, m) = multi_domain%subdomains(n, m)%dx ** 2.0_8 / sqrt(dt) / 100.0_8
    end do
  end do
  call diffusion%init(sbp42_2, coefs, multi_domain)

  call create_timescheme(timescheme, state, 'rk4')
  call create_timescheme(explicit_Euler, state, 'rk4')
  call set_two_cyclones_initial_conditions(state, multi_domain, H_MEAN)

  do t = 0, Nt
    if(mod(t,nzap)==0) then
        irec = t/nzap+1
        print *, "write t=", t,irec
        !call write_field(state%u%subfields(1, 1), multi_domain%subdomains(1, 1), './data/two_cyclones_u.dat', irec)
        !call write_field(state%v%subfields(1, 1), multi_domain%subdomains(1, 1), './data/two_cyclones_v.dat', irec)
        !call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/two_cyclones_h.dat', irec)
        call calc_curl(curl, state%u, state%v, multi_domain, sbp42, sbp42)
        if (num_sub_x > 1 .or. num_sub_y > 1) then
          call write_field(curl%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test6_96_11_curl.dat', irec)
          call write_field(curl%subfields(1, 2), multi_domain%subdomains(1, 2), './data/test6_96_12_curl.dat', irec)
          call write_field(curl%subfields(1, 3), multi_domain%subdomains(1, 3), './data/test6_96_13_curl.dat', irec)
          call write_field(curl%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test6_96_21_curl.dat', irec)
          call write_field(curl%subfields(2, 2), multi_domain%subdomains(2, 2), './data/test6_96_22_curl.dat', irec)
          call write_field(curl%subfields(2, 3), multi_domain%subdomains(2, 3), './data/test6_96_23_curl.dat', irec)
          call write_field(curl%subfields(3, 1), multi_domain%subdomains(3, 1), './data/test6_96_31_curl.dat', irec)
          call write_field(curl%subfields(3, 2), multi_domain%subdomains(3, 2), './data/test6_96_32_curl.dat', irec)
          call write_field(curl%subfields(3, 3), multi_domain%subdomains(3, 3), './data/test6_96_33_curl.dat', irec)

          call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test6_96_11_h.dat', irec)
          call write_field(state%h%subfields(1, 2), multi_domain%subdomains(1, 2), './data/test6_96_12_h.dat', irec)
          call write_field(state%h%subfields(1, 3), multi_domain%subdomains(1, 3), './data/test6_96_13_h.dat', irec)
          call write_field(state%h%subfields(2, 1), multi_domain%subdomains(2, 1), './data/test6_96_21_h.dat', irec)
          call write_field(state%h%subfields(2, 2), multi_domain%subdomains(2, 2), './data/test6_96_22_h.dat', irec)
          call write_field(state%h%subfields(2, 3), multi_domain%subdomains(2, 3), './data/test6_96_23_h.dat', irec)
          call write_field(state%h%subfields(3, 1), multi_domain%subdomains(3, 1), './data/test6_96_31_h.dat', irec)
          call write_field(state%h%subfields(3, 2), multi_domain%subdomains(3, 2), './data/test6_96_32_h.dat', irec)
          call write_field(state%h%subfields(3, 3), multi_domain%subdomains(3, 3), './data/test6_96_33_h.dat', irec)
        else
          call write_field(curl%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test6_96_curl.dat', irec)
          call write_field(state%h%subfields(1, 1), multi_domain%subdomains(1, 1), './data/test6_96_h.dat', irec)
        end if
    end if
    call timescheme%step(state, op, multi_domain, dt)
    call explicit_Euler%step(state, diffusion, multi_domain, dt)
  end do

end program test_6_two_cyclones
