module initial_conditions_mod
  use field_mod,            only : field_t
  use domain_mod,           only : domain_t
  use multi_domain_mod,     only : multi_domain_t
  use multi_grid_field_mod, only : multi_grid_field_t
  use stvec_mod,            only : stvec_t
  use stvec_swe_mod,        only : stvec_swe_t
  use const_mod,            only : Earth_grav, Earth_radii, pcori, pi
implicit none

contains

  subroutine swm_rotor_velocity(state, multi_domain, velocity_scale)
  !Sets one elongated "Gaussian hill" in the middle of each subdivision subdomain,
  !used in conjunction with the net transfer operator
    type(stvec_swe_t),    intent(inout) :: state
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: velocity_scale

    real(kind=8)    :: dx, dy
    integer(kind=8) :: i, j, n, m

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
            state%u%subfields(n, m)%f(i, j) =  velocity_scale * sin((multi_domain%subdomains(n, m)%x(i)) / (2.0_8 * Earth_radii)) * cos((multi_domain%subdomains(n, m)%y(j)) / (2.0_8 * Earth_radii))
            state%v%subfields(n, m)%f(i, j) = -velocity_scale * cos((multi_domain%subdomains(n, m)%x(i)) / (2.0_8 * Earth_radii)) * sin((multi_domain%subdomains(n, m)%y(j)) / (2.0_8 * Earth_radii))
          end do
        end do
      end do
    end do

  end subroutine swm_rotor_velocity



  subroutine swm_gaussian_hill(state, multi_domain, h_mean, kx, ky, one_hill)
  !Sets the Gaussian slide, the scale is adjusted by the parameters kx and ky
    type(stvec_swe_t),    intent(inout) :: state
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean, kx, ky
    integer(kind=2),      intent(in)    :: one_hill

    real(kind=8)    :: dx, dy
    integer(kind=8) :: i, j, n, m

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
            if (one_hill.eq.1) then
              dx = (multi_domain%subdomains(n, m)%x(i) - 0.8_8 * (maxval(multi_domain%global_domain%x) - minval(multi_domain%global_domain%x)) / 2.0_8) / multi_domain%global_domain%x(multi_domain%global_domain%ie)
              dy = (multi_domain%subdomains(n, m)%y(j) - 0.8_8 * (maxval(multi_domain%global_domain%y) - minval(multi_domain%global_domain%y)) / 2.0_8) / multi_domain%global_domain%y(multi_domain%global_domain%je)
              state%h%subfields(n, m)%f(i, j) = h_mean + 0.1_8 * h_mean * exp(- kx * (dx ** 2.0_8) - ky * (dy ** 2.0_8))
            else
              dx = (multi_domain%subdomains(n, m)%x(i) - (maxval(multi_domain%subdomains(n, m)%x) + minval(multi_domain%subdomains(n, m)%x)) / 2.0_8) / multi_domain%global_domain%x(multi_domain%global_domain%ie)
              dy = (multi_domain%subdomains(n, m)%y(j) - (maxval(multi_domain%subdomains(n, m)%y) + minval(multi_domain%subdomains(n, m)%y)) / 2.0_8) / multi_domain%global_domain%y(multi_domain%global_domain%je)
              state%h%subfields(n, m)%f(i, j) = h_mean + 0.1_8 * h_mean * exp(- kx * (dx ** 2.0_8) - ky * (dy ** 2.0_8))
            end if
          end do
        end do
      end do
    end do

  end subroutine swm_gaussian_hill



  subroutine swm_geostrophic_balance(state, multi_domain, h_mean, scale_h, scale_sigma)
  !Geostrophic balanced vortex of characteristic depth scale_h
  !and characteristic scale scale_sigma.
    type(stvec_swe_t),    intent(inout) :: state
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean, scale_h, scale_sigma

    type(multi_grid_field_t) :: dhdr, vtan
    real(kind=8)             :: dx, dy
    integer(kind=8)          :: i, j, n, m

    call dhdr%init(multi_domain)
    call vtan%init(multi_domain)

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
            dx = (multi_domain%subdomains(n, m)%x(i) - 0.98_8 * (maxval(multi_domain%global_domain%x) - minval(multi_domain%global_domain%x)) / 2.0_8)
            dy = (multi_domain%subdomains(n, m)%y(j) - (maxval(multi_domain%global_domain%y) - minval(multi_domain%global_domain%y)) / 2.0_8)

            dhdr%subfields(n, m)%f(i, j) = h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8) * scale_h * 2.0_8 * sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma ** 2.0_8
            vtan%subfields(n, m)%f(i, j) = (- sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori + sqrt((sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori) ** 2.0_8 + 4.0_8 * Earth_grav * sqrt(dx ** 2.0_8 + dy ** 2.0_8) * dhdr%subfields(n, m)%f(i, j))) / 2.0_8

            state%h%subfields(n, m)%f(i, j) = h_mean - (h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8)) * scale_h
            state%u%subfields(n, m)%f(i, j) = - vtan%subfields(n, m)%f(i, j) * cos(pi / 2.0_8 - atan2(dy, dx))
            state%v%subfields(n, m)%f(i, j) =   vtan%subfields(n, m)%f(i, j) * sin(pi / 2.0_8 - atan2(dy, dx))
          end do
        end do
      end do
    end do

  end subroutine swm_geostrophic_balance



  subroutine swm_geostrophic_cyclone(state, multi_domain, h_mean, scale_h, scale_sigma, h0_const, u0_const, v0_const, field_type)
  !Geostrophic balanced cyclone in a constant velocity field
  !of characteristic depth scale_h and characteristic scale scale_sigma.
    type(stvec_swe_t),    intent(inout) :: state
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean, scale_h, scale_sigma, u0_const, v0_const, h0_const
    character(len=*),     intent(in)    :: field_type

    type(multi_grid_field_t) :: dhdr, vtan, h0, u0, v0
    real(kind=8)             :: dx, dy
    integer(kind=8)          :: i, j, n, m

    call dhdr%init(multi_domain)
    call vtan%init(multi_domain)
    call h0%init(multi_domain)
    call u0%init(multi_domain)
    call v0%init(multi_domain)

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
            dx = (multi_domain%subdomains(n, m)%x(i) - 0.6_8 * (maxval(multi_domain%global_domain%x) - minval(multi_domain%global_domain%x)) / 2.0_8)
            dy = (multi_domain%subdomains(n, m)%y(j) - (maxval(multi_domain%global_domain%y) - minval(multi_domain%global_domain%y)) / 2.0_8)

            dhdr%subfields(n, m)%f(i, j) = h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8) * scale_h * 2.0_8 * sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma ** 2.0_8
            vtan%subfields(n, m)%f(i, j) = (- sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori + sqrt((sqrt(dx ** 2.0_8 + dy ** 2.0_8) * pcori) ** 2.0_8 + 4.0_8 * Earth_grav * sqrt(dx ** 2.0_8 + dy ** 2.0_8) * dhdr%subfields(n, m)%f(i, j))) / 2.0_8

            if (field_type == 'const_periodic') then
              h0%subfields(n, m)%f(i, j) = pcori  / Earth_grav * cos((multi_domain%subdomains(n, m)%y(j) - 10.0_8 ** 7.0_8) / (6371.22_8 * 1000.0_8)) * u0_const * 6371.22_8 * 1000.0_8
              u0%subfields(n, m)%f(i, j) = sin((multi_domain%subdomains(n, m)%y(j) - 10.0_8 ** 7.0_8) / (6371.22_8 * 1000.0_8)) * u0_const
              v0%subfields(n, m)%f(i, j) = v0_const
            else if (field_type == 'const_discontinuum') then
              h0%subfields(n, m)%f(i, j) = h0_const
              u0%subfields(n, m)%f(i, j) = u0_const
              v0%subfields(n, m)%f(i, j) = v0_const
            else
              print *, 'Undefine field type'
            end if

            state%h%subfields(n, m)%f(i, j) = h_mean - (h_mean * exp( - (sqrt(dx ** 2.0_8 + dy ** 2.0_8) / scale_sigma) ** 2.0_8)) * scale_h + h0%subfields(n, m)%f(i, j)
            state%u%subfields(n, m)%f(i, j) = - vtan%subfields(n, m)%f(i, j) * cos(pi / 2.0_8 - atan2(dy, dx)) + u0%subfields(n, m)%f(i, j)
            state%v%subfields(n, m)%f(i, j) =   vtan%subfields(n, m)%f(i, j) * sin(pi / 2.0_8 - atan2(dy, dx)) + v0%subfields(n, m)%f(i, j)
          end do
        end do
      end do
    end do

  end subroutine swm_geostrophic_cyclone

  subroutine swm_baratropic_instability(state, multi_domain, h_mean, u0)

    type(stvec_swe_t),    intent(inout) :: state
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean, u0

    integer(kind=8) :: i, j, n, m, k
    real(kind=8)    :: c, lx, ly
    type(multi_grid_field_t) :: d1, d2, h_pret

    c = 1000.0_8

    call state%h%create_similar(h_pret)
    call state%h%create_similar(d1)
    call state%h%create_similar(d2)

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je

            state%u%subfields(n, m)%f(i, j) = u0 * (sin(2.0_8 * pi * multi_domain%subdomains(n, m)%y(j) / abs(multi_domain%global_domain%ye - multi_domain%global_domain%ys)) ** 81.0_8)
            state%h%subfields(n, m)%f(i, j) = h_mean !phone field for h

            do k = multi_domain%subdomains(n, m)%js, j
              state%h%subfields(n, m)%f(i, j) = state%h%subfields(n, m)%f(i, j) - pcori / Earth_grav * state%u%subfields(n, m)%f(i, k) * multi_domain%subdomains(n, m)%dy
            end do

            lx = abs(multi_domain%global_domain%xe - multi_domain%global_domain%xs)
            ly = abs(multi_domain%global_domain%ye - multi_domain%global_domain%ys)

            d1%subfields(n, m)%f(i, j) = ((multi_domain%subdomains(n, m)%x(i) - 0.85_8 * lx) / lx) ** 2.0_8 + ((multi_domain%subdomains(n, m)%y(j) - 0.75_8 * ly) / ly) ** 2.0_8
            d2%subfields(n, m)%f(i, j) = ((multi_domain%subdomains(n, m)%x(i) - 0.15_8 * lx) / lx) ** 2.0_8 + ((multi_domain%subdomains(n, m)%y(j) - 0.25_8 * ly) / ly) ** 2.0_8
            h_pret%subfields(n, m)%f(i, j) = 0.0_8!0.01_8 * h_mean * (exp(- c * d1%subfields(n, m)%f(i, j)) + exp(- c * d2%subfields(n, m)%f(i, j)))

            state%h%subfields(n, m)%f(i, j) = state%h%subfields(n, m)%f(i, j) + h_pret%subfields(n, m)%f(i, j)

          end do
        end do
      end do
    end do

  end subroutine swm_baratropic_instability

  subroutine set_two_cyclones_initial_conditions(state, multi_domain, h_mean)

    type(stvec_swe_t),    intent(inout) :: state
    type(multi_domain_t), intent(in)    :: multi_domain
    real(kind=8),         intent(in)    :: h_mean

    real(kind=8), parameter :: c1 = 0.3398_8, c2 = 5.377e-4
    real(kind=8), parameter :: v0 = 71.521_8, f0 = 0.62e-4
    real(kind=8), parameter :: r0 = 100e3_8

    integer(kind=8) :: i, j, n, m, k
    real(kind=8)    :: cx1, cy1, s1, ex1, ey1, v1
    real(kind=8)    :: cx2, cy2, s2, ex2, ey2, v2
    type(multi_grid_field_t) :: d1, d2, h_pret

    cy1 = 0.5_8*(multi_domain%global_domain%ye + multi_domain%global_domain%ys)
    cx1 = 0.5_8*(multi_domain%global_domain%xe + multi_domain%global_domain%xs)+2e5_8
    cy2 = cy1
    cx2 = cx1 - 4e5_8

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
          do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je

              s1 = sqrt((multi_domain%subdomains(n, m)%x(i)-cx1)**2+ &
                        (multi_domain%subdomains(n, m)%y(j)-cy1)**2)
              ex1 = (multi_domain%subdomains(n, m)%x(i)-cx1) / max(s1,1e-14)
              ey1 = (multi_domain%subdomains(n, m)%y(j)-cy1) / max(s1,1e-14)
              s1 = s1 / r0

              s2 = sqrt((multi_domain%subdomains(n, m)%x(i)-cx2)**2+ &
                        (multi_domain%subdomains(n, m)%y(j)-cy2)**2)
              ex2 = (multi_domain%subdomains(n, m)%x(i)-cx2) / max(s2,1e-14)
              ey2 = (multi_domain%subdomains(n, m)%y(j)-cy2) / max(s2,1e-14)
              s2 = s2 / r0

              v1 = cyclone_v(s1); v2 = cyclone_v(s2)

              state%u%subfields(n, m)%f(i, j) =-v1*ey1-v2*ey2
              state%v%subfields(n, m)%f(i, j) = v1*ex1+v2*ex2
              state%h%subfields(n, m)%f(i, j) = h_mean+integrate_h(s1)+integrate_h(s2)

          end do
        end do
      end do
    end do

    contains
        real(kind=8) function cyclone_v(s) result(v)
            real(kind=8), intent(in) :: s
            v = v0*s*(1.0_8+3*c2/c1*s**4) / (1.0_8+c1*s**2+c2*s**6)**2
        end function
        real(kind=8) function cyclone_dghds(s) result(dhds)
            real(kind=8), intent(in) :: s
            real(kind=8) :: vt
            vt = cyclone_v(s)
            dhds = r0*f0*vt+vt**2/s
        end function
        real(kind=8) function integrate_h(s) result(h)
            real(kind=8) :: s
            real(kind=8), parameter :: dgh0 = 3228.906700369434_8
            real(kind=8), parameter :: x1 = 1.0_8/3.0_8*sqrt(5.0_8+2.0_8*sqrt(10.0_8/7.0_8))
            real(kind=8), parameter :: x2 = 1.0_8/3.0_8*sqrt(5.0_8-2.0_8*sqrt(10.0_8/7.0_8))
            real(kind=8), parameter :: x5(5) = [-x1,-x2,0.0_8,x2,x1]
            real(kind=8), parameter :: w1 = (322.0_8-13.0_8*sqrt(70.0_8))/900.0_8
            real(kind=8), parameter :: w2 = (322.0_8+13.0_8*sqrt(70.0_8))/900.0_8
            real(kind=8), parameter :: w5(5) = [w1,w2,128.0_8/225.0_8,w2,w1]
            real(kind=8) :: a, b, x
            integer(kind=4) :: iq, iseg

            h = -dgh0
            do iseg= 1, floor(s)+1
                a = (iseg-1.0_8); b = min(1.0_8*iseg,s)
                do iq = 1, 5
                    x = 0.5_8*(b-a)*x5(iq)+0.5_8*(a+b)
                    h = h+w5(iq)*cyclone_dghds(x)*0.5_8*(b-a)
                end do
            end do
            h = h / Earth_grav
        end function
  end subroutine set_two_cyclones_initial_conditions

end module initial_conditions_mod
