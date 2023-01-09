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

  subroutine set_swm_rotor_velocity(state, multi_domain, velocity_scale)
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

  end subroutine set_swm_rotor_velocity



  subroutine set_swm_gaussian_hill(state, multi_domain, h_mean, kx, ky, one_hill)
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

  end subroutine set_swm_gaussian_hill



  subroutine set_swm_geostrophic_balance(state, multi_domain, h_mean, scale_h, scale_sigma)
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

  end subroutine set_swm_geostrophic_balance



  subroutine set_swm_geostrophic_cyclone(state, multi_domain, h_mean, scale_h, scale_sigma, h0_const, u0_const, v0_const, field_type)
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

  end subroutine set_swm_geostrophic_cyclone

  subroutine set_swm_baratropic_instability(state, multi_domain, h_mean, u0)

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

            if (j <= multi_domain%subdomains(n, m)%je / 2) then
              do k = multi_domain%subdomains(n, m)%js, j
                state%h%subfields(n, m)%f(i, j) = state%h%subfields(n, m)%f(i, j) - pcori / Earth_grav * state%u%subfields(n, m)%f(i, k) * multi_domain%subdomains(n, m)%dy
              end do
            else
              state%h%subfields(n, m)%f(i, j) = state%h%subfields(n, m)%f(i, multi_domain%subdomains(n, m)%je - j)
            end if

            lx = abs(multi_domain%global_domain%xe - multi_domain%global_domain%xs)
            ly = abs(multi_domain%global_domain%ye - multi_domain%global_domain%ys)

            d1%subfields(n, m)%f(i, j) = ((multi_domain%subdomains(n, m)%x(i) - 0.85_8 * lx) / lx) ** 2.0_8 + ((multi_domain%subdomains(n, m)%y(j) - 0.75_8 * ly) / ly) ** 2.0_8
            d2%subfields(n, m)%f(i, j) = ((multi_domain%subdomains(n, m)%x(i) - 0.15_8 * lx) / lx) ** 2.0_8 + ((multi_domain%subdomains(n, m)%y(j) - 0.25_8 * ly) / ly) ** 2.0_8
            h_pret%subfields(n, m)%f(i, j) = 0.01_8 * h_mean * (exp(- c * d1%subfields(n, m)%f(i, j)) + exp(- c * d2%subfields(n, m)%f(i, j)))

            state%h%subfields(n, m)%f(i, j) = state%h%subfields(n, m)%f(i, j) + h_pret%subfields(n, m)%f(i, j)

          end do
        end do
      end do
    end do

  end subroutine set_swm_baratropic_instability

end module initial_conditions_mod
