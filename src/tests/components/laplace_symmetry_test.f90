program laplace_symmetry_test
use field_mod,                     only : field_t
use domain_mod,                    only : domain_t
use multi_grid_field_mod,          only : multi_grid_field_t
use multi_domain_mod,              only : multi_domain_t
use sbp_differential_operator_mod, only : sbp21_2_t, sbp42_2_t
use SAT_mod,                       only : sbp_SAT_penalty_two_block_diffusion
use laplacian_mod,                 only : calc_laplacian
use norm_operators_mod,            only : weigh_with_Gram21_1d, weigh_with_Gram42_1d

implicit none

    type(multi_domain_t)         :: multi_domain
    type(domain_t)               :: global_domain
    type(multi_grid_field_t)     :: u, v, lap_u, lap_v, Hv
    type(sbp21_2_t)              :: sbp21
    type(sbp42_2_t)              :: sbp42
    integer(kind=4), allocatable :: deg(:, :)
    integer(kind=4)              :: n, m, i, j, numx = 1, numy = 1
    real(kind=8), allocatable    :: coefs(:, :)
    real(kind=8)                 :: dot1, dot2

    allocate(deg(1:numx, 1:numy))
    allocate(coefs(1:numx, 1:numy))
    deg(1, 1) = 1
    coefs(1, 1) = 1.0_8

    call global_domain%init(0.0_8, 10.0_8, 0, 10, 0.0_8, 10.0_8, 0, 10)
    call multi_domain%init(global_domain, numx, numy, deg)
    call u%init(multi_domain)
    call v%init(multi_domain)
    call Hv%init(multi_domain)
    call lap_u%init(multi_domain)
    call lap_v%init(multi_domain)

    sbp21%name = 'sbp21_2'
    sbp42%name = 'sbp42_2'

    do m = 1, multi_domain%num_sub_y
    do n = 1, multi_domain%num_sub_x
        do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
            call random_number(u%subfields(n, m)%f(i, j))
            call random_number(v%subfields(n, m)%f(i, j))
            Hv%subfields(n, m)%f(i, j) = v%subfields(n, m)%f(i, j)
        end do 
        end do
    end do
    end do

    call weigh_with_Gram42_1d(Hv, multi_domain, 'x') 
    call weigh_with_Gram42_1d(Hv, multi_domain, 'y') 
    call calc_laplacian(lap_u, u, multi_domain, coefs, sbp42) 
    call calc_laplacian(lap_v, v, multi_domain, coefs, sbp42)
    call weigh_with_Gram42_1d(lap_v, multi_domain, 'x') 
    call weigh_with_Gram42_1d(lap_v, multi_domain, 'y') 

    dot1 = 0.0_8
    dot2 = 0.0_8
    do m = 1, multi_domain%num_sub_y
    do n = 1, multi_domain%num_sub_x
        do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
        do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
            dot1 = dot1 + lap_u%subfields(n, m)%f(i, j) * Hv%subfields(n, m)%f(i, j)
            dot2 = dot2 + lap_v%subfields(n, m)%f(i, j) * u%subfields(n, m)%f(i, j)
        end do 
        end do
    end do
    end do

    print *, sqrt(dot1) - sqrt(dot2)


end program laplace_symmetry_test