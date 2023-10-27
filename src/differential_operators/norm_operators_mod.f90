module norm_operators_mod
    use multi_domain_mod, only     : multi_domain_t
    use multi_grid_field_mod, only : multi_grid_field_t
implicit none

contains
    subroutine weigh_with_Gram21_1d(in, domains, direction)
        type(multi_grid_field_t), intent(inout) :: in
        type(multi_domain_t),     intent(in)    :: domains
        character(len=1),         intent(in)    :: direction
        integer(kind = 8) :: n, m, i, j
        
        if (direction == 'x') then
        do m = 1, domains%num_sub_y
        do n = 1, domains%num_sub_x
            do j = domains%subdomains(n, m)%js, domains%subdomains(n, m)%je
                in%subfields(n, m)%f(domains%subdomains(n, m)%is, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%is, j) * domains%subdomains(n, m)%dx / 2.0_8
                do i = domains%subdomains(n, m)%is + 1, domains%subdomains(n, m)%ie - 1
                    in%subfields(n, m)%f(i, j) = in%subfields(n, m)%f(i, j) * domains%subdomains(n, m)%dx
                end do
                in%subfields(n, m)%f(domains%subdomains(n, m)%ie, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%ie, j) * domains%subdomains(n, m)%dx / 2.0_8
            end do
        end do
        end do

        else if (direction == 'y') then
        do m = 1, domains%num_sub_y
        do n = 1, domains%num_sub_x
            do i = domains%subdomains(n, m)%is, domains%subdomains(n, m)%ie
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%js) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%js) * domains%subdomains(n, m)%dy / 2.0_8
                do j = domains%subdomains(n, m)%js + 1, domains%subdomains(n, m)%je - 1
                    in%subfields(n, m)%f(i, j) = in%subfields(n, m)%f(i, j) * domains%subdomains(n, m)%dy
                end do
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%je) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%je) * domains%subdomains(n, m)%dy / 2.0_8
            end do
        end do
        end do
        end if
    end subroutine weigh_with_Gram21_1d

    subroutine weigh_with_Gram42_1d(in, domains, direction)
        type(multi_grid_field_t), intent(inout) :: in
        type(multi_domain_t),     intent(in)    :: domains
        character(len=1),         intent(in)    :: direction
        integer(kind=8) :: n, m, i, j
        
        if (direction == 'x') then
        do m = 1, domains%num_sub_y
        do n = 1, domains%num_sub_x
            do j = domains%subdomains(n, m)%js, domains%subdomains(n, m)%je
                in%subfields(n, m)%f(domains%subdomains(n, m)%is, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%is, j) * domains%subdomains(n, m)%dy * 17.0_8 / 48.0_8
                in%subfields(n, m)%f(domains%subdomains(n, m)%is + 1, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%is + 1, j) * domains%subdomains(n, m)%dy * 59.0_8 / 48.0_8
                in%subfields(n, m)%f(domains%subdomains(n, m)%is + 2, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%is + 2, j) * domains%subdomains(n, m)%dy * 43.0_8 / 48.0_8
                in%subfields(n, m)%f(domains%subdomains(n, m)%is + 3, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%is + 3, j) * domains%subdomains(n, m)%dy * 49.0_8 / 48.0_8
                do i = domains%subdomains(n, m)%is + 4, domains%subdomains(n, m)%ie - 4
                    in%subfields(n, m)%f(i, j) = in%subfields(n, m)%f(i, j) * domains%subdomains(n, m)%dx
                end do
                in%subfields(n, m)%f(domains%subdomains(n, m)%ie, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%ie, j) * domains%subdomains(n, m)%dy * 17.0_8 / 48.0_8
                in%subfields(n, m)%f(domains%subdomains(n, m)%ie - 1, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%ie - 1, j) * domains%subdomains(n, m)%dy * 59.0_8 / 48.0_8
                in%subfields(n, m)%f(domains%subdomains(n, m)%ie - 2, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%ie - 2, j) * domains%subdomains(n, m)%dy * 43.0_8 / 48.0_8
                in%subfields(n, m)%f(domains%subdomains(n, m)%ie - 3, j) = in%subfields(n, m)%f(domains%subdomains(n, m)%ie - 3, j) * domains%subdomains(n, m)%dy * 49.0_8 / 48.0_8
            end do
        end do
        end do

        else if (direction == 'y') then
        do m = 1, domains%num_sub_y
        do n = 1, domains%num_sub_x
            do i = domains%subdomains(n, m)%is, domains%subdomains(n, m)%ie
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%js) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%js) * domains%subdomains(n, m)%dy * 17.0_8 / 48.0_8
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%js + 1) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%js + 1) * domains%subdomains(n, m)%dy * 59.0_8 / 48.0_8
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%js + 2) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%js + 2) * domains%subdomains(n, m)%dy * 43.0_8 / 48.0_8
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%js + 3) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%js + 3) * domains%subdomains(n, m)%dy * 49.0_8 / 48.0_8
                do j = domains%subdomains(n, m)%js + 4, domains%subdomains(n, m)%je - 4
                    in%subfields(n, m)%f(i, j) = in%subfields(n, m)%f(i, j) * domains%subdomains(n, m)%dy
                end do
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%je) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%je) * domains%subdomains(n, m)%dy * 17.0_8 / 48.0_8
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%je - 1) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%je - 1) * domains%subdomains(n, m)%dy * 59.0_8 / 48.0_8
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%je - 2) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%je - 2) * domains%subdomains(n, m)%dy * 43.0_8 / 48.0_8
                in%subfields(n, m)%f(i, domains%subdomains(n, m)%je - 3) = in%subfields(n, m)%f(i, domains%subdomains(n, m)%je - 3) * domains%subdomains(n, m)%dy * 49.0_8 / 48.0_8
            end do
        end do
        end do
        end if
    end subroutine weigh_with_Gram42_1d
end module norm_operators_mod