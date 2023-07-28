program interpolation_test
use multi_grid_field_mod, only : multi_grid_field_t
use multi_domain_mod,     only : multi_domain_t
use domain_mod,           only : domain_t
use field_mod,            only : field_t
use cell_integration_mod, only : integrate_multi_medium_rectangles, integrate_medium_rectangles
use cell_interpolation_mod, only : interpolate_cell_centers, interpolate_cell_faces, interpolate_faces_centers
use const_mod,            only: pi
implicit none

	type(field_t)            :: in_field, out_field, cos_field
	type(multi_grid_field_t), allocatable :: diff(:), cent_diff(:)
  	type(domain_t)           :: domain
  	type(multi_grid_field_t) :: inf, outf
  	type(multi_domain_t)     :: multi_domain
  	integer        :: i, j, n, m
  	integer(kind=4), allocatable :: deg(:, :)

  	call domain%init(0.0_8, 128.0_8, 0, 64, 0.0_8, 128.0_8, 0, 64)
  	call in_field%init(0, 64, 0, 64)
  	call cos_field%init(0, 64, 0, 64)

  	allocate(diff(1:1))
  	allocate(cent_diff(1:1))

  	allocate(deg(1:1, 1:1))
  	deg(1, 1) = 1

  	call multi_domain%init(domain, 1, 1, deg)
  	call inf%init(multi_domain)

  	do m = 1, multi_domain%num_sub_y
  		do n = 1, multi_domain%num_sub_x
  			do j = 0, multi_domain%subdomains(n, m)%ny
    			do i = 0, multi_domain%subdomains(n, m)%nx
      				inf%subfields(n, m)%f(i, j) = sin(multi_domain%subdomains(n, m)%x(i))
    			end do
  			end do
  		end do
  	end do

  	outf = interpolate_cell_centers(inf, multi_domain)

  	diff = interpolate_cell_faces(inf, outf, multi_domain)

  	cent_diff = interpolate_faces_centers(diff, multi_domain)

  	print *, cent_diff(1)%subfields(1,1)%f

end program interpolation_test