module cell_interpolation_mod
use multi_domain_mod,     only : multi_domain_t
use multi_grid_field_mod, only : multi_grid_field_t
use laplacian_mod,        only : modify_cell_FV_laplacian, modify_faces_FV_laplacian
use field_mod,			  only : field_t
use domain_mod,           only : domain_t
implicit none

contains

	function interpolate_cell_centers(cell_volume, multi_domain) result(centers_cell)

		class(multi_grid_field_t), allocatable :: centers_cell
		type(multi_grid_field_t),  intent(in)  :: cell_volume
		type(multi_domain_t),      intent(in)  :: multi_domain

		integer(kind=8) :: i, j, n, m

		allocate(centers_cell)
		call cell_volume%create_similar(centers_cell)

		do m = 1, multi_domain%num_sub_y
			do n = 1, multi_domain%num_sub_x
				do j = centers_cell%subfields(n, m)%js, centers_cell%subfields(n, m)%je
					do i = centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%ie
						centers_cell%subfields(n, m)%f(i, j) = cell_volume%subfields(n, m)%f(i, j)
					end do
				end do
			end do 
		end do

		call modify_cell_FV_laplacian(multi_domain%subdomains(1, 1)%dx ** 2 / 24.0_8, centers_cell, cell_volume, multi_domain)

	end function interpolate_cell_centers

	function interpolate_cell_faces(cell_volume, cell_centers, multi_domain) result(cell_faces_volume)

		class(multi_grid_field_t),  allocatable :: cell_faces_volume(:)
		type(multi_grid_field_t),  intent(in)  :: cell_volume, cell_centers
		type(multi_domain_t),      intent(in)  :: multi_domain

		integer(kind=8) :: i, j, n, m

		allocate(cell_faces_volume(1:2))
		call cell_centers%create_similar(cell_faces_volume(1))
		call cell_centers%create_similar(cell_faces_volume(2))

		do m = 1, multi_domain%num_sub_y
		do n = 1, multi_domain%num_sub_x
			do j = multi_domain%subdomains(n, m)%js + 1, multi_domain%subdomains(n, m)%je - 2
			do i = multi_domain%subdomains(n, m)%is + 1, multi_domain%subdomains(n, m)%ie - 2

				cell_faces_volume(1)%subfields(n, m)%f(i, j) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, j) + cell_volume%subfields(n, m)%f(i + 1, j)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i - 1, j) + cell_volume%subfields(n, m)%f(i + 2, j))
				cell_faces_volume(2)%subfields(n, m)%f(i, j) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, j) + cell_volume%subfields(n, m)%f(i, j + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, j - 1) + cell_volume%subfields(n, m)%f(i, j + 2))

				! y-direction
				cell_faces_volume(2)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js)     = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js + 2))
				cell_faces_volume(2)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je)     = cell_faces_volume(2)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js)
				cell_faces_volume(2)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je - 1) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je - 2) + cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js + 1))

				! x-direction
				cell_faces_volume(1)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js)         = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(i + 1, cell_faces_volume(1)%subfields(n, m)%js)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i - 1, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(i + 2, cell_faces_volume(1)%subfields(n, m)%js))
				cell_faces_volume(1)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je)         = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je) + cell_volume%subfields(n, m)%f(i + 1, cell_faces_volume(1)%subfields(n, m)%je)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i - 1, cell_faces_volume(1)%subfields(n, m)%je) + cell_volume%subfields(n, m)%f(i + 2, cell_faces_volume(1)%subfields(n, m)%je))
				cell_faces_volume(1)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je - 1)     = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(i + 1, cell_faces_volume(1)%subfields(n, m)%je - 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(i - 1, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(i + 2, cell_faces_volume(1)%subfields(n, m)%je - 1))
			end do

				! x-direction
				cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j)     = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 1, j)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 2, j))
				cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, j)     = cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j)
				cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, j)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 2, j) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 1, j))

				!y-direction
				cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j)      = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j + 2))
				cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, j)      = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, j) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, j + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, j - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, j + 2))
				cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j)  = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, j + 2))
			end do

			! x-direction
			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js)     = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 1, cell_faces_volume(1)%subfields(n, m)%js)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 2, cell_faces_volume(1)%subfields(n, m)%js))
			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 1, cell_faces_volume(1)%subfields(n, m)%je - 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 2, cell_faces_volume(1)%subfields(n, m)%je - 1))
			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je)     = cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js)

			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%js)     = cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js)
			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%je - 1) = cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1)
			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%je)     = cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je)

			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js)     = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%js)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 2, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 1, cell_faces_volume(1)%subfields(n, m)%js))
			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je)     = cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js)
			cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je - 1) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%je - 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 2, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is + 1, cell_faces_volume(1)%subfields(n, m)%je - 1))
			
			! y-direction
			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js)     = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js + 2))
			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js + 2))
			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%js)     = cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js)

			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je)     = cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js)    
			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je) = cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js)
			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%je)     = cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%js)

			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1)     = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js + 2))
			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je - 1) = 7.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js + 1)) - 1.0_8 / 12.0_8 * (cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%je - 1) + cell_volume%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie - 1, cell_faces_volume(1)%subfields(n, m)%js + 2))
			cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%ie, cell_faces_volume(1)%subfields(n, m)%je - 1)     = cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1)    
		end do
		end do

	end function interpolate_cell_faces

	function interpolate_faces_centers(cell_faces_volume, multi_domain) result(cell_faces_center)
		class(multi_grid_field_t), allocatable              :: cell_faces_center(:)
		type(multi_grid_field_t),  dimension(2), intent(in) :: cell_faces_volume
		type(multi_domain_t),                    intent(in) :: multi_domain

		integer(kind=8) :: i, j, n, m

		allocate(cell_faces_center(1:2))
		call cell_faces_volume(1)%create_similar(cell_faces_center(1))
		call cell_faces_volume(2)%create_similar(cell_faces_center(2))

		do m = 1, multi_domain%num_sub_y
		do n = 1, multi_domain%num_sub_x
			do j = multi_domain%subdomains(n, m)%js, multi_domain%subdomains(n, m)%je
			do i = multi_domain%subdomains(n, m)%is, multi_domain%subdomains(n, m)%ie
				cell_faces_center(1)%subfields(n, m)%f(i ,j) = cell_faces_volume(1)%subfields(n, m)%f(i ,j)
				cell_faces_center(2)%subfields(n, m)%f(i ,j) = cell_faces_volume(2)%subfields(n, m)%f(i ,j)
			end do
			end do
		end do
		end do

		call modify_faces_FV_laplacian(-1.0_8 *  multi_domain%subdomains(1, 1)%dx ** 2.0_8 / 24.0_8, cell_faces_center, cell_faces_volume, multi_domain)
	end function interpolate_faces_centers
	
end module cell_interpolation_mod