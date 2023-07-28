module cell_integration_mod
use multi_domain_mod,     only : multi_domain_t
use multi_grid_field_mod, only : multi_grid_field_t
use field_mod,			  only : field_t
use domain_mod,           only : domain_t
implicit none

contains

	function integrate_medium_rectangles(in, domain) result(out)
		class(field_t), allocatable :: out
		type(field_t),  intent(in)  :: in
		type(domain_t), intent(in)  :: domain

		integer(kind=8) :: i, j

		allocate(out)
		call in%create_similar(out)

		do j = out%js, out%je
			do i = out%is, out%ie
				out%f(i, j) = in%f(i, j) * (domain%dx * domain%dy)
			end do
		end do

	end function integrate_medium_rectangles

	function integrate_multi_medium_rectangles(in, multi_domain) result(out)
		class(multi_grid_field_t), allocatable :: out
		type(multi_grid_field_t),  intent(in)  :: in
		type(multi_domain_t),      intent(in)  :: multi_domain

		integer(kind=8) :: n, m

		allocate(out)
		call in%create_similar(out)

		do m = 1, multi_domain%num_sub_y
			do n = 1, multi_domain%num_sub_x
				out%subfields(n, m) = integrate_medium_rectangles(in%subfields(n, m), multi_domain%subdomains(n, m))
			end do
		end do
	end function integrate_multi_medium_rectangles

end module cell_integration_mod