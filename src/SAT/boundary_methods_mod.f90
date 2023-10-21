module boundary_methods_mod
use field_mod,                          only : field_t
use domain_mod,							only : domain_t
use multi_grid_field_mod,               only : multi_grid_field_t
implicit none

contains

	subroutine apply_sbp21_2_boundary_method(out, in, domain, direction, sten)
		type(field_t),     intent(inout) :: out
  		type(field_t),     intent(in)    :: in
  		type(domain_t),	   intent(in)    :: domain
  		character(len=1),  intent(in)	 :: direction
  		character(len=1),  intent(in)    :: sten

  		integer(kind=8) :: i

  		if (direction == 'x') then
  			if (sten == 's') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (3.0_8 * in%f(in%is, i) - 4.0_8 * in%f(in%is + 1, i) + in%f(in%is + 2, i)) / 2.0_8 / domain%dx
	  			end do
	  		else if (sten == 'e') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (3.0_8 * in%f(in%ie, i) - 4.0_8 * in%f(in%ie - 1, i) + in%f(in%ie - 2, i)) / 2.0_8 / domain%dx
	  			end do
	  		end if
	  	else if (direction == 'y') then
	  		if (sten == 's') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (3.0_8 * in%f(i, in%js) - 4.0_8 * in%f(i, in%js + 1) + in%f(i, in%js + 2)) / 2.0_8 / domain%dy
	  			end do
	  		else if (sten == 'e') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (3.0_8 * in%f(i, in%je) - 4.0_8 * in%f(i, in%je - 1) + in%f(i, in%je - 2)) / 2.0_8 / domain%dy
	  			end do
	  		end if
	  	end if

	end subroutine apply_sbp21_2_boundary_method

	subroutine apply_sbp42_2_boundary_method(out, in, domain, direction, sten)
		type(field_t),     intent(inout) :: out
  		type(field_t),     intent(in)    :: in
  		type(domain_t),	   intent(in)    :: domain
  		character(len=1),  intent(in)	 :: direction
  		character(len=1),  intent(in)    :: sten

  		integer(kind=8) :: i

  		if (direction == 'x') then
  			if (sten == 's') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (11.0_8 * in%f(in%is, i) - 18.0_8 * in%f(in%is + 1, i) + 9.0_8 * in%f(in%is + 2, i) - 2.0_8 * in%f(in%is + 3, i) ) / 6.0_8 / domain%dx
	  			end do
	  		else if (sten == 'e') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (11.0_8 * in%f(in%ie, i) - 18.0_8 * in%f(in%ie - 1, i) + 9.0_8 * in%f(in%ie - 2, i) - 2.0_8 * in%f(in%ie - 3, i)) / 6.0_8 / domain%dx
	  			end do
	  		end if
	  	else if (direction == 'y') then
	  		if (sten == 's') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (11.0_8 * in%f(i, in%js) - 18.0_8 * in%f(i, in%js + 1) + 9.0_8 * in%f(i, in%js + 2) - 2.0_8 * in%f(i, in%js + 3) ) / 6.0_8 / domain%dy
	  			end do
	  		else if (sten == 'e') then
	  			do i = out%is, out%ie
	  				out%f(i, out%js) = (11.0_8 * in%f(i, in%je) - 18.0_8 * in%f(i, in%je - 1) + 9.0_8 * in%f(i, in%je - 2) - 2.0_8 * in%f(i, in%je - 3)) / 6.0_8 / domain%dy
	  			end do
	  		end if
	  	end if

	end subroutine apply_sbp42_2_boundary_method

	subroutine apply_sbp21_2_boundary_method_second_terms(out, in, domain, direction, sten)
		type(field_t),     intent(inout) :: out
  		type(field_t),     intent(in)    :: in
  		type(domain_t),	   intent(in)    :: domain
  		character(len=1),  intent(in)	 :: direction
  		character(len=1),  intent(in)	 :: sten

  		real(kind=8) :: dx
  		integer(kind=8) :: i

  		if (direction == 'x') then
  			dx = domain%dx
  		else if (direction == 'y') then
  			dx = domain%dy
  		end if

  		if (sten == 's') then
  			do i = out%is, out%ie
  				out%f(i, out%js)     =  3.0_8 * in%f(i, 0) / 2.0_8 / (dx ** 2.0_8) * 2.0_8
  				out%f(i, out%js + 1) = -4.0_8 * in%f(i, 0) / 2.0_8 / (dx ** 2.0_8)
  				out%f(i, out%js + 2) =  1.0_8 * in%f(i, 0) / 2.0_8 / (dx ** 2.0_8)
  				out%f(i, out%js + 3) =  0.0_8
  			end do
  		else if (sten == 'e') then
  			do i = out%is, out%ie
  				out%f(i, out%js)     =  3.0_8 * in%f(i, 0) / 2.0_8 / (dx ** 2.0_8) * 2.0_8
  				out%f(i, out%js + 1) = -4.0_8 * in%f(i, 0) / 2.0_8 / (dx ** 2.0_8)
  				out%f(i, out%js + 2) =  1.0_8 * in%f(i, 0) / 2.0_8 / (dx ** 2.0_8)
  				out%f(i, out%js + 3) =  0.0_8
  			end do
  		end if

	end subroutine apply_sbp21_2_boundary_method_second_terms

	subroutine apply_sbp42_2_boundary_method_second_terms(out, in, domain, direction, sten)
		type(field_t),     intent(inout) :: out
  		type(field_t),     intent(in)    :: in
  		type(domain_t),	   intent(in)    :: domain
  		character(len=1),  intent(in)	 :: direction
  		character(len=1),  intent(in)	 :: sten

  		real(kind=8) :: dx
  		integer(kind=8) :: i

  		if (direction == 'x') then
  			dx = domain%dx
  		else if (direction == 'y') then
  			dx = domain%dy
  		end if

  		if (sten == 's') then
  			do i = out%is, out%ie
  				out%f(i, out%js)     =  11.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 17.0_8
  				out%f(i, out%js + 1) = -18.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 59.0_8
  				out%f(i, out%js + 2) =   9.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 43.0_8
  				out%f(i, out%js + 3) = - 2.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 49.0_8
  			end do
  		else if (sten == 'e') then
  			do i = out%is, out%ie
  				out%f(i, out%js)     =  11.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 17.0_8
  				out%f(i, out%js + 1) = -18.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 59.0_8
  				out%f(i, out%js + 2) =   9.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 43.0_8
  				out%f(i, out%js + 3) =  -2.0_8 * in%f(i, 0) * 48.0_8 / 6.0_8 / (dx ** 2.0_8) / 49.0_8
  			end do
  		end if

	end subroutine apply_sbp42_2_boundary_method_second_terms

end module boundary_methods_mod
