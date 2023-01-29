module boundary_methods_mod
use field_mod,                          only : field_t
use domain_mod,							only : domain_t
use multi_grid_field_mod,               only : multi_grid_field_t
implicit none

contains

	subroutine apply_sbp21_2_boundary_method(outs, oute, in, domain, direction)
		type(field_t),     intent(inout) :: outs, oute
  		type(field_t),     intent(in)    :: in
  		type(domain_t),	   intent(in)    :: domain
  		character(len=1),  intent(in)	 :: direction

  		integer(kind=8) :: i

  		if (direction == 'x') then
	  		do i = outs%is, outs%ie
	  			outs%f(i, outs%js) = (3.0_8 * in%f(in%is, i) - 4.0_8 * in%f(in%is + 1, i) + in%f(in%is + 2, i)) / 2.0_8 / domain%dx
	  		end do
	  		do i = oute%is, oute%ie
	  			oute%f(i, outs%js) = (3.0_8 * in%f(in%ie, i) - 4.0_8 * in%f(in%ie - 1, i) + in%f(in%ie - 2, i)) / 2.0_8 / domain%dx
	  		end do
	  	else if (direction == 'y') then
	  		do i = outs%is, outs%ie
	  			outs%f(i, outs%js) = (3.0_8 * in%f(i, in%js) - 4.0_8 * in%f(i, in%js + 1) + in%f(i, in%js + 2)) / 2.0_8 / domain%dy
	  		end do
	  		do i = oute%is, oute%ie
	  			oute%f(i, outs%js) = (3.0_8 * in%f(i, in%je) - 4.0_8 * in%f(i, in%je - 1) + in%f(i, in%je - 2)) / 2.0_8 / domain%dy
	  		end do
	  	end if
  		
	end subroutine apply_sbp21_2_boundary_method

	subroutine apply_sbp42_2_boundary_method(outs, oute, in, domain, direction)
		type(field_t),     intent(inout) :: outs, oute
  		type(field_t),     intent(in)    :: in
  		type(domain_t),	   intent(in)    :: domain
  		character(len=1),  intent(in)	 :: direction

  		integer(kind=8) :: i

  		if (direction == 'x') then
	  		do i = outs%is, outs%ie
	  			outs%f(i, outs%js) = (11.0_8 * in%f(in%is, i) - 18.0_8 * in%f(in%is + 1, i) + 9.0_8 * in%f(in%is + 2, i) - 2.0_8 * in%f(in%is + 3, i) ) / 6.0_8 / domain%dx
	  		end do
	  		do i = oute%is, oute%ie
	  			oute%f(i, outs%js) = (11.0_8 * in%f(in%ie, i) - 18.0_8 * in%f(in%ie - 1, i) + 9.0_8 * in%f(in%ie - 2, i) - 2.0_8 * in%f(in%ie - 3, i)) / 6.0_8 / domain%dx
	  		end do
	  	else if (direction == 'y') then
	  		do i = outs%is, outs%ie
	  			outs%f(i, outs%js) = (11.0_8 * in%f(i, in%js) - 18.0_8 * in%f(i, in%js + 1) + 9.0_8 * in%f(i, in%js + 2) - 2.0_8 * in%f(i, in%js + 3) ) / 6.0_8 / domain%dy
	  		end do
	  		do i = oute%is, oute%ie
	  			oute%f(i, outs%js) = (11.0_8 * in%f(i, in%je) - 18.0_8 * in%f(i, in%je - 1) + 9.0_8 * in%f(i, in%je - 2) - 2.0_8 * in%f(i, in%je - 3)) / 6.0_8 / domain%dy
	  		end do
	  	end if
  		
	end subroutine apply_sbp42_2_boundary_method

end module boundary_methods_mod