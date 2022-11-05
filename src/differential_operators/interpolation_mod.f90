module interpolation_mod
use field_mod,            only : field_t
use multi_grid_field_mod, only : multi_grid_field_t
implicit none

contains

subroutine interp_1d_sbp21_2to1_ratio(in, out, direction)

  type(field_t),     intent(inout) :: in, out
  character(len=11), intent(in)    :: direction
  integer(kind=8) :: i, j

  if (direction == 'coarse2fine') then

      do i = in%is, in%ie - 1
        out%f(2 * i, out%is) = in%f(i, in%is)
        out%f(2 * i + 1, out%is) = (in%f(i, in%is) + in%f(i + 1, in%is)) / 2.0_8
      end do
      out%f(out%ie, out%is) = in%f(in%ie, in%is)

  else if (direction == 'fine2coarse') then

      do i = out%is + 1, out%ie - 1
        out%f(i, out%is) = (in%f(2 * i - 1, in%is) + 2.0_8 * in%f(2 * i, in%is) + in%f(2 * i + 1, in%is)) / 4.0_8
      end do
      out%f(out%is, out%is) = (in%f(in%ie - 1, in%is) + 2.0_8 * in%f(in%is, in%is)  + in%f(in%is + 1, in%is)) / 4.0_8
      out%f(out%ie, out%is) = (in%f(in%ie - 1, in%is) + 2.0_8 * in%f(in%ie, in%is)  + in%f(in%is + 1, in%is)) / 4.0_8

  end if

end subroutine interp_1d_sbp21_2to1_ratio

end module interpolation_mod
