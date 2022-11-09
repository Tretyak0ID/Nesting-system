module interpolation_mod
use field_mod,            only : field_t
use multi_grid_field_mod, only : multi_grid_field_t
implicit none

contains

subroutine identity(in, out, direction)

  type(field_t),     intent(inout) :: out
  type(field_t),     intent(in)    :: in
  character(len=*), intent(in)     :: direction
  integer(kind=8) :: i

    do i = in%is, in%ie
      out%f(i, out%js) = in%f(i, in%js)
    end do

end subroutine identity

subroutine interp_1d_sbp21_2to1_ratio(in, out, direction)

  type(field_t),     intent(inout) :: out
  type(field_t),     intent(in)    :: in
  character(len=*), intent(in)     :: direction
  integer(kind=8) :: i, j

  if (direction == 'coarse2fine') then

      do i = in%is, in%ie - 1
        out%f(2 * i, out%js)     = in%f(i, in%js)
        out%f(2 * i + 1, out%js) = (in%f(i, in%js) + in%f(i + 1, in%js)) / 2.0_8
      end do
      out%f(out%ie, out%js) = in%f(in%ie, in%js)

  else if (direction == 'fine2coarse') then

      do i = out%is + 1, out%ie - 1
        out%f(i, out%js) = (in%f(2 * i - 1, in%js) + 2.0_8 * in%f(2 * i, in%js) + in%f(2 * i + 1, in%js)) / 4.0_8
      end do
      out%f(out%is, out%js) = (in%f(in%ie - 1, in%js) + 2.0_8 * in%f(in%is, in%js)  + in%f(in%is + 1, in%js)) / 4.0_8
      out%f(out%ie, out%js) = (in%f(in%ie - 1, in%js) + 2.0_8 * in%f(in%ie, in%js)  + in%f(in%is + 1, in%js)) / 4.0_8

  end if

end subroutine interp_1d_sbp21_2to1_ratio

subroutine interp_1d_sbp42_2to1_ratio(in, out, direction)
  type(field_t),     intent(inout) :: out
  type(field_t),     intent(in)    :: in
  character(len=*), intent(in)     :: direction
  integer(kind=8) :: i, j

  if (direction == 'coarse2fine') then

      do i = in%is + 1, in%ie - 2
          out%f(2 * i, out%js)     = in%f(i, in%js)
          out%f(2 * i + 1, out%js) = (- in%f(i - 1, in%js) + 9.0_8 * in%f(i, in%js) + 9.0_8 * in%f(i + 1, in%js) - in%f(i + 2, in%js)) / 16.0_8
      end do
      out%f(out%is, out%js)     = in%f(in%is, in%js)
      out%f(out%is + 1, out%js) = (- in%f(in%ie - 1, in%js) + 9.0_8 * in%f(in%is, in%js) + 9.0_8 * in%f(in%is + 1, in%js) - in%f(in%is + 2, in%js)) / 16.0_8
      out%f(out%ie - 1, out%js) = (- in%f(in%ie - 2, in%js) + 9.0_8 * in%f(in%ie - 1, in%js) + 9.0_8 * in%f(in%ie, in%js) - in%f(in%is + 1, in%js)) / 16.0_8
      out%f(out%ie - 2, out%js) = in%f(in%ie - 1, in%js)
      out%f(out%ie, out%js)     = in%f(in%ie, in%js)

  else if (direction == 'fine2coarse') then

      do i = out%is + 3, out%ie - 2
          out%f(i, out%js) = (-in%f(2 * i - 3, in%js) +9 * in%f(2*i-1, in%js) + 16*in%f(2*i, in%js)+ 9*in%f(2*i+1, in%js) - in%f(2*i+3, in%js)) / 32.0_8
      end do
      out%f(out%is, out%js)     = (- in%f(in%ie - 3, in%js) + 9.0_8 * in%f(in%ie - 1, in%js) + 16.0_8 * in%f(in%is + 0, in%js) + 9.0_8 * in%f(in%is + 1, in%js) - in%f(in%is + 3, in%js)) / 32.0_8
      out%f(out%is + 1, out%js) = (- in%f(in%ie - 1, in%js) + 9.0_8 * in%f(in%is + 1, in%js) + 16.0_8 * in%f(in%is + 2, in%js) + 9.0_8 * in%f(in%is + 3, in%js) - in%f(in%is + 5, in%js)) / 32.0_8
      out%f(out%is + 2, out%js) = (- in%f(in%is + 1, in%js) + 9.0_8 * in%f(in%is + 3, in%js) + 16.0_8 * in%f(in%is + 4, in%js) + 9.0_8 * in%f(in%is + 5, in%js) - in%f(in%is + 7, in%js)) / 32.0_8
      out%f(out%ie - 2, out%js) = (- in%f(in%ie - 7, in%js) + 9.0_8 * in%f(in%ie - 5, in%js) + 16.0_8 * in%f(in%ie - 4, in%js) + 9.0_8 * in%f(in%ie - 3, in%js) - in%f(in%ie - 1, in%js)) / 32.0_8
      out%f(out%ie - 1, out%js) = (- in%f(in%ie - 5, in%js) + 9.0_8 * in%f(in%ie - 3, in%js) + 16.0_8 * in%f(in%ie - 2, in%js) + 9.0_8 * in%f(in%ie - 1, in%js) - in%f(in%is + 1, in%js)) / 32.0_8
      out%f(out%ie, out%js)     = (- in%f(in%ie - 3, in%js) + 9.0_8 * in%f(in%ie - 1, in%js) + 16.0_8 * in%f(in%ie - 0, in%js) + 9.0_8 * in%f(in%is + 1, in%js) - in%f(in%is + 3, in%js)) / 32.0_8

  end if
end subroutine

end module interpolation_mod
