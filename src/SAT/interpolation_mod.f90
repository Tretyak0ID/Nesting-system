module interpolation_mod
use field_mod,                          only : field_t
use multi_grid_field_mod,               only : multi_grid_field_t
implicit none

contains

subroutine interp_identity(in, out, direction)

  type(field_t),     intent(inout) :: out
  type(field_t),     intent(in)    :: in
  character(len=*), intent(in)     :: direction
  integer(kind=8) :: i

    do i = in%is, in%ie
      out%f(i, out%js) = in%f(i, in%js)
    end do

end subroutine interp_identity

subroutine interp_MC2order_2to1ratio(in, out, direction)

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
      out%f(out%is, out%js) = (in%f(in%is, in%js) + in%f(in%is + 1, in%js)) / 2.0_8
      out%f(out%ie, out%js) = (in%f(in%ie - 1, in%js) + in%f(in%ie, in%js)) / 2.0_8

  end if

end subroutine interp_MC2order_2to1ratio

subroutine interp_MC2order_2to1ratio_periodic(in, out, direction)

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

end subroutine interp_MC2order_2to1ratio_periodic

subroutine interp_MC4order_2to1ratio_periodic(in, out, direction)
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
end subroutine interp_MC4order_2to1ratio_periodic

subroutine interp_MC4order_2to1ratio(in, out, direction)
  type(field_t),     intent(inout) :: out
  type(field_t),     intent(in)    :: in
  character(len=*), intent(in)     :: direction
  integer(kind=8) :: i, j

  if (direction == 'coarse2fine') then

      do i = in%is + 5, in%ie - 5
          out%f(2 * i, out%js)     = - 3.0_8 / 128.0_8 * in%f(i - 2, in%js) + 3.0_8 / 32.0_8 * in%f(i - 1, in%js) + 55.0_8 / 64.0_8 * in%f(i, in%js) + 3.0_8 / 32.0_8 * in%f(i + 1, in%js) - 3.0_8 / 128.0_8 * in%f(i + 2, in%js)
          out%f(2 * i + 1, out%js) = (- in%f(i - 1, in%js) + 9.0_8 * in%f(i, in%js) + 9.0_8 * in%f(i + 1, in%js) - in%f(i + 2, in%js)) / 16.0_8
      end do
      out%f(out%is, out%js)      = 603.0_8 / 641.0_8 * in%f(in%is, in%js) + 129.0_8 / 1088.0_8 * in%f(in%is + 1, in%js) - 129.0_8 / 2176.0_8 * in%f(in%is + 2, in%js) 
      out%f(out%is + 1, out%js)  = 429.0_8 / 944.0_8 * in%f(in%is, in%js) + 279.0_8 / 472.0_8 * in%f(in%is + 1, in%js) - 43.0_8 / 944.0_8 * in%f(in%is + 2, in%js)
      out%f(out%is + 2, out%js)  = 111.0_8 / 2752.0_8 * in%f(in%is, in%js) + 956.0_8 / 1071.0_8 * in%f(in%is + 1, in%js) + 3.0_8 / 32.0_8 * in%f(in%is + 2, in%js) -147.0_8 / 5504.0_8 * in%f(in%is + 3, in%js)
      out%f(out%is + 3, out%js)  = - 103.0_8 / 784.0_8 * in%f(in%is, in%js) + 549.0_8 / 784.0_8 * in%f(in%is + 1, in%js) + 387.0_8 / 784.0_8 * in%f(in%is + 2, in%js) - 1.0_8 / 16.0_8 * in%f(in%is + 3, in%js)
      out%f(out%is + 4, out%js)  = - 335.0_8 / 3072.0_8 * in%f(in%is, in%js) + 205.0_8 / 768.0_8 * in%f(in%is + 1, in%js) + 378.0_8 / 491.0_8 * in%f(in%is + 2, in%js) + 49.0_8 / 512.0_8 * in%f(in%is + 3, in%js) - 3.0_8 / 128.0_8 * in%f(in%is + 4, in%js)
      out%f(out%is + 5, out%js)  = -9.0_8 / 256.0_8 * in%f(in%is, in%js) + 5.0_8 / 256.0_8 * in%f(in%is + 1, in%js) + 129.0_8 / 256.0_8 * in%f(in%is + 2, in%js) + 147.0_8 / 256.0_8 * in%f(in%is + 3, in%js) - 1.0_8 / 16.0_8 * in%f(in%is + 4, in%js)
      out%f(out%is + 6, out%js)  = 5.0_8 / 192.0_8 * in%f(in%is, in%js) - 59.0_8 / 1024.0_8 * in%f(in%is + 1, in%js) + 43.0_8 / 512.0_8 * in%f(in%is + 2, in%js) + 722.0_8 / 823.0_8 * in%f(in%is + 3, in%js) + 3.0_8 / 32.0_8 * in%f(in%is + 4, in%js) - 3.0_8 / 128.0_8 * in%f(in%is + 5, in%js)
      out%f(out%is + 7, out%js)  = 23.0_8 / 768.0_8 * in%f(in%is, in%js) - 37.0_8 / 768.0_8 * in%f(in%is + 1, in%js) - 43.0_8 / 768.0_8 * in%f(in%is + 2, in%js) + 147.0_8 / 256.0_8 * in%f(in%is + 3, in%js) + 9.0_8 / 16.0_8 * in%f(in%is + 4, in%js) - 1.0_8 / 16.0_8 * in%f(in%is + 5, in%js)
      out%f(out%is + 8, out%js)  = 13.0_8 / 2048.0_8 * in%f(in%is, in%js) - 11.0_8 / 1024.0_8 * in%f(in%is + 1, in%js) - 43.0_8 / 2048.0_8 * in%f(in%is + 2, in%js) + 49.0_8 / 512.0_8 * in%f(in%is + 3, in%js) + 55.0_8 / 64.0_8 * in%f(in%is + 4, in%js) + 3.0_8 / 32.0_8 * in%f(in%is + 5, in%js) - 3.0_8 / 128.0_8 * in%f(in%is + 6, in%js)
      out%f(out%is + 9, out%js)  = - 1.0_8 / 384.0_8 * in%f(in%is, in%js) + 1.0_8/256.0_8 * in%f(in%is + 1, in%js) + 0.0_8 * in%f(in%is + 2, in%js) - 49.0_8/768.0_8 * in%f(in%is + 3, in%js) + 9.0_8/16.0_8 * in%f(in%is + 4, in%js) + 9.0_8 / 16.0_8 * in%f(in%is + 5, in%js) - 1.0_8/16.0_8 * in%f(in%is + 6, in%js)
      out%f(out%is + 10, out%js) = - 1.0_8/1024.0_8 * in%f(in%is, in%js) + 3.0_8/2048.0_8 * in%f(in%is + 1, in%js) + 0.0_8 * in%f(in%is + 2, in%js) - 49.0_8/2048.0_8 * in%f(in%is + 3, in%js) + 3.0_8/32.0_8 * in%f(in%is + 4, in%js) + 55.0_8/64.0_8 * in%f(in%is + 5, in%js) + 3.0_8/32.0_8 * in%f(in%is + 6, in%js) - 3.0_8/128.0_8 * in%f(in%is + 7, in%js)
      out%f(out%ie - 10, out%js) = -3.0_8/128.0_8 * in%f(in%ie - 7, in%js) + 3.0_8/32.0_8 * in%f(in%ie - 6, in%js) + 55.0_8/64.0_8 * in%f(in%ie - 5, in%js) + 3.0_8/32.0_8 * in%f(in%ie - 4, in%js) - 49.0_8/2048.0_8 * in%f(in%ie - 3, in%js) + 0.0_8 * in%f(in%ie - 2, in%js) + 3.0_8/2048.0_8 * in%f(in%ie - 1, in%js) - 1.0_8/1024.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 9, out%js) = -1.0_8/16.0_8 * in%f(in%ie - 6, in%js) + 9.0_8/16.0_8 * in%f(in%ie - 5, in%js) + 9.0_8/16.0_8 * in%f(in%ie - 4, in%js) - 49.0_8/768.0_8 * in%f(in%ie - 3, in%js) + 0.0_8 * in%f(in%ie - 2, in%js) + 1.0_8/256.0_8 * in%f(in%ie - 1, in%js) - 1.0_8/384.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 8, out%js) = -3.0_8/128.0_8 * in%f(in%ie - 6, in%js) + 3.0_8/32.0_8 * in%f(in%ie - 5, in%js) + 55.0_8/64.0_8 * in%f(in%ie - 4, in%js) + 49.0_8/512.0_8 * in%f(in%ie - 3, in%js) - 43.0_8/2048.0_8 * in%f(in%ie - 2, in%js) - 11.0_8/1024.0_8 * in%f(in%ie - 1, in%js) + 13.0_8/2048.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 7, out%js) = -1.0_8/16.0_8 * in%f(in%ie - 5, in%js) +  9.0_8/16.0_8 * in%f(in%ie - 4, in%js) + 147.0_8/256.0_8 * in%f(in%ie - 3, in%js) - 43.0_8/768.0_8 * in%f(in%ie - 2, in%js) - 37.0_8/768.0_8 * in%f(in%ie - 1, in%js) + 23.0_8/768.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 6, out%js) = -3.0_8/128.0_8 * in%f(in%ie - 5, in%js) + 3.0_8/32.0_8 * in%f(in%ie - 4, in%js) + 722.0_8/823.0_8 * in%f(in%ie - 3, in%js) + 43.0_8/512.0_8 * in%f(in%ie - 2, in%js) - 59.0_8/1024.0_8 * in%f(in%ie - 1, in%js) + 5.0_8/192.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 5, out%js) = -1.0_8/16.0_8 * in%f(in%ie - 4, in%js) + 147.0_8/256.0_8 * in%f(in%ie - 3, in%js) + 129.0_8/256.0_8 * in%f(in%ie - 2, in%js) + 5.0_8/256.0_8 * in%f(in%ie - 1, in%js) - 9.0_8/256.0_8 * in%f(in%ie, in%js) 
      out%f(out%ie - 4, out%js) = -3.0_8/128.0_8 * in%f(in%ie - 4, in%js) + 49.0_8/512.0_8 * in%f(in%ie - 3, in%js) + 378.0_8/491.0_8 * in%f(in%ie - 2, in%js) + 205.0_8/768.0_8 * in%f(in%ie - 1, in%js) - 335.0_8/3072.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 3, out%js) = -1.0_8/16.0_8 * in%f(in%ie -3 , in%js) + 387.0_8/784.0_8 * in%f(in%ie - 2, in%js) + 549.0_8/784.0_8 * in%f(in%ie - 1, in%js) -103.0_8/784.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 2, out%js) = -147.0_8/5504.0_8 * in%f(in%ie - 3, in%js) + 3.0_8/32.0_8 * in%f(in%ie - 2, in%js) + 956.0_8/1071.0_8 * in%f(in%ie - 1, in%js) + 111.0_8/2752.0_8 * in%f(in%ie, in%js)
      out%f(out%ie - 1, out%js) = -43.0_8 / 944.0_8 * in%f(in%ie - 2, in%js) + 279.0_8 / 472.0_8 * in%f(in%ie - 1, in%js) + 429.0_8 / 944.0_8 * in%f(in%ie, in%js) 
      out%f(out%ie, out%js)     = -129.0_8 / 2176.0_8 * in%f(in%ie - 2, in%js) + 129.0_8 / 1088.0_8 * in%f(in%ie - 1, in%js) + 603.0_8 / 641.0_8 * in%f(in%ie, in%js)

  else if (direction == 'fine2coarse') then

      do i = out%is + 2, out%ie - 2
          out%f(i, out%js) = - 3.0_8 / 256.0_8 * in%f(2*i - 4, in%js) - 1.0_8 / 32.0_8 * in%f(2*i - 3, in%js) + 3.0_8 / 64.0_8 * in%f(2*i - 2, in%js) + 9.0_8 / 32.0_8 * in%f(2*i - 1, in%js) + 55.0_8 / 128.0_8 * in%f(2*i, in%js) - 3.0_8 / 256.0_8 * in%f(2*i + 4, in%js) - 1.0_8 / 32.0_8 * in%f(2*i + 3, in%js) + 3.0_8 / 64.0_8 * in%f(2*i + 2, in%js) + 9.0_8 / 32.0_8 * in%f(2*i + 1, in%js)
      end do


      out%f(out%is, out%js)     = 722.0_8/1535.0_8 * in%f(in%is, in%js)   +    429.0_8/544.0_8 * in%f(in%is + 1, in%js)     +   111.0_8/2176.0_8 * in%f(in%is + 2, in%js)      -103.0_8/544.0_8 * in%f(in%is + 3, in%js)       -335.0_8/2176.0_8 * in%f(in%is + 4, in%js)       -27.0_8/544.0_8 * in%f(in%is + 5, in%js)     +     5.0_8/136.0_8 * in%f(in%is + 6, in%js)    +     23.0_8/544.0_8 * in%f(in%is + 7, in%js)    +     39.0_8/4352.0_8 * in%f(in%is + 8, in%js)        -1.0_8/272.0_8 * in%f(in%is + 9, in%js)         -3.0_8/2176.0_8 * in%f(in%is + 10, in%js)
      out%f(out%is + 1, out%js) = 129.0_8/7552.0_8 * in%f(in%is, in%js)    +   279.0_8/944.0_8 * in%f(in%is + 1, in%js)     +   848.0_8/2607.0_8 * in%f(in%is + 2, in%js)   +    549.0_8/1888.0_8 * in%f(in%is + 3, in%js)    +   205.0_8/1888.0_8 * in%f(in%is + 4, in%js)    +    15.0_8/1888.0_8 * in%f(in%is + 5, in%js)        -3.0_8/128.0_8 * in%f(in%is + 6, in%js)        -37.0_8/1888.0_8 * in%f(in%is + 7, in%js)       -33.0_8/7552.0_8 * in%f(in%is + 8, in%js)     +    3.0_8/1888.0_8 * in%f(in%is + 9, in%js)     +    9.0_8/15104.0_8 * in%f(in%is + 10, in%js)
      out%f(out%ie - 1, out%js) = 9.0_8/15104.0_8 * in%f(in%ie - 10, in%js)    +    3.0_8/1888.0_8 * in%f(in%ie - 9, in%js)       -33.0_8/7552.0_8 * in%f(in%ie - 8, in%js)       -37.0_8/1888.0_8 * in%f(in%ie - 7, in%js)        -3.0_8/128.0_8 * in%f(in%ie - 6, in%js)      +   15.0_8/1888.0_8 * in%f(in%ie - 5, in%js)     +  205.0_8/1888.0_8 * in%f(in%ie - 4, in%js)   +    549.0_8/1888.0_8 * in%f(in%ie - 3, in%js)   +    848.0_8/2607.0_8 * in%f(in%ie - 2, in%js)   +    279.0_8/944.0_8 * in%f(in%ie - 1, in%js)    +    129.0_8/7552.0_8 * in%f(in%ie, in%js)
      out%f(out%ie, out%js)     = -3.0_8/2176.0_8 * in%f(in%ie - 10, in%js)        -1.0_8/272.0_8 * in%f(in%ie - 9, in%js)     +    39.0_8/4352.0_8 * in%f(in%ie - 8, in%js)    +    23.0_8/544.0_8 * in%f(in%ie - 7, in%js)    +     5.0_8/136.0_8 * in%f(in%ie - 6, in%js)        -27.0_8/544.0_8 * in%f(in%ie - 5, in%js)       -335.0_8/2176.0_8 * in%f(in%ie - 4, in%js)      -103.0_8/544.0_8 * in%f(in%ie - 3, in%js)    +   111.0_8/2176.0_8 * in%f(in%ie - 2, in%js)    +  429.0_8/544.0_8 * in%f(in%ie - 1, in%js)    +   722.0_8/1535.0_8 * in%f(in%ie, in%js)

  end if
end subroutine interp_MC4order_2to1ratio

end module interpolation_mod
