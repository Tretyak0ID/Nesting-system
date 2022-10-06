module sbp_operators_mod
use field_mod,  only: field_t
use domain_mod, only: domain_t
implicit none

contains

  subroutine sbp21(out, in, direction, domain)
    class    (field_t),  intent(in)  :: in
    class    (domain_t), intent(in)  :: domain
    class    (field_t),  intent(out) :: out
    character(len=1),    intent(in)  :: direction
    integer  (kind=8)                :: i, j

    call out%init(0, domain%nx, 0, domain%ny)

    if (direction == 'x') then

      do i = 1, domain%nx - 1
        do j = 0, domain%ny
          out%f(i, j) = (in%f(i + 1, j) - in%f(i - 1, j)) / 2.0_8 / domain%dx
        end do
      end do

      do j = 0, domain%ny
        out%f(domain%nx, j) = (in%f(domain%nx, j) - in%f(domain%nx - 1, j)) / domain%dx
        out%f(0, j)         = (in%f(1, j) - in%f(0, j)) / domain%dx
      end do

    else if (direction == 'y') then

      do i = 0, domain%nx
        do j = 1, domain%ny - 1
          out%f(i, j) = (in%f(i, j + 1) - in%f(i, j - 1)) / 2.0_8 / domain%dy
        end do
      end do

      do i = 0, domain%nx
        out%f(i, domain%ny) = (in%f(i, domain%ny) - in%f(i, domain%ny - 1)) / domain%dy
        out%f(i, 0)         = (in%f(i, 1) - in%f(i, 0)) / domain%dy
      end do

    else
      print *, 'Error in diff_sbp21. Wrong direction value'
    end if

  end subroutine sbp21

end module sbp_operators_mod
