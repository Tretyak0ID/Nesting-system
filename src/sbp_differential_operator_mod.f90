module sbp_differential_operator_mod
use differential_operator_mod, only : differential_operator_t
use field_mod,  only: field_t
use domain_mod, only: domain_t
implicit none

  type, extends(differential_operator_t) :: sbp21_t

  contains

    procedure :: apply => apply_sbp21

  end type sbp21_t

  ! type, extends(differential_operator_t) :: sbp42_t
  !
  ! contains
  !
  !   procedure :: apply => apply_sbp42
  !
  ! end type sbp_42_t

contains

  subroutine apply_sbp21(this, out, in, domain, direction)

    class    (sbp21_t),  intent(in)  :: this
    type     (field_t),  intent(out) :: out
    type     (field_t),  intent(in)  :: in
    type     (domain_t), intent(in)  :: domain
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

  end subroutine apply_sbp21

end module sbp_differential_operator_mod
