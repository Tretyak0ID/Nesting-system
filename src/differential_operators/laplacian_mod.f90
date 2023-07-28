module laplacian_mod
use multi_domain_mod,          only : multi_domain_t
use multi_grid_field_mod,      only : multi_grid_field_t
use differential_operator_mod, only : differential_operator_t
use SAT_mod,                   only : sbp_SAT_penalty_two_block_diffusion
implicit none

contains

  subroutine calc_laplacian(out, in, multi_domain, coefs, diff2_op)
    type (multi_grid_field_t),      intent(inout) :: out
    type (multi_grid_field_t),      intent(in)    :: in
    type (multi_domain_t),          intent(in)    :: multi_domain
    real (kind=8), allocatable,     intent(in)    :: coefs(:, :)
    class(differential_operator_t), intent(in)    :: diff2_op

    type(multi_grid_field_t)  :: dinx, diny
    integer(kind=4) :: n, m

    call dinx%init(multi_domain)
    call diny%init(multi_domain)

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        call diff2_op%apply(dinx%subfields(n, m), in%subfields(n, m), multi_domain%subdomains(n, m), 'x')
        call diff2_op%apply(diny%subfields(n, m), in%subfields(n, m), multi_domain%subdomains(n, m), 'y')
        call out%subfields(n, m)%assign(coefs(n, m), dinx%subfields(n, m), coefs(n, m), diny%subfields(n, m), multi_domain%subdomains(n, m))
      end do
    end do

    !call out%assign(coefs(1, 1), dinx, coefs(1, 1), diny, multi_domain)

    call sbp_SAT_penalty_two_block_diffusion(out, in, multi_domain, coefs, 'x', diff2_op%name)
    call sbp_SAT_penalty_two_block_diffusion(out, in, multi_domain, coefs, 'y', diff2_op%name)

  end subroutine calc_laplacian

  subroutine modify_cell_FV_laplacian(scalar, centers_cell, cell_volume, multi_domain)

    type(multi_grid_field_t), intent(inout) :: centers_cell
    type(multi_grid_field_t), intent(in)    :: cell_volume
    type(multi_domain_t),     intent(in)    :: multi_domain
    real(kind=8),             intent(in)    :: scalar

    integer(kind=8) :: i, j, n, m

    do m = 1, multi_domain%num_sub_y
    do n = 1, multi_domain%num_sub_x
      do j = multi_domain%subdomains(n, m)%js + 1, multi_domain%subdomains(n, m)%je - 1
      do i = multi_domain%subdomains(n, m)%is + 1, multi_domain%subdomains(n, m)%ie - 1
        centers_cell%subfields(n, m)%f(i, j) = centers_cell%subfields(n, m)%f(i, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_volume%subfields(n, m)%f(i - 1, j) - 2.0_8 * cell_volume%subfields(n, m)%f(i, j) + cell_volume%subfields(n, m)%f(i + 1, j)))
        centers_cell%subfields(n, m)%f(i, j) = centers_cell%subfields(n, m)%f(i, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_volume%subfields(n, m)%f(i, j - 1) - 2.0_8 * cell_volume%subfields(n, m)%f(i, j) + cell_volume%subfields(n, m)%f(i, j + 1)))

        centers_cell%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js) = centers_cell%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_volume%subfields(n, m)%f(i - 1, centers_cell%subfields(n, m)%js) - 2.0_8 * cell_volume%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(i + 1, centers_cell%subfields(n, m)%js)))
        centers_cell%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js) = centers_cell%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_volume%subfields(n, m)%f(i, centers_cell%subfields(n, m)%je - 1) - 2.0_8 * cell_volume%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js + 1)))
        centers_cell%subfields(n, m)%f(i, centers_cell%subfields(n, m)%je) = centers_cell%subfields(n, m)%f(i, centers_cell%subfields(n, m)%js)
      end do
        centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%ie - 1, j) - 2.0_8 * cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j) + cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is + 1, j)))
        centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j - 1) - 2.0_8 * cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j) + cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j + 1)))
        centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%ie, j) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, j)
      end do
      centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%ie - 1, centers_cell%subfields(n, m)%js) - 2.0_8 * cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is + 1, centers_cell%subfields(n, m)%js)))
      centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%je - 1) - 2.0_8 * cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js) + cell_volume%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js + 1)))
      centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%ie, centers_cell%subfields(n, m)%js) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js)
      centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%je) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js)
      centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%ie, centers_cell%subfields(n, m)%je) = centers_cell%subfields(n, m)%f(centers_cell%subfields(n, m)%is, centers_cell%subfields(n, m)%js)
    end do
    end do

  end subroutine modify_cell_FV_laplacian

  subroutine modify_faces_FV_laplacian(scalar, cell_faces_center, cell_faces_volume, multi_domain)
    type(multi_grid_field_t), intent(inout), dimension(2) :: cell_faces_center
    type(multi_grid_field_t), intent(in),    dimension(2) :: cell_faces_volume
    type(multi_domain_t),     intent(in)                  :: multi_domain
    real(kind=8),             intent(in)                  :: scalar

    integer(kind=8) :: n, m, i, j

    do m = 1, multi_domain%num_sub_y
    do n = 1, multi_domain%num_sub_x
      do j = multi_domain%subdomains(n, m)%js + 1, multi_domain%subdomains(n, m)%je - 1
      do i = multi_domain%subdomains(n, m)%is + 1, multi_domain%subdomains(n, m)%ie - 1
        cell_faces_center(2)%subfields(n, m)%f(i, j) = cell_faces_center(2)%subfields(n, m)%f(i, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_faces_volume(2)%subfields(n, m)%f(i - 1, j) - 2.0_8 * cell_faces_volume(2)%subfields(n, m)%f(i, j) + cell_faces_volume(2)%subfields(n, m)%f(i + 1, j)))
        cell_faces_center(1)%subfields(n, m)%f(i, j) = cell_faces_center(1)%subfields(n, m)%f(i, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_faces_volume(1)%subfields(n, m)%f(i, j - 1) - 2.0_8 * cell_faces_volume(1)%subfields(n, m)%f(i, j) + cell_faces_volume(1)%subfields(n, m)%f(i, j + 1)))

        cell_faces_center(2)%subfields(n, m)%f(i, cell_faces_center(2)%subfields(n, m)%js) = cell_faces_center(2)%subfields(n, m)%f(i, cell_faces_center(2)%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_faces_volume(2)%subfields(n, m)%f(i - 1, cell_faces_volume(2)%subfields(n, m)%js) - 2.0_8 * cell_faces_volume(2)%subfields(n, m)%f(i, cell_faces_volume(2)%subfields(n, m)%js) + cell_faces_volume(2)%subfields(n, m)%f(i + 1, cell_faces_volume(2)%subfields(n, m)%js)))
        cell_faces_center(1)%subfields(n, m)%f(i, cell_faces_center(1)%subfields(n, m)%js) = cell_faces_center(1)%subfields(n, m)%f(i, cell_faces_center(1)%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_faces_volume(1)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%je - 1) - 2.0_8 * cell_faces_volume(1)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js) + cell_faces_volume(1)%subfields(n, m)%f(i, cell_faces_volume(1)%subfields(n, m)%js + 1)))
        cell_faces_center(2)%subfields(n, m)%f(i, cell_faces_center(2)%subfields(n, m)%je) = cell_faces_center(2)%subfields(n, m)%f(i, cell_faces_center(2)%subfields(n, m)%js)
        cell_faces_center(1)%subfields(n, m)%f(i, cell_faces_center(1)%subfields(n, m)%je) = cell_faces_center(1)%subfields(n, m)%f(i, cell_faces_center(1)%subfields(n, m)%js)
      end do
        cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, j) = cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(2)%subfields(n, m)%ie - 1, j) - 2.0_8 * cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(2)%subfields(n, m)%is, j) + cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(2)%subfields(n, m)%is + 1, j)))
        cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, j) = cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, j) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j - 1) - 2.0_8 * cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j) + cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, j + 1)))
        cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%ie, j) = cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, j)
        cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%ie, j) = cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, j)
      end do
      cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, cell_faces_center(2)%subfields(n, m)%js) = cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, cell_faces_center(2)%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dx ** 2.0_8 * (cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(2)%subfields(n, m)%ie - 1, cell_faces_volume(2)%subfields(n, m)%js) - 2.0_8 * cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(2)%subfields(n, m)%is, cell_faces_volume(2)%subfields(n, m)%js) + cell_faces_volume(2)%subfields(n, m)%f(cell_faces_volume(2)%subfields(n, m)%is + 1, cell_faces_volume(2)%subfields(n, m)%js)))
      cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, cell_faces_center(1)%subfields(n, m)%js) = cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, cell_faces_center(1)%subfields(n, m)%js) + scalar * (1.0_8 / multi_domain%subdomains(n, m)%dy ** 2.0_8 * (cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%je - 1) - 2.0_8 * cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js) + cell_faces_volume(1)%subfields(n, m)%f(cell_faces_volume(1)%subfields(n, m)%is, cell_faces_volume(1)%subfields(n, m)%js + 1)))
      cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%ie, cell_faces_center(2)%subfields(n, m)%js) = cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, cell_faces_center(2)%subfields(n, m)%js)
      cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, cell_faces_center(2)%subfields(n, m)%je) = cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, cell_faces_center(2)%subfields(n, m)%js)
      cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%ie, cell_faces_center(2)%subfields(n, m)%je) = cell_faces_center(2)%subfields(n, m)%f(cell_faces_center(2)%subfields(n, m)%is, cell_faces_center(2)%subfields(n, m)%js)
      cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%ie, cell_faces_center(1)%subfields(n, m)%js) = cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, cell_faces_center(1)%subfields(n, m)%js)
      cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, cell_faces_center(1)%subfields(n, m)%je) = cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, cell_faces_center(1)%subfields(n, m)%js)
      cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%ie, cell_faces_center(1)%subfields(n, m)%je) = cell_faces_center(1)%subfields(n, m)%f(cell_faces_center(1)%subfields(n, m)%is, cell_faces_center(1)%subfields(n, m)%js)
    end do
    end do

  end subroutine modify_faces_FV_laplacian

end module laplacian_mod
