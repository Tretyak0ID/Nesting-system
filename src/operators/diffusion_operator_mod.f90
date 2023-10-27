module diffusion_operator_mod
  use operator_mod,              only: operator_t
  use differential_operator_mod, only: differential_operator_t
  use multi_domain_mod,          only: multi_domain_t
  use stvec_mod,                 only: stvec_t
  use stvec_swe_mod,             only: stvec_swe_t
  use multi_domain_mod,          only: multi_domain_t
  use multi_grid_field_mod,      only: multi_grid_field_t
  use laplacian_mod,             only: calc_laplacian
implicit none


  type, public, extends(operator_t) :: diffusion_operator_t
  !Diffusion operator
    class(differential_operator_t), allocatable :: diff2_op
    real(kind=8) , allocatable                  :: coefs(:, :)

    type(multi_grid_field_t) :: lap_u_buff, lap_v_buff, lap_h_buff, lap_u, lap_v, lap_h

  contains

    procedure :: init => init_diffusion
    procedure :: apply => apply_diffusion

  end type diffusion_operator_t


contains

  subroutine init_diffusion(this, diff2_op, coefs, multi_domain)

    class(diffusion_operator_t),    intent(inout) :: this
    class(differential_operator_t), intent(in)    :: diff2_op
    class(multi_domain_t),          intent(in)    :: multi_domain
    real(kind=8), allocatable,      intent(in)    :: coefs(:, :)

    integer(kind=4) :: n, m

    this%diff2_op = diff2_op

    allocate(this%coefs(1 : multi_domain%num_sub_x, 1 : multi_domain%num_sub_y))

    do n = 1, multi_domain%num_sub_x
      do m = 1, multi_domain%num_sub_y
        this%coefs(n, m) = coefs(n, m)
      end do
    end do

    call this%lap_u%init(multi_domain)
    call this%lap_v%init(multi_domain)
    call this%lap_h%init(multi_domain)
    call this%lap_u_buff%init(multi_domain)
    call this%lap_v_buff%init(multi_domain)
    call this%lap_h_buff%init(multi_domain)

  end subroutine

  subroutine apply_diffusion(this, out, in, multi_domain)

    class (diffusion_operator_t), intent(inout) :: this
    class (stvec_t),              intent(inout) :: out
    class (stvec_t),              intent(inout) :: in
    type  (multi_domain_t),       intent(in)    :: multi_domain

    select type (out)
    class is (stvec_swe_t)
      select type(in)
      class is (stvec_swe_t)

        call calc_laplacian(this%lap_u_buff, in%u, multi_domain, this%coefs, this%diff2_op)
        call calc_laplacian(this%lap_v_buff, in%v, multi_domain, this%coefs, this%diff2_op)
        call calc_laplacian(this%lap_h_buff, in%h, multi_domain, this%coefs, this%diff2_op)
        call calc_laplacian(this%lap_u, this%lap_u_buff, multi_domain, this%coefs, this%diff2_op)
        call calc_laplacian(this%lap_v, this%lap_v_buff, multi_domain, this%coefs, this%diff2_op)
        call calc_laplacian(this%lap_h, this%lap_h_buff, multi_domain, this%coefs, this%diff2_op)

        call out%u%assign(-1.0_8, this%lap_u, multi_domain)
        call out%v%assign(-1.0_8, this%lap_v, multi_domain)
        call out%h%assign(-1.0_8, this%lap_h, multi_domain)

      class default
      end select
    class default
    end select

  end subroutine apply_diffusion

end module diffusion_operator_mod
