module swe_advective_operator_mod
use stvec_swe_mod,             only: stvec_swe_t
use stvec_mod,                 only: stvec_t
use operator_mod,              only: operator_t
use differential_operator_mod, only: differential_operator_t
use field_mod,                 only: field_t
use mesh_mod,                only: mesh_t
use grad_mod,                  only: calc_grad
use div_mod,                   only: calc_div
use const_mod,                 only: Earth_grav, pcori
implicit none

  type, extends(operator_t) :: swe_advective_operator_t

    class(differential_operator_t) :: diff_opx
    class(differential_operator_t) :: diff_opy

  contains

    procedure :: apply => apply_swe_advective

  end type swe_advective_operator_t

contains

  subroutine apply_swe_advective(this, out, in, mesh)

    class (swe_advective_operator_t), intent(in)    :: this
    class (stvec_t),                  intent(inout) :: out
    class (stvec_t),                  intent(in)    :: in
    type  (mesh_t),                 intent(in)    :: mesh

    type(field_t) :: gx, gy, div, dudx, dudy, dvdx, dvdy, hu, hv
    integer       :: i, j

    call gx%init_on_mesh(mesh)
    call gy%init_on_mesh(mesh)
    call div%init_on_mesh(mesh)
    call dudx%init_on_mesh(mesh)
    call dudy%init_on_mesh(mesh)
    call dvdx%init_on_mesh(mesh)
    call dvdy%init_on_mesh(mesh)
    call hu%init_on_mesh(mesh)
    call hv%init_on_mesh(mesh)

    do i = 0, mesh%nx
      do j = 0, mesh%ny
        hu%f(i, j) = in%h%f(i, j) * in%u%f(i, j)
        hv%f(i, j) = in%h%f(i, j) * in%v%f(i, j)
      end do
    end do

    call calc_grad(gx, gy, in%h, mesh, this%diff_opx, this%diff_opy)
    call calc_div(div, hu, hv, mesh, this%diff_opx, this%diff_opy)
    call calc_grad(dudx, dudy, in%u, mesh, this%diff_opx, this%diff_opy)
    call calc_grad(dvdx, dvdy, in%v, mesh, this%diff_opx, this%diff_opy)

    do i = 0, mesh%nx
      do j = 0, mesh%ny
        out%u%f(i, j) = - Earth_grav * gx%f(i, j) + pcori * gy%f(i, j) - in%u%f(i, j) * dudx%f(i, j) + in%v%f(i, j) * dudy%f(i, j)
        out%v%f(i, j) = - Earth_grav * gx%f(i, j) - pcori * gy%f(i, j) - in%u%f(i, j) * dvdx%f(i, j) + in%v%f(i, j) * dvdy%f(i, j)
        out%h%f(i, j) = - div%f(i, j)
      end do
    end do

    ! out%u = sum(sum(mult_const(- Earth_grav, gx), mult_const(  pcori, in%v)), mult_const(-1.0_8, sum(mult(in%u, dudx), mult(in%v, dudy))))
    ! out%v = sum(sum(mult_const(- Earth_grav, gy), mult_const(- pcori, in%u)), mult_const(-1.0_8, sum(mult(in%u, dvdx), mult(in%v, dvdy))))
    ! out%h = mult_const(-1.0_8, div)

  end subroutine apply_swe_advective

end module swe_advective_operator_mod
