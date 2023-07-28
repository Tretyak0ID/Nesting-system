program sbp_operators_test
use sbp_differential_operator_mod, only: sbp21_t, sbp42_t
use field_mod,         only: field_t
use domain_mod,        only: domain_t
use const_mod,         only: pi
use read_write_mod,    only: write_field

  type(field_t)  :: in_field, out_field, cos_field
  type(domain_t) :: domain
  type(sbp42_t)  :: sbp42
  integer        :: i, j

  N = 64
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))

  !-------------------------------------------------
N = 128
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))

  !-------------------------------------------------
  N = 256
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))

  !-------------------------------------------------
  N = 512
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))

  !-------------------------------------------------
  N = 1024
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))

  !-------------------------------------------------
  N = 2048
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))

  !-------------------------------------------------
  N = 4096
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))

  !-------------------------------------------------
  N = 4096 * 2
  call domain%init(0.0_8, 2.0_8 * pi, 0, N, 0.0_8, 2.0_8 * pi, 0, N)
  call in_field%init(0, N, 0, N)
  call out_field%init(0, N, 0, N)
  call cos_field%init(0, N, 0, N)

  do j = 0, domain%ny
    do i = 0, domain%nx
      in_field%f(i, j) = sin(domain%x(i))
      cos_field%f(i, j) = cos(domain%x(i))
    end do
  end do

  call sbp42%apply(out_field, in_field, domain, 'x')

  print *, maxval(out_field%f(3:N-2, 3:N-2) - cos_field%f(3:N-2, 3:N-2))



end program sbp_operators_test
