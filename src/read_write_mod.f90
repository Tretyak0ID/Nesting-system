module read_write_mod
use field_mod,  only : field_t
use domain_mod, only : domain_t
implicit none

contains

  subroutine write_field(f, domain, file_name, rec_num)

    type(field_t),               intent(inout) :: f
    type(domain_t),              intent(in)    :: domain
    character(len=*),            intent(in)    :: file_name
    integer(kind=4), optional,   intent(in)    :: rec_num

    integer(kind=4) :: irec, reclen

    irec = 1
    if (present(rec_num)) irec = rec_num

    ! number of points to write
    reclen = (domain%nx+1)*(domain%ny+1)

    open(unit=111, file = trim(file_name), access="direct", recl = reclen)

    write(111, rec = irec) real(f%f,4)

    close(111)

  end subroutine write_field

end module read_write_mod
