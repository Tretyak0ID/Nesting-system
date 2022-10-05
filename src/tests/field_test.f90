program field_test
use field_mod, only: field_t

  type(field_t) :: field
  call field%init(-2, 2, -2, 2)

  field%f(1,-1) = 10.0

  print *, field%f

end program field_test
