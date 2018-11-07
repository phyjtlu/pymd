! File hello.f90
subroutine foo (a)
implicit none
integer a
print*, "Hello from Fortran foo!"
print*, "a=", a
end

subroutine bar (b)
implicit none
real(8),intent(in) :: b(3)
print*, "Hello from Fortran bar !"
print*, b
end
