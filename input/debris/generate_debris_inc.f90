program main
  implicit none
  !***********************************
  ! Setting
  !***********************************
  real(16), parameter :: Rd = 1.0d-6
  real(16), parameter :: Ad = 3.141592653589793d-12
  real(16), parameter :: Md = 1.1309733552923d-14
  real(16), parameter :: height = 800.0d0
  integer, parameter  :: d_inc   = 30
  real(16), parameter :: Omega = 120.0d0
  !************************************
  real(16) :: X0, V0
  integer, parameter :: fo1 = 10
  integer :: count
  real(16) :: pi
  real(16) :: deg2rad
  real(16) :: x_omega
  integer :: i_inc, i_Omega, i_nu
  integer :: i, j
  real(16) :: A1(3,3), B2(3,3), C3(3,3), D(3,3), E(3,3)
  real(16) :: X(3), V(3)
  real(16) :: inc

  X0 = (6378.142d0 + height)/6378.142d0
  V0 = 1.0d0/sqrt(X0)

  count   = 0
  pi      = acos(-1.0d0)
  deg2rad = pi/180.0d0

  x_omega = omega*deg2rad
  i_inc   = 180/d_inc
  
  print*, "The Total Number of Debris =", i_inc + 1

  open(fo1, file = 'debris_new.dat', status = 'unknown')

  X(1) = X0
  X(2) = 0.0d0
  X(3) = 0.0d0
  
  V(1) = 0.0d0
  V(2) = V0
  V(3) = 0.0d0
     
  write(fo1,'(3ES24.15)') Rd, Ad, Md
  write(fo1,'(3ES24.15)') X(1), X(2), X(3)
  write(fo1,'(3ES24.15)') V(1), V(2), V(3)
  write(fo1,*)
  count = count + 1
  
  do i = 1, i_inc-1, 1
     inc   = i*d_inc*deg2rad
     
     X(1) = X0*cos(x_omega)
     X(2) = X0*sin(x_omega)
     X(3) = 0.0d0
     
     V(1) = -V0*cos(inc)*sin(x_omega)
     V(2) = V0*cos(inc)*cos(x_omega)
     V(3) = V0*sin(inc)

     if(abs(X(1)) < 1.0d-15)then
        X(1) = 0.0d0
     end if
     if(abs(X(2)) < 1.0d-15)then
        X(2) = 0.0d0
     end if
     if(abs(X(3)) < 1.0d-15)then
        X(3) = 0.0d0
     end if
     if(abs(V(1)) < 1.0d-15)then
        V(1) = 0.0d0
     end if
     if(abs(V(2)) < 1.0d-15)then
        V(2) = 0.0d0
     end if
     if(abs(V(3)) < 1.0d-15)then
        V(3) = 0.0d0
     end if
     
     write(fo1,'(3ES24.15)') Rd, Ad, Md
     write(fo1,'(3ES24.15)') X(1), X(2), X(3)
     write(fo1,'(3ES24.15)') V(1), V(2), V(3)
     write(fo1,*)
     
     count = count + 1
  end do
  
  X(1) = X0
  X(2) = 0.0d0
  X(3) = 0.0d0
  
  V(1) = 0.0d0
  V(2) = -V0
  V(3) = 0.0d0
  
  write(fo1,'(3ES24.15)') Rd, Ad, Md
  write(fo1,'(3ES24.15)') X(1), X(2), X(3)
  write(fo1,'(3ES24.15)') V(1), V(2), V(3)
  write(fo1,*)
  count = count + 1
  
  print*, count
  close(fo1)
end program main
