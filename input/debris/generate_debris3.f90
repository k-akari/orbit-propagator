program main
  implicit none
  !***********************************
  ! Setting
  !***********************************
  real(16), parameter :: Rd = 1.0d-6
  real(16), parameter :: Ad = 3.141592653589793d-12
  real(16), parameter :: Md = 1.1309733552923d-14
  real(16), parameter :: height = 800.0d0
  integer, parameter :: d_inc   = 5
  integer, parameter :: d_Omega = 5
  !************************************
  real(16) :: X0, V0
  real(8), parameter :: r0 = 6378.142d0
  real(8) :: vv0
  real(8), parameter :: mu_e = 3.986004356d14
  integer, parameter :: fo1 = 10
  integer :: count
  real(16) :: pi
  real(16) :: deg2rad
  integer :: i_inc, i_Omega, i_nu
  integer :: i, j
  real(16) :: A1(3,3), B2(3,3), C3(3,3), D(3,3), E(3,3)
  real(16) :: X(3), V(3)
  real(16) :: inc, Omega, nu

  X0 = (r0 + height)/r0
  V0 = 1.0d0/sqrt(X0)
  vv0 = sqrt(mu_e/(r0*1.0d3))

  count   = 0
  pi      = acos(-1.0d0)
  deg2rad = pi/180.0d0

  i_inc   = 180/d_inc
  i_Omega = 360/d_Omega
  
  print*, "The Total Number of Debris =", (i_inc-1)*i_Omega + 2

  open(fo1, file = 'debris_new.dat', status = 'unknown')

  X(1) = X0
  X(2) = 0.0d0
  X(3) = 0.0d0

  V(1) = 0.0d0
  V(2) = V0
  V(3) = 0.0d0

  write(fo1,'(3ES24.15)') Rd, Ad, Md
  write(fo1,'(3ES24.15)') X(1)*r0*1.0d3, X(2)*r0*1.0d3, X(3)*r0*1.0d3
  write(fo1,'(3ES24.15)') V(1)*vv0, V(2)*vv0, V(3)*vv0
  write(fo1,*)
  count = count + 1
  
  do i = 1, i_inc-1, 1
     do j = 1, i_omega, 1
        inc   = i*d_inc*deg2rad
        omega = (j-1)*d_omega*deg2rad
        X(1) = X0*cos(omega)
        X(2) = X0*sin(omega)
        X(3) = 0.0d0
        
        V(1) = -V0*cos(inc)*sin(omega)
        V(2) = V0*cos(inc)*cos(omega)
        V(3) = V0*sin(inc)
        write(fo1,'(3ES24.15)') Rd, Ad, Md
        write(fo1,'(3ES24.15)') X(1)*r0*1.0d3, X(2)*r0*1.0d3, X(3)*r0*1.0d3
        write(fo1,'(3ES24.15)') V(1)*vv0, V(2)*vv0, V(3)*vv0
        write(fo1,*)

        count = count + 1
     end do
  end do

  X(1) = X0
  X(2) = 0.0d0
  X(3) = 0.0d0

  V(1) = 0.0d0
  V(2) = -V0
  V(3) = 0.0d0

  write(fo1,'(3ES24.15)') Rd, Ad, Md
  write(fo1,'(3ES24.15)') X(1)*r0*1.0d3, X(2)*r0*1.0d3, X(3)*r0*1.0d3
  write(fo1,'(3ES24.15)') V(1)*vv0, V(2)*vv0, V(3)*vv0
  write(fo1,*)
  count = count + 1
  
  print*, count
  close(fo1)
end program main
