program main
  implicit none
  integer :: i
  real(8) :: r_min, r_max, r_inc
  integer :: n_debris
  real(8), allocatable :: r_debris(:)
  real(8), allocatable :: A_debris(:)
  real(8), allocatable :: M_debris(:)
  real(8) :: X_GEI_initial(3)
  real(8) :: V_GEI_initial(3)
  real(8) :: pi
  real(8), parameter :: rho_al = 2.7d3

  open(10,file='debris_new.dat',status='unknown')

  X_GEI_initial(1) = 1.125428377104179d0
  X_GEI_initial(2) = 0.0d0
  X_GEI_initial(3) = 0.0d0
  V_GEI_initial(1) = 0.0d0
  V_GEI_initial(2) = 0.0d0
  V_GEI_initial(3) = 0.942629591581051d0

  pi = acos(-1.0d0)
  write(6,*) "Input MINIMUM radius of debris (m)"
  read *, r_min
  write(6,*) "Input MAX radius of debris (m)"
  read *, r_max
  write(6,*) "What the size is increased by (m)?"
  read *, r_inc
  n_debris = int((r_max - r_min)/r_inc) + 1
  allocate(r_debris(n_debris))
  allocate(A_debris(n_debris))
  allocate(M_debris(n_debris))
  print*, n_debris

  r_debris(1) = r_min
  do i = 2, n_debris
     r_debris(i) = r_debris(1) + (i-1)*r_inc
  end do

  A_debris(:) = pi*r_debris(:)**2
  M_debris(:) = 4.0d0*pi*r_debris(:)**3*rho_al/3.0d0

  do i = 1, n_debris
     write(10,'(3ES24.15)') r_debris(i), A_debris(i), M_debris(i)
     write(10,'(3ES24.15)') X_GEI_initial(1), X_GEI_initial(2), X_GEI_initial(3)
     write(10,'(3ES24.15)') V_GEI_initial(1), V_GEI_initial(2), V_GEI_initial(3)
     write(10,*)
  end do
  
  close(10)
end program main

