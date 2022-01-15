program main
  implicit none
  !***********************************
  ! Setting
  !***********************************
  real(16), parameter :: Xinitial(3) &
       = (/1.125428377104179d0, 0.0d0, 0.0d0/)
  real(16), parameter :: Vinitial(3) &
       = (/0.0d0, 0.942629591581051d0, 0.0d0/)
  real(16), parameter :: Rd = 1.0d-6
  real(16), parameter :: Ad = 3.141592653589793d-12
  real(16), parameter :: Md = 1.1309733552923d-14
  integer, parameter :: d_inc   = 5
  integer, parameter :: d_Omega = 5
  integer, parameter :: d_nu    = 90
  !************************************
  integer, parameter :: fo1 = 10
  integer :: count
  real(16) :: pi
  real(16) :: deg2rad
  integer :: i_inc, i_Omega, i_nu
  integer :: ii, jj, kk
  real(16) :: A1(3,3), B2(3,3), C3(3,3), D(3,3), E(3,3)
  real(16) :: X(3), V(3)
  real(16) :: inc, Omega, nu

  count   = 0
  pi      = acos(-1.0d0)
  deg2rad = pi/180.0d0

  i_inc   = 180/d_inc
  i_Omega = 360/d_Omega
  i_nu    = 360/d_nu

  print*, "The Total Number of Debris =", (i_inc-1)*i_Omega*i_nu+2*4

  open(fo1, file = 'debris_new.dat', status = 'unknown')

  inc = 0.0d0
  Omega = 0.0d0
  call set_matrixX(inc, B2(:,:))
  call set_matrixZ(Omega, A1(:,:))
  do kk = 1, i_nu
     nu = d_nu*(kk-1)*deg2rad
     call set_matrixZ(nu, C3(:,:))

     call multiple_matrices(A1,B2,D)
     call multiple_matrices(D,C3,E)
     call multiple_matrix_vector(E,Xinitial,X)
     call multiple_matrix_vector(E,Vinitial,V)

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

     count = count + 1
     write(fo1,'(3ES24.15)') Rd, Ad, Md
     write(fo1,'(3ES24.15)') X(1), X(2), X(3)
     write(fo1,'(3ES24.15)') V(1), V(2), V(3)
     write(fo1,*)
  end do
  
  
  do ii = 1, i_inc-1
     inc = d_inc*ii*deg2rad
     call set_matrixX(inc, B2(:,:))
     do jj = 1, i_Omega
        Omega = d_Omega*(jj-1)*deg2rad
        call set_matrixZ(Omega, A1(:,:))
        do kk = 1, i_nu
           nu = d_nu*(kk-1)*deg2rad
           call set_matrixZ(nu, C3(:,:))

           call multiple_matrices(A1,B2,D)
           call multiple_matrices(D,C3,E)
           call multiple_matrix_vector(E,Xinitial,X)
           call multiple_matrix_vector(E,Vinitial,V)

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

           count = count + 1
           write(fo1,'(3ES24.15)') Rd, Ad, Md
           write(fo1,'(3ES24.15)') X(1), X(2), X(3)
           write(fo1,'(3ES24.15)') V(1), V(2), V(3)
           write(fo1,*)
        end do
     end do
  end do

  inc = 180.0d0*deg2rad
  Omega = 0.0d0
  call set_matrixX(inc, B2(:,:))
  call set_matrixZ(Omega, A1(:,:))
  do kk = 1, i_nu
     nu = d_nu*(kk-1)*deg2rad
     call set_matrixZ(nu, C3(:,:))

     call multiple_matrices(A1,B2,D)
     call multiple_matrices(D,C3,E)
     call multiple_matrix_vector(E,Xinitial,X)
     call multiple_matrix_vector(E,Vinitial,V)

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

     count = count + 1
     write(fo1,'(3ES24.15)') Rd, Ad, Md
     write(fo1,'(3ES24.15)') X(1), X(2), X(3)
     write(fo1,'(3ES24.15)') V(1), V(2), V(3)
     write(fo1,*)
  end do

  print*, count
  close(fo1)
contains
  subroutine set_matrixX(theta, A)
    real(16), intent(in)  :: theta
    real(16), intent(out) :: A(3,3)

    A(1,1) = 1.0d0
    A(2,1) = 0.0d0
    A(3,1) = 0.0d0

    A(1,2) = 0.0d0
    A(2,2) = cos(theta)
    A(3,2) = sin(theta)

    A(1,3) = 0.0d0
    A(2,3) = -sin(theta)
    A(3,3) = cos(theta)

    return
  end subroutine set_matrixX


  subroutine set_matrixZ(theta, A)
    real(16), intent(in)  :: theta
    real(16), intent(out) :: A(3,3)

    A(1,1) = cos(theta)
    A(2,1) = sin(theta)
    A(3,1) = 0.0d0

    A(1,2) = -sin(theta)
    A(2,2) = cos(theta)
    A(3,2) = 0.0d0

    A(1,3) = 0.0d0
    A(2,3) = 0.0d0
    A(3,3) = 1.0d0

    return
  end subroutine set_matrixZ


  subroutine multiple_matrices(A,B,C)
    real(16), intent(in)  :: A(3,3), B(3,3)
    real(16), intent(out) :: C(3,3)
    integer :: i, j, k

    C(:,:) = 0.0d0
    
    do i = 1, 3
       do j = 1, 3
          do k = 1, 3
             C(j,i) = C(j,i) + A(j,k)*B(k,i)
          end do
       end do
    end do

    return
  end subroutine multiple_matrices


  subroutine multiple_matrix_vector(A,Vin,Vout)
    real(16), intent(in) :: A(3,3), Vin(3)
    real(16), intent(out) :: Vout(3)
    integer :: i, j

    Vout(:) = 0.0d0
    
    do i = 1, 3
       do j = 1, 3
          Vout(i) = Vout(i) + A(i,j)*Vin(j)
       end do
    end do

    return
  end subroutine multiple_matrix_vector
end program main


















































































  
