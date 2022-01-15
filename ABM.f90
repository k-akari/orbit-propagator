module ABM_mod
  use global_mod, only : &
       set_clear, &
       X, V, K, XdeGEI, VdeGEI, Atotal, &
       limit_hGEI, i_max, i_step, dt, &
       N_ABM, N_debris, &
       Cd, Ad, Md, &
       keisuu_airdrag, keisuu_ph
  use accel_mod, only : &
       renew_acceleration1, &
       renew_acceleration2
  use subprog_mod, only : &
       factorial
  use time_mod, only : &
       renew_clock, &
       set_initial_time, &
       set_time1
  use air_drag_mod, only : &
       set_airdrag_constants
  use sun_mod, only : &
       set_solar_radiation_pressure_constants
  use conversion_coordinates_mod, only : &
       Car2Kep
  use output_mod, only : &
       open_file, &
       close_file, &
       output
  use input_mod, only : &
       set_debris
  implicit none
  real(16), save, allocatable :: AB(:)
  real(16), save, allocatable :: AM(:)
  real(8),  save, allocatable :: VV(:,:)
  real(8),  save, allocatable :: AA(:,:)
contains
  !**********************************************************
  ! This module program can solve the equation of motion
  ! by using predictor-corrector method of Nth-order
  ! Adams-Bashforth-Moulton method.
  !
  ! An order of this calculation scheme can be varied.
  ! Coefficients required for this scheme of Nth order should
  ! be calculated before the numerical integration starts.
  !
  ! Nth order Adams-Bashforth-Moulton method requires the
  ! initial values plus N-1 previous steps to predict the
  ! next step values. Therefore, another numerical integration
  ! scheme is required for starting. In this program, 8th order
  ! Runge-Kutta method is employed to prepare the past values.
  ! 
  ! Author : Keisuke Akari
  ! Date   : November 23, 2016
  !**********************************************************
  subroutine start_orbital_propagation(ID_debris)
    !---- args
    integer, intent(in) :: ID_debris
    !---- vars
    !---- body    
    call open_file(ID_debris)
    call set_debris(ID_debris)
    call set_clear
    call set_initial_time
    call set_time1
    call set_airdrag_constants(Cd, Ad, Md, keisuu_airdrag)
    call set_solar_radiation_pressure_constants(Ad, Md, keisuu_ph)
    call start_Adams_Bashforth_Moulton
    call close_file
    
    print*, ID_debris, '/', N_debris

    return
  end subroutine start_orbital_propagation

  
  subroutine start_Adams_Bashforth_Moulton
    !---- args
    !---- vars
    integer :: i
    integer :: j
    !---- body
    call start_Runge_Kutta(N_ABM)

    do i = 1, i_step-N_ABM, 1
       call Adams_Bashforth_Moulton(N_ABM)
       if(sqrt(X(1)**2 + X(2)**2 + X(3)**2) <= limit_hGEI)then
          go to 10
       end if
    end do
    call output

    do i = 1, i_max-1, i_step
       do j = 1, i_step, 1
          call Adams_Bashforth_Moulton(N_ABM)
          if(sqrt(X(1)**2 + X(2)**2 + X(3)**2) <= limit_hGEI)then
             go to 10
          end if
       end do
       call output
    end do

10  call output
    
    return
  end subroutine start_Adams_Bashforth_Moulton

  
  subroutine start_Runge_Kutta(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    integer :: i
    integer :: j
    !---- body
    VV(:,:) = 0.0d0
    AA(:,:) = 0.0d0
    
    do i = 1, N, 1
       call Runge_Kutta8
       
       do j = 0, N-2, 1
          VV(:,j) = VV(:,j+1)
          AA(:,j) = AA(:,j+1)
       end do

       XdeGEI(:) = X(:)
       VdeGEI(:) = V(:)
       call renew_acceleration2

       VV(:,N-1) = V(:)
       AA(:,N-1) = Atotal(:)
    end do

    return
  end subroutine start_Runge_Kutta
    
       
  subroutine Adams_Bashforth_Moulton(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    integer  :: i
    real(16) :: V_AB(3), A_AB(3), V_AM(3), A_AM(3)
    real(8)  :: V_AM1(3), V_AM2(3)
    real(8)  :: A_AM1(3), A_AM2(3)
    !---- body
    V_AB(:)  = 0.0d0
    A_AB(:)  = 0.0d0
    V_AM(:)  = 0.0d0
    A_AM(:)  = 0.0d0
    V_AM1(:) = 0.0d0
    A_AM1(:) = 0.0d0
    V_AM2(:) = 0.0d0
    A_AM2(:) = 0.0d0
    
    call renew_clock(dt)

    !*********************************
    ! Nth order Adams-Bashforth Method
    !*********************************
    do i = 0, N-1, 1
       V_AB(:) = V_AB(:) + VV(:,i)*AB(i)
       A_AB(:) = A_AB(:) + AA(:,i)*AB(i)
    end do

    XdeGEI(:) = X(:) + dt*dble(V_AB(:))
    VdeGEI(:) = V(:) + dt*dble(A_AB(:))
    call renew_acceleration1

    !*******************************
    ! Nth order Adams-Moulton Method
    !*******************************
    do i = 0, N-2, 1
       V_AM(:) = V_AM(:) + VV(:,i+1)*AM(i)
       A_AM(:) = A_AM(:) + AA(:,i+1)*AM(i)
    end do
    V_AM1(:) = dble(V_AM(:) + AM(N-1)*VdeGEI(:))
    A_AM1(:) = dble(A_AM(:) + AM(N-1)*Atotal(:))

    XdeGEI(:) = X(:) + dt*V_AM1(:)
    VdeGEI(:) = V(:) + dt*A_AM1(:)
    call renew_acceleration2

    V_AM2(:) = dble(V_AM(:) + AM(N-1)*VdeGEI(:))
    A_AM2(:) = dble(A_AM(:) + AM(N-1)*Atotal(:))

    X(:) = X(:) + dt*V_AM2(:)
    V(:) = V(:) + dt*A_AM2(:)
    call renew_acceleration2

    !**************************************************
    ! Renew the past values of the positions, velocity,
    ! and accelerations.
    !**************************************************
    do i = 0, N-2, 1
       VV(:,i) = VV(:,i+1)
       AA(:,i) = AA(:,i+1)
    end do
    VV(:,N-1) = V(:)
    AA(:,N-1) = Atotal(:)

    return
  end subroutine Adams_Bashforth_Moulton


  subroutine Runge_Kutta8
    !---- args
    !---- vars
    real(8) :: K1(3),K2(3),K3(3),K4(3),K5(3)
    real(8) :: K6(3),K7(3),K8(3),K9(3),K10(3),K11(3)
    real(8) :: L1(3),L2(3),L3(3),L4(3),L5(3)
    real(8) :: L6(3),L7(3),L8(3),L9(3),L10(3),L11(3)
    !---- body
    XdeGEI = X
    VdeGEI = V
    call renew_acceleration1
    K1 = VdeGEI*dt
    L1 = Atotal*dt

    call renew_clock(0.5d0*dt)
    XdeGEI = X + 0.5d0*K1
    VdeGEI = V + 0.5d0*L1
    call renew_acceleration1
    K2 = VdeGEI*dt
    L2 = Atotal*dt

    XdeGEI = X + 0.25d0*K1 + 0.25d0*K2
    VdeGEI = V + 0.25d0*L1 + 0.25d0*L2
    call renew_acceleration2
    K3 = VdeGEI*dt
    L3 = Atotal*dt
    
    call renew_clock(-sqrt(21.0d0)/14.0d0*dt)
    XdeGEI = X + K1/7.0d0 &
         + (-7.0d0 + 3.0d0*sqrt(21.0d0))/98.0d0*K2 &
         + (21.0d0 - 5.0d0*sqrt(21.0d0))/49.0d0*K3
    VdeGEI = V + L1/7.0d0 &
         + (-7.0d0 + 3.0d0*sqrt(21.0d0))/98.0d0*L2 &
         + (21.0d0 - 5.0d0*sqrt(21.0d0))/49.0d0*L3
    call renew_acceleration1
    K4 = VdeGEI*dt
    L4 = Atotal*dt

    XdeGEI = X + (11.0d0 - sqrt(21.0d0))/84.0d0*K1 &
         + (18.0d0 - 4.0d0*sqrt(21.0d0))/63.0d0*K3 &
         + (21.0d0 + sqrt(21.0d0))/252.0d0*K4
    VdeGEI = V + (11.0d0 - sqrt(21.0d0))/84.0d0*L1 &
         + (18.0d0 - 4.0d0*sqrt(21.0d0))/63.0d0*L3 &
         + (21.0d0 + sqrt(21.0d0))/252.0d0*L4
    call renew_acceleration2
    K5 = VdeGEI*dt
    L5 = Atotal*dt

    call renew_clock(sqrt(21.0d0)/14.0d0*dt)
    XdeGEI = X + (5.0d0 - sqrt(21.0d0))/48.0d0*K1 &
         + (9.0d0 - sqrt(21.0d0))/36.0d0*K3 &
         + (-231.0d0 - 14.0d0*sqrt(21.0d0))/360.0d0*K4 &
         + (63.0d0 + 7.0d0*sqrt(21.0d0))/80.0d0*K5
    VdeGEI = V + (5.0d0 - sqrt(21.0d0))/48.0d0*L1 &
         + (9.0d0 - sqrt(21.0d0))/36.0d0*L3 &
         + (-231.0d0 - 14.0d0*sqrt(21.0d0))/360.0d0*L4 &
         + (63.0d0 + 7.0d0*sqrt(21.0d0))/80.0d0*L5  
    call renew_acceleration1
    K6 = VdeGEI*dt
    L6 = Atotal*dt

    call renew_clock(sqrt(21.0d0)/14.0d0*dt)
    XdeGEI = X + (10.0d0 + sqrt(21.0d0))/42.0d0*K1 &
         + (-432.0d0 - 92.0d0*sqrt(21.0d0))/315.0d0*K3 &
         + (633.0d0 + 145.0d0*sqrt(21.0d0))/90.0d0*K4 &
         + (-504.0d0 - 115.0d0*sqrt(21.0d0))/70.0d0*K5 &
         + (63.0d0 + 13.0d0*sqrt(21.0d0))/35.0d0*K6
    VdeGEI = V + (10.0d0 + sqrt(21.0d0))/42.0d0*L1 &
         + (-432.0d0 - 92.0d0*sqrt(21.0d0))/315.0d0*L3 &
         + (633.0d0 + 145.0d0*sqrt(21.0d0))/90.0d0*L4 &
         + (-504.0d0 - 115.0d0*sqrt(21.0d0))/70.0d0*L5 &
         + (63.0d0 + 13.0d0*sqrt(21.0d0))/35.0d0*L6         
    call renew_acceleration1
    K7 = VdeGEI*dt
    L7 = Atotal*dt

    XdeGEI = X + K1/14.0d0 &
         + (14.0d0 + 3.0d0*sqrt(21.0d0))/126.0d0*K5 &
         + (13.0d0 + 3.0d0*sqrt(21.0d0))/63.0d0*K6 + K7/9.0d0
    VdeGEI = V + L1/14.0d0 &
         + (14.0d0 + 3.0d0*sqrt(21.0d0))/126.0d0*L5 &
         + (13.0d0 + 3.0d0*sqrt(21.0d0))/63.0d0*L6 + L7/9.0d0
    call renew_acceleration2
    K8 = VdeGEI*dt
    L8 = Atotal*dt

    call renew_clock(-sqrt(21.0d0)/14.0d0*dt)
    XdeGEI = X + K1/32.0d0 &
         + (91.0d0 + 21.0d0*sqrt(21.0d0))/576.0d0*K5 &
         + 11.0d0/72.0d0*K6 &
         + (-385.0d0 + 75.0d0*sqrt(21.0d0))/1152.0d0*K7 &
         + (63.0d0 - 13.0d0*sqrt(21.0d0))/128.0d0*K8
    VdeGEI = V + L1/32.0d0 &
         + (91.0d0 + 21.0d0*sqrt(21.0d0))/576.0d0*L5 &
         + 11.0d0/72.0d0*L6 &
         + (-385.0d0 + 75.0d0*sqrt(21.0d0))/1152.0d0*L7 &
         + (63.0d0 - 13.0d0*sqrt(21.0d0))/128.0d0*L8    
    call renew_acceleration1
    K9 = VdeGEI*dt
    L9 = Atotal*dt

    call renew_clock(-sqrt(21.0d0)/14.0d0*dt)
    XdeGEI = X + K1/14.0d0 + K5/9.0d0 &
         + (-733.0d0 + 147.0d0*sqrt(21.0d0))/2205.0d0*K6 &
         + (515.0d0 - 111.0d0*sqrt(21.0d0))/504.0d0*K7 &
         + (-51.0d0 + 11.0d0*sqrt(21.0d0))/56.0d0*K8 &
         + (132.0d0 - 28.0d0*sqrt(21.0d0))/245.0d0*K9
    VdeGEI = V + L1/14.0d0 + L5/9.0d0 &
         + (-733.0d0 + 147.0d0*sqrt(21.0d0))/2205.0d0*L6 &
         + (515.0d0 - 111.0d0*sqrt(21.0d0))/504.0d0*L7 &
         + (-51.0d0 + 11.0d0*sqrt(21.0d0))/56.0d0*L8 &
         + (132.0d0 - 28.0d0*sqrt(21.0d0))/245.0d0*L9
    call renew_acceleration1
    K10 = VdeGEI*dt
    L10 = Atotal*dt

    call renew_clock((0.5d0 + sqrt(21.0d0)/14.0d0)*dt)
    XdeGEI = X + (-6.0d0 - sqrt(21.0d0))*7.0d0/18.0d0*K5 &
         + (-18.0d0 - 28.0d0*sqrt(21.0d0))/45.0d0*K6 &
         + (-273.0d0 + 53.0d0*sqrt(21.0d0))/72.0d0*K7 &
         + (301.0d0 - 53.0d0*sqrt(21.0d0))/72.0d0*K8 &
         + (1.0d0 + sqrt(21.0d0))*28.0d0/45.0d0*K9 &
         + (7.0d0 + sqrt(21.0d0))*7.0d0/18.0d0*K10
    VdeGEI = V + (-6.0d0 - sqrt(21.0d0))*7.0d0/18.0d0*L5 &
         + (-18.0d0 - 28.0d0*sqrt(21.0d0))/45.0d0*L6 &
         + (-273.0d0 + 53.0d0*sqrt(21.0d0))/72.0d0*L7 &
         + (301.0d0 - 53.0d0*sqrt(21.0d0))/72.0d0*L8 &
         + (1.0d0 + sqrt(21.0d0))*28.0d0/45.0d0*L9 &
         + (7.0d0 + sqrt(21.0d0))*7.0d0/18.0d0*L10
    call renew_acceleration1
    K11 = VdeGEI*dt
    L11 = Atotal*dt

    X = X + (9.0d0*K1 + 49.0d0*K8 + 64.0d0*K9 + &
         49.0d0*K10 + 9.0d0*K11)/180.0d0
    V = V + (9.0d0*L1 + 49.0d0*L8 + 64.0d0*L9 + &
         49.0d0*L10 + 9.0d0*L11)/180.0d0

    return
  end subroutine Runge_Kutta8


  subroutine set_ABM_coefficients(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    !---- body
    call set_array(N)
    call set_AB_coefficients(N)
    call set_AM_coefficients(N)

    return
  end subroutine set_ABM_coefficients


  subroutine set_array(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    !---- body
    allocate(VV(1:3,0:N-1))
    allocate(AA(1:3,0:N-1))

    VV(:,:) = 0.0d0
    AA(:,:) = 0.0d0

    return
  end subroutine set_array

  
  subroutine set_AB_coefficients(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    integer :: i, j, m
    !---- body
    real(16), allocatable :: AB2(:,:,:)
    real(16), allocatable :: AB3(:,:)
    real(16), allocatable :: AB4(:)
    real(16), allocatable :: AB5(:)
    real(16), allocatable :: AB6(:)
    real(16), allocatable :: AB7(:)

    allocate(AB2(0:N-1,0:N-1,2:N))
    allocate(AB3(0:N,0:N-1))
    allocate(AB4(0:N-1))
    allocate(AB5(0:N-1))    
    allocate(AB6(0:N-1))
    allocate(AB7(0:N-1))
    allocate(AB(0:N-1))

    AB2(:,:,:) = 0.0d0
    AB3(:,:)   = 0.0d0
    AB4(:)     = 0.0d0
    AB5(:)     = 0.0d0
    AB6(:)     = 0.0d0
    AB7(:)     = 0.0d0
    AB(:)      = 0.0d0

    AB2(0,0,2) = 0.0d0
    AB2(1,0,2) = 1.0d0
    AB2(0,1,2) = 1.0d0
    AB2(1,1,2) = 1.0d0
    
    do j = 3, N, 1
       AB6(:) = AB2(:,0,j-1)
       call multiply_polynomials(AB6, j-2, AB7)
       AB2(:,0,j) = AB7(:)
       do i = 1, j-1, 1
          AB6(:) = AB2(:,i-1,j-1)
          call multiply_polynomials(AB6, j-1, AB7)
          AB2(:,i,j) = AB7(:)
       end do
    end do

    do m = 0, N-1, 1
       call integrate_polynomial(AB2(:,m,N), AB3(:,m))
    end do

    do m = 0, N-1, 1
       AB4(m) = factorial(m) * factorial(N-m-1) * (-1)**(N+m+1)
    end do

    do m = 0, N-1, 1
       AB5(m) = 0.0d0
       do i = 0, N, 1
          AB5(m) = AB5(m) + AB3(i,m)
       end do
    end do

    print*, ""
    print*,'     **********************************************'
    print*,N,'th order Adams-Bashforth Coefficients'
    print*,'     **********************************************'
    do m = 0, N-1, 1
       AB(m) = AB5(m) / AB4(m)
       print*,'    ',AB(m)
    end do
    print*, ""

    deallocate(AB2)
    deallocate(AB3)
    deallocate(AB4)
    deallocate(AB5)
    deallocate(AB6)
    deallocate(AB7)

    return
  end subroutine set_AB_coefficients
  

  subroutine set_AM_coefficients(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    integer :: i, j, m
    real(16), allocatable :: AM2(:,:,:)
    real(16), allocatable :: AM3(:,:)
    real(16), allocatable :: AM4(:)
    real(16), allocatable :: AM5(:)
    real(16), allocatable :: AM6(:)
    real(16), allocatable :: AM7(:)
    !---- body
    allocate(AM2(0:N-1,0:N-1,2:N))
    allocate(AM3(0:N,0:N-1))
    allocate(AM4(0:N-1))
    allocate(AM5(0:N-1))
    allocate(AM6(0:N-1))
    allocate(AM7(0:N-1))
    allocate(AM(0:N-1))

    AM2(:,:,:) = 0.0d0
    AM3(:,:)   = 0.0d0
    AM4(:)     = 0.0d0
    AM5(:)     = 0.0d0
    AM6(:)     = 0.0d0
    AM7(:)     = 0.0d0
    AM(:)      = 0.0d0
    
    AM2(0,0,2) = -1.0d0
    AM2(1,0,2) = 1.0d0
    AM2(0,1,2) = 0.0d0
    AM2(1,1,2) = 1.0d0

    do j = 3, N, 1
       AM6(:) = AM2(:,0,j-1)
       call multiply_polynomials(AM6, j-3, AM7)
       AM2(:,0,j) = AM7(:)
       do i = 1, j-1, 1
          AM6(:) = AM2(:,i-1,j-1)
          call multiply_polynomials(AM6, j-2, AM7)
          AM2(:,i,j) = AM7(:)
       end do
    end do

    do m = 0, N-1, 1
       call integrate_polynomial(AM2(:,m,N), AM3(:,m))
    end do
    
    do m = 0, N-1, 1
       AM4(m) = factorial(N-m-1)*factorial(m)*(-1.0d0)**(N-m-1)
    end do

    do m = 0, N-1, 1
       AM5(m) = 0.0d0
       do i = 0, N, 1
          AM5(m) = AM5(m) + AM3(i,m)
       end do
    end do

    print*, ""
    print*,'     **********************************************'
    print*, N,'th order Adams-Moulton Coefficients'
    print*,'     **********************************************'
    do m = 0, N-1, 1
       AM(m) = AM5(m) / AM4(m)
       print*,'    ',AM(m)
    end do
    print*, ""
    
    deallocate(AM2)
    deallocate(AM3)
    deallocate(AM4)
    deallocate(AM5)
    deallocate(AM6)
    deallocate(AM7)

    return
  end subroutine set_AM_coefficients


  subroutine multiply_polynomials(A,M,C)
    !******************************************************
    ! Input  : A = A(0) + A(1)*x + ... + A(n-1)*x^(n-1)
    !          the number of components of input array (N)
    !          multiplcation of the polynomial (x + M)
    ! output : C = C(0) + C(1)*x + ... + C(n-1)*x^(n-1)
    !******************************************************
    !---- args
    integer, intent(in)  :: M
    real(16), intent(in)  :: A(0:)
    real(16), intent(out) :: C(0:)
    !---- vars
    integer :: i, N
    !---- body
    N = size(A)

    C(0) = dble(M)*A(0)
    do i = 1, N-1
       C(i) = A(i-1) + dble(M)*A(i)
    end do

    return
  end subroutine multiply_polynomials


  subroutine integrate_polynomial(A,B)
    !**************************************************
    ! Input  : A = A(0) + A(1)*x + ... + A(n)*x^n
    ! Output : B = B(0) + B(1)*x + ... + B(n+1)*x^(n+1)
    !**************************************************
    !---- args
    real(16), intent(in)  :: A(0:)
    real(16), intent(out) :: B(0:)
    !---- vars
    integer :: i
    !---- body
    B(0) = 0.0d0
    
    do i = 1, size(A)
       B(i) = A(i-1)/dble(i)
    end do
  end subroutine integrate_polynomial
end module ABM_Mod
