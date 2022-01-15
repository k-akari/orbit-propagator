module subprog_mod
  implicit none
contains
  subroutine interpolation(Xi1,Xi2,Xi3,T,T2,interval,Xo)
    !************************************************************
    ! This subroutine predict the true value from 3 tabular data.
    ! Use this subroutine to determine the celestial positions.
    !
    ! Ref    : Jean Meeus, Astronomical Algorithms
    !              - Second Edition - p23
    ! Author : Keisuke Akari
    ! Date   : November 3, 2016
    !************************************************************
    !---- args
    real(8), intent(in)  :: Xi1(3), Xi2(3), Xi3(3)
    real(8), intent(in)  :: T, T2, interval
    real(8), intent(out) :: Xo(3)
    !---- vars
    real(8) :: n, A(3), B(3), C(3)
    !---- body
    n = (T - T2)/interval
    A = Xi2 - Xi1
    B = Xi3 - Xi2
    C = B - A
    Xo = Xi2 + 0.5d0*n*(A + B + n*C)

    return
  end subroutine interpolation

  
  subroutine angle(A,B,theta)
    !**********************************************************
    ! This subroutine calculates the angle between two vectors.
    !**********************************************************
    !---- args
    real(8), intent(in)  :: A(3),B(3)
    real(8), intent(out) :: theta
    !---- vars
    real(8) :: nolmA
    real(8) :: nolmB
    real(8) :: AdotB
    real(8) :: cos_theta
    !---- body
    nolmA = sqrt(A(1)**2 + A(2)**2 + A(3)**2)
    nolmB = sqrt(B(1)**2 + B(2)**2 + B(3)**2)
    
    if(nolmA == 0.0d0 .or. nolmB == 0.0d0)then
       print*, 'Error : cannot calculate an angle between two vectors'
       stop
    end if
    
    AdotB = (A(1)*B(1)) + (A(2)*B(2)) + (A(3)*B(3))

    cos_theta = AdotB/(nolmA*nolmB)
    if(cos_theta >= 1.0d0)then
       cos_theta = 0.999999999999999d0
    else if(cos_theta <= -1.0d0)then
       cos_theta = -0.999999999999999d0
    end if
    
    theta = acos(cos_theta)

    return
  end subroutine angle
  
  
  subroutine dot_product(A, B, AdotB)
    !---- args
    real(8), intent(in)  :: A(3), B(3)
    real(8), intent(out) :: AdotB
    !---- vars
    !---- body
    AdotB = (A(1)*B(1)) + (A(2)*B(2)) + (A(3)*B(3))

    return
  end subroutine dot_product
  
  
  subroutine cross_product(A, B, AcrossB)
    !---- args
    real(8), intent(in)  :: A(3), B(3)
    real(8), intent(out) :: AcrossB(3)
    !---- vars
    !---- body
    AcrossB(1) = A(2)*B(3) - A(3)*B(2)
    AcrossB(2) = A(3)*B(1) - A(1)*B(3)
    AcrossB(3) = A(1)*B(2) - A(2)*B(1)

    return
  end subroutine cross_product


  subroutine calc_unit_vector(A, B)
    !---- args
    real(8), intent(in) :: A(3)
    real(8), intent(out) :: B(3)
    !---- vars
    real(8) :: nolmA
    !---- body
    nolmA = sqrt(A(1)**2 + A(2)**2 + A(3)**2)
    B = A / nolmA

    return
  end subroutine calc_unit_vector
  
  
  function dot(A,B)
    !********************************************
    ! This function calculates the inner product.
    !********************************************
    !---- args
    real(8), intent(in) :: A(3),B(3)
    !---- vars
    real(8) :: dot
    !---- body
    dot = (A(1)*B(1)) + (A(2)*B(2)) + (A(3)*B(3))
    
    return
  end function dot
  

  function cross(A,B)
    !********************************************
    ! This function calculates the cross product.
    !********************************************
    !---- args
    real(8), intent(in) :: A(3), B(3)
    !---- vars
    real(8) :: cross(3)
    !---- body
    cross(1) = A(2)*B(3) - A(3)*B(2)
    cross(2) = A(3)*B(1) - A(1)*B(3)
    cross(3) = A(1)*B(2) - A(2)*B(1)
    
    return
  end function cross


  function unit_vector(A)
    !*****************************************
    ! This function calculate the unit vector.
    !*****************************************
    !---- args
    real(8), intent(in) :: A(3)
    !---- vars
    real(8) :: unit_vector(3)
    real(8) :: nolmA
    !---- body
    nolmA = sqrt(A(1)**2 + A(2)**2 + A(3)**2)
    if(nolmA == 0.0d0)then
       unit_vector = 0.0d0
    else
       unit_vector = A / nolmA
    end if
    
    return
  end function unit_vector


  function nolm(A)
    !************************************************
    ! This function returns the nolm of input vector.
    !***********************************************
    !---- args
    real(8), intent(in) :: A(3)
    !---- vars
    real(8) :: nolm
    !---- body
    nolm = sqrt(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
    
    return
  end function nolm


  function signal(A)
    !***********************************************
    ! This function returns the sign of input value.
    !***********************************************
    !---- args
    real(8), intent(in) :: A
    !---- vars
    real(8) :: signal
    !---- body
    if(A >= 0.0d0)then
       signal = 1.0d0
    else if(A < 0.0d0)then
       signal = -1.0d0
    end if

    return
  end function signal


  recursive real(16) function factorial(n) result(ret)
    !*************************************************************
    ! Warning :
    ! When calculating factorial values, set double precision type
    ! as function type. Otherwise, incorrect values are returned 
    ! if argument is large enough. This error is due to the inner
    ! expression of the computer.
    !*************************************************************
    !---- args
    integer, intent(in) :: n
    !---- vars
    !---- body
    if(n >= 1)then
       ret = n * factorial(n-1)
    else if(N == 0)then
       ret = 1
    end if
    
    return
  end function factorial


  function permutation(n,k)
    !---- args
    integer, intent(in) :: n, k
    !---- vars
    real(16)            :: permutation
    integer             :: i
    !---- body
    permutation = 1.0d0
    do i = k+1, n
       permutation = permutation*i
    end do

    return
  end function permutation
end module subprog_mod
