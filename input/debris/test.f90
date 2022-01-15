program main
  implicit none
  integer, parameter :: fi1 = 10
  integer, parameter :: fo1 = 20
  integer :: i
  integer :: count
  real(8) :: a, b, c
  real(8) :: X(3), V(3), K(6)
  real(8) :: pi, rad2deg

  pi = acos(-1.0d0)
  rad2deg = 180.0d0/pi
  
  open(fi1, file = "debris_new.dat", status = "old")
  open(fo1, file = "check.dat", status = "unknown")

  count = 0
  do
     read(fi1,*,end=99)
     read(fi1,*) a, b, c
     read(fi1,*) a, b, c
     read(fi1,*)
     count = count + 1
  end do
99 rewind(fi1)
  
  do i = 1, count
     read(fi1,*)
     read(fi1,*) X(1), X(2), X(3)
     read(fi1,*) V(1), V(2), V(3)
     read(fi1,*)
     call Car2Kep(X, V, K)
     write(fo1,*) K(1)*rad2deg, K(2)*rad2deg, K(3), K(4)*6378.142d0
  end do

  close(fi1)
  close(fo1)
contains
  subroutine Car2Kep(X, V, K)
    !*****************************************************************
    ! This subroutine convert Cartesian Elements to Keplerian Elements 
    !
    ! Input  : normalized state vector in Cartesian Coordinate (X,V)
    ! Output : normalized Keplerian Elements (K)
    !            K(1) : the orbit inclination (rad)
    !            K(2) : the longitude of the ascending node (rad)
    !            K(3) : the orbit eccentricity    
    !            K(4) : the semi-major radius (RE)                       
    !            K(5) : the the argument of periapsis (rad)
    !            K(6) : the true anomaly (rad)
    ! Ref    : Fundamentals of Astrodynamics and Applications
    !                        - Third Edition - p
    ! Author : Keisuke Akari
    ! Date   : September 1, 2016
    !*****************************************************************
    !---- args
    real(8), intent(in)  :: X(3),V(3)
    real(8), intent(out) :: K(6)
    !---- vars
    real(8), parameter   :: z(1:3) = (/0.0d0, 0.0d0, 1.0d0/)
    real(8) :: h(3)
    real(8) :: n(3)
    real(8) :: e(3), unit_e(3)
    real(8) :: r, xi, dotXV
    !---- body
    r = nolm(X)
    dotXV = dot(X,V)    
    h = cross(X,V)
    n = cross(z,h)
    e = cross(V,h) - unit_vector(X)
    unit_e = unit_vector(e)
    
    K(1) = acos(h(3)/nolm(h))
    K(2) = acos(n(1)/nolm(n))
    if(n(2) < 0.0d0)then
       K(2) = 2.0d0*pi - K(2)
    end if
    
    K(3) = nolm(e)
    K(4) = r/(2.0d0 - r*dot(V,V))
    K(5) = acos(dot(unit_vector(n), unit_e))
    if(e(3) < 0.0d0)then
       K(5) = 2.0d0*pi - K(5)
    end if
    
    K(6) = acos(dot(unit_e, unit_vector(X)))
    if(dotXV < 0.0d0)then
       K(6) = 2.0d0*pi - K(6)
    end if

    return
  end subroutine Car2Kep

  
  function dot(A,B)
    !********************************************
    ! This function calculates the inner product.
    !********************************************
    real(8), intent(in) :: A(3),B(3)
    real(8) :: dot

    dot = (A(1)*B(1)) + (A(2)*B(2)) + (A(3)*B(3))
    
    return
  end function dot
  

  function cross(A,B)
    !********************************************
    ! This function calculates the cross product.
    !********************************************
    real(8), intent(in) :: A(3), B(3)
    real(8) :: cross(3)

    cross(1) = A(2)*B(3) - A(3)*B(2)
    cross(2) = A(3)*B(1) - A(1)*B(3)
    cross(3) = A(1)*B(2) - A(2)*B(1)
    
    return
  end function cross


  function unit_vector(A)
    !*****************************************
    ! This function calculate the unit vector.
    !*****************************************
    real(8), intent(in) :: A(3)
    real(8) :: unit_vector(3)
    real(8) :: nolmA

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
    real(8), intent(in) :: A(3)
    real(8) :: nolm

    nolm = sqrt(A(1)*A(1) + A(2)*A(2) + A(3)*A(3))
    
    return
  end function nolm
end program main
