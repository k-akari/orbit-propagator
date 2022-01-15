module geopotential_mod
  use global_mod, only : &
       C_geo, S_geo, fi1
  use subprog_mod, only : permutation
  implicit none
  real(8), save, allocatable :: Legendre(:,:)
  real(8), save, allocatable :: C_normalized(:,:)
  real(8), save, allocatable :: S_normalized(:,:)
contains
  subroutine Attraction_of_the_earth(Xde, Rde_inv3, Aout)
    !---- args
    real(8), intent(in)  :: Xde(3)
    real(8), intent(in)  :: Rde_inv3
    real(8), intent(out) :: Aout(3)
    !---- vars
    !---- body
    Aout(:) = - Xde(:) * Rde_inv3

    return
  end subroutine Attraction_of_the_earth


  subroutine calc_RdeGEI_inverse(N, Rde, Rde_inv)
    !---- args
    integer, intent(in)  :: N
    real(8), intent(in)  :: Rde
    real(8), intent(out) :: Rde_inv(:)
    !---- vars
    integer :: l
    !---- body
    if(N >= 3)then
       Rde_inv(1) = 1.0d0 / Rde
       do l = 2, N
          Rde_inv(l) = Rde_inv(1) * Rde_inv(l-1)
       end do
    else if(0 <= N .and. N <= 2)then
       Rde_inv(2) = Rde_inv(1)*Rde_inv(1)
       Rde_inv(3) = Rde_inv(1)*Rde_inv(2)
    else
       print*, 'Error : Input order of geopotential is invalid'
       stop
    end if

    return
  end subroutine calc_RdeGEI_inverse

  
  subroutine Geopotential_Acceleration &
       (N_geo, XdeGEI, RdeGEI_inv, lat, lon, Ageo)
    !************************************************************************
    ! This subroutine is the core to calculate the geopotential perturbation.
    !************************************************************************
    !---- args
    integer, intent(in)  :: N_geo
    real(8), intent(in)  :: XdeGEI(3)
    real(8), intent(in)  :: RdeGEI_inv(:)
    real(8), intent(in)  :: lat
    real(8), intent(in)  :: lon
    real(8), intent(out) :: Ageo(3)
    !---- vars
    real(8) :: dUdR
    real(8) :: dUdPHI
    real(8) :: dUdLAMBDA
    !---- body
    !*******************************************************************
    ! First, expand associated Legendre functions up to specified order.
    !*******************************************************************
    call Associated_Legendre_functions(sin(lat),N_geo)
    !**************************************************************
    ! Next, calculate the acceleration in the spherical expression.
    !**************************************************************
    call perturbations_derived_from_geopotential &
         (N_geo, lat, lon, RdeGEI_inv, dUdR, dUdPHI, dUdLAMBDA)
    !*****************************************************************
    ! At last, Convert the acceleration into expression that refers to
    ! Geocentric Equatorial Inertial coordinate.
    !*****************************************************************
    call Sph2Car_acceleration &
         (XdeGEI, RdeGEI_inv, dUdR, dUdPHI, dUdLAMBDA, Ageo)

    return
  end subroutine Geopotential_Acceleration

  
  subroutine Associated_Legendre_functions(y,N)
    !**************************************************************
    ! This subroutine solve a recurrsion formula of associated
    ! Legendre functions numerically.
    ! Input   : Order of expansion terms(N)
    !           argument of associated Legendre functions(y)
    ! Output  : the values of associated Legendre functions
    ! Ref     : Arfken & Weber, Mathematical Methods for Physicists
    !                        -Sixth Edition- p775
    ! Author  : Keisuke Akari
    ! Date    : November 12, 2016
    !**************************************************************
    !---- args
    real(8), intent(in)  :: y
    integer, intent(in)  :: N
    !---- vars
    integer              :: m, l
    real(8)              :: constant
    !---- body
    do m = 0, N-2
       constant         = sqrt((1.0d0 - y**2)**(m))
       Legendre(m,m)    = permutation(2*m,m)*constant/dble(2**(m))
       Legendre(m+1,m)  = permutation(2*m+2,m+1)*constant*y/dble(2**(m+1))

       do l = m+2, N
          Legendre(l,m) = dble((2*l-1)*y*Legendre(l-1,m) - (l+m-1)*Legendre(l-2,m))/dble(l - m)
       end do
    end do
    
    constant            = sqrt((1.0d0 - y**2)**(N-1))
    Legendre(N-1,N-1)   = permutation(2*N-2,N-1)*constant/dble(2**(N-1))
    Legendre(N,N-1)     = permutation(2*N,N)*constant*y/dble(2**(N))
    constant            = sqrt((1.0d0 - y**2)**(N))
    Legendre(N,N)       = permutation(2*N,N)*constant/dble(2**(N))

    return
  end subroutine Associated_Legendre_functions


  subroutine perturbations_derived_from_geopotential &
       (N, lat, lon, RdeGEI_inv, dUdR, dUdPHI, dUdLAMBDA)
    !**************************************************************************
    ! This subroutine calculates perturbations derived from geopotential
    ! including high order expansion terms. Perturbations are calculated
    ! in spherical coordinates here.
    !
    ! Input   : position of a debri (X)
    !           (in geocentric spherical coordinate)
    !           Order of expansion term (N)
    ! Output  : Perturbations of geopotential
    !           (in geocentric spherical coordinate)
    ! Warning : Longitude should be measured positively east from the
    !           Greenwich meridian [rad]
    ! Ref     : David.A.Vallado, Fundamentals of Astrodynamics and Applications
    !                           - Third Edition- p548
    ! Author  : Keisuke Akari
    ! Date    : November 12, 2016
    !**************************************************************************
    !---- args
    integer, intent(in)  :: N
    real(8), intent(in)  :: lon
    real(8), intent(in)  :: lat
    real(8), intent(in)  :: RdeGEI_inv(:)
    real(8), intent(out) :: dUdR
    real(8), intent(out) :: dUdPHI
    real(8), intent(out) :: dUdLAMBDA
    !---- vars
    integer              :: l
    integer              :: m
    real(8)              :: constant(4)
    !---- body
    dUdR      = 0.0d0
    dUdPHI    = 0.0d0
    dUdLAMBDA = 0.0d0
    
    !****************************
    ! Sumation of Non-Zonal terms
    !****************************
    do l = N, 2, -1
       do m = l, 1, -1
          constant(1) = sin(m*lon)
          constant(2) = cos(m*lon)
          constant(3) = C_geo(m,l)*constant(2) + S_geo(m,l)*constant(1)

          dUdR      = dUdR      &
               + RdeGEI_inv(l)*(l+1)*Legendre(l,m)*constant(3)
          dUdPHI    = dUdPHI    &
               + RdeGEI_inv(l)*constant(3) &
               * (Legendre(l,m+1) - m*tan(lat)*Legendre(l,m))
          dUdLAMBDA = dUdLAMBDA &
               + RdeGEI_inv(l)*m*Legendre(l,m) &
               * (S_geo(m,l)*constant(2) - C_geo(m,l)*constant(1))
       end do
    end do
    
    !************************
    ! Sumation of Zonal terms
    !************************
    do l = N, 2, -1
       constant(4) = C_geo(0,l) * RdeGEI_inv(l)
       dUdR        = dUdR   + constant(4)*(l+1)*Legendre(l,0)
       dUdPHI      = dUdPHI + constant(4)*Legendre(l,1)
    end do
    
    !**************************
    ! Multiply the coefficients
    !**************************
    dUdR        = - RdeGEI_inv(2)*dUdR
    dUdPHI      =   RdeGEI_inv(1)*dUdPHI
    dUdLAMBDA   =   RdeGEI_inv(1)*dUdLAMBDA

    return
  end subroutine perturbations_derived_from_geopotential


  subroutine Sph2Car_acceleration &
       (XdeGEI, RdeGEI_inv, dUdR, dUdPHI, dUdLAMBDA, Aout)
    !---- args
    real(8), intent(in)  :: XdeGEI(3)
    real(8), intent(in)  :: RdeGEI_inv(:)
    real(8), intent(in)  :: dUdR
    real(8), intent(in)  :: dUdPHI
    real(8), intent(in)  :: dUdLAMBDA
    real(8), intent(out) :: Aout(3)
    !---- vars
    real(8)              :: constant(6)
    !---- body
    constant(1) = dUdR*RdeGEI_inv(1)
    constant(2) = 1.0d0/(XdeGEI(1)**2 + XdeGEI(2)**2)
    constant(3) = sqrt(constant(2))
    constant(4) = constant(3)*XdeGEI(3)*dUdPHI*RdeGEI_inv(2)
    constant(5) = constant(1) - constant(4)
    constant(6) = constant(2)*dUdLAMBDA
    
    Aout(1) = constant(5)*XdeGEI(1) - constant(6)*XdeGEI(2)
    Aout(2) = constant(5)*XdeGEI(2) + constant(6)*XdeGEI(1)
    Aout(3) = constant(1)*XdeGEI(3) + RdeGEI_inv(2)*dUdPHI/constant(3)

    return
  end subroutine Sph2Car_acceleration


  subroutine set_order_of_geopotential_polynomials(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    !---- body
    allocate(Legendre(0:N,0:N))
    allocate(C_normalized(0:N,0:N))
    allocate(S_normalized(0:N,0:N))        
    allocate(C_geo(0:N,0:N))
    allocate(S_geo(0:N,0:N))

    return
  end subroutine set_order_of_geopotential_polynomials


  subroutine input_geopotential_coefficients(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    integer             :: l, m, i, j
    !---- body
    open(fi1, file = 'input/EGM2008.dat',   status = 'old')
    
    do l = 2, N
       do m = 0, l
          read(fi1,*) i, j , C_normalized(m,l), S_normalized(m,l)
          if((i .ne. l) .or. (j .ne. m))then
             print*, 'Error : cannot read geopotential coefficients correctly'
             stop
          end if
       end do
    end do

    close(fi1)
    
    return
  end subroutine input_geopotential_coefficients


  subroutine set_geopotential_coefficients(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    integer             :: l, m
    real(8)             :: PI_Legendre
    !---- body
    do l = 2, N
       do m = 0, l
          if(m == 0)then
             PI_Legendre = sqrt(dble(permutation(l+m,l-m))/dble(2*l+1))
             C_geo(m,l) = C_normalized(m,l)/PI_Legendre
             S_geo(m,l) = S_normalized(m,l)/PI_Legendre
          else
             PI_Legendre = sqrt(dble(permutation(l+m,l-m))/dble(4*l+2))
             C_geo(m,l) = C_normalized(m,l)/PI_Legendre
             S_geo(m,l) = S_normalized(m,l)/PI_Legendre
          end if
       end do
    end do

    return
  end subroutine set_geopotential_coefficients
end module geopotential_mod
