module conversion_coordinates_mod
  use global_mod, only : &
       r0, pi2, pi, sec2rad
  use subprog_mod, only : &
       dot, cross, unit_vector, nolm, signal
  implicit none
contains
  subroutine precession(Vi,sze,cze,sz,cz,st,ct,Vo)
    !***********************************************************
    ! This subroutine can calculate direction cosine 
    ! considering the precession.
    !
    ! Input  : Position Vector (referred to J2000.0)
    !          Julian century (referred to J2000.0)
    ! Output : Mean Position Vector
    ! Ref    : 長沢 工, 天体の位置計算 - 増補版 - p49~56,228,229
    !          Jean Meeus, Astronomical Algorithms
    !              - Second Edition - p134,135
    ! Author : Keisuke Akari
    ! Date   : November 1, 2016
    !***********************************************************
    !---- args
    !---- vars
    real(8), intent(in)  :: Vi(3)
    real(8), intent(in)  :: sze,cze,sz,cz,st,ct
    real(8), intent(out) :: Vo(3)
    !---- body
    Vo(1) = (-sz*sze + cz*ct*cze)*Vi(1) &
         + (-sz*cze - cz*ct*sze)*Vi(2) &
         - cz*sze*Vi(3)
    Vo(2) = (-cz*sze + sz*ct*cze)*Vi(1) &
         + (cz*cze - sz*ct*sze)*Vi(2) &
         - sz*sze*Vi(3)
    Vo(3) = sz*cze*Vi(1) - st*sze*Vi(2) + ct*Vi(3)

    return
  end subroutine precession


  subroutine set_precession_matrix(T,sze,cze,sz,cz,st,ct)
    !***********************************************************
    ! This subroutine can calculate direction cosine 
    ! considering the precession.
    !
    ! Input  : Position Vector (referred to J2000.0)
    !          Julian century (referred to J2000.0)
    ! Output : Mean Position Vector
    ! Ref    : 長沢 工, 天体の位置計算 - 増補版 - p49~56,228,229
    !          Jean Meeus, Astronomical Algorithms
    !              - Second Edition - p134,135
    ! Author : Keisuke Akari
    ! Date   : November 1, 2016
    !***********************************************************
    !---- args
    real(8), intent(in)  :: T
    real(8), intent(out) :: sze,cze,sz,cz,st,ct
    !---- vars
    real(8) :: T2, T3
    real(8) :: zeta, z, theta
    !---- body
    T2 = T*T
    T3 = T*T2

    zeta  = (2306.2181d0*T + 3.0188d-1*T2 + 1.7998d-2*T3)*sec2rad
    z     = (2306.2181d0*T + 1.09468d0*T2 + 1.8203d-2*T3)*sec2rad
    theta = (2004.3109d0*T - 4.2665d-1*T2 - 4.1833d-2*T3)*sec2rad
    
    sze = sin(zeta)
    cze = cos(zeta)
    sz  = sin(z)
    cz  = cos(z)
    st  = sin(theta)
    ct  = cos(theta)

    return
  end subroutine set_precession_matrix


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
    real(8) :: r, dotXV
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


  subroutine Car2Sph(X, dis, dec, RA)
    !*********************************************************
    ! This subroutine cconverts cartesian coordinte into
    ! spherical coordinate.
    !
    ! Input   : Cartesian components (X,Y,Z)
    ! Output  : Spherical components
    !             dis = distance from the original point
    !             dec = geocentric latitude (-pi/2 to pi/2)
    !             RA  = Right Ascension, measured positively
    !                   east from the X-axis. (0 to 2pi)
    ! Warning : An unit of 'alt' depends on an unit of inputs
    ! Author  : Keisuke Akari
    ! Date    : November 3, 2016
    !*********************************************************
    !---- args
    real(8), intent(in)  :: X(3)
    real(8), intent(out) :: dis
    real(8), intent(out) :: dec
    real(8), intent(out) :: RA
    !---- vars
    real(8) :: R_XY
    real(8) :: aaa
    !---- body
    R_XY = X(1)**2 + X(2)**2
    aaa  = X(1)/sqrt(R_XY)
    if(aaa >= 1.0d0)then
       aaa = 0.999999999999999d0
    else if(aaa <= -1.0d0)then
       aaa = -0.999999999999999d0
    end if
       
    dis = sqrt(R_XY + X(3)**2)
    dec = atan(X(3)/sqrt(R_XY))
    RA = acos(aaa)
    if(X(2) < 0.0d0)then
       RA = pi2 - RA
    end if

    return
  end subroutine Car2Sph


  subroutine Car2NTW(XX,VV,Acar,Antw)
    !************************************************************
    ! Conversion from cartesian coordinate system into
    ! satellite coordinate system.
    ! Input  : XX(3)
    !          Position Vector in Cartesian coordinates
    !          VV(3)
    !          Velocity Vector in Cartesian coordinates
    !          Acar(3)
    !          Acceleration Vector in Cartesian coordinates
    ! Output : Antw(3)
    !          Acceleration Vector in Frenet coordinate system.
    !          The N axis is normal to the velocity vector.
    !          The T axis is tangential to the orbit and always
    !          points to the velocity vector.
    !          The W axis is normal to the orbital plane.
    ! Ref    : David. A. Vallado, Fundamentals of Astrodynamics
    !          and Applications, Third Edition, p163-165
    ! Date   : April 30, 2017
    ! Author : Keisuke Akari
    !************************************************************
    !---- args
    real(8), intent(in)  :: XX(3)
    real(8), intent(in)  :: VV(3)
    real(8), intent(in)  :: Acar(3)
    real(8), intent(out) :: Antw(3)
    !---- vars
    real(8) :: XXX(3)
    real(8) :: YYY(3)
    real(8) :: ZZZ(3)
    real(8) :: XcrossV(3)
    !---- body
    YYY(:)     = VV(:)/nolm(VV)
    XcrossV(:) = cross(XX,VV)
    ZZZ(:)     = XcrossV/nolm(XcrossV)
    XXX(:)     = cross(YYY,ZZZ)

    Antw(1) = dot(Acar,XXX) ! N : Radial component
    Antw(2) = dot(Acar,YYY) ! T : In-track component
    Antw(3) = dot(Acar,ZZZ) ! W : Cross-track component
    
    return
  end subroutine Car2NTW


  subroutine Car2RSW(XX,VV,Acar,Arsw)
    !***********************************************************
    ! Conversion from cartesian coordinate system into
    ! satellite coordinate system.
    ! Input  : XX(3)
    !          Position Vector in Cartesian coordinates
    !          VV(3)
    !          Velocity Vector in Cartesian coordinates
    !          Acar(3)
    !          Acceleration Vector in Cartesian coordinates
    ! Output : Arsw(3)
    !          Acceleration Vector in Gaussian coordinate system
    !          The R axis always points from the Earth's center
    !          along the radius vector toward the satellite.
    !          The S axis points in the direction (but not
    !          necessarily parallel to) the velocity vector and
    !          is perpendicular to the raius vector.
    !          The W axis is normal to the orbital plane. 
    ! Ref    : David. A. Vallado, Fundamentals of Astrodynamics
    !          and Applications, Third Edition, p163-165
    ! Date   : April 30, 2017 (modified on Nov. 24, 2017)
    ! Author : Keisuke Akari
    !***********************************************************
    !---- args
    real(8), intent(in)  :: XX(3)
    real(8), intent(in)  :: VV(3)
    real(8), intent(in)  :: Acar(3)
    real(8), intent(out) :: Arsw(3)
    !---- vars
    real(8) :: XXX(3) !dummy
    real(8) :: YYY(3)
    real(8) :: ZZZ(3)
    real(8) :: XcrossV(3)
    !---- body
    XXX(:)     = XX(:)/nolm(XX)
    XcrossV(:) = cross(XX,VV)
    ZZZ(:)     = XcrossV/nolm(XcrossV)
    YYY(:)     = cross(ZZZ,XXX)

    Arsw(1) = dot(Acar,XXX) ! R : Radial component
    Arsw(2) = dot(Acar,YYY) ! S : Along-track component
    Arsw(3) = dot(Acar,ZZZ) ! W : Cross-track component
    
    return
  end subroutine Car2RSW
  
    
  subroutine GEI2ECEF(XXX, YYY, theta)
    !---- args
    real(8), intent(in)  :: XXX(3), theta
    real(8), intent(out) :: YYY(3)
    !---- vars
    real(8) :: sth, cth
    !---- body
    sth = sin(theta)
    cth = cos(theta)

    YYY(1) = cth*XXX(1) - sth*XXX(2)
    YYY(2) = sth*XXX(1) + cth*XXX(2)
    YYY(3) = XXX(3)

    return
  end subroutine GEI2ECEF
  

  subroutine ECEF2LatLon(XX,phi,lambda,h)
    !**************************************************************
    ! This subroutine converts ECEF(Earth Centered, Earth Fixed)
    ! coordinate(X,Y,Z) into Geodetic coordinate(phi,lambda,h).
    !
    ! Input  : position vector in ECEF coordinates  (X,Y,Z [km])
    ! Output : position vector (Geodetic coordinate)
    !            h      : Height above the reference ellipsoid [km]
    !            phi    : Geodetic latitude, measured positively
    !                     north from the equator (-pi/2 to pi/2)
    !            lambda : Longitude, measured positively east
    !                     from the Greenwich meridian (0 to 2pi)
    ! Ref    : Fundamentals of Astrodynamics and Applications
    !                 - Third Edition - p180,181
    ! Author : Keisuke Akari
    ! Date   : Novenber 1, 2016
    !**************************************************************
    !---- args
    real(8), intent(in)  :: XX(3)
    real(8), intent(out) :: phi,lambda,h
    !---- vars
    real(8) :: r_sat, a, b, E, alpha
    real(8) :: F,     P, Q, D, nu
    real(8) :: G,     t
    !---- body
    r_sat = sqrt(XX(1)**2 + XX(2)**2)
    a     = r0 * 1.0d-3
    b     = 6356.75160056d0*signal(XX(3))
    E     = (b*XX(3) - (a*a - b*b)) / (a*r_sat)
    
    if(XX(2)/r_sat >= 0.0d0)then
       alpha = acos(XX(1)/r_sat)
    else
       alpha = pi2 - acos(XX(1)/r_sat)
    end if
    lambda = alpha

    F = (b*XX(3) + (a*a - b*b)) / (a*r_sat)
    P = 4.0d0*(E*F + 1.0d0)/3.0d0
    Q = 2.0d0*(E*E - F*F)
    D = P**3 + Q*Q

    if(D >= 0.0d0)then
       nu = (sqrt(D) - Q)**(1.0d0/3.0d0) - (sqrt(D) + Q)**(1.0d0/3.0d0)
    else if(D < 0.0d0)then
       nu = 2.0d0*sqrt(-P)*cos(acos(Q/(P*sqrt(-P)))/3.0d0)
    end if
    
    G = 0.5d0*(sqrt(E*E + nu) + E)
    t = sqrt(G*G + (F - nu*G)/(2.0d0*G - E)) - G
    
    phi = atan(a*(1.0d0 - t*t)/(2.0d0*b*t))

    h = (r_sat - a*t)*cos(phi) + (XX(3) - b)*sin(phi)

    return
  end subroutine ECEF2LatLon


  subroutine calc_longitude(XX, YY, lon)
    !---- args
    real(4), intent(in)  :: XX, YY
    real(4), intent(out) :: lon
    !---- vars
    real(8) :: aaa
    !---- body
    aaa = dble(XX/sqrt(XX*XX + YY*YY))
    if(aaa >= 1.0d0)then
       aaa = 0.999999999999999d0
    else if(aaa <= -1.0d0)then
       aaa = -0.999999999999999d0
    end if
    
    lon = real(acos(aaa))
    if(YY < 0.0d0)then
       lon = pi2 - lon
    end if

    return
  end subroutine calc_longitude
end module conversion_coordinates_mod
