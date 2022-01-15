module em_force_mod
  use global_mod, only  : &
       r0, r0_geopack, omega_e, pi, pi2, rad2deg, rad2lt
  use weimer_mod, only  : EpotVal
  use subprog_mod, only : cross_product, dot, nolm, unit_vector
  implicit none
contains
  subroutine IGRF_magnetic_field(YEFR_in, alt_gd_in, colat_gd_in, &
       lon_gd_in, colat_gc_in, lon_gc_in, Bout)
    !**************************************************
    ! This subroutine calculates magnetic field vector.
    ! Local spherical components of magneric field is
    ! calculated in 12th IGRF package, so conversion
    ! of coordinates is required.
    ! Subroutine, BSPCAR_08 is converter of coordinates
    ! from local spherical one into cartesian one.
    !
    ! Author : Keisuke Akari
    ! Date   : December 1, 2016
    !**************************************************
    !---- args
    real(8), intent(in)  :: YEFR_in
    real(8), intent(in)  :: alt_gd_in
    real(8), intent(in)  :: colat_gd_in
    real(8), intent(in)  :: lon_gd_in
    real(8), intent(in)  :: colat_gc_in
    real(8), intent(in)  :: lon_gc_in
    real(8), intent(out) :: Bout(3)
    !********************************************************
    ! YEFR_in     : a fraction of a year
    ! alt_gd_in   : height above the reference ellipsoid (km)
    ! colat_gd_in : geodetic colatitude (0 - pi [rad])
    ! lon_gd_in   : geodetic longitude  (0 - 2pi [rad])
    ! colat_gc_in : geocentric colatitude (0 - pi [rad])
    ! lon_gc_in   : geocentric longitude, measured positively
    !               east from the direction of vernal equinox
    ! Bout(3)     : Magnetic Flux Density Vector (T)
    !********************************************************
    !---- vars
    real(8) :: Br
    real(8) :: Btheta
    real(8) :: Bphi
    real(4) :: BnGEI(3)
    real(8) :: Bn_nolm
    !---- body
    call igrf12syn(0, YEFR_in, 1, alt_gd_in, colat_gd_in*rad2deg, &
         lon_gd_in*rad2deg, Btheta, Bphi, Br, Bn_nolm)
    !****************************************************
    ! Warning : Downward direction is positive for the
    !           vertical component of geomagnetic fields.
    !****************************************************
    call BSPCAR_08(real(colat_gc_in), real(lon_gc_in), real(-Br), &
         real(-Btheta), real(Bphi), BnGEI(1), BnGEI(2), BnGEI(3))

    Bout(:) = dble(BnGEI(:))*1.0d-9

    return
  end subroutine IGRF_magnetic_field


  subroutine dipole_magnetic_field(XGEO, BGEO)
    !****************************************************
    ! This subroutine calculate the magnetic flux density
    ! vector at the debris' position. To find the value,
    ! dipole magnetic field model is employed here.
    !
    ! Input  : Position Vector of Debris [Normalized]
    ! Output : The Magnetic Flux Density Vector [T]
    ! Ref    : 人工衛星の力学と制御ハンドブック, p288-293
    ! Author : Keisuke Akari
    ! Date   : 6 July, 2017
    !****************************************************
    !---- args
    real(8), intent(in)  :: XGEO(3) !地球固定座標系
    real(8), intent(out) :: BGEO(3) !地球固定座標系
    !---- vars
    ! d(3) = unit vector of magnetic polar
    real(8), parameter :: d(3) = (/0.0d0,0.0d0,-1.0d0/)
    real(8) :: vec(3)
    real(8) :: vec_unit(3)
    real(8) :: R
    !---- body
    R = nolm(XGEO)
    vec_unit(:) = XGEO(:)/R
    vec(:)      = 3.0d0*dot(d(:),vec_unit(:))*vec_unit(:) - d(:)
    BGEO(:)     = 3.0257d-5*(6.3712d6/(r0*R))**3*vec(:)

    return
  end subroutine dipole_magnetic_field


  subroutine tilted_dipole_magnetic_field(XGEO, MJDN, DAFR, BGEO)
    !****************************************************
    ! This subroutine calculate the magnetic flux density
    ! vector at the debris' position. To find the value,
    ! dipole magnetic field model is employed here.
    !
    ! Warning : 一時的なプログラムであるから, 2005年から
    !           2020年までしか適応していない
    ! Input   : Position Vector of Debris [Normalized]
    ! Output  : The Magnetic Flux Density Vector [T]
    ! Ref     : 人工衛星の力学と制御ハンドブック, p288-293
    ! Author  : Keisuke Akari
    ! Date    : 6 July, 2017
    !****************************************************
    !---- args
    real(8), intent(in)  :: XGEO(3) !地球固定座標系
    integer, intent(in)  :: MJDN
    real(8), intent(in)  :: DAFR
    real(8), intent(out) :: BGEO(3) !地球固定座標系
    !---- vars
    real(8) :: g10_2005 = -29554.63d0 !IGRF coefficients
    real(8) :: g11_2005 = -1669.05d0
    real(8) :: h11_2005 = 5077.99d0
    real(8) :: g10_dot_2005 = 8.8d0   !Secular variation per year
    real(8) :: g11_dot_2005 = 10.8d0
    real(8) :: h11_dot_2005 = -21.3d0

    real(8) :: g10_2010 = -29496.57d0
    real(8) :: g11_2010 = -1586.42d0
    real(8) :: h11_2010 = 4944.26d0
    real(8) :: g10_dot_2010 = 11.4d0
    real(8) :: g11_dot_2010 = 16.7d0
    real(8) :: h11_dot_2010 = -28.8d0

    real(8) :: g10_2015 = -29442.0d0
    real(8) :: g11_2015 = -1501.0d0
    real(8) :: h11_2015 = 4797.1d0
    real(8) :: g10_dot_2015 = 10.3d0
    real(8) :: g11_dot_2015 = 18.1d0
    real(8) :: h11_dot_2015 =-26.6d0

    real(8) :: YY
    real(8) :: g11
    real(8) :: g10
    real(8) :: h11
    real(8) :: d(3)
    real(8) :: d_unit(3)
    real(8) :: vec(3)
    real(8) :: vec_unit(3)
    real(8) :: R
    !---- body
    if(53371 <= MJDN .and. MJDN < 55197)then !2005年~2010年
       YY = (MJDN - 53371.0d0 + DAFR/86400.0d0)/365.25 !2005年1月1日0時0分からの経過年数
    
       g11 = g11_2005 + g11_dot_2005*YY
       h11 = h11_2005 + h11_dot_2005*YY
       g10 = g10_2005 + g10_dot_2005*YY
    else if(55197 <= MJDN .and. MJDN < 57023)then !2010年~2015年
       YY = (MJDN - 55197.0d0 + DAFR/86400.0d0)/365.25 !2010年1月1日0時0分からの経過年数
    
       g11 = g11_2010 + g11_dot_2010*YY
       h11 = h11_2010 + h11_dot_2010*YY
       g10 = g10_2010 + g10_dot_2010*YY
    else if(57023 <= MJDN .and. MJDN < 58849)then !2015年~2020年
       YY = (MJDN - 57023.0d0 + DAFR/86400.0d0)/365.25 !2015年1月1日0時0分からの経過年数
    
       g11 = g11_2015 + g11_dot_2015*YY
       h11 = h11_2015 + h11_dot_2015*YY
       g10 = g10_2015 + g10_dot_2015*YY
    else
       print*, "Error : Input time is invalid &
            to calculate the tilted dipole magnetic field."
       stop
    end if
    
    d(1) = g11
    d(2) = h11
    d(3) = g10

    d_unit(:) = unit_vector(d(:))
    
    R = nolm(XGEO)
    vec_unit(:) = XGEO(:)/R
    vec(:)      = 3.0d0*dot(d_unit(:),vec_unit(:))*vec_unit(:) - d_unit(:)
    BGEO(:)     = 3.0257d-5*(6.3712d6/(r0*R))**3*vec(:)

    return
  end subroutine tilted_dipole_magnetic_field
  

  subroutine corotation_electric_field(Xin,Bin,Eout)
    !***********************************************
    ! This subroutine sets corotation electric field
    !***********************************************
    !---- args
    real(8), intent(in)  :: Xin(3)
    real(8), intent(in)  :: Bin(3)
    real(8), intent(out) :: Eout(3)
    !*****************************************
    ! Xin  : Position Vector of the Debris (m)
    ! Bin  : Magnetic Flux Density Vector (T)
    ! Eout : Corotation Electric Field (V/m)
    !*****************************************
    !---- vars
    real(8) :: omegar(3)
    !---- body
    call cross_product(Xin, omega_e, omegar)
    call cross_product(omegar, Bin, Eout)

    return
  end subroutine corotation_electric_field


  subroutine set_convection_electric_field_coefficients(Din,XGRID,YGRID,ZGRID)
    !---- args
    real(8), intent(in)  :: Din
    real(8), intent(out) :: XGRID(3)
    real(8), intent(out) :: YGRID(3)
    real(8), intent(out) :: ZGRID(3)
    !*******************************
    ! Din : Grid length (Normalized)
    !*******************************
    !---- vars
    !---- body
    XGRID(1) = Din
    XGRID(2) = 0.0d0
    XGRID(3) = 0.0d0
    
    YGRID(1) = 0.0d0
    YGRID(2) = Din
    YGRID(3) = 0.0d0
    
    ZGRID(1) = 0.0d0
    ZGRID(2) = 0.0d0
    ZGRID(3) = Din

    return
  end subroutine set_convection_electric_field_coefficients
  
  
  subroutine convection_electric_field(Xin,DXYZ,XGRID,YGRID,ZGRID,Eout)
    !********************************************************
    ! This subroutine sets convection electric field without
    ! tracing a footprint along geomagnetic line.
    ! Therefore, no dependency of an altitude on Earth'
    ! surface electric potential is assumed.
    ! In order to convert coordinates, Geopack2008 is used in
    ! this program, and Weimer model is employed for the
    ! calculation of Earth' surface electric potential.
    ! Be careful about variable type.
    ! Both package cannot treat double precision real number.
    !
    ! Author : Keisuke Akari
    ! Date   : December 3, 2016
    !********************************************************
    !---- args
    real(8), intent(in)  :: Xin(3)
    real(8), intent(in)  :: DXYZ
    real(8), intent(in)  :: XGRID(3)
    real(8), intent(in)  :: YGRID(3)
    real(8), intent(in)  :: ZGRID(3)
    real(4), intent(out) :: Eout(3)
    !---- vars
    integer :: i
    real(4) :: XGEO(3)
    real(4) :: XGSW(3)
    real(4) :: XSM(3)
    real(4) :: MLAT
    real(4) :: MLT
    real(4) :: Xplus(3)
    real(4) :: Xminus(3)
    real(4) :: Yplus(3)
    real(4) :: Yminus(3)
    real(4) :: Zplus(3)
    real(4) :: Zminus(3)    
    real(4) :: Eplus(3)
    real(4) :: Eminus(3)
    !---- body
    Xplus(:)  = (Xin(:) + XGRID(:))*r0_geopack
    Xminus(:) = (Xin(:) - XGRID(:))*r0_geopack
    Yplus(:)  = (Xin(:) + YGRID(:))*r0_geopack
    Yminus(:) = (Xin(:) - YGRID(:))*r0_geopack
    Zplus(:)  = (Xin(:) + ZGRID(:))*r0_geopack
    Zminus(:) = (Xin(:) - ZGRID(:))*r0_geopack

    call GEIGEO_08(Xplus(1),Xplus(2),Xplus(3),XGEO(1),XGEO(2),XGEO(3),1)
    call GEOGSW_08(XGEO(1),XGEO(2),XGEO(3),XGSW(1),XGSW(2),XGSW(3),1)
    call SMGSW_08(XSM(1),XSM(2),XSM(3),XGSW(1),XGSW(2),XGSW(3),-1)
    call SM2MLATMLT(XSM, MLAT, MLT)
    Eplus(1) = EpotVal(abs(MLAT), MLT)

    call GEIGEO_08(Xminus(1),Xminus(2),Xminus(3),XGEO(1),XGEO(2),XGEO(3),1)
    call GEOGSW_08(XGEO(1),XGEO(2),XGEO(3),XGSW(1),XGSW(2),XGSW(3),1)
    call SMGSW_08(XSM(1),XSM(2),XSM(3),XGSW(1),XGSW(2),XGSW(3),-1)
    call SM2MLATMLT(XSM, MLAT, MLT)
    Eminus(1) = EpotVal(abs(MLAT), MLT)
    
    call GEIGEO_08(Yplus(1),Yplus(2),Yplus(3),XGEO(1),XGEO(2),XGEO(3),1)
    call GEOGSW_08(XGEO(1),XGEO(2),XGEO(3),XGSW(1),XGSW(2),XGSW(3),1)
    call SMGSW_08(XSM(1),XSM(2),XSM(3),XGSW(1),XGSW(2),XGSW(3),-1)
    call SM2MLATMLT(XSM, MLAT, MLT)
    Eplus(2) = EpotVal(abs(MLAT), MLT)
    
    call GEIGEO_08(Yminus(1),Yminus(2),Yminus(3),XGEO(1),XGEO(2),XGEO(3),1)
    call GEOGSW_08(XGEO(1),XGEO(2),XGEO(3),XGSW(1),XGSW(2),XGSW(3),1)
    call SMGSW_08(XSM(1),XSM(2),XSM(3),XGSW(1),XGSW(2),XGSW(3),-1)
    call SM2MLATMLT(XSM, MLAT, MLT)
    Eminus(2) = EpotVal(abs(MLAT), MLT)
    
    call GEIGEO_08(Zplus(1),Zplus(2),Zplus(3),XGEO(1),XGEO(2),XGEO(3),1)
    call GEOGSW_08(XGEO(1),XGEO(2),XGEO(3),XGSW(1),XGSW(2),XGSW(3),1)
    call SMGSW_08(XSM(1),XSM(2),XSM(3),XGSW(1),XGSW(2),XGSW(3),-1)
    call SM2MLATMLT(XSM, MLAT, MLT)
    Eplus(3) = EpotVal(abs(MLAT), MLT)
    
    call GEIGEO_08(Zminus(1),Zminus(2),Zminus(3),XGEO(1),XGEO(2),XGEO(3),1)
    call GEOGSW_08(XGEO(1),XGEO(2),XGEO(3),XGSW(1),XGSW(2),XGSW(3),1)
    call SMGSW_08(XSM(1),XSM(2),XSM(3),XGSW(1),XGSW(2),XGSW(3),-1)
    call SM2MLATMLT(XSM, MLAT, MLT)
    Eminus(3) = EpotVal(abs(MLAT), MLT)

    do i = 1, 3
       if(Eplus(i) == 0.0d0 .or. Eminus(i) == 0.0d0)then
          Eout(:) = 0.0d0
          return
       end if
    end do
    
    Eout(:) = -1.0d3*((Eplus(:) - Eminus(:)) / (2.0d0 * DXYZ))

    return
  end subroutine convection_electric_field


  subroutine SM2MLATMLT(XX, mlat, mlt)
    !********************************************************
    ! This subroutine calculates geomagnetic latitude and
    ! geomagnetic local time from Solar Magnetic coordinates
    ! (SM). In SM coordinates, the Z-axis is chosen parallel
    ! to the north magnetic pole and the Y-axis perpendicular
    ! to the Earth-Sun line towards dusk. 
    ! Input  : Geomagnetic coordinate
    ! Output : Geomagnetic latitude   [degree]
    !          Geomagnetic local time [hour]
    !          (MLT = 12 at noon, MLT = 0,24 at midnight)
    ! Author : Keisuke Akari
    ! Date   : December 1, 2016
    !******************************************************
    !---- args
    real(4), intent(in)  :: XX(3)
    real(4), intent(out) :: mlat, mlt
    !---- vars
    real(4) :: Rxy
    real(4) :: theta
    !---- body
    Rxy = sqrt(XX(1)*XX(1) + XX(2)*XX(2))
    
    if(Rxy == 0.0d0)then
       mlat = 0.0d0
       mlt  = 0.0d0
       return
    end if
    
    mlat = atan2(XX(3),Rxy)*rad2deg

    theta = acos(XX(1)/Rxy)
    
    if(XX(2) < 0.0d0)then
       theta = pi2 - theta
    end if

    if(0.0d0 <= theta .and. theta < pi)then
       mlt = (theta + pi)*rad2lt
    else if(pi <= theta .and. theta < pi2)then
       mlt = (theta - pi)*rad2lt
    end if
    
    return
  end subroutine SM2MLATMLT
end module em_force_mod
