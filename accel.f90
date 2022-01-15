module accel_mod
  use global_mod, only : &
       r0, v0, g0, t0, eps0, pi, pi2, pi_half, &
       charge_d, Rd, Md, Ad, Ud, rad2deg,&
       UT_year, UT_DoY, UT_hour, UT_min, UT_sec, time, &
       UT_DAFR, UT_MJDN, UT_YEFR, TD_JD, D_aphelion, &
       VGSEX, VGSEY, VGSEZ, &
       alt_gd, lat_gc, lat_gd, lon_gc, lon_gd, colat_gc, colat_gd, &
       pi_half, XdeGEO, N_geo, &
       BGEI, BGEO, VcrossB, Ecn, Ecr, &
       Aair, Aph, Amoon, Asun, Ageo, Aearth, Atotal, &
       N_air, &
       RA_sun, dec_sun, &
       XdeltaGEI, YdeltaGEI, ZdeltaGEI, deltaXYZ, &
       moon_interval, index_moon, sun_interval, index_sun, &
       XdeGEI, RdeGEI, RdeGEI_inv, VdeGEI, &
       XmeGEI, RmeGEI, RmeGEI_inv, &
       XdmGEI, RdmGEI, RdmGEI_inv, &
       XseGEI, RseGEI, RseGEI_inv, &
       XdsGEI, RdsGEI, RdsGEI_inv, &
       keisuu_airdrag, keisuu_sun, keisuu_moon, keisuu_ph, &
       Rsun, air_frag
  use time_mod, only : JD2CalenderDate
  use geopotential_mod, only : &
       calc_RdeGEI_inverse, &
       Attraction_of_the_earth, &
       Geopotential_Acceleration
  use air_drag_mod, only : air_drag
  use moon_mod, only : &
       calc_moon_position, &
       distance_from_the_moon, &
       attraction_of_the_moon
  use sun_mod, only : &
       calc_sun_position, &
       distance_from_the_sun1, &
       distance_from_the_sun2, &
       sun_attraction, &
       solar_radiation_pressure
  use em_force_mod, only : &
       IGRF_magnetic_field, &
       dipole_magnetic_field, &
       tilted_dipole_magnetic_field, &
       corotation_electric_field, &
       convection_electric_field
  use subprog_mod, only : cross_product, nolm
  use conversion_coordinates_mod, only : &
       calc_longitude, &
       Car2Sph
  use Charge, only : Renew_Charge
  use charge_solve_mod, only : solve_current_balance
  implicit none
contains
  !*******************************************************
  ! When time is renewed, call renew_acceleration1.
  ! Otherwise, call renew_acceleration2.
  ! The difference of these two subroutines is whether
  ! the positions of celestial bodies are renewed, or not.
  !*******************************************************
  subroutine renew_acceleration1
    !---- args
    !---- vars
    real(4) :: BB(3)
    real(8) :: dhour
    integer :: year, month, day, mmdd, shadow_function    
    !---- body
    !********************************************************
    ! Prepare debri's position expressed in some coordinates.
    !********************************************************
    ! When time is renewd, RECALC_08 should be called before
    ! conversion of coordinates.  
    !********************************************************
    call RECALC_08(UT_year, UT_DoY, UT_hour, UT_min, &
         int(UT_sec + time*t0), VGSEX, VGSEY, VGSEZ)
    
    !***********************************************************
    !'Car2Sph' converts cartesian components into sphrical ones.
    !***********************************************************
    call Car2Sph(XdeGEI, RdeGEI, lat_gc, lon_gc)
    colat_gc = pi_half - lat_gc

    !*************************************************************
    ! 'GEODGEO_08' converts geocentric spherical components
    ! into geodetic spherical ones. This part calculates
    ! geodetic altitude and latitude from geocentric counterparts.
    !*************************************************************
    call GEODGEO_08 &
         (alt_gd, lat_gd, real(RdeGEI*r0*1.0d-3),real(colat_gc), -1)
    if(pi_half <= lat_gd)then
       lat_gd = real(pi_half - 1.0d-6)
    else if(lat_gd <= -pi_half)then
       lat_gd = real(-pi_half + 1.0d-6)
    end if
    
    colat_gd = pi_half - dble(lat_gd)

    !***************************************************************
    ! 'GEIGEO_08' converts geocentric eqatorial cartesian components
    ! into geodetic cartesian ones. This part calculates the geodetic
    ! components at the position in question.
    !***************************************************************
    call GEIGEO_08(real(XdeGEI(1)), real(XdeGEI(2)), real(XdeGEI(3)), &
         XdeGEO(1), XdeGEO(2), XdeGEO(3), 1)
    call calc_longitude(XdeGEO(1), XdeGEO(2), lon_gd)
    if(lon_gd >= pi2)then
       lon_gd = real(lon_gd - pi2)
    end if
    !*******************************************************
    call calc_RdeGEI_inverse(N_geo, RdeGEI, RdeGEI_inv)

    call calc_sun_position(TD_JD, sun_interval, index_sun, XseGEI)
    call distance_from_the_sun1(XseGEI, RseGEI, RseGEI_inv)
    call distance_from_the_sun2 &
         (XdeGEI, XseGEI, XdsGEI, RdsGEI, RdsGEI_inv)
    call Car2Sph(XseGEI, RseGEI, dec_sun, RA_sun)   

    call calc_moon_position(TD_JD, moon_interval, index_moon, XmeGEI)
    call distance_from_the_moon &
         (XdeGEI,XmeGEI,XdmGEI,RdmGEI,RdmGEI_inv,RmeGEI,RmeGEI_inv)
    
    !***************************************************************
    ! Perturbations derived from the Earth's electromagnetic fields.
    !***************************************************************
    call sun_attraction &
         (keisuu_sun, XdsGEI, RdsGEI_inv(3), &
         XseGEI, RseGEI_inv(3), Asun)
    call solar_radiation_pressure(keisuu_ph, D_aphelion, Rsun, &
         XdeGEI, XdsGEI, XseGEI, RdeGEI, RdsGEI, &
         RdeGEI_inv(1), RdsGEI_inv(1), Aph, shadow_function)

    !****************************************************
    ! See 'Makefile' in order to understand how to switch
    ! perturbing force considered in this model.
    !****************************************************
#ifdef change_in_charge    
    call Renew_Charge(dble(alt_gd), lat_gc*rad2deg, lon_gc*rad2deg, time*t0, charge_d)
#endif

#ifdef charge_solve
    call JD2CalenderDate(UT_MJDN+UT_DAFR/86400.0d0+2400000.5d0, &
         year,month,day)
    mmdd  = month*100+day
    dhour = UT_DAFR/3600.0d0
    call solve_current_balance(4.0d0*Ad,v0*nolm(VdeGEI(:)),dble(alt_gd), &
         dble(lat_gd*rad2deg),dble(lon_gd*rad2deg), &
         year,mmdd,dhour,shadow_function,Ud)!input shadow_function; 0=sunlit,2=ecllipse
    charge_d = 4.0d0 * pi * eps0 * Rd * Ud
#endif
    
#ifdef IGRF
       call IGRF_magnetic_field(UT_YEFR, dble(alt_gd), dble(colat_gd), &
            dble(lon_gd), colat_gc, lon_gc, BGEI)
#endif

#ifdef ideal_dipole
       call dipole_magnetic_field(dble(XdeGEO), BGEO)
       call GEIGEO_08(BB(1), BB(2), BB(3), &
            real(BGEO(1)), real(BGEO(2)), real(BGEO(3)), -1)
       BGEI(:) = BB(:)
#endif

#ifdef tilted_dipole
       call tilted_dipole_magnetic_field(dble(XdeGEO), UT_MJDN, UT_DAFR, BGEO)
       call GEIGEO_08(BB(1), BB(2), BB(3), &
            real(BGEO(1)), real(BGEO(2)), real(BGEO(3)), -1)
       BGEI(:) = BB(:)
#endif
       
#ifdef mag
       call cross_product((v0*VdeGEI),BGEI,VcrossB)
#endif
       
#ifdef ecor
       call corotation_electric_field((r0*XdeGEI),BGEI,Ecr)
#endif
       
#ifdef econv
       call convection_electric_field &
            (XdeGEI,deltaXYZ,XdeltaGEI,YdeltaGEI,ZdeltaGEI,Ecn)
#endif
    
    !************************************************
    ! Perturbations derived from the dynamical force.
    !************************************************
    call Attraction_of_the_earth(XdeGEI, RdeGEI_inv(3), Aearth)
    call Geopotential_Acceleration &
         (N_geo, XdeGEI, RdeGEI_inv, lat_gc, dble(lon_gd), Ageo)
#ifdef drag    
    call Air_Drag(N_air, UT_MJDN, UT_DAFR, &
         alt_gd*1.0d3, dble(lat_gd), lat_gc, lon_gc, RA_sun, dec_sun, &
         keisuu_airdrag, XdeGEI, VdeGEI, air_frag, Aair)
#endif
    
    call Attraction_of_the_moon &
         (keisuu_moon,XdmGEI,RdmGEI_inv(3),XmeGEI,RmeGEI_inv(3),Amoon)

    !***********************************
    ! Sum of all kinds of perturbations.
    !***********************************
    call add_up_to_acceleration

    return
  end subroutine renew_acceleration1


  subroutine renew_acceleration2
    !---- args
    !---- vars
    real(4) :: BB(3)
    real(8) :: dhour
    integer :: year, month, day, mmdd, shadow_function
    !---- body    
    !*******************************************************
    ! Prepare debri's position expressed in some coordinates
    !*******************************************************    
    call Car2Sph(XdeGEI, RdeGEI, lat_gc, lon_gc)
    colat_gc = pi_half - lat_gc
    call GEODGEO_08 &
         (alt_gd, lat_gd, real(RdeGEI*r0*1.0d-3), real(colat_gc), -1)
    if(pi_half <= lat_gd)then
       lat_gd = real(pi_half - 1.0d-6)
    else if(lat_gd <= -pi_half)then
       lat_gd = real(-pi_half + 1.0d-6)
    end if
    
    call GEIGEO_08(real(XdeGEI(1)), real(XdeGEI(2)), real(XdeGEI(3)), &
         XdeGEO(1), XdeGEO(2), XdeGEO(3), 1)
    call calc_longitude(XdeGEO(1), XdeGEO(2), lon_gd)
    if(lon_gd >= pi2)then
       lon_gd = real(lon_gd - pi2)
    end if    
    !*******************************************************
    call calc_RdeGEI_inverse(N_geo, RdeGEI, RdeGEI_inv)

    call distance_from_the_sun2 &
         (XdeGEI, XseGEI, XdsGEI, RdsGEI, RdsGEI_inv)
    call distance_from_the_moon &
         (XdeGEI,XmeGEI,XdmGEI,RdmGEI,RdmGEI_inv,RmeGEI,RmeGEI_inv)
    
    !***************************************************************
    ! Perturbations derived from the Earth's electromagnetic fields.
    !***************************************************************
    call sun_attraction &
         (keisuu_sun, XdsGEI, RdsGEI_inv(3), &
         XseGEI, RseGEI_inv(3), Asun)
    call solar_radiation_pressure(keisuu_ph, D_aphelion, Rsun, &
         XdeGEI, XdsGEI, XseGEI, RdeGEI, RdsGEI, &
         RdeGEI_inv(1), RdsGEI_inv(1), Aph, shadow_function)

    !****************************************************
    ! See 'Makefile' in order to understand how to switch
    ! perturbing force considered in this model.
    !****************************************************
#ifdef change_in_charge    
    call Renew_Charge(dble(alt_gd), lat_gc*rad2deg, lon_gc*rad2deg, time*t0, charge_d)
#endif

#ifdef charge_solve
    call JD2CalenderDate(UT_MJDN+UT_DAFR/86400.0d0+2400000.5d0, &
         year,month,day)
    mmdd  = month*100+day
    dhour = UT_DAFR/3600.0d0
    call solve_current_balance(4.0d0*Ad,v0*nolm(VdeGEI(:)),dble(alt_gd), &
         dble(lat_gd*rad2deg),dble(lon_gd*rad2deg), &
         year,mmdd,dhour,shadow_function,Ud)!input shadow_function; 0=sunlit,2=ecllipse
    charge_d = 4.0d0 * pi * eps0 * Rd * Ud
#endif    

#ifdef IGRF
    call IGRF_magnetic_field(UT_YEFR, dble(alt_gd), dble(colat_gd), &
         dble(lon_gd), colat_gc, lon_gc, BGEI)
#endif

#ifdef ideal_dipole
       call dipole_magnetic_field(dble(XdeGEO), BGEO)
       call GEIGEO_08(BB(1), BB(2), BB(3), &
            real(BGEO(1)), real(BGEO(2)), real(BGEO(3)), -1)
       BGEI(:) = BB(:)
#endif

#ifdef tilted_dipole
       call tilted_dipole_magnetic_field(dble(XdeGEO), UT_MJDN, UT_DAFR, BGEO)
       call GEIGEO_08(BB(1), BB(2), BB(3), &
            real(BGEO(1)), real(BGEO(2)), real(BGEO(3)), -1)
       BGEI(:) = BB(:)
#endif
    
#ifdef mag
    call cross_product((v0*VdeGEI), BGEI, VcrossB)
#endif
    
#ifdef ecor
    call corotation_electric_field((r0*XdeGEI),BGEI,Ecr)
#endif
    
#ifdef econv
    call convection_electric_field &
         (XdeGEI,deltaXYZ,XdeltaGEI,YdeltaGEI,ZdeltaGEI,Ecn)
#endif
    
    !*************************************************
    ! Perturbations derived from the mechanical force.
    !*************************************************
    call Attraction_of_the_earth(XdeGEI, RdeGEI_inv(3), Aearth)
    call Geopotential_Acceleration &
         (N_geo, XdeGEI, RdeGEI_inv, lat_gc, dble(lon_gd), Ageo)
#ifdef drag    
    call Air_Drag(N_air, UT_MJDN, UT_DAFR, &
         alt_gd*1.0d3, dble(lat_gd), lat_gc, lon_gc, RA_sun, dec_sun, &
         keisuu_airdrag, XdeGEI, VdeGEI, air_frag, Aair)
#endif    
    call Attraction_of_the_moon &
         (keisuu_moon,XdmGEI,RdmGEI_inv(3),XmeGEI,RmeGEI_inv(3),Amoon)
    
    !***********************************
    ! Sum of all kinds of perturbations.
    !***********************************
    call add_up_to_acceleration

    return
  end subroutine renew_acceleration2

  
  subroutine add_up_to_acceleration
    Atotal(:) = charge_d*(Ecn(:) + Ecr(:) + VcrossB(:))/(Md*g0) &
         + Aair(:) + Aph(:) + Amoon(:) + Asun(:) + Ageo(:)  &
         + Aearth

    return
  end subroutine add_up_to_acceleration
end module accel_mod
