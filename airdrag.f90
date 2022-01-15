module air_drag_mod
  !**************************************************
  ! This is the module for calculation of atmospheric
  ! drag. In this library, Jacchia-Roberts model is
  ! employed to know the atmospheric density at the
  ! point in question.
  !**************************************************
  use global_mod,          only : r0, t0, omega_e
  use jacchia_roberts_mod, only : RSDAMO, SFDJ70
  use Exponential_Model,   only : Exponential_Density
  implicit none
contains
  subroutine set_airdrag_constants(Cin, Ain, Min, Cout)
    !---- args
    real(8), intent(in)  :: Cin
    real(8), intent(in)  :: Ain
    real(8), intent(in)  :: Min
    real(8), intent(out) :: Cout
    !***********************************************
    ! Cin  : Drag Coefficient [Dimensionless]       
    ! Ain  : Reference Area [m^2]                   
    ! Min  : Mass [kg]                              
    ! Cout : Coefficient for Calculation of Air Drag
    !***********************************************
    !---- vars
    !---- body
    Cout = -0.5d0*Cin*Ain*r0/Min

    return
  end subroutine set_airdrag_constants
  
  
  subroutine Air_Drag(Nin, MJDN_in, DAFR_in, alt_gd, lat_gd, &
       lat_gc, lon_gc, RA_sun_in, dec_sun_in, Cin, &
       Xin, Vin, flag, Aout)
    !*******************************************************
    ! This subroutine is the main program to calculate an
    ! atmospheric drag.
    ! First, call 'RDYMOS' to find the atmospheric density
    ! of the point in question.
    ! Next, calculate the relative velocity of a debri to
    ! the atmosphere.
    ! Finally, call 'calc_air_drag' so that calculate the
    ! perturbation.
    !
    ! Warning :
    !  In this program, data published by INPE is employed
    !  for an averaged daily observed solar flux defined
    !  by Jacchia.
    !  The data held in the file have a limit of time span :
    !  Janualy 1, 1958 to December 31, 2008.
    !  Therefore, in this simulation, time should be gone
    !  back when referencing the date which is not described
    !  in the file.
    !
    !  Author : Keisuke Akari
    !  Date   : November 30, 2016
    !*******************************************************
    !---- args
    integer, intent(in)    :: Nin
    integer, intent(in)    :: MJDN_in
    real(8), intent(in)    :: DAFR_in
    real(8), intent(in)    :: alt_gd
    real(8), intent(in)    :: lat_gd
    real(8), intent(in)    :: lat_gc
    real(8), intent(in)    :: lon_gc    
    real(8), intent(in)    :: RA_sun_in
    real(8), intent(in)    :: dec_sun_in
    real(8), intent(in)    :: Cin
    real(8), intent(in)    :: Xin(3)
    real(8), intent(in)    :: Vin(3)
    integer, intent(inout) :: flag
    real(8), intent(out)   :: Aout(3)
    !***************************************************************
    ! Nin        : Latest Julian Day Number recorded in the input
    !              file which contains data of the Sun's activities.
    ! MJDN_in    : Modified Julian Day Number (UT)
    ! DAFR_in    : a fraction of a day in sec (UT) [sec]
    ! alt_gd     : Geodetic Altitude of the Debris [m]
    ! lat_gd     : Geodetic Latitude [rad]
    ! alt_gc     : Geocentric Altitude [m]
    ! lon_gc     : Right Ascension of the Debris [rad]
    ! lat_gc     : Geocentric Latitude of the Debris [rad]
    ! RA_sun_in  : Right Ascension of the Sun [rad]
    ! dec_sun_in : Declination of the Sun [rad]
    ! Cin        : Coefficient for calculation of Atmospheric Drag
    ! Xin(3)     : Position Vector of the Debris (Normalized)
    ! Vin(3)     : Velocity Vector of the Debris (Normalized)
    ! Aout(3)    : Acceleration Vector of the Atmosphric Drag
    !              (Normalized)
    !***************************************************************
    !---- vars
    real(8) :: alt_gc
    real(8) :: VdaGEI(3)
    real(8) :: rho_air
    !---- body
#ifdef JR_model
    alt_gc = alt_gd/cos(lat_gd - lat_gc)
    call Jacchia_Roberts_Model(Nin, MJDN_in, DAFR_in, &
       alt_gc, lon_gc, lat_gc, RA_sun_in, dec_sun_in, flag, rho_air)
#endif
#ifdef Expo_model
    call Exponential_Density(alt_gd, rho_air)
#endif    
    call calc_relative_velocity_to_the_atmosphere(Vin, Xin, VdaGEI)
    call calc_air_drag(Cin,VdaGEI,rho_air,Aout)

    return
  end subroutine Air_Drag


  subroutine Jacchia_Roberts_Model(Nin, MJDN_in, DAFR_in, &
       alt_in, lon_in, lat_in, RA_sun_in, dec_sun_in, flag, rho_air)
    !--- args
    integer, intent(in)    :: Nin
    integer, intent(in)    :: MJDN_in
    real(8), intent(in)    :: DAFR_in
    real(8), intent(in)    :: alt_in
    real(8), intent(in)    :: lon_in
    real(8), intent(in)    :: lat_in
    real(8), intent(in)    :: RA_sun_in
    real(8), intent(in)    :: dec_sun_in
    integer, intent(inout) :: flag    
    real(8), intent(out)   :: rho_air
    !--- vars
    real(8) :: UT_MJDN1950
    real(8) :: RJUD
    real(8) :: OUTR
    real(8) :: SF(3)
    real(8) :: SA(3)
    real(8) :: SU(2)
    real(8) :: TE(2)
    real(8) :: ADRHO(6)
    real(8) :: WMOL    
    !--- body
    UT_MJDN1950 = MJDN_in - 33282.0d0
    
    if(MJDN_in < Nin)then
       RJUD = UT_MJDN1950
    else if(MJDN_in >= Nin)then
       RJUD = UT_MJDN1950 - 8035.0d0 ! go back to 22 years
    end if
    
    call SFDJ70(RJUD, DAFR_in, flag, SF, OUTR)

    SA(1) = lon_in
    SA(2) = lat_in
    SA(3) = alt_in

    SU(1) = RA_sun_in
    SU(2) = dec_sun_in
    
    call RSDAMO(SA, SU, SF, UT_MJDN1950, DAFR_in, &
         TE, ADRHO, WMOL, rho_air)

    return
  end subroutine Jacchia_Roberts_Model
  
  
  subroutine calc_relative_velocity_to_the_atmosphere(Vin, Xin, Vout)
    !***********************************************************
    ! This subroutine calcuates the relative velocity of a debri
    ! to the Earth's atmosphere.
    !***********************************************************
    !---- args
    real(8), intent(in)  :: Vin(3)
    real(8), intent(in)  :: Xin(3)
    real(8), intent(out) :: Vout(3)
    !***************************************************
    ! Vin  : Velocity Vector of the Debris
    ! Xin  : Position Vector of the Debris
    ! Vout : Rerative Velocity Vector to the Atomosphere
    !***************************************************
    !---- vars
    !---- body
    Vout(1) = Vin(1) + t0*omega_e(3)*Xin(2)
    Vout(2) = Vin(2) - t0*omega_e(3)*Xin(1)  
    Vout(3) = Vin(3)

    return
  end subroutine calc_relative_velocity_to_the_atmosphere
  
  
  subroutine calc_air_drag(Cin,Vin,rho,Aout)
    !---- args
    real(8), intent(in)  :: Cin
    real(8), intent(in)  :: Vin(3)
    real(8), intent(in)  :: rho
    real(8), intent(out) :: Aout(3)
    !---- vars
    real(8) :: Vin_nolm
    !----body
    Vin_nolm = sqrt(Vin(1)**2 + Vin(2)**2 + Vin(3)**2)    
    Aout(:) = Cin * rho * Vin_nolm * Vin(:)

    return
  end subroutine calc_air_drag
end module air_drag_mod
