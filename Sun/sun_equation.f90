module sun_mod
  use global_mod
  use subprog_mod
  implicit none
contains
  subroutine set_sun_attraction_constants
    keisuu1_sun = mu_s / mu_e
  end subroutine set_sun_attraction_constants
  
  
  subroutine set_solar_radiation_pressure_constants
    keisuu_ph = (2.0d0*Ad*WE*gammaE*Rse_mean**2)/(c*g0*Md*r0**2)
  end subroutine set_solar_radiation_pressure_constants

  
  subroutine sun_position
    real(8) :: lambda_M
    real(8) :: M_sun
    real(8) :: M2_sun
    real(8) :: lambda_ecliptic
    real(8) :: R_sun
    real(8) :: epsilon_sun
    
    lambda_M = 280.46d0 + 36000.771d0*T_J2000
    M_sun = (357.5277233d0*T_J2000)*deg2rad
    M2_sun = 2.0d0*M_sun
    lambda_ecliptic = (lambda_M + 1.914666471d0*sin(M_sun) &
         + 0.019994643d0*sin(M2_sun))*deg2rad
    R_sun = 1.000140612d0 - 0.016708617d0*cos(M_sun) &
         - 1.39589d-4*cos(M2_sun)
    epsilon_sun = (23.439291d0 - 0.0130042*T_J2000)*deg2rad

    XseGEI(1) = R_sun*cos(lambda_ecliptic)*AU2RE
    XseGEI(2) = R_sun*cos(epsilon_sun)*sin(lambda_ecliptic)*AU2RE
    XseGEI(3) = R_sun*sin(epsilon_sun)*sin(lambda_ecliptic)*AU2RE

    call Car2Sph(XseGEI,SSPHE)
  end subroutine sun_position

  
  subroutine distance_from_the_sun
    RseGEI = SSPHE(1)
    RseGEI_imverse = 1.0d0 / RseGEI
    XdsGEI = XdeGEI - XseGEI    
    RdsGEI = sqrt(XdsGEI(1)**2 + XdsGEI(2)**2 + XdsGEI(3)**2)
    RdsGEI_imverse = 1.0d0 / RdsGEI
  end subroutine distance_from_the_sun

  
  subroutine sun_attraction
    Asun = -keisuu1_sun * (RdsGEI_imverse**3 * XdsGEI + RseGEI_imverse**3 * XseGEI)
  end subroutine sun_attraction


  subroutine solar_radiation_pressure
    if(switch_shadow >= 1)then
       Aph = 0.0d0
       return
    end if

    Aph = keisuu_ph * RdsGEI_imverse**3 * XdsGEI
  end subroutine solar_radiation_pressure


  subroutine cylindrical_shadow
    !*****************************************************************!
    ! This subroutine determine if the flying object is in            !
    ! the shadow of the Earth or not. If the Sun is infinitely        !
    ! far from the Earth, so the light rays can be assumed parallel,  !
    ! producing a cylindrical Earth shadow of radius RE.              !
    ! This program use this assumpsion.                               !
    !                                                                 !
    ! 'switch_shadow = 0' : out of the shadow of the Earth            !
    ! 'switch_shadow = 1' : in the shadow of the Earth                !
    !*****************************************************************!
    real(8) :: zeta
    
    call angle(dble(-XseGEI),XdeGEI,zeta)
    if(zeta < asin(1.0d0/RdeGEI))then
       switch_shadow = 1
    else
       switch_shadow = 0
    end if
  end subroutine cylindrical_shadow
  

  subroutine shadow
    !************************************************************!
    ! This subroutine determine if the flying object is in       !
    ! the umbral of penumbral regions. The umbra is totally      !
    ! eclipsed by the Earth; the penumbra is only partially      !
    ! obscured by the Earth.                                     !
    !                                                            !
    ! Inputs : the position vectors of the object and the Sun    !
    ! Outputs : set the state variable, shadow                   !
    !           'switch_shadow = 1' : in penumbra                !
    !           'switch_shadow = 2' : in umbra                   !
    !                                                            !
    ! Reherence : Fundamentals of Astrodynamics and Applications !
    !                          ---Third Edition---               !
    !************************************************************!
    real(8), parameter :: alpha_pen = 4.695009689825d-3 ! [rad]
    real(8), parameter :: alpha_umb = 4.60974407268d-3  ! [rad]
    real(8) :: zeta
    real(8) :: SAT_horiz
    real(8) :: SAT_vert
    real(8) :: PEN_vert
    real(8) :: UMB_vert
    
    call angle(dble(-XseGEI),XdeGEI,zeta)
    
    SAT_horiz = RdeGEI*cos(zeta)
    SAT_vert = RdeGEI*sin(zeta)
    PEN_vert = RseGEI + tan(alpha_pen)*SAT_horiz
    if(SAT_vert <= PEN_vert)then
       switch_shadow = 1   
       UMB_vert = RseGEI - tan(alpha_umb)*SAT_horiz
    else if(SAT_vert <= UMB_vert)then
       switch_shadow = 2   
    else
       switch_shadow = 0   
    end if
  end subroutine shadow
end module sun_mod
