module moon_mod
  use global_mod
  implicit none
contains
  !**************************************************************!
  ! This subroutine calculate the position of the Moon using     !
  ! less precise formula from the astronoical almanac and Green. !
  ! This formula is based on Brown's theory that consists of     !
  ! a series expansion involving hundreds of terms.              !
  ! However,formula in the astronomical almanac includes only    !
  ! the few terms.                                               !
  !                                                              !
  ! Reference : Fundamentals of Astrodynamics and Applications   !
  !                          ---Third Edition---                 !
  !**************************************************************!
  subroutine set_moon_attraction_constants
      keisuu_moon = mu_m / mu_e
  end subroutine set_moon_attraction_constants

    
  subroutine moon_position    
    lambda_ecliptic = (218.32d0 + 481267.883d0*T_J2000 &
         + 6.29d0*sin(deg2rad*(134.9d0 + 477198.85d0*T_J2000)) &
         - 1.27d0*sin(deg2rad*(259.2d0 - 413335.38d0*T_J2000)) &
         + 0.66d0*sin(deg2rad*(235.7d0 + 890534.23d0*T_J2000)) &
         + 0.21d0*sin(deg2rad*(269.9d0 + 954397.70d0*T_J2000)) &
         - 0.19d0*sin(deg2rad*(357.5d0 + 35999.05d0*T_J2000)) &
         - 0.11d0*sin(deg2rad*(186.6d0 + 966404.05d0*T_J2000)))*deg2rad 
    
    phi_ecliptic = (5.13d0*sin(deg2rad*(93.3d0 + 483202.03d0*T_J2000)) &
         + 0.28d0*sin(deg2rad*(228.2d0 + 960400.87d0*T_J2000)) &
         - 0.28d0*sin(deg2rad*(318.3d0 + 6003.18d0*T_J2000)) &
         - 0.17d0*sin(deg2rad*(217.6d0 - 407332.20d0*T_J2000)))*deg2rad
    
    gamma_ecliptic = (0.9508d0 &
         + 0.0518d0*cos(deg2rad*(134.9d0 + 477198.85d0*T_J2000)) &
         + 0.0095d0*cos(deg2rad*(259.2d0 - 413335.38d0*T_J2000)) &
         + 0.0078d0*cos(deg2rad*(235.7d0 + 890534.23d0*T_J2000)) &
         + 0.0028d0*cos(deg2rad*(269.9d0 + 954397.70d0*T_J2000)))*deg2rad
    
    epsilon_ecliptic = (23.439291d0 - 0.0130042*T_J2000 &
         - 1.64d-7*T_J2000**2 + 5.04d-7*T_J2000**3)*deg2rad
    
    r_ecliptic = 1.0d0/sin(gamma_ecliptic)
    
    XmeGEI(1) = r_ecliptic*cos(phi_ecliptic)*cos(lambda_ecliptic)
    XmeGEI(2) = r_ecliptic*(cos(epsilon_ecliptic)*cos(phi_ecliptic)*sin(lambda_ecliptic) &
         - sin(epsilon_ecliptic)*sin(phi_ecliptic))
    XmeGEI(3) = r_ecliptic*(sin(epsilon_ecliptic)*cos(phi_ecliptic)*sin(lambda_ecliptic) &
         + cos(epsilon_ecliptic)*sin(phi_ecliptic))
  end subroutine moon_position
  
  
  subroutine distance_from_the_moon
    XdmGEI = XdeGEI - XmeGEI
    RmeGEI = sqrt(XmeGEI(1)**2 + XmeGEI(2)**2 + XmeGEI(3)**2)
    RdmGEI = sqrt(XdmGEI(1)**2 + XdmGEI(2)**2 + XdmGEI(3)**2)
  end subroutine distance_from_the_moon
  
  
  subroutine attraction_of_the_moon
    Amoon = -keisuu_moon * (XdmGEI/RdmGEI**3 + XmeGEI/RmeGEI**3)
  end subroutine attraction_of_the_moon
end module moon_mod
