module global_mod
  use charge_solve_mod, only : initialize_IRI
  implicit none
  !=======Physical Constants=======!
  real(8), save :: G              ! Newtonian Constant of Gravity
  real(8), save :: mu_e           ! Gravitational parameter of the Earth : G*Me [m^3/s^2]
  real(8), save :: mu_s           ! Gravitational parameter of the Sun 
  real(8), save :: mu_m           ! Gravitational parameter of the Moon
  real(8), save :: Rsun           ! Radius of the Sun [m]
  real(8), save :: omega_e(3)     ! the Earth's mean rotation rate on its axis [rad/sec] 
  real(8), save :: eps0           ! Dielectric Constant of vacuum
  real(8), save :: c              ! The Speed of Light
  real(8), save :: Rse_mean       ! Mean Distance between the Earth and the Sun
  real(8), save :: Rme_mean       ! Mean Distance between the Earth and the Moon
  real(8), save :: r0             ! Radius of the Earth
  real(8), save :: g0             ! Gravitational acceleration [m/s^2]
  real(8), save :: v0             ! Normalization Constant of velocity [m/s]
  real(8), save :: t0             ! Normalization Constant of time [sec]
  real(8), save :: pi             ! Circular Constant
  real(8), save :: pi2            ! 2*pi
  real(8), save :: pi_half        ! pi/2
  real(8), save :: RE2km          ! convert an unit from RE to km
  real(8), save :: km2RE          ! convert an unit from km to RE
  real(8), save :: r0_geopack     ! convert an unit to use Geopack-2008
  real(8), save :: rad2deg        ! unit converter : radian into degree
  real(8), save :: rad2lt         ! unit converter : radian into localtime [hour]
  real(8), save :: deg2rad        ! unit converter : degree into radian
  real(8), save :: sec2rad        ! unit converter : second (hour angle) into radian
  real(8), save :: min2rad        ! unit converter : minute (hour angle) into radian
  real(8), save :: AU             ! Astronomical Unit [m]
  real(8), save :: AU2RE          ! unit converter : AU into RE
  real(8), save :: charge_d              ! Charge Amount of the debri
  real(8), save :: into_day       ! unit converter : normalized time to day
  
  !=======Setting=======!
  !*******************************************************
  ! Unit number 20 to 49 for the moon's position data file
  ! Unit number 50 to 79 for the sun's position data file
  !*******************************************************
  integer, parameter :: fi1     = 10
  integer, parameter :: fi_deb  = 11
  integer, parameter :: fi_flux = 12
  integer, parameter :: fo1     = 100
  integer, parameter :: fo2     = 101
  integer, parameter :: fo3     = 102
  integer, parameter :: fo4     = 103
  integer, parameter :: fo5     = 104
  integer, parameter :: fo6     = 105
  integer, parameter :: fo7     = 106
  integer, parameter :: fo8     = 107
  integer, parameter :: fo9     = 108
  integer, parameter :: fo10    = 109
  integer, parameter :: fo11    = 110
  integer, parameter :: fo12    = 120
  integer, parameter :: fo13    = 130
  integer, parameter :: fo14    = 140
  integer, parameter :: fo15    = 150
  integer, parameter :: fi_moon = 24
  integer, parameter :: fi_sun  = 25
  integer, save      :: N_ABM            ! the order of Adams-Bashforth-Moulton
  integer, save      :: N_debris         ! the number of debris
  integer, save      :: i_step
  integer(8), save   :: i_max
  real(8), allocatable, save :: table_debris(:,:)
  real(8), save      :: dt_imverse
  real(8), save      :: dt               ! time step
  real(16), save     :: simulation_span  ! simulation span in days [day]
  real(8), save      :: simulation_year 
  real(8), save      :: simulation_day
  real(8), save      :: limit_h
  real(8), save      :: limit_hGEI
  real(8), save      :: deltaXYZ         ! Stepsize of Earth's electric potential
  real(4), save      :: VGSEX, VGSEY, VGSEZ
  
  !=======Time=======!
  !**********************************************************************
  ! Table as below shows what each variable means
  !
  ! time      : elapsed time from the epoch
  ! TD_year, TD_month, TD_day, TD_hour, TD_min, TD_sec ... Dynamical Time
  ! TD_JD     : Julian Day in Dynamical Time
  ! TD_DY     : Day of the Year in Dynamical Time scale
  ! TD_DAFR   : a fraction of a day in sec (0 to 86400)
  !
  ! UT_year, UT_month, UT_day, UT_hour, UT_min, UT_sec ... Universal Time   
  ! UT_JD     : Julian Day in Universal Time
  !             (epoch is January 1, 4713 B.C, 12:00)  
  ! UT_MJD    : Modified Julian Day in Universal Time
  ! UT_MJDN   : Modified Julian Day Number
  ! UT_DAFR   : a fraction of a day in sec (0 to 86400)
  ! UT_YEFR   : a fraction of a year in year
  ! UT_DoY    : the number of days elapsed from the beginning of the year
  ! UT_JC2000 : Julian Century in Universal Time elapsed from J2000.0
  ! UT_JC1950 : Julian Century in Dynamical Time elapsed from B1950.0
  ! T_GST     : Greenwich sidereal time at a given time UT
  ! T_GMST_0h : Greenwich Mean Sidereal Time of 0h UT at a given date(rad)
  !**********************************************************************    
  integer, save :: UT_year
  integer, save :: UT_month
  integer, save :: UT_day
  integer, save :: UT_hour
  integer, save :: UT_min
  real(8), save :: UT_sec
  real(8), save :: TD_JD
  real(8), save :: TD_JC2000
  real(8), save :: UT_JD   
  real(8), save :: UT_MJD
  integer, save :: UT_MJDN
  real(8), save :: UT_DAFR
  real(8), save :: UT_YEFR
  integer, save :: UT_DoY
  real(8), save :: UT_JC2000
  real(8), save :: UT_MJDN1950
  real(8), save :: UT_aphelion
  real(8), save :: T_GMST_0h
  real(8), save :: T_GST
  integer, save :: TD_MJDN
  real(8), save :: TD_DAFR
  real(8), save :: delta_T
  real(8), save :: time
  real(8), save :: FOY
  
  !=======Data of a debri=======!
  real(8), save :: XdeGEI_initial(3)
  real(8), save :: VdeGEI_initial(3)
  real(8), save :: Rd             ! Characteristic length of a debri [m]
  real(8), save :: Ad             ! Cross-sectional area of a debri [m^2]
  real(8), save :: Md             ! Mass of a debri [kg]
  real(8), save :: Ud             ! Surface potential of a debri
  
  !=======Position, Velocity, Acceleration=======!
  real(8), save :: X(3)           ! Position of a debri [normalized]
  real(8), save :: V(3)           ! Velocity of a debri [normalized]
  real(8), save :: K(6)           ! Keplerian Elements of a debri
  !*******Warning*****************************
  ! Make sure the direction of vectors
  ! before calculations of some purterbations.
  !*******************************************
  real(8), save :: XdeGEI(3)      ! state vector of a debri referred to J2000.0 (normalized)
  real(4), save :: XdeGEO(3)      ! position vector of a debri in the geodetic coordinate (normalized)
  real(8), save :: VdeGEI(3)      ! Velocity of a debri (normalized)  
  real(4), save :: alt_gd         ! height above the reference ellipsoid (km)
  real(4), save :: lat_gd         ! geodetic latitude   (-pi/2 - pi/2 [rad])
  real(8), save :: lat_gc         ! geocentric latitude (-pi/2 - pi/2 [rad])
  real(8), save :: colat_gc       ! geocentric colatitude (0 - pi [rad])
  real(8), save :: colat_gd       ! geodetic colatitude (0 - pi [rad])
  real(4), save :: lon_gd         ! geodetic longitude  (0 - 2pi [rad])
  real(8), save :: lon_gc         ! geocentric longitude, measured positively east from the direction of vernal equinox
  real(8), save :: XseGEI(3)      ! position vector of the Sun       (normalized)
  real(8), save :: XdsGEI(3)      ! position of a debri from the Sun (normalized)
  real(8), save :: dec_sun        ! Declination of the Sun (-pi/2 to pi/2 [rad])
  real(8), save :: RA_sun         ! Right Ascension of the Sun (0 to 2pi [rad])
  real(8), save :: XmeGEI(3)      ! Mean Moon position vector referred to J2000.0 (normalized)
  real(8), save :: XdmGEI(3)      ! Mean Position vector of a debri from the Moon (normalized)
  real(8), save :: RdeGEI         ! Distance between barycenter of the Earthand a debri (normalized)
  real(8), save, allocatable :: RdeGEI_inv(:) ! Inverse of RdeGEI (normalized)
  real(8), save :: RdsGEI         ! Distance between the barycenter of the Sun and a debri (normalized)
  real(8), save :: RdsGEI_inv(3)  ! Inverse of RdsGEI (normalized)
  real(8), save :: RseGEI         ! Distance of the barycenter between the Earth and the Sun (normalized)
  real(8), save :: RseGEI_inv(3)  ! Inverse of RseGEI (normalized)
  real(8), save :: RmeGEI         ! Distance of the barycenter between the Earth and the Moon (normalized)
  real(8), save :: RdmGEI         ! Distance between the barycenter of the Moon and a debri (normalized)
  real(8), save :: RmeGEI_inv(3)  ! Inverse of RmeGEI (normalized)
  real(8), save :: RdmGEI_inv(3)  ! Inverse of RdmGEI (normalized)
  real(8), save :: Atotal(3)      ! Total Perturbation acting on a debri (normalized)
  real(8), save :: Aearth(3)      ! Acceleration derived from the gravity of the Earth (normalized)
  real(8), save :: Ageo(3)        ! Force derived from high-order geopotential (normalized)
  real(8), save :: Aair(3)        ! Acceleration derived from the atmospheric drag (normalized)
  real(8), save :: Amoon(3)       ! Acceleration derived from the attraction of the Moon (normalized)
  real(8), save :: Asun(3)        ! Acceleration derived from the attraction of the Sun (normalized)
  real(8), save :: Aph(3)         ! Acceleration derived from the solar radiation pressure (normalized)
  real(8), save :: Ecr(3)         ! Corotation electric field vector [N/C]
  real(4), save :: Ecn(3)         ! Convection electric field vector [N/C)
  real(8), save :: VcrossB(3)     ! V cross B, V : Velocity of a debri, B : Magnetic field vector
  real(8), save :: BGEI(3)        ! Magnetic field vector in GEI coordinate [T]
  real(8), save :: BGEO(3)        ! Magnetic field vector in GEO coordinate [T]
  
  !*************
  ! Geopotential
  !*************
  real(8), save, allocatable :: C_geo(:,:)
  real(8), save, allocatable :: S_geo(:,:)
  integer, save :: N_geo
  
  !*********
  ! Air Drag
  !*********
  integer, save :: N_air
  real(8), save :: keisuu_airdrag
  real(8), save :: Cd             ! Coefficient of the Air Drag
  integer, save :: air_frag

  !*****
  ! Moon
  !*****
  real(8), save :: keisuu_moon
  real(8), save :: moon_interval  
  character(40), allocatable, save :: file_moon(:) ! ephemeris file name
  integer, save :: N_file ! the number of ephemeris data file
  integer, save :: N_data_moon ! the number of ephemeris data
  integer, save :: index_moon  
  integer, save :: initial_index_moon
  real(8), allocatable, save :: TD_JD_moon(:)
  real(8), allocatable, save :: posi_table_moon(:,:)  

  !***************
  ! Solar Pressure
  !***************
  real(8), save :: keisuu_ph
  real(8), save :: eps_r
  real(8), save :: D_aphelion

  !****
  ! Sun
  !****
  real(8), save :: keisuu_sun
  integer, save :: N_data_sun ! the number of ephemeris data
  real(8), save :: sun_interval
  integer, save :: index_sun
  integer, save :: initial_index_sun
  character(40), allocatable, save :: file_sun(:) ! ephemeris file name
  real(8), allocatable, save :: TD_JD_sun(:)
  real(8), allocatable, save :: posi_table_sun(:,:)  

  !*********************************
  ! Convection Electric Field Vector
  !*********************************
  real(8), save :: XdeltaGEI(3)
  real(8), save :: YdeltaGEI(3)
  real(8), save :: ZdeltaGEI(3)
contains    
  subroutine set_physical_constants
    !---- args
    !---- vars
    !---- body
    g0         = mu_e / r0**2
    v0         = sqrt(r0 * g0)
    t0         = sqrt(r0 / g0)
    pi         = acos(-1.0d0)
    pi2        = 2.0d0*pi
    pi_half    = 0.5d0*pi
    RE2km      = r0*1.0d-3
    km2RE      = 1.0d3/r0
    r0_geopack = r0/6.3712d6
    rad2deg    = 180.0d0 / pi
    rad2lt     = 12.0d0 / pi
    deg2rad    = pi / 180.0d0
    sec2rad    = deg2rad / 3600.0d0
    AU2RE      = AU / r0
    
    into_day        = t0 / 86400.0d0
    dt              = 1.0d0 / dt_imverse
    simulation_span = 365.0d0*simulation_year + simulation_day
    i_max           = int((86400.0d0*simulation_span)/(dt*t0))
    limit_hGEI      = (limit_h*1.0d3 + r0) / r0

    call initialize_IRI

    return
  end subroutine set_physical_constants

  
  subroutine set_clear
    !***********************************************
    ! Unit of the charge amount is Coulomb [C = A*s]
    !***********************************************
    !---- args
    !---- vars
    !---- body    
    charge_d          = 4.0d0 * pi * eps0 * Rd * Ud
    
    X(:)       = XdeGEI_initial(:)
    V(:)       = VdeGEI_initial(:)

    index_sun  = initial_index_sun
    index_moon = initial_index_moon
    air_frag   = 0
    
    time       = 0.0d0

    alt_gd        = 0.0d0
    lat_gd        = 0.0d0
    lat_gc        = 0.0d0
    colat_gc      = 0.0d0
    colat_gd      = 0.0d0
    lon_gd        = 0.0d0
    lon_gc        = 0.0d0

    dec_sun       = 0.0d0
    RA_sun        = 0.0d0

    VdeGEI(:)     = 0.0d0
    XdeGEI(:)     = 0.0d0
    RdeGEI        = 0.0d0
    RdeGEI_inv(:) = 0.0d0

    XdsGEI(:)     = 0.0d0    
    RdsGEI        = 0.0d0
    RdsGEI_inv(:) = 0.0d0
    
    XseGEI(:)     = 0.0d0    
    RseGEI        = 0.0d0
    RseGEI_inv(:) = 0.0d0

    XmeGEI(:)     = 0.0d0    
    RmeGEI        = 0.0d0
    RmeGEI_inv(:) = 0.0d0
    
    XdmGEI(:)     = 0.0d0    
    RdmGEI        = 0.0d0
    RdmGEI_inv(:) = 0.0d0
    
    Atotal(:)     = 0.0d0
    Aearth(:)     = 0.0d0
    Ageo(:)       = 0.0d0
    Aair(:)       = 0.0d0
    Amoon(:)      = 0.0d0
    Asun(:)       = 0.0d0
    Aph(:)        = 0.0d0
    Ecr(:)        = 0.0d0
    Ecn(:)        = 0.0d0
    VcrossB(:)    = 0.0d0
    BGEI(:)       = 0.0d0
    
    return
  end subroutine set_clear
end module global_mod
