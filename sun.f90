module sun_mod
  use global_mod, only : &
       pi, pi_half, c, g0, r0, eps_r, fi_sun, &
       TD_JD_sun, posi_table_sun
  use subprog_mod, only : &
       angle, interpolation, nolm
  use conversion_coordinates_mod, only : &
       precession, &
       set_precession_matrix, &
       Car2Sph
  implicit none
contains
  subroutine set_sun_attraction_constants(mu_s, mu_e, Cout)
    !---- args
    real(8), intent(in)  :: mu_s ! Gravity Constant of the Sun
    real(8), intent(in)  :: mu_e ! Gravity Constant of the Earth
    real(8), intent(out) :: Cout
    !---- vars
    !---- body
    Cout = mu_s / mu_e

    return
  end subroutine set_sun_attraction_constants
  
    
  subroutine set_solar_radiation_pressure_constants(Ad, Md, Cout)
    !---- args
    real(8), intent(in)  :: Ad ! Reference Area
    real(8), intent(in)  :: Md ! Mass
    real(8), intent(out) :: Cout
    !---- vars
    !---- body
    Cout = (Ad*(1.0d0 + eps_r))/(c*g0*Md)

    return
  end subroutine set_solar_radiation_pressure_constants


  subroutine load_sun_ephemeris_table &
       (Tin, num_file, filename, num_data, interval)
    !**********************************
    ! Store all ephemeris data in array
    !**********************************
    !--- args
    real(8), intent(in)          :: Tin
    integer, intent(in)          :: num_file
    character(len=*), intent(in) :: filename(1:num_file)
    integer, intent(out)         :: num_data
    real(8), intent(out)          :: interval
    !--- vars
    integer :: idx
    integer :: iost
    integer :: i
    real(8) :: sze,cze,sz,cz,st,ct
    real(8) :: XseKM(3)
    real(8) :: JD_sun
    real(8) :: sun_interval1
    real(8) :: sun_interval2    
    !--- body
    call set_precession_matrix(Tin,sze,cze,sz,cz,st,ct)
    call count_sun_ephemeris_file(num_file, filename(:), num_data)
    allocate(TD_JD_sun(num_data))
    allocate(posi_table_sun(1:3,1:num_data))

    idx = 1
    ! read ephemeris data recorded in each file
    do i = 1, num_file
       open(fi_sun, file = filename(i), status = 'old')
       rewind(fi_sun)
       do while(.true.)
          read(fi_sun, *, iostat = iost) JD_sun, XseKM(1), XseKM(2), XseKM(3)
          if(iost /= 0)then
             exit ! finish counting if all data are read
          end if
          ! add to entry arrangements
          TD_JD_sun(idx) = JD_sun
          call precession(XseKM(1:3)*1.0d3/r0,sze,cze,sz,cz,st,ct,&
               posi_table_sun(1:3, idx))
          idx = idx + 1
       end do
       close(fi_sun)
    end do

    print*, "*****************************************************"
    print*, "Conversion of the Sun's ephemeris is completed,"
    print*, "considering the Precession of the Earth."
    print*, "-----------------------------------------------------"
    print*, ""
    print*, "Julian Century refered to J2000.0(DT)"
    print*, Tin
    print*, ""
    print*, "*****************************************************"        

    ! check if time interval is a constant
    do i = 1, num_file, 2
       sun_interval1 = TD_JD_sun(i+1) - TD_JD_sun(i)
       sun_interval2 = TD_JD_sun(i+2) - TD_JD_sun(i+1)
       if(sun_interval1 - sun_interval2 > 1.0d-5)then
          print*, "Error : Time Interval of the Sun's Ephemeris is not a Constant"
          print*, "Stop at sun.f90"
          stop
       end if
    end do
    interval = TD_JD_sun(2) - TD_JD_sun(1)

    return
  end subroutine load_sun_ephemeris_table


  subroutine count_sun_ephemeris_file(num_file, file_name, count)
    !***************************************
    ! Count the sum number of ephemeris data
    !***************************************
    !--- args
    integer,            intent(in)    :: num_file
    character(len = *), intent(in)    :: file_name(:)
    integer,            intent(inout) :: count
    !--- vars
    integer :: iost
    integer :: num
    integer :: i
    real(8) :: buff
    !--- body
    
    num = 0
    do i = 1, num_file
       ! count the number of tabular row recorded in each file
       open(fi_sun, file = file_name(i), status = 'old')
       rewind(fi_sun)
       do while(.true.)
          read(fi_sun, '(A255)', iostat = iost) buff
          if(iost /= 0)then
             exit
          end if
          num = num + 1
       end do
       close(fi_sun)
    end do
    count = num

    return
  end subroutine count_sun_ephemeris_file

  
  subroutine find_sun_ephemeris_JD(Tin, num_data, index)
  !************************************************************
  ! Find the closest index to the given JD from ephemeris table
  !************************************************************    
    !--- args
    real(8), intent(in)  :: Tin
    integer, intent(in)  :: num_data
    integer, intent(out) :: index
    !--- vars
    integer :: id_min
    integer :: id_max
    integer :: id_mid
    real(8) :: jd_mid
    !--- body

    ! check the input range
    if(Tin < TD_JD_sun(1))then
       print*, 'Error : Ephemeris Table is INVALID for the Simulation Epoch'
       stop
    end if
    if(Tin > TD_JD_sun(num_data))then
       print*, 'Error : Ephemeris Table is INVALID for the Simulation Epoch'
       stop
    end if

    id_min = 1
    id_max = num_data

    ! find the closest index to given JD by half deviding method
    do while(id_min < id_max)
       id_mid = id_min + (id_max - id_min)/2
       jd_mid = TD_JD_sun(id_mid)
       if(jd_mid == Tin)then
          index  = id_mid
          exit
       else if(id_max - id_min == 1)then ! case of achieving the purpose
          index  = id_min
          exit
       else if(Tin < jd_mid)then
          id_max = id_mid
       else if(Tin > jd_mid)then
          id_min = id_mid
       else
          print*, 'Error : Input the unknown number on TD_JD'
          print*, 'Stop at sun.f90'
          stop
       end if
    end do

    return
  end subroutine find_sun_ephemeris_JD


  subroutine calc_sun_position(Tin, interval, index, Xout)
    !*********************************************
    ! Choose the 3 tabular data points, and decide
    ! the position vector of the Sun
    !*********************************************
    !---- args
    real(8), intent(in)    :: Tin      ! Julian Day (TD)
    real(8), intent(in)    :: interval ! Ephemeris Data Interval
    integer, intent(inout) :: index    !
    real(8), intent(out)   :: Xout(3)  ! Position Vector of the Sun
    !---- vars
    real(8) :: XseGEI3(3), XseGEI2(3), XseGEI1(3)
    !---- body
    if(abs(Tin - TD_JD_sun(index)) &
         > abs(Tin - TD_JD_sun(index+1)))then
       index = index + 1
    end if
    
    XseGEI3(:) = posi_table_sun(:,index+1)
    XseGEI2(:) = posi_table_sun(:,index)
    XseGEI1(:) = posi_table_sun(:,index-1)
    
    call interpolation(XseGEI1(:),XseGEI2(:),XseGEI3(:), &
         Tin, TD_JD_sun(index), interval, Xout(:))

    return
  end subroutine calc_sun_position    

  
  subroutine distance_from_the_sun1(Xse, Rse, Rse_inv)
    !---- args
    real(8), intent(in)  :: Xse(3)
    real(8), intent(out) :: Rse
    real(8), intent(out) :: Rse_inv(3)
    !---- vars
    !---- body
    Rse        = nolm(Xse)
    Rse_inv(1) = 1.0d0 / Rse
    Rse_inv(2) = Rse_inv(1) * Rse_inv(1)
    Rse_inv(3) = Rse_inv(1) * Rse_inv(2)
        
    return
  end subroutine distance_from_the_sun1


  subroutine distance_from_the_sun2 &
       (Xde, Xse, Xds, Rds, Rds_inv)
    !---- args
    real(8), intent(in)  :: Xde(3)
    real(8), intent(in)  :: Xse(3)
    real(8), intent(out) :: Xds(3)
    real(8), intent(out) :: Rds
    real(8), intent(out) :: Rds_inv(3)
    !---- vars
    !---- body
    Xds(:)     = Xde(:) - Xse(:)
    Rds        = nolm(Xds)
    Rds_inv(1) = 1.0d0 / Rds
    Rds_inv(2) = Rds_inv(1)*Rds_inv(1)
    Rds_inv(3) = Rds_inv(1)*Rds_inv(2)
    
    return
  end subroutine distance_from_the_sun2

  
  subroutine sun_attraction &
       (Cin, Xds, Rds_inv3, Xse, Rse_inv3, Aout)
    !---- args
    real(8), intent(in) :: Cin
    real(8), intent(in) :: Xds(3)
    real(8), intent(in) :: Rds_inv3
    real(8), intent(in) :: Xse(3)
    real(8), intent(in) :: Rse_inv3
    real(8), intent(out) :: Aout(3)
    !---- vars
    !---- body
    Aout(:) = -Cin*(Rds_inv3*Xds(:) + Rse_inv3*Xse(:))

    return
  end subroutine sun_attraction


  subroutine solar_radiation_pressure &
       (Cin, D_aph, Rsun, Xde, Xds, Xse, Rde, Rds, &
       Rde_inv1, Rds_inv1, Aout,idx)
    !***********************************************************************
    ! This subroutine calculates the solar radiation pressure which varies
    ! with the amount of Sun's shadowing ratio by the Earth.
    !
    ! Input  : Parameter 'switch_shadow' which denote the shadow region
    !          Debri's position referred to the Sun (XdsGEI)
    !          Distance between the Sun and a debri (RdsGEI)
    ! Output : The acceleration caused by the solar radiation pressure (Aph)
    ! Ref    :                Toshimichi Otsubo et al
    !          Error Control of Numerical Integration in SLR Analysis
    !                            Software CONCERTO
    ! Author : Keisuke Akari
    ! Date   : Novenber 28, 2016
    !***********************************************************************
    !---- args
    real(8), intent(in)  :: Cin
    real(8), intent(in)  :: D_aph
    real(8), intent(in)  :: Rsun
    real(8), intent(in)  :: Xde(3)
    real(8), intent(in)  :: Xds(3)
    real(8), intent(in)  :: Xse(3)
    real(8), intent(in)  :: Rde
    real(8), intent(in)  :: Rds
    real(8), intent(in)  :: Rde_inv1
    real(8), intent(in)  :: Rds_inv1
    real(8), intent(out) :: Aout(3)
    integer, intent(out) :: idx
    !*************************************
    ! Cin       :
    ! D_aph     :
    ! Rsun      :
    ! Xde(3)    :
    ! Xds(3)    :
    ! Xse(3)    :
    ! Rde       :
    ! Rds       :
    ! Rde_inv1  :
    ! Rds_inv1  :
    ! id_shadow :
    ! Aout(3)   :
    ! idx       : shadow_function
    !             0: sunlit
    !             1: 
    !             2: ecllipse
    !*************************************
    !---- vars
    real(8) :: SF
    real(8) :: shadow_ratio
    !---- body
    call set_solar_flux(D_aph, SF)
#ifdef cylindrical_shadow_model    
    call cylindrical_shadow(Xse, Xde, Rde_inv1, idx)
#endif
#ifdef conical_shadow_model
    call conical_shadow(Xse, Xde, Rde, idx)
#endif    
    
    if(idx == 0)then
       Aout(:) = SF * Cin * Rds_inv1 * Xds(:)
       return
    else if(idx == 2)then
       Aout(:) = 0.0d0
       return
    else if(idx == 1)then
       call calc_shadowing_ratio &
            (Rsun, Xds, Rds, Xde, Rde, shadow_ratio)
       Aout(:) = shadow_ratio * SF * Cin * Rds_inv1 * Xds(:)

       return
    end if
  end subroutine solar_radiation_pressure


  subroutine set_solar_flux(D_aph, SF)
    !---- args
    real(8), intent(in)  :: D_aph 
    real(8), intent(out) :: SF
    !*********************************************************
    ! D_aph : Rotation Angle mearsured from the Aphelion [rad]
    ! SF    : Solar Flux []
    !*********************************************************
    !---- vars
    !---- body
    SF = 1.358d3/(1.004d0 + 3.34d-2*cos(D_aph))

    return
  end subroutine set_solar_flux


  subroutine cylindrical_shadow(Xse, Xde, Rde_inv, idx)
    !***************************************************************
    ! This subroutine determine if the flying object is in            
    ! the shadow of the Earth or not. If the Sun is infinitely        
    ! far from the Earth, so the light rays can be assumed parallel,  
    ! producing a cylindrical Earth shadow of radius RE.              
    ! This program use this assumpsion.                               
    !                                                                 
    ! Input  : Positions of the debri and the Sun (XdeGEI, XseGEI)
    ! Output : set the parameter, 'switch_shadow'
    !          'idx = 0' : sunlit
    !          'idx = 2' : eclipse
    ! Ref    :                 David.A.Vallado
    !          Fundamentals of Astrodynamics and Applications
    !                       -Third Edition-, p304
    ! Author : Keisuke Akari
    ! Date   : October 25, 2016
    !***************************************************************
    !---- args
    real(8), intent(in)  :: Xse(3)
    real(8), intent(in)  :: Xde(3)
    real(8), intent(in)  :: Rde_inv
    integer, intent(out) :: idx
    !******************************************************
    ! Xse(3)  : Sun's Position Vector from the Earth [RE]
    ! Xde(3)  : Debri's Position Vector from the Earth [RE]
    ! Rde     : Distance between Debris and the Earth [RE]
    ! Rde_inv : Inverse of Rde
    ! idx     : State Parameter
    !******************************************************
    !---- vars
    real(8) :: zeta
    !---- body
    call angle(dble(-Xse),Xde,zeta)
    if(zeta < asin(Rde_inv))then
       idx = 2
    else
       idx = 0
    end if

    return
  end subroutine cylindrical_shadow
  

  subroutine conical_shadow(Xse, Xde, Rde, idx)
    !*********************************************************
    ! This subroutine determine if the flying object is in       
    ! the umbral or penumbral regions. The umbra is totally      
    ! eclipsed by the Earth; the penumbra is only partially      
    ! obscured by the Earth.                                     
    !                                                            
    ! Input  : the position vectors of the object and the Sun   
    ! Output : set the state parameter, 'switch_shadow'
    !           'switch_shadow' = 0 : sunlit 
    !           'switch_shadow' = 1 : penumbra                
    !           'switch_shadow' = 2 : umbra                   
    ! Ref    :                 David.A.Vallado
    !          Fundamentals of Astrodynamics and Applications 
    !                     -Third Edition-, p300-303
    ! Author : Keisuke Akari
    ! Date   : October 25, 2016
    !********************************************************
    !---- args
    real(8), intent(in)  :: Xse(3)
    real(8), intent(in)  :: Xde(3)
    real(8), intent(in)  :: Rde
    integer, intent(out) :: idx
    !*****************************************************
    ! Xse(3) : Sun's Position Vector from the Earth [RE]
    ! Xde(3) : Debris' Position Vector from the Earth [RE]
    ! Rde    : Distance between Debris and the Earth [RE]
    ! idx    : State Parameter
    !*****************************************************
    !---- vars
    real(8), parameter :: alpha_pen = 4.695009689825d-3 ! [rad]
    real(8), parameter :: alpha_umb = 4.60974407268d-3  ! [rad]
    real(8) :: zeta
    real(8) :: SAT_horiz
    real(8) :: SAT_vert
    real(8) :: PEN_vert
    real(8) :: UMB_vert
    !---- body
    call angle(-Xse,Xde,zeta)
    if(zeta >= pi_half)then
       idx = 0
       return
    end if
    
    SAT_horiz = Rde*cos(zeta)
    SAT_vert  = Rde*sin(zeta)
    PEN_vert  = 1.0d0 + tan(alpha_pen)*SAT_horiz
    if(SAT_vert <= PEN_vert)then
       idx = 1
       UMB_vert = 1.0d0 - tan(alpha_umb)*SAT_horiz
       if(SAT_vert <= UMB_vert)then
          idx = 2
       end if
    else
       idx = 0
    end if

    return
  end subroutine conical_shadow


  subroutine calc_shadowing_ratio(radius, Xds, Rds, Xde, Rde, ratio)
    !************************************************************
    ! Warning : 
    !   This subroutine sometimes output NAN because arguments of 
    !   acos become larger than 1 or less than -1.
    !************************************************************
    !---- args
    real(8), intent(in)  :: radius 
    real(8), intent(in)  :: Xds(3) 
    real(8), intent(in)  :: Rds
    real(8), intent(in)  :: Xde(3)
    real(8), intent(in)  :: Rde
    real(8), intent(out) :: ratio
    !*******************************************************
    ! radius : Radius of the Sun []
    ! Xds(3) : Position Vector of Debris from the Sun [RE]
    ! Rds    : Distance between Debris and the Sun [RE]
    ! Xde(3) : Position Vector of Debris from the Earth [RE]
    ! Rde    : Distance between Debris and the Earth [RE]
    ! ratio  : Sun's Visible Ratio from Debris
    !*******************************************************
    !---- vars
    real(8) :: alpha
    real(8) :: beta
    real(8) :: gamma
    real(8) :: delta
    real(8) :: theta1
    real(8) :: theta2
    !---- body
    alpha   = asin(radius/(r0*Rds))
    beta    = asin(1.0d0/Rde)
    call angle(-Xde, -Xds, gamma)
    delta   = (beta/alpha)**2

    theta1  = acos((alpha**2 + gamma**2 - beta**2)/(2.0d0*alpha*gamma))
    theta2  = acos((beta**2  + gamma**2 - alpha**2)/(2.0d0*beta*gamma))
    
    ratio = 1.0d0 - (theta1 + theta2*delta - gamma/alpha*sin(theta1))/pi

    return
  end subroutine calc_shadowing_ratio
end module sun_mod
