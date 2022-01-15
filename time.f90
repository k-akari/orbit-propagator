module time_mod
  use global_mod, only : &
       pi2, sec2rad, t0, time, &
       UT_DAFR, UT_MJDN,  UT_JD,  UT_JC2000, UT_aphelion, UT_MJD, &
       TD_DAFR, TD_MJDN,  TD_JD,  TD_JC2000, D_aphelion, &
       UT_year, UT_month, UT_day, UT_hour,   UT_min, UT_sec, UT_DoY, &
       UT_YEFR, delta_T, T_GMST_0h, T_GST, FOY
  implicit none
contains
  !**********************************************************************
  ! Table as below shows what each variable means
  !
  ! time        : elapsed time from the simulation epoch (normalized)
  ! TD_year, TD_month, TD_day, TD_hour, TD_min, TD_sec ... Dynamical Time
  ! TD_DY       : Day of the Year in Dynamical Time scale
  ! TD_DAFR     : a fraction of a day in sec (0 to 86400)
  !
  ! UT_year, UT_month, UT_day, UT_hour, UT_min, UT_sec ... Universal Time   
  ! UT_JD       : Julian Day in Universal Time
  !             (epoch is January 1, 4713 B.C, 12:00)  
  ! UT_MJD      : Modified Julian Day in Universal Time
  ! UT_MJDN     : Modified Julian Day Number
  ! UT_DAFR     : a fraction of a day in sec (0 to 86400)
  ! UT_YEFR     : Date described in year of Anno Domini
  ! UT_JC2000   : Julian Century in Universal Time elapsed from J2000.0
  ! UT_JC1950   : Julian Century in Dynamical Time elapsed from B1950.0
  ! T_GMST_0h   : The sidereal time at Greenwich at 0h of a given date
  ! T_GST       : Greenwich Mean Sidereal Time at the moment [rad]
  ! UT_aphelion : Elapsed day number from the aphelion of the Earth.
  ! D_aphelion  : Rotated angle of the Earth, measured positively east
  !               from the aphelion.
  !**********************************************************************
  subroutine renew_clock(t)
    !---- args
    real(8), intent(in) :: t
    !---- vars
    !---- body
    time = time + t

    UT_DAFR = UT_DAFR + t*t0
    if(86400.0d0 <= UT_DAFR)then
       UT_DAFR = UT_DAFR - 86400.0d0
       UT_MJDN = UT_MJDN + 1.0d0
       call set_time2(1)
    else if(UT_DAFR < 0.0d0)then
       UT_DAFR = UT_DAFR + 86400.0d0
       UT_MJDN = UT_MJDN - 1.0d0
       call set_time2(-1)
    else if(0.0d0 <= UT_DAFR .and. UT_DAFR < 86400.0d0)then
       call set_time1
    end if
    
    TD_DAFR = TD_DAFR + t*t0
    if(TD_DAFR >= 86400.0d0)then
       TD_DAFR = TD_DAFR - 86400.0d0
       TD_MJDN = TD_MJDN + 1
    else if(TD_DAFR < 0.0d0)then
       TD_DAFR = TD_DAFR + 86400.0d0
       TD_MJDN = TD_MJDN - 1
    end if
    
    return
  end subroutine renew_clock


  subroutine set_time1
    !---- args
    !---- vars
    !---- body
    UT_MJD      = UT_MJDN + UT_DAFR/86400.0d0
    UT_JD       = UT_MJD + 2400000.5d0
    TD_JD       = UT_JD + delta_T/86400.0d0
    UT_JC2000   = (UT_MJD - 51544.5d0)   / 36525.0d0
    TD_JC2000   = (TD_JD  - 2451545.0d0) / 36525.0d0

    return
  end subroutine set_time1


  subroutine set_time2(frag)
    !---- args
    integer, intent(in) :: frag
    !---- vars
    !---- body
    UT_MJD      = UT_MJDN + UT_DAFR/86400.0d0
    UT_JD       = UT_MJD + 2400000.5d0
    TD_JD       = UT_JD + delta_T/86400.0d0
    if(frag == 1)then
       UT_YEFR     = UT_YEFR + FOY
       UT_aphelion = UT_aphelion + 1       
    else if(frag == -1)then
       UT_YEFR = UT_YEFR - FOY
       UT_aphelion = UT_aphelion - 1              
    end if
    UT_JC2000   = (UT_MJD - 51544.5d0)   / 36525.0d0
    TD_JC2000   = (TD_JD  - 2451545.0d0) / 36525.0d0
    D_aphelion = pi2 * UT_aphelion / 365.0d0

    return
  end subroutine set_time2

  
  subroutine set_initial_time
    !---- args
    !---- vars
    integer :: Leap
    !---- body
    FOY = 1.0d0/365.25d0
    call calc_Julian_Day &
         (UT_year, UT_month, UT_day, UT_hour, UT_min, UT_sec, UT_JD)
    UT_MJD  = UT_JD - 2400000.5d0
    UT_MJDN = int(UT_MJD)
    UT_DAFR = UT_hour*3600.0d0 + UT_min*60.0d0 + UT_sec
    
    call YMD2Day_of_the_Year(UT_year,UT_month,UT_day,UT_DoY,Leap)
    UT_aphelion = UT_DoY - 155.0d0
    D_aphelion = pi2 * UT_aphelion / 365.0d0    
    
    if(Leap == 1)then
       UT_YEFR = UT_year + (UT_DoY + UT_DAFR/86400.0d0)/366.0d0
    else if(Leap == 2)then
       UT_YEFR = UT_year + (UT_DoY + UT_DAFR/86400.0d0)/365.0d0
    else
       print*, 'Error occured in "set_initial_time" of "time.f90"'
       print*, 'Error : Input Calender Date is invalid'
       stop
    end if

    call MJD2JC(UT_MJD, UT_JC2000)

    call set_deltaT(UT_year,UT_month,delta_T)

    TD_DAFR = UT_DAFR + delta_T
    TD_MJDN = UT_MJDN
    if(TD_DAFR >= 86400.0d0)then
       TD_DAFR = TD_DAFR - 86400.0d0
       TD_MJDN = TD_MJDN + 1
    end if
    TD_JD      = UT_JD + delta_T/86400.0d0    
    
    return
  end subroutine Set_Initial_Time


  subroutine calc_Julian_Day(year, month, day, hour, min, sec, JD)
    !****************************************************************************
    ! This subroutine calculates Julian Centuries from J2000.0
    ! that reference a particular time scale.
    !
    ! Input     : year, month , day, hour, min, sec
    ! Output    : Julian Day
    !
    ! Reference : David A Vallado, Fundamentals of Astrodynamics and Applications
    !                          - Third Edition - p189
    ! Author    : Keisuke Akari
    ! Date      : November 1, 2017
    !****************************************************************************    
    !--- args
    integer, intent(in)  :: year
    integer, intent(in)  :: month
    integer, intent(in)  :: day
    integer, intent(in)  :: hour
    integer, intent(in)  :: min
    real(8), intent(in)  :: sec
    real(8), intent(out) :: JD
    !--- vars
    !--- body
    JD = + ((sec/60.0d0 + min)/60.0d0 + hour)/24.0d0 + day &
         + 367*year - int((7*(year + int((month + 9)/12)))/4) &
         + int(275*month/9) + 1721013.5d0
    
    return
  end subroutine calc_Julian_Day
  

  subroutine JD2JC(JD,JC)
    !*************************************************************************
    ! This subroutine calculates Julian Centuries from J2000.0
    ! that reference a particular time scale.
    !
    ! Input  : Julian Day (JD)
    ! Output : Julian Century from J2000.0 (JC)
    !
    ! Ref    : David A Vallado, Fundamentals of Astrodynamics and Applications
    !                          - Third Edition - p190,191
    ! Author : Keisuke Akari
    ! Date   : November 1, 2016
    !*************************************************************************
    !---- args
    real(8), intent(in)  :: JD
    real(8), intent(out) :: JC
    !---- vars
    !---- body
    JC = (JD - 2451545.0d0) / 36525.0d0

    return
  end subroutine JD2JC


  subroutine MJD2JC(MJD,JC)
    !*************************************************************************
    ! This subroutine calculates Julian Centuries from J2000.0
    ! that reference a particular time scale.
    !
    ! Input  : Modified Julian Day (MJD)
    ! Output : Julian Century from J2000.0 (JC)
    !
    ! Ref    : David A Vallado, Fundamentals of Astrodynamics and Applications
    !                          - Third Edition - p190,191
    ! Author : Keisuke Akari
    ! Date   : November 1, 2016
    !*************************************************************************
    !---- args
    real(8), intent(in)  :: MJD
    real(8), intent(out) :: JC
    !---- vars
    !---- body
    JC = (MJD - 51544.5d0) / 36525.0d0    

    return
  end subroutine MJD2JC
  

  subroutine set_deltaT(year,month,delta_T)
    !*****************************************************************
    ! This subroutine calculate the value of the difference between
    ! Universal Time(UT) and Dynamical Time(TD).
    ! Although the exact value can be deduced only from observations,
    ! polynomial expressions for delta_T is employed in this program.
    ! In addition, the value of delta_T is assumed to be as a constant
    ! during this simulation of orbit propagation, though the value of
    ! delta_T depends on time.
    !
    ! Input  : Universal Time (year, month)
    ! Output : delta_T = TD - UT
    ! Ref    : NASA Eclipse Web Site
    !          http://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    ! Author : Keisuke Akari
    ! Date   : October 30, 2016
    !*****************************************************************
    !---- args
    integer, intent(in)  :: year
    integer, intent(in)  :: month
    real(8), intent(out) :: delta_T
    !---- vars
    real(8) :: y, T
    !---- body
    y = year + (month - 0.5d0)/12.0d0
    T = y - 2000.0d0
    
    if((year >= 1986) .and. (year < 2005))then
       delta_T = 63.86d0 + 3.345d-1*T - 6.0374d-2*T**2 &
            + 1.7275d-3*T**3 + 6.51814d-4*T**4 &
            + 2.373599d-5*T**5
    else if((year >= 2005) .and. (year < 2050))then
       delta_T = 62.92d0 + 0.32217d0*T + 5.589d-3*T**2
    else
       print*, 'Error occured in "set_deltaT" of "time.f90"'
       print*, 'Error : "Cannot calculate delta_T (= TD - UT)"'
       stop
    end if

    return
  end subroutine set_deltaT
  

  subroutine calc_GMST_0h(julian_century_2000, GMST_0h)
    !****************************************************************
    ! This subroutine calculates Greenwich Mean Sidereal Time at
    ! 0h UT of a given date.
    ! To achive the purpose, the modified Julian Century elapsed 
    ! from the epoch J2000.0 is required.
    ! In addition, UT1 time scale should be employed as an input.
    !
    ! Input     : Julian Century refered to J2000.0
    ! Output    : Greenwich Mean Sidereal Time at 00:00(UT) (GMST_0h)
    ! Reference : Fundamentals of Astrodynamics and Applications
    !                    -Third Edition-, p191~196
    ! Author    : Keisuke Akari
    ! Date      : 29 October, 2016 (modified on 1 Nov. 2017)
    !****************************************************************
    !---- args
    real(8), intent(in)  :: julian_century_2000
    real(8), intent(out) :: GMST_0h
    !---- vars
    !---- body
    GMST_0h = (67310.54841d0 &
         + (876600.0d0 + 8640184.812866d0)*julian_century_2000 &
         + 0.093104d0*julian_century_2000**2 &
         - 6.2d-6*julian_century_2000**3)*sec2rad
    if(GMST_0h < 0.0d0)then
       do while(GMST_0h < 0.0d0)
          GMST_0h = GMST_0h + pi2
       end do
    else if(GMST_0h >= pi2)then
       do while(GMST_0h >= pi2)
          GMST_0h = GMST_0h - pi2
       end do
    end if

    return
  end subroutine calc_GMST_0h


  subroutine calc_GST(GMST_0h, DAFR, GMST)
    !**********************************************************
    ! This subroutine calculates Greenwich sidereal time
    ! at a giben time. 
    !
    ! Input     : GMST_0h
    !                Greenwich Mean Sidereal Time at 0h UT
    !             DAFR
    !                A fraction of the day (UT)
    ! Output    : GMST
    !                Greenwich Mean Sidereal Time at the moment
    ! Reference : 長沢　工, 天体の位置計算 -増補版-, p23~25
    ! Author    : Keisuke Akari
    ! Date      : November 23, 2016 (modified on 1 Nov. 2017)
    !**********************************************************
    !---- args
    real(8), intent(in)  :: GMST_0h
    real(8), intent(in)  :: DAFR
    real(8), intent(out) :: GMST
    !---- vars
    !---- body
    GMST = GMST_0h + pi2*DAFR/86400.0d0
    if(GMST >= pi2)then
       GMST = GMST - pi2
    end if

    return
  end subroutine calc_GST  


  subroutine JD2CalenderDate(JD,year,month,day)
    !*****************************************************
    ! This subroutine calculate the Calender Date from JD.
    !
    ! Input  : Julian Date
    ! Output : year, month, day
    ! Ref    : Jean Meeus, Astronomical Algorithms
    !               - Second Edition -, p63
    ! Author : Keisuke Akari
    ! Date   : November 2, 2016
    !*****************************************************
    !---- args
    real(8), intent(in)  :: JD
    integer, intent(out) :: year, month, day
    !---- vars
    integer :: Z, A, alpha, B, C, D, E
    real(8) :: Y, F
    !---- body
    Y = JD + 0.5d0
    Z = int(Y)
    F = Y - Z
    
    if(Z < 2299161)then
       A = Z
    else if(Z >= 2299161)then
       alpha = int((Z - 1867216.25d0)/36524.25d0)
       A = Z + 1 + alpha - int(alpha/4)
    end if

    B = A + 1524
    C = int((B - 122.1d0)/365.25d0)
    D = int(365.25d0*C)
    E = int((B - D)/30.6001d0)
    
    day = B - D - int(30.6001d0*E) + F
    
    if(E < 14)then
       month = E - 1
    else if(E == 14 .or. E == 15)then
       month = E - 13
    end if

    if(month > 2)then
       year = C - 4716
    else if(month == 1 .or. month == 2)then
       year = C - 4715
    end if

    return
  end subroutine JD2CalenderDate


  subroutine YMD2Day_of_the_Year(year, month, day, DoY, Leap)
    !************************************************************
    ! This subroutine calculates the number of days elapsed from
    ! the beginning of the year.
    !
    ! Input  : year, month, day
    ! Output : Day of the Year(DoY)
    !          Leap Year Frag
    !          (for a leap year 'Leap' = 1, otherwise 'Leap' = 2)
    ! Ref    : Jean Meeus, Astronomical Algorithms
    !               - Second Edition - , p65
    ! Author : Keisuke Akari
    ! Date   : November 2, 2016
    !************************************************************
    !---- args
    integer, intent(in)  :: year, month, day
    integer, intent(out) :: DoY, Leap
    !---- vars
    !---- body
    if(year/4 == 0)then
       Leap = 1   ! for a leap year
    else
       Leap = 2   ! for a common year
    end if

    DoY = int(275*month/9) - Leap*int((month + 9)/12) + day - 30

    return
  end subroutine YMD2Day_of_the_Year
  
    
  subroutine Day_of_the_Year2YMD(N, year, month, day)
    !******************************************************
    ! This subroutine converts the day number N, which is
    ! the day elapsed from the beginning of the year into
    ! the corresponding date.
    !
    ! Input  : day of the year(N), year
    ! Output : month, day
    ! Ref    : Jean Meeus, Astronomical Algorithm
    !              - Second Edition -, p66
    ! Author : Keisuke Akari
    ! Date   : November 2, 2016
    !******************************************************
    !---- args
    integer, intent(in)  :: N, year
    integer, intent(out) :: month, day
    !---- vars
    integer :: K
    !---- body
    if(year/4 == 0)then
       K = 1   ! in the case of a leap year
    else
       K = 2   ! in the case of a common year
    end if

    if(N < 32)then
       month = 1
    else
       month = int(9.0d0*(K + N)/275.0d0 + 0.98d0)
    end if

    day = N - int(275*month/9) + K*int((month + 9)/12) + 30

    return
  end subroutine Day_of_the_Year2YMD  
end module time_mod
