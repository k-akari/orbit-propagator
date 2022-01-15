module Charge
  implicit none
  real(8), save, allocatable :: table_altitude(:)
  real(8), save, allocatable :: table_declination(:)
  real(8), save, allocatable :: table_right_ascension(:)
  real(8), save, allocatable :: spline_coefficients(:,:,:)
  integer, parameter :: n_hour = 288
  integer, parameter :: n_alt = 5
  integer, save :: n_dec
  integer, save :: n_RA
  real(8), parameter :: alt_interval = 50.0d0 ![km]
  real(8), parameter :: dec_interval = 2.0d0  ![deg.]
  real(8), parameter :: RA_interval  = 5.0d0  ![deg.]
  integer, parameter :: fi_charge = 999
  integer, save :: n_RA12
  integer, save :: n_RA14
  integer, save :: n_RA34
  integer, save :: n_dec12
contains
  subroutine Input_Charge_Database
    !--- args
    !--- vars
    integer :: i
    integer :: j
    integer :: k
    integer :: n
    integer :: m
    integer :: l
    !--- body
    n_dec = nint(180.0d0/dec_interval) + 1
    n_RA  = nint(360.0d0/RA_interval)
    
    n_dec12 = (n_dec-1)/2+1
    n_RA12  = n_RA/2+1
    n_RA14  = n_RA/4+1
    n_RA34  = 3*n_RA/4+1

    allocate(table_altitude(1:n_alt))
    allocate(table_declination(1:n_dec))
    allocate(table_right_ascension(1:n_RA))
    allocate(spline_coefficients(1:n_alt*n_dec*n_RA, 0:4, 1:n_hour))
    
    open(fi_charge, file='input/charge_database.dat', status='old')
    n = 0
    do i = 1, n_alt, 1
       do j = 1, n_dec, 1
          do k = 1, n_RA, 1
             n = n + 1
             read(fi_charge, *) table_altitude(i), table_declination(j), &
                  table_right_ascension(k), ((spline_coefficients(n,m,l), m=0,4), l=1,n_hour-1)
          end do
       end do
    end do

    print*, "*******************************************"
    print*, "Reading 'charge_database.dat' is completed."
    print*, "*******************************************"

    return
  end subroutine Input_Charge_Database

  
  subroutine Renew_Charge &
       (altitude, declination, right_ascension, time, charge)
    !--- args
    real(8), intent(in)  :: altitude        ![km]
    real(8), intent(in)  :: declination     ![deg.]
    real(8), intent(in)  :: right_ascension ![deg.]
    real(8), intent(in)  :: time   !Elapsed Time from Simulation Epoch [sec.]
    real(8), intent(out) :: charge ![C=A・s]
    !--- vars
    integer :: id_position
    !--- body
    call Identify_Position &
         (altitude, declination, right_ascension, id_position)

    call Interpolation_Time &
         (time, id_position, charge)

    return
  end subroutine Renew_Charge

  
  subroutine Identify_Position &
       (altitude, declination, right_ascension, id_position)
    !--- args
    real(8), intent(in)  :: altitude
    real(8), intent(in)  :: declination
    real(8), intent(in)  :: right_ascension
    integer, intent(out) :: id_position
    !--- vars
    integer :: id_altitude
    integer :: id_declination
    integer :: id_RA
    !--- body
    call Identify_Altitude(altitude, id_altitude)
    call Identify_Declination(declination, id_declination)
    call Identify_Right_Ascension(right_ascension, id_RA)

    id_position = id_RA + (id_declination-1)*n_RA + (id_altitude-1)*n_dec*n_RA

    return
  end subroutine Identify_Position

  
  subroutine Identify_Altitude(altitude, id_altitude)
    !--- args
    real(8), intent(in)  :: altitude ![km]
    integer, intent(out) :: id_altitude
    !--- vars
    integer :: i
    !--- body
    if(table_altitude(n_alt) - 0.5d0*alt_interval <= altitude)then
       id_altitude = n_alt
       return
    end if
    
    do i = n_alt-1, 2, -1
       if((table_altitude(i) - 0.5d0*alt_interval <= altitude) .and. &
            (altitude < table_altitude(i) + 0.5d0*alt_interval))then
          id_altitude = i

          return
       end if
    end do

    if(altitude < table_altitude(1) + 0.5d0*alt_interval)then
       id_altitude = 1
       return
    end if
    
    print*, "Error : 'Identify_Altitude' in 'charge.f90'"
    print*, "Error : input altitude is invalid"
    print*, "altitude =", altitude
    stop
    
  end subroutine Identify_Altitude

  
  subroutine Identify_Declination(declination, id_declination)
    !--- args
    real(8), intent(in)  :: declination ![degree]
    integer, intent(out) :: id_declination
    !--- vars
    integer :: i
    !--- body
    if((0.0d0 <= declination) .and. (declination <= 90.0d0))then
       do i = n_dec12, n_dec, 1
          if((table_declination(i) - 0.5d0*dec_interval <= declination) .and. &
               (declination < table_declination(i) + 0.5d0*dec_interval))then
             id_declination = i

             return
          end if
       end do
    else if((-90.0d0 <= declination) .and. (declination <= 0.0d0))then
       do i = n_dec12, 1, -1
          if((table_declination(i) - 0.5d0*dec_interval <= declination) .and. &
               (declination < table_declination(i) + 0.5d0*dec_interval))then
             id_declination = i

             return
          end if
       end do
    else
       print*, "Error : 'Identify_Declination' in 'charge.f90'"
       print*, "Error : input declination is invalid"
       print*, "Declination =", declination
       
       stop
    end if
    
    return
  end subroutine Identify_Declination

  
  subroutine Identify_Right_Ascension(right_ascension, id_RA)
    !--- args
    real(8), intent(in)  :: right_ascension ![degree]
    integer, intent(out) :: id_RA
    !--- vars
    integer :: i
    !--- body
    if((table_right_ascension(1) <= right_ascension) .and. &
         (right_ascension <= table_right_ascension(n_RA12)))then
       if((table_right_ascension(1) <= right_ascension) .and. &
            (right_ascension <= table_right_ascension(n_RA14)))then
          do i = 1, n_RA14, 1
             if((table_right_ascension(i) - 0.5d0*RA_interval <= right_ascension) .and. &
                  (right_ascension < table_right_ascension(i) + 0.5d0*RA_interval))then
                id_RA = i
                
                return
             end if
          end do
       else if((table_right_ascension(n_RA14) < right_ascension) .and. &
            (right_ascension <= table_right_ascension(n_RA12)))then
          do i = n_RA14, n_RA12, 1
             if((table_right_ascension(i) - 0.5d0*RA_interval <= right_ascension) .and. &
                  (right_ascension < table_right_ascension(i) + 0.5d0*RA_interval))then
                id_RA = i
                
                return
             end if
          end do          
       end if
    else if((table_right_ascension(n_RA12) < right_ascension) .and. &
         (right_ascension <= 360.0d0))then
       if((table_right_ascension(n_RA12) <= right_ascension) .and. &
            (right_ascension <= table_right_ascension(n_RA34)))then
          do i = n_RA12, n_RA34, 1
             if((table_right_ascension(i) - 0.5d0*RA_interval <= right_ascension) .and. &
                  (right_ascension < table_right_ascension(i) + 0.5d0*RA_interval))then
                id_RA = i
                
                return
             end if
          end do          
       else if((table_right_ascension(n_RA34) < right_ascension) .and. &
            (right_ascension <= 360.0d0))then
          do i = n_RA34, n_RA, 1
             if((table_right_ascension(i) - 0.5d0*RA_interval <= right_ascension) .and. &
                  (right_ascension < table_right_ascension(i) + 0.5d0*RA_interval))then
                id_RA = i
                
                return
             end if
          end do
          if((table_right_ascension(n_RA) + 0.5d0*RA_interval <= right_ascension) .and. &
               (right_ascension < table_right_ascension(n_RA) + RA_interval))then
             id_RA = 1

             return
          end if
       end if       
    else
        
       print*, "Error : 'Identify_Right_ascension' in 'charge.f90'"
       print*, "Error : input right ascension is invalid"
       print*, "Right Ascension =", right_ascension
       
       stop
    end if
    
    return
  end subroutine Identify_Right_Ascension


  subroutine Interpolation_Time(time, id_position, charge)
    !--- args
    real(8), intent(in)  :: time ![sec.]
    integer, intent(in)  :: id_position
    real(8), intent(out) :: charge ![C]
    !--- vars
    real(8) :: hour
    !--- body
    hour = time/3600.0d0
    charge = Poly3(hour, spline_coefficients(id_position, :, int(hour)+1))

    return
  end subroutine Interpolation_Time

  
  function Poly3(x,c)
    ! returns value of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
    !--- args
    real(8)             :: poly3
    real(8), intent(in) :: x
    ! point at which to evaluate polynomial
    real(8), intent(in) :: c(0:4)
    ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
    !--- vars
    real(8) dx
    !--- body
    dx=x-c(0)
    poly3=c(1)+c(2)*dx+c(3)*dx**2+c(4)*dx**3

    return
  end function Poly3
end module Charge

  
