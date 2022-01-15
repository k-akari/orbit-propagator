module moon_mod
  use global_mod, only : &
       r0, fi_moon, TD_JD_moon, posi_table_moon
  use conversion_coordinates_mod, only : &
       set_precession_matrix, &
       precession
  use subprog_mod, only : &
       interpolation, &
       nolm
  implicit none
contains
  subroutine set_moon_attraction_constants(mu_m, mu_e, Cout)
    !---- args
    real(8), intent(in)  :: mu_m ! Coefficient for calculation of Moon's attraction
    real(8), intent(in)  :: mu_e ! Gravity constant of the Moon
    real(8), intent(out) :: Cout ! Gravity constant of the Earth
    !---- vars
    !---- body
    Cout = mu_m / mu_e
    
    return
  end subroutine set_moon_attraction_constants


  subroutine load_moon_ephemeris_table &
       (Tin, num_file, filename, num_data, interval)
    !**********************************
    ! Store all ephemeris data in array
    !**********************************
    !--- args
    real(8), intent(in)          :: Tin
    integer, intent(in)          :: num_file
    character(len=*), intent(in) :: filename(1:num_file)
    integer, intent(out)         :: num_data
    real(8), intent(out)         :: interval
    !***************************************************
    ! Tin      :
    ! num      : The Number of Input Ephemeris Data File
    ! filename : Ephemeris File Name
    !***************************************************
    !--- vars
    integer :: idx
    integer :: iost
    integer :: i
    real(8) :: sze,cze,sz,cz,st,ct
    real(8) :: XseKM(3)
    real(8) :: JD_moon
    real(8) :: moon_interval1
    real(8) :: moon_interval2    
    !--- body
    call set_precession_matrix(Tin,sze,cze,sz,cz,st,ct)
    call count_moon_ephemeris_file(num_file, filename(:), num_data)
    allocate(TD_JD_moon(num_data))
    allocate(posi_table_moon(1:3,1:num_data))

    idx = 1
    ! read ephemeris data recorded in each file
    do i = 1, num_file
       open(fi_moon, file = filename(i), status = 'old')
       rewind(fi_moon)
       do while(.true.)
          read(fi_moon, *, iostat = iost) JD_moon, XseKM(1), XseKM(2), XseKM(3)
          if(iost /= 0)then
             exit ! finish counting if all data are read
          end if
          ! add to entry arrangements
          TD_JD_moon(idx) = JD_moon
          call precession(XseKM(1:3)*1.0d3/r0,sze,cze,sz,cz,st,ct,&
               posi_table_moon(1:3, idx))
          idx = idx + 1
       end do
       close(fi_moon)
    end do

    print*, "*****************************************************"
    print*, "Conversion of the Moon's ephemeris is completed,"
    print*, "considering the Precession of the Earth."
    print*, "-----------------------------------------------------"
    print*, ""
    print*, "Julian Century refered to J2000.0(DT)"
    print*, Tin
    print*, "" 
    print*, "*****************************************************"        

    ! check if time interval is a constant
    do i = 1, num_data, 2
       moon_interval1 = TD_JD_moon(i+1) - TD_JD_moon(i)
       moon_interval2 = TD_JD_moon(i+2) - TD_JD_moon(i+1)
       if(moon_interval1 - moon_interval2 > 1.0d-5)then
          print*, "Error : Time Interval of the Moon's Ephemeris is not a Constant"
          print*, "Data No.", i
          print*, "Stop at moon.f90"
          stop
       end if
    end do
    interval = TD_JD_moon(2) - TD_JD_moon(1)

    return
  end subroutine load_moon_ephemeris_table


  subroutine count_moon_ephemeris_file(num_file, file_name, count)
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
       open(fi_moon, file = file_name(i), status = 'old')
       rewind(fi_moon)
       do while(.true.)
          read(fi_moon, '(A255)', iostat = iost) buff
          if(iost /= 0)then
             exit
          end if
          num = num + 1
       end do
       close(fi_moon)
    end do
    count = num

    return
  end subroutine count_moon_ephemeris_file

  
  subroutine find_moon_ephemeris_JD(Tin,num_data,index)
  !************************************************************
  ! Find the closest index to the given JD from ephemeris table
  !************************************************************    
    !--- args
    real(8), intent(in)  :: Tin      ! Julian Day (TD)
    integer, intent(in)  :: num_data ! the number of ephemeris data
    integer, intent(out) :: index    ! initial_index_moon
    !--- vars
    integer :: id_min, id_max, id_mid
    real(8) :: jd_mid
    !--- body

    ! check the input range
    if(Tin < TD_JD_moon(1))then
       print*, 'Error : Ephemeris Table is INVALID for the Simulation Epoch'
       stop
    end if
    if(Tin > TD_JD_moon(num_data))then
       print*, 'Error : Ephemeris Table is INVALID for the Simulation Epoch'
       stop
    end if

    id_min = 1
    id_max = num_data

    ! find the closest index to given JD by half deviding method
    do while(id_min < id_max)
       id_mid = id_min + (id_max - id_min)/2
       jd_mid = TD_JD_moon(id_mid)
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
          print*, 'Stop at moon.f90'
          stop
       end if
    end do

    return
  end subroutine find_moon_ephemeris_JD


  subroutine calc_moon_position(Tin, interval, index, Xout)
    !*********************************************
    ! Choose the 3 tabular data points, and decide
    ! the position vector of the Moon
    !*********************************************
    !---- args
    real(8), intent(in)     :: Tin      ! Julian Day (TD)
    real(8), intent(in)     :: interval
    integer, intent(inout)  :: index    ! index of the Moon
    real(8), intent(out)    :: Xout(3)  ! Position Vector of the Moon
    !---- vars
    real(8) :: XmeGEI3(3), XmeGEI2(3), XmeGEI1(3)
    !---- body
    if(abs(Tin - TD_JD_moon(index)) &
         > abs(Tin - TD_JD_moon(index+1)))then
       index = index + 1
    end if
    
    XmeGEI3(:) = posi_table_moon(:,index+1)
    XmeGEI2(:) = posi_table_moon(:,index)
    XmeGEI1(:) = posi_table_moon(:,index-1)
    
    call interpolation(XmeGEI1(:),XmeGEI2(:),XmeGEI3(:), &
         Tin, TD_JD_moon(index), interval, Xout(:))

    return
  end subroutine calc_moon_position    
  
    
  subroutine distance_from_the_moon &
       (Xde, Xme, Xdm, Rdm, Rdm_inv, Rme, Rme_inv)
    !---- args
    real(8), intent(in)  :: Xde(3)
    real(8), intent(in)  :: Xme(3)
    real(8), intent(out) :: Xdm(3)
    real(8), intent(out) :: Rdm
    real(8), intent(out) :: Rdm_inv(3)
    real(8), intent(out) :: Rme
    real(8), intent(out) :: Rme_inv(3)
    !---- vars
    !---- body
    Xdm(:)     = Xde(:) - Xme(:)
    Rdm        = nolm(Xdm)
    Rdm_inv(1) = 1.0d0/Rdm
    Rdm_inv(2) = Rdm_inv(1)*Rdm_inv(1)
    Rdm_inv(3) = Rdm_inv(1)*Rdm_inv(2)
    
    Rme        = nolm(Xme)    
    Rme_inv(1) = 1.0d0/Rme
    Rme_inv(2) = Rme_inv(1)*Rme_inv(1)
    Rme_inv(3) = Rme_inv(1)*Rme_inv(2)

    return
  end subroutine distance_from_the_moon
  
  
  subroutine attraction_of_the_moon &
       (Cin, Xdm, Rdm_inv3, Xme, Rme_inv3, Aout)
    !---- args
    real(8), intent(in)  :: Cin
    real(8), intent(in)  :: Xdm(3)
    real(8), intent(in)  :: Rdm_inv3
    real(8), intent(in)  :: Xme(3)
    real(8), intent(in)  :: Rme_inv3
    real(8), intent(out) :: Aout(3)
    !---- vars
    !---- body
    Aout(:) = -Cin*(Xdm(:)*Rdm_inv3 + Xme(:)*Rme_inv3)

    return
  end subroutine attraction_of_the_moon
end module moon_mod
