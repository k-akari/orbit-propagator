module input_mod
  use global_mod, only : &
       G, mu_e, mu_S, mu_m, r0, v0, Rsun, AU, &
       eps0, c, omega_e, eps_r, Cd, &
       dt_imverse, simulation_year, simulation_day, i_step, &
       N_ABM, N_debris, N_geo, N_air, Ud, &
       UT_year, UT_month, UT_day, UT_hour, UT_min, UT_sec, &
       limit_h, deltaXYZ, VGSEX, VGSEY, VGSEZ, &
       N_file, file_sun, file_moon, table_debris, &
       fi1, fi_deb, &
       Rd, Md, Ad, XdeGEI_initial, VdeGEI_initial
  implicit none
  character(60) :: debris_file
contains
  subroutine input_physical_constants
    !--- args
    !--- vars
    character(128) :: text
    !--- body
    open(fi1, file = 'input/constants.dat', status = 'old')
    
    read(fi1,*) text
    read(fi1,*) text
    read(fi1,*) text    

    read(fi1,*) text, G
    read(fi1,*) text, mu_e
    read(fi1,*) text, mu_s
    read(fi1,*) text, mu_m
    read(fi1,*) text, r0
    read(fi1,*) text, Rsun
    read(fi1,*) text, AU
    read(fi1,*) text, eps0
    read(fi1,*) text, c
    read(fi1,*) text, omega_e(1), omega_e(2), omega_e(3)
    read(fi1,*) text
    read(fi1,*) text, eps_r
    read(fi1,*) text
    read(fi1,*) text, Cd

    close(fi1)
    
    return
  end subroutine input_physical_constants
  

  subroutine input_simulation_condition
    !--- args
    !--- vars
    integer        :: i
    character(128) :: text
    character(20)  :: name1
    character(20)  :: name2
    character(20)  :: name3
    character(20)  :: name4
    !--- body
    open(fi1, file = 'input/setting.dat',   status = 'old')
    
    read(fi1,*) text
    read(fi1,*) text

    read(fi1,*) text, debris_file    
    read(fi1,*) text, dt_imverse
    read(fi1,*) text, simulation_year
    read(fi1,*) text, simulation_day
    read(fi1,*) text, i_step

    read(fi1,*) text, N_ABM
    read(fi1,*) text, N_debris
    read(fi1,*) text, N_geo
    read(fi1,*) text, Ud
    read(fi1,*) text, UT_year, UT_month, UT_day, UT_hour, UT_min, UT_sec
    read(fi1,*) text, limit_h
    read(fi1,*) text, deltaXYZ
    read(fi1,*) text, VGSEX, VGSEY, VGSEZ
    read(fi1,*) text, N_air

    read(fi1,*) text, N_file
    
    allocate(file_sun(1:N_file))
    allocate(file_moon(1:N_file))
    do i = 1, N_file
       read(fi1,*) name1, name2, name3, name4
       file_sun(i)  = trim(adjustl(name1))//'/'//trim(adjustl(name2))
       file_moon(i) = trim(adjustl(name3))//'/'//trim(adjustl(name4))
    end do    
    close(fi1)

    print*, "*****************************************************"    
    print*, "This simulation uses Geocentric Equatorial Coordinate"
    print*, "referred to the date(UT) as shown below."
    print*, "-----------------------------------------------------"    
    print*, "       Year       ","    Month      ","Day       "
    print*, ""
    print*, UT_year, UT_month, UT_day
    print*, ""
    print*, "-----------------------------------------------------"    
    print*, "       Hour       ","    Min        ","Sec       "
    print*, ""
    print*, UT_hour, UT_min, UT_sec
    print*, ""
    print*, "*****************************************************"
    
    return
  end subroutine input_simulation_condition

  
  subroutine input_data_of_debris(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    integer :: i
    !---- body
    open(fi_deb, file = debris_file, status = 'old')
    
    allocate(table_debris(9,N))

    do i = 1, 6
       read(fi_deb,*)
    end do
    !*******************************************************************
    ! For an input file containing objects data for orbital calculation,
    ! data should be put in order as follows.
    !*******************************************************************
    !   Radius [m], projected cross-sectional area[m^2], mass [kg]
    !   X [m],    Y [m],    Z [m]
    !   Vx [m/s], Vy [m/s], Vz [m/s]
    !*****************************************************************
    
    do i = 1, N
       read(fi_deb,*) table_debris(1,i), table_debris(2,i), table_debris(3,i)
       read(fi_deb,*) table_debris(4,i), table_debris(5,i), table_debris(6,i)
       read(fi_deb,*) table_debris(7,i), table_debris(8,i), table_debris(9,i)
       read(fi_deb,*)
    end do
    close(fi_deb)

    return
  end subroutine input_data_of_debris
  
  subroutine set_debris(N)
    !---- args
    integer, intent(in) :: N
    !---- vars
    !---- body
    Rd                = table_debris(1,N)
    Ad                = table_debris(2,N)
    Md                = table_debris(3,N)
    XdeGEI_initial(1) = table_debris(4,N)/r0
    XdeGEI_initial(2) = table_debris(5,N)/r0
    XdeGEI_initial(3) = table_debris(6,N)/r0
    VdeGEI_initial(1) = table_debris(7,N)/v0
    VdeGEI_initial(2) = table_debris(8,N)/v0
    VdeGEI_initial(3) = table_debris(9,N)/v0

    return
  end subroutine set_debris
end module input_mod
