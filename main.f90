program main
  use global_mod, only : set_physical_constants, &
       N_file,  N_debris,  N_ABM,  N_geo,  N_data_moon,  N_data_sun, &
       deltaXYZ,  XdeltaGEI,  YdeltaGEI,  ZdeltaGEI, &
       r0,  mu_m,  mu_e,  mu_s, &
       keisuu_moon,  keisuu_sun, &
       TD_JC2000,  TD_JD,  moon_interval,  sun_interval, &
       file_moon,  file_sun, &
       initial_index_moon,  initial_index_sun, &
       RdeGEI_inv
  use input_mod, only : &
       input_physical_constants, &
       input_simulation_condition, &
       input_data_of_debris
  use ABM_mod, only : &
       set_ABM_coefficients, &
       start_orbital_propagation
  use time_mod, only : &
       set_time1, &
       set_initial_time
  use geopotential_mod, only : &
       set_order_of_geopotential_polynomials, &
       input_geopotential_coefficients, &
       set_geopotential_coefficients
  use sun_mod, only : &
       set_sun_attraction_constants, &
       load_sun_ephemeris_table, &
       find_sun_ephemeris_JD
  use moon_mod, only : &
       set_moon_attraction_constants, &
       load_moon_ephemeris_table, &
       find_moon_ephemeris_JD
  use em_force_mod, only : &
       set_convection_electric_field_coefficients
  use weimer_mod, only : SetModel
  use Charge, only : Input_Charge_Database
  implicit none
  !---- vars
  integer :: i
  real(4) :: tt1,tt2
  !---- body
  call set_timer
  call input_physical_constants
  call input_simulation_condition
  call input_data_of_debris(N_debris)
  
#ifdef change_in_charge
  call Input_Charge_Database
#endif
  
  call Set_Initial_Time
  call set_time1

  call set_physical_constants  
  call set_ABM_coefficients(N_ABM)  
  call set_order_of_geopotential_polynomials(N_geo)
  allocate(RdeGEI_inv(1:N_geo))
  call input_geopotential_coefficients(N_geo)
  call set_geopotential_coefficients(N_geo)

  !******************************************************************  
  ! Moon
  call set_moon_attraction_constants(mu_m, mu_e, keisuu_moon)
  call load_moon_ephemeris_table &
       (TD_JC2000, N_file, file_moon, N_data_moon, moon_interval)
  call find_moon_ephemeris_JD(TD_JD, N_data_moon, initial_index_moon)
  !*******************************************************************

  !*********************************************************************
  ! Sun
  call set_sun_attraction_constants(mu_s, mu_e, keisuu_sun)
  call load_sun_ephemeris_table &
       (TD_JC2000, N_file, file_sun, N_data_sun, sun_interval)
  call find_sun_ephemeris_JD(TD_JD, N_data_sun, initial_index_sun)
  !*********************************************************************
  
  call set_convection_electric_field_coefficients &
       (deltaXYZ/r0, XdeltaGEI, YdeltaGEI, ZdeltaGEI)
  
  !******* Ordinary ********************************
  call SetModel(90.0,3.0,0.0,350.0,10.0,0.0,.false.)
  !*************************************************
  
  !******* Substorm ******************************
  !call SetModel(180.,15.,0.0,600.,10.,0.,.false.)
  !***********************************************

  do i = 1, N_debris
     call start_orbital_propagation(i)
  end do

  call finish_timer
contains
  subroutine set_timer
    !---- args
    !---- vars
    !---- body
    call cpu_time(tt1)

    return
  end subroutine set_timer


  subroutine finish_timer
    !---- args
    !---- vars
    !---- body    
    call cpu_time(tt2)
    write(*,*) 'CPU time =', (tt2-tt1)
    
    return
  end subroutine finish_timer
end program main
