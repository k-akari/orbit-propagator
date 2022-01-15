module charge_solve_mod
  implicit none
  real(8), parameter :: e    = 1.60217657d-19
  real(8), parameter :: kB   = 1.3806488d-23
  real(8), parameter :: M_e  = 9.10938356d-31
  real(8), parameter :: M_H  = 1.672621284d-27 !M_electron * 1836.152
  real(8), parameter :: M_He = 6.690485136d-27 ! 4  * M_proton
  real(8), parameter :: M_O  = 2.676194054d-26 ! 16 * M_proton
  real(8), parameter :: M_O2 = 5.352388109d-26 ! 32 * M_proton (16 * M_protonでは?)
  real(8), parameter :: pi = 3.14159265358979323846

  integer, save :: JF(50), iut, jmag
contains
  subroutine initialize_IRI
    !---- args
    !---- vars
    integer :: i
    !---- body
    do i = 1, 50, 1
       JF(i) = 1
    end do

    JF(4)  = 0
    JF(5)  = 0
    JF(6)  = 0
    JF(21) = 0
    JF(22) = 0 ! give the actual ion density
    JF(23) = 0
    JF(28) = 0
    JF(29) = 0
    JF(30) = 0
    JF(33) = 0 ! auroral boundary model
    JF(34) = 0 ! message off
    JF(35) = 0

    ! Check Options
    iut    = 1 !iut=1: universal time, iut=0: local time
    jmag   = 0 !geomagnetic
    JF(33) = 1 !auroral boundary
    JF(35) = 1 !foE storm
    JF(36) = 0 !hmF2 foF2 storm
    JF(37) = 0 !topside foF2 storm

    call read_ig_rz()
    call readapf107()
    
    return
  end subroutine initialize_IRI

  
  subroutine solve_current_balance &
       (surface_area, vd, height, latitude, longitude, &
       year, mmdd, dhour, inp_shadow_function, surface_potential)
    !---- args
    real(8), intent(in)    :: surface_area,vd,height,latitude,longitude,dhour
    integer, intent(in)    :: year, mmdd, inp_shadow_function
    real(8), intent(out)   :: surface_potential
    !---- vars
    integer :: i, shadow_function
    real(8) :: fx, fx_prime, pre_surface_potential
    real(8) :: Ie, Iph, IH_thermal, IO_thermal, IO2_thermal, IHe_thermal
    real(8) :: IH_drift, IO_drift, IO2_drift, IHe_drift
    real(8) :: Ie_prime, Iph_prime, IH_prime, IO_prime, IO2_prime, IHe_prime
    real(8) :: N_e, N_H, N_O, N_O2, N_He
    real(8) :: T_e, T_i, v_e, v_O, v_H, v_He, v_O2
    real(8) :: N(5), T(2), v(5)        
    !---- body
    !conversion
    if(inp_shadow_function == 2)then
       shadow_function = 0
    else if(inp_shadow_function == 0)then
       shadow_function = 1
    else
       print*, "Error: charge_solver.f90 solve_current_balance"
       stop
    end if
    
    !Plasma Parameter
    call Plasma_Parameter &
         (height,latitude,longitude,year,mmdd,dhour,N,T,v)

    N_e  = N(1)
    N_O  = N(2)
    N_H  = N(3)
    N_He = N(4)
    N_O2 = N(5)
    
    T_e  = T(1)
    T_i  = T(2)
    
    v_e  = v(1)
    v_O  = v(2)
    v_H  = v(3)
    v_He = v(4)
    v_O2 = v(5)

    
    !初期化
    surface_potential     = -50.0d0
    pre_surface_potential = 50.0d0
    fx       = 0.0d0
    fx_prime = 0.0d0

    !Newton-Rapson法による数値解法
    do i = 0, 1000, 1
       ! electron current
       Ie  = get_electron_current(N_e, T_e, v_e, surface_area, surface_potential)
       Ie_prime  = get_electron_current_prime(T_e, N_e, v_e, surface_area, surface_potential)

       ! photoelectron current
       Iph = get_photoelectron_current(surface_area, shadow_function)
       Iph_prime = 0.0d0

       ! H+ ion current
       if(N_H > 0.0d0)then
          IH_thermal  = get_thermal_ion_current(N_H, T_i, v_H, surface_area, surface_potential)
          IH_drift    = get_drift_ion_current(surface_area, vd, N_H)
          IH_prime    = get_ion_current_prime(T_i, N_H, v_H, surface_area, surface_potential)
       else
          IH_thermal = 0.0d0
          IH_drift   = 0.0d0
          IH_prime   = 0.0d0
       end if
       
       ! O+ ion current
       if(N_O > 0.0d0)then
          IO_thermal  = get_thermal_ion_current(N_O, T_i, v_O, surface_area, surface_potential)
          IO_drift    = get_drift_ion_current(surface_area, vd, N_O)
          IO_prime    = get_ion_current_prime(T_i, N_O, v_O, surface_area, surface_potential)
       else
          IO_thermal  = 0.0d0
          IO_drift    = 0.0d0
          IO_prime    = 0.0d0
       end if
       
       ! O2+ ion current
       if(N_O2 > 0.0d0)then
          IO2_thermal = get_thermal_ion_current(N_O2, T_i, v_O2, surface_area, surface_potential)
          IO2_drift   = get_drift_ion_current(surface_area, vd, N_O2)
          IO2_prime   = get_ion_current_prime(T_i, N_O2, v_O2, surface_area, surface_potential)
       else
          IO2_thermal = 0.0d0
          IO2_drift   = 0.0d0
          IO2_prime   = 0.0d0
       end if
       
       ! He+ ion current
       if(N_He > 0.0d0)then
          IHe_thermal = get_thermal_ion_current(N_He, T_i, v_He, surface_area, surface_potential)
          IHe_drift   = get_drift_ion_current(surface_area, vd, N_He)
          IHe_prime   = get_ion_current_prime(T_i, N_He, v_He, surface_area, surface_potential)
       else
          IHe_thermal = 0.0d0
          IHe_drift   = 0.0d0
          IHe_prime   = 0.0d0
       end if
       
       fx       = Ie - IH_thermal - IO_thermal - IO2_thermal - IHe_thermal &
            - IH_drift - IO_drift - IO2_drift - IHe_drift -Iph
       fx_prime = Ie_prime - IH_prime - IO_prime - IO2_prime - IHe_prime - Iph_prime
       
       pre_surface_potential = surface_potential
       surface_potential     = surface_potential - fx/fx_prime
       
       if(abs(pre_surface_potential - surface_potential) <= 1.0d-7)then
          exit
       end if
    end do

    return
  end subroutine solve_current_balance


  subroutine Plasma_Parameter &
       (height,latitude,longitude,year,mmdd,dhour,N,T,v)
    !---- args
    real(8), intent(in)  :: height, latitude, longitude, dhour
    integer, intent(in)  :: year, mmdd
    real(8), intent(out) :: N(1:5), T(1:2), v(1:5)
    !---- vars
    real(4) :: OUTF(1:20,1:1000), OARR(1:100)
    real(8) :: xhour
    integer :: i
    !---- body
    !プラズマパラメータの取得
    call Check_Inputs(height,latitude,longitude,year,mmdd,dhour)
    xhour = dhour+iut*25.0d0
    OARR(:) = -1.0
    call IRI_SUB(JF,jmag,real(latitude),real(longitude),year, &
         mmdd,real(xhour),real(height),real(height),1.0,OUTF,OARR)
    !Universal Timeを用いる場合は25.0を足す

    !Density
    N(1) = dble(OUTF(1,1)) !Electron Density
    N(2) = dble(OUTF(5,1)) !O+ Density
    N(3) = dble(OUTF(6,1)) !H+ Density
    N(4) = dble(OUTF(7,1)) !He+ Density
    N(5) = dble(OUTF(8,1)) !O2+ Density

    !Temperature
    T(1)   = Temperature2Electronvolts(dble(OUTF(4,1))) !Electron Temperature[eV]
    T(2)   = Temperature2Electronvolts(dble(OUTF(3,1))) !Ion Temperature[eV]
    
    !Molecular Velocity
    v(1) = get_molecular_velocity(T(1),M_e)  !Electron Velocity
    v(2) = get_molecular_velocity(T(2),M_O)  !O+ ion Velocity
    v(3) = get_molecular_velocity(T(2),M_H)  !H+ ion Velocity
    v(4) = get_molecular_velocity(T(2),M_He) !He+ ion Velocity
    v(5) = get_molecular_velocity(T(2),M_O2) !O2+ ion Velocity    

    return
  end subroutine Plasma_Parameter
  

  subroutine Check_Inputs(height,latitude,longitude,year,mmdd,dhour)
    !---- args
    real(8), intent(in) :: height, latitude, longitude, dhour
    integer, intent(in) :: year, mmdd
    !---- vars
    !---- body
    if(height < 80.0d0)then
       print*, 'Error: height must be higher than 80km.'
       print*, 'Error: charge_solver.f90, Check_Inputs'
       stop
    end if
    
    if(latitude < -90.0d0 .or. 90.0d0 < latitude)then
       print*, 'Error: latitude must be in the range from -90 deg. to 90 deg.'
       print*, 'Error: charge_solver.f90, Check_Inputs'
       stop
    end if
    
    if(longitude < -180.0d0 .or. 360.0d0 < longitude)then
       print*, 'Error: longitude must be in the range from 0 deg. to 360 deg.'
       print*, '       ( if -180 to 0, it will be converted to 0 to 360.'
       print*, 'Error: charge_solver.f90, Check_Inputs'
       stop
    end if

    if(2019 < year)then
       print*, 'Error: year must be lower than 2020.'
       print*, 'Error: charge_solver.f90, Check_Inputs'
       stop
    end if

    if(mmdd < -365.0d0 .or. 1231 < mmdd)then
       print*, 'Error: mmdd must be in the range from -365 to 1231.'
       print*, 'Error: charge_solver.f90, Check_Inputs'
       stop
    end if

    if(24.0d0 < dhour)then
       print*, 'Error: dhour must be lower than 24.0'
       print*, 'Error: charge_solver.f90, Check_Inputs'
       stop
    end if

    return
  end subroutine Check_Inputs


  function Temperature2Electronvolts(temp)
    !---- args
    real(8), intent(in) :: temp
    real(8)             :: Temperature2Electronvolts
    !---- vars
    !---- body
    Temperature2Electronvolts = temp * kB / e

    return
  end function Temperature2Electronvolts
  
    
  function get_electron_current(Ne, Te, ve, surface_area, surface_potential)
    !---- args
    real(8), intent(in) :: surface_area, Ne, Te, ve, surface_potential
    real(8)             :: get_electron_current    
    !---- vars
    real(8) :: Ie0, fact
    !---- body
    Ie0 = surface_area * Ne * e * ve / (2 * sqrt(pi))
    
    if(surface_potential < 0.0d0)then
       fact = exp(surface_potential/Te)
    else if(surface_potential >= 0.0d0)then
       fact = 1.0d0 + surface_potential/Te
    else
       print*, 'Error: charge_solver.f90, get_electron_current'
       stop
    end if

    get_electron_current = Ie0 * fact
    
    return
  end function get_electron_current


  function get_thermal_ion_current(Ni, Ti, vi, surface_area, surface_potential)
    !---- args
    real(8), intent(in) :: Ni, Ti, vi, surface_area, surface_potential
    real(8)             :: get_thermal_ion_current    
    !---- vars
    real(8) :: Ii0, fact, Ii_thermal, Ii_drift
    !---- body
    Ii0 = surface_area * Ni * e * vi / (2.0d0 * sqrt(pi))
    
    if(surface_potential > 0.0d0)then
       fact = exp(-surface_potential/Ti)
    else if(surface_potential <= 0.0d0)then
       fact = 1.0d0 - surface_potential/Ti
    else
       print*, 'Error: charge_solver.f90, get_ion_current'
       stop
    end if

    get_thermal_ion_current = Ii0 * fact  ! Thermal Ion Current
    
    return
  end function get_thermal_ion_current


  function get_drift_ion_current(surface_area,vd,Ni)
    !---- args
    real(8), intent(in) :: surface_area, vd, Ni
    !---- vars
    real(8)             :: get_drift_ion_current
    real(8)             :: projected_surface_area
    !---- body
    projected_surface_area = 0.25d0 * surface_area
    get_drift_ion_current = projected_surface_area * Ni * vd * e ! Drift Current

    return
  end function get_drift_ion_current
  

  function get_photoelectron_current(surface_area, shadow_function)
    !---- args
    integer, intent(in)  :: shadow_function
    real(8), intent(in)  :: surface_area
    real(8)              :: get_photoelectron_current
    !---- vars
    real(8), parameter :: Jph = 3.0d-5
    real(8)            :: projected_surface_area
    !---- body
    projected_surface_area = 0.25d0 * surface_area
    get_photoelectron_current = Jph * projected_surface_area * shadow_function

    return
  end function get_photoelectron_current


  function get_electron_current_prime(Te, Ne, ve, surface_area, surface_potential)
    !---- args
    real(8), intent(in) :: surface_area, Te, Ne, ve, surface_potential
    real(8)             :: get_electron_current_prime
    !---- vars
    real(8) :: Ie0
    !---- body
    Ie0 = surface_area * Ne * e * ve / (2.0d0 * sqrt(pi))

    if(surface_potential >= 0.0d0)then
       get_electron_current_prime = Ie0 / Te
    else if(surface_potential < 0.0d0)then
       get_electron_current_prime = &
            get_electron_current(Ne,Te,ve,surface_area,surface_potential) / Te
    else
       print*, 'Error: charge_solver.f90, get_electron_current_prime'
       stop
    end if

    return
  end function get_electron_current_prime


  function get_ion_current_prime(Ti, Ni, vi, surface_area, surface_potential)
    !---- args
    real(8), intent(in) :: surface_area, Ti, Ni, vi, surface_potential
    real(8)             :: get_ion_current_prime
    !---- vars
    real(8) :: Ii0
    !---- body
    Ii0 = surface_area * Ni * e * vi / (2.0d0 * sqrt(pi))

    if(surface_potential <= 0.0d0)then
       get_ion_current_prime = -Ii0 / Ti
    else if(surface_potential > 0.0d0)then
       get_ion_current_prime = &
            - get_thermal_ion_current(Ni,Ti,vi,surface_area,surface_potential) / Ti
    else
       print*, 'Error: charge_solver.f90, get_ion_current_prime'
       stop
    end if

    return
  end function get_ion_current_prime


  function get_molecular_velocity(T,M)
    !---- args
    real(8), intent(in) :: T, M
    real(8)             :: get_molecular_velocity
    !---- vars
    !---- body

    get_molecular_velocity = sqrt(2.0d0 * e * T / M)

    return
  end function get_molecular_velocity
end module charge_solve_mod
