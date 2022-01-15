module output_mod
  use global_mod, only : &
       Aearth, Ageo, Aair, Aph, Amoon, Asun, VcrossB, Ecr, Ecn, &
       g0, r0, t0, v0, rad2deg, &
       alt_gd, lon_gd, lat_gd, &
       lon_gc, lat_gc, &
       Ud, charge_d, Md, time, X, V, K, &
       fo1,  fo2,  fo3,  fo4,  fo5, &
       fo6,  fo7,  fo8,  fo9,  fo10, &
       fo11, fo12, fo13, fo14, fo15
  use subprog_mod, only : nolm
  use conversion_coordinates_mod, only : Car2NTW, Car2RSW, Car2Kep
  use time_mod, only : JD2CalenderDate, UT_MJDN,UT_DAFR
  implicit none
contains
  subroutine output
    !---- args
    !---- vars
    real(8) :: output_day
    real(8) :: R
    real(8) :: alt_gc
    integer :: year, month, day, mmdd
    real(8) :: dhour
    !---- body
    !******************************************************************
    !   Car2Kep
    !******************************************************************
    !   Input  : X(1:3) - Position vector of debris [RE]
    !            V(1:3) - Velocity vector of debris [normalized]                       
    !   Output : Keplerian Elements
    !            K(1) - Orbital Inclination [rad]　
    !            K(2) - Right Ascension of the Ascending Node [rad]　
    !            K(3) - Eccentricity
    !            K(4) - Semi-major Axis[RE]　
    !            K(5) - Argument of Periapsis [rad]
    !            K(6) - True Anomaly [rad]
    !******************************************************************

    !******************************************************************
    !   Car2RSW
    !******************************************************************
    !   Input  : X(1:3)    - Position vector of debris [RE]
    !            V(1:3)    - Velocity vector of debris [normalized]
    !            A(1:3)    - Acceleration vector of debris [normalized]
    !   Output : Arsw(1:3) - Acceleration vector
    !                        in Gaussian coordinate system [normalized]
    !******************************************************************

    !******************************************************************
    !   Car2NTW
    !******************************************************************
    !   Input  : X(1:3)    - Position vector of debris [RE]
    !            V(1:3)    - Velocity vector of debris [normalized]
    !            A(1:3)    - Acceleration vector of debris [normalized]
    !   Output : Antw(1:3) - Acceleration vector
    !                        in Frenet coordinate system [normalized]
    !******************************************************************
    
    output_day = time*t0/86400.0d0
    R = nolm(X)
    
#ifdef output_x    
    write(fo1,'(5ES16.8)')  output_day, &
         X(1)*r0*1.0d-3, X(2)*r0*1.0d-3, X(3)*r0*1.0d-3, r0*1.0d-3*R
#endif

#ifdef output_v    
    write(fo2,'(5ES16.8)')  output_day, &
         V(1)*v0*1.0d-3, V(2)*v0*1.0d-3, V(3)*v0*1.0d-3, v0*1.0d-3*nolm(V(:))
#endif

#ifdef output_oe    
    call Car2Kep(X,V,K)
    write(fo3,'(7ES16.8)')  output_day, &
         K(1)*rad2deg, K(2)*rad2deg, K(3), K(4)*r0*1.0d-3, K(5)*rad2deg, K(6)*rad2deg
#endif

#ifdef output_sph_gd
    write(fo4,'(4ES16.8)')  output_day, alt_gd, lat_gd*rad2deg, lon_gd*rad2deg
#endif

#ifdef output_sph_gc
    alt_gc = R*r0*1.0d-3
    write(fo5,'(4ES16.8)')  output_day, alt_gc, lat_gc*rad2deg, lon_gc*rad2deg
#endif

#ifdef output_a_geo
    write(fo6,'(5ES16.8)')  output_day, &
         g0*Ageo(1),  g0*Ageo(2),  g0*Ageo(3),  g0*nolm(Ageo)
#endif    

#ifdef output_a_drag    
    write(fo7,'(5ES16.8)')  output_day, &
         g0*Aair(1),  g0*Aair(2),  g0*Aair(3),  g0*nolm(Aair)
#endif

#ifdef output_a_srp    
    write(fo8,'(5ES16.8)')  output_day, &
         g0*Aph(1),   g0*Aph(2),   g0*Aph(3),   g0*nolm(Aph)
#endif

#ifdef output_a_moon    
    write(fo9,'(5ES16.8)')  output_day, &
         g0*Amoon(1), g0*Amoon(2), g0*Amoon(3), g0*nolm(Amoon)
#endif

#ifdef output_a_sun    
    write(fo10,'(5ES16.8)') output_day, &
         g0*Asun(1),  g0*Asun(2),  g0*Asun(3),  g0*nolm(Asun)
#endif

#ifdef output_a_mag    
    write(fo11,'(5ES16.8)') output_day, charge_d*VcrossB(1)/Md, &
         charge_d*VcrossB(2)/Md, charge_d*VcrossB(3)/Md, g0*nolm(VcrossB)/Md
#endif

#ifdef output_a_ecr    
    write(fo12,'(5ES16.8)') output_day, charge_d*Ecr(1)/Md,     &
         charge_d*Ecr(2)/Md,     charge_d*Ecr(3)/Md,     g0*nolm(Ecr)/Md
#endif

#ifdef output_a_ecn    
    write(fo13,'(5ES16.8)') output_day, charge_d*Ecn(1)/Md,     &
         charge_d*Ecn(2)/Md,     charge_d*Ecn(3)/Md,     g0*nolm(dble(Ecn))/Md
#endif

#ifdef output_charge_amount
    write(fo14,'(3ES16.8)') output_day, Ud, charge_d
#endif
    
#ifdef output_reserve2
    
#endif    
    
    return
  end subroutine output  

  
  subroutine open_file(id_debris)
    !---- args
    integer, intent(in) :: id_debris
    !---- vars
    character(1)  :: name0    
    character(2)  :: name00
    character(3)  :: name000
    character(4)  :: name0000
    character(5)  :: name00000
    character(30) :: &
         name001, name002, name003, name004, name005, &
         name006, name007, name008, name009, name010, &
         name011, name012, name013, name014, name015
    character(:), allocatable :: &
         name1,  name2,  name3,  name4,  name5, &
         name6,  name7,  name8,  name9,  name10, &
         name11, name12, name13, name14, name15
    integer :: n(15)
    !---- body
    ! Set the file name
    name001 = 'x'
    n(1)    = len_trim(name001) + 11
    name002 = 'v'
    n(2)    = len_trim(name002) + 11
    name003 = 'oe'
    n(3)    = len_trim(name003) + 11    
    name004 = 'sph_geodetic'
    n(4)    = len_trim(name004) + 11
    name005 = 'sph_geocentric'
    n(5)    = len_trim(name005) + 11
    name006 = 'A_geo'
    n(6)    = len_trim(name006) + 11
    name007 = 'A_drag'
    n(7)    = len_trim(name007) + 11
    name008 = 'A_SRP'
    n(8)    = len_trim(name008) + 11
    name009 = 'A_moon'
    n(9)    = len_trim(name009) + 11
    name010 = 'A_sun'
    n(10)   = len_trim(name010) + 11
    name011 = 'A_mag'
    n(11)   = len_trim(name011) + 11
    name012 = 'A_ecr'
    n(12)   = len_trim(name012) + 11
    name013 = 'A_ecn'
    n(13)   = len_trim(name013) + 11
    name014 = 'charge'
    n(14)   = len_trim(name014) + 11
    name015 = 'reserve2'
    n(15)   = len_trim(name015) + 11    

    ! Assign the file number corresponding the debris' ID
    ! "file name" + "file number"
    if(id_debris <= 9)then
       write(name0,'(i1.1)') id_debris

       n(:) = n(:) + 1
       
       allocate(character(n(1)) ::name1)
       allocate(character(n(2)) ::name2)
       allocate(character(n(3)) ::name3)
       allocate(character(n(4)) ::name4)
       allocate(character(n(5)) ::name5)
       allocate(character(n(6)) ::name6)
       allocate(character(n(7)) ::name7)
       allocate(character(n(8)) ::name8)
       allocate(character(n(9)) ::name9)
       allocate(character(n(10))::name10)
       allocate(character(n(11))::name11)
       allocate(character(n(12))::name12)
       allocate(character(n(13))::name13)
       allocate(character(n(14))::name14)
       allocate(character(n(15))::name15)
       
       name1  = 'result/'//trim(name001)//name0//'.dat'
       name2  = 'result/'//trim(name002)//name0//'.dat'
       name3  = 'result/'//trim(name003)//name0//'.dat'
       name4  = 'result/'//trim(name004)//name0//'.dat'
       name5  = 'result/'//trim(name005)//name0//'.dat'
       name6  = 'result/'//trim(name006)//name0//'.dat'
       name7  = 'result/'//trim(name007)//name0//'.dat'
       name8  = 'result/'//trim(name008)//name0//'.dat'
       name9  = 'result/'//trim(name009)//name0//'.dat'
       name10 = 'result/'//trim(name010)//name0//'.dat'
       name11 = 'result/'//trim(name011)//name0//'.dat'
       name12 = 'result/'//trim(name012)//name0//'.dat'
       name13 = 'result/'//trim(name013)//name0//'.dat'
       name14 = 'result/'//trim(name014)//name0//'.dat'
       name15 = 'result/'//trim(name015)//name0//'.dat'
    else if(10 <= id_debris .and. id_debris <= 99)then
       write(name00,'(i2.2)') id_debris

       n(:) = n(:) + 2
       
       allocate(character(n(1)) ::name1)
       allocate(character(n(2)) ::name2)
       allocate(character(n(3)) ::name3)
       allocate(character(n(4)) ::name4)
       allocate(character(n(5)) ::name5)
       allocate(character(n(6)) ::name6)
       allocate(character(n(7)) ::name7)
       allocate(character(n(8)) ::name8)
       allocate(character(n(9)) ::name9)
       allocate(character(n(10))::name10)
       allocate(character(n(11))::name11)
       allocate(character(n(12))::name12)
       allocate(character(n(13))::name13)
       allocate(character(n(14))::name14)
       allocate(character(n(15))::name15)
       
       name1  = 'result/'//trim(name001)//name00//'.dat'
       name2  = 'result/'//trim(name002)//name00//'.dat'
       name3  = 'result/'//trim(name003)//name00//'.dat'
       name4  = 'result/'//trim(name004)//name00//'.dat'
       name5  = 'result/'//trim(name005)//name00//'.dat'
       name6  = 'result/'//trim(name006)//name00//'.dat'
       name7  = 'result/'//trim(name007)//name00//'.dat'
       name8  = 'result/'//trim(name008)//name00//'.dat'
       name9  = 'result/'//trim(name009)//name00//'.dat'
       name10 = 'result/'//trim(name010)//name00//'.dat'
       name11 = 'result/'//trim(name011)//name00//'.dat'
       name12 = 'result/'//trim(name012)//name00//'.dat'
       name13 = 'result/'//trim(name013)//name00//'.dat'
       name14 = 'result/'//trim(name014)//name00//'.dat'
       name15 = 'result/'//trim(name015)//name00//'.dat'
    else if(100 <= id_debris .and. id_debris <= 999)then
       write(name000,'(i3.3)') id_debris
       
       n(:) = n(:) + 3
       
       allocate(character(n(1)) ::name1)
       allocate(character(n(2)) ::name2)
       allocate(character(n(3)) ::name3)
       allocate(character(n(4)) ::name4)
       allocate(character(n(5)) ::name5)
       allocate(character(n(6)) ::name6)
       allocate(character(n(7)) ::name7)
       allocate(character(n(8)) ::name8)
       allocate(character(n(9)) ::name9)
       allocate(character(n(10))::name10)
       allocate(character(n(11))::name11)
       allocate(character(n(12))::name12)
       allocate(character(n(13))::name13)
       allocate(character(n(14))::name14)
       allocate(character(n(15))::name15)
       
       name1  = 'result/'//trim(name001)//name000//'.dat'
       name2  = 'result/'//trim(name002)//name000//'.dat'
       name3  = 'result/'//trim(name003)//name000//'.dat'
       name4  = 'result/'//trim(name004)//name000//'.dat'
       name5  = 'result/'//trim(name005)//name000//'.dat'
       name6  = 'result/'//trim(name006)//name000//'.dat'
       name7  = 'result/'//trim(name007)//name000//'.dat'
       name8  = 'result/'//trim(name008)//name000//'.dat'
       name9  = 'result/'//trim(name009)//name000//'.dat'
       name10 = 'result/'//trim(name010)//name000//'.dat'
       name11 = 'result/'//trim(name011)//name000//'.dat'
       name12 = 'result/'//trim(name012)//name000//'.dat'
       name13 = 'result/'//trim(name013)//name000//'.dat'
       name14 = 'result/'//trim(name014)//name000//'.dat'
       name15 = 'result/'//trim(name015)//name000//'.dat'
    else if(1000 <= id_debris .and. id_debris <= 9999)then
       write(name0000,'(i4.4)') id_debris
       
       n(:) = n(:) + 4
       
       allocate(character(n(1)) ::name1)
       allocate(character(n(2)) ::name2)
       allocate(character(n(3)) ::name3)
       allocate(character(n(4)) ::name4)
       allocate(character(n(5)) ::name5)
       allocate(character(n(6)) ::name6)
       allocate(character(n(7)) ::name7)
       allocate(character(n(8)) ::name8)
       allocate(character(n(9)) ::name9)
       allocate(character(n(10))::name10)
       allocate(character(n(11))::name11)
       allocate(character(n(12))::name12)
       allocate(character(n(13))::name13)
       allocate(character(n(14))::name14)
       allocate(character(n(15))::name15)
       
       name1  = 'result/'//trim(name001)//name0000//'.dat'
       name2  = 'result/'//trim(name002)//name0000//'.dat'
       name3  = 'result/'//trim(name003)//name0000//'.dat'
       name4  = 'result/'//trim(name004)//name0000//'.dat'
       name5  = 'result/'//trim(name005)//name0000//'.dat'
       name6  = 'result/'//trim(name006)//name0000//'.dat'
       name7  = 'result/'//trim(name007)//name0000//'.dat'
       name8  = 'result/'//trim(name008)//name0000//'.dat'
       name9  = 'result/'//trim(name009)//name0000//'.dat'
       name10 = 'result/'//trim(name010)//name0000//'.dat'
       name11 = 'result/'//trim(name011)//name0000//'.dat'
       name12 = 'result/'//trim(name012)//name0000//'.dat'
       name13 = 'result/'//trim(name013)//name0000//'.dat'
       name14 = 'result/'//trim(name014)//name0000//'.dat'
       name15 = 'result/'//trim(name015)//name0000//'.dat'
    else if(10000 <= id_debris .and. id_debris <= 99999)then
       write(name00000,'(i5.5)') id_debris
       
       n(:) = n(:) + 5
       
       allocate(character(n(1)) ::name1)
       allocate(character(n(2)) ::name2)
       allocate(character(n(3)) ::name3)
       allocate(character(n(4)) ::name4)
       allocate(character(n(5)) ::name5)
       allocate(character(n(6)) ::name6)
       allocate(character(n(7)) ::name7)
       allocate(character(n(8)) ::name8)
       allocate(character(n(9)) ::name9)
       allocate(character(n(10))::name10)
       allocate(character(n(11))::name11)
       allocate(character(n(12))::name12)
       allocate(character(n(13))::name13)
       allocate(character(n(14))::name14)
       allocate(character(n(15))::name15)
       
       name1  = 'result/'//trim(name001)//name00000//'.dat'
       name2  = 'result/'//trim(name002)//name00000//'.dat'
       name3  = 'result/'//trim(name003)//name00000//'.dat'
       name4  = 'result/'//trim(name004)//name00000//'.dat'
       name5  = 'result/'//trim(name005)//name00000//'.dat'
       name6  = 'result/'//trim(name006)//name00000//'.dat'
       name7  = 'result/'//trim(name007)//name00000//'.dat'
       name8  = 'result/'//trim(name008)//name00000//'.dat'
       name9  = 'result/'//trim(name009)//name00000//'.dat'
       name10 = 'result/'//trim(name010)//name00000//'.dat'
       name11 = 'result/'//trim(name011)//name00000//'.dat'
       name12 = 'result/'//trim(name012)//name00000//'.dat'
       name13 = 'result/'//trim(name013)//name00000//'.dat'
       name14 = 'result/'//trim(name014)//name00000//'.dat'
       name15 = 'result/'//trim(name015)//name00000//'.dat'       
    end if

#ifdef output_x
    open(fo1,  file = name1,  status = 'unknown')
#endif
#ifdef output_v    
    open(fo2,  file = name2,  status = 'unknown')
#endif
#ifdef output_oe    
    open(fo3,  file = name3,  status = 'unknown')
#endif
#ifdef output_sph_gd    
    open(fo4,  file = name4,  status = 'unknown')
#endif
#ifdef output_sph_gc    
    open(fo5,  file = name5,  status = 'unknown')
#endif
#ifdef output_a_geo    
    open(fo6,  file = name6,  status = 'unknown')
#endif
#ifdef output_a_drag    
    open(fo7,  file = name7,  status = 'unknown')
#endif
#ifdef output_a_srp    
    open(fo8,  file = name8,  status = 'unknown')
#endif
#ifdef output_a_moon    
    open(fo9,  file = name9,  status = 'unknown')
#endif
#ifdef output_a_sun    
    open(fo10, file = name10, status = 'unknown')
#endif
#ifdef output_a_mag    
    open(fo11, file = name11, status = 'unknown')
#endif
#ifdef output_a_ecr    
    open(fo12, file = name12, status = 'unknown')
#endif
#ifdef output_a_ecn    
    open(fo13, file = name13, status = 'unknown')
#endif
#ifdef output_charge_amount
    open(fo14, file = name14, status = 'unknown')
#endif
#ifdef output_reserve2 
    open(fo15, file = name15, status = 'unknown')
#endif        
    
    return
  end subroutine open_file
  

  subroutine close_file
    !---- args
    !---- vars
    !---- body
#ifdef output_x    
    close(fo1)
#endif
#ifdef output_v    
    close(fo2)
#endif
#ifdef output_oe    
    close(fo3)
#endif
#ifdef output_sph_gd    
    close(fo4)
#endif
#ifdef output_sph_gc    
    close(fo5)
#endif
#ifdef output_a_geo    
    close(fo6)
#endif
#ifdef output_a_drag    
    close(fo7)
#endif
#ifdef output_a_srp    
    close(fo8)
#endif
#ifdef output_a_moon    
    close(fo9)
#endif
#ifdef output_a_sun    
    close(fo10)
#endif
#ifdef output_a_mag    
    close(fo11)
#endif
#ifdef output_a_ecr    
    close(fo12)
#endif
#ifdef output_a_ecn    
    close(fo13)
#endif
#ifdef output_charge_amount
    close(fo14)
#endif
#ifdef output_reserve2
    close(fo15)
#endif        
    
    return
  end subroutine close_file     
end module output_mod
