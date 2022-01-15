module Exponential_Model
  implicit none
  real(8), parameter :: nominal_density(1:28) = &
       (/ 3.019d-15,  5.245d-15,  1.170d-14,  3.614d-14, &
          1.454d-13,  6.967d-13,  1.585d-12,  3.725d-12, &
          9.518d-12,  2.418d-11,  7.248d-11,  2.789d-10, &
          5.464d-10,  2.070d-9,   3.845d-9,   8.484d-9, &
          2.438d-8,   9.661d-8,   5.297d-7,   3.396d-6, &
          1.905d-5,   8.77d-5,    3.206d-4,   1.057d-3, &
          3.972d-3,   1.774d-2,   3.899d-2,   1.225d0 /)
  real(8), parameter :: base_altitude(1:28) = &
       (/ 1.0d6,      9.0d5,      8.0d5,      7.0d5, &
          6.0d5,      5.0d5,      4.5d5,      4.0d5, &
          3.5d5,      3.0d5,      2.5d5,      2.0d5, &
          1.8d5,      1.5d5,      1.4d5,      1.3d5, &
          1.2d5,      1.1d5,      1.0d5,      9.0d4, &
          8.0d4,      7.0d4,      6.0d4,      5.0d4, &
          4.0d4,      3.0d4,      2.5d4,      0.0d0 /)
  real(8), parameter :: scale_height(1:28) = &
       (/ 2.6800d5,   1.8105d5,   1.2464d5,   8.8667d4, &
          7.1835d4,   6.3822d4,   6.0828d4,   5.8515d4, &
          5.3298d4,   5.3628d4,   4.5546d4,   3.7105d4, &
          2.9740d4,   2.2523d4,   1.6149d4,   1.2636d4, &
          9.473d3,    7.263d3,    5.877d3,    5.382d3, &
          5.799d3,    6.549d3,    7.714d3,    8.382d3, &
          7.554d3,    6.682d3,    6.349d3,    7.249d3 /)
contains
  !****************************Exponential Model******************************!
  !  This simple, static model assumes the density of the atmosphere decays   !
  ! exponentially with increasing altitude. It also assumes a spherically     !
  ! symmetrical distribution of particles.                                    !
  !  Although a very simple approach, this method yields moderate results for !
  ! general studies. Source : Wertz, 1978:820, which uses the U.S Standard    !
  ! Atmosphere for 0km, CIRA-72 for 25-500km, and CIRA-72 with exospheric     !
  ! temperature, exospheric temperature, T_infinite = 1000K for 500-1000km.   !
  ! The scale heights have been adjusted to maintain a piecewise-continuous   !
  ! formulation of the density.                                               !
  !                                                                           !
  ! Input     : geodetic altitude [m]                                         !
  ! Output    : atmospheric_density [kg/m^3]                                  !
  ! Author    : Keisuke Akari                                                 !
  ! Date      : October 9, 2017                                               !
  ! Reference : Fundamentals of Astrodynamics and Applications, pp562-564     !
  !                          ---Third Edition---                              !
  !***************************************************************************!
  Subroutine Exponential_Density(altitude, atmospheric_density)
    !--- args
    real(8), intent(in)  :: altitude            ![m]
    real(8), intent(out) :: atmospheric_density ![kg/m^3]
    !--- vars
    integer :: i
    !--- body
    if(altitude < base_altitude(28))then
       print*, 'Error in "exponential.f90" : input altitude is incorrect.'
       stop
    else if(base_altitude(1) <= altitude)then
       atmospheric_density &
            = nominal_density(1)*exp((base_altitude(1) - altitude)/scale_height(1))
    else
       do i = 1, 27, 1
          if((base_altitude(i+1) <= altitude).and.(altitude < base_altitude(i)))then
             atmospheric_density &
                  = nominal_density(i+1)*exp((base_altitude(i+1) - altitude)/scale_height(i+1))
             exit
          end if
       end do
    end if    
  end Subroutine Exponential_Density
end module Exponential_Model
