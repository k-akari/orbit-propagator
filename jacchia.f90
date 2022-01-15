module Jacchia_Roberts_mod
  use global_mod, only : fi_flux
  implicit none
contains
!***********************************************************
!***********************************************************
!     LIBRARY DENSITY
!***********************************************************
!-----------------------------------------------------------
!
!     The library DENSITY includes several routines to
!     compute the high atmospheric properties, using the
!     Robert's version of the Jacchia 1970 model.
!
!-----------------------------------------------------------

  subroutine RSDAMO (SA,SU,SF,RJUD,DAFR,TE,AD,WMOL,RHOD)
!***********************************************************    
!    INPUTS:
!       SA(1)    RIGHT ASCENTION OF THE POINT, IN RADIANS.
!       SA(2)    DECLINATION (GEOCENTRIC  LATITUDE)  OF  THE
!                POINT, IN RADIANS (-PI TO PI).
!       SA(3)    GEOCENTRIC ALTITUDE IN METERS, BETWEEN  THE
!                RANGE 110,000-2,000,000.
!       SU(1)    RIGHT ASCENTION OF THE SUN AT THE DATE,  IN
!                RADIANS (0 TO 2*PI).
!       SU(2)    SUN DECLINATION IN RADIANS (-PI TO PI).
!       SF(1)    DAILY OBSERVED SOLAR FLUX AT  10.7  CM,  AT
!                THE  TIME  1.71  DAYS  EARLIER,  IN   1E-22
!                W/M/M/HZ.
!       SF(2)    AVERAGED DAILY OBSERVED FLUX  AS DEFINED BY
!                JACCHIA, IN 1E-22 W/M/M/HZ.
!       SF(3)    3-HOURLY PLANETARY GEOMAGNETIC INDEX KP, AT
!                THE TIME 0.279 DAYS EARLIER.
!       RJUD     MODIFIED JULIAN  DATE
!       DAFR     TIME (UT) OF THE DAY, IN SECONDS.    
!       GSTI     GREENWICH SIDEREAL TIME, IN RADIANS, AT THE
!                TIME DAFR OF THE DATE RJUD (0 TO 2*PI).(NOT
!                USED. FOR COMPATIBILITY PURPOSE WITH  OTHER
!                MODELS ONLY)
!                ＊このサブルーチンではGSTIは使われていなかっ
!                  たので、消去した。(2017 6/14 Akari)    
!
! OUTPUTS:
!       TE(1)    EXOSPHERIC TEMPERATURE ABOVE THE  POINT  AS
!                DEFINED BY JACCHIA'S 70 MODEL, IN KELVIN.
!       TE(2)    LOCAL  TEMPERATURE  AROUND  THE  POINT,  IN
!                KELVIN.
!       AD(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
!       AD(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
!       AD(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
!       AD(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
!       AD(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
!       AD(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
!       WMOL     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE  AT
!                THE POINT, IN KG/KGMOL.
!       RHOD     MEAN MASS DENSITY OF THE ATMOSPHERE AT  THE
!                POINT, IN KG/M/M/M.
!***********************************************************
    real(8), intent(in)  :: SA(3), SU(2), SF(3), RJUD, DAFR
    real(8), intent(out) :: TE(2), AD(6), WMOL,  RHOD
    real(8) :: AMJD, AL(6)

    AMJD = RJUD + 33282.0d0 + DAFR/86400.0d0

    call DYJRMO (AMJD,SU,SA,SF,TE,AL,WMOL,RHOD)
  
    AD(1) = AL(3)
    AD(2) = AL(4)
    AD(3) = AL(1)
    AD(4) = AL(2)
    AD(5) = AL(5)
    AD(6) = AL(6)
    
    return
  end subroutine RSDAMO


  subroutine DYJRMO(DJM,SUN,SAT,GEO,TEMP,DN,AMW,DENS)
!*******************************************************
! INPUTS : DJM ... MODIFIED JULIAN DATE DJM=JD-2400000.5
!          SUN(1). RIGHT ASCENSION OF SUN (RAD)
!          SUN(2). DECLINATION     OF SUN (RAD)
!          SAT(1). RIGHT ASCENSION OF THE POINT (RAD)
!          SAT(2). DECLINATION     OF THE POINT (RAD)
!          SAT(3). ALTITUDE        OF THE POINT (M)
!          GEO(1). 10.7 CM SOLAR FLUX,IN UNITS OF
!                  1.E-22 WATTS M**2/HERTZ , FOR A
!                  TABULAR TIME 1.71 DAYS EARLIER
!          GEO(2). 10.7 CM SOLAR FLUX AVERAGED OVER
!                  FOUR SOLAR ROTATIONS,CENTERED ON
!                  THE PRESENT TIME
!          GEO(3). GEOMAGNETIC PLANETARY THREE HOUR
!                  RANGE  INDEX "KP" FOR A TABULAR
!                  TIME 0.279 DAYS EARLIER
!
! OUTPUTS : TEMP(1)...EXOSPHERIC TEMPERATURE (KELVIN)
!           TEMP(2)...LOCAL TEMPERATURE      (KELVIN)
!           DN(1) ... LOG10 OF N2 DENSITY NUMBER (M**-3)
!           DN(2) ... LOG10 OF A  DENSITY NUMBER
!           DN(3) ... LOG10 OF HE DENSITY NUMBER
!           DN(4) ... LOG10 OF O2 DENSITY NUMBER
!           DN(5) ... LOG10 OF O  DENSITY NUMBER
!           DN(6) ... LOG10 OF H  DENSITY NUMBER
!           AMW   ... MEAN MOLECULAR WEIGHT (KG/KGMOL)
!           DENS  ... ATMOSPHERIC DENSITY (KG/M**3)
!******************************************************
    real(8), intent(in)  :: DJM,     SUN(2), SAT(3), GEO(3)
    real(8), intent(out) :: TEMP(2), DN(6),  AMW,    DENS
    real(8), parameter   :: PIV2   = 6.2831853d0
    real(8), parameter   :: PIV4   = 12.566371d0
    real(8), parameter   :: PID4   = 0.78539816d0
    real(8), parameter   :: CONS25 = 0.35355339d0
    real(8) :: SUN1,  SUN2,  SAT1,  SAT2,  SAT3,   FS,    FSM   
    real(8) :: PK,    TSUBC, ETA,   THETA, H,      TAU,   S
    real(8) :: DF,    TSUBL, EXPKP, DTG18, DTG20,  DLR20, F
    real(8) :: DLRGM, DTG,   TINF,  TZ,    CAPPHI, GDFT,  FDFZ
    real(8) :: DLRSA, DLRSL, DLR,   DLHE,  D1,     D2,    D3
    real(8) :: D4,    D5,    D6,    SUMN,  SUMNM
!
!      PIV2 = 2 * PI
!      PIV4 = 4 * PI
!      PID4 = PI / 4
!      CONS25 = SIN (PI/4) **3
!
    SUN1 = SUN(1)
    SUN2 = SUN(2)
    SAT1 = SAT(1)
    SAT2 = SAT(2)
    SAT3 = 1.0d-3*SAT(3)
    FS   = GEO(1)
    FSM  = GEO(2)
    PK   = GEO(3)
!
!       MINIMUM NIGHT-TIME TEMPERATURE OF THE GLOBAL
!       EXOSPHERIC TEMPERATURE DISTRIBUTION WHEN THE
!       GEOMAGNETIC ACTIVITY INDEX KP = 0
!       EQUATION 14J
!
    TSUBC = 379.0d0 + 3.24d0*FSM + 1.3d0*(FS - FSM)
!
!       EQUATION 15J
!
    ETA   = 0.5d0*abs(SAT2-SUN2)
    THETA = 0.5d0*abs(SAT2+SUN2)
!
!       EQUATION 16J
!
    H   = SAT1 - SUN1
    TAU = H - 0.64577182d0 + 0.10471976d0*sin(H + 0.75049158d0)
!
!       EXOSPHERIC TEMPERATURE TSUBL WITHOUT CORRECTION
!       FOR GEOMAGNETIC ACTIVITY
!       EQUATION 17J
!
    S     = sin(THETA)**2.2d0
    DF    = S + (cos(ETA)**2.2d0 - S) * abs(cos(0.5d0*TAU))**3
    TSUBL = TSUBC*(1.0d0 + 0.3d0*DF)
!
!       EQUATION 18J
!
    EXPKP  = exp(PK)
    DTG18  = 28.0d0*PK + 0.03d0*EXPKP
!
!       EQUATION 20J
!
    DTG20 = 14.0d0*PK + 0.02d0*EXPKP
    DLR20 = 0.012d0*PK + 1.2d-05*EXPKP
!
!       THE FOLLOWING STATEMENTS EFFECT A CONTINUOUS
!       TRANSITION FROM EQ. 20J AT HEIGHTS WELL BELOW
!       350 KM TO EQ. 18J AT HEIGHTS WELL ABOVE
!       350 KM .
!
    F     = 0.5d0*(tanh(0.04d0*(SAT3 - 350.0d0)) + 1.0d0)
    DLRGM = DLR20 * (1.0d0 - F)
    DTG   = DTG20 * (1.0d0 - F) + DTG18 * F
!
!       EXOSPHERIC TEMPERATURE
!
    TINF  = TSUBL + DTG
!
!   STATIC MODEL OF JACCHIA-ROBERTS FOR THE
!   ATMOSPHERIC DENSITY
!
    CALL STJRMO(TINF,SAT3,TZ,DN)
!
!   EQ. 23J   PHASE OF THE SEMI-ANNUAL VARIATION
!
    CAPPHI = mod((DJM - 36204.0d0)/365.2422d0, 1.0d0)
!
!   EQ. 22J
!
    TAU = CAPPHI + 0.09544d0*((0.5d0 &
         + 0.5d0*sin(PIV2*CAPPHI + 6.035d0))**1.650d0 - 0.5d0)
    GDFT = 0.02835d0 + 0.3817d0*(1.0d0 &
         + 0.4671d0*sin(PIV2*TAU + 4.137d0))*sin(PIV4*TAU + 4.259d0)
    FDFZ = (5.876d-07*SAT3**2.331d0 + 6.328d-2)*exp(-2.868d-03*SAT3)
!
!   EQ. 21J  SEMI-ANNUAL VARIATION
!
    DLRSA = FDFZ*GDFT
!
!   EQ. 24J  SEASONAL-LATITUDINAL VARIATION OF THE
!            LOWER THERMOSPHERE
!
    DLRSL = 0.014d0*(SAT3 - 90.0d0) &
         * exp(-0.0013d0*(SAT3 - 90.0d0)**2) &
         * sign(1.0d0, SAT2)*sin(PIV2*CAPPHI + 1.72d0) &
         * sin(SAT2)**2
!
!   SUM THE CORRECTIONS AND APPLY TO THE
!   NUMBER DENSITIES
!
    DLR = DLRGM + DLRSA + DLRSL
    DN(1) = DN(1) +  DLR
    DN(2) = DN(2) +  DLR
    DN(3) = DN(3) +  DLR
    DN(4) = DN(4) +  DLR
    DN(5) = DN(5) +  DLR
    DN(6) = DN(6) +  DLR
!
!   EQ. 25J  SEASONAL-LATITUDINAL VARIATION
!            OF HELIUM
!
    DLHE = 0.65d0*abs(SUN2 / 0.4091609d0) &
         * (sin(PID4 - 0.5d0*SAT2*sign(1.0d0, SUN2))**3 - CONS25)
    DN(3)   = DN(3) + DLHE
!
!  COMPUTE DENSITY AND MEAN MOLECULAR WEIGHT
!
    D1 = 10.0d0**DN(1)
    D2 = 10.0d0**DN(2)
    D3 = 10.0d0**DN(3)
    D4 = 10.0d0**DN(4)
    D5 = 10.0d0**DN(5)
    D6 = 10.0d0**DN(6)
    SUMN  = D1 + D2 + D3 + D4 + D5 + D6
    SUMNM = 28.0134d0*D1 + 39.9480d0*D2 + 4.0026d0*D3 &
         + 31.9988d0*D4 + 15.9994d0*D5 + 1.00797d0*D6
    AMW   = SUMNM/SUMN
    DENS  = SUMNM/6.02257d+26
    TEMP(1) = TINF
    TEMP(2) = TZ

    return
  end subroutine DYJRMO


  subroutine STJRMO(TINF,SAT3,TZ,DN)
!****************************************************
!  INPUTS : TINF...EXOSPHERIC TEMPERATURE (KELVIN)
!           SAT3...ALTITUDE               (KM)
!
!  OUTPUTS : TZ  ... LOCAL TEMPERATURE    (KELVIN)
!            DN  ... LOG10 OF DENSITY NUMBERS (M**-3)
!            DN(1) ... N2
!            DN(2) ... A
!            DN(3) ... HE
!            DN(4) ... O2
!            DN(5) ... O
!            DN(6) ... H
!****************************************************
    real(8), intent(in)  :: TINF, SAT3
    real(8), intent(out) :: TZ,   DN(6)
!-----
    if(SAT3 > 125.0d0)then
       call STJR03(TINF,SAT3,TZ,DN)
       return
    else if(SAT3 > 100.0d0)then
       call STJR02(TINF,SAT3,TZ,DN)
       return
    else if(SAT3 < 90.0d0)then
       print*, 'A flying object is reentering the atmosphere'
       print*, 'Altitude is', SAT3
       stop
    end if
    call STJR01(TINF,SAT3,TZ,DN)
    
    return
  end subroutine STJRMO


  subroutine STJR01(TINF,SAT3,TL2,AL10N)
!*****************************************************    
!  INPUTS : TINF ... EXOSPHERIC TEMPERATURE (KELVIN)
!           SAT3 ... ALTITUDE               (KM)
!
!  OUTPUTS : TL2 ... LOCAL TEMPERATURE (KELVIN)
!            AL10N . ALOG10 OF DENSITY NUMBERS (M**-3)
!
!            AL10N(1) ... N2
!            AL10N(2) ... A
!            AL10N(3) ... HE
!            AL10N(4) ... O2
!            AL10N(5) ... O
!            AL10N(6) ... H
!*****************************************************
    real(8), intent(in)  :: TINF, SAT3
    real(8), intent(out) :: TL2,  AL10N(6)
!
!   RA = POLAR EARTH RADIUS (KM)
!   WT = WEIGHTS FOR THE NEWTON-COTES
!        FIVE POINT QUADRATURE FORMULAE
!
    integer :: i, j, N
    real(8), parameter :: RA = 6356.766d0
    real(8), parameter :: WT(1:5) = &
         (/ 0.31111111d0,  1.4222222d0, 0.53333333d0, &
            1.4222222d0,   0.31111111d0 /)
    real(8), parameter :: AM1   = 28.82678d0
    real(8), parameter :: TL1   = 183.0d0
    real(8), parameter :: FACT1 = 0.12027444181d0
    real(8) :: TX,  GX, AL,   ZR,   ZEND, SUM2
    real(8) :: AIN, Z,  DZ,   SUM1, ZD,   AM2
    real(8) :: DZX, GZ, DENS, ANM,  AN,   FACT2
    
    TX    = 371.668d0 + 0.0518806d0*TINF &
              - 294.3503d0*exp(-0.00216222d0*TINF)
    GX    = 0.054285714d0*(TX - 183.0d0)
    AL    = log(SAT3/90.0d0)
    N     = int(AL/0.050d0) + 1
    ZR    = exp(AL/DFLOAT(N))
    ZEND  = 90.0d0
    SUM2  = 0.0d0
    AIN   = AM1*9.534750028d0/TL1

    do i = 1, N
       Z = ZEND
       ZEND = ZR*Z
       DZ = 0.25d0*(ZEND-Z)
       SUM1 = 0.31111111d0*AIN
       do j = 2, 5
          Z = Z + DZ
!
!       MOLE!ULAR WEIGHT FOR Z BETWEEN 90 KM AND
!       100 KM . ACCORDING TO JACCHIA 1971,EQ.1J
!
          ZD  = Z - 90.0d0
          AM2 = 28.82678d0 &
                  - 7.40066d-2 * ZD &
                  + ZD*(-1.19407d-2*ZD &
                  + ZD*( 4.51103d-4*ZD &
                  + ZD*(-8.21895d-6*ZD &
                  + ZD*( 1.07561d-5*ZD &
                  - 6.97444d-7*ZD*ZD ))))
!
!       TEMPERATURE FOR Z BETWEEN 90 AND 100 KM
!       EQ. 5R
!
          DZX  = Z - 125.0d0
          TL2  = TX &
                 + ((-9.8204695d-6*DZX &
                 - 7.3039742d-4)*DZX*DZX &
                 + 1.0d0)*DZX*GX
          GZ   = 9.80665d0*(RA/(RA + Z))**2
          AIN  = AM2*GZ/TL2
          SUM1 = SUM1 + WT(J)*AIN
       end do
       SUM2 = SUM2 + DZ*SUM1
    end do
   
    DENS  = 3.46d-6*AM2*TL1*exp(-FACT1*SUM2)/AM1/TL2
    ANM   = 6.02257d26*DENS
    AN    = ANM/AM2
    FACT2 = ANM/28.960d0
    
    AL10N(1) = log10(0.78110d0*FACT2)
    AL10N(2) = log10(9.3432d-3*FACT2)
    AL10N(3) = log10(6.1471d-6*FACT2)
    AL10N(4) = log10(1.20955d0*FACT2-AN)
    AL10N(5) = log10(2.0d0*(AN - FACT2))
    AL10N(6) = AL10N(5)-15.0d0
   
    return
  end subroutine STJR01


  subroutine STJR02(TINF,SAT3,TZ,DN)
    real(8), intent(in)  :: TINF, SAT3
    real(8), intent(out) :: TZ,   DN(6)
!-----
!       R ... UNIVERSAL GAS CONSTANT(JOULES/K MOLE)
!       RA... POLAR EARTH RADIUS    (KM)
!       RAS.. RA**2                 (KM**2)
!
    integer :: i
    real(8), parameter :: R   = 8.31432d0
    real(8), parameter :: RA  = 6356.766d0
    real(8), parameter :: RAS = 4.04084739788d7
    real(8) :: TSUBX, TXMT0, GSUBX,  SKSF,   C0A
    real(8) :: TEMP,  R1,    R2,     PZ1,    PZ2
    real(8) :: DPZ1,  DPZ2,  R1N,    R2N,    SOMA
    real(8) :: PROD,  DIFE,  X,      X2Y2,   H2
    real(8) :: H3,    H4,    UR1H2,  UR2H3,  WR1
    real(8) :: WR2,   VRA,   CX,     DE100,  T100
    real(8) :: Z,     AM100, DEAVOG, D1,     D2
    real(8) :: D3,    D4,    D5,     UR1,    UR2
    real(8) :: Q1,    Q2,    Q3,     Q4,     Q5
    real(8) :: Q6,    DZX,   AUX,    Y,      AUX1
    real(8) :: AUX2,  F3,    F4,     T100TZ, SKSF34
!
!       DENSITY ANALYTICALLY CALCULATED
!
!         EQ. 9J = EQ. 2R
!
    TSUBX = 371.6678d0 + 0.0518806d0*TINF &
         - 294.3503d0*exp(-0.00216222d0*TINF)
!
!       EQ. 11J
!
    TXMT0 = TSUBX - 183.0d0
    GSUBX = 0.054285714d0*TXMT0
!
!       VALUE OF SMALL K <= SK  AND SMALL F <= SF
!
    SKSF  = 9.80665d0/(R*TXMT0)*1500625.0d0*RAS/0.8d0
!
!       VALUE OF  C0* <= C0A FOR COMPOSING THE
!       FOURTH DEGREE POLYNOMIAL
!
    C0A   =  -87783750.0d0 + 274614375.0d0/TXMT0
!
!       NEWTON-RAPHSON PROCEDURE FOR OBTAINING
!       THE TWO REAL ROOTS OF THE QUARTIC
!       POLYNOMIAL C4*P(Z) , EQ. 10R
!
!               INITIAL GUESSES
!
    TEMP = (TSUBX-300.0d0)/200.0d0
    R1   = 167.77d0 - 3.35d0*TEMP
    R2   =  57.34d0 + 7.95d0*TEMP
    
    do i = 1, 7
       PZ1 = C0A + 3542400.0d0*R1 &
            + R1*(R1*(-52687.5d0 + 340.5d0*R1 - 0.8d0*R1*R1))
       PZ2 = C0A + 3542400.0d0*R2 &
            + R2*(R2*(-52687.5d0 + 340.5d0*R2 - 0.8d0*R2*R2))
       DPZ1 = 3542400.0d0 - 105375.0d0*R1 &
            + R1*(1021.5d0*R1 - 3.2d0*R1*R1)
       DPZ2 = 3542400.0d0 - 105375.0d0*R2 &
            + R2*(1021.5d0*R2 - 3.2d0*R2*R2)
       R1N = R1 - PZ1/DPZ1
       R2N = R2 - PZ2/DPZ2
       if(abs(R1N - R1) < 1.0d-7 .and. abs(R2N - R2) < 1.0d-7)then
          exit
       end if
       R1 = R1N
       R2 = R2N
    end do
    
    R1 = R1N
    R2 = R2N
!
!       COMPLEX ROOTS OR X & X**2+Y**2
!
    SOMA = R1 + R2
    PROD = R1 * R2
    DIFE = R1 - R2
    X    = -0.5d0*(SOMA - 425.625d0)
    X2Y2 = -C0A/(0.8d0*PROD)
!
!       CALCULATE U(R1),U(R2),W(R1),W(R2),CX(CAPITAL X),
!                 AND V(-RA)
!
!  EXPRESSION OF W CORRECTED ACCORDING TO GSFC (NASA,1976)
!
    H2 = R1 + RA
    H3 = R2 + RA
    H4 = RAS + 2.0d0*X*RA + X2Y2
    
    UR1H2 = H2*(R1*R1 - 2.0d0*X*R1 + X2Y2)*DIFE
    UR2H3 = H3*(R2*R2 - 2.0d0*X*R2 + X2Y2)*DIFE
    WR1 = RA + X2Y2/R1
    WR2 = RA + X2Y2/R2
    VRA = H4 * H2 * H3
    CX  = -H4 - H4

    DE100 = (((((0.7026942D-32*TINF*TINF &
         -0.7734110D-28*TINF)*TINF &
         +0.3727894D-24*TINF)*TINF &
         -0.1021474D-20*TINF)*TINF &
         +0.1711735D-17*TINF)*TINF &
         -0.1833490D-14*TINF &
         +0.1985549D-10)
    T100   = TSUBX - 0.94585589d0*TXMT0
    Z      = SAT3
!
!      NUMBER OF PARTICLES PER M**3 AT 100 KM
!
!         D  ... TOTAL NUMBER
!         D1 ... N2 NITROGEN
!         D2 ... AR ARGON
!         D3 ... HE HELIUM
!         D4 ... O2 DIATOMYC OXYGEN
!         D5 ... O  MONOATOMYC OXYGEN
!
    AM100 = 27.6396281382d0
    DEAVOG = DE100*6.02257d29
    D1 = 0.78110d0   * DEAVOG
    D2 = 0.0093432d0 * DEAVOG
    D3 = 6.1471d-6   * DEAVOG
    D4 = (1.20955d0 - 28.96d0/AM100)*DEAVOG
    D5 = 2.0d0*(28.96d0 - AM100)/AM100*DEAVOG
!
!      Q(I) PARAMETERS
!
    UR1 = H2*UR1H2
    UR2 = H3*UR2H3
    Q2  =  1.0d0/UR1
    Q3  = -1.0d0/UR2
    Q5  =  1.0d0/VRA
    Q4  = ( 1.0d0/(PROD*RA) &
         +(RA-X2Y2/RA)/VRA &
         +WR1/UR1H2 &
         -WR2/UR2H3)/CX
    Q6 = -Q5 &
         -2.0d0*(X + RA)*Q4 &
         +1.0d0/UR2H3 &
         -1.0d0/UR1H2
    Q1 = -Q4-Q4-Q3-Q2
!
!       TEMPERATURE FOR Z BETWEEN 100 AND 125 KM
!       EQ. 5R
!
    DZX = Z - 125.0d0
    TZ  = TSUBX +((-9.8204695d-6*DZX &
         - 7.3039742d-4)*DZX*DZX + 1.0d0)*DZX*GSUBX
    
    AUX = Z - 100.0d0
    Y   = SQRT(X2Y2 - X*X)
    AUX1 = Z + RA
    AUX2 = RA + 100.0d0
    
    F3 = DLOG(AUX1/AUX2)*Q1 &
         + log((Z-R1)/(100.0d0 - R1))*Q2 &
         + log((Z-R2)/(100.0d0 - R2))*Q3 &
         + log((Z*Z - 2.0d0*X*Z + X2Y2) &
         /(10000.0d0 - 200.0d0*X + X2Y2))*Q4 
    
    F4 = Q5*AUX/(AUX1*AUX2) &
         + Q6/Y*atan(Y*AUX/(X2Y2 + 100.0d0*Z - (100.0d0 + Z)*X))
!
!      DENSITY NUMBERS D(I) : N2,AR,HE,O2,O,H
!      EQ. 20R
!
    T100TZ = T100/TZ
    SKSF34 = SKSF*(F3+F4)
    
    DN(1) = log10(D1*T100TZ*exp(28.0134d0*SKSF34))
    DN(2) = log10(D2*T100TZ*exp(39.9480d0*SKSF34))
    DN(3) = log10(D3*T100TZ**0.62*exp(4.0026d0*SKSF34))
    DN(4) = log10(D4*T100TZ*exp(31.9988d0*SKSF34))
    DN(5) = log10(D5*T100TZ*exp(15.9994d0*SKSF34))
    
    return
  end subroutine STJR02


  subroutine STJR03(TINF,SAT3,TZ,DN)
    real(8), intent(in)  :: TINF, SAT3
    real(8), intent(out) :: TZ,   DN(6)
!
!-----
!
!       R....UNIVERSAL GAS CONSTANT (JOULES/K MOLE)
!       RA...POLAR EARTH RADIUS     (KM)
!       RAS..RA**2                  (KM**2)
!
    real(8), parameter :: R   = 8.31432d0
    real(8), parameter :: RA  = 6356.766d0
    real(8), parameter :: RAS = 4.04084739788d7
    real(8) :: TSUBX, D1,    D2,    D3,    D4
    real(8) :: D5,    Z,     AL,    TXMT0, TETX
    real(8) :: AUX,   AUX1,  AUX2,  A1,    A2
    real(8) :: A1A2A, TZ500, H500,  D6
!
!       DENSITY ANALYTICALLY CALCULATED
!
!         EQ. 9J = EQ. 2R
!
    TSUBX = 371.6678 + 0.0518806*TINF &
         -294.3503*exp(-0.00216222*TINF)

    D1 = ((((-0.2296182d-19*TINF*TINF &
         + 0.1969715d-15*TINF)*TINF &
         - 0.7139785d-12*TINF)*TINF &
         + 0.1420228d-8*TINF)*TINF &
         - 0.1677341d-5*TINF)*TINF &
         + 0.1186783d-2*TINF &
         + 0.1093155d2
    D2 = ((((-0.4837461d-19*TINF*TINF &
         + 0.4127600d-15*TINF)*TINF &
         - 0.1481702d-11*TINF)*TINF &
         + 0.2909714d-8*TINF)*TINF &
         - 0.3391366d-5*TINF)*TINF &
         + 0.2382822d-2*TINF &
         + 0.8049405d1
    D3 = (((-0.1270838d-16*TINF*TINF &
         + 0.9451989d-13*TINF)*TINF &
         - 0.2894886d-09*TINF)*TINF &
         + 0.4694319d-6*TINF)*TINF &
         - 0.4383486d-3*TINF &
         + 0.7646886d1
    D4 = ((((-0.3131808d-19*TINF*TINF &
         + 0.2698450d-15*TINF)*TINF &
         - 0.9782183d-12*TINF)*TINF &
         + 0.1938454d-8*TINF)*TINF &
         - 0.2274761d-5*TINF)*TINF &
         + 0.1600311d-2*TINF &
         + 0.9924237d1
    D5 = (((0.5116298d-17*TINF*TINF &
         - 0.3490739d-13*TINF)*TINF &
         + 0.9239354d-10*TINF)*TINF &
         - 0.1165003d-6*TINF)*TINF &
         + 0.6118742d-4*TINF &
         + 0.1097083d2
    Z   = SAT3
    AL = ((0.2462708d-9*TINF*TINF &
         - 0.1252487d-5*TINF)*TINF &
         + 0.1579202d-2*TINF)*TINF &
         + 0.2341230d1*TINF &
         + 0.1031445d5
!
!      TEMPERATURE PROFILE EQ. 23R
!
    TXMT0 = TSUBX - 183.0d0
    TETX = TINF - TSUBX
    TZ   = TETX*exp(-TXMT0/TETX*(Z - 125.0d0)/35.0d0*AL/(Z + RA))
!
!      PARAMETERS G(I) : N2,AR,HE,O2,O EQ. 25'R
!
    AUX = 9.80665d0*RAS/(R*AL*TINF)*TETX/TXMT0*35.0d0/6481.766d0
    
    AUX1  = TSUBX/(TINF - TZ)
    AUX2  = TZ / TETX
    A1    = log10(AUX1)
    A2    = log10(AUX2)
    A1A2A = (A1 + A2)*AUX
!
!        DENSITY NUMBERS D(I) EQ.25R
!
    D1 = D1 +        A1 + 28.0134d0*A1A2A
    D2 = D2 +        A1 + 39.9480  *A1A2A
    D3 = D3 + 0.62d0*A1 + 4.0026d0 *A1A2A
    D4 = D4 +        A1 + 31.9988d0*A1A2A
    D5 = D5 +        A1 + 15.9994d0*A1A2A
!
!       CALCULATE TZ(500) EQ.23R
!       T(Z) ALREADY CALCULATED
!
    TZ500 = TINF &
         - TETX*exp(-TXMT0/TETX*(375.0d0/35.0d0*AL/(500.0d0 + RA)))
!
!       INCLUSION OF HYDROGEN
!
!       DENSITY NUMBER FROM EQS. 26R,27R
!
    AUX1 = log10(TINF)
    H500 = 73.13d0-(39.4d0 - 5.5d0*AUX1)*AUX1
    A1   = log10(TZ500/(TINF - TZ))
    A2   = log10(TZ/(TINF-TZ500))
    D6   = H500 + A1 + 1.00797d0*AUX*(A1 + A2)
!
!        LOAD ALOG10 OF DENSITY NUMBERS
!        IN M**-3
!
    DN(1) = D1 + 6.0d0
    DN(2) = D2 + 6.0d0
    DN(3) = D3 + 6.0d0
    DN(4) = D4 + 6.0d0
    DN(5) = D5 + 6.0d0
    DN(6) = D6 + 6.0d0
    TZ    = TINF - TZ
    
    return
  end subroutine STJR03
  

  subroutine SFLUXD (RJUD, frag, SD, OUTR)
!**************************************************************
! PURPOSE:
!       THE SUBROUTINE SFLUXD RETURNS WITH  THE  SOLAR  FLUX
!                                                -      ----
!       AND GEO-PHYSICAL DATA AT THE DATE RJUD.
!                        -
! INPUTS:
!       RJUD  MODIFIED JULIAN DATE  (IF  OUT  OF  RANGE  THE 
!             ROUTINE WILL RETURN A NON-ZERO VALUE IN OUTR.
!			(REAL*8)
!
! OUTPUTS:
!       SD(1) MODIFIED JULIAN DATE (1950.0) OF THE INPUT RECORD
!       SD(2) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 0-3 UTC.
!       SD(3) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 3-6 UTC.
!       SD(4) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 6-9 UTC.
!       SD(5) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 9-12 UTC.
!       SD(6) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 12-15 UTC.
!       SD(7) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 15-18 UTC.
!       SD(8) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 18-21 UTC.
!       SD(9) THREE-HOURLY GEOMAGNETIC PLANETARY INDICE  KP,
!             FOR 21-24 UTC.
!       SD(10) DAILY SUMM OF  PLANETARY  GEOMAGNETIC  INDICES
!             KP
!       SD(11) DAILY EQUIVALENT PLANETARY AMPLITUDE AP
!       SD(12) DAILY OBSERVED SOLAR FLUX AT  2800  MHZ  (10.7
!             CM) IN 1E-22 W/M/M/HZ.
!       SD(13) ADJUSTED TO 1 AU DAILY SOLAR  FLUX,  IN  1E-22 
!             W/M/M/HZ.
!       SD(14) AVERAGED DAILY OBSERVED SOLAR FLUX, AS  DEFI-
!             NED BY JACCHIA (1).
!       SD(15) ARITMETIC MEAN OF DAILY OBSERVED SOLAR  FLUX,
!             AS DEFINED BY JACCHIA (2).
!
!       OUTR  OUTPUT FLAG:
!             0  - NORMAL EXECUTION.
!             -1.- MISSING DATA FILE.
!             1. - MODIFIED JULIAN DATE LESS THAN THE  FIRST
!                  DATE OF THE DATA FILE.
!             2. - MODIFIED JULIAN  DATE  GREATER  THAN  THE 
!                  LAST DATE OF THE DATA FILE.
!
! OBS:
!       THIS ROUTINE USES THE SOLAR  FLUX  AND   GEOMAGNETIC
!       INDEX DATA FILE FLUX_MEAN.DAT
!
! REFERENCES:
!   (1) JACCHIA, L. G.  "THERMOSPHERIC TEMPERATURE,  DENSITY
!       AND COMPOSITIONS: NEW MODELS."  CAMBRIDGE,  MA,  SAO
!       1977. (SAO SPECIAL REPORT NO 375)
!
!   (2) JACCHIA, L. G. ATMOSPHERIC MODELS IN THE REGION FROM
!       110 TO 2000 KM. IN COMMITTEE ON SPACE RESEARCH (COS-
!       PAR). "CIRA 1972." BERLIM, AKADEMIC-VERLAG, 1972.
!
!
! AUTHORS:
!       VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
!
! DATE:
!       JAN. 1985              V. 1.0
!       SEP  2011              V. 2.0
!***********************************************************
    real(8), intent(in)    :: RJUD
    integer, intent(inout) :: frag
    real(8), intent(out)   :: SD(15), OUTR
    integer, save :: TIDA, TEDA
    integer :: i, INDI
    integer :: NATI, NATF, DATI, DATF

    if(frag == 0)then
       inquire(file = 'flux_mean.dat')
       open(fi_flux, file = 'input/flux_mean.dat', access = 'direct', &
            recl = 49, form = 'formatted', status='old')
       read(fi_flux, '(I5,1X,I5)', rec = 1) NATI,NATF
       DATI   = NATI
       DATF   = NATF
       TIDA   = DATI + 2
       TEDA   = DATF

       frag   = 1
    end if

    OUTR    = 0.0d0
    if(RJUD < TIDA)then
       OUTR   = 1.0d0
       do i = 1, 11
          SD(i) = 0.0d0
       end do
       return
    end if
    
    if(RJUD > TEDA)then
       OUTR   = 2.0d0
       do i = 1, 11
          SD(i) = 0.0d0
       end do
       return
    end if
    
    INDI    = RJUD - TIDA + 4

    read(fi_flux, '(F5.0,8F2.0,2F3.0,4F5.1)', rec = INDI) SD
    do i = 2, 10
       SD(i) = 0.1d0*SD(i)
    end do

    return
  end subroutine SFLUXD


  subroutine SFDJ70(RJUD, DAFR, frag, SF, OUTR)
!*************************************************************
! INPUTS:
!       RJUD  MODIFIED JULIAN DATE 
!             referred to 1950.0
!			(REAL*8)
!       DAFR  TIME (UT) OF THE DAY, IN SECONDS (0 TO 86400).
!			(REAL*8)
!
! OUTPUTS:
!       SF(1)    DAILY OBSERVED SOLAR FLUX AT  10.7  CM,  AT
!                THE  TIME  1.71  DAYS  EARLIER,  IN   1E-22
!                W/M/M/HZ.
!       SF(2)    AVERAGED DAILY OBSERVED SOLAR FLUX, AS  DEFI-
!                NED BY JACCHIA (1).
!       SF(3)    3-HOURLY PLANETARY GEOMAGNETIC INDEX KP, AT
!                THE TIME 0.279 DAYS EARLIER.
!
!       OUTR  OUTPUT FLAG:
!             0  - NORMAL EXECUTION.
!             -1.- DATA FILE NOT ON DISK.
!             1. - MODIFIED JULIAN DATE LESS THAN THE  FIRST
!                  DATE OF THE DATA FILE.
!             2. - MODIFIED JULIAN  DATE  GREATER  THAN  THE 
!                  LAST DATE OF THE DATA FILE.
!
! OBS:
!       THIS ROUTINE USES THE SOLAR  FLUX  AND   GEOMAGNETIC
!       INDEX DATA FILE FLUX_MEAN.DAT
!***********************************************************
    real(8), intent(in)    :: RJUD,  DAFR
    integer, intent(inout) :: frag
    real(8), intent(out)   :: SF(3), OUTR
    real(8) :: TAUO, SD(15)
    integer :: ND 

    if(DAFR < 61344.0d0)then
       call SFLUXD(RJUD-2, frag, SD, OUTR)

       SF(1)   = SD(12)
       SF(2)   = SD(15)
       TAUO    = DAFR/3600.0d0 - 6.696d0 ! = 0.279 days
       if (TAUO < 0.0d0)then
          call SFLUXD(RJUD-1, frag, SD, OUTR)
          ND   = int((TAUO + 24)/3) + 2
       else
          call SFLUXD(RJUD, frag, SD, OUTR)
          ND   = int(TAUO/3) + 2
       end if

       SF(3)   = SD(ND)
    else
       call SFLUXD(RJUD-1, frag, SD, OUTR)

       SF(1)   = SD(12)
       SF(2)   = SD(15)
       TAUO    = DAFR/3600.0d0 - 6.696d0
       
       if(TAUO < 0.0d0)then
          ND   = int((TAUO + 24)/3) + 2
       else
          call SFLUXD(RJUD, frag, SD, OUTR)
          ND   = int(TAUO/3) + 2
       end if
       SF(3)   = SD(ND)
    end if

    return
  end subroutine SFDJ70
end module Jacchia_Roberts_mod
