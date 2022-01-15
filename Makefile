#*************************************************************************
# Before input the 'make' command, select the kinds of perturbing force
# and the environmental model in consideration of the orbital calculation,
# and choose parameters you want to output.
#*************************************************************************

debug   = off
#**************************************************************************
#   select the perturbing force in consideration of the orbital calcultaion
#**************************************************************************
airdrag = on   # Atmospheric Drag
	       # on: Considering atmospheric drag.
	       #     In LEO region (< 2000 km),
	       #     atmospheric drag should be included.
mag     = on   # Perturbing force derived from the geomagnetic field.
	       # (ignoring the rotation of geomagnetic field)
	       # (if you want to consider the Lorentz force,
	       #  both flags of "mag" and "ecor" should be turned on.)
	       # on: Considering.
	       #     This effects are remarkable for HHMR or small objects.
ecor    = on   # Perturbing force derived from the corotational electric field.
	       # (this is the effect of the rotation of the geomagnetic field.)
	       # on: Considering.
	       #     This effects are remarkable for HHMR or small objects.
econv   = off  # Perturbing force derived from the convection electric field.
	       # on: Considering.
	       #     This effects are remarkable for HHMR or small objects.

#**************************************************************************
#   select the environmental model for the simulation
#**************************************************************************
charge  = 0      # 0:  Assuming Charge as a Constance
	         # 1:  Considering charge variation during orbital motion.
	         #     Charge amount is calculated from the charge database.
		 # 2:  Considering charge variation during orbital motion.
		 #     Charge amount and orbital motion are calculated
		 #     at the same time.
mag_model = 3    # Select the geomagnetic field model
	         # 0:  Exclude
	         # 1:  IGRF
	         # 2:  Ideal Dipole (ignoring tilts of geomagnetic field)
	         # 3:  Tilted Dipole
air_model = 1    # Select the atmospheric density model
	         # 0:  Exponential Model
	         # 1:  Jacchia-Roberts Model
	         # 2:  NLRMSISE00 Model
shadow_model = 1 # Select the Earth's Shadow Model
		 # 1:  Cylindrical Shadow Model
	         # 2:  Conical Shadow Model

#**************************************************************************
#   select output parameters
#**************************************************************************
out_x      = yes  # Position Vector
		  # 1: time [day]
		  # 2: x [km]
		  # 3: y [km]
		  # 4: z [km]
		  # 5: geocentric altitude [km]
out_v      = yes  # Velocity Vector
		  # 1: time [day]
		  # 2: Vx [km/s]
		  # 3: Vy [km/s]
		  # 4: Vz [km/s]
		  # 5: velocity [km/s]
out_oe     = yes  # Orbital Elements
		  # 1: time [day]
		  # 2: orbital inclination [deg.],
		  # 3: right ascension of the ascending node [deg.]
		  # 4: eccentricity
		  # 5: semi-major axis [km]
		  # 6: argument of perigee[deg.]
		  # 7: true anomaly [deg.]
out_sph_gd = no  # Geodetic Spherical Components
		  # 1: time [day]
		  # 2: geodetic altitude [km]
		  # 3: geodetic latitude [deg.]
		  # 4: geodetic longitude [deg.]
out_sph_gc = no  # Geocentric Spherical Components
		  # 1: time [day]
		  # 2: geocentric distance from the center of the Earth [km]
		  # 3: declination [deg.]
		  # 4: right ascension [deg.]
out_a_geo  = no  # Perturbing acceleration derived from
	 	  # heigher harmonics of the Earth gravity
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_a_drag = no  # Perturbing acceleration derived from
		  # atmospheric drag
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_a_srp  = no  # Perturbing force derived from
		  # solar radiation pressure
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_a_moon = no  # Perturbing force derived from
		  # 3rd body effect of the Moon
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_a_sun  = no  # Perturbing force derived from
		  # 3rd body effect of the Sun
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_a_mag  = yes  # Perturbing force derived from
		  # geomagnetic feild
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_a_ecr  = no  # Perturbing force derived from
		  # corotation electri field
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_a_ecn  = no  # Perturbing force derived from
		  # convection electric field
		  # 1: time [day]
		  # 2: Ax [m/s^2]
		  # 3: Ay [m/s^2]
		  # 4: Az [m/s^2]
out_charge = yes  # Charge Amount [C]
out_reserve2 = no

#****************************************
# Definition for specific project
#****************************************
TARGET = go
OBJECTS =  irisub.o          irifun.o         iriflip.o       iritec.o \
	   iridreg.o         igrf.o           cira.o          geopack.o \
	   charge_solver.o   global.o         charge.o        weimer.o \
	   igrf12.o          sub.o            exponential.o   jacchia.o \
	   time.o            coordinates.o    em.o            input.o \
	   output.o          geopotential.o   airdrag.o       moon.o \
	   sun.o             accel.o          ABM.o           main.o
COMMON_MOD = charge_solver.mod \
	     global.mod        charge.mod    weimer.mod   sub.mod \
	     exponential.mod   jacchia.mod   time.mod     coordinates.mod \
	     em.mod            input.mod     output.mod   geopotential.mod \
	     airdrag.mod       moon.mod      sun.mod      accel.mod \
	     ABM.mod


#****************************************
# Default rule for Fortran
#****************************************
##############
# The compiler
FC = gfortran
##############

#################
# Compiler Option
ifeq "$(strip $(debug))" "on"
  # For Debugging
  FFLAGS = -Wall -fbounds-check -fbacktrace
#####################################################################
#   デバッグ用のオプション
#####################################################################
# -warn all
# 全てのコンパイル時警告メッセージを有効
#
# -warn declarations
# 暗黙の型宣言を警告
#
# -check uninit
# 初期化されていない変数を検出
#
# -traceback -g
# デバッグのために必要な情報をオブジェクトファイルに埋め込みます。
# Segmentation Faultなどのエラー終了時にエラーの発生箇所を
# 表示するために必要です。
#
# -CB -traceback -g
# 実行時に配列の領域外参照を検出します。
# 3つのオプションを同時に指定してください。
#
# -fpe{0|1|3} -traceback -g
# 浮動小数点演算の例外処理をどのように扱うかを指定します。
# 浮動小数点演算の例外処理は次の4つです。
# オーバーフロー、アンダーフロー、ゼロ割り、不正な処理(Invalid)
#
# -fpe0 -traceback -g
# オーバーフロー、ゼロ割り、不正な処理を検出するとエラー終了
# アンダーフローはゼロに置き換えて処理続行
#
# -fpe1 -traceback -g
# オーバーフロー、ゼロ割り、不正な処理を検出しても処理続行
# アンダーフローはゼロに置き換えて処理続行
#
# -r8
# real/compelx型で宣言された変数を
# real*8/complex*16型の変数として取り扱います。
#
# -i8
# integer型で宣言された変数をinteger*8型の変数として取り扱います。
######################################################################

else
  FFLAGS = -O1 -cpp
######################################################################
#   Perturbing Force Considered in the Orbital Calculation
######################################################################
  ifeq "$(strip $(airdrag))" "on"
    FFLAGS += -Ddrag
    ifeq "$(strip $(air_model))" "0"
	FFLAGS += -DExpo_model
    endif
    ifeq "$(strip $(air_model))" "1"
	FFLAGS += -DJR_model
    endif
    ifeq "$(strip $(air_model))" "2"
	FFLAGS += -DMSIS_model
    endif
  endif

  ifeq "$(strip $(mag))" "on"
    FFLAGS += -Dmag
  endif

  ifeq "$(strip $(ecor))" "on"
    FFLAGS += -Decor
  endif

  ifeq "$(strip $(econv))" "on"
    FFLAGS += -Deconv
  endif

######################################################################
#   The Environmental Model
######################################################################
  ifeq "$(strip $(shadow_model))" "1"
    FFLAGS += -Dcylindrical_shadow_model
  endif

  ifeq "$(strip $(shadow_model))" "2"
    FFLAGS += -Dconical_shadow_model
  endif

  ifeq "$(strip $(charge))" "1"
    FFLAGS += -Dchange_in_charge
  endif

  ifeq "$(strip $(charge))" "2"
    FFLAGS += -Dcharge_solve
  endif

  ifeq "$(strip $(mag_model))" "1"
    FFLAGS += -DIGRF
  endif

  ifeq "$(strip $(mag_model))" "2"
    FFLAGS += -Dideal_dipole
  endif

  ifeq "$(strip $(mag_model))" "3"
    FFLAGS += -Dtilted_dipole
  endif

######################################################################
#   Output Parameters
######################################################################
  ifeq "$(strip $(out_x))" "yes"
    FFLAGS += -Doutput_x
  endif
  ifeq "$(strip $(out_v))" "yes"
    FFLAGS += -Doutput_v
  endif
  ifeq "$(strip $(out_oe))" "yes"
    FFLAGS += -Doutput_oe
  endif
  ifeq "$(strip $(out_sph_gd))" "yes"
    FFLAGS += -Doutput_sph_gd
  endif
  ifeq "$(strip $(out_sph_gc))" "yes"
    FFLAGS += -Doutput_sph_gc
  endif
  ifeq "$(strip $(out_a_geo))" "yes"
    FFLAGS += -Doutput_a_geo
  endif
  ifeq "$(strip $(out_a_drag))" "yes"
    FFLAGS += -Doutput_a_drag
  endif
  ifeq "$(strip $(out_a_srp))" "yes"
    FFLAGS += -Doutput_a_srp
  endif
  ifeq "$(strip $(out_a_moon))" "yes"
    FFLAGS += -Doutput_a_moon
  endif
  ifeq "$(strip $(out_a_sun))" "yes"
    FFLAGS += -Doutput_a_sun
  endif
  ifeq "$(strip $(out_a_mag))" "yes"
    FFLAGS += -Doutput_a_mag
  endif
  ifeq "$(strip $(out_a_ecr))" "yes"
    FFLAGS += -Doutput_a_ecr
  endif
  ifeq "$(strip $(out_a_ecn))" "yes"
    FFLAGS += -Doutput_a_ecn
  endif
  ifeq "$(strip $(out_charge))" "yes"
    FFLAGS += -Doutput_charge_amount
  endif
  ifeq "$(strip $(out_reserve2))" "yes"
    FFLAGS += -Doutput_reserve2
  endif
endif

######################################################################

######################################################################
# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o : %.f90
	${FC} ${FFLAGS} -c $<

%.o : %.f
	${FC} ${FFLAGS} -c $<

%.o : %.for
	${FC} ${FFLAGS} -c $<
######################################################################

###################################################################
# This has the advantage of not recompiling unnecessarily,
# but if the .mod file is missing then make will be unable to build
# at all.
%.mod : %.f90 %.o
	@:
###################################################################

########################################
# "make" builds the Target
${TARGET} : ${OBJECTS}
	${FC} ${FFLAGS} -o $@ ${OBJECTS}
########################################

################################################
# Utility targets
.PHONY : clean
clean :
	rm -f *.o *.mod  go fort.* *__genmod.f90
################################################
