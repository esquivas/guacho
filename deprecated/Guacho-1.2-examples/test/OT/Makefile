#---------------------------------------------------
#   for Guacho 3D
#---------------------------------------------------
#   Name of the executable
PROGRAM=OT

#  Search Path, object files are searched in the order below
VPATH = . : ../src

#   Choice of compiler: ifort and gfortran are tested
#   if MPI enabled, do not use mpif90, but rather the base compiler
#   of the mpi build (see 'mpif90 --version')
#COMPILER= gfortran
COMPILER= ifort

#   Compiler flags, make sure they are compatible with the previous selection
### ifortran ###
#FLAGS= -O3 -cpp -i_dynamic -mcmodel=large
#FLAGS= -O3 -cpp  -mcmodel=large -vec-report0 -traceback -check
#FLAGS= -O3 -cpp -vec-report0
FLAGS= -O3 -cpp -g -traceback -check bounds -warn all

### gfortran ###
# DEBUG
#FLAGS = -cpp -O3 -fdefault-double-8  
#FLAGS = -cpp -ggdb -O0 -fdefault-double-8 -pedantic -Wall -Waliasing \
        -Wsurprising -Wconversion -Wunderflow -ffpe-trap=invalid,zero,overflow \
        -fbounds-check -fimplicit-none -fstack-protector-all

#---------------------------------------------------
#   compilation time parameters (Y=on, N=off)
#   carefull all of the following is case sensitive
#---------------------------------------------------
#   MPI paralalelization (Y/N)
MPIP = Y
#   Double Precision (Y/N Default is single)
DOUBLEP= Y

#   Enable Passive MHD (Y/N)
#   Includes induction eq with no back reaction to the flow
#   This is compatible with HD solvers (i.e. HLL, HLLC)   
PMHD = N

#   Enable MHD (Y/N)
#   This is compatible with MHD solvers (HLLE, HLLD, under construction)
MHD = Y

#   Solver: 
#   HD  Solvers: HLL  (too difussive), or HLLC
#   MHD Solvers: HLLE (too difussive), or HLLD (under construction)
SOLVER = HLLD

#   Type of output
OUTBIN = Y
OUTVTK = N
OUTSILO= N

#   IF silo was selected make sure the libraries are in place,
#   adjust the following line for that purpose
#FLAGS += -I/usr/local/silo/include -I/usr/local/szip/include -I/usr/local/hdf5/include/
#LINKFLAGS = -L/usr/local/szip/lib -L/usr/local/hdf5/lib/ -L/usr/local/silo/lib/
#LINKFLAGS += -lsiloh5  -lhdf5_fortran -lsz -lz -lstdc++ 

#   additional equations (i.e. passive scalars) Y/N?
PASSIVES = N

#   Equation of state used to compute T and P
#   ADIABATIC     : Does not modify P, and T=(P/rho)*Tempsc
#   SINGLE_SPECIE : Uses only n (e.g. to use with tabulated cooling curves)
#   H_RATE        : Using n_HI and n_HII
#   CHEM          : Enables a full chemical network
EOS = ADIABATIC

#   Type of cooling (choose only one)
#   NONE: Turns off the cooling
#   H   : Single parametrized cooling function (ionization fraction and T)
#   BBC : Cooling function of Benjamin, Benson and Cox (2003)
#   DMC : coronal eq. (tabulated) from Dalgarno & Mc Cray (1972)
#   CHI : From table(s) generated with Chianti
#   CHEM : enables cooling from a full chemical network
COOLING = NONE

#   boundary conditions (OUTFLOW, CLOSED, PERIODIC)
#   choose only one per boundary
#   + OTHERB (user defined in user_mod.f90) if needed)
LEFTX   = PERIODIC
RIGHTX  = PERIODIC
BOTTOMY = PERIODIC
TOPY    = PERIODIC
INZ     = PERIODIC
OUTZ    = PERIODIC
OTHERB  = N

#   choice of slope limiter, available limiters are:
#   limiter =-1: no average, 0 : no limiter, 1: Minmod,
#   2: Van Leer, 3: Van Albada,4: UMIST, 5: Woodward
#   6: Superbee
LIMITER = 1

#   Thermal conduction?
#   options are:
#   OFF : No thermal conduction
#   ISO : Isotropic thermal conduction
#   ANISOTROPIC : Anisotropic thermal conduction (requires B field)
THERMAL_COND = N

#   Enable Saturation in the thermal conduction (Y/N)
TC_SATURATION = N

#   Enable diffuse radiation (Y/N) 
RADDIFF = N

# Implementing C2ray
C2ray = N

#   Include gravity (Y/N) (from point sources)   
GRAV = N

#   Include radiative pressure (Y/N)
RADPRES = N

#   Include terms proportional to DIV B (powell et al. 1999) (Y/N)
EIGHT_WAVE = Y

#   Flux Interpolated Central Diff. scheme (Toth 2000) (Y/N)
CT = N

#
CEXCHANGE = N

# Activate the write of the DivB in the output
DIVB = Y

#   Enable the N-Body module (Y/N) not implemented yet!
#   NBODY = N

#####################################################
# There should be no need to modify below this line #
#####################################################

OBJECTS = \
constants.o  \
parameters.o \
globals.o	\
orzag_tang.o \
user_mod.o \
hydro_core.o \
difrad.o \
cooling_h.o	\
cooling_dmc.o	\
cooling_chi.o	\
thermal_cond.o \
init.o \
Out_Silo_Module.o \
Out_VTK_Module.o \
Out_BIN_Module.o \
output.o \
boundaries.o	\
hll.o 	\
hllc.o 	\
hlle.o 	\
hlld.o 	\
sources.o \
hydro_solver.o \
main.o 

#  for projection alog a LOS
OBJECTSCOLDENS = \
constants.o \
parameters.o 	\
globals.o	\
exoplanet.o \
user_mod.o \
hydro_core.o \
difrad.o \
init.o \
coldens.o

#  For the Lyman Alpha Tau calculation
OBJECTSLYAT = \
constants.o \
parameters.o 	\
globals.o	\
exoplanet.o \
hydro_core.o \
lyman_alpha_tau.o\

#---------------------------------------------------
# Set flags
ifeq ($(DOUBLEP),Y) 
FLAGS += -DDOUBLEP
ifeq ($(COMPILER),ifort)
FLAGS += -r8
endif
ifeq ($(COMPILER),gfortran)
FLAGS += -fdefault-real-8
endif
endif
ifeq ($(SOLVER),HLL)
FLAGS += -DHLL
endif
ifeq ($(SOLVER),HLLC)
FLAGS += -DHLLC
endif
ifeq ($(SOLVER),HLLE)
FLAGS += -DHLLE
endif
ifeq ($(SOLVER),HLLD)
FLAGS += -DHLLD
endif
ifeq ($(OUTBIN),Y)
FLAGS += -DOUTBIN
endif
ifeq ($(OUTVTK),Y)
FLAGS += -DOUTVTK
endif
ifeq ($(OUTSILO),Y)
FLAGS += -DOUTSILO
endif
ifeq ($(PASSIVES),Y)
FLAGS += -DPASSIVES
endif
ifeq ($(EOS),ADIABATIC)
FLAGS += -DEOS_ADIABATIC
endif
ifeq ($(EOS),SINGLE_SPECIE)
FLAGS += -DEOS_SINGLE_SPECIE
endif
ifeq ($(EOS),H_RATE)
FLAGS += -DEOS_H_RATE
endif
ifeq ($(EOS),CHEM)
FLAGS += -DEOS_CHEM
endif
ifeq ($(COOLING),NONE)
FLAGS += -DNO_COOL
endif
ifeq ($(COOLING),H)
FLAGS += -DCOOLINGH
endif
ifeq ($(COOLING),DMC)
FLAGS += -DCOOLINGDMC
endif
ifeq ($(COOLING),CHI)
FLAGS += -DCOOLINGCHI
endif
ifeq ($(COOLING),BBC)
FLAGS += -DCOOLINGBBC
endif
ifeq ($(COOLING),CHEM)
FLAGS += -DCOOLINGCHEM
endif
ifeq ($(LEFTX),PERIODIC)
FLAGS += -DPERIODX
endif
ifeq ($(BOTTOMY),PERIODIC)
FLAGS += -DPERIODY
endif
ifeq ($(INZ),PERIODIC)
FLAGS += -DPERIODZ
endif
ifeq ($(LEFTX),CLOSED)
FLAGS += -DREFXL
endif
ifeq ($(RIGHTX),CLOSED)
FLAGS += -DREFXR
endif
ifeq ($(BOTTOMY),CLOSED)
FLAGS += -DREFYB
endif
ifeq ($(TOPY),CLOSED)
FLAGS += -DREFYT
endif
ifeq ($(INZ),CLOSED)
FLAGS += -DREFZI
endif
ifeq ($(OUTZ),CLOSED)
FLAGS += -DREFZO
endif
ifeq ($(LEFTX),OUTFLOW)
FLAGS += -DOUTFXL
endif
ifeq ($(RIGHTX),OUTFLOW)
FLAGS += -DOUTFXR
endif
ifeq ($(TOPY),OUTFLOW)
FLAGS += -DOUTFYT
endif
ifeq ($(BOTTOMY),OUTFLOW)
FLAGS += -DOUTFYB
endif
ifeq ($(INZ),OUTFLOW)
FLAGS += -DOUTFZI
endif
ifeq ($(OUTZ),OUTFLOW)
FLAGS += -DOUTFZO
endif
ifeq ($(OTHERB),Y)
FLAGS += -DOTHERB
endif
FLAGS += -DLIMITER=$(LIMITER)
ifeq ($(MPIP),Y)
FLAGS += -DMPIP
COMPILER = mpif90
endif
ifeq ($(RADDIFF),Y)
FLAGS += -DRADDIFF
endif
ifeq ($(RADDIFF_GLOBAL),Y)
FLAGS += -DRADDIFF_GLOBAL
endif
ifeq ($(NBODY),Y)
FLAGS += -DNBODY
endif
ifeq ($(THERMAL_COND),ISO)
FLAGS += -DTHERMAL_COND
FLAGS += -DTC_ISOTROPIC
endif
ifeq ($(THERMAL_COND),ANISOTROPIC)
FLAGS += -DTHERMAL_COND
FLAGS += -DTC_ANISOTROPIC
endif
ifeq ($(TC_SATURATION),Y)
FLAGS += -DSAT_COND
endif
ifeq ($(GRAV),Y)
FLAGS += -DGRAV
endif
ifeq ($(RADPRES),Y)
FLAGS += -DRADPRES
endif
ifeq ($(PMHD),Y)
FLAGS += -DPMHD
endif
ifeq ($(MHD),Y)
FLAGS += -DMHD
endif
ifeq ($(EIGHT_WAVE),Y)
FLAGS += -DEIGHT_WAVE
endif
ifeq ($(CT),Y)
FLAGS += -DCT
endif
ifeq ($(CEXCHANGE),Y)
FLAGS += -DCEXCHANGE
endif
ifeq ($(DIVB),Y)
FLAGS += -DDIVB
endif
ifeq ($(C2ray),Y)
FLAGS += -DC2ray
endif
#---------------------------------------------------
# Compilation rules
$(PROGRAM)  : prebuild ${OBJECTS}
	@echo Linking object files ...
	@echo Linking flags used: $(LINKFLAGS)
	@$(COMPILER) $(FLAGS) $(LINKFLAGS)  $(OBJECTS) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
  ${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
  ${SRC_PATH}/C2Ray/*.mod
	@echo "Done! (`date`)"	

coldens : prebuild ${OBJECTSCOLDENS}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTSCOLDENS) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
	${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
	${SRC_PATH}/C2Ray/*.mod
	@echo "Done! (`date`)"	

lyman_alpha_tau : prebuild ${OBJECTSLYAT}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTSLYAT) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
	${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
	${SRC_PATH}/C2Ray/*.mod
	@echo "Done! (`date`)"	

prebuild :
	@echo "Guacho build started `date`"
	@echo Compiling flags used: $(FLAGS)
%.o:%.f95
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:%.f90
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:%.f
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:%.F95
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:%.F90
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:%.F
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@

# C2ray routines and modules
%.o:src/C2ray/%.f95
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:src/C2Ray/%.f90
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:.src/C2Ray/%.f
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:src/C2Ray/%.F95
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:/src/C2Ray/%.F90
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@
%.o:src//C2Ray/%.F
	@echo Compiling $^ ...
	@$(COMPILER) $(FLAGS) -c $< -o $@

clean :
	rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
	  ${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
		${SRC_PATH}/C2Ray/*.mod
	rm -rf *genmod.f90 ${SRC_PATH}/*genmod.f90
	rm -f $(PROGRAM) lyman_alpha_tau *.out *.err 

cleanall :
	rm -f $(PROGRAM) lyman_alpha_tau *.out *.err
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
		${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
		${SRC_PATH}/C2Ray/*.modrm -f SS/BIN/*
	rm -f ${USR_PATH}/SILO/*.root
	rm -f ${USR_PATH}/SILO/BLOCKS/*
