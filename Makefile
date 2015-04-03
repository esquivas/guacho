#---------------------------------------------------
#   for Guacho 3D
#---------------------------------------------------
#   Name of the executable
PROGRAM=Guacho

#   Choice of compiler: ifort and gfortran are tested
#   if MPI enabled, do not use mpif90, but rather the base compiler
#   of the mpi build (see 'mpif90 --version')
#COMPILER= gfortran
COMPILER= ifort

#   Compiler flags, make sure they are compatible with 
#   the previous selection
#FLAGS= -O3 -cpp -i_dynamic -mcmodel=large
#FLAGS= -O3 -cpp  -mcmodel=large -vec-report0 -traceback -check
#FLAGS= -O3 -cpp -vec-report0
FLAGS= -O3 -cpp 
# -g -traceback -check bounds -warn all

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
MHD = N

#   Solver: 
#   HD  Solvers: HLL  (too difussive), or HLLC
#   MHD Solvers: HLLE (too difussive), or HLLD (under construction)
SOLVER = HLLC

#   Type of output
OUTDAT = N
OUTBIN = Y
OUTVTK = N
OUTSILO= N

#   IF silo was selected make sure the libraries are in place,
#   adjust the following line for that purpose
#FLAGS += -I/usr/local/silo/include -I/usr/local/szip/include -I/usr/local/hdf5/include/
#LINKFLAGS = -L/usr/local/szip/lib -L/usr/local/hdf5/lib/ -L/usr/local/silo/lib/
#LINKFLAGS += -lsiloh5  -lhdf5_fortran -lsz -lz -lstdc++ 

#   additional equations (i.e. passive scalars)?
PASSIVES = Y

#   Equation of state used to compute T and P
#   ADIABATIC     : Does not modify P, and T=(P/rho)*Tempsc
#   SINGLE_SPECIE : Uses only n (e.g. to use with tabulated cooling curves)
#   H_RATE        : Using n_HI and n_HII
#   MULTI_SPECIES : not implemented yet
EOS = H_RATE

#   Type of cooling (choose only one)
#   NONE: Turns off the cooling
#   H   : Single parametrized cooling function (ionization fraction and T)
#   BBC : Cooling function of Benjamin, Benson and Cox (2003)
#   DMC : coronal eq. (tabulated) from Dalgarno & Mc Cray (1972)
#   CHI : From table(s) generated with Chianti
COOLING = H

#   boundary conditions (OUTFLOW, CLOSED, PERIODIC)
#   choose only one per boundary
#   + OTHERB (user defined in user_mod.f90) if needed)
LEFTX   = OUTFLOW
RIGHTX  = OUTFLOW
BOTTOMY = OUTFLOW
TOPY    = OUTFLOW
INZ     = OUTFLOW
OUTZ    = OUTFLOW
OTHERB  = Y

#   choice of slope limiter, available limiters are:
#   limiter =-1: no average, 0 : no limiter, 1: Minmod,
#   2: Van Leer, 3: Van Albada,4: UMIST, 5: Woodward
#   6: Superbee
LIMITER = 1

#   Enable (isotropic) thermal conduction (Y/N)
THERMAL_COND = N

#   Enable diffuse radiation (Y/N) 
RADDIFF = Y

# Implementing C2ray
C2ray = N

#   Include gravity (Y/N) (from point sources)   
GRAV = Y

#   Include radiative pressure (y/N)
RADPRES = N

#   Include terms proportional to DIV B (powell et al. 1999) (y/N)
EIGHT_WAVE = N
#
CEXCHANGE = N

#   Enable the N-Body module (Y/N) not implemented yet!
#   NBODY = N

#####################################################
# There should be no need to modify below this line #
#####################################################

OBJECTS = \
./src/constants.o  \
./src/parameters.o \
./src/globals.o	\
./src/exoplanet.o \
./src/hydro_core.o \
./src/difrad.o \
./src/cooling_h.o	\
./src/cooling_dmc.o	\
./src/cooling_chi.o	\
./src/thermal_cond.o \
./src/user_mod.o \
./src/init.o \
./src/Out_Silo_Module.o \
./src/output.o \
./src/boundaries.o	\
./src/hll.o 	\
./src/hllc.o 	\
./src/hlle.o 	\
./src/hlld.o 	\
./src/sources.o \
./src/hydro_solver.o \
./src/main.o 
#./src/cexchange.o \

#  for projection alog a LOS
OBJECTSCOLDENS = \
./src/constants.o \
./src/parameters.o 	\
./src/globals.o	\
./src/exoplanet.o \
./src/hydro_core.o \
./src/difrad.o \
./src/user_mod.o \
./src/init.o \
./src/coldens.o

#  For the Lyman Alpha Tau calculation
OBJECTSLYAT = \
./src/constants.o \
./src/parameters.o 	\
./src/globals.o	\
./src/hydro_core.o \
./src/lyman_alpha_tau.o
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
ifeq ($(OUTDAT),Y)
FLAGS += -DOUTDAT
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
ifeq ($(EOS),MULTI_SPECIES)
FLAGS += -DEOS_MULTI_SPECIES
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
ifeq ($(THERMAL_COND),Y)
FLAGS += -DTHERMAL_COND
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
ifeq ($(CEXCHANGE),Y)
FLAGS += -DCEXCHANGE
endif
ifeq ($(C2ray),Y)
FLAGS += -DC2ray
endif
#---------------------------------------------------
# Compilation rules
$(PROGRAM)  : prebuild ${OBJECTS}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS)  $(OBJECTS) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod src/*.o src/*.mod src/C2Ray/*.o src/C2Ray/*.mod
	@echo "Done! (`date`)"	

coldens : prebuild ${OBJECTSCOLDENS}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTSCOLDENS) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod src/*.o src/*.mod src/C2Ray/*.o src/C2Ray/*.mod
	@echo "Done! (`date`)"	

lyman_alpha_tau : prebuild ${OBJECTSLYAT}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTSLYAT) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod src/*.o src/*.mod src/C2Ray/*.o src/C2Ray/*.mod
	@echo "Done! (`date`)"	

prebuild :
	@echo "Guacho build started `date`"

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
	rm -f *.o *.mod src/*.o src/*.mod src/C2ray/*.o scr/C2Ray/*.mod 
	rm -rf *genmod.f90 src/*genmod.f90
	rm -f $(PROGRAM) lyman_alpha_tau *.out *.err 

cleanall :
	rm -f $(PROGRAM) lyman_alpha_tau *.out *.err
	rm -f *.o *.mod src/*.o src/*.mod src/*/*.o src/*/*.mod
	rm -f SS/BIN/*
	rm -f SILO/*.root
	rm -f SILO/BLOCKS/*
