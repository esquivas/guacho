#---------------------------------------------------
#   for Guacho 3D
#---------------------------------------------------
#   Name of the executable
PROGRAM=guacho

#  Seaerh Path, object files are searched in the order below
VPATH = . : ../src

#   Choice of compiler: ifort and gfortran are the only implemented
#   if MPI enabled, do not use mpif90, but rather the base compiler
#   of the mpi build (e.g. 'mpif90 --version')
COMPILER= ifort
#COMPILER= gfortran

#   Compiler flags, make sure they are compatible with the previous selection
### ifortran ###
 FLAGS= -O3 -cpp -mcmodel=large -g -traceback -check bounds -warn all
#FLAGS= -O3 -cpp -mcmodel=large -traceback -check
#FLAGS= -O3 -cpp
#FLAGS= -O3 -cpp -g -traceback -check bounds -warn all

### gfortran ###
#FLAGS = -cpp -O3
# For DEBUG, suggestions:
#FLAGS = -cpp -Og -ggdb -pedantic -Wall -Waliasing -fbacktrace -Wsurprising  \
        -Wconversion -Wunderflow -ffpe-trap=invalid,zero,overflow           \
				-fbounds-check -fimplicit-none -fstack-protector-all

#---------------------------------------------------
#   compilation time parameters (Y=on, N=off)
#   carefull all of the following is case sensitive
#---------------------------------------------------
#   MPI paralalelization (Y/N)
MPIP = Y

#   Include Magnetic field (Y/N) , (either full MHD or passive B)
BFIELD = Y

#   Include passive scalars (Y/N)
PASSIVES = N

#   Double Precision (Y/N Default is single)
DOUBLEP= Y

#   Enable  Silo Output (requires external libraries)
SILO = N
#   IF silo was selected make sure the libraries are in place,
#   adjust the following line for that purpose
#FLAGS += -I/usr/local/silo/include -I/usr/local/szip/include -I/usr/local/hdf5/include/
#LINKFLAGS = -L/usr/local/szip/lib -L/usr/local/hdf5/lib/ -L/usr/local/silo/lib/
#LINKFLAGS += -lsiloh5  -lhdf5_fortran -lsz -lz -lstdc++

# Additional USER modules to compile (they can only call constants,
# parameters and globals)
MODULES_USER= \
snr.o
#nbody_mod2.o \
#starwind.o

#####################################################
# There should be no need to modify below this line #
#####################################################

MODULES_MAIN = \
constants.o  \
parameters.o \
globals.o \
utils.o

OBJECTS_MAIN=\
lmp_module.o \
user_mod.o \
self_gravity.o \
network.o \
flux_cd_module.o \
hydro_core.o \
linear_system.o \
chemistry.o \
difrad.o \
cooling_h.o \
cooling_dmc.o	\
cooling_chi.o \
cooling_schure.o \
thermal_cond.o \
init.o \
Out_Silo_Module.o \
Out_VTK_Module.o \
Out_BIN_Module.o \
output.o \
boundaries.o  \
hll.o 	\
hllc.o 	\
hlle.o 	\
hlld.o 	\
hlle_split_all.o \
sources.o \
hydro_solver.o \
main.o

#  for projection alog a LOS
COLDENS = \
user_mod.o \
network.o \
hydro_core.o \
difrad.o \
coldens.o

#  For the Lyman Alpha Tau calculation
LYAT = \
hydro_core.o \
lyman_alpha_tau.o  #   not updated to Guacho V1.3

#  to get stokes parameters with the Lagrangian particles
STOKES_LP = \
lmp_module.o \
hydro_core.o \
stokes_lp.o

# Join object lists
OBJECTS_ALL       = ${MODULES_MAIN} ${MODULES_USER} ${OBJECTS_MAIN}
OBJECTS_NB        = ${MODULES_MAIN} ${NBODY}
OBJECTS_COLDENS   = ${MODULES_MAIN} ${MODULES_USER} ${COLDENS}
OBJECTS_LYAT      = ${MODULES_MAIN} ${MODULES_USER} ${LYAT}
OBJECTS_STOKES_LP = ${MODULES_MAIN} ${MODULES_USER} ${STOKES_LP}
#---------------------------------------------------
# Set flags
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

ifeq ($(MPIP),Y)
FLAGS += -DMPIP
COMPILER = mpif90
endif

ifeq ($(PASSIVES),Y)
FLAGS += -DPASSIVES
endif

ifeq ($(BFIELD),Y)
FLAGS += -DBFIELD
endif

ifeq ($(SILO),Y)
FLAGS += -DOUT_SILO
endif
#---------------------------------------------------
# Compilation rules
$(PROGRAM)  : prebuild ${OBJECTS_ALL}
	@echo Linking object files ...
	@echo Linking flags used: $(LINKFLAGS)
	@$(COMPILER) $(FLAGS) $(LINKFLAGS)  $(OBJECTS_ALL) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
  ${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
  ${SRC_PATH}/C2Ray/*.mod
	@echo "Done! (`date`)"

test_nbody : prebuild ${OBJECTS_NB}
		@echo Linking object files ...
		@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTS_NB) -o $@
		@echo Cleaning up ...
		@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
		${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
		${SRC_PATH}/C2Ray/*.mod
		@echo "Done! (`date`)"

coldens : prebuild ${OBJECTS_COLDENS}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTS_COLDENS) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
	${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
	${SRC_PATH}/C2Ray/*.mod
	@echo "Done! (`date`)"

lyman_alpha_tau : prebuild ${OBJECTS_LYAT}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTS_LYAT) -o $@
	@echo Cleaning up ...
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
	${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
	${SRC_PATH}/C2Ray/*.mod
	@echo "Done! (`date`)"

stokes_lp : prebuild ${OBJECTS_STOKES_LP}
	@echo Linking object files ...
	@$(COMPILER) $(FLAGS) $(LINKFLAGS) $(OBJECTS_STOKES_LP) -o $@
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
	rm -f $(PROGRAM) lyman_alpha_tau coldens stokes_lp test_nbody *.out *.err

cleanall :
	rm -f $(PROGRAM) lyman_alpha_tau coldens stokes_lp test_nbody *.out *.err
	@rm -f *.o *.mod ${SRC_PATH}/*.o  ${SRC_PATH}/*.mod \
		${USR_PATH}/*.o ${USR_PATH}/*.mod ${SRC_PATH}/C2Ray/*.o  \
		${SRC_PATH}/C2Ray/*.modrm -f SS/BIN/*
	rm -f ${USR_PATH}/SILO/*.root
	rm -f ${USR_PATH}/SILO/BLOCKS/*
