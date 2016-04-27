#---------------------------------------------------
#   for Guacho 3D
#---------------------------------------------------
#   Name of the executable
PROGRAM= grav

#  Seaerh Path, object files are searched in the order below
VPATH = . : ../src

#   Choice of compiler: ifort and gfortran are the only implemented
#   if MPI enabled, do not use mpif90, but rather the base compiler
#   of the mpi build (e.g. 'mpif90 --version')
COMPILER= gfortran
#COMPILER= ifort

#   Compiler flags, make sure they are compatible with the previous selection 
### ifortran ###
#FLAGS= -O3 -cpp -i_dynamic -mcmodel=large
#FLAGS= -O3 -cpp  -mcmodel=large -vec-report0 -traceback -check
#FLAGS= -O3 -cpp -vec-report0
#FLAGS= -O3 -cpp -g -traceback -check bounds -warn all

### gfortran ###
# DEBUG
FLAGS = -cpp -O3 
#FLAGS = -cpp -ggdb -O0 -fdefault-double-8 -pedantic -Wall -Waliasing \
#        -Wsurprising -Wconversion -Wunderflow -ffpe-trap=invalid,zero,overflow \
#        -fbounds-check -fimplicit-none -fstack-protector-all

#---------------------------------------------------
#   compilation time parameters (Y=on, N=off)
#   carefull all of the following is case sensitive
#---------------------------------------------------
#   MPI paralalelization (Y/N)
MPIP = Y
#   Double Precision (Y/N Default is single)
DOUBLEP= Y

#   Enable  Silo Output (requires external libraries)
SILO = N
#   IF silo was selected make sure the libraries are in place,
#   adjust the following line for that purpose
#FLAGS += -I/usr/local/silo/include -I/usr/local/szip/include -I/usr/local/hdf5/include/
#LINKFLAGS = -L/usr/local/szip/lib -L/usr/local/hdf5/lib/ -L/usr/local/silo/lib/
#LINKFLAGS += -lsiloh5  -lhdf5_fortran -lsz -lz -lstdc++ 

#####################################################
# There should be no need to modify below this line #
#####################################################

OBJECTS = \
constants.o  \
parameters.o \
globals.o	\
exoplanet.o \
jet.o \
user_mod.o \
network.o \
field_cd_module.o \
hydro_core.o \
linear_system.o \
chemistry.o \
difrad.o \
cooling_h.o \
cooling_dmc.o	\
cooling_chi.o \
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
sources.o \
hydro_solver.o \
main.o \

#  for projection alog a LOS
OBJECTSCOLDENS = \
constants.o \
parameters.o 	\
globals.o	\
jet.o \
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

ifeq ($(MPIP),Y)
FLAGS += -DMPIP
COMPILER = mpif90
endif

ifeq ($(SILO),Y)
FLAGS += -DOUT_SILO
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
