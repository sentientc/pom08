#
# sbPOM makefile
#

#-----------------------------------------------------------------------
# Settings that depend on the system and the compiler
#-----------------------------------------------------------------------
# Set macros
CPP = cpp -P
FC = ifort 
LD = ifort
CLEAN = rm
# Set libraries and include files
NETCDFINC = -I/aracbox/lib/netcdf/4.1.2/intel_12/include
NETCDFLIB = -L/aracbox/lib/netcdf/4.1.2/intel_12/lib

#FFLAGS = -O3 -mcmodel large -shared-intel -fp-model precise -assume byterecl $(NETCDFINC)
FFLAGS = -check bounds -warn interface -g -traceback -mcmodel large -shared-intel -fp-model precise -assume byterecl $(NETCDFINC)

#LIBS = -O3 -mcmodel large -shared-intel -fp-model precise -assume byterecl $(NETCDFLIB) -lnetcdf -lnetcdff -limf -lm 
LIBS = -check bounds -warn interface -g -traceback -mcmodel large -shared-intel -fp-model precise -assume byterecl $(NETCDFLIB) -lnetcdf -lnetcdff -limf -lm 
#-----------------------------------------------------------------------
# Set the executable
#-----------------------------------------------------------------------
BIN = pom08.exe  #yoyo


#-----------------------------------------------------------------------
# Define source directory
#-----------------------------------------------------------------------
SRCDIR = pom

#-----------------------------------------------------------------------
# Define objects
#-----------------------------------------------------------------------
OBJS = pom08_iosub.o    \
       module_time.o    \
       interp.o         \
       gridgen.o        \
       wave.o           \
       wind.o           \
       pom08.o
#       interp.o         
VPATH = $(SRCDIR)

#-----------------------------------------------------------------------
# Set implicit rules for compilation
#-----------------------------------------------------------------------
%.o: %.f 
	@echo
	$(FC) -c $(FFLAGS)  $<

#-----------------------------------------------------------------------
# Set implicit rules for dependencies
#-----------------------------------------------------------------------
%.f: %.F
	@echo
	$(CPP) $(FFLAGS) $< > $*.f

#%.o:%.f90 
#	@echo
#	$(FC) -c $(FFLAGS)  $<
#-----------------------------------------------------------------------
# Create the executable
#-----------------------------------------------------------------------
$(BIN): $(OBJS)
	@echo
	$(LD) $(FFLAGS) -o $(BIN) $(OBJS) $(LIBS)

#-----------------------------------------------------------------------
# Cleaning target
#-----------------------------------------------------------------------
clean:
	@rm -f *.o *.mod *.il *.exe  *genmod.f90
