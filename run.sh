#!/bin/bash
rundir=`pwd`
cd $rundir
mkdir out
### original compile option
#g77 -o pom2k ./pom08.f
### add -lnetcdff to include netcdf library
### should be usable in most linux with gfortran or g77 install
#gfortran -o pom08 ./pom08.f -lnetcdff -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal
#gfortran -o pom08 ./pom08.f -lnetcdff
####enabling netcdf.1. using ifort as the library is compiled with ifort. 2. sequence of include , soource and library matters and this is the working sequence
### used in atop cluster
#-fpe0 so when n/0 and other situation occur the program will stop, it is easier to debug with -traceback option
ifort -o pom08.x -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
#ifort -fpe0 -o pom08 -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
#ifort -check bounds -align all -o pom2k -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08_WithnetcdfOutput.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
#ifort -check bounds -warn interface -g -traceback -align all -o pom2k -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08_WithnetcdfOutput.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
./pom08.x > out/pom08.log 
