#!/bin/bash -xe
#rundir=`pwd`
scnm='e1'
wrkdir='/wrk/simon/fortran/pom2k/latest/pom08'
rundir='/archive7/simon/pom_wave/pac11/'$scnm
gridf='/archive5/simon/gridf/pac11/out.30km/grid'
windf='/archive7/simon/mpiPOM/global/exp001'
\rm -rf $rundir
mkdir -p $rundir
cd $rundir
mkdir out
mkdir pom 
\cp -f $wrkdir/pom/* pom/
\cp -f $wrkdir/makefile .
ln -sf $wrkdir/specavs .
\cp -f $wrkdir/pom08.c .
echo '-----------Running '$scnm'-------------'
sed -i s/'      time_start = '.*/'      time_start = '\''2001-01-01 00:00:00'\'/ pom/pom08.f
sed -i s/'      days = '.*/'      days = '2/ pom/pom08.f
if [[ $scnm == *"e0" ]]; then
 sed -i s/'      iproblem = '.*/'      iproblem = '1/ pom/pom08.f
# sed -i s/'      mode='.*/'      mode=0'/ pom/pom08.f
 sed -i s/'    $ (im='.*/'    $ (im=24        ,jm=11            ,kb=21)'/ pom08.c
 sed -i s/'      days = '.*/'      days = '10/ pom/pom08.f
elif [[ $scnm == *e1 ]]; then
 mkdir in
 ln -sf $gridf in/
 ln -sf $windf in/wind
 sed -i s/'      time_start = '.*/'      time_start = '\''2002-07-10 00:00:00'\'/ pom/pom08.f
 sed -i s/'      days = '.*/'      days = 35'/ pom/pom08.f
 sed -i s/'      iproblem = '.*/'      iproblem = '6/ pom/pom08.f
 sed -i s/'      mode='.*/'      mode=0'/ pom/pom08.f
 sed -i s/'      dte=2'.*/'      dte=30'/ pom/pom08.f
 sed -i s/'      prtd2='.*/'      prtd2=125e-3'/ pom/pom08.f
 sed -i s/'    $ (im='.*/'    $ (im=650        ,jm=452            ,kb=21)'/ pom08.c
elif [[ $scnm == *e2 ]]; then
 mkdir in
 ln -sf $gridf in/
 ln -sf $windf in/wind
 sed -i s/'      time_start = '.*/'      time_start = '\''2002-07-10 00:00:00'\'/ pom/pom08.f
 sed -i s/'      days = '.*/'      days = 25'/ pom/pom08.f
 sed -i s/'      iproblem = '.*/'      iproblem = '6/ pom/pom08.f
 sed -i s/'      mode='.*/'      mode=0'/ pom/pom08.f
 sed -i s/'      dte=2'.*/'      dte=30'/ pom/pom08.f
 sed -i s/'      prtd2='.*/'      prtd2=1e0'/ pom/pom08.f
 sed -i s/'    $ (im='.*/'    $ (im=650        ,jm=452            ,kb=21)'/ pom08.c
elif [[ $scnm == *e3 ]]; then
 mkdir in
 ln -sf $gridf in/
 ln -sf $windf in/wind
 sed -i s/'      time_start = '.*/'      time_start = '\''2002-07-10 00:00:00'\'/ pom/pom08.f
 sed -i s/'      days = '.*/'      days = 25'/ pom/pom08.f
 sed -i s/'      iproblem = '.*/'      iproblem = '6/ pom/pom08.f
 sed -i s/'      mode='.*/'      mode=0'/ pom/pom08.f
 sed -i s/'      dte=2'.*/'      dte=15'/ pom/pom08.f
 sed -i s/'      prtd2='.*/'      prtd2=1e0'/ pom/pom08.f
 sed -i s/'    $ (im='.*/'    $ (im=650        ,jm=452            ,kb=21)'/ pom08.c
elif [[ $scnm == *e4 ]]; then
 mkdir in
 ln -sf $gridf in/
 ln -sf $windf in/wind
 sed -i s/'      iproblem = '.*/'      iproblem = '6/ pom/pom08.f
 sed -i s/'      mode='.*/'      mode=0'/ pom/pom08.f
# sed -i s/'      dte=2'.*/'      dte=0.05'/ pom/pom08.f
 sed -i s/'    $ (im='.*/'    $ (im=650        ,jm=452            ,kb=21)'/ pom08.c
fi
#export LD_LIBRARY_PATH=/aracbox/mpi/openmpi/1.5.4/intel_12/lib:/aracbox/lib/netcdf/4.1.2/intel_12/lib:/aracbox/intel/composerxe-2011.5.220/ipp/lib:/aracbox/intel/composerxe-2011.5.220/mkl/lib/intel64:/aracbox/intel/composerxe-2011.5.220/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21:/aracbox/intel/composerxe-2011.5.220/compiler/lib/intel64
export LD_LIBRARY_PATH=/aracbox/lib/netcdf/4.1.2/intel_12/lib:/aracbox/intel/composerxe-2011.5.220/ipp/lib:/aracbox/intel/composerxe-2011.5.220/mkl/lib/intel64:/aracbox/intel/composerxe-2011.5.220/tbb/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21:/aracbox/intel/composerxe-2011.5.220/compiler/lib/intel64
#/aracbox/Modules/3.2.9/bin/modulecmd bash load intel_poe
#module load intel_poe
make clean
make
./pom08.exe > out/pom08.log 
#C -- seamount      (iproblem=1)
#c     $ (im=24        ,jm=11            ,kb=21)
#c -- file2ic (iproblem=3)
#C    $ (im=41          ,jm=61          ,kb=16)
#C -- nybight  (iproblem=4)
#C    $ (im=40          ,jm=40          ,kb=21)
#c    $ (im=200         ,jm=200         ,kb=21)
#c    $ (im=401        ,jm=251          ,kb=25) 
#C -- read data  (iproblem=6)
#echo '      parameter (iproblem=1)' >> pom08.c
### original compile option
#g77 -o pom2k ./pom08.f
### add -lnetcdff to include netcdf library
### should be usable in most linux with gfortran or g77 install
#gfortran -o pom08 ./pom08.f -lnetcdff -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,precision,denormal
#gfortran -o pom08 ./pom08.f -lnetcdff
####enabling netcdf.1. using ifort as the library is compiled with ifort. 2. sequence of include , soource and library matters and this is the working sequence
### used in atop cluster
#-fpe0 so when n/0 and other situation occur the program will stop, it is easier to debug with -traceback option
#ifort -o pom08.x -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
#ifort -fpe0 -o pom08 -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
#ifort -check bounds -align all -o pom2k -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08_WithnetcdfOutput.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
#ifort -check bounds -warn interface -g -traceback -align all -o pom2k -I/aracbox/lib/netcdf/4.1.2/intel_12/include ./pom08_WithnetcdfOutput.f -L/aracbox/lib/netcdf/4.1.2/intel_12/lib -lnetcdf -lnetcdff -limf -lm
