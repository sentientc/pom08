# README #

* Configuration and Dependencies  
	1. compiler  
		+ make sure there is ifort install  
		+ enable compiler option in run.sh   
	2. netcdf library  
		+ make sure there is netcdf library installed in the system  
		+ modify the include path for netcdf include accordingly. usually in /usr/include/netcdf.inc  
                + netcdf library is best to be compiled in the fortran used  
* Run the program with run.sh
	1. run the model with run script(execute run.sh)
	2. after execution, 3 files are generated. pom08(the executable), out/pom08.nc(netcdf output), and out/pom08.log(model output messages)
	3. results are in the out folder pom08.nc and pom08.log

* Dependencies  
	+ netcdf library version at least 4.1.2  

# History #
* hosting source code for pom08 per request from Leo Oey 
* the followings are concatenated from pomchanges and moreREADME from http://shoni2.princeton.edu/ftp/glm/
* see the pom08.f source for more comments
As of Feb. 28, 2008, I have made changes based on symetry
analysis. Thus, for a straight channel with jm=11, properties
should be symmetrical about j=6. They are now to 3 or 4
significant places if cor(i,j)=0. 


For mode=3, vertical u(i,j,k) profiles are noisy, but not
so for cor(i,j)=1.e-4 or when aam(i,j,k) (a la Smagorinsky)
is replaced by aam(i,j,k)=aam_init. 

 
In deep regions where q2's are very small, length scales are
meaningless.

The directory, PROBLEMS, contain various problem definitions
so that if one wishes to purge pom08.f of definitons in favor
of one's own problem, they can  be stored in the directory.

 Herewith is a list of changes to pom08.f or pom08.c

2/4/08
  In the initialization of s.r. wave, change tht_wnd=0 to tht_wnd(i,j)=0.

  Change data statement cdmax/.0025 to cdmax/.005 The former is too
      restrctive; e.g., for laboratory scale waves.
 
  In s.r. specaves add zeta,zetf to save list; e.g.,
      save init,kpDd,F1d,F2d,F3d,F4d,zeta,zetf
2/13/08
  In the do 8000 loop, -tpx0(i,j) and -tpy0(i,j) has been included in with
(wusurf(i,j)-wubot(i,j)) and (wvsurf(i,j)-wubot(i,j)) respectively.
2/18/08
  In s.r. advave, + Sxx, +Sxy, +Syy  should be removed.
2/28/08
  Changes in s.r.advct, s.r.advave, profu, profv.profq  and s.r. wave and 
its supporting s.r.'s. Replace these s.r.'s OR the whole pom08.f 
3/7/08
  In s.r. seamount, there is an obvious syntax error in defining ele(j)
and elw(j).
4/28/08
  Some changes in s.r. profq. Best to swap the s.r.
1/23/10
  Some changes in s.r. cur2wave and s.r. advct. Best to swap the s.r.'s
1/1/11
  In s.r. wave, change save init to save dtw2,dtinv
  in s.r. advct, call wave2cur is redundant and has been commented out.
2/25/12
  Stokes drift has been removed from velocities in profu, profv and profq.
  Stokes drift is governed by wave properties and should not be mixed
  in profu, profv or participate in TKE production.
3/8/12
  A number of changes in profq particularly in calculating prod. Best to
swap old profq for the new.
8/16/12
  specavs has been chnged so that that which was a surface delta function
is now distributed. S.R. wave2cur has been changed slightly; there are
2 sign changes and the additional delta function has been removed. It's
distributed version is now incorporated in Sxx and Syy.

Feb. 12, 2014, changes were made in recognition that the
term that was a delta function in the stress radiation is
now distributed across the water column. See papers in the
previous directory. Thus the declaration of Srad0 in S.R.
cur2wave should be deleted  and the minus sign, just above,
after Srad(i,j)= should also be deleted. Srad0 can be
deleted (or not) everywhere else in the code.

