C***********************************************************************
C
C     NetCDF subroutines for pom08.f
C
C     For version number and date, see start of subroutine write_netcdf
C
C     Add wusurf,wvsurf (May/2008; lyo@princeton.edu)!lyo:_20080507:
C     Add km,kh,q2,q2l (Apr/2008; zhibins@princeton.edu)
C     Add WAD variables (Aug/2007; lyo@princeton.edu)
C
C***********************************************************************
C
      subroutine def_var_netcdf(ncid,name,nvdims,vdims,varid,
     $                          n_long_name,long_name,n_units,units,
     $                          n_coords,coords,lcoords,
     $                          nf_float,nf_noerr)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Defines a netCDF variable and its attributes.       *
C *                                                                    *
C *                ncid ........ the netCDF I.D.                       *
C *                name ........ the variable name                     *
C *                nvdims ...... the number of dimensions of name      *
C *                vdims ....... a vector of length at least nvdims,   *
C *                              containing the dimension I.D.s        *
C *                varid ....... the variable I.D. returned            *
C *                n_long_name . the number of characters in long_name *
C *                long_name ... the long name for the variable        *
C *                n_units ..... the number of characters in units     *
C *                units ....... the units of the variable             *
C *                n_coords .... the number of characters in coords    *
C *                              (if applicable)                       *
C *                coords ...... the names of variables for the        *
C *                              "coordinates" attribute               *
C *                              (if applicable)                       *
C *                lcoords ..... .true. if a "coordinates" attribute   *
C *                              required, otherwise .false.           *
c *                nf_float .... an integer defining the netCDF        *
C *                              variable type                         *
C *                nf_noerr .... an integer defining the netCDF "no    *
C *                              error" status                         *
C *                                                                    *
C *                (nf_float and nf_noerr are declared and defined in  *
C *                 the file netcdf.inc "included" in the subroutine   *
C *                 write_netcdf.)                                     *
C *                                                                    *
C **********************************************************************
C
      integer vdims(4)
C
      integer ncid,nf_noerr,nvdims,varid,n_long_name,n_units,n_coords
C
      logical lcoords
C
      character*(*) name,long_name,units,coords
C
      integer status
C
      status=nf_def_var(ncid,name,nf_float,nvdims,vdims,varid)
C
      call handle_netcdf_error('nf_def_var          ',
     $                         status,nf_noerr)
C
C     write(6,1) name,varid
C   1 format('Variable ID returned by nf_def_var for variable ',
C    $       a7,'  = ',i5)
C
      status=nf_put_att_text(ncid,varid,'long_name',
     $                       n_long_name,long_name)
      call handle_netcdf_error('nf_put_att_text     ',
     $                         status,nf_noerr)
C
      status=nf_put_att_text(ncid,varid,'units',
     $                       n_units,units)
      call handle_netcdf_error('nf_put_att_text     ',
     $                         status,nf_noerr)
C
C     Add coordinates attribute, if necessary:
C
      if(lcoords) then
C
        status=nf_put_att_text(ncid,varid,'coordinates',
     $                         n_coords,coords)
        call handle_netcdf_error('nf_put_att_text     ',
     $                           status,nf_noerr)
C
      endif
C
      return
C
      end
C
      subroutine handle_netcdf_error(routine,status,nf_noerr)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Checks for netCDF error.                            *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      integer status,nf_noerr
C
      character*20 routine
C
      character*80 nf_strerror
C
      if (status.ne.nf_noerr) then
C
        write(6,1) routine,nf_strerror(status)
    1   format('NetCDF routine ',a20,' terminated with error:'/2x,a80)
        stop
C
      else
C
        return
C
      endif
C
      end
C
      subroutine write_netcdf(netcdf_file,option)
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Initialises, writes to and closed netCDF file.      *
C *                                                                    *
C *                netcdf_file ....... name of output file             *
C *                option ............ 1 to initialise file            *
C *                                    2 to write a set of data        *
C *                                    3 to close file                 *
C *                                                                    *
C *                In order to include another variable:               *
C *                                                                    *
C *                1. Within the loop "if(option.eq.1) then":          *
C *                                                                    *
C *                   Insert appropriate "call def_var_netcdf(.....)", *
C *                   followed by "status=nf_put_att_<type>(.....)"    *
C *                   statements for any additional attributes. Each   *
C *                   of these latter statements should be followed by *
C *                   a "call handle_netcdf_error(.....)" statement.   *
C *                                                                    *
C *                2. Within the loop "if(option.eq.2) then":          *
C *                                                                    *
C *                   Insert appropriate                               *
C *                   "status=nf_put_var_<type>(.....)" statement,     *
C *                   followed by a "call handle_netcdf_error(.....)"  *
C *                   statement.                                       *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom08.c'
C
C     Filename of netCDF "include" file:
C
C      include '/usr/local/include/netcdf.inc'  !GFDL:
c      include '/usr/include/netcdf.inc'        !Princeton U:!sc: usable
c      in most linux installation with gfortran or g77 installed
      include '/aracbox/lib/netcdf/4.1.2/intel_12/include/netcdf.inc'        !sc:netcdf.inc on ATOP cluster 
C
      integer option
C
      character*120 netcdf_file
C
      integer count(4),start(4),vdims(4)
C
      integer dum_varid,dvm_varid
      integer dx_varid,dy_varid
      integer east_c_varid,east_e_varid,east_u_varid,east_v_varid
      integer elb_varid
      integer fsm_varid
      integer h_varid
      integer iout
      integer ncid
      integer north_c_varid,north_e_varid,north_u_varid,north_v_varid
      integer rho_varid,rmean_varid,rot_varid
      integer status,s_varid
      integer time_dimid,time_varid,t_varid
      integer uab_varid,u_varid
      integer vab_varid,v_varid
      integer w_varid
      integer x_dimid
      integer y_dimid
      integer z_dimid,z_varid,zz_varid
      integer wetmask_varid,d_varid,elbmsl_varid !lyo:!wad:
c      integer i,j  !lyo:!wad:!sc:already defined in pom08.c
      integer km_varid,kh_varid,q2_varid,q2l_varid  !zb:EX1-1A
      integer wusurf_varid,wvsurf_varid !lyo:EX1-3A!lyo:_20080507:
C!sc:wave variables out:start
      integer cg_varid,thtav_varid,Hs_varid,ent_varid,sigp_varid
     & ,kthav_varid,tpx0_varid,tpy0_varid,ttx0_varid,tty0_varid
C!sc:wave variables out:end
C
      save dum_varid,dvm_varid
      save dx_varid,dy_varid
      save east_c_varid,east_e_varid,east_u_varid,east_v_varid
      save elb_varid
      save fsm_varid
      save h_varid
      save iout
      save ncid
      save north_c_varid,north_e_varid,north_u_varid,north_v_varid
      save rho_varid,rmean_varid,rot_varid
      save s_varid
      save time_dimid,time_varid,t_varid
      save uab_varid,u_varid
      save vab_varid,v_varid
      save w_varid
      save x_dimid
      save y_dimid
      save z_dimid,z_varid,zz_varid
      save wetmask_varid,d_varid,elbmsl_varid !lyo:!wad:
      save km_varid,kh_varid,q2_varid,q2l_varid  !zb:EX1-1A
      save wusurf_varid,wvsurf_varid  !lyo:EX1-3A!lyo:_20080507:
C!sc:wave variables out:start
      save cg_varid,thtav_varid,Hs_varid,ent_varid,sigp_varid
     & ,kthav_varid,tpx0_varid,tpy0_varid,ttx0_varid,tty0_varid
C!sc:wave variables out:end
C
C***********************************************************************
C
C     source_n should agree with source defined in pom08.f. 
C
      character*40 source_n
      parameter(source_n='pom08  2008-04-19')
C
C!sc:remove source check!no idea what this does
c      if(source.ne.source_n) then
c
c        write(6,4)
c    4   format(/'Incompatible versions of program and include files ',
c     $          '..... program terminated; in pom08.n'/)
c        stop
C
c      endif
C
C***********************************************************************
C
      if(option.eq.1) then
C
C     Initialise netCDF ouput:
C
C     Initialise time index:
C
      iout=0
C
C     Create netcdf file:
C
        status=nf_create(netcdf_file,nf_clobber,ncid)
        call handle_netcdf_error('nf_create           ',
     $                           status,nf_noerr)
C       write(6,1) ncid
C   1   format('ncid returned by nf_create = ',i5)
C
C     Define global attributes:
C
        status=nf_put_att_text(ncid,nf_global,'source',
     $                         40,source)
        call handle_netcdf_error('nf_put_att_text     ',
     $                           status,nf_noerr)
C
        status=nf_put_att_text(ncid,nf_global,'title',
     $                         40,title)
        call handle_netcdf_error('nf_put_att_text     ',
     $                           status,nf_noerr)
C
C     Define dimensions:
C
        status=nf_def_dim(ncid,'time',nf_unlimited,time_dimid)
        call handle_netcdf_error('nf_def_dim          ',
     $                           status,nf_noerr)
C       write(6,2) 'time',time_dimid
C   2   format('Dimension ID returned by nf_def_dim for variable ',
C    $         a4,' = ',i5)
C
        status=nf_def_dim(ncid,'z',kb,z_dimid)
        call handle_netcdf_error('nf_def_dim          ',
     $                           status,nf_noerr)
C       write(6,2) 'z',z_dimid
C
        status=nf_def_dim(ncid,'y',jm,y_dimid)
        call handle_netcdf_error('nf_def_dim          ',
     $                           status,nf_noerr)
C       write(6,2) 'y',y_dimid
C
        status=nf_def_dim(ncid,'x',im,x_dimid)
        call handle_netcdf_error('nf_def_dim          ',
     $                           status,nf_noerr)
C       write(6,2) 'x',x_dimid
C
C     Define variables and their attributes:
C
        vdims(1)=time_dimid
C
        call def_var_netcdf(ncid,'time',1,vdims,time_varid,
     $                      4,'time',37,'days since '//time_start,
     $                      1,' ',.false.,
     $                      nf_float,nf_noerr)
C
        vdims(1)=z_dimid
C
        call def_var_netcdf(ncid,'z',1,vdims,z_varid,
     $                      18,'sigma of cell face',11,'sigma_level',
     $                      1,' ',.false.,
     $                      nf_float,nf_noerr)
        status=nf_put_att_text(ncid,z_varid,'standard_name',
     $                         22,'ocean_sigma_coordinate')
        call handle_netcdf_error('nf_put_att_text     ',
     $                           status,nf_noerr)
        status=nf_put_att_text(ncid,z_varid,'formula_terms',
     $                         26,'sigma: z eta: elb depth: h')
        call handle_netcdf_error('nf_put_att_text     ',
     $                           status,nf_noerr)
C
        call def_var_netcdf(ncid,'zz',1,vdims,zz_varid,
     $                      20,'sigma of cell centre',11,'sigma_level',
     $                      1,' ',.false.,
     $                      nf_float,nf_noerr)
        status=nf_put_att_text(ncid,zz_varid,'standard_name',
     $                         22,'ocean_sigma_coordinate')
        call handle_netcdf_error('nf_put_att_text     ',
     $                           status,nf_noerr)
        status=nf_put_att_text(ncid,zz_varid,'formula_terms',
     $                         27,'sigma: zz eta: elb depth: h')
        call handle_netcdf_error('nf_put_att_text     ',
     $                           status,nf_noerr)
C
        vdims(1)=x_dimid
        vdims(2)=y_dimid
C
        call def_var_netcdf(ncid,'dx',2,vdims,dx_varid,
     $                      19,'grid increment in x',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'dy',2,vdims,dy_varid,
     $                      19,'grid increment in y',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'east_u',2,vdims,east_u_varid,
     $                      19,'easting of u-points',5,'metre',
     $                      14,'east_u north_u',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'east_v',2,vdims,east_v_varid,
     $                      19,'easting of v-points',5,'metre',
     $                      14,'east_v north_v',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'east_e',2,vdims,east_e_varid,
     $                      27,'easting of elevation points',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'east_c',2,vdims,east_c_varid,
     $                      23,'easting of cell corners',5,'metre',
     $                      14,'east_c north_c',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'north_u',2,vdims,north_u_varid,
     $                      20,'northing of u-points',5,'metre',
     $                      14,'east_u north_u',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'north_v',2,vdims,north_v_varid,
     $                      20,'northing of v-points',5,'metre',
     $                      14,'east_v north_v',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'north_e',2,vdims,north_e_varid,
     $                      28,'northing of elevation points',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'north_c',2,vdims,north_c_varid,
     $                      24,'northing of cell corners',5,'metre',
     $                      14,'east_c north_c',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'rot',2,vdims,rot_varid,
     $               34,'Rotation angle of x-axis wrt. east',6,'degree',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'h',2,vdims,h_varid,
     $                      23,'undisturbed water depth',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'fsm',2,vdims,fsm_varid,
     $                      17,'free surface mask',13,'dimensionless',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'dum',2,vdims,dum_varid,
     $                      15,'u-velocity mask',13,'dimensionless',
     $                      14,'east_u north_u',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'dvm',2,vdims,dvm_varid,
     $                      15,'v-velocity mask',13,'dimensionless',
     $                      14,'east_v north_v',.true.,
     $                      nf_float,nf_noerr)
C
        vdims(1)=x_dimid
        vdims(2)=y_dimid
        vdims(3)=z_dimid
C
        call def_var_netcdf(ncid,'rmean',3,vdims,rmean_varid,
     $                25,'horizontally-averaged rho',13,'dimensionless',
     $                      17,'east_e north_e zz',.true.,
     $                      nf_float,nf_noerr)
C
        vdims(1)=x_dimid
        vdims(2)=y_dimid
        vdims(3)=time_dimid
C
        call def_var_netcdf(ncid,'uab',3,vdims,uab_varid,
     $                      16,'depth-averaged u',9,'metre/sec',
     $                      14,'east_u north_u',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'vab',3,vdims,vab_varid,
     $                      16,'depth-averaged v',9,'metre/sec',
     $                      14,'east_v north_v',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'elb',3,vdims,elb_varid,
     $                      17,'surface elevation',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
!lyo:!wad:beg.
        call def_var_netcdf(ncid,'wetmask',3,vdims,wetmask_varid,
     $                       7,'wetmask',13,'dimensionless',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'d',3,vdims,d_varid,
     $                      17,'Total Water Depth',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'elbmsl',3,vdims,elbmsl_varid,
     $                      22,'surf.elevation wrt MSL',5,'metre',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
!lyo:!wad:end.
!lyo:!lyo:_20080507:
        call def_var_netcdf(ncid,'xwindstrs',3,vdims,wusurf_varid,
     $                       21,'Kinematic Xwindstress',7,'(m/s)^2',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'ywindstrs',3,vdims,wvsurf_varid,
     $                       21,'Kinematic Ywindstress',7,'(m/s)^2',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C
C!sc:wave variables out:start
        call def_var_netcdf(ncid,'thtav',3,vdims,thtav_varid,
     $              34,'average wave propagation direction',5,'(rad)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'ent',3,vdims,ent_varid,
     $          36,'en averaged over all wave directions',9,'(m+3 s-2)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'Hs',3,vdims,Hs_varid,
     $          23,'significant wave height',3,'(m)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'sigp',3,vdims,sigp_varid,
     $          35,'peak frequency of wind-driven waves',5,'(s-1)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'kthav',3,vdims,kthav_varid,
     $          34,'wave number, averaged on en(i,j,m)',5,'(s-1)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'tpx0',3,vdims,tpx0_varid,
     $          32,'wave pressure stress at surface',9,'(m+2 s-2)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'tpy0',3,vdims,tpy0_varid,
     $          32,'wave pressure stress at surface',9,'(m+2 s-2)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'ttx0',3,vdims,ttx0_varid,
     $          25,'turbulence surface stress',9,'(m+2 s-2)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
        call def_var_netcdf(ncid,'tty0',3,vdims,tty0_varid,
     $          25,'turbulence surface stress',9,'(m+2 s-2)',
     $                      14,'east_e north_e',.true.,
     $                      nf_float,nf_noerr)
C!sc:wave variables out:end
C
        vdims(1)=x_dimid
        vdims(2)=y_dimid
        vdims(3)=z_dimid
        vdims(4)=time_dimid
C
        call def_var_netcdf(ncid,'u',4,vdims,u_varid,
     $                      10,'x-velocity',9,'metre/sec',
     $                      17,'east_u north_u zz',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'v',4,vdims,v_varid,
     $                      10,'y-velocity',9,'metre/sec',
     $                      17,'east_v north_v zz',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'w',4,vdims,w_varid,
     $                      10,'z-velocity',9,'metre/sec',
     $                      16,'east_e north_e z',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'t',4,vdims,t_varid,
     $                      21,'potential temperature',1,'K',
     $                      17,'east_e north_e zz',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'s',4,vdims,s_varid,
     $                      23,'salinity x rho / rhoref',3,'PSS',
     $                      17,'east_e north_e zz',.true.,
     $                      nf_float,nf_noerr)
C
        call def_var_netcdf(ncid,'rho',4,vdims,rho_varid,
     $                    21,'(density-1000)/rhoref',13,'dimensionless',
     $                      17,'east_e north_e zz',.true.,
     $                      nf_float,nf_noerr)
C    !zb:EX1-1A
        call def_var_netcdf(ncid,'km',4,vdims,km_varid,
     $                    28,'vertical kinematic viscosity',4,'m2/s',
     $                      16,'east_e north_e z',.true.,
     $                      nf_float,nf_noerr)
C    !zb:EX1-1A
        call def_var_netcdf(ncid,'kh',4,vdims,kh_varid,
     $                    20,'vertical diffusivity',4,'m2/s',
     $                      16,'east_e north_e z',.true.,
     $                      nf_float,nf_noerr)
C    !zb:EX1-1A
        call def_var_netcdf(ncid,'q2',4,vdims,q2_varid,
     $                    15,'twice turb kinetic energy',5,'m2/s2',
     $                      16,'east_e north_e z',.true.,
     $                      nf_float,nf_noerr)
C    !zb:EX1-1A
        call def_var_netcdf(ncid,'q2l',4,vdims,q2l_varid,
     $                    26,'q2 times turb length scale',5,'m3/s2',
     $                      16,'east_e north_e z',.true.,
     $                      nf_float,nf_noerr)
C
C!sc:wave variables out:start
        call def_var_netcdf(ncid,'cg',4,vdims,cg_varid,
     $                32,'spectrally averaged group speeds',5,'m/s-1',
     $                      16,'east_e north_e z',.true.,
     $                      nf_float,nf_noerr)
C!sc:wave variables out:end

C     End definitions:
C
        status=nf_enddef(ncid)
        call handle_netcdf_error('nf_enddef           ',
     $                           status,nf_noerr)
C
C     Write initial data:
C
        status=nf_put_var_real(ncid,z_varid,z)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,zz_varid,zz)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,dx_varid,dx)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,dy_varid,dy)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,east_u_varid,east_u)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,east_v_varid,east_v)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,east_e_varid,east_e)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,east_c_varid,east_c)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,north_u_varid,north_u)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,north_v_varid,north_v)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,north_e_varid,north_e)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,north_c_varid,north_c)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,rot_varid,rot)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,h_varid,h)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,fsm_varid,fsm)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,dum_varid,dum)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,dvm_varid,dvm)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
        status=nf_put_var_real(ncid,rmean_varid,rmean)
        call handle_netcdf_error('nf_put_var_real     ',
     $                           status,nf_noerr)
C
C-----------------------------------------------------------------------
C
      else if(option.eq.2) then
C
C     Write a set of data:
C
        iout=iout+1
        start(1)=iout
        count(1)=1
C
        status=nf_put_vara_real(ncid,time_varid,start,count,time)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        start(1)=1
        start(2)=1
        start(3)=iout
        count(1)=im
        count(2)=jm
        count(3)=1
C
      write(6,*)'nc put uab'
        status=nf_put_vara_real(ncid,uab_varid,start,count,uab)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        status=nf_put_vara_real(ncid,vab_varid,start,count,vab)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        status=nf_put_vara_real(ncid,elb_varid,start,count,elb)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
!lyo:_20091229:!pom08_2d.n:tide_only for nwatl:
!       if(nwad.ne.0) then
!
!lyo:!wad:beg.
        status=nf_put_vara_real(ncid,wetmask_varid,start,count,wetmask)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
          do j=1,jm; do i=1,im
             tps(i,j)=d(i,j)*fsm(i,j)
             enddo; enddo
        status=nf_put_vara_real(ncid,d_varid,start,count,tps)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
          do j=1,jm; do i=1,im
             tps(i,j)=(elb(i,j)+hhi)*wetmask(i,j)
             enddo; enddo
        status=nf_put_vara_real(ncid,elbmsl_varid,start,count,tps)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
!lyo:!wad:end.
!lyo:!lyo:_20080507:
          do j=1,jm; do i=1,im
             tps(i,j)=-wusurf(i,j)
             enddo; enddo
        status=nf_put_vara_real(ncid,wusurf_varid,start,count,tps)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
          do j=1,jm; do i=1,im
             tps(i,j)=-wvsurf(i,j)
             enddo; enddo
        status=nf_put_vara_real(ncid,wvsurf_varid,start,count,tps)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C!sc:wave variables out:start
      write(6,*)'nc put wave'
        status=nf_put_vara_real(ncid,thtav_varid,start,count,thtav)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
        status=nf_put_vara_real(ncid,Hs_varid,start,count,Hs)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave 1'
        status=nf_put_vara_real(ncid,ent_varid,start,count,ent)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave 2'
        status=nf_put_vara_real(ncid,sigp_varid,start,count,sigp)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave 3'
        status=nf_put_vara_real(ncid,kthav_varid,start,count,kthav)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave 4'
        status=nf_put_vara_real(ncid,tpx0_varid,start,count,tpx0)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave 5'
        status=nf_put_vara_real(ncid,tpy0_varid,start,count,tpy0)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave 6'
        status=nf_put_vara_real(ncid,ttx0_varid,start,count,ttx0)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave 7'
        status=nf_put_vara_real(ncid,tty0_varid,start,count,tty0)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
      write(6,*)'nc put wave end'
C!sc:wave variables out:end
c
!lyo:_20091229:!pom08_2d.n:tide_only for nwatl:
!	endif !        if(nwad.ne.0) then
!
!pom08_2d.n:
        if(mode.ne.2) then
!
        start(1)=1
        start(2)=1
        start(3)=1
        start(4)=iout
        count(1)=im
        count(2)=jm
        count(3)=kb
        count(4)=1
C
        status=nf_put_vara_real(ncid,u_varid,start,count,u)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        status=nf_put_vara_real(ncid,v_varid,start,count,v)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        status=nf_put_vara_real(ncid,w_varid,start,count,wr) !lyo:_20080706:20080726:change "w" to "wr"
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        status=nf_put_vara_real(ncid,t_varid,start,count,t)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        status=nf_put_vara_real(ncid,s_varid,start,count,s)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C
        status=nf_put_vara_real(ncid,rho_varid,start,count,rho)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C    !zb:EX1-1A
        status=nf_put_vara_real(ncid,km_varid,start,count,km)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C    !zb:EX1-1A
        status=nf_put_vara_real(ncid,kh_varid,start,count,kh)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C    !zb:EX1-1A
        status=nf_put_vara_real(ncid,q2_varid,start,count,q2)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C    !zb:EX1-1A
        status=nf_put_vara_real(ncid,q2l_varid,start,count,q2l)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C!sc:wave variables out:start
        status=nf_put_vara_real(ncid,cg_varid,start,count,cg)
        call handle_netcdf_error('nf_put_vara_real    ',
     $                           status,nf_noerr)
C!sc:wave variables out:end
C
C
C-----------------------------------------------------------------------
C
!pom08_2d.n:
	endif !if(mode.ne.2) ...
!
      else if(option.eq.3) then
C
C     Close file:
C
        status=nf_close(ncid)
        call handle_netcdf_error('nf_close            ',
     $                           status,nf_noerr)
C
C-----------------------------------------------------------------------
C
      else
C
        write(6,3)
    3   format(/'Invalid option for subroutine write_netcdf ..... ',
     $          'program terminated'/)
        stop
C
      endif
C
C-----------------------------------------------------------------------
C
      return
C
      end
C
!_______________________________________________________________________
      subroutine read_grid
      implicit none
      include 'pom08.c'
      include '/aracbox/lib/netcdf/4.1.2/intel_12/include/netcdf.inc'        !sc:netcdf.inc on ATOP cluster 
      character*120 netcdf_grid_file 
      integer z_varid,zz_varid,dx_varid,dy_varid,east_c_varid,
     $        east_e_varid,east_u_varid,east_v_varid,north_c_varid,
     $        north_e_varid,north_u_varid,north_v_varid,rot_varid,
     $        h_varid,fsm_varid,dum_varid,dvm_varid
      integer ncid,status
      integer vdims(2)
      integer start(2),edge(2)
! open netcdf file
      write(netcdf_grid_file,'(a)') "in/grid/grid.nc"
      status=nf_open(netcdf_grid_file,nf_nowrite,ncid)
      call handle_netcdf_error('nf_open: '//netcdf_grid_file,
     $                          status,nf_noerr)
! get variables
      status=nf_inq_varid(ncid,'z',z_varid)
      call handle_netcdf_error('nf_inq_varid: z     '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'zz',zz_varid)
      call handle_netcdf_error('nf_inq_varid: zz    '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dx',dx_varid)
      call handle_netcdf_error('nf_inq_varid: dx    '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dy',dy_varid)
      call handle_netcdf_error('nf_inq_varid: dy    '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_u',east_u_varid)
      call handle_netcdf_error('nf_inq_varid: east_u'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_v',east_v_varid)
      call handle_netcdf_error('nf_inq_varid: east_v'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_e',east_e_varid)
      call handle_netcdf_error('nf_inq_varid: east_e'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'east_c',east_c_varid)
      call handle_netcdf_error('nf_inq_varid: east_c'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_u',north_u_varid)
      call handle_netcdf_error('nf_inq_varid:north_u'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_v',north_v_varid)
      call handle_netcdf_error('nf_inq_varid:north_v'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_e',north_e_varid)
      call handle_netcdf_error('nf_inq_varid:north_e'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'north_c',north_c_varid)
      call handle_netcdf_error('nf_inq_varid:north_c'
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'rot',rot_varid)
      call handle_netcdf_error('nf_inq_varid: rot   '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'h',h_varid)
      call handle_netcdf_error('nf_inq_varid: h     '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'fsm',fsm_varid)
      call handle_netcdf_error('nf_inq_varid: fsm   '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dum',dum_varid)
      call handle_netcdf_error('nf_inq_varid: dum   '
     $                          ,status,nf_noerr)
      status=nf_inq_varid(ncid,'dvm',dvm_varid)
      call handle_netcdf_error('nf_inq_varid: dvm   '
     $                          ,status,nf_noerr)
! get data
      start(1)=1
      edge(1)=kb
      status=nf_get_vara_real(ncid,z_varid,start,edge,z)
      call handle_netcdf_error('nf_get_vara_real    ',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,zz_varid,start,edge,zz)
      call handle_netcdf_error('nf_get_vara_real    ',
     $                          status,nf_noerr)
      start(1)=1
      start(2)=1
      edge(1)=im
      edge(2)=jm
      status=nf_get_vara_real(ncid,dx_varid,start,edge,dx)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,dy_varid,start,edge,dy)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,east_u_varid,start,edge,
     $                               east_u)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,east_v_varid,start,edge,
     $                               east_v)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,east_e_varid,start,edge,
     $                               east_e)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,east_c_varid,start,edge,
     $                               east_c)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,north_u_varid,start,edge,
     $                               north_u)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,north_v_varid,start,edge,
     $                               north_v)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,north_e_varid,start,edge,
     $                               north_e)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,north_c_varid,start,edge,
     $                               north_c)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,rot_varid,start,edge,rot)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,h_varid,start,edge,h)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,fsm_varid,start,edge,fsm)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,dum_varid,start,edge,dum)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
      status=nf_get_vara_real(ncid,dvm_varid,start,edge,dvm)
      call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
! close file:
      status=nf_close(ncid)
      call handle_netcdf_error('nf_close: grid      ',status,nf_noerr)
      return
      end
!_______________________________________________________________________
      subroutine read_wind(wfile,ti,wu,wv)
      implicit none
      include 'pom08.c'
      include '/aracbox/lib/netcdf/4.1.2/intel_12/include/netcdf.inc'        !sc:netcdf.inc on ATOP cluster 
      character*16 wfile 
      integer wu_varid,wv_varid
      integer ncid,status
      integer vdims(3)
      integer start(3),edge(3)
      integer :: ti
      real :: wu(im,jm),wv(im,jm)
      logical :: lexist
      inquire(file='in/'//trim(wfile),
     $  exist=lexist)
! open netcdf file
      if (lexist)then
       status=nf_open('in/'//trim(wfile),nf_nowrite,ncid)
       call handle_netcdf_error('open'//wfile,
     $                          status,nf_noerr)
! get variables
       status=nf_inq_varid(ncid,'uwsrf',wu_varid)
       call handle_netcdf_error('nf_inq_varid: uwsrf '
     $                          ,status,nf_noerr)
       status=nf_inq_varid(ncid,'vwsurf',wv_varid)
       call handle_netcdf_error('nf_inq_varid: vwsrf '
     $                          ,status,nf_noerr)
! get data
       start(1)=1
       start(2)=1
       start(3)=ti
       edge(1)=im
       edge(2)=jm
       edge(3)=1
       status=nf_get_vara_real(ncid,wu_varid,start,edge,
     $                               wu)
       call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
       status=nf_get_vara_real(ncid,wv_varid,start,edge,
     $                               wv)
       call handle_netcdf_error('nf_get_vara_real_all',
     $                          status,nf_noerr)
! close file:
       status=nf_close(ncid)
       call handle_netcdf_error('nf_close: grid      ',status,nf_noerr)
      endif
      return
      end
!_______________________________________________________________________
C     End of source code
C
C-----------------------------------------------------------------------
C
