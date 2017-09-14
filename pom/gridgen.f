
c      open(11, file='model_grid',form='formatted')  !used in s.r. nybight

      subroutine gridgen
C **********************************************************************
C *                                                                    *
C * FUNCTION    :  Sets up for NY BIght                                *
C *                                                                    *
C **********************************************************************
C
      implicit none
C
      include 'pom08.c'
      
C
      integer list,n,ii,jj,id,jd,imm,jmm
      parameter(imm=im,jmm=jm)
      real latd,longd,re,vel
      real dxd,dyd,hd,rotd
      real delxat,delxon,delyat,delyon
c      real hmm(imm,jmm),dxm(imm,jmm),dym(imm,jmm),rotm(imm,jmm)
c      real corm(imm,jmm),latm(imm,jmm),longm(imm,jmm) 
      
c      write(6,'(''im,jm,imm,jmm='',4i5)') im,jm,imm,jmm 
c      id=1; jd=1
c      if(im.eq.40) then; id=5; jd=5; endif   ! im,jm set in pom08.c              
c      write(6,'(''id,jd='',2i5)') id,jd 
C
      re=6365.e3            ! earth's readius
c      do j=1,jmm
c        do i=1,imm
c          hmm(i,j)=1.0
c          dxm(i,j)=1. 
c          dym(i,j)=1. 
c          corm(i,j)=1.e-4
c          latm(i,j)=0.
c          longm(i,j)=0.
c        enddo  
c      enddo
c!sc: read grid.nc and using north_e and east_e as long and lat
      call read_grid
      lat=north_e
      long=east_e     
c!sc: the rest is covered in read_grid
C  N.B. There are no j=1 or i=1 data in model_grid
c      do n=1,8; read(11,'(2x)');  enddo
c      do list=9,156794
c        read(11,'(2i5,4f10.2,2f10.5)') i,j,dxd,dyd,hd,rotd,latd,longd
c       write(6,'(2i5,4f10.2,2f10.5)') i,j,dxd,dyd,hd,rotd,latd,longd
c        if(i.le.imm) then      ! retain only i=1 to imm       
c        dxm(i,j)=dxd  
c        dym(i,j)=dyd
c        hmm(i,j)=hd
c        rotm(i,j)=rotd  
c        latm(i,j)=latd
c        longm(i,j)=longd
C      if(i.eq.4.and.j.eq.2) write(6,'(''delyat'',2i5,3e12.3)')
C    &    i,j,latm(4,2),longm(4,2),hmm(4,2)
c        endif
c      enddo

c      do jj=1,jmm,jd
c        do ii=1,imm,id      !retain every id,jd points
c          if(id.eq.5) then
c            i=1+ii/id
c            j=1+jj/jd
c          else
c            i=ii/id
c            j=jj/jd
c          endif
c          h(i,j)=hmm(ii,jj) 
c          rot(i,j)=rotm(ii,jj) 
c          lat(i,j)=latm(ii,jj)
c          long(i,j)=longm(ii,jj)
c        enddo
c      enddo

          
      do j=1,jm
        do i=1,im
c!sc:testing
c          if(lat(i,j).ne.0.) 
c     &    cor(i,j)=4.*pi*sin(lat(i,j)*pi/180.)/86400.
c          if(lat(i,j).ne.0.) 
          cor(i,j)=4.*pi*sin(lat(i,j)*pi/180.)/86400.
        end do
      end do
      dum(1,:)=dum(2,:)
! calculate areas of "t" and "s" cells
      do j=1,jm
        do i=1,im
          art(i,j)=dx(i,j)*dy(i,j)
        end do
      end do

! calculate areas of "u" and "v" cells
      do j=2,jm
        do i=2,im
          aru(i,j)=.25e0*(dx(i,j)+dx(i-1,j))*(dy(i,j)+dy(i-1,j))
          arv(i,j)=.25e0*(dx(i,j)+dx(i,j-1))*(dy(i,j)+dy(i,j-1))
        end do
      end do
C
c      do i=1,im
c        h(i,1)=h(i,2)
c        h(i,jm)=1.                 
c      enddo
c      do j=1,jm
c        h(1,j)=h(2,j)
c        h(im,j)=h(im-1,j)
c      enddo
      call areas_masks         ! obtain fsm,dum,dvm
C
C     Adjust bottom topography so that cell to cell variations
C     in h do not exceed parameter slmax:
C
C     if(slmax.lt.1.e0) call slpmax
C
      hmax=100.        ! used in profu and profv
C
c      do j=2,jm-1
c        do i=2,im-1

c         if(lat(i+1,j).gt.0..and.lat(i-1,j).gt.0.) then
c          delxon= (long(i+1,j)-long(i-1,j))
c     &                *sin(lat(i,j)*pi/180.)
c          delxat= (lat(i+1,j)-lat(i-1,j))
c          dx(i,j)=0.5*re*pi/180.*sqrt(delxon**2+delxat**2)  
c         endif
  
c         if(lat(i,j+1).gt.0..and.lat(i,j-1).gt.0.) then 
c          delyat=(lat(i,j+1)-lat(i,j-1))
c          delyon=(long(i,j+1)-long(i,j-1))
c     &           *sin(lat(i,j)*pi/180.)
c          dy(i,j)=0.5*re*pi/180.*sqrt(delyon**2+delyat**2)        
c         endif
c        enddo
c      enddo
c      do j=2,jm
c        dx(2,j)=dx(3,j)
c        dx(1,j)=dx(2,j)
c        dy(2,j)=dy(3,j)
c        dy(1,j)=dy(2,j)
c        dy(im,j)=dy(im-1,j)
c        dx(im,j)=dx(im-1,j)
c        do i=1,im
c         if(dx(i,j).lt.1.) dx(i,j)=dx(i,j-1)
c         if(dy(i,j).lt.1.) dy(i,j)=dy(i,j-1)
c        enddo
c      enddo
c      do i=1,im
c      dx(i,2)=dx(i,3)
c      dx(i,1)=dx(i,2)
c      dy(i,2)=dy(i,3)
c      dy(i,1)=dy(i,2)
c      enddo

c      call areas_masks      ! obtain cell areas
C     Set initial conditions:
C
c!sc:vel is used for ub
      vel=0.1e0 ! current velocity

      do k=1,kbm1
!sc:add dz and dzz calculation
        dz(k)=z(k)-z(k+1)
        dzz(k)=zz(k)-zz(k+1)
        do j=1,jm
          do i=1,im
            tb(i,j,k)=5.e0+15.e0*exp(zz(k)*h(i,j)/1000.e0)-tbias
            sb(i,j,k)=35.e0-sbias
            tclim(i,j,k)=tb(i,j,k)
            sclim(i,j,k)=sb(i,j,k)
            ub(i,j,k)=vel*dum(i,j)
          end do
        end do
      end do
c!sc: added according to mpipom read grid
      dz(kb) = dz(kb-1) !=0. !lyo:20110202 =0 is dangerous but checked
      dzz(kb)=dzz(kb-1) !=0. !thro' code - and found to be ok.
C
C     Initialise uab and vab as desired
C     (NOTE that these have already been initialised to zero in the
C     main program):
c!sc:vel was used before this
c      vel=0.0
C
      do j=1,jm
        do i=1,im
          uab(i,j)=vel*dum(i,j)
        end do
      end do
C
C     Set surface boundary conditions, e_atmos, vflux, wusurf,
C     wvsurf, wtsurf, wssurf and swrad, as necessary
C     (NOTE:
C      1. These have all been initialised to zero in the main program.
C      2. The temperature and salinity of inflowing water must be
C         defined relative to tbias and sbias.):
C
C
C     Initialise elb, etb, dt and aam2d:
C
      do j=1,jm
        do i=1,im
          elb(i,j)=-e_atmos(i,j)
          etb(i,j)=-e_atmos(i,j)
          dt(i,j)=h(i,j)-e_atmos(i,j)
          aam2d(i,j)=aam(i,j,1)
        end do
      end do
C
      call dens(sb,tb,rho)
C
C     Generated horizontally averaged density field (in this
C     application, the initial condition for density is a function
C     of z (the vertical cartesian coordinate) -- when this is not
C     so, make sure that rmean has been area averaged BEFORE transfer
C     to sigma coordinates):
C
      do k=1,kbm1
        do j=1,jm
          do i=1,im
            rmean(i,j,k)=rho(i,j,k)
          end do
        end do
      end do
C
C     Set lateral boundary conditions, for use in subroutine bcond
C     (in the seamount problem, the east and west boundaries are open,
C     while the south and north boundaries are closed through the
C     specification of the masks fsm, dum and dvm):
C
      rfe=1.e0
      rfw=1.e0
      rfn=1.e0
      rfs=1.e0
C
      do j=2,jmm1
        uabw(j)=uab(2,j)
        uabe(j)=uab(imm1,j)

c!sc:turn on according mpipom setting
C
C     Set geostrophically conditioned elevations at the boundaries:
C
       ele(j)=ele(j-1)-cor(imm1,j)*uab(imm1,j)/grav*dy(imm1,j-1)
       elw(j)=elw(j-1)-cor(2,j)*uab(2,j)/grav*dy(2,j-1)
      end do
C
C     Adjust boundary elevations so that they are zero in the middle
C     of the channel:
C
C
C     Set thermodynamic boundary conditions (for the seamount
C     problem, and other possible applications, lateral thermodynamic
C     boundary conditions are set equal to the initial conditions and
C     are held constant thereafter - users may, of course, create
C     variable boundary conditions):
C
      do k=1,kbm1
C
        do j=1,jm
          tbe(j,k)=tb(im,j,k)
          tbw(j,k)=tb(1,j,k)
          sbe(j,k)=sb(im,j,k)
          sbw(j,k)=sb(1,j,k)
        end do
C
        do i=1,im
          tbn(i,k)=tb(i,jm,k)
          tbs(i,k)=tb(i,1,k)
          sbn(i,k)=sb(i,jm,k)
          sbs(i,k)=sb(i,1,k)
        end do
C
      end do
C
      return
      end
C
