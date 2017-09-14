C     Begin wave source code
C
C-----------------------------------------------------------------
C                  BEGINNING OF WAVE CODE             !WAVE
C-----------------------------------------------------------------

      subroutine wave(dtw,dth,nwave)        
C------------------------------------------------------------
C Subroutine to calculate wave energy density as a function of
C x, y and theta. Grid is orthogonal curvilinear. For rectilinear
C grid, simply set dx(i,j) or dy(i,j) = constant
C
C  Input:
C    z(k),dz(k) = vertical coordinate and spacing (non-d.)
C    zz(k),dzz(k) = middle cell coordinate and spacing
C    dx(i,j), dy(i,j) = spatial grid cell increments located ati
C        cell center (m)
C    h(i,j) = bottom depth
C    d(i,j) = water column depth (m)
C    fsm(i,j) = 1 for water cell, = 0 for land cell
C    u10x(i,j),u10y(i,j) = wind vector at 10 m height (m s-1)
C    u(i,j,k),v(i,j,k) = current vector (m s-1)
c    nitera = # of iterations in s.r. stepmpd

C  Output wave variables:
C    enb(i,j,m), en(i,j,m), enf(i,j,m) = past, present,future
C        energy densities  (m+3 s-2)
C    ent(i,j) en averaged over all wave directions (m+3 s-2)
C    Hs(i,j) = significant wave height = 4*sqrt(ent/grav) (m)
C    thtav(i,j) = average wave propagation direction (rad)
C    Sxx,Sxy,Syy = wave radiation stress for momentum equation
C       related to wave velocity moments (m+3 s-2)
C    (tpx0,tpy0) = wave pressure stress at surface (m+2 s-2)
C    (tpx,tpy)=(tpx0,tpy0)*tpzdist =  wave pressure below the
C       surface (m+2 s-2)

C  Internal wave variables (some could also be output)
C    theta(m) = wave directional coordinate located at cell edge (rad) 
C    cs(m),sn(m) = cosine and sine of theta
C    thetam(m) = wave directional coordinate located at cell center 
C    csm(m),snm(m) = cosine and sine of thetam
C    sigth(i,j,m)=theta dependent frequency, swell and/or wind
C          waves (s-1)
C    sigthav(i,j)=frequency, averaged on en(i,j,m)   (s-1)
C    cth(i,j,m)=theta dependent phase speed (m s-1)
C    kth(i,j,m)= theta dependent wave number   (m-1)
C    kthav(i,j)=wave number, averaged on en(i,j,m)   (s-1)
C    kthavD(i,j)= average wave number *d(i,j) 
C    cg(i,j,m) = spectrally averaged group speeds located at
C         cell center (m s-1)
C    sigp(i,j)=peak frequency of wind-driven waves (s-1)
C    cp(i,j) = peak frequency phase speed from sigp (m s-1)
C    kp(i,j) = peak frequency wave number from sigp (m-1)
C    kpD(i,j) = kp(i,j)*d(i,j)  (non-d)
C    cgx(i,j,m), cgy(i,j,m) = propagation vector in x, y
C         directions located at cell edge (m s-1)
C    ud(i,j),vd(i,j) = Dopler velocities (m s-1)
C    cthth(i,j,m) = propagation speed in the theta direction 
C        located at cell edge (rad s-1)
C    u10(i,j) = wind speed
C    ageinv(i,j)=u10*sigp/grav (= 0.83 for fully developed
C        waves) (non-d)
C    fspr(i,j)=spreading function for swin (non-d)
C    beta = constant in fspr; higher beta = narrower spread
C        but recommended value = 2,2 (non-d)
C    swin(i,j,m) = wave energy source term (m3 s-3)
C    sdis(i,j,m) = surface wave energy dissipation (m3 s-3)
C    bdis1(i,j,m) = bottom friction dissipation (m3 s-3)
C    bdis2(i,j,m) = bottom depth-induced dissipation (m3 s-3)
C    gama = empirical constant in formulation of bdis2 (non-d)
C    Sradx(i,j),Srady(i,j),Sradxy(i,j),Srad(i,j) = wave radiation
C        stress for wave energy equation related to wave velocity 
C        moments (s-1)
C    Fcg(i,jm), Fct(i,j,m) = correction for cg and cth to  
C        account for spectral averaging (non-d)
C    Cin, swcon = constants in swin curve fit (non-d)
C    adis, bdis = empirically determined constants in sdis term
C        (non-d)
C    swfct = factor to reduce dissipation for swells (non-d)
C    dif = horizontal advective diffusion coeficient (see note
C        in x-y step) (non-d)
C------------------------------------------------------------
      implicit none
C
      include 'pom08.c' 
C
      integer 
     &    nitera,mbig,mploc,mmloc,loop,
     &    imax,jmax,init 
      real
     &    dxinv,dyinv,dx2inv,dy2inv,Hsmax,
     &    diff,summ,fsprm,tht_ops,fb,Tp,bb,Hsbc,
     &    cgxm,cgym,sorc,adv,sor,cgnd,fcoshinv,fsinhinv,
     &    snmav,csmav,Tsig,Cf,Dcoeff,cdmax
      real                    !  2-d variables
     &    usdif(im,jm),vsdif(im,jm),
     &    Sradx(im,jm),Sradxy(im,jm),Srady(im,jm),Srad(im,jm),
     &    Srad0(im,jm)
      real                    !  3-d variables
     &    xflux(im,jm,mm),yflux(im,jm,mm),tflux(im,jm,mm),
     &    uak(im,jm,mm),
     &    Fcg(im,jm,mm),Fct(im,jm,mm),tpsw(im,jm,mm),
     &    tpsw2(im,jm,mm),tpsw3(im,jm,mm),tpsw4(im,jm,mm),
     &    tpsw5(im,jm,mm),
     &    westang,eastang,southang,northang
      real smoth,beta,Cin,swcon,adis,bdis,swfct,dif,gama
      save dtw2
      data smoth/.20/
      data beta/2.2/                ! spreading const. 
      data Cin/370./                ! Sin const.
      data swcon/0.33/              ! Sin const.
      data adis/0.925/              ! Sdis const.
      data bdis/0.18e-4/            ! Sdis const.
      data swfct/0.2/               ! factor to decrease swell diss.
      data nitera/3/                ! # of iterations in MPDATA
      data dif/0.5/                 ! diffusion coeficient 
                                    ! in x-y step section
      data gama/0.70/               ! depth-induced breaking const.
      data cdmax/.0025/             ! max. drag coefficient
      data initwave/0/

C
C------------------------------------------------------------------
                      if(iint.eq.1) then
C------------------------------------------------------------------
      rewind 10
      open(10, file='specavs',form='formatted')
      dtw2=2.*dtw
      dthinv=1./dth                ! inverse angle increment
      write(6,'(''wave parameters'')')
      write(6,'('' Cin ='',f8.1,''  adis ='',f8.3,'' bdis ='',e12.4)')
     &      Cin,adis,bdis
      write(6,'('' beta ='',f10.3)') beta
      write(6,'('' nitera ='',i5 )') nitera     
      write(6,'('' dif ='',f8.3)') dif
      write(6,'('' cdmax ='',f8.3)') cdmax   
      do j=1,jm
        do i=1,im-1
          aren(i,j)=dx(i,j)*dy(i,j)
          areninv(i,j)=1./aren(i,j)
        enddo
      enddo
      write(6,'(5x,''   m      theta       cs       sn      thetam'',
     &''      csm       snm'')')
      do m=1,mm
        thetam(m)=dth*float(m-mm/2)
        theta(m)=dth*float(m-mm/2)-0.5*dth
        cs(m)=cos(theta(m))
        sn(m)=sin(theta(m))
        csm(m)=cos(thetam(m))
        snm(m)=sin(thetam(m))
        write (6,'(5x,i5,6f10.3)') m, theta(m),cs(m),sn(m),
     &        thetam(m),csm(m),snm(m)
      enddo
C
C Open boundary conditions. (For a closed basin, not used)
C Add north conditions as needed.
      Hsbc=0.2         !Sig. wave height at boundaries
      Tp=10            !wave period at boundaries
      westang=0.       !western wave propagation angle; incoming
C     westang=pi       !western wave propagation angle; outgoing
      eastang=0.       !eastern wave propagation angle; outgoing 
C     eastang=pi       !eastern wave propagation angle; incoming   
      southang=pi/2.   !southern wave propagation angle; incoming
      do m=1,mm
        do j=1,jm
          enw(j,m)=grav*(Hsbc/4.)**2*fsprm(beta,m,mm,westang,thetam(m))
          enw(j,m)=max(.0001,enw(j,m))*fsm(1,j)
          ene(j,m)=grav*(Hsbc/4.)**2*fsprm(beta,m,mm,eastang,thetam(m))
          ene(j,m)=max(.0001,ene(j,m))*fsm(im,j) 
          sigthe(j,m)=2.*pi/Tp
          sigthw(j,m)=2.*pi/Tp
        enddo
      enddo
      do m=1,mm
        do i=1,im
          ens(i,m)=grav*(1.00/4.)**2*fsprm(beta,m,mm,southang,thetam(m))
          ens(i,m)=max(.0001,ens(i,m))*fsm(i,1)   
          sigths(i,m)=2.*pi/Tp
        enddo
      enddo
                         endif
C-------------------------Initialize----------------------------- 
                      if(initwave.eq.0) then
C----------------------------------------------------------------   
      do k=1,kb
        F1(k)=0.
        F2(k)=0.
        F3(k)=0.
        F4(k)=0.
      enddo
      do j=1,jm    
        do i=1,im
          ud(i,j)=0.
          vd(i,j)=0.
          sdistot(i,j)=0.
          bdistot(i,j)=0.
          swintot(i,j)=0.
          thtav(i,j)=0.
          qdis(i,j)=0.5
          do m=1,mm
            cth(i,j,m)=10.
            kth(i,j,m)=.1
            kthD(i,j,m)=1.
            Fcg(i,j,m)=0.5
            Fct(i,j,m)=1. 
            bdis1(i,j,m)=0.
            bdis2(i,j,m)=0.
          enddo
          ageinv(i,j)=5.
          sigp(i,j)=1.5
          ustp(i,j)=0.
          cp(i,j)=grav/sigp(i,j)
          kp(i,j)=sigp(i,j)/cp(i,j)
          kpD(i,j)=kp(i,j)*d(i,j)
          kthavD(i,j)=kthD(i,j,1)
          Sradx(i,j)=0.
          Sradxy(i,j)=0.
          Srady(i,j)=0.
          Srad(i,j)=0.
          ud(i,j)=0.
          vd(i,j)=0.
          tht_wnd(i,j)=0.
          do k=1,kb
            Sxx(i,j,k)=0.
            Sxy(i,j,k)=0.
            Syy(i,j,k)=0.
          enddo
        enddo
      enddo
C
      do m=1,mm
        do j=1,jm
          do i=1,im
            sigthb(i,j,m)=sigp(i,j)
            sigth(i,j,m)=sigthb(i,j,m)
            sigthf(i,j,m)=sigth(i,j,m)
            cth(i,j,m)=grav/sigth(i,j,m)
            kth(i,j,m)=sigth(i,j,m)**2/grav
            kthD(i,j,m)=kth(i,j,m)*d(i,j)
C  Set up intitial enb, en, enf
            Hs(i,j)=0.25
c!sc: change hs on land to 0
            if (fsm(i,j).eq.1)then
             Hs(i,j)=0.
            end if  
            enb(i,j,m)=grav*(Hs(i,j)/4.)**2
     &                *fsprm(beta,m,mm,tht_wnd(i,j),thetam(m))
            enb(i,j,m)=max(.0001,enb(i,j,m))*fsm(i,j)
            en(i,j,m)=enb(i,j,m)
            enf(i,j,m)=en(i,j,m)
C
            xflux(i,j,m)=0.
            yflux(i,j,m)=0.
            tflux(i,j,m)=0.
            fspr(i,j,m)=fsprm(beta,m,mm,tht_wnd(i,j),thetam(m))
          enddo
        enddo
      enddo
      do j=1,jm
        do i=1,im
          ent(i,j)=summ(enb,dth,i,j,im,jm,mm)
        enddo
      enddo
C
C
      call prxy('tht_wnd',0.,tht_wnd,im,iskp,jm,jskp,0.)
      call prxy('Hs     ',0.,Hs     ,im,iskp,jm,jskp,0.)
      call prxz('fspr',0.,fspr,im,iskp,jm,mm,2,jm/2,jm-2,0.,d,thetam)
c     call prxz('enf ',0.,enf ,im,iskp,jm,mm,2,100,199,0.,d,thetam)
c     call prxy('ent',0.,ent,im,iskp,jm,jskp,0.)
C---------------------------------------------------------------- 
                   initwave=1
                   endif
C---------------------------------------------------------------- 
      call wave2cur (dth)  !WAVE radiation:    

      dtw2=2.*dtw
c      write(*,*)'xxx2',u10x(2,2),u(2,2,1),u10y(2,2),v(2,2,1)
      do j=2,jm-1
        do i=2,im-1
          u10(i,j)=sqrt(u10x(i,j)**2+u10y(i,j)**2)
          usdif(i,j)=(u10x(i,j)-u(i,j,1))*fsm(i,j) ! Used in s.r. stress
          vsdif(i,j)=(u10y(i,j)-v(i,j,1))*fsm(i,j)
          if(u10(i,j).gt.0.0001)
     &      tht_wnd(i,j)=acos(u10x (i,j)/(u10(i,j)))*sign(1.,u10y (i,j))
          if(u10 (i,j).gt.4.00 ) then
            do m=1,mm
              fspr(i,j,m)=fsprm(beta,m,mm,tht_wnd(i,j),thetam(m))
     &                    *fsm(i,j)
            enddo
          else
            do m=1,mm
              fspr(i,j,m)=0.
            enddo
          endif
        enddo
      enddo
      call stress (iint,nwave,Hs,ageinv,usdif,vsdif,cd,tpx0,tpy0,
     &                 cdmax,ttx0,tty0,ustw,ustp,im,jm,fsm)
C  Calculate sigp for wind driven portion of enf
      do j=1,jm
        do i=1,im
c         ageinv(i,j)=1.
          if (fsm(i,j).gt.0.) then
            if(u10(i,j).gt.4.0) then
              enwnd(i,j)= 0.001
              do m=2,mm-1
                if(fspr(i,j,m).gt.0.10)
     &               enwnd(i,j)=enwnd(i,j)+en(i,j,m)*dth
              enddo
              Dcoeff=.0022*(0.7+.12*sqrt(enwnd(i,j)))
              sigp(i,j) =(grav/u10(i,j))
     &          *(Dcoeff*u10(i,j)**4/(grav*enwnd(i,j)))**.303
              call wvfcns(sigp(i,j),d(i,j),kp(i,j),cp(i,j))   
              ageinv(i,j)=u10(i,j)*sigp(i,j)/grav       
            endif
          endif
        enddo
      enddo
c     call prxz('enf ',time,enf ,im,iskp,jm,mm,2,111,jm-2, 0.,d,thetam)
c     call prxy('ent',time,ent,im,iskp,jm,jskp  ,0. )
c     call prxy('sigp',time,sigp,im,iskp,jm,jskp  ,0. )
c     call prxy('kp  ',time,kp  ,im,iskp,jm,jskp  ,0. )
c     call prxy('cp  ',time,cp  ,im,iskp,jm,jskp  ,0. )
c     call prxy('cg ',time,cg ,im,iskp,jm,jskp,0. )
C---------------- Calculate propagation speeds-----------------------
      do j=1,jm
        do i=1,im  
          if(fsm(i,j).gt.0.) then
          ent(i,j)=summ(en,dth,i,j,im,jm,mm)
          if(ent(i,j).lt.0.) then
          write(6,'(''entlt0'',4i5,e12.3)') iint,nwave,i,j,ent(i,j)
          stop
          endif
          Hs(i,j)=0.5*(4.*sqrt(ent(i,j)/grav)+Hs(i,j))
C  Find kth and cth 
            do m=2,mm-1
              call wvfcns(sigthb(i,j,m),d(i,j),kth(i,j,m),cth(i,j,m))  
              kthD(i,j,m)=kth(i,j,m)*d(i,j)
C  Find spectrally averaged corrections to cg and cthth
C  cg used in x-y step; cthth used in theta step
              call speeds(kthD(i,j,m),Fcg(i,j,m),Fct(i,j,m))
              cg(i,j,m)=Fcg(i,j,m)*cth(i,j,m)
            enddo
          endif
        enddo
      enddo
C
C-------------------------------------------------------------
      do j=2,jm-1  
        do i=2,im-1  
          if(fsm(i,j).eq.0.) then
            tps(i,j)=(fsm(i+1,j)+fsm(i-1,j)+fsm(i,j+1)+fsm(i,j-1))
            do m=1,mm
              cg(i,j,m)=
     &          fsm(i+1,j)*cg(i+1,j,m)+fsm(i-1,j)*cg(i-1,j,m)
     &         +fsm(i,j+1)*cg(i,j+1,m)+fsm(i,j-1)*cg(i,j-1,m)
            enddo
            if(tps(i,j).gt.0) then
              do m=1,mm
                cg(i,j,m)=cg(i,j,m)/tps(i,j)  
              enddo
            endif
          endif
        enddo
      enddo
C---------------------------------------------------------------
C   The solution is split between x-y steps, theta steps and
C   source steps.
C ---------------------------------------------------------------
C                     x-y step for wave energy
C ---------------------------------------------------------------
C  If dif=0.2, then diff = 0.2 everywhere, but, on open water
C  cells imediately next to land cells, diff = 0.5 for local upwind
C  B.C. If dif=0.5, then diff=0.5 everywhere.
C  N.B. not a great difference if diff = 0.5 everywhere.
      do m=2,mm-1
        do j=2,jm
          do i=2,im
            cgx(i,j,m)=0.5*(cg(i,j,m)+cg(i-1,j,m))*csm(m)    
     &                +0.5*(ud(i,j)+ud(i-1,j))
            cgy(i,j,m)=0.5*(cg(i,j,m)+cg(i,j-1,m))*snm(m)   
     &                +0.5*(vd(i,j)+vd(i,j-1))
            diff=0.5*(1.-dum(i,j))+dif*dum(i,j)
            xflux(i,j,m)= 0.5*cgx(i,j,m)*(enb(i,j,m)+enb(i-1,j,m))
     &           -diff*abs(cgx(i,j,m))*(enb(i,j,m)-enb(i-1,j,m))
            diff=0.5*(1.-dvm(i,j))+dif*dvm(i,j)
            yflux(i,j,m)= 0.5*cgy(i,j,m)*(enb(i,j,m)+enb(i,j-1,m))
     &           -diff*abs(cgy(i,j,m))*(enb(i,j,m)-enb(i,j-1,m))
            xflux(i,j,m)=0.5*(dy(i,j)+dy(i-1,j))*xflux(i,j,m) !*dum(i,j)
            yflux(i,j,m)=0.5*(dx(i,j)+dx(i,j-1))*yflux(i,j,m) !*dvm(i,j)
          enddo
        enddo
      enddo
      
      do m=2,mm-1
        do j=2,jm-1
          do i=2,im-1
             enf(i,j,m)=enb(i,j,m)
     &        -(xflux(i+1,j,m)-xflux(i,j,m)
     &         +yflux(i,j+1,m)-yflux(i,j,m))*areninv(i,j)*dtw2
            enf(i,j,m)=enf(i,j,m)*fsm(i,j)   
          enddo
        enddo
      enddo
C-------------------Essential  Cyclic B.C.s on m ---------------
      do j=1,jm
        do i=1,im
        sigth(i,j,1)=sigth(i,j,mm-1)
        sigth(i,j,mm)=sigth(i,j,2)
        enf(i,j,1)=enf(i,j,mm-1)
        enf(i,j,mm)=enf(i,j,2)
        enddo
      enddo
C -- -------------------------------------------------------------
C                  theta step (MDATA) for wave energy
C ---------------------------------------------------------------
      do m=1,mm
        do j=1,jm
          do i=1,im 
            uak(i,j,m)=cs(m)*ud(i,j)+sn(m)*vd(i,j)
          enddo
        enddo
      enddo
      do j=2,jm-1
        do i=2,im-1 
          if(fsm(i,j).eq.1.) then
            dx2inv=1./(dx(i+1,j)+dx(i-1,j))
            dy2inv=1./(dy(i,j+1)+dy(i,j-1))
            do m=1,mm
              cthth(i,j,m)=                    ! speed in theta coordinate
     &          Fct(i,j,m)*grav/(2.*cth(i,j,m))*fcoshinv(kthD(i,j,m))**2
     &       *( sn(m)*(d(i+1,j)-d(i-1,j))*dx2inv*fsm(i+1,j)*fsm(i-1,j)
     &         -cs(m)*(d(i,j+1)-d(i,j-1))*dy2inv*fsm(i,j+1)*fsm(i,j-1))
     &         +sn(m)*(uak(i+1,j,m)-uak(i-1,j,m))*dx2inv
     &                  *fsm(i+1,j)*fsm(i-1,j)
     &         -cs(m)*(uak(i,j+1,m)-uak(i,j-1,m))*dy2inv
     &                  *fsm(i,j+1)*fsm(i,j-1)
c Protect against refraction CFL violation, but, if invoked, may cause 
C errors in shoaling waters
            cthth(i,j,m)=min(cthth(i,j,m),0.99*dth/dtw2)
          enddo
         endif
        enddo
      enddo
      do m=1,mm
        do j=1,jm
          cthth(im,j,m)=cthth(im-1,j,m)
          cthth(1,j,m)=cthth(2,j,m)
         enddo   
         do i=1,im
          cthth(i,jm,m)=cthth(i,jm-1,m)
          cthth(i,1,m)=cthth(i,2,m)
         enddo
       enddo
C
       call stepmpd(enb,en,enf,cthth,im,jm,mm,dtw2,dthinv,nitera)
C ---------------------------------------------------------------
C                  Source step for wave energy
C ---------------------------------------------------------------
      do i=2,im-1
        do j=2,jm-1
          if(fsm(i,j).eq.1.) then
C         Calculate Hm and qdis for later use in depth-induced breaking
            Hm(i,j)=gama*d(i,j)
            fb=8.0*ent(i,j)/grav/Hm(i,j)**2
c           fb=min(fb,1.0)     ! necessary ?
c!sc:testing limiting qdis 
c            if(abs((qdis(i,j)-1)/fb).le.40.) then
            
            qdis(i,j)=exp((qdis(i,j)-1)/fb)
C
            do m=2,mm-1
C Create neg. "tail" to oppose wave travelling opposit to wind. 
C Const. 0.4 from Donelan 1999.
              tht_ops=tht_wnd(i,j)+pi
              tpsw(i,j,m)=fsprm(beta,m,mm,tht_wnd(i,j),thetam(m))
     &             -0.4*fsprm(beta,m,mm,tht_ops     ,thetam(m))
              swin(i,j,m)=Cin *exp(-swcon *ageinv(i,j))*fsm(i,j)
              swin(i,j,m)=max(swin(i,j,m),10. )*ustp(i,j)**3*tpsw(i,j,m)
C Surface dissipation. bdis is used for wind driven portion of en.
C Reduced swfct*bdis is used for swell portion.
              bb=bdis                       
              if(abs(swin(i,j,m)).lt.1.e-6) bb=swfct*bdis  
              sdis(i,j,m)=adis*swin(i,j,m)+bb*sigth(i,j,m)*en(i,j,m)
              enf(i,j,m)=enf(i,j,m)+(swin(i,j,m)-sdis(i,j,m))*dtw2
C Bottom dissipation, first due to friction, second to depth-
C induced breaking a la Battjes and Janssen 1978.
             bdis1(i,j,m)=.003*(sigthb(i,j,m)
     &             *sqrt(2.*enb(i,j,m)/grav)*fsinhinv(kthD(i,j,m)))**3
             bdis2(i,j,m)=(en(i,j,m)/ent(i,j))*grav*sigth(i,j,m)
 
     &                /(8.*pi)*qdis(i,j)*(Hm(i,j))**2
             enf(i,j,m)=enf(i,j,m)-(bdis1(i,j,m)+bdis2(i,j,m))*dtw2
             enf(i,j,m)=max(.0001,enf(i,j,m))*fsm(i,j)
           enddo
c           endif
c!sc:testing limiting qdis 
C
         endif
        enddo
      enddo
       
      if(iint.eq.iprint) then
        write(6,'(''enb  ='',10e12.3)') (enb (12,6,m),m=1,26)
        write(6,'(''sigthb='',10e12.3)') (sigthb(12,6,m),m=1,26)
        write(6,'(''kthD  ='',10e12.3)') (kthD  (12,6,m),m=1,26)
        write(6,'(''bdis1='',10e12.3)') (bdis1(12,6,m),m=1,26)
        write(6,'(''bdis2='',10e12.3)') (bdis2(12,6,m),m=1,26)
        write(6,'(''qdis ='',10e12.3)') (qdis (i,6),i=1,24)
        write(6,'(''Hm   ='',10e12.3)') (Hm   (i,6),i=1,24)
      endif
C ---------------------------------------------------------------
C  Efficiency could be improved (slightly) for mode=0 or mode=2
C  calculations by adding code which would bypass using u and v.
      if(mode.eq.2) then
        do k=1,kb
          do j=1,jm
            do i=1,im
              u(i,j,k)=ua(i,j)
              v(i,j,k)=va(i,j)
            enddo
          enddo
        enddo
      endif
C  Include interactive radiation stress/current terms
       call cur2wave(im,jm,kb,dx,dy,fsm,z,dz,u,v,kthavD, 
     &       Sradx,Sradxy,Srady,Srad,ud,vd,iint)
      do i=3,im-2
        do j=3,jm-2
          do m=2,mm-1
C Eq.(14) has D instead of kth*D 
            sorc= kth(i,j,m)*dt(i,j)*en(i,j,m)
     &          * (  csm(m)*csm(m)*Sradx(i,j)
     &              +csm(m)*snm(m)*Sradxy(i,j)
     &              +snm(m)*snm(m)*Srady(i,j) 
     &              +              Srad(i,j)  )
            sorc=sorc*fsm(i-1,j)*fsm(i+1,j)*fsm(i,j-1)*fsm(i,j+1)
            enf(i,j,m)=enf(i,j,m)-sorc*dtw2
          enddo
        enddo
      enddo 
C ---------------------------------------------------------------
C                     x-y step for frequency
C ---------------------------------------------------------------
      do m=1,mm-1
        do j=2,jm-1
          do i=2,im-1
            cgxm=cg(i,j,m)*csm(m)+ud(i,j)
            cgym=cg(i,j,m)*snm(m)+vd(i,j)
            adv=(0.5*cgxm*(sigthb(i+1,j,m)-sigthb(i-1,j,m))
     &           -0.5*abs(cgxm)
     &            *(sigthb(i+1,j,m)-2.*sigthb(i,j,m)+sigthb(i-1,j,m)) )
     &                  /dx(i,j)
     &         +(0.5*cgym*(sigthb(i,j+1,m)-sigthb(i,j-1,m))
     &           -0.5*abs(cgym)
     &            *(sigthb(i,j+1,m)-2.*sigthb(i,j,m)+sigthb(i,j-1,m)) )
     &                  /dy(i,j)
            sigthf(i,j,m)=sigthb(i,j,m)-adv*fsm(i,j)*dtw2
!---- Nudge in wind driven sigp only when fspr > 0. --------------
            Cf= sqrt(fspr(i,j,m))*dtw2*sigthb(i,j,m)  
            sigthf(i,j,m)=sigthf(i,j,m)/(1.+Cf)+sigp(i,j)*Cf/(1.+Cf)   
          enddo
        enddo
      enddo
C ---------------------------------------------------------------
C                  Source step for frequency     
C ---------------------------------------------------------------
C N.B. Effect of time derivative of dt neglected for now.
      do i=2,im-1
        do j=2,jm-1
          dx2inv=1./(dx(i+1,j)+dx(i-1,j))
          dy2inv=1./(dy(i,j+1)+dy(i,j-1))
          do m=2,mm-1
            sor=-cg(i,j,m)*kth(i,j,m)*0.5*(
     &          (csm(m)*csm(m))*(ud(i+1,j)-ud(i-1,j))*dx2inv   
     &         +(csm(m)*snm(m))*(vd(i+1,j)-vd(i-1,j))*dx2inv    
     &         +(csm(m)*snm(m))*(ud(i,j+1)-ud(i,j-1))*dy2inv    
     &         +(snm(m)*snm(m))*(vd(i,j+1)-vd(i,j-1))*dy2inv   )
     &         +(sigth(i,j,m)/d(i,j))*(Fcg(i,j,m)-0.5)
     &           *( ud(i,j)*(d(i+1,j)-d(i-1,j))*dx2inv  
     &              +vd(i,j)*(d(i,j+1)-d(i,j-1))*dy2inv   )
           sor=sor*fsm(i-1,j)*fsm(i+1,j)*fsm(i,j-1)*fsm(i,j+1)
           sigthf(i,j,m)=sigthf(i,j,m)+sor*dtw2
           sigthf(i,j,m)=max(sigthf(i,j,m),0.1)
          enddo
        enddo
      enddo
      do j=2,jm-1
        do m=2,mm-1
          sigthf(1,j,m)=sigthf(2,j,m)
          sigthf(im,j,m)=sigthf(im-1,j,m)
        enddo
      enddo
      do i=2,im-1
        do m=2,mm-1
          sigthf(i,1,m)=sigthf(i,2,m)
          sigthf(i,jm,m)=sigthf(i,jm-1,m)
        enddo
      enddo
C--------------------------------------------------------------
C                  Boundary Condtions
C   Make sure that fsm=1 on open boundaries and
C   fsm=0 for closed boundaries.
C   For solid boundaries waves should pass through the boundaries
C   assuming they are thus dissipated. Otherwise, one of the
C   following 2 B.C.s can be implemeneted (uncommented) depending
C   on the desired physics.
C
C---------------------Reflection at side walls -----------
c     do i=1,im
c       do m=1,12
c         enf(i,2,m)=enf(i,3,26-m)
c         enf(i,jm-1,26-m)=enf(i,jm-2,m)
c       enddo
c     enddo
C--- For cyclic B.C.s on the north-south boundaries ----------
C--- to create j-independent problem -------------------------
c      do m=2,mm
c      do i=1,im
c        enf(i,2,m)=enf(i,jm-2,m)        !cyclic B.C.s
c        enf(i,jm-1,m)=enf(i,3,m)
c       enddo
c     enddo
C--- For open B.C.s on the eastern boundary ------------------
      do j=1,jm
         do m=1,mm
           cgnd=cgx(im,j,m)/(dx(im,j)+dx(im-1,j))*dtw
           enf(im,j,m)=en(im,j,m)
     &       -(cgnd+abs(cgnd))*(en(im,j,m)-en(im-1,j,m))
     &       -(cgnd-abs(cgnd))*(ene(j,m)-en(im,j,m))
         enddo
       enddo
C--- For open B.C.s on the western boundary -------------------
      do j=1,jm
         do m=1,mm
           cgnd=cgx(2,j,m)/(dx(2,j)+dx(1,j))*dtw
           enf(1,j,m)=en(1,j,m)
     &       -(cgnd+abs(cgnd))*(en(1,j,m)-enw(j,m))
     &       -(cgnd-abs(cgnd))*(en(2,j,m)-en(1,j,m))
         enddo
       enddo
C--- For open B.C.s on the southern boundary -------------------
      do i=1,im 
         do m=1,mm
           cgnd=cgy(i,2,m)/(dy(i,2)+dy(i,1))*dtw
           enf(i,1,m)=en(i,1,m)
     &       -(cgnd+abs(cgnd))*(en(i,1,m)-ens(i,m))
     &       -(cgnd-abs(cgnd))*(en(i,2,m)-en(i,1,m))
         enddo
       enddo
C------------------------------------------------------------
C      Mask out land values and apply Asselin filter          
C------------------------------------------------------------
      do m=1,mm
        do j=1,jm
          do i=1,im
            enf(i,j,m)=enf(i,j,m)*fsm(i,j)
            en(i,j,m)=en(i,j,m)
     &        +0.5*smoth*(enb(i,j,m)-2.*en(i,j,m)+enf(i,j,m))
            sigth(i,j,m)=sigth(i,j,m)
     &        +0.5*smoth*(sigthb(i,j,m)-2.*sigth(i,j,m)+sigthf(i,j,m))
             en(i,j,m)=max(.0001,en(i,j,m))*fsm(i,j)
             enf(i,j,m)=max(.0001,enf(i,j,m))*fsm(i,j)
          enddo
        enddo
      enddo
C------ Reset time sequence -----------------------------------
      do m=1,mm
        do j=1,jm
          do i=1,im
            enb(i,j,m)=en(i,j,m) 
            en(i,j,m)=enf(i,j,m)
            sigthb(i,j,m)=sigth(i,j,m) 
            sigth(i,j,m)=sigthf(i,j,m)
          enddo
        enddo
      enddo
c ----Calculate total dissipation for use in the tke equation ---
      do j=1,jm
        do i=1,im
          sdistot(i,j)=summ(sdis,dth,i,j,im,jm,mm)*fsm(i,j)  
          do m=1,mm
            tpsw(i,j,m)=bdis1(i,j,m)+bdis2(i,j,m)
          enddo
          bdistot(i,j)=summ(tpsw,dth,i,j,im,jm,mm)*fsm(i,j) 
        enddo
      enddo
C------------------------------------------------------------
C     Determine some diagnostic wave properties
C------------------------------------------------------------
      do j=1,jm
        do i=1,im
                        if(fsm(i,j).eq.1.) then
C  Determine average wave propagation angle
          do m=1,mm
            tpsw(i,j,m)=en(i,j,m)*csm(m)  
          enddo
          csmav=summ(tpsw,dth,i,j,im,jm,mm)
          do m=1,mm
            tpsw(i,j,m)=en(i,j,m)*snm(m)  
          enddo
          snmav=summ(tpsw,dth,i,j,im,jm,mm)
          if(csmav.eq.0.) csmav=csmav+1.e-7
          thtav(i,j)=atan(snmav/csmav)
          if (snmav.gt.0..and.csmav.lt.0.) thtav(i,j)=thtav(i,j)+pi
          if (snmav.lt.0..and.csmav.lt.0.) thtav(i,j)=thtav(i,j)-pi
          swintot(i,j)=summ(swin,dth,i,j,im,jm,mm)
C Determine average frequency and wave number 
          do m=1,mm
            tpsw(i,j,m)=en(i,j,m)*sigth(i,j,m)
          enddo
          ent(i,j)=summ(en,dth,i,j,im,jm,mm)
          sigthav(i,j)=summ(tpsw ,dth,i,j,im,jm,mm)/ent(i,j)
          do m=1,mm
            tpsw(i,j,m)=en(i,j,m)*kth(i,j,m)
          enddo
          kthav(i,j)=summ(tpsw ,dth,i,j,im,jm,mm)/ent(i,j)
          kthavD(i,j)=kthav(i,j)*dt(i,j)
                            endif
        enddo
      enddo
      return
      end
C--------------------------------------------------------------
       subroutine botdrag (im,jm,fsm,
     &     kp,cp,ent,wubot,wvbot,d,zzkbm1,z0b,cbcmin,cbcmax,cbc)
C--------------------------------------------------------------
C Provides wave influenced bottom fricion. Derived from
C Nielson, P., 1992, Coastal bottom boundary layers and sediment
C transport. World Scientific.
C This aspect needs more basic research.
C 
C  Input:
C    fsm(I,j) = 1 for water cell, = 0 for land cell
C    cp(i,j) = peak frequency phase speed (m s-1)
C    kp(i,j) = peak frequency wave number
C    ent(i,j) = wave energy averaged over all wave angles
C    wubot(i,j),wvbot(i,j) = bottom stress vector as determned
C      by bottom current (at zz(kb-1)) and cbc found here
C    zzkbm1 = zz(kb-1) (non-d.)
C    d(i,j) = water column depth (m)
C    zob = roughness parameter (m)
C    cbcmin,cbcmax = limits on cbs (non-d.)    
C
C  Output:
C    cbc(i,j) = bottom drag coefficient (non-d.)
C--------------------------------------------------------------
      implicit none
      integer i,j,im,jm
      real kappa,grav,fsm(im,jm),d(im,jm)
      real z0b,z0a,zzkbm1,cbcmin,cbcmax,cbc(im,jm)
      real uboscil,utau2,fsinhinv
      real kp(im,jm),cp(im,jm),ent(im,jm),
     &     wubot(im,jm),wvbot(im,jm)
      data kappa/0.4/,grav/9.807/
      do j=1,jm
        do i=1,im
          if (fsm(i,j).eq.1.) then
c!sc:fsinhinv(kp(i,j)*d(i,j)) was 0
            uboscil=cp(i,j)*kp(i,j)*sqrt(2.*ent(i,j)/grav)
     &               /(fsinhinv(kp(i,j)*d(i,j))+1.e-12)
c!sc:prevent utau2 become 0
            utau2=sqrt(wubot(i,j)**2+wvbot(i,j)**2+1.e-12)
c            utau2=sqrt(wubot(i,j)**2+wvbot(i,j)**2)
            z0a=z0b*(1.+0.05*uboscil**2/utau2)
            cbc(i,j)=(kappa/(log(1.+(1.0+zzkbm1)*d(i,j)/z0a)+1.e-12))**2
            cbc(i,j)=max(cbcmin,cbc(i,j))
C
C     If the following is invoked, then it is probable that the wrong
C     choice of z0b or vertical spacing has been made:
C
            cbc(i,j)=min(cbcmax,cbc(i,j))
C!sc:alternative cbc setting
c!sc:          if ((wubot(i,j).le.1.e-12).and.(wvbot(i,j).le.1.e-12))then
c!sc:            cbc(i,j)=cbcmin
c!sc:          endif
          endif
        enddo
      enddo
      return
      end

C--------------------------------------------------------------------
      subroutine cur2wave(im,jm,kb,dx,dy,fsm,z,dz,u,v,kthavD,
     &      Sradx,Sradxy,Srady,Srad,ud,vd,iint)
C--------------------------------------------------------------------
C                  Current to wave interaction terms
C     Calculates vertical averages using spectrally averaged functions
C     provided by subroutine specavs.
C     Terms described by Eq. (14) in Mellor et al (2008)
C--------------------------------------------------------------------
      implicit none
C   Input       
      integer i,j,k,im,jm,kb,iint
      real kthavD(im,jm),fsm(im,jm) 
      real dx(im,jm),dy(im,jm),z(kb),dz(kb),u(im,jm,kb),v(im,jm,kb)
C   Output
      real  Sradx(im,jm),Sradxy(im,jm),Srady(im,jm),Srad(im,jm),
     &      Srad0(im,jm),ud(im,jm),vd(im,jm)
C   Internal
      real sumk 
      real F1(kb),F2(kb),F3(kb),F4(kb),
     &    tp1(kb),tp2(kb),tp3(kb),
     &    gradxu(kb),gradxv(kb),gradyu(kb),gradyv(kb),
     &    dxiv,dyiv
      do j=3,jm-2
        do i=3,im-2  
          call specavs (kthavD(i,j),F1,F2,F3,F4,z,kb,i,j,iint)
          dxiv=1./dx(i,j)   
          dyiv=1./dy(i,j)   
          do k=1,kb-1
            tp1(k)=0.25*(u(i,j,k)+u(i+1,j,k))*(F3(k)+F3(k+1))*dz(k)
            tp2(k)=0.25*(v(i,j,k)+v(i,j+1,k))*(F3(k)+F3(k+1))*dz(k)
            gradxu(k)=(u(i+1,j,k)-u(i,j,k))*dxiv*fsm(i+1,j)*fsm(i,j)    
            gradxv(k)=0.25*(v(i+1,j+1,k)-v(i-1,j,k)
     &                +v(i+1,j,k)-v(i-1,j+1,k))*dxiv  ! corrected 07/07/09
     &               *fsm(i+1,j+1)*fsm(i-1,j)*fsm(i+1,j)*fsm(i-1,j+1)
            gradyu(k)=0.25*(u(i+1,j+1,k)-u(i,j-1,k)
     &                +u(i,j+1,k)-u(i+1,j-1,k))*dyiv  ! corrected 07/07/09
     &               *fsm(i+1,j+1)*fsm(i,j-1)*fsm(i,j+1)*fsm(i+1,j-1)
            gradyv(k)=(v(i,j+1,k)-v(i,j,k))*dyiv*fsm(i,j+1)*fsm(i,j)  
          enddo
          ud(i,j)=sumk(tp1,1.,kb)    ! Dopler velocity
          vd(i,j)=sumk(tp2,1.,kb)    !   ''     ''
C
          do k=1,kb-1
            tp1(k)=0.5*(F1(k)+F1(k+1))*gradxu(k)*dz(k)
            tp2(k)=0.5*(F1(k)+F1(k+1))*(gradxv(k)+gradyu(k))*dz(k)
            tp3(k)=0.5*(F1(k)+F1(k+1))*gradyv(k)*dz(k)
          enddo
          Sradx(i,j)=sumk(tp1,1.,kb)
          Sradxy(i,j)=sumk(tp2,1.,kb)
          Srady(i,j)=sumk(tp3,1.,kb)
          do k=1,kb-1
            tp1(k)=0.5*(F2(k)+F2(k+1))*gradxu(k)*dz(k)
            tp2(k)=0.5*(F2(k)+F2(k+1))*gradyv(k)*dz(k)
          enddo
          Srad(i,j)=-(sumk(tp1,1.,kb)+sumk(tp2,1.,kb))
        enddo
      enddo
      do j=1,jm
        ud(im-1,j)=ud(im-2,j)
        ud(im,j)=ud(im-1,j)
        ud(2,j)=ud(3,j)
        ud(1,j)=ud(2,j)
        vd(im-1,j)=vd(im-2,j)
        vd(im,j)=vd(im-1,j)
        vd(2,j)=vd(3,j)
        vd(1,j)=vd(2,j)
      enddo
      do i=1,im
        ud(i,jm-1)=ud(i,jm-2)
        ud(i,jm)=ud(i,jm-1)
        ud(i,2)=ud(i,3)
        ud(i,1)=ud(i,2)
        vd(i,jm-1)=vd(i,jm-2)
        vd(i,jm)=vd(i,jm-1)
        vd(i,2)=vd(i,3)
        vd(i,1)=vd(i,2)
      enddo
      return
      end

C----------------------------------------------------------------
      subroutine specavs (kpD,F1,F2,F3,F4,z,kb,i,j,iint)
C       This is a subroutine to deliver spectrally averaged F1 - F4.
C       averaged. 
C       F1-F4 are spectrally averaged components of wave radiation
C       stresses. 
C       N.B. This subroutine uses the file, specavs.
      implicit none
      integer me,ne
      parameter (me=17,ne=21)
      integer i,j,k,m,n,kb,init,iint
      integer m1,n1
      real kpD,kpDd(me)
      real F1d(me,ne),F2d(me,ne),F3d(me,ne),F4d(me,ne) 
      real F1i(ne),F2i(ne),F3i(ne),F4i(ne)
      real F1(kb),F2(kb),F3(kb),F4(kb)
      real zeta(ne),zetf(ne),z(kb) 
      save init,kpDd,F1d,F2d,F3d,F4d,zeta,zetf
      data init/0/
C----------------------------------------------------------------
                       if(init.eq.0) then
C----------------------------------------------------------------
        do m=1,16
          kpDd(m)=0.2*float(m-1)
        enddo
        write(6,'(''Read spectrally averaged functions from specavs'')')
        read(10,'(2x)')
        write(6,'(/,''   F1'')')
        do n=1,ne 
          read(10,'(f6.2,15(f6.3),2f8.3)')
     &              zeta(n),(F1d(m,n),m=2,me),zetf(n)
          F1d(1,n)=1.0/0.1    ! avoid infinity at kpD=0
          write(6,'(f6.2,f6.2,15(f6.3),2f8.3)')
     &              zeta(n),(F1d(m,n),m=1,me),zetf(n)
        enddo
        read(10,'(2x)')
        write(6,'(/,''   F2'')')
        do n=1,ne 
          read(10,'(f6.2,15(f6.3),2f8.3)')
     &              zeta(n),(F2d(m,n),m=2,me),zetf(n)
          F2d(1,n)=0.               ! asymptote for kpD=0
         write(6,'(f6.2,16(f6.3),2f8.3)')
     &              zeta(n),(F2d(m,n),m=1,me),zetf(n)
        enddo
        read(10,'(2x)')
        write(6,'(/,''   F3'')')
        do n=1,ne 
          read(10,'(f6.2,15(f6.3),2f8.3)')
     &              zeta(n),(F3d(m,n),m=2,me),zetf(n)
          F3d(1,n)=1.5+zeta(n)    
         write(6,'(f6.2,16(f6.3),2f8.3)')
     &              zeta(n),(F3d(m,n),m=1,me),zetf(n)
        enddo
        read(10,'(2x)')
        write(6,'(/,''   F4'')')
        do n=1,ne 
          read(10,'(f6.2,15(f6.3),2f8.3)')
     &              zeta(n),(F4d(m,n),m=2,me),zetf(n)
          F4d(1,n)=1.+zeta(n)
         write(6,'(f6.2,16(f6.3),2f8.3)')
     &              zeta(n),(F4d(m,n),m=1,me),zetf(n)

       enddo 
C----------------------------------------------------------------
                          init=1 

                          endif  
               
C----------------------------------------------------------------
C  Solve for F1, F2, F3, F4
C----------------------------------------------------------------
      m=kpD/0.2+1
c      if (m.lt.16) then
c!sc:m can be less than 1
      if (m.lt.16.and.m.gt.0) then
C Interpolate w.r.t. kpD =< 3. (0 < m< 17)
        do n=1,ne
          F1i(n)=F1d(m,n)+(F1d(m+1,n)-F1d(m,n))*(kpD-kpDd(m))/0.2
          F2i(n)=F2d(m,n)+(F2d(m+1,n)-F2d(m,n))*(kpD-kpDd(m))/0.2
          F3i(n)=F3d(m,n)+(F3d(m+1,n)-F3d(m,n))*(kpD-kpDd(m))/0.2
          F4i(n)=F4d(m,n)+(F4d(m+1,n)-F4d(m,n))*(kpD-kpDd(m))/0.2
        enddo
C Interpolate w.r.t. z
        do k=1,kb-1
          n=-z(k)/.05+1
          F1(k)=F1i(n)+(F1i(n+1)-F1i(n))*(-z(k)+zeta(n))/.05
          F2(k)=F2i(n)+(F2i(n+1)-F2i(n))*(-z(k)+zeta(n))/.05
          F3(k)=F3i(n)+(F3i(n+1)-F3i(n))*(-z(k)+zeta(n))/.05
          F4(k)=F4i(n)+(F4i(n+1)-F4i(n))*(-z(k)+zeta(n))/.05
        enddo
        F1(kb)=F1i(ne)
        F2(kb)=F2i(ne)
        F3(kb)=F3i(ne)
        F4(kb)=F4i(ne)
      else
C For m=17, kpD=100 and F1 is scaled different from m<17
C Use similarity for large kpD > 3. (m = 17)
C e.g., F1(kpD, z) = (kpD/100)*F1(100, kpD*z/100)
        do k=1,kb
          F1(k)=0.
          F2(k)=0.
          F3(k)=0.
          F4(k)=0.
          n=1-z(k)*kpD/.1
c          if(n.lt.ne) then     !asumptotes for kpD > 3 
c!sc:n can be less than 1
          if(n.lt.ne.and.n.gt.0) then
            F1(k)=F1d(17,n)+(F1d(17,n+1)-F1d(17,n))/.001
     &              *(-z(k)*kpD/100.+zetf(n))
            F2(k)=F1(k)
            F3(k)=F3d(17,n)+(F3d(17,n+1)-F3d(17,n))/.001
     &              *(-z(k)*kpD/100.+zetf(n))
            F3(k)=F3(k)*kpD*0.01
            F4(k)=F4d(17,n)+(F4d(17,n+1)-F4d(17,n))/.001
     &              *(-z(k)*kpD/100.+zetf(n))
          else if(n.le.0) then
            write(*,*)'xxx_wave_specavs',n,z(k),kpD
          endif  
        enddo
      endif
      return  
      end

C----------------------------------------------------------------
       subroutine speeds(kpD,Fcg,Fct)
C       This is a subroutine to deliver spectrally averaged Fcg, Fct
C       Fcg = ratio of cg to cp where cp is spectral peak phase speed.
C       Fct = is factor in cthth (theta speed) to make it a spectrally
C       averaged. Fcg and Fct are functions of kpD.
       implicit none
       integer k
       real kpD,kpDd(17),Fcg,Fct
       real Fcgd(17),Fctd(17),c_cp(17)
c!sc:dum may conflict with common variable 
c, dum(17)
       data kpDd/
     &    0.0,   0.2,   0.4,   0.6,   0.8,   1.0,   1.2,   1.4,   1.6,
     &    1.8,   2.0,   2.2,   2.4,   2.6,   2.8,   3.0, 100.0/
       data c_cp/
     &  1.000, 1.005, 1.001, 0.968, 0.951, 0.936, 0.924, 0.915, 0.910,
     &  0.906, 0.905, 0.904, 0.904, 0.904, 0.905, 0.905, 0.910/
       data Fcgd/
     &  1.000, 0.984, 0.928, 0.834, 0.756, 0.689, 0.635, 0.592, 0.560,
     &  0.535, 0.516, 0.502, 0.491, 0.483, 0.477, 0.472, 0.472/
c!sc:dum may conflict with common variable 
c       data dum/
c     &  1.000, 0.987, 0.950, 0.897, 0.837, 0.776, 0.720, 0.671, 0.631,
c     &  0.598, 0.573, 0.554, 0.539, 0.529, 0.521, 0.515, 0.500/
       data Fctd/
     &  1.000, 0.984, 0.954, 0.931, 0.892, 0.858, 0.830, 0.809, 0.795,
     &  0.787, 0.784, 0.786, 0.791, 0.801, 0.814, 0.829, 0.829/

        k=kpD/0.2+1
c!sc: too many output
c        if(kpD.lt.1) call printall
c!sc: k anomaly
c      if (k.lt.16) then
      if (k.lt.16.and.k.gt.0) then
        Fcg=Fcgd(k)+(Fcgd(k+1)-Fcgd(k))*(kpD-kpDd(k))/0.2
        Fct=Fctd(k)+(Fctd(k+1)-Fctd(k))*(kpD-kpDd(k))/0.2
      else
        if(k.le.0) then
          write(*,*)'xxx_wave_speeds',k,kpD
        endif
        k=16
        Fcg=Fcgd(k)+(Fcgd(k+1)-Fcgd(k))*(kpD-kpDd(k))/97.
        Fct=Fctd(k)+(Fctd(k+1)-Fctd(k))*(kpD-kpDd(k))/97.
      endif
      return
      end

C------------------------------------------------------------------
      subroutine stress (iint,nwave,Hs,ageinv,usdif,vsdif,cd,tpx0,tpy0,
     &                     cdmax,ttx0,tty0,ustw,ustp,im,jm,fsm)
C-----------------------------------------------------------------
C  Input: usdif,vsdif = difference between 10m wind velocity and
C                       surface velocity
C         Hs = significant wave height
C         ageinv = inverse wave age
C  Output: tpx0,tpy0 = wind pressure wave slope correlation =
C          wind induced surface stress
C          ttx0,tty0 = turbulence surface stress (non-zero only
C          for low winds (< about 4 m s-1)
C              Surface roughness due to Donelan (1990)
C-----------------------------------------------------------------        
C      implicit none
      integer i,j,im,jm
      real u1,u2,z0w,z0t,cdt,cdw,grav,pi,r,rhalf 
      real Hs(im,jm),ageinv(im,jm),usdif(im,jm),vsdif(im,jm)  
      real cd(im,jm),tpx0(im,jm),tpy0(im,jm),ustw(im,jm),ustp(im,jm) 
      real ttx0(im,jm),tty0(im,jm),fsm(im,jm)
      real kappa,nu
      data nu/1.8e-6/,r/0.0011630/,rhalf/.0034100/,kappa/0.41/
      data z0w/1.e-6/,z0t/1.e-6/,grav/9.807/,pi/3.1415926/
      do j=2,jm-1
        do i=2,im-1
          if(fsm(i,j).gt.0.) then
c          write(*,*)'xxx1',usdif(i,j),vsdif(i,j),i,j
          u2=usdif(i,j)**2+vsdif(i,j)**2
          u1=sqrt(u2)
          z0w=1.38e-4*Hs(i,j)*(ageinv(i,j))**2.66+1.e-5   !Donelan
          z0t=0.18*nu/(ustw(i,j) +.000001)
          cdw=(kappa/log(10./z0w))**2
          cdw=min(cdw,cdmax)
          cdt=(kappa/log(10./z0t))**2
C The transition from smooth surface turbulent flow and friction
C drag to a wave surface and form drag is abrupt; this probably
C needs improvement but nevertheless gives a continuos Cd vs U10
C result that agrees with data.
          if(cdw.gt.cdt) then
            tpx0(i,j)=r*cdw*u1*usdif(i,j)
            tpy0(i,j)=r*cdw*u1*vsdif(i,j)
            ustp(i,j)=sqrt(sqrt(tpx0(i,j)**2+tpy0(i,j)**2))
            ttx0(i,j)=0.
            tty0(i,j)=0.
            cdt=0.      
          else
            ttx0(i,j)=r*cdt*u1*usdif(i,j)
            tty0(i,j)=r*cdt*u1*vsdif(i,j)
            tpx0(i,j)=0.
            tpy0(i,j)=0.
            cdw=0.
          endif
          cd(i,j)=cdw+cdt
          ustw(i,j)=sqrt(r*(cdw+cdt)*u2)
          endif
        enddo
      enddo
      return
      end

C----------------------------------------------------------------
      subroutine stepmpd(fb,f,ff,c,im,jm,mm,dtw2,dthinv,nitera)
C**********************************************************************
C *                                                                    *
C * FUNCTION    :  Integrates conservative scalar equations.           *
C *                                                                    *
C *                This is a first-order upstream scheme, which        *
C *                reduces implicit diffusion using the Smolarkiewicz  *
C *                iterative upstream scheme with an antidiffusive     *
C *                velocity.                                           *
C *                                                                    *
C *                It is based on the subroutines of Gianmaria Sannino *
C *                (Inter-university Computing Consortium, Rome, Italy)*
C *                and Vincenzo Artale (Italian National Agency for    *
C *                New Technology and Environment, Rome, Italy),       *
C *                downloaded from the POM FTP site on 1 Nov. 2001.    *
C *                                                                    *
C *                                                                    *
C *                Reference:                                          *
C *                                                                    *
C *                Smolarkiewicz, P.K.; A fully multidimensional       *
C *                  positive definite advection transport algorithm   *
C *                  with small implicit diffusion, Journal of         *
C *                  Computational Physics, 54, 325-362, 1984.         *
C *                                                                    *
C **********************************************************************
C
      implicit none
      integer itera,nitera,i,j,m,im,jm,mm
      real sw,dtw2,dthinv,value_min,epsilon
      real fbmem(im,jm,mm),ffmem(im,jm,mm),c(im,jm,mm)
      real flux(im,jm,mm),speed(im,jm,mm)
      real fb(im,jm,mm),f(im,jm,mm),ff(im,jm,mm)
      real mol,abs_1,abs_2
      real udx,u2dt,vdy,v2dt,wdz,w2dt 
      value_min=1.e-9; epsilon=1.e-14; sw=1.0
c!sc:testing larger value min
c      value_min=1.e-3; epsilon=1.e-14; sw=1.0
C
      do m=1,mm
        do j=1,jm
          do i=1,im
            speed(i,j,m)=c(i,j,m)
            fbmem(i,j,m)=fb(i,j,m)
            ffmem(i,j,m)=ff(i,j,m)
          end do
        end do
      end do
c      write(*,*)'xxx3a',speed(141,178,2),dtw2,dthinv,c(141,178,2)
C
      do itera=1,nitera
C
        do m=2,mm-1
          do j=1,jm
            do i=1,im
              flux(i,j,m)=0.5e0
     $                      *((speed(i,j,m)+abs(speed(i,j,m)))
     $                       *fbmem(i,j,m-1)+
     $                        (speed(i,j,m)-abs(speed(i,j,m)))
     $                       *fbmem(i,j,m))
            end do
          end do
        end do
        do j=1,jm
          do i=1,im
            flux(i,j,1)=flux(i,j,mm-1)
            flux(i,j,mm)=flux(i,j,2)
          end do
        end do
C
C     Add net advective fluxes and step forward in time:
C
        do m=1,mm-1
          do j=1,jm
            do i=1,im
              ff(i,j,m)=ffmem(i,j,m)
     $             +(flux(i,j,m)-flux(i,j,m+1))*dthinv*dtw2
            end do
          end do
        end do
C *                Calculate the antidiffusive speed used to           *
C *                reduce the numerical diffusion associated with the  *
C *                upstream differencing scheme.                       *
C *                                                                    *
        do m=2,mm-1
          do j=1,jm
            do i=1,im
            if(ff(i,j,m).lt.value_min.or.
     $         ff(i,j,m-1).lt.value_min) then
              speed(i,j,m)=0.e0
            else
              wdz=abs(speed(i,j,m))
c      write(*,*)'xxx3',speed(i,j,m),dtw2,dthinv,i,j,m,c(i,j,m)
              w2dt=speed(i,j,m)*speed(i,j,m)*dtw2*dthinv              
              mol=(ff(i,j,m)-ff(i,j,m-1))
     $             /(ff(i,j,m)+ff(i,j,m-1)+epsilon)
              speed(i,j,m)=(wdz-w2dt)*mol*sw
              abs_1=abs(wdz)
              abs_2=abs(w2dt)
              if(abs_1.lt.abs_2)speed(i,j,m)=0.e0
            end if
          end do
        end do
      end do
C
        do m=1,mm
          do j=1,jm
            do i=1,im
              ffmem(i,j,m)=ff(i,j,m)
              fbmem(i,j,m)=ff(i,j,m)
            end do
          end do
        end do
C
C     End of Smolarkiewicz scheme
      end do
      return
      end

C-----------------------------------------------------------------
      subroutine wave2cur  (dth)
C  This routine supplies Sxx,Sxy.Syy for insertion into the
C  into the momentum equation when waves are coupled with 
C  the circulation model. 

      implicit none
      include 'pom08.c'  
      real tpxx(mm),tpxy(mm),tpyy(mm),tp(mm),summ1
      real stpxx,stpxy,stpyy,stp,dxinv,dyinv
      real F1m(kb),F2m(kb),us(im,jm,kb)
C        if(iint.eq.3.and.nwave.eq.2) stop
      do j=2,jm-1
        do i=2,im-1  
          if(fsm(i,j).gt.0.) then
          do m=1,mm
            tpxx(m)=kth(i,j,m)*en(i,j,m)*csm(m)*csm(m)         
            tpxy(m)=kth(i,j,m)*en(i,j,m)*csm(m)*snm(m)
            tpyy(m)=kth(i,j,m)*en(i,j,m)*snm(m)*snm(m)
            tp(m)=kth(i,j,m)*en(i,j,m)
          enddo
            stpxx=summ1(tpxx,dth,mm)
            stpxy=summ1(tpxy,dth,mm)
            stpyy=summ1(tpyy,dth,mm)
            stp  =summ1(tp  ,dth,mm)
            call specavs (kthavD(i,j),F1,F2,F3,F4,z,kb,i,j,iint)
          do k=1,kb-1
            F1m(k)=0.5*(F1(k)+F1(k+1))
            F2m(k)=0.5*(F2(k)+F2(k+1))
            Sxx(i,j,k)=stpxx*F1m(k)+stp*F2m(k)
            Sxy(i,j,k)=stpxy*F1m(k)
            Syy(i,j,k)=stpyy*F1m(k)+stp*F2m(k)
            tpzdist(i,j,k)=F4(k)     ! subsurface momentum
                                     ! tranfer after multiplication
                                     ! by (tpx0,tpy0) after do 9000
C Calculate Stokes drift
c!sc:denominator term have to changce to be 0 so 1.e-12 is added
          us(i,j,k)=ent(i,j)*kthav(i,j)/(sigthav(i,j)*dt(i,j)+1.e-12)
     &        *(F4(k)-F4(k+1))/dz(k)
          ust(i,j,k)=us(i,j,k)*cos(thtav(i,j))
          vst(i,j,k)=us(i,j,k)*sin(thtav(i,j))
          enddo

c         if(i.eq.7.and.j.eq.6.and.iint.eq.180) then
c           do k=1,kb
c           write(6,'(''wave2'',i5,4e12.3)')k,kthav(i,j),kthavD(i,j),
c    &         F1(k),F2(k)
c           write(6,'(''wave2'',i5,6e12.3)')
c    $         k,F1m(k),F2m(k),stpxx,en(i,j,13),ent(i,j),Sxx(i,j,k)
c           enddo
c         endif
        endif
        enddo
      enddo
      return
      end

C --------------------------------------------------------------
      subroutine wvfcns(sig,d,wvnu,c)   
      implicit none
      real geeinv
      parameter(geeinv=1./9.807)
C
C This lookup table inputs sig and d and outputs wvnu, c and cg  
C using the dispersion relation
C
C sig = frequency (s-1)
C d = water column depth (m) 
C wvnu = wave number (m-1)
C c = phase speed (m s-1)
C cg = group speed (m s-1)  
C f1 = sig**2*d/gee
C f2 = wvnu*d
C f3 = cg/c
C
C ------------------------------------------------------------------
      integer k
      real sig,d,wvnu,c,cg,dknd,cgnd,fqnd
      real f1(81),f2(81),f3(81)
      data f1/ 
     &  0.0000,  0.0500,  0.1000,  0.1500,  0.2000,  0.2500,  0.3000,
     &  0.3500,  0.4000,  0.4500,  0.5000,  0.5500,  0.6000,  0.6500,
     &  0.7000,  0.7500,  0.8000,  0.8500,  0.9000,  0.9500,  1.0000,
     &  1.0500,  1.1000,  1.1500,  1.2000,  1.2500,  1.3000,  1.3500,
     &  1.4000,  1.4500,  1.5000,  1.5500,  1.6000,  1.6500,  1.7000,
     &  1.7500,  1.8000,  1.8500,  1.9000,  1.9500,  2.0000,  2.0500,
     &  2.1000,  2.1500,  2.2000,  2.2500,  2.3000,  2.3500,  2.4000,
     &  2.4500,  2.5000,  2.5500,  2.6000,  2.6500,  2.7000,  2.7500,
     &  2.8000,  2.8500,  2.9000,  2.9500,  3.0000,  3.0500,  3.1000,
     &  3.1500,  3.2000,  3.2500,  3.3000,  3.3500,  3.4000,  3.4500,
     &  3.5000,  3.5500,  3.6000,  3.6500,  3.7000,  3.7500,  3.8000,
     &  3.8500,  3.9000,  3.9500,  4.0000/
      data f2/
     &  0.0000,  0.2255,  0.3216,  0.3973,  0.4627,  0.5218,  0.5767,
     &  0.6284,  0.6778,  0.7255,  0.7717,  0.8168,  0.8611,  0.9046,
     &  0.9476,  0.9902,  1.0324,  1.0744,  1.1163,  1.1580,  1.1997,
     &  1.2414,  1.2831,  1.3249,  1.3668,  1.4088,  1.4511,  1.4934,
     &  1.5360,  1.5788,  1.6218,  1.6651,  1.7085,  1.7523,  1.7962,
     &  1.8405,  1.8849,  1.9297,  1.9746,  2.0199,  2.0653,  2.1110,
     &  2.1569,  2.2031,  2.2495,  2.2960,  2.3428,  2.3898,  2.4369,
     &  2.4843,  2.5318,  2.5795,  2.6273,  2.6752,  2.7233,  2.7716,
     &  2.8199,  2.8684,  2.9170,  2.9657,  3.0144,  3.0633,  3.1123,
     &  3.1613,  3.2104,  3.2596,  3.3088,  3.3581,  3.4074,  3.4568,
     &  3.5063,  3.5558,  3.6053,  3.6548,  3.7044,  3.7541,  3.8037,
     &  3.8500,  3.9000,  3.9500,  4.0000/
      data f3/
     &  1.0000,  0.9834,  0.9671,  0.9510,  0.9352,  0.9196,  0.9042,
     &  0.8891,  0.8743,  0.8598,  0.8455,  0.8315,  0.8179,  0.8045,
     &  0.7914,  0.7786,  0.7662,  0.7541,  0.7422,  0.7308,  0.7196,
     &  0.7088,  0.6983,  0.6882,  0.6784,  0.6689,  0.6598,  0.6511,
     &  0.6426,  0.6345,  0.6268,  0.6193,  0.6122,  0.6054,  0.5990,
     &  0.5928,  0.5870,  0.5814,  0.5761,  0.5711,  0.5664,  0.5619,
     &  0.5577,  0.5538,  0.5500,  0.5465,  0.5432,  0.5401,  0.5373,
     &  0.5345,  0.5320,  0.5297,  0.5274,  0.5254,  0.5235,  0.5217,
     &  0.5200,  0.5185,  0.5171,  0.5157,  0.5145,  0.5134,  0.5123,
     &  0.5114,  0.5104,  0.5096,  0.5088,  0.5081,  0.5075,  0.5069,
     &  0.5063,  0.5058,  0.5053,  0.5049,  0.5045,  0.5041,  0.5038,
     &  0.5000,  0.5000,  0.5000,  0.5000/
      fqnd=sig*sig*d*geeinv
      k=min(1+fqnd*20.,80.)  
        dknd=f2(k)+(f2(k+1)-f2(k))*(fqnd-f1(k))*20.
        cgnd=f3(k)+(f3(k+1)-f3(k))*(fqnd-f1(k))*20.
        if(fqnd.lt.0.2) dknd=sqrt(fqnd)
      wvnu=dknd/d
      c=sig/wvnu
      cg=cgnd*c
      return
      end
C--------------------------------------------------------------------
      function fsprm(beta,m,mm,tht_wnd,thetam)
      save tht,fdat,dtht,ifirst
      dimension tht(35),fdat(35) 
      data pi/3.141592654/
      data ifirst/0/
      if(ifirst.eq.0) then
        write(6,'(''spreading function'',/,
     &            ''    m     tht       fdat'')')
       dtht=pi/float(32)
        do k=1,35 
          fdat(k)=0.
          tht(k)=dtht*(k-1)
          if(k.lt.17) fdat(k)=0.5*beta/cosh(beta*tht(k))**2
          write(6,'(i5,2f10.4)') k,tht(k),fdat(k)
        enddo
      endif
      ifirst=1
C Create spreading function around the wind direction
      fsprm=0.
      thtrel=abs(thetam-tht_wnd)
      if(thtrel.gt.pi) thtrel=abs(thtrel-2.*pi)  
      kth=thtrel/dtht+1  
c     if(kth.le.14)
        fsprm=fdat(kth)
     &    +(fdat(kth+1)-fdat(kth))*(thtrel-tht(kth))/dtht   
      return
      end
     
C------------------------------------------------------------
      function fcoshinv(x)
      dimension xd(15),fdat(15)
      data xd/0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7./
      data fdat/1.0,0.887,0.648,0.425,0.265,0.163,0.099,0.060,
     &         0.037,0.022,0.013,0.008,0.005,0.003,0.002/
      fcoshinv=0.0 
      if(x.lt.7.0) then
        ixd=int(2.0*x)+1
        fcoshinv=fdat(ixd)+(fdat(ixd+1)-fdat(ixd))*(x-xd(ixd))*2.0
       endif
       return
       end

C----------------------------------------------------------------
      function fsinhinv(x)
      dimension xd(15),fdat(15)
      data xd/0.,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7./
      data fdat/1.919,1.919,0.851,0.470,0.276,0.165,0.100,0.060,
     &         0.037,0.022,0.013,0.008,0.005,0.003,0.002/
      fsinhinv=0.0 
      if(x.lt.7.0) then
        ixd=int(2.0*x)+1
        fsinhinv=fdat(ixd)+(fdat(ixd+1)-fdat(ixd))*(x-xd(ixd))*2.0
       endif
       return
       end
       
C----------------------------------------------------------------
      function mploc(a,b,i,j,im,jm,mm)
      dimension a(im,jm,mm),b(im,jm,mm)
      amx=0.
      mmx=1
      do m=2,mm
        if(b(i,j,m).gt.0.10) then
          if(a(i,j,m).gt.amx) then
            amx=a(i,j,m)
            mmx=m
          endif
        endif
      enddo
      mploc=mmx
      return
      end
      
C----------------------------------------------------------------
      function mmloc(a,i,j,im,jm,mm)
      dimension a(im,jm,mm)
      amx=0.
      mmx=1
      do m=2,mm
          if(a(i,j,m).gt.amx) then
            amx=a(i,j,m)
            mmx=m
          endif
c       if(i.eq.16.and.j.eq.31)then
c      write(6,'(''mmloc'',4i5,2e12.3)') i,j,m,mmx,amx,a(i,j,m)
c       endif
      enddo
      mmloc=mmx
      return
      end

C----------------------------------------------------------------
      function summ(a,dpi,i,j,im,jm,mm)
      dimension a(im,jm,mm)
      summ=0.
      do m=2,mm-1 
         summ=summ+a(i,j,m)*dpi  
c     if(i.eq.3.and.j.eq.23) write(6,'(i5,e12.3)')m,summ
      enddo
      return
      end
C----------------------------------------------------------------
      function summ1(a,dpi,mm)
      dimension a(mm)
      summ1=0.
      do m=2,mm-1
         summ1=summ1+a(m)*dpi
      enddo
      return
      end
C----------------------------------------------------------------
      function sumk(a,dpi,kb)
      dimension a(kb)
      sumk=0.
      do k=2,kb-1 
         sumk=sumk+a(k)*dpi  
      enddo
      return
      end

C----------------------------------------------------------------------
C                        END of WAVE CODE             !WAVE
C----------------------------------------------------------------------
