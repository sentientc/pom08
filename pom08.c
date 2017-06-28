C***********************************************************************
C
C     Common blocks for pom08.f
C
C***********************************************************************
C
C     Array sizes:
C
      integer
     $  i,j,k,
     $  im             ,imm1           ,imm2           ,jm             ,
     $  jmm1           ,jmm2           ,kb             ,kbm1           ,
     $  kbm2           ,mm             ,i1,i2,i3,j1,j2,j3,j4,k1,k2,k3   
C
C***********************************************************************
C
C     Set size of problem here:
C
C
      parameter
C -- seamount      (iproblem=1)
     $ (im=24        ,jm=11            ,kb=21)
c -- file2ic (iproblem=3)
C    $ (im=41          ,jm=61          ,kb=16)
C -- nybight  (iproblem=4)
C    $ (im=40          ,jm=40          ,kb=21)
c    $ (im=200         ,jm=200         ,kb=21)
c    $ (im=401        ,jm=251          ,kb=25) 

C    wave parameter mm = Number of angle segment + 2
      parameter (mm=24+2)
C***********************************************************************
C                                     
      parameter
     $ (imm1=im-1      ,imm2=im-2      ,jmm1=jm-1      ,jmm2=jm-2      ,
     $  kbm1=kb-1      ,kbm2=kb-2      )
C
C-----------------------------------------------------------------------
C
C     Scalars:
C
      real
     $  alpha          ,dte            ,dti            ,dti2           ,  
     $  grav           ,hmax           ,kappa          ,pi             ,
     $  ramp           ,rfe            ,rfn            ,rfs            ,
     $  rfw            ,rhoref         ,sbias          ,slmax          ,
     $  small          ,tbias          ,time           ,tprni          ,
     $  umol           ,vmax           

      integer
     $  iint           ,iprint         ,iskp           ,jskp           ,
     $  kl1            ,kl2            ,mode           ,modew          ,
     $  ntp            ,initwave,
     $  iperx          ,ipery          ,iproblem !lyo:_20080419:!lyo:_20080507:
C
      common/blkcon/ 
     $  alpha          ,dte            ,dti            ,dti2           ,
     $  grav           ,hmax           ,kappa          ,pi             ,
     $  ramp           ,rfe            ,rfn            ,rfs            ,
     $  rfw            ,rhoref         ,sbias          ,slmax          ,
     $  small          ,tbias          ,time           ,tprni          ,
     $  umol           ,vmax           ,
     $  iint           ,iprint         ,iskp           ,jskp           ,
     $  kl1            ,kl2            ,mode           ,modew          ,
     $  ntp            ,initwave,
     $  iperx          ,ipery          ,iproblem !lyo:_20080419:!lyo:_20080507:
C
C-----------------------------------------------------------------------
C
C     1-D arrays:
C
      real
     $  dz             ,dzz            ,z              ,zz
C
      common/blk1d/ 
     $  dz(kb)         ,dzz(kb)        ,z(kb)          ,zz(kb) 
C
C-----------------------------------------------------------------------
C
C     2-D arrays:
C
      real
     $  lat            ,long           , 
     $  aam2d          ,advua          ,advva          ,adx2d          ,
     $  ady2d          ,art            ,aru            ,arv            ,
     $  cbc            ,cor            ,d              ,drx2d          ,
     $  dry2d          ,dt             ,dum            ,dvm            ,
     $  dx             ,dy             ,east_c         ,east_e         ,
     $  east_u         ,east_v         ,e_atmos        ,egb            ,
     $  egf            ,el             ,elb            ,elf            ,
     $  et             ,etb            ,etf            ,fluxua         ,
     $  fluxva         ,fsm            ,h              ,north_c        ,
     $  north_e        ,north_u        ,north_v        ,psi            ,
     $  rot            ,ssurf          ,swrad          ,vfluxb         ,
     $  tps            ,tsurf          ,ua             ,vfluxf         ,
     $  uab            ,uaf            ,utb            ,utf            ,
     $  va             ,vab            ,vaf            ,
     $  vtb            ,vtf            ,wssurf         ,wtsurf         ,
     $  wubot          ,wusurf         ,wvbot          ,wvsurf         ,
     $  u10            ,u10x           ,u10y
C
      common/blk2d/
     $  lat(im,jm)     ,long(im,jm)    , 
     $  aam2d(im,jm)   ,advua(im,jm)   ,advva(im,jm)   ,adx2d(im,jm)   ,
     $  ady2d(im,jm)   ,art(im,jm)     ,aru(im,jm)     ,arv(im,jm)     ,
     $  cbc(im,jm)     ,cor(im,jm)     ,d(im,jm)       ,drx2d(im,jm)   ,
     $  dry2d(im,jm)   ,dt(im,jm)      ,dum(im,jm)     ,dvm(im,jm)     ,
     $  dx(im,jm)      ,dy(im,jm)      ,east_c(im,jm)  ,east_e(im,jm)  ,
     $  east_u(im,jm)  ,east_v(im,jm)  ,e_atmos(im,jm) ,egb(im,jm)     ,
     $  egf(im,jm)     ,el(im,jm)      ,elb(im,jm)     ,elf(im,jm)     ,
     $  et(im,jm)      ,etb(im,jm)     ,etf(im,jm)     ,fluxua(im,jm)  ,
     $  fluxva(im,jm)  ,fsm(im,jm)     ,h(im,jm)       ,north_c(im,jm) ,
     $  north_e(im,jm) ,north_u(im,jm) ,north_v(im,jm) ,psi(im,jm)     ,
     $  rot(im,jm)     ,ssurf(im,jm)   ,swrad(im,jm)   ,vfluxb(im,jm)  ,
     $  tps(im,jm)     ,tsurf(im,jm)   ,ua(im,jm)      ,vfluxf(im,jm)  ,
     $  uab(im,jm)     ,uaf(im,jm)     ,utb(im,jm)     ,utf(im,jm)     ,
     $  va(im,jm)      ,vab(im,jm)     ,vaf(im,jm)     ,
     $  vtb(im,jm)     ,vtf(im,jm)     ,wssurf(im,jm)  ,wtsurf(im,jm)  ,
     $  wubot(im,jm)   ,wusurf(im,jm)  ,wvbot(im,jm)   ,wvsurf(im,jm)  ,
     $  u10(im,jm)     ,u10x(im,jm)    ,u10y(im,jm)
C
C-----------------------------------------------------------------------
C
C     3-D arrays:
C
      real 
     $  aam            ,advx           ,advy           ,a              ,
     $  c              ,drhox          ,drhoy          ,dtef           ,
     $  ee             ,gg             ,kh             ,km             ,
     $  kq             ,l              ,q2b            ,q2             ,
     $  q2lb           ,q2l            ,rho            ,rmean          ,
     $  sb             ,sclim          ,s              ,tb             ,
     $  tclim          ,t              ,ub             ,uf             ,
     $  u              ,vb             ,vf             ,v              ,
     $  w              ,wr             ,zflux          ,rclim
!lyo:_20080602:add wr  !lyo:_20081226:add rclim
C
      common/blk3d/
     $  aam(im,jm,kb)  ,advx(im,jm,kb) ,advy(im,jm,kb) ,a(im,jm,kb)    ,
     $  c(im,jm,kb)    ,drhox(im,jm,kb),drhoy(im,jm,kb),dtef(im,jm,kb) ,
     $  ee(im,jm,kb)   ,gg(im,jm,kb)   ,kh(im,jm,kb)   ,km(im,jm,kb)   ,
     $  kq(im,jm,kb)   ,l(im,jm,kb)    ,q2b(im,jm,kb)  ,q2(im,jm,kb)   ,
     $  q2lb(im,jm,kb) ,q2l(im,jm,kb)  ,rho(im,jm,kb)  ,rmean(im,jm,kb),
     $  sb(im,jm,kb)   ,sclim(im,jm,kb),s(im,jm,kb)    ,tb(im,jm,kb)   ,
     $  tclim(im,jm,kb),t(im,jm,kb)    ,ub(im,jm,kb)   ,uf(im,jm,kb)   ,
     $  u(im,jm,kb)    ,vb(im,jm,kb)   ,vf(im,jm,kb)   ,v(im,jm,kb)    ,
     $  w(im,jm,kb)    ,wr(im,jm,kb)   ,zflux(im,jm,kb),rclim(im,jm,kb)
!lyo:_20080602:add wr  !lyo:_20081226:add rclim
C
C-----------------------------------------------------------------------
C
C     1 and 2-D boundary value arrays:
C
      real
     $  ele            ,eln            ,els            ,elw            ,
     $  sbe            ,sbn            ,sbs            ,sbw            ,
     $  tbe            ,tbn            ,tbs            ,tbw            ,
     $  uabe           ,uabw           ,ube            ,ubw            ,
     $  vabn           ,vabs           ,vbn            ,vbs       
C
      common/bdry/
     $  ele(jm)        ,eln(im)        ,els(im)        ,elw(jm)        ,
     $  sbe(jm,kb)     ,sbn(im,kb)     ,sbs(im,kb)     ,sbw(jm,kb)     ,
     $  tbe(jm,kb)     ,tbn(im,kb)     ,tbs(im,kb)     ,tbw(jm,kb)     ,
     $  uabe(jm)       ,uabw(jm)       ,ube(jm,kb)     ,ubw(jm,kb)     ,
     $  vabn(im)       ,vabs(im)       ,vbn(im,kb)     ,vbs(im,kb)
C
C-----------------------------------------------------------------------
C
C     Character variables:
C
      character*26
     $  time_start
C
      character*40
     $  source,title
C
      common/blkchar/
     $  time_start     ,source         ,title
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C   Common variables associated with the wave subroutines
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      integer
     &    m,nwave,isplitw
      real 
     &    dth,dthinv,dtw,dtw2,dtww,ent,enwnd,thtav,thtpws,Hs,
     &    aren,areninv,x,y,
     &    cs,sn,theta,thetam,csm,snm,
     &    enw,ens,ene,enn,
     &    sigthw,sigths,sigthe,sigthn,
     &    sigthf,sigth,sigthb,sigthav,
     &    cth,cg,kth,kthD,kthav,kthavD,
     &    kpD,cgx,cgy,
     &    age,ageinv,tht_wnd,
     &    sigp,sigpws,ustw,cd,
     &    cp,kp, cthth,
     &    enb,en,enf,
     &    qdis,Hm,bdis1,bdis2,
     &    sdis,swin,fspr,
     &    swintot,sdistot,bdistot,ustp,
     &    tps2d,ud,vd,
C variables used by circulation model
     &    tpx0,tpy0,ttx0,tty0,
     &    tpx,tpy,tpzdist,ttx,tty,
     &    ust,vst,
     &    ucb,ucf,vcb,vcf,
     &    Sxx,Sxy,Syy,
     &    F1,F2,F3,F4
C 
      common/blkwave/ 
     &   dthinv,
     &    ent(im,jm),enwnd(im,jm),thtav(im,jm),thtpws(im,jm),Hs(im,jm),
     &    aren(im,jm),areninv(im,jm),x(im),y(jm),
     &    cs(mm),sn(mm),theta(mm),thetam(mm),csm(mm),snm(mm),
     &    enw(jm,mm),ens(im,mm),ene(jm,mm),enn(im,mm),
     &    sigthw(jm,mm),sigths(im,mm),sigthe(jm,mm),sigthn(im,mm),
     &    sigthf(im,jm,mm),sigth(im,jm,mm),sigthb(im,jm,mm),
     &    sigthav(im,jm),
     &    cth(im,jm,mm),cg(im,jm,mm),kth(im,jm,mm),kthD(im,jm,mm),
     &    kpD(im,jm),kthav(im,jm),kthavD(im,jm),
     &    cgx(im,jm,mm),cgy(im,jm,mm),
     &    age(im,jm),ageinv(im,jm),tht_wnd(im,jm),
     &    sigp(im,jm),sigpws(im,jm),ustw(im,jm),cd(im,jm),
     &    cp(im,jm),kp(im,jm), cthth(im,jm,mm),
     &    enb(im,jm,mm),en(im,jm,mm),enf(im,jm,mm),
     &    qdis(im,jm),Hm(im,jm),bdis1(im,jm,mm),bdis2(im,jm,mm),
     &    sdis(im,jm,mm),swin(im,jm,mm),fspr(im,jm,mm),
     &    swintot(im,jm),sdistot(im,jm),bdistot(im,jm),
     &    ustp(im,jm),
     &    tps2d(im,jm),ud(im,jm),vd(im,jm),
C variables used by circulation model
     &    tpx0(im,jm),tpy0(im,jm),ttx0(im,jm),tty0(im,jm),
     &    tpx(im,jm,kb),tpy(im,jm,kb),tpzdist(im,jm,kb),
     &    ttx(im,jm,kb),tty(im,jm,kb),
     &    ust(im,jm,kb),vst(im,jm,kb),
     &    ucb(im,jm,kb),ucf(im,jm,kb),vcb(im,jm,kb),vcf(im,jm,kb),
     &    Sxx(im,jm,kb),Sxy(im,jm,kb),Syy(im,jm,kb),
     &    F1(kb),F2(kb),F3(kb),F4(kb)


!                                                                      !
!----------------------------------------------------------------------!
!lyo:!wad:begins:                                                      !
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!   Common variables associated with WAD                               !
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
      integer nwad
!
      real 
     $    alph0,hhi0,tidamp,hco,
     $    zsh,
     $    hmin,hc,hhi,
     $    dtrat,cwetrlx1,cwetrlx2 !Relax params. for when river water 
                                  !is discharged onto dry cells
!
      common/blkconwad/ 
     $  alph0          ,hhi0           ,tidamp         ,hco            ,
     $  zsh            ,hmin           ,hc             ,hhi            ,
     $  dtrat          ,cwetrlx1       ,cwetrlx2       ,nwad
!
      real
     $  wriv           ,dwet           ,wetmask        ,wmarsh         ,
     $  ddf            ,dd0            ,ddb            ,ssurfwet       ,
     $  hkeep
!
      common/blk2dwad/    !ALL but wriv are used in main prog. only:
     $    wriv(im,jm),    !used in main & subr.advt!not implemented yet
     $    dwet(im,jm),wetmask(im,jm),wmarsh(im,jm),
     $    ddf(IM,JM),dd0(IM,JM),ddb(IM,JM),
     $    hkeep(im,jm),ssurfwet(im,jm) !used for river specifications
!
      real
     $  pdens
!
      common/blk3dwad/
     $    pdens(im,jm,kb)              !used in dens & profq
!                                                                      !
!lyo:!wad:ends:                                                        !
!----------------------------------------------------------------------!
!                                                                      !
C-----------------------------------------------------------------------
C
C     End of common blocks
C
C-----------------------------------------------------------------------
C
