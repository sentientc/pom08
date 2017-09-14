      subroutine get_wind
      
       use module_time
c      use interp
       implicit none
C
       include 'pom08.c'
  
       type(date)  :: mtime,mtimes(2),mtimei
       character*19 :: tnm
       character*62 :: wfile 
       character*19 :: tdate
       integer :: ti,tis(2)
c       real :: wu(im,jm,2),wv(im,jm,2)
c       real :: wu1(im,jm),wv1(im,jm)
       real :: td(2),tt,r1,r2
       mtime=str2date(time_start)
       mtime=mtime+int(time*86400)
       if (time.le.dti/86400) then
        tis(:)=-1
       endif 
c       tnm=date2str(mtime)
c       mtime=str2date(tnm)
c       hr=mtime%hour+(mtime%min+mtime%sec/60)/60
       do i=1,2
        mtimei=mtime
c!sc:none merged wind(each nc file is 6 hourly wind data a day)
c        ti=int(hr/6)+(i-1)
c        if (ti.eq.4) then
c         mtimei=mtime+86400
c         ti=0
c        endif
c        tdate=date2str(mtimei)
c        td(i)=abs(mtime-(str2date(tdate(1:10)//'_00_00_00')+ti*6*3600))
c        write( wfile, '( a7,"_",i4.4,2i2.2,".nc" )' )
c     $  'mergew1', mtimei%year, mtimei%month, mtimei%day     
c!sc:merged wind(each nc file is 6 hourly wind data for a year)
        tdate=date2str(mtimei)
        ti=int((mtimei-str2date(tdate(1:4)//'_01_01_00_00_00'))/3600/6)
        ti=ti+(i-1)
        write(wfile,'(i4,"-01-01/out/SRF.globalexp001."
     &  ,i4,"-01-01_00_",i4,"-12-31_18.nc")' )
     &   mtimei%year,mtimei%year,mtimei%year
        if (tis(i).ne.ti) then
         tis(i)=ti
         call read_wind(wfile,ti+1,wu(:,:,i),wv(:,:,i),tt)
         mtimes(i)=str2date('1950-01-01-00-00-00')
         mtimes(i)=mtimes(i)+int(tt*86400)
        else
        endif
        td(i)=abs(mtime-mtimes(i))
       enddo
       r1=td(2)/(td(1)+td(2))
       r2=td(1)/(td(1)+td(2))
c       do j=1,jm
c         do i=1,im
c         u10x(i,j)=wu(i,j,1)*r1+wu(i,j,2)*r2
c         u10y(i,j)=wv(i,j,1)*r1+wv(i,j,2)*r2
c        enddo
c       enddo
      u10x(:,:)=r1*wu(:,:,1)+r2*wu(:,:,2)
      u10y(:,:)=r1*wv(:,:,1)+r2*wv(:,:,2)
      end
      subroutine get_time(dtnow,treft,rt)

      use module_time
      include 'pom08.c'
      type(date)  :: rtime,reft
      real ::rt
      real ::dtnow
      character ::treft*19
      rtime=str2date(time_start)
      rtime=rtime+int(dtnow*86400)
c      reft=str2date('1950-01-01-00-00-00')
      reft=str2date(treft)
      rt=(rtime-reft)/86400
      return
      end
