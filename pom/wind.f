      subroutine get_wind(mtime,wu1,wv1)
      
       use module_time
c      use interp
       implicit none
C
       include 'pom08.c'
  
       type(date)  :: mtime,mtime0,mtimei
       character*8 :: tnm
       character*16 :: wfile 
       character*19 :: tdate
       integer :: ti
       real :: wu(im,jm,2),wv(im,jm,2)
       real :: wu1(im,jm),wv1(im,jm)
       real :: hr,td(2)
       tnm=date2str(mtime)
       hr=mtime%hour+(mtime%min+mtime%sec/60)/60
       do i=1,2
        mtimei=mtime
        ti=int(hr/6)+(i-1)
        if (ti.eq.4) then
         mtimei=mtime+86400
         ti=0
        endif
        tdate=date2str(mtimei)
        td(i)=abs(mtime-(str2date(tdate(1:10)//'_00_00_00')+ti*6*3600))
        write( wfile, '( a7,"_",i4.4,2i2.2,".nc" )' )
     $  'mergew1', mtimei%year, mtimei%month, mtimei%day     
        call read_wind(wfile,ti+1,wu(:,:,i),wv(:,:,i))
       enddo
       wu1=(1-td(1)/(td(1)+td(2)))*wu(:,:,1)
     & +(1-td(2)/(td(1)+td(2)))*wu(:,:,2)
       wv1=(1-td(1)/(td(1)+td(2)))*wv(:,:,1)
     & +(1-td(2)/(td(1)+td(2)))*wv(:,:,2)
      end
