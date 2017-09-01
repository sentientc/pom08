      module  module_time
      
      implicit none
      
!      include 'pom.h'
      
      type date
      integer :: year, month, day, hour, min, sec
      end type date

      interface operator(+)
      module procedure add_date
      end interface

      interface operator(-)
      module procedure dif_date
      end interface
  
      interface operator(>=)
      module procedure ge_date
      end interface

      interface operator(>)
      module procedure gt_date
      end interface

      contains

!=====================================================================
! whether d1 is equal to or earlier than d2 
!---------------------------------------------------------------------
      logical function ge_date( d1, d2 )


      type(date), intent(in) :: d1, d2
      

      ge_date = lge( date2str( d1 ), date2str( d2 ) )
      

      end function ge_date
!=====================================================================


!=====================================================================
! whether d1 is earlier than d2 
!---------------------------------------------------------------------
      logical function gt_date( d1, d2 )


      type(date), intent(in) :: d1, d2
      

      gt_date = lgt( date2str( d1 ), date2str( d2 ) )
      

      end function gt_date
!=====================================================================


!=====================================================================
! char --> type
!--------------------------------------------------------
      type(date) function str2date( str ) result(d)
 
      character*19 :: str
           
      read(str, '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                    
     $     d%year, d%month, d%day, d%hour, d%min, d%sec
      
      
      end function str2date
!=====================================================================


!=====================================================================
! type --> char
!--------------------------------------------------------      
      character(19) function date2str( d ) result(str)
 
      type(date) :: d
           
      write(str, 
     $     '(i4.4,"-",i2.2,"-",i2.2,"_",i2.2,"_",i2.2,"_",i2.2)' )
     $     d%year, d%month, d%day, d%hour, d%min, d%sec
      
      
      end function date2str
!=====================================================================     


!=====================================================================     
! add dt_in [sec] to d
!---------------------------------------------------------------------
      type(date) function add_date( d, dt_in ) result(d_result)


      type(date), intent(in) :: d
      integer, intent(in) :: dt_in

      integer :: mday(0:12) = (/31, 31, 28, 31, 30, 31, 30,
     $                         31, 31, 30, 31, 30, 31/)

     
!     dt_in must be positive.
      if  ( dt_in .lt. 0 ) then
         write(*,'(a,i)') 
     $  " Eror. module_time : dt_in must be positive. dt_in = ", dt_in
         stop
      endif

!     check leap year. 
                             
      if( ( mod( d%year, 4 ) == 0 .and. mod( d%year, 100 ) /= 0 ) .or.   
     $                             mod( d%year, 400 ) == 0 ) then
         mday(2) = 29
      else
         mday(2) = 28
      endif

       
!     seconds
         d_result%sec = d%sec + dt_in
         d_result%min = d_result%sec / 60
         d_result%sec = mod( d_result%sec, 60 )

!     minutes
         d_result%min  = d_result%min + d%min
         d_result%hour = d_result%min / 60
         d_result%min  = mod( d_result%min, 60 )
       
!     hour
         d_result%hour = d_result%hour + d%hour
         d_result%day  = d_result%hour / 24
         d_result%hour = mod( d_result%hour, 24 )
       
!!     day
!         d_result%day = d_result%day + d%day
!         d_result%month = ( d_result%day - 1 ) / mday( d%month )
!         if ( d_result%month .gt. 1 ) then 
!            write( *, * ) 
!     $     'Error: module_time,  d_result%month > 1.',
!     $           d_result%month
!            stop
!         end if
!         d_result%day = mod( d_result%day, mday( d%month ) )
!         if( d_result%day == 0 ) then
!            d_result%day = mday( d%month )
!         endif
!       
!!     month
!         d_result%month = d_result%month + d%month
!         d_result%year = d%year + ( d_result%month - 1 ) / 12
!         d_result%month = mod( d_result%month, 12 )
!         if( d_result%month == 0 ) then
!            d_result%month = 12
!         endif
 

         d_result%day = d_result%day + d%day
         d_result%month = d%month
         d_result%year = d%year
         do while ( d_result%day > mday( d_result%month ) ) 
            d_result%day = d_result%day - mday( d_result%month )
            d_result%month = d_result%month + 1
            if ( d_result%month == 13 ) then
               d_result%year = d_result%year + 1 
               d_result%month = 1 
!     check leap year. 
               if( ( mod( d_result%year, 4 ) == 0 
     $              .and. mod( d_result%year, 100 ) /= 0 ) .or.   
     $              mod( d_result%year, 400 ) == 0 ) then
                  mday(2) = 29
               else
                  mday(2) = 28
               endif
            endif
         enddo
            
                  
         end function add_date
!=====================================================================     


!=====================================================================     
! calculate interval [sec] from d1 to d2 ( d1 > d2 )
!---------------------------------------------------------------------
      integer function dif_date( d1, d2 )
      
      
      type(date), intent(in) :: d1, d2
      
      integer :: mday(12) = (/ 31, 28, 31, 30, 31, 30,        
     $                         31, 31, 30, 31, 30, 31 /)
      
      integer :: day1, day2, sec1, sec2, iy

!     Number of days since 1 January d2%year.

      mday(2) = 28 + inc_leap(d1%year)
      day1 = sum( mday(1:d1%month) ) - mday(d1%month) + d1%day
      mday(2) = 28 + inc_leap(d2%year)
      day2 = sum( mday(1:d2%month) ) - mday(d2%month) + d2%day

      do iy = d2%year, d1%year - 1
         day1 = day1 + 365 + inc_leap(iy)
      end do

!     Number of seconds since 0h in each day.

      sec1 = 3600 * d1%hour + 60 * d1%min + d1%sec
      sec2 = 3600 * d2%hour + 60 * d2%min + d2%sec

!     sum      

      dif_date = 86400 * (day1 - day2) + (sec1 - sec2) 


      end function dif_date
!=====================================================================     


!=====================================================================     
!                          
!---------------------------------------------------------------------            
      integer function inc_leap( year )

      integer :: year
      
      if ( ( mod(year, 4) == 0 .and. mod(year, 100) /= 0 ) .or.  
     $     mod(year, 400) == 0 ) then
         inc_leap = 1
      else
         inc_leap = 0
      end if
      
      
      end function inc_leap
!=====================================================================     

      end module
