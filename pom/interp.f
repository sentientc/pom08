! ----------------------------------------------------------------
      subroutine interp_wind(uwnd_out, vwnd_out, 
     $     lon_wind,lat_wind,uwnd_wind_in,vwnd_wind_in)

      implicit none
      
      include 'pom08.c'
      
      real, dimension(im,jm) :: uwnd_out, vwnd_out
      real, dimension(1439,602) :: 
     $     uwnd_wind_in, vwnd_wind_in
      real, dimension(1439) :: lon_wind
      real, dimension(602) :: lat_wind      

      integer :: i_tmp, j_tmp 
      real, parameter :: dx1 = 0.25
      real :: lon_tmp, lat_tmp, uwnd_tmp, vwnd_tmp
      real, dimension(4) :: x1, y1, ff
      real :: bilinear_interp, bilinear_interp0


      do j=1,jm
         do i=1,im

!            if ( lon ( i, j ) .lt. 0. ) then
!               lon_tmp = lon( i, j ) + 360. 
!               i_tmp = int( ( lon_tmp ) / dx ) + 1
!            else 
!               lon_tmp = lon( i, j )
!               i_tmp = int( lon( i, j ) / dx ) + 1
!            endif         

!            lat_tmp = lat( i, j )            
!            j_tmp = int( ( lat( i, j ) + 78.3750 )/ dx ) + 1

             lon_tmp = east_e( i, j ) 
!            i_tmp = int( (lon_tmp - lonl) / dx ) + 1
            i_tmp = int( (lon_tmp -  lon_wind(1)) / dx1 ) + 1
            lat_tmp = north_e( i, j ) 
!            j_tmp = int( (lat_tmp - latb) / dx ) + 1           
            j_tmp = int( (lat_tmp - lat_wind(1)) / dx1 ) + 1           

!     surrounding 4 points

            x1(1) = lon_wind( i_tmp )
            x1(2) = lon_wind( i_tmp )
            x1(3) = lon_wind( i_tmp + 1 ) 
            x1(4) = lon_wind( i_tmp + 1 )
            
            y1(1) = lat_wind( j_tmp )
            y1(2) = lat_wind( j_tmp + 1 )
            y1(3) = lat_wind( j_tmp + 1 )
            y1(4) = lat_wind( j_tmp )

!     interpolation for u wind velocity

            ff(1) = uwnd_wind_in( i_tmp    , j_tmp     )
            ff(2) = uwnd_wind_in( i_tmp    , j_tmp + 1 )
            ff(3) = uwnd_wind_in( i_tmp + 1, j_tmp + 1 )
            ff(4) = uwnd_wind_in( i_tmp + 1, j_tmp     )
  
            uwnd_tmp 
     $           = bilinear_interp0( x1, y1, ff, lon_tmp, lat_tmp )


!     interpolation for v wind velocity

            ff(1) = vwnd_wind_in( i_tmp    , j_tmp     )
            ff(2) = vwnd_wind_in( i_tmp    , j_tmp + 1 )
            ff(3) = vwnd_wind_in( i_tmp + 1, j_tmp + 1 )
            ff(4) = vwnd_wind_in( i_tmp + 1, j_tmp     )
  
            vwnd_tmp 
     $           = bilinear_interp0( x1, y1, ff, lon_tmp, lat_tmp )
            

!     rotation 
            uwnd_out( i, j ) 
     $           =   uwnd_tmp * cos( rot( i, j ) * pi / 180.0 )
     $             + vwnd_tmp * sin( rot( i, j ) * pi / 180.0 ) 

            vwnd_out( i, j ) 
     $           = - uwnd_tmp * sin( rot( i, j ) * pi / 180.0 )
     $             + vwnd_tmp * cos( rot( i, j ) * pi / 180.0 ) 


      goto 100
      write(*,*)
      write(*,*)'---------from interp_wind subr.------------------'
      write(*,*) 'i,j=',i,j
      write(*,*) 'lon_tmp,lat_tmp=',lon_tmp,lat_tmp
      write(*,*) 'i_tmp,j_tmp=',i_tmp,j_tmp
      write(*,*) 'x(1:4)', x1(1:4) 
      write(*,*) 'y(1:4)', y1(1:4) 
      write(*,*) 'rot(i,j)',rot( i, j )
      write(*,*) 'cos,sin='
     $,cos( rot( i, j ) * pi / 180.0 ),sin( rot( i, j ) * pi / 180.0 )
      write(*,*) 'uf(1:4)' 
     $  ,uwnd_wind_in( i_tmp    , j_tmp     )
     $  ,uwnd_wind_in( i_tmp    , j_tmp + 1 )
     $  ,uwnd_wind_in( i_tmp + 1, j_tmp + 1 )
     $  ,uwnd_wind_in( i_tmp + 1, j_tmp     )
      write(*,*) 'vf(1:4)' 
     $  ,vwnd_wind_in( i_tmp    , j_tmp     )
     $  ,vwnd_wind_in( i_tmp    , j_tmp + 1 )
     $  ,vwnd_wind_in( i_tmp + 1, j_tmp + 1 )
     $  ,vwnd_wind_in( i_tmp + 1, j_tmp     )
      write(*,*) 'uwnd_tmp,uwnd_out(i,j)=', uwnd_tmp,uwnd_out(i,j)
      write(*,*) 'vwnd_tmp,vwnd_out(i,j)=', vwnd_tmp,vwnd_out(i,j)
      write(*,*)
      write(*,*)'-------------------------------------------------'
100   continue
         enddo
      enddo

      return

      end

! --------------------------------------------------
      real function bilinear_interp( x_ref, y_ref, f_ref
     $     , x_tar, y_tar ) result( var_out )


      implicit none

      real, dimension(4), intent(in) :: x_ref, y_ref, f_ref
      real, intent(in) :: x_tar, y_tar
      
      real :: a1 ,a2, a3, a4, b1, b2, b3, b4
      real :: aa, bb, cc
      real :: ss, tt

      ss  = ( x_tar - x_ref(1) ) / ( x_ref(4) - x_ref(1) ) 
      tt  = ( y_tar - y_ref(1) ) / ( y_ref(2) - y_ref(1) ) 


      var_out = f_ref(1) * ( 1.0 - tt ) * ( 1.0 - ss ) +
     $          f_ref(2) *         tt   * ( 1.0 - ss ) +
     $          f_ref(3) *         tt   *         ss   +
     $          f_ref(4) * ( 1.0 - tt ) *         ss

      return 
      
      end

!     --------------------------------------------------
      real function bilinear_interp0( x_ref, y_ref, f_ref
     $     , x_tar, y_tar ) result( var_out )

      implicit none

      real, dimension(4), intent(in) :: x_ref, y_ref, f_ref
      real, intent(in) :: x_tar, y_tar
      
      real :: a1 ,a2, a3, a4, b1, b2, b3, b4
      real :: aa, bb, cc
      real :: ss, tt

      a1 =   x_ref(1) - x_ref(2) + x_ref(3) - x_ref(4)
      a2 = - x_ref(1) + x_ref(4)
      a3 = - x_ref(1) + x_ref(2)
      a4 =   x_ref(1) - x_tar
      b1 =   y_ref(1) - y_ref(2) + y_ref(3) - y_ref(4)
      b2 = - y_ref(1) + y_ref(4)
      b3 = - y_ref(1) + y_ref(2)
      b4 =   y_ref(1) - y_tar

      aa =   a3 * b1 - a1 * b3
      bb =   b2 * a3 + b1 * a4 - a1 * b4 - a2 * b3
      cc = - a2 * b4 + a4 * b2

      if ( abs( aa*cc ) .gt. 0.002*bb**2 ) then 
         tt = ( - bb - sqrt( bb**2 - 4.0 * aa * cc ) ) 
     $        / ( 2.0 * aa )
      else
         tt = cc / abs( bb )
      endif

      aa =   a2 * b1 - a1 * b2
      bb =   b3 * a2 + b1 * a4 - a1 * b4 - a3 * b2
      cc = - a3 * b4 + a4 * b3

      if ( abs( aa * cc ) .gt. 0.002*bb**2 ) then
         ss = ( - bb - sqrt( bb**2 - 4.0 * aa * cc ) )
     $        / ( 2.0 * aa )
      else
         ss = - cc / abs( bb )
      endif

      var_out = f_ref(1) * ( 1.0 - tt ) * ( 1.0 - ss ) +
     $          f_ref(2) *         tt   * ( 1.0 - ss ) +
     $          f_ref(3) *         tt   *         ss   +
     $          f_ref(4) * ( 1.0 - tt ) *         ss

      return 
      
      end

!-------------------------------------------------------
