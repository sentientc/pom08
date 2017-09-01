      module interp

      implicit none
      
      private

      public :: interp_init,interp_mask_2d,bilinear_coef_mask
     $          ,interp_mask_3d

      include 'pom.h'

      character*16 :: infile, gridfile

!      integer, parameter :: x_division = 2, y_division = 2 
      integer, parameter :: is = 1,  js = 1
      real,    parameter :: bound = 0.0
      contains
!################################################################################
!-------------------------------------------------------------------------------- 
      subroutine interp_init  
! read in original grids
!--------------------------------------------------------------------------------

       implicit none
       integer i,j,k


! read grid
        call read_grid_pnetcdf0
       return
       end subroutine interp_init  
!################################################################################
!--------------------------------------------------------------------------------

      subroutine interp_mask_2d(
     $    var_coarse, flag_uv, alon_fine, alat_fine,
     $    mask_fine, var_fine)

!--------------------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------------------


      real, dimension(im_coarse, jm_coarse) ::  var_coarse

      real, dimension(im,jm) :: alon_fine, alat_fine, mask_fine

      integer, intent(in) :: flag_uv
      real, intent(out)   :: var_fine(im, jm)

!--------------------------------------------------------------------------------

      integer :: i, j, i1, j1, flag_interp, ii
      real :: x, y, x1, x2, x3, x4, y1, y2, y3, y4, 
     $        m1, m2, m3, m4, a1, a2, a3, a4

!      integer :: ie1, je1, im_fine1, jm_fine1

      real, dimension(im_coarse, jm_coarse) :: 
     $          alon_coarse1, alat_coarse1
      real, dimension(im,jm) :: alon_fine1, alat_fine1

!--------------------------------------------------------------------------------

! --- For e-points.

       if (flag_uv == 0) then
         alon_coarse1 = alon_coarse;
         alat_coarse1 = alat_coarse
         alon_fine1 = alon_fine;
         alat_fine1 = alat_fine; 
       endif

! --- For u-points.

       if (flag_uv == 1) then

!------------------------

       alon_coarse1(2:im_coarse, 1:jm_coarse) =            
     $    (alon_coarse(2:im_coarse  , 1:jm_coarse) 
     $    + alon_coarse(1:im_coarse-1, 1:jm_coarse) ) * 0.5
       alon_coarse1(1, :) = 2.*alon_coarse1(2, :) - alon_coarse1(3, :)

!------------------------

       alat_coarse1(2:im_coarse, 1:jm_coarse) =             
     $    (alat_coarse(2:im_coarse  , 1:jm_coarse) 
     $    + alat_coarse(1:im_coarse-1, 1:jm_coarse) ) * 0.5
      alat_coarse1(1, :) = 2.*alat_coarse1(2, :) - alat_coarse1(3, :)

!------------------------

      alon_fine1(2:im, 1:jm) = ( alon_fine(2:im  , 1:jm)       
     $                         + alon_fine(1:im-1, 1:jm) ) * 0.5
      alon_fine1(1, :) = 2.*alon_fine1(2, :) - alon_fine1(3, :)

!------------------------

      alat_fine1(2:im, 1:jm) = ( alat_fine(2:im  , 1:jm)       
     $                         + alat_fine(1:im-1, 1:jm) ) * 0.5
      alat_fine1(1, :) = 2.*alat_fine1(2, :) - alat_fine1(3, :)

      endif

! --- For v-points.

      if (flag_uv == 2) then

!------------------------

      alon_coarse1(1:im_coarse, 2:jm_coarse) =           
     $  ( alon_coarse(1:im_coarse, 2:jm_coarse) 
     $  + alon_coarse(1:im_coarse, 1:jm_coarse-1) ) * 0.5
      alon_coarse1(:, 1) = 2.*alon_coarse1(:, 2) - alon_coarse1(:, 3)

!------------------------

       alat_coarse1(1:im_coarse, 2:jm_coarse) =           
     $  ( alat_coarse(1:im_coarse, 2:jm_coarse) 
     $  + alat_coarse(1:im_coarse, 1:jm_coarse-1) ) * 0.5
       alat_coarse1(:, 1) = 2.*alat_coarse1(:, 2) - alat_coarse1(:, 3)

!------------------------

       alon_fine1(1:im, 2:jm) = ( alon_fine(1:im, 2:jm) 
     $                       + alon_fine(1:im, 1:jm-1) ) * 0.5
       alon_fine1(:, 1) = 2.*alon_fine1(:, 2) - alon_fine1(:, 3)

!------------------------

       alat_fine1(1:im, 2:jm) = ( alat_fine(1:im, 2:jm) 
     $                       + alat_fine(1:im, 1:jm-1) ) * 0.5
       alat_fine1(:, 1) = 2.*alat_fine1(:, 2) - alat_fine1(:, 3)

       endif

!--------------------------------------------------------------------------------

!  if (flag_uv /= 0) then

!    write(32)flag_uv, im, jm, im_fine, jm_fine,                &
!             alon_coarse, alat_coarse, alon_coarse1, alat_coarse1, &
!             alon_fine, alat_fine, alon_fine1, alat_fine1

!  endif

!--------------------------------------------------------------------------------

       call exchange2d_mpi(alon_coarse1,im_coarse,jm_coarse)
       call exchange2d_mpi(alat_coarse1,im_coarse,jm_coarse)
       call exchange2d_mpi(alon_fine1,im,jm)
       call exchange2d_mpi(alat_fine1,im,jm)

       do j = 1, jm
          do i = 1, im

             flag_interp = -1; ii = 1

             i1 = float(i-1)/x_division + is - 0.5
             j1 = float(j-1)/y_division + js - 0.5
             if(flag_uv == 2)i1 = (i-1)/x_division + is
             if(flag_uv == 1)j1 = (j-1)/y_division + js

             if (ii == 2) then

                if (flag_interp == 0) cycle

                select case(flag_interp)
                case (1, 2, 3)
                   i1 = i1 - 1
                case (5, 6, 7)
                   i1 = i1 + 1
                case default
                   i1 = i1
                end select

                select case(flag_interp)
                case (3, 4, 5)
                   j1 = j1 - 1
                case (1, 8, 7)
                   j1 = j1 + 1
                case default
                   i1 = j1
                end select

             endif

!         if(i1 <= 0 .or. j1 <= 0 .or. i1 >= im_coarse 
!     $         .or. j1 >= jm_coarse ) then
!             var_fine(i, j) = bound
!             cycle
!         endif

         if (i1 <= 0) i1 = 1
         if (j1 <= 0) j1 = 1
         if (i1 >= im_coarse) i1 = im_coarse - 1
         if (j1 >= jm_coarse) j1 = jm_coarse - 1

         x1 = alon_coarse1(i1  , j1  )
         x2 = alon_coarse1(i1  , j1+1)
         x3 = alon_coarse1(i1+1, j1+1)
         x4 = alon_coarse1(i1+1, j1  )
         
         y1 = alat_coarse1(i1  , j1  )
         y2 = alat_coarse1(i1  , j1+1)
         y3 = alat_coarse1(i1+1, j1+1)
         y4 = alat_coarse1(i1+1, j1  )

         m1 = mask_coarse(i1  , j1  )
         m2 = mask_coarse(i1  , j1+1)
         m3 = mask_coarse(i1+1, j1+1)
         m4 = mask_coarse(i1+1, j1  )
         
         if(flag_uv.ne.0) then
            m1=1;m2=1;m3=1;m4=1
         endif

         x = alon_fine1(i, j)
         y = alat_fine1(i, j)

         call bilinear_coef_mask(x, y, x1, x2, x3, x4, y1, y2, y3, y4,  
     $          m1, m2, m3, m4, a1, a2, a3, a4, flag_interp)

!    if(a1*a2*a3*a4 /= 0.0) then
!    if (flag_interp /= 0 .and. a1*a2*a3*a4 /= 0.0) then
!         if (flag_interp /= 0 ) then

!            write(10, *) 'iji1j1 = ', i, j, i1, j1, flag_uv
!            write(10, *) 'a = ', a1, a2, a3, a4, flag_interp

!         endif
         

         var_fine(i, j) = a1 * var_coarse(i1  , j1  )  
     $              + a2 * var_coarse(i1  , j1+1)  
     $              + a3 * var_coarse(i1+1, j1+1)  
     $              + a4 * var_coarse(i1+1, j1  )


!         var_fine(i, j) = var_fine(i, j) * mask_fine(i, j)

!   enddo
         enddo
       enddo

!      if(my_task == master_task) then

!       write(10, *) 'interp_wind:' 
!       write(10, *) 'wind_coarse(1,1)', var_coarse(1,1)
!       write(10, *) 'wind_coarse(im,jm)', var_coarse(im_coarse,jm_coarse)
!       write(10, *) ' '
!       write(10, *) 'wind_fine(1,1)', var_fine(1,1)
!       write(10, *) '(100,100)',var_fine(100,100)
!       write(10,*) 'alon,alat',alon_coarse(100,100),alat_coarse(100,100)
!       write(10,*) 'mask',mask_coarse(100,100) 
!       write(10, *) '(im,jm)',var_fine(im,jm)
!       write(10, *) '(2,1),(1,2)', var_fine(2,1),var_fine(1,2) 
!       write(10, *) ' '
!      endif  

!--------------------------------------------------------------------------------

       return

!--------------------------------------------------------------------------------

       end subroutine interp_mask_2d

!--------------------------------------------------------------------------------
      subroutine interp_mask_3d(
     $    var_coarse, flag_uv, alon_fine, alat_fine,
     $    mask_fine, var_fine)

!--------------------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------------------


      real, dimension(im_coarse, jm_coarse, kb) ::  var_coarse

      real, dimension(im,jm) :: alon_fine, alat_fine, mask_fine

      integer, intent(in) :: flag_uv
      real, intent(out)   :: var_fine(im, jm, kb)

!--------------------------------------------------------------------------------
      integer :: k
!--------------------------------------------------------------------------------

       do k = 1, kb
        call interp_mask_2d(var_coarse(:, :, k), flag_uv, 
     $   alon_fine,  alat_fine,  mask_fine,  var_fine(:, :, k))
       enddo

!--------------------------------------------------------------------------------

       return

!--------------------------------------------------------------------------------

       end subroutine interp_mask_3d


!--------------------------------------------------------------------------------
!################################################################################
!--------------------------------------------------------------------------------

       subroutine bilinear_coef_mask(x, y, x1, x2, x3, x4, y1, y2, y3,  
     $                     y4, m1, m2, m3, m4, w1, w2, w3, w4, flag   )

!--------------------------------------------------------------------------------

         implicit none

!--------------------------------------------------------------------------------

         real, intent(in) :: x, y, x1, x2, x3, x4, y1, y2, y3, y4, 
     $                       m1, m2, m3, m4
         real, intent(out) :: w1, w2, w3, w4
         integer, intent(out) :: flag

!--------------------------------------------------------------------------------

         real :: t, s, a1, a2, a3, a4, b1, b2, b3, b4, A, B, C, 
     $           tot_w, tot_m

!--------------------------------------------------------------------------------

! --- Pre-process for input.

         a1 = x1 - x2 + x3 - x4
         a2 = x4 - x1
         a3 = x2 - x1
         a4 = x1 - x

         b1 = y1 - y2 + y3 - y4
         b2 = y4 - y1
         b3 = y2 - y1
         b4 = y1 - y

!--------------------------------------------------------------------------------

! --- Solve for t;

         A = a3 * b1 - a1 * b3
         B = a4 * b1 + a3 * b2 - a1 * b4 - a2 * b3
         C = a4 * b2 - a2 * b4

         if(A /= 0) then
!  if(abs(A*C > 0.002*B**2))then
            t = ( - B - sqrt(B**2 - 4*A*C) ) / (2*A)
         else
            t = - C / B
         endif

!--------------------------------------------------------------------------------

! --- Solve for s

         A = a2 * b1 - a1 * b2
         B = a4 * b1 + a2 * b3 - a1 * b4 - a3 * b2
         C = a4 * b3 - a3 * b4

!  if(A == 0) then
!    s = - C / B
!  elseif(abs(A*C > 0.002*B**2))then
!    s = ( - B + sqrt(B**2 - 4*A*C) ) / (2*A)
!  else
!    s = 0.0
!  endif
         if(A /= 0) then
!  if(abs(A*C > 0.002*B**2))then
            s = ( - B + sqrt(B**2 - 4*A*C) ) / (2*A)
         else
            s = - C / B
         endif

!--------------------------------------------------------------------------------

! --- For w1, w2, w3, w4

         w1 = (1.-t) * (1.-s);
         w2 = t * (1.-s);
         w3 = s * t;
         w4 = (1.-t) * s;

!--------------------------------------------------------------------------------

! --- Adjust for normalization.
!
! *** Location od Flags
!  -----------------------------
!      1  |   8   |  7
!  -------+-------+----
!      2  |   0   |  6
!  -------+-------+-----
!      3  |   4   |  5
!  -----------------------------

         flag = 0

         if(w1 < 0 .and. w3 < 0) then
            flag = 1; w1 = 0; w2 = 1.0; w3 = 0; w4 = 0
         endif
         if(w3 < 0 .and. w4 < 0) then
            flag = 2; w3 = 0; w4 = 0
         endif
         if(w2 < 0 .and. w4 < 0) then
            flag = 3; w1 = 1.0; w2 = 0; w3 = 0; w4 = 0
         endif
         if(w2 < 0 .and. w3 < 0) then
            flag = 4; w2 = 0; w3 = 0
         endif
         if(w1 < 0 .and. w3 < 0) then
            flag = 5; w1 = 1.0; w2 = 0; w3 = 0; w4 = 1.0
         endif
         if(w1 < 0 .and. w2 < 0) then
            flag = 6; w1 = 0; w2 = 0
         endif
         if(w2 < 0 .and. w4 < 0) then
            flag = 7; w1 = 1.0; w2 = 0; w3 = 1.0; w4 = 0
         endif
         if(w1 < 0 .and. w4 < 0) then
            flag = 8; w1 = 0; w4 = 0
         endif

!--------------------------------------------------------------------------------

! --- Adjust for normalization.

         w1 = w1 * m1; w2 = w2 * m2; w3 = w3 * m3; w4 = w4 * m4;
         tot_w = w1 + w2 + w3 + w4;
         if(tot_w > 0 .and. tot_w /= 1) then
            w1 = w1 / tot_w; w2 = w2 / tot_w; 
            w3 = w3 / tot_w; w4 = w4 / tot_w
         endif

!  if (m1+m2+m3+m4 /= 0) then
!    write(20, *)'v ', var_coarse(i1  , j1  ), var_coarse(i1  , j1+1), var_coarse(i1+1, j1+1), var_coarse(i1+1, j1  )
!    write(20, *)'a ', a1, a2, a3, a4 
!    write(20, *)'w=[', w1, w2, w3, w4, '];' 
!    write(20, *)'m=[', m1, m2, m3, m4, '];' 
!    write(20, *)'x=[', x1, x2, x3, x4, '];' 
!    write(20, *)'y=[', y1, y2, y3, y4 , '];'
!    write(20, *)'xy=[', x, y, '];' !, i, j, i1, j1 
!    write(20, *)'st=[', s, t, '];' !, i, j, i1, j1 
!  endif

         return

!--------------------------------------------------------------------------------

         tot_m = m1 + m2 + m3 + m4; 

!         if(tot_m < 4) then
!            w1 = 0; w2 = 0; w3 = 0; w4 = 0
!         endif

!--------------------------------------------------------------------------------

         return

         
       end subroutine bilinear_coef_mask

!--------------------------------------------------------------------------------
       end module interp
  
