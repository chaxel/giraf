!-------------------------------------------------------------
      integer function minint(a,b)
      
      integer :: a,b
      
      if (a.le.b) then
        minint=a
      else
        minint=b
      end if	 
      
      return
         
      end function
!-------------------------------------------------------------
      subroutine interp3d(nxi,nyi,nzi,nxo,nyo,nzo,wx,wy,ix,iy,in3d,out3d)
      
      implicit none
      
      integer :: i, j, k
      integer :: nxi,nyi,nzi      
      integer :: nxo,nyo,nzo
      integer :: ix(nxo,nyo),iy(nxo,nyo)
      real    :: wx(nxo,nyo),wy(nxo,nyo)      
      real    :: in3d (nxi,nyi,nzi)
      real    :: out3d(nxo,nyo,nzo)
            
      do i=1,nxo
       do j=1,nyo
        do k=1,nzo
         out3d(i ,j, k ) =(1.-wx(i ,j ))*(1.-wy(i ,j ))* in3d   (ix(i ,j )  ,iy(i ,j )  , k) &
        	        +     wx(i ,j ) *(1.-wy(i ,j ))* in3d   (ix(i ,j )+1,iy(i ,j )  , k) &
        	        + (1.-wx(i ,j ))*    wy(i ,j ) * in3d   (ix(i ,j )  ,iy(i ,j )+1, k) &
        	        +     wx(i ,j ) *    wy(i ,j ) * in3d   (ix(i ,j )+1,iy(i ,j )+1, k)        
        end do
       end do
      end do
	 
      end subroutine interp3d
!-------------------------------------------------------------
      subroutine interp2d(nxi,nyi,nxo,nyo,wx,wy,ix,iy,in2d,out2d)
      
      implicit none
      
      integer :: i, j
      integer :: nxi,nyi,nxo,nyo
      integer :: ix(nxo,nyo),iy(nxo,nyo)
      real    :: wx(nxo,nyo),wy(nxo,nyo)      
      real    :: in2d (nxi,nyi)
      real    :: out2d(nxo,nyo)
      
      do i=1,nxo
       do j=1,nyo
         out2d(i ,j ) =	(1.-wx(i ,j ))*(1.-wy(i ,j ))* in2d   (ix(i ,j )  ,iy(i ,j )  ) &
        	      +     wx(i ,j ) *(1.-wy(i ,j ))* in2d   (ix(i ,j )+1,iy(i ,j )  ) &
        	      + (1.-wx(i ,j ))*    wy(i ,j ) * in2d   (ix(i ,j )  ,iy(i ,j )+1) &
        	      +     wx(i ,j ) *    wy(i ,j ) * in2d   (ix(i ,j )+1,iy(i ,j )+1)	     
        end do
       end do
	 
      end subroutine interp2d      
!-------------------------------------------------------------      
      subroutine recopy3d(nxi,nyi,nzi,nxo,nyo,nzo,ix,iy,in3d,out3d)
      
      implicit none
      
      integer :: i, j, k
      integer :: nxi,nyi,nzi      
      integer :: nxo,nyo,nzo
      integer :: ix(nxo,nyo),iy(nxo,nyo)    
      real    :: in3d (nxi,nyi,nzi)
      real    :: out3d(nxo,nyo,nzo)
            
      do i=1,nxo
       do j=1,nyo
        do k=1,nzo
         out3d(i ,j, k ) = in3d(ix(i,j),iy(i,j),k)      
        end do
       end do
      end do
	 
      end subroutine recopy3d
!-------------------------------------------------------------
      subroutine recopy2d(nxi,nyi,nxo,nyo,ix,iy,in2d,out2d)
      
      implicit none
      
      integer :: i, j
      integer :: nxi,nyi,nxo,nyo
      integer :: ix(nxo,nyo),iy(nxo,nyo)     
      real    :: in2d (nxi,nyi)
      real    :: out2d(nxo,nyo)
      
      do i=1,nxo
       do j=1,nyo
         out2d(i ,j ) =	in2d(ix(i,j),iy(i,j))    
       end do
      end do
	 
      end subroutine recopy2d      
!-------------------------------------------------------------
      subroutine smoothing2d(nx,ny,in2d,out2d,ratio)
      
      !moyenne glissante d'une grille haute resolution vers une grille haute resolution
      
      implicit none
      
      integer :: i1, j1
      integer :: nx,ny
      integer :: ix(nx,ny),iy(nx,ny)     
      real    :: in2d (nx,ny)
      real    :: out2d(nx,ny)
      
      integer :: imin, jmin
      integer :: imax, jmax          
      
      integer :: i2, j2
      
      integer :: ratio
      
      integer :: radius
            
      real :: nval
      
      out2d = 0.
      
      radius = int((ratio-1.+.1)/2.)
            
      do i1=1,nx
       do j1=1,ny
         
	 imin = i1 - radius
	 jmin = j1 - radius
	 imax = i1 + radius
	 jmax = j1 + radius
	 
	 if (imin.lt. 1)imin=1
	 if (jmin.lt. 1)jmin=1
	 if (imax.gt.nx)imax=nx
	 if (jmax.gt.ny)jmax=ny	 
	 
	 nval = 0.
	 
	 do i2 = imin, imax
	 do j2 = jmin, jmax	 	 
	   nval = nval + 1.
           out2d(i1 ,j1 ) = out2d(i1 ,j1 ) + in2d(i2 ,j2 )	 
	 end do
	 end do
	 	 
	 out2d(i1 ,j1 ) = out2d(i1 ,j1 ) / nval
	 	 
       end do
      end do
	 
      end subroutine smoothing2d 
!-------------------------------------------------------------    
      subroutine somme_var(nvar,nxi,nyi,nxo,nyo,ix,iy,in2d,out2d)
      
      ! somme les valeurs d'une grille fine vers une grille grossière
      
      implicit none
      
      
      integer :: i, j, io, jo
      integer :: nvar, nxi,nyi,nxo,nyo
      integer :: ix(nxi,nyi),iy(nxi,nyi)     
      real    :: in2d (nvar,nxi,nyi)
      real    :: out2d(nvar,nxo,nyo)      
      
      out2d = 0.
      do i=1,nxi
       do j=1,nyi
         io = ix(i,j)
         jo = iy(i,j)	 
         out2d(:, io, jo) = out2d(:,io, jo) + in2d(:,i,j)    
       end do
      end do
	 
      end subroutine somme_var       
!-------------------------------------------------------------    
      subroutine moyenne_2d(nxi,nyi,nxo,nyo,ix,iy,in2d,out2d)

      ! moyenne les valeurs d'une grille fine vers une grille grossière
      
      implicit none
      
      integer :: i, j, io, jo
      integer :: nxi,nyi,nxo,nyo
      integer :: ix(nxi,nyi),iy(nxi,nyi)     
      real    :: in2d (nxi,nyi)
      real    :: out2d(nxo,nyo)
      real    :: nval(nxo,nyo)
      
      ! somme
      nval = 0.
      out2d = 0.
      do i=1,nxi
       do j=1,nyi
         io = ix(i,j)
         jo = iy(i,j)	 
         out2d(io, jo) = out2d(io, jo) + in2d(i,j)  
	 nval(io,jo) =   nval(io,jo) + 1.
       end do
      end do

      ! moyenne
      do io=1,nxo
       do jo=1,nyo
       if (nval(io,jo).gt.0.) then
        out2d(io,jo)  = out2d(io,jo) / nval(io,jo)
       else
        out2d(io,jo)  = 0.0
       end if  
       end do
      end do       
            
      
      end subroutine moyenne_2d
!-------------------------------------------------------------        
