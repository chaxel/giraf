!--------------------------------------------------------------------        
      subroutine calcul_domaine_reduit(    &
      nxc,nyc,nzc,lonc,latc,nx2,ny2,nz2,lon2,lat2, &
      nx0,ny0,nz0,dx1min,dx1max,dy1min,dy1max,idebug)
      
      implicit none
      
      ! entrees
      integer :: nxc,nyc,nzc
      integer :: nx2,ny2,nz2
      real    :: lonc(nxc,nyc),latc(nxc,nyc)
      real    :: lon2(nx2,ny2),lat2(nx2,ny2)
      logical :: idebug
      
      ! sorties
      integer :: nx0,ny0,nz0,dx1min,dx1max,dy1min,dy1max  
      
      ! local 
      integer :: i1,j1
      
      dx1min = 1      
      dx1max = nxc
      dy1min = 1      
      dy1max = nyc  
          
      do i1=1,nxc      
       do j1=1,nyc         
        if ( ( lonc(i1,j1) .gt. maxval( lon2 ) ) .and. ( i1 .lt. dx1max ) )dx1max = i1	
        if ( ( lonc(i1,j1) .lt. minval( lon2 ) ) .and. ( i1 .gt. dx1min ) )dx1min = i1		
        if ( ( latc(i1,j1) .gt. maxval( lat2 ) ) .and. ( j1 .lt. dy1max ) )dy1max = j1
        if ( ( latc(i1,j1) .lt. minval( lat2 ) ) .and. ( j1 .gt. dy1min ) )dy1min = j1			 	
       end do     
      end do   
      
      ! On joue la securite en ajoutant 2 mailles
      dx1min = max( 1, dx1min - 2 )       
      dy1min = max( 1, dy1min - 2 )      
      dx1max = dx1max + 2
      if ( dx1max.ge. nxc ) dx1max = nxc
      dy1max = dy1max + 2
      if ( dy1max.ge. nyc ) dy1max = nyc           
      
      if (idebug)write(*,*)'Domaine réduit'
      if (idebug)write(*,'("dx1min, dx1max, dy1min, dy1max",4I4)') dx1min, dx1max, dy1min, dy1max
      
      if (idebug)write(*,*) 'Domaine sortie'       
      if (idebug)write(*,*) 'lon2=',minval(lon2),maxval(lon2)
      if (idebug)write(*,*) 'lat2=',minval(lat2),maxval(lat2)  
      if (idebug)write(*,*) 'Domaine entree brut'       
      if (idebug)write(*,*) 'lonc=',minval(lonc),maxval(lonc)
      if (idebug)write(*,*) 'latc=',minval(latc),maxval(latc)
      if (idebug)write(*,*) 'Domaine entree reduit'                          
      if (idebug)write(*,*) 'lon0=',minval(lonc(dx1min:dx1max,dy1min:dy1max)),maxval(lonc(dx1min:dx1max,dy1min:dy1max))
      if (idebug)write(*,*) 'lat0=',minval(latc(dx1min:dx1max,dy1min:dy1max)),maxval(latc(dx1min:dx1max,dy1min:dy1max))         
      nx0 = dx1max - dx1min + 1     
      ny0 = dy1max - dy1min + 1
      nz0 = nzc 
      
      end subroutine  calcul_domaine_reduit
!--------------------------------------------------------------------
      subroutine calcul_correspondance
     
      use params
      
      implicit none
        
      integer :: i0, j0, i2, j2
      real    :: rsearch, rtemp
      
      open(unit=10,file='WCROSS',status='unknown')
      
      do i2=1,nx2      
        do j2=1,ny2
		  
          ix2(i2,j2) = 0	  
          iy2(i2,j2) = 0	
	  
	  rsearch = 1.E30 !distGrid_dg * 5
	  
          do i0=1,nx0
            do j0=1,ny0
	      rtemp = ((lon_g0(i0,j0)*100.-lon_g2(i2,j2)*100.)**2+(lat_g0(i0,j0)*100.-lat_g2(i2,j2)*100.)**2)**.5
	      if (rtemp.lt.rsearch) then
                ix2(i2,j2) = i0
                iy2(i2,j2) = j0
	        rsearch=rtemp	          
	      end if 	  
	    end do
          end do 
	       
          if (ix2(i2,j2)*iy2(i2,j2).eq.0) then
            write(*,*) 'Point en dehors du domaine [lon lat]:',lon_g2(i2,j2),lat_g2(i2,j2)	    
            stop
	  else
	    write(10,'(2I6,2F16.5,2I7)') i2,j2,lon_g2(i2,j2),lat_g2(i2,j2),ix2(i2,j2),iy2(i2,j2)	  
          end if
	  	  	
	end do
	
      end do
      
      close(10)
      
      if (idebug)write(*,*) '-> Ecrit WCROSS'

      end subroutine  
!--------------------------------------------------------------------
      subroutine calcul_poids

      use params
      
      implicit none

      ! LOCAL
      integer :: i1, j1, i2, j2
           
      call wcal2d(nx0,ny0,lon_g0,lat_g0,nx2,ny2,lon_g2,lat_g2,ix1,iy1,wx1,wy1)
 
      open(unit=10,file='WDOT',status='unknown')
      do i2=1,nx2
        do j2=1,ny2    
	  write(10,*) ix1(i2,j2),iy1(i2,j2),wx1(i2,j2),wy1(i2,j2)
	end do
      end do		
      close(10)
      
      if (idebug)write(*,*) '-> Ecrit WDOT'      
 
      end subroutine  calcul_poids
!--------------------------------------------------------------------
