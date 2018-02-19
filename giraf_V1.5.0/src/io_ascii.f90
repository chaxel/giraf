!--------------------------------------------------------------------        
      subroutine read_mailles_region_ra_1km
      
      use params 
      
      implicit none 
      
      ! LOCAL
      integer :: i1, j1, i2, j2, ix, iy          

! Mailles region
!	 allocate(rregion(nx2,ny2))	
!	 allocate(idept_region(nx2,ny2,99))		 
!	 allocate(izone_region(nx2,ny2,13))	 	 	 	  	 
!	 allocate(pop_region(nx2,ny2))	

! Fichier disponible ? OUI DONNEES REGIONALES DISPONIBLES
      if ( ifregion ) then
         region1km=.true.
         if (idebug)write(*,*) 'Lit mailles [x_utm,y_utm,pop] region Rhone-Alpes dans fichier: '//trim(fregion)
         
	 allocate(iregion(nx2,ny2))
	 allocate(rregion(nx2,ny2))	
	 allocate(idept_region(nx2,ny2,99))		 
	 allocate(izone_region(nx2,ny2,13))	 	 	 	  	 
	 allocate(pop_region(nx2,ny2))	 
	 	 
         do i2=1,nx2
           do j2=1,ny2	  	
	      iregion(i2,j2) = .false.
	      rregion(i2,j2) = 0.
	      pop_region(i2,j2) = 0.	      
	      idept_region(i2,j2,:) = 0.
	      izone_region(i2,j2,:) = 0.	      	      	      
	    end do
          end do 	 
	 
	 open( unit=10,file=trim(fregion),status='old' )
	 
         do i2=1,nx2
           do j2=1,ny2	  	
	      read(10,*,end=999) i1, j1, pop !, zone_region, dept_region
	      ix =  (i1- 553500)/1000 + 1
	      iy =  (j1-4887500)/1000 + 1
	      if ( (ix .ge.1).and.(iy .ge.1).and.(ix .le.nx2).and.(iy .le.ny2) ) then 
	        iregion(ix,iy) = .true.
	        rregion(ix,iy) = 1.
		i1 = dept_region
		j1 = zone_region - 100
		pop_region(ix, iy) = pop
		!write(*,*) ix, iy, i1, j1, pop	
		if ((i1.ge.1).and.(i1.le.99))idept_region(ix, iy, i1) = 1.
		if ((j1.ge.1).and.(j1.le.13))izone_region(ix, iy, j1) = 1.		
	      end if        
	    end do
          end do   
	  
999      continue	  	 	 
	 close(10)
	 
! Fichier disponible ? NON PAS DE DONNEES REGIONALES
      else
        if (idebug)write(*,*) 'Ne lit pas les mailles region'   	   
        region1km=.false.
	iregion=.true.
	rregion=1.
	idept_region=0
	izone_region=0
	pop_region=0
	
      end if
      
      end subroutine
!--------------------------------------------------------------------                   
      subroutine read_ascii_ndim(file_ascii,nx,ny)
      
      implicit none  

      ! local 
      character(len=*) :: file_ascii
      integer :: ligne
      integer, parameter :: nmax = 10000000
      real    :: x2, x1
      real    :: y2, y1
      integer :: nx, ny  
      
      open(unit=10,file=file_ascii, status='old' )
      
      x1 = -1.E30
      y1 = -1.E30  
      nx = 0 
      ny = 0          
      do ligne=1, nmax
         x2 = x1
	 y2 = y1
         read(10,*,end=99) x1,y1
	 if (x1.lt.x2.and.nx.eq.0) nx = ligne - 1
      end do
      
      99 continue
      
      ny = (ligne - 1)/nx
      
      close(10)

      
      end subroutine    
!-----------------------------------------------------------------------------      
      subroutine read_chimere_ground_ascii(file_ascii, nx, ny, lon, lat, landuse, topo, idebug )

      use netcdf
      use typesizes
      
      implicit none  
      
      character(len=*) :: file_ascii      
      integer :: nx, ny
      real    :: lon(nx,ny)
      real    :: lat(nx,ny)
      real    :: landuse(nx,ny)
      real    :: topo(nx,ny)
      integer :: dxmin, dymin 
      
      integer :: i,j
      logical :: idebug     
       
      ! Grille CHIMERE
      if (idebug)write(*,*) 'Lit ASCII '//trim(file_ascii)
      open(unit=10,file=file_ascii, status='old' )      
      
      do j=1,ny
        do i=1, nx
	  read(10,*) lon(i,j), lat(i,j), topo(i,j)
	end do
      end do      	  
  
      end subroutine            
!-----------------------------------------------------------------------------      
      subroutine read_ncol_ascii(file_ascii, ncol, nx, ny, valeur, idebug )

      use netcdf
      use typesizes
      
      implicit none  
      
      character(len=*) :: file_ascii      
      integer :: ncol, nx, ny
      real    :: valeur(ncol,nx,ny)
      real    :: topo(nx,ny)
      integer :: dxmin, dymin 
      
      integer :: i,j
      logical :: idebug     
       
      ! Grille CHIMERE
      if (idebug)write(*,*) 'Lit '//trim(file_ascii)
      open(unit=10,file=file_ascii, status='old' )      
      
      do j=1,ny
        do i=1, nx
	  read(10,*,end=99) valeur(:,i,j)
	end do
      end do
      
      go to 98
      
99    continue
      write(*,*) 'ERREUR de lecture de '//trim(file_ascii)//' a la ligne',i+(j-1)*nx
      write(*,*) 'Le fichier doit comporter un nombre de colonnes egal a',ncol
      write(*,*) 'Le fichier doit comporter un nombre de lignes egal a',nx*ny        
      stop
      
98    continue            	  
  
      end subroutine 
!-----------------------------------------------------------------------------      
      subroutine read_x_y_ncol_ascii(file_ascii, ncol, nx, ny, x, y, valeur, idebug )

      use netcdf
      use typesizes
      
      implicit none  
      
      character(len=*) :: file_ascii      
      integer :: ncol, nx, ny
      real    :: valeur(ncol,nx,ny)
      real    :: x(nx,ny)
      real    :: y(nx,ny)      
      integer :: dxmin, dymin 
      
      integer :: i,j
      logical :: idebug     
       
      ! Grille CHIMERE
      if (idebug)write(*,*) 'Lit '//trim(file_ascii)
      open(unit=10,file=file_ascii, status='old' )      
      
      do j=1,ny
        do i=1, nx
	  read(10,*,end=99) x(i,j),y(i,j),valeur(:,i,j)
	end do
      end do
      
      go to 98
      
99    continue
      write(*,*) 'ERREUR de lecture de '//trim(file_ascii)//' a la ligne',i+(j-1)*nx
      write(*,*) 'Le fichier doit comporter un nombre de colonnes (X,Y,N) egal a',ncol
      write(*,*) 'Le fichier doit comporter un nombre de lignes egal a',nx*ny        
      stop
      
98    continue            	  
  
      end subroutine 
