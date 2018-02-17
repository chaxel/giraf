subroutine wcal2d(nx_coa,ny_coa,lon_coa,lat_coa,nx_fin,ny_fin,lon_fin,lat_fin,ix,iy,wx,wy) 

  implicit none
  
  

  !*****************************************************************************************
  ! subroutine arguments
  integer,intent(in) :: nx_coa
  integer,intent(in) :: ny_coa
  integer,intent(in) :: nx_fin
  integer,intent(in) :: ny_fin 
  integer,intent(out) :: ix(nx_fin,ny_fin) 
  integer,intent(out) :: iy(nx_fin,ny_fin) 
  real,   intent(out) :: wx(nx_fin,ny_fin)
  real,   intent(out) :: wy(nx_fin,ny_fin)

  ! parameters
  integer,parameter :: ntry=100 

  ! local variables
  integer :: nsec
  integer :: nx,ny,i,j,ntx,nty
  integer :: iflag,ierr
  real :: a,a0,a1,b,b0,b1,d2
  real :: d2min
  
  double precision :: x, y, xx, yy

  real :: lon_coa(nx_coa,ny_coa),lat_coa(nx_coa,ny_coa) ! coordonnees grille large
  real :: lon_fin(nx_fin,ny_fin),lat_fin(nx_fin,ny_fin) ! coordonnees grille fine  

  double precision :: x_coa(nx_coa,ny_coa),y_coa(nx_coa,ny_coa) ! coordonnees grille large
  double precision :: x_fin(nx_fin,ny_fin),y_fin(nx_fin,ny_fin) ! coordonnees grille fine  

  ! PROJECTION (UTILE dans version 1.54)
  double precision  ::  lon2r1, lat2r1, lon2r2, lat2r2
  double precision  :: x1, y1, h1, x2, y2,  h2
  integer           :: gs1
  integer           :: gs2 
  integer, parameter:: projtype  = 1 ! Lambert 2 etendu - NTF
  double precision, parameter :: utmzone = 31.
  double precision,parameter  :: d2r=.0174532925199  ! Degrees	->Radians
  double precision,parameter  :: r2d=57.2957795131  ! Radians	->Degrees
  double precision, parameter :: projunits = 1. !metres
  
  gs1 = 1 !WGS84
  if (projtype.eq.0)gs2 = 1 !geographique - WGS84
  if (projtype.eq.1)gs2 = 3 !Lambert 2 etendu - NTF
  if (projtype.eq.2)gs2 = 1 !Transverse Mercator (UTM) - WGS84  
  if (projtype.eq.3)gs2 = 4 !Lambert 93 - GRS80

!Projection
  do i=1,nx_coa
    do j=1,ny_coa
      x1 = lon_coa(i,j) * d2r
      y1 = lat_coa(i,j) * d2r
   !PB a cette etape avec domaine non-carre (utilise la projection)
!     ****COMPUTE TRANSFORMATION****
!     PROJ1
!    write(*,*)'projectinv2(gs1,x1,y1,zone,x2,y2,proj1,div)'
!    call projectinv2(gs1,x1,y1,zone,x2,y2,proj1,div)
!    write(*,*) gs1,x1,y1,zone,x2,y2,proj1,div
!   CHANGE GEOID    
    h1 = 0.
    call georef(gs1,x1,y1,h1,gs2,x1,y1,h1)
    !write(*,*) gs1,x2,y2,h1,gs2,x2,y2,h2

    call projectfor2(gs2,x1,y1,x_coa(i,j) ,y_coa(i,j) ,utmzone,projtype,projunits)
    !write(*,*) gs2,x2,y2,x3,y3,zone,proj2,div   
    end do
  end do

  do i=1,nx_fin
    do j=1,ny_fin
      x1 = lon_fin(i,j) * d2r
      y1 = lat_fin(i,j) * d2r    
   !PB a cette etape avec domaine non-carre (utilise la projection)
!     ****COMPUTE TRANSFORMATION****
!     PROJ1
!    write(*,*)'projectinv2(gs1,x1,y1,zone,x2,y2,proj1,div)'
!    call projectinv2(gs1,x1,y1,zone,x2,y2,proj1,div)
!    write(*,*) gs1,x1,y1,zone,x2,y2,proj1,div
!   CHANGE GEOID    
    h1 = 0.
    call georef(gs1,x1,y1,h1,gs2,x1,y1,h1)
    !write(*,*) gs1,x2,y2,h1,gs2,x2,y2,h2

    call projectfor2(gs2,x1,y1,x_fin(i,j) ,y_fin(i,j) ,utmzone,projtype,projunits)
    !write(*,*) gs2,x2,y2,x3,y3,zone,proj2,div   
    end do
  end do
  !*****************************************************************************************

  do i=1,nx_fin
    do j=1,ny_fin
    
      iflag = 1 

      ! Search whether the point is within one of the polygons   
      ! decale de 1.e-5 pour combler les erreurs d'interpolation            
      x = x_fin(i,j) + 1. !metre
      y = y_fin(i,j) + 1. !metre
     
      do nx=1,nx_coa-1 
        do ny=1,ny_coa-1

           nsec = 0 
	   
           if(  (y_coa(nx+1,ny+0).ne.y_coa(nx+0,ny+0))       .and.  &
                (y.le.max(y_coa(nx+1,ny+0),y_coa(nx+0,ny+0))).and.  &
                (y.ge.min(y_coa(nx+1,ny+0),y_coa(nx+0,ny+0))) )  then
              a = (y_coa(nx+1,ny+0)-y) / (y_coa(nx+1,ny+0)-y_coa(nx+0,ny+0))                        
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(x_coa(nx+0,ny+0)-x) + (1-a)*(x_coa(nx+1,ny+0)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(  (y_coa(nx+1,ny+1).ne.y_coa(nx+1,ny+0))       .and.  & 
                (y.le.max(y_coa(nx+1,ny+1),y_coa(nx+1,ny+0))).and.  &
                (y.ge.min(y_coa(nx+1,ny+1),y_coa(nx+1,ny+0))) ) then
              a = (y_coa(nx+1,ny+1)-y) / (y_coa(nx+1,ny+1)-y_coa(nx+1,ny+0))
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(x_coa(nx+1,ny+0)-x) + (1-a)*(x_coa(nx+1,ny+1)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(  (y_coa(nx+0,ny+1).ne.y_coa(nx+1,ny+1))       .and.  &
                (y.le.max(y_coa(nx+0,ny+1),y_coa(nx+1,ny+1))).and.  &
                (y.ge.min(y_coa(nx+0,ny+1),y_coa(nx+1,ny+1))) ) then
              a = (y_coa(nx+0,ny+1)-y) / (y_coa(nx+0,ny+1)-y_coa(nx+1,ny+1))                        
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(x_coa(nx+1,ny+1)-x) + (1-a)*(x_coa(nx+0,ny+1)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(  (y_coa(nx+0,ny+0).ne.y_coa(nx+0,ny+1))       .and.  &
                (y.le.max(y_coa(nx+0,ny+0),y_coa(nx+0,ny+1))).and.  &
                (y.ge.min(y_coa(nx+0,ny+0),y_coa(nx+0,ny+1))) ) then
              a = (y_coa(nx+0,ny+0)-y) / (y_coa(nx+0,ny+0)-y_coa(nx+0,ny+1))                        
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(x_coa(nx+0,ny+1)-x) + (1-a)*(x_coa(nx+0,ny+0)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(x.eq.x_coa(nx+0,ny+0).and.y.eq.y_coa(nx+0,ny+0)) then
              nsec = 1
           endif
           if(x.eq.x_coa(nx+1,ny+0).and.y.eq.y_coa(nx+1,ny+0)) then
              nsec = 1
           endif
           if(x.eq.x_coa(nx+1,ny+1).and.y.eq.y_coa(nx+1,ny+1)) then
              nsec = 1
           endif
           if(x.eq.x_coa(nx+0,ny+1).and.y.eq.y_coa(nx+0,ny+1)) then
              nsec = 1
           endif
	   
           if(nsec.eq.1) then 
              ix(i,j) = nx 
              iy(i,j) = ny 
              go to 1111 
           endif

        enddo
      enddo
      write(*,'("* POINT : ",2I5,2F14.3," OUTSIDE DOMAIN (nsec=",I2,")")')i,j,x,y,nsec
!      print *,'nsec=',nsec
      print *,'x_coa=',minval(x_coa),'->',maxval(x_coa) 
      print *,'y_coc=',minval(y_coa),'->',maxval(y_coa) 
      print *,'lon_coa=',minval(lon_coa),'->',maxval(lon_coa) 
      print *,'lat_coc=',minval(lat_coa),'->',maxval(lat_coa)	    
      print *,'a=',a,'b=',b
      !ix(i) = nx
      !iy(i) = ny
      iflag = 0 
1111 continue 

      !  Search for weights                                                   

      if(iflag.eq.1) then 
        d2min = 1e20 
        nx = ix(i,j) 
        ny = iy(i,j) 
        do ntx=0,ntry 
           a0 = float(ntx)/ntry 
           a1 = 1. - a0 
           do nty=0,ntry 
              b0 = float(nty)/ntry 
              b1 = 1. - b0 
              xx =     a0*b0*x_coa(nx+1,ny+1) &
                   & + a1*b0*x_coa(nx  ,ny+1) &
                   & + a0*b1*x_coa(nx+1,ny  ) &
                   & + a1*b1*x_coa(nx  ,ny  )                               
              yy =     a0*b0*y_coa(nx+1,ny+1) &
                   & + a1*b0*y_coa(nx  ,ny+1) &
                   & + a0*b1*y_coa(nx+1,ny  ) &
                   & + a1*b1*y_coa(nx  ,ny  )                               
              d2 = (x-xx)*(x-xx)+(y-yy)*(y-yy) 
              if(d2.lt.d2min) then 
                 d2min = d2 
                 a  = a0 
                 b  = b0 
              endif
           enddo
        enddo
        wx(i,j)=a
        wy(i,j)=b
      else
        ix(i,j)=0
        iy(i,j)=0
        wx(i,j)=0
        wy(i,j)=0
      endif


    enddo
  enddo  

END subroutine wcal2d

