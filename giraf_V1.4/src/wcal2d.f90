subroutine wcal2d(nxx,nyy,ix,iy,wx,wy,xc,yc,nx2,ny2,x2,y2) 

  implicit none

  !*****************************************************************************************
  ! subroutine arguments
  integer,intent(in) :: nxx
  integer,intent(in) :: nyy
  integer,intent(in) :: nx2
  integer,intent(in) :: ny2  
  integer,intent(out) :: ix(nx2,ny2) 
  integer,intent(out) :: iy(nx2,ny2) 
  real,   intent(out) :: wx(nx2,ny2)
  real,   intent(out) :: wy(nx2,ny2)

  ! parameters
  integer,parameter :: ntry=100 

  ! local variables
  integer :: nsec
  integer :: nx,ny,i,j,ntx,nty
  integer :: iflag,ierr
  real :: x,y,xx,yy,a,a0,a1,b,b0,b1,d2
  real :: d2min

  real :: xc(nxx,nyy),yc(nxx,nyy) ! coordonnees grille large
  real :: x2(nx2,ny2),y2(nx2,ny2) ! coordonnees grille fine  

  !*****************************************************************************************

  do i=1,nx2
    do j=1,ny2
    
      iflag = 1 

      ! Search whether the point is within one of the polygons   
      ! decale de 1.e-5 pour combler les erreurs d'interpolation            
      x = x2(i,j) + 1.e-5
      y = y2(i,j) + 1.e-5
     
      do nx=1,nxx-1 
        do ny=1,nyy-1 
           nsec = 0 
	   
           if(  (yc(nx+1,ny+0).ne.yc(nx+0,ny+0))       .and.  &
                (y.le.max(yc(nx+1,ny+0),yc(nx+0,ny+0))).and.  &
                (y.ge.min(yc(nx+1,ny+0),yc(nx+0,ny+0))) )  then
              a = (yc(nx+1,ny+0)-y) / (yc(nx+1,ny+0)-yc(nx+0,ny+0))                        
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(xc(nx+0,ny+0)-x) + (1-a)*(xc(nx+1,ny+0)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(  (yc(nx+1,ny+1).ne.yc(nx+1,ny+0))       .and.  & 
                (y.le.max(yc(nx+1,ny+1),yc(nx+1,ny+0))).and.  &
                (y.ge.min(yc(nx+1,ny+1),yc(nx+1,ny+0))) ) then
              a = (yc(nx+1,ny+1)-y) / (yc(nx+1,ny+1)-yc(nx+1,ny+0))
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(xc(nx+1,ny+0)-x) + (1-a)*(xc(nx+1,ny+1)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(  (yc(nx+0,ny+1).ne.yc(nx+1,ny+1))       .and.  &
                (y.le.max(yc(nx+0,ny+1),yc(nx+1,ny+1))).and.  &
                (y.ge.min(yc(nx+0,ny+1),yc(nx+1,ny+1))) ) then
              a = (yc(nx+0,ny+1)-y) / (yc(nx+0,ny+1)-yc(nx+1,ny+1))                        
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(xc(nx+1,ny+1)-x) + (1-a)*(xc(nx+0,ny+1)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(  (yc(nx+0,ny+0).ne.yc(nx+0,ny+1))       .and.  &
                (y.le.max(yc(nx+0,ny+0),yc(nx+0,ny+1))).and.  &
                (y.ge.min(yc(nx+0,ny+0),yc(nx+0,ny+1))) ) then
              a = (yc(nx+0,ny+0)-y) / (yc(nx+0,ny+0)-yc(nx+0,ny+1))                        
              if(a.ge.0.and.a.lt.1) then 
                 b =    a *(xc(nx+0,ny+1)-x) + (1-a)*(xc(nx+0,ny+0)-x)
                 if(b.gt.0) nsec = nsec + 1 
              endif
           end if

           if(x.eq.xc(nx+0,ny+0).and.y.eq.yc(nx+0,ny+0)) then
              nsec = 1
           endif
           if(x.eq.xc(nx+1,ny+0).and.y.eq.yc(nx+1,ny+0)) then
              nsec = 1
           endif
           if(x.eq.xc(nx+1,ny+1).and.y.eq.yc(nx+1,ny+1)) then
              nsec = 1
           endif
           if(x.eq.xc(nx+0,ny+1).and.y.eq.yc(nx+0,ny+1)) then
              nsec = 1
           endif
	   
           if(nsec.eq.1) then 
              ix(i,j) = nx 
              iy(i,j) = ny 
              go to 1111 
           endif

        enddo
      enddo
      print *,'* POINT : ',i,j,x,y,' OUTSIDE DOMAIN'
      print *,'nsec=',nsec
      print *,'xc=',minval(xc),'->',maxval(xc) 
      print *,'yc=',minval(yc),'->',maxval(yc) 
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
              xx =     a0*b0*xc(nx+1,ny+1) &
                   & + a1*b0*xc(nx  ,ny+1) &
                   & + a0*b1*xc(nx+1,ny  ) &
                   & + a1*b1*xc(nx  ,ny  )                               
              yy =     a0*b0*yc(nx+1,ny+1) &
                   & + a1*b0*yc(nx  ,ny+1) &
                   & + a0*b1*yc(nx+1,ny  ) &
                   & + a1*b1*yc(nx  ,ny  )                               
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

