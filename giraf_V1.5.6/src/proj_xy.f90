subroutine proj_xy(proj,nx_fin,ny_fin,lon_fin,lat_fin) 

  implicit none
    
  !*****************************************************************************************
  ! projection des fichiers d'entree
  character(len=12) :: proj
  
  integer :: nx_fin,ny_fin

  real :: lon_fin(nx_fin,ny_fin),lat_fin(nx_fin,ny_fin) ! coordonnees grille fine  

  ! PROJECTION (UTILE dans version > 1.54)
  integer :: i, j
  double precision  :: x1, y1, h1, x2, y2,  h2
  integer           :: gs1
  integer           :: gs2  
  integer          :: projtype 
  double precision :: utmzone 
  double precision,parameter  :: d2r=.0174532925199  ! Degrees	->Radians
  double precision,parameter  :: r2d=57.2957795131  ! Radians	->Degrees
  double precision, parameter :: projunits = 1. !metres

!Projection des entrées
  projtype=1 ! par defaut Lambert 2 etendu - NTF
  utmzone=31.
  if (index(proj,'utm')  .ne.0)projtype=2
  if (index(proj,'utm31').ne.0)utmzone=31.
  if (index(proj,'utm32').ne.0)utmzone=32.
  if (index(proj,'l2')   .ne.0)projtype=1
  
  gs1 = 1 !WGS84
  if (projtype.eq.0)gs2 = 1 !geographique - WGS84
  if (projtype.eq.1)gs2 = 3 !Lambert 2 etendu - NTF
  if (projtype.eq.2)gs2 = 1 !Transverse Mercator (UTM) - WGS84  
  if (projtype.eq.3)gs2 = 4 !Lambert 93 - GRS80

  do i=1,nx_fin
    do j=1,ny_fin
 
  !Besoin de transformation
  if ((lon_fin(i,j).gt.180).or.(lat_fin(i,j).gt.90.)) then !    
    
    x2=lon_fin(i,j)
    y2=lat_fin(i,j)
    
    call projectinv2(gs2,x2 ,y2 ,utmzone, x1, y1, projtype,projunits)	    
    h1 = 0.
    call georef(gs1,x1,y1,h1,gs2,x1,y1,h1)
    lon_fin(i,j) = x1 * r2d
    lat_fin(i,j) = y1 * r2d   
  end if 
  
    end do
  end do 
  

END subroutine proj_xy

