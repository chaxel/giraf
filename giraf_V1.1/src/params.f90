
      module params

      ! PARAMETRES PROGRAMMES
      character(len=20), parameter :: exe_str='giraf.exe'
      real             , parameter :: version = 1.10
      integer          , parameter :: fnstrl=256
      
      ! PARAMETRES SUR CALCUL
      logical :: iraf_alti =.true.
      logical :: iraf_chim =.true.
      
      ! MODE DEBOGAGE     
      logical :: idebug ! en entree par defaut .false.

      ! PRESENCE DES FICHIERS D'ENTREE ?
      logical :: iout1           ! obligatoire
      logical :: iground2,iemis2 ! obligatoire si rafinement
      logical :: iground1,iemis1 ! facultatif (recalculé)
      logical :: ifregion        ! facultatif (recalculé)      
      logical :: iout2,imet      ! completement facultatif

      ! FICHIERS d'ENTREE AU FORMAT NETCDF ?
      logical  :: inc_topo, inc_emis

      ! parametres grille
      real, parameter  :: distGrid_g2 = 1.	  !Taille de la maille de la grille FINE (par defaut 1 km) 
      
      ! parametres emissions/variables    
      integer, parameter    :: nvar = 9   ! utilise 3 polluants gaz O3, NO, NO2          PMs :  PM10  pPPM  p10HNO3
      character(len=fnstrl) :: var_name(nvar)        = (/ 'O3' ,   'NO'  ,  'NO2' ,     'PM10',     'PM25',     'pPPM',    'pHNO3', 'temp'  ,  'wins'  /)
      character(len=fnstrl) :: var_unit(nvar)        = (/ 'ppb',   'ppb' ,  'ppb' ,'microg/m3','microg/m3','microg/m3','microg/m3', 'degC'  ,  'm/s'   /)     
      character(len=fnstrl) :: var_name_cdf(nvar)   
      integer, parameter    :: voisinage_xy_var(nvar)= (/    10,       10,      10,         10, 	10,          0,         10,     10  ,      10  /) ! taille des demi-voisinages (1 km)

      logical               :: iinput_var(nvar)      ! la variable est disponible dans le NetCDF ? 
      logical               :: icorrect_z_var(nvar)  = (/ .true.,  .true.,  .true.,      .true.,    .true.,    .false.,     .true.,  .true. ,   .true. /) ! corrige O3/NO2/PM10 avec profils Z ?    
      logical               :: ioutput_var(nvar)     = (/ .true.,  .true.,  .true.,      .true.,    .true.,    .false.,     .true.,  .true. ,   .true. /) ! ecrit la variables dans NetCDF ?
       ! DEBUG
      logical               :: ioutput_dz_var(nvar)  = (/ .true.,  .true.,  .true.,      .true.,   .false.,    .false.,     .true.,   .true.,   .true. /) ! ecrit la variables dans NetCDF (DEBUG mode) ?
      logical               :: ioutput_de_var(nvar)  = (/ .true.,  .true.,  .true.,      .true.,   .false.,    .false.,    .false.,  .false.,   .false./) ! ecrit la variables dans NetCDF (DEBUG mode) ?

      logical               :: ioutput_krig_var(nvar)=.false. ! (/ .false., .false., .false., .false., .false.,.false., .false., .false./) ! creer une variable *_KRIG dans NetCDF ?         
   
      integer, parameter    :: nvar_emis  = 2     ! utilise 3 polluants NOX, VOC, PM10        
      character(len=fnstrl) :: var_emis_name(nvar_emis)  = (/'NOX','PM10'/)  ! ne pas modifier  
      integer, parameter    :: nfile_emis = 1     ! utilise 1 fichiers d entree d'emissions 
      integer, parameter    :: nsnap1 = 11        ! nombre de classe snap niveau 1 (par defaut 11)
      ! Tableau de reparition par SNAP (fraction d'emission dans le niveau 1 du modele CHIMERE)
      !                        S01   S02   S03   S04   S05   S06   S07   S08   S09   S10   S11 
      real, parameter       :: tbl_ratio_emi_snap(nsnap1) = (/ 0.00, 0.20, 0.20, 0.20, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.00 /)                 
     
      integer, parameter    :: nvar_ppm   = 3     ! utilise 3 especes primaires PPM        
      character(len=fnstrl) :: var_ppm_name_cdf(nvar_ppm)    != (/'p10PPM','p10OCAR','p10BCAR'/)  ! ne pas modifier       
       
      ! CORRECTIONS                                      						 
      logical, parameter :: correct_no    = .true.
      logical, parameter :: correct_o3    = .false.
      logical, parameter :: correct_no2   = .false.
      logical, parameter :: correct_pm10  = .true.         
      
      ! INDICES DES POLLUANTS POUR LA CHIMIE
      integer :: chim_ino
      integer :: chim_ino2
      integer :: chim_io3
      integer :: chim_ipm10
      integer :: chim_ippm
      integer :: chim_iphno3
      integer :: met_itemp
      
      ! INDICES DES ESPECES POUR LA METEO
      integer :: met_iwinm
      integer :: met_iwinz
      integer :: met_iwins
      
      ! INDICES DES POLLUANTS POUR LES EMISSINS
      integer :: emis_icovnm
      integer :: emis_inox
      integer :: emis_ipm10   
      
      ! NOM DES FICHIERS          
      character(len=fnstrl) :: fout1, fdep1, fmet1, fground1, fregion
      character(len=fnstrl) :: femis1(nfile_emis)

      character(len=fnstrl) :: fout2, fground2
      character(len=fnstrl) :: femis2(nfile_emis)      

      ! Grille CHIMERE 
      integer :: nxc,nyc,nzc    ! domaine CHIMERE
      integer :: nx0,ny0,nz0    ! domaine CHIMERE coarse  
      integer :: nx1,ny1,nz1    ! domaine CHIMERE raffine
      integer :: nx2,ny2,nz2    ! domaine FINE
      integer :: dx1min, dx1max  ! bornes X du domaine reduit
      integer :: dy1min, dy1max  ! bornes Y du domaine reduit
      
      ! LOGICAL SUR DONNEES CHIMERE                      
      logical :: ipm, ippm, iwin
      logical :: imeteo
      logical :: iinterp

      ! PRINT profils
      logical :: iprint 
      integer :: ixp,iyp     !coordonnee du point d'affichage   
      integer :: itp = 1     !coordonnee du point d'affichage                                
      
      ! PHYSIQUE
      REAL :: es , qs
      REAL , PARAMETER :: svp1        =     0.6112
      REAL , PARAMETER :: eps         =     0.622      
      REAL , PARAMETER :: svp2        =    17.67
      REAL , PARAMETER :: svp3        =    29.65
      REAL , PARAMETER :: svpt0       =   273.15
      REAL , PARAMETER :: rayon_km    =   6340.0
      
      ! METEO
      real, parameter   :: t0k = 273.15       

      ! PROJECTION
      double precision  :: h1, h2, lon2r1, lat2r1, lon2r2, lat2r2, x2, y2, projunits
      integer           :: projtype=2 ! UTM
      double precision  :: utmzone=31
      integer           :: gs1, gs2
      double precision  :: d2r	 ! Degrees   ->Radians	 
      double precision  :: r2d	 ! Radians   ->Degrees	 	     
      parameter( d2r=.0174532925199	   )
      parameter( r2d=57.2957795131	   ) 

      ! Grille CHIMERE brute
      real,allocatable :: lon_gc(:,:)
      real,allocatable :: lat_gc(:,:)
      real,allocatable :: topo_gc(:,:)
      real,allocatable :: landuse_gc(:,:)
            
      ! Grille COARSE 
      real,allocatable :: lon_g0(:,:)
      real,allocatable :: lat_g0(:,:)
      real,allocatable :: topo_g0(:,:)
      real,allocatable :: landuse_g0(:,:)
      real,allocatable :: hght_g0(:,:)
      real,allocatable :: tem2_g0(:,:) 
      real,allocatable :: w10m_g0(:,:)            
      real,allocatable :: alti_agl_g0(:,:,:) !3D 
                        
      ! Grille FINE interpolee      
      real,allocatable :: lon_g1(:,:)
      real,allocatable :: lat_g1(:,:)
      real,allocatable :: topo_g1(:,:)
      real,allocatable :: landuse_g1(:,:)

      ! variables meteo        
      real,allocatable    :: hght_g1(:,:)
      real,allocatable    :: tem2_g1(:,:)
      real,allocatable    :: w10m_g1(:,:)      
      real,allocatable    :: alti_asl_g1(:,:,:)
      real,allocatable    :: alti_agl_g1(:,:,:) 
      real,allocatable    :: temp_g1(:,:,:)
      integer,allocatable :: iz_alti_g1(:,:,:)
      real,allocatable    :: i_alti_g1(:,:,:)        
               
      ! Grille FINE rafinée
      real,allocatable    :: lon_g2(:,:)
      real,allocatable    :: lat_g2(:,:)
      real,allocatable    :: topo_g2(:,:)
      real,allocatable    :: landuse_g2(:,:)
      real,allocatable    :: hght_g2(:,:)      
      real,allocatable    :: tem2_g2(:,:)
      real,allocatable    :: w10m_g2(:,:)                    

      ! especes gazeuses (pour l'instant O3 seulement)
      real,allocatable :: val_g0_var (:,:,:,:)
      real,allocatable :: val_g1_var (:,:,:,:)          
      real,allocatable :: val_g2_var (:,:,:,:)   

      ! FDMS
      real,allocatable :: pm10_nonvolat_g2(:,:,:)
      real,allocatable :: ecart_fdms_g2(:,:,:)

      ! Profils 1D + 3D   
      real,allocatable :: alti_1D(:)             ! Profil d'altitude pour le calcul (orthonorme)
      real,allocatable :: profil_1D_var(:,:)     ! Profil moyen espece N                
      real,allocatable :: profil_3D_var(:,:,:,:) ! Profil moyen d'O3      
      real,allocatable :: dvar_dz_loc(:,:,:)     ! gradients local altitude
      real,allocatable :: dvar_de_loc(:,:,:)     ! gradients local traitements chimiques 
  
      ! interpolation grille CHIMERE -> grille fine (wcal1.f90)
      integer,allocatable :: ix2 (:,:)
      integer,allocatable :: iy2 (:,:)      
      integer,allocatable :: ix1 (:,:)
      integer,allocatable :: iy1 (:,:)       
      real   ,allocatable :: wx1 (:,:)
      real   ,allocatable :: wy1 (:,:)        
          
      ! General inputs     
      integer :: numSteps	     !Step number
      integer :: numFrTimes	     !Numbers of time steps
      integer :: timeStringLen       !Duration de one time step (in minutes)
      real    :: distGrid_dg	     !Distance de grid cells (in degrees)
      real    :: distGrid_km	     !Distance de grid cells (in km)  
        
      integer :: Year	     !Four-digit year of obs time
      integer :: Mnth	     !Month of obs time
      integer :: Day	     !Day of obs time
      integer :: Hour	     !Hour of obs time
      integer :: Min	     !Minute of obs time
      integer :: Sec	     !Second of obs time 
      integer :: nday        !number of day in a month
      integer :: it          !timestep of observations in mmoutcdf file  

      ! MAILLES REGION
      ! information sur la region rhone-alpes      
      logical :: region1km
      logical,allocatable :: iregion(:,:)
      real,allocatable    :: rregion(:,:)        
      real                :: dept_region, zone_region,pop
      real,allocatable    :: izone_region(:,:,:) 
      real,allocatable    :: idept_region(:,:,:) !1,7,26,38,73,74,69,42           
      real,allocatable    :: pop_region(:,:)         

      ! EMISSIONS 1=NOX, 2=VOC, 3=PM10
      real,allocatable :: emis_g0(:,:,:)        ! (id_polluant,x,y)
      real,allocatable :: emis_g2(:,:,:)        ! (id_polluant,x,y)
      real,allocatable :: ratio_emis_g2(:,:,:)  ! rapport des emissions maille fine/maille chimere (id_polluant,x,y)
      real,allocatable :: ratio_emis_g0(:,:,:)  ! rapport des emissions maille fine/maille chimere (id_polluant,x,y)
      real,allocatable :: nratio_emis_g0(:,:,:) ! rapport des emissions maille fine/maille chimere (id_polluant,x,y)

      ! NETCDF
      integer :: met1FileID, out1FileID, out2FileID
  
      end module params
