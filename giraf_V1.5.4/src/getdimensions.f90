      subroutine getdimensions(fileName,numLons,numLats,numLays,&
                      numFrTimes,timeStringLen,&
                      distGrid,startYear,startMnth,startDay,&
		      startHour,startMin,startSec)

      use netcdf
      use typesizes
		
      implicit none
      		      
      ! General  
      character(len=256):: fileName   !Name of the input file              
      integer :: numLats	     !Dimension of grid in Y direction
      integer :: numLons	     !Dimension of grid in X direction
      integer :: numLays	     !Number of sigma-P levels
      integer :: numSteps	     !Step number
      integer :: numFrTimes	     !Numbers of time steps
      integer :: timeStringLen       !Duration of one time step (in minutes)
      real    :: distGrid	     !Distance of grid cells (in meters)
      integer :: startYear	     !Four-digit year of start time
      integer :: startMnth	     !Month of start time
      integer :: startDay	     !Day of start time
      integer :: startHour	     !Hour of start time
      integer :: startMin	     !Minute of start time
      integer :: startSec	     !Second of start time
      character(len=33) :: tunits
      		      
      !netcdf related
      integer :: ncFileID, &
                 latDimID, lonDimID, zDimID, frTimeDimID, timeDimID, &
                 latVarID, lonVarID, zVarID,frTimeVarID, scalarVarID, xvarID 
		 
      real,allocatable :: xCoord(:,:) ! X coordinates of the grid 
      real,allocatable :: tCoord(:)           	 		 
		 
      integer :: n
      
      character(len=19) :: datestr       		      

    write(*,*) 'Read existing CDF '//trim(fileName) 
     
    call check(nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID))       

    call check(nf90_inq_dimid(ncFileID, 'west_east'    , lonDimID))     
    call check(nf90_inq_dimid(ncFileID, 'south_north'  , latDimID))
    call check(nf90_inq_dimid(ncFileID, 'bottom_top'  , zDimID))  
    call check(nf90_inq_dimid(ncFileID, 'Time', frTimeDimID))  
                   
    call check(nf90_Inquire_Dimension(ncFileID, latDimID   , len=numLats))
    call check(nf90_Inquire_Dimension(ncFileID, lonDimID   , len=numLons))    
    call check(nf90_Inquire_Dimension(ncFileID, zDimID     , len=numLays))            
    call check(nf90_Inquire_Dimension(ncFileID, frTimeDimID, len= numFrTimes))    

! nombre de pas de temps
      call check(nf90_inq_dimid(ncFileID, 'Time', frTimeDimID)) 
      call check(nf90_Inquire_Dimension(ncFileID, frTimeDimID, len= numFrTimes))
      write(*,*) 'Trouve pas de temps : ',numFrTimes   

      ! date
      write(*,*) 'Lit date de CHIMERE'             	      
      call check(nf90_inq_varid(ncFileID, "Times", xVarID))    
      call check(nf90_get_var(ncFileID, xVarID, datestr, start = (/1,1/) ))
     	      
      read(datestr,'(I4,"-",I2,"-",I2,"_",I2,":",I2,":",I2)')&
     	     startYear,startMnth,startDay,startHour,startMin,startSec
      write(*,*) 'Trouve date initiale : '//datestr

    allocate(xcoord(numLons,numLats))
    
    call check(nf90_inq_varid(ncFileID, 'lon', lonVarID))
    call check(nf90_get_var(ncFileID, lonVarID, xcoord, start = (/ 1,1/) ))
    distgrid= xcoord(2,1)-xcoord(1,1)
    
    write(*,*) 'distgrid(degrees)=',distgrid

    deallocate(xcoord)
        
    call check(nf90_close(ncFileID)) 
               
      end subroutine
