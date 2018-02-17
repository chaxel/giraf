      module params_cdf

      integer           :: ncFileID, xVarID
      character(len=19) :: datestr 
      integer           :: ivar   
       
      end module params_cdf
!---------------------------------------------------------------------------------------------------------------------   
      subroutine check_cdf_var
      
      use netcdf
      use typesizes

      use params
      use params_cdf
      
      implicit none
      
      if (idebug)write(*,*) 'Verifie la presence des variables dans CHIMERE'    
      call check(nf90_open(fout1, nf90_nowrite, ncFileID)) 
          
      ! initialisation
      iinput_var(nvar)=.false.
      var_name_cdf = var_name
      
      do ivar=1, nvar
	if (nf90_inq_varid(ncFileID, trim(var_name(ivar)), xVarID).eq.nf90_noerr) then
	  var_name_cdf(ivar)=var_name(ivar)
	  iinput_var(ivar)=.true.
	end if
      end do
        	 
      do ivar=1, nvar	                       
        if (trim(var_name(ivar)).eq.'NO'    )chim_ino    = ivar
        if (trim(var_name(ivar)).eq.'NO2'   )chim_ino2   = ivar
        if (trim(var_name(ivar)).eq.'O3'    )chim_io3    = ivar
        if (trim(var_name(ivar)).eq.'PM10'  )chim_ipm10  = ivar
        if (trim(var_name(ivar)).eq.'pPPM'  )chim_ippm   = ivar
        if (trim(var_name(ivar)).eq.'pHNO3' )chim_iphno3 = ivar
        if (trim(var_name(ivar)).eq.'temp'  )met_itemp   = ivar
        if (trim(var_name(ivar)).eq.'wins'  )met_iwins   = ivar				
      end do
      
      if (chim_ippm.ne.0) var_name_cdf(chim_ippm) = var_name(chim_ippm)
      
      if (chim_iphno3.ne.0) then
      ! Regles speciale pHNO3      
       if (nf90_inq_varid(ncFileID, 'p10HNO3', xVarID).eq.nf90_noerr) then
	var_name_cdf(chim_iphno3)='p10HNO3'
	iinput_var(chim_iphno3)=.true.	
	var_ppm_name_cdf(:)    = (/'p10PPM','p10OCAR','p10BCAR'/)
       end if		 
      ! CHIMERE standard
       if (nf90_inq_varid(ncFileID, 'pHNO3', xVarID).eq.nf90_noerr) then
        var_name_cdf(chim_iphno3)='pHNO3'
	iinput_var(chim_iphno3)=.true.
	var_ppm_name_cdf(:)    = (/'pPPM','pOCAR','pBCAR'/)
       end if 
      end if
      
      ! Regles speciale temperature   
      if (met_itemp.ne.0) then   
       if (nf90_inq_varid(ncFileID, 'temp', xVarID).eq.nf90_noerr) then
	var_name_cdf(met_itemp)='temp'
	iinput_var(met_itemp)=.true.	
       end if		  
      end if
      
      ! Regles speciale vent
      if (met_iwins.ne.0) then
        if ( (nf90_inq_varid(ncFileID, 'winm', xVarID).eq.nf90_noerr).and.&
             (nf90_inq_varid(ncFileID, 'winz', xVarID).eq.nf90_noerr) )then
	  iinput_var(met_iwins)=.true.	
        end if      		  
      end if
                                             
      ! redefini les sorties
      do ivar=1, nvar
            
        if (.not.iinput_var(ivar)) then         
          icorrect_z_var(ivar)=.false. 
          ioutput_var(ivar)   =.false.
          ioutput_dz_var(ivar)=.false.
          ioutput_de_var(ivar)=.false.
        end if   	
       
      end do            

      call check(nf90_close( ncFileID))

      end subroutine check_cdf_var     
!--------------------------------------------------------------------     
      subroutine check_cdf
      
      use netcdf
      use typesizes

      use params
      use params_cdf     
        
      implicit none 
      
    ! grille CHIMERE
      if (idebug)write(*,*) 'Verifie la presence de PM10 dans CHIMERE'    
      call check(nf90_open(fout1, nf90_nowrite, ncFileID)) 
      if (nf90_inq_varid(ncFileID, 'PM10', xVarID).eq.nf90_noerr) then
        if (idebug)write(*,*) '> AEROSOLS DISPONIBLES'
	ipm=.true.
      else
        if (idebug)write(*,*) '> AEROSOLS INDISPONIBLES'      
        ipm=.false.      
      end if  
      !vent
      if (idebug)write(*,*) 'Verifie la presence de winm,winz dans CHIMERE'       
        if ( (nf90_inq_varid(ncFileID, 'winm', xVarID).eq.nf90_noerr).and.&
             (nf90_inq_varid(ncFileID, 'winz', xVarID).eq.nf90_noerr) )then   
        if (idebug)write(*,*) '> VENT/TEMP. DISPONIBLE'
	iwin = .true.
      else     
        if (idebug)write(*,*) '> VENT/TEMP. INDISPONIBLE'  
	iwin = .false.	     
      end if         
      call check(nf90_close( ncFileID))
                   
      ! METEO
      if (idebug)write(*,*) 'Verifie la presence du fichier meteo. CHIMERE'        
      if ( nf90_open(fmet1, nf90_nowrite, ncFileID).eq.nf90_noerr) then	      
        if (idebug)write(*,*) '> METEO DISPONIBLE'
        imeteo=.true.
        call check(nf90_close( ncFileID))	
      else
        if (idebug)write(*,*) '> METEO INDISPONIBLE'      
        imeteo=.false.	 
      end if

      end subroutine check_cdf     
!--------------------------------------------------------------------                   
      subroutine read_cdf_var2d(file_cdf, var_cdf, var_nx, var_ny, dxmin, dymin, var2d)

      use netcdf
      use typesizes
      
      implicit none  
      
      ! local 
      character(len=*) :: file_cdf
      character(len=*) :: var_cdf
      integer :: var_nx, var_ny
      real    :: var2d(var_nx, var_ny)
      integer :: xVarID, ncFileID
      integer :: dxmin, dymin      
      
      integer :: start2(2), count2(2)
      
      start2 = (/dxmin , dymin  /)
      count2 = (/var_nx, var_ny /)          
            
      call check(nf90_open(file_cdf, nf90_nowrite, ncFileID))     
      call check(nf90_inq_varid(ncFileID, trim(adjustl(var_cdf)) , xVarID))    
      call check(nf90_get_var(ncFileID, xVarID, var2d,  start = start2, count = count2 ))     
      call check(nf90_close(ncFileID))
      
      end subroutine
!--------------------------------------------------------------------                   
      subroutine read_cdf_ndim(file_cdf,dim_cdf,ndim_cdf)

      use netcdf
      use typesizes
      
      implicit none  

      ! local 
      character(len=*) :: file_cdf
      character(len=*) :: dim_cdf
      integer :: ndim_cdf
      integer :: xDimID, ncFileID    
      
      !write(*,*) trim(dim_cdf)
      call check(nf90_open(file_cdf, nf90_nowrite, ncFileID))
      call check(nf90_inq_dimid(ncFileID, trim(dim_cdf)  , xDimID))  
      call check(nf90_Inquire_Dimension(ncFileID, xDimID   , len=ndim_cdf))    
      call check(nf90_close(ncFileID))
      
      !write(*,*) 'dim '//trim(adjustl(dim_cdf))//'=', ndim_cdf
      
      end subroutine
       
!--------------------------------------------------------------------                   
      subroutine read_cdf_chimere_date(file_cdf,time_cdf,Year,Mnth,Day,Hour,Min,Sec)

      use netcdf
      use typesizes
      
      implicit none  
      
      ! local 
      character(len=*) :: file_cdf
      character(len=*) :: time_cdf
      
      integer :: xVarID, ncFileID  
      
      character(len=19) :: datestr      
      
      integer :: Year		!Four-digit year of start time
      integer :: Mnth		!Month of start time
      integer :: Day		!Day of start time
      integer :: Hour		!Hour of start time
      integer :: Min		!Minute of start time
      integer :: Sec		!Second of start time	      
        
      
      !write(*,*) 'Lit date de CHIMERE'
      call check(nf90_open(file_cdf, nf90_nowrite, ncFileID))      
      call check(nf90_inq_varid(ncFileID, trim(time_cdf), xVarID))    
      call check(nf90_get_var(ncFileID, xVarID, datestr, start = (/1,1/) ))
      call check(nf90_close(ncFileID))
               	      
      read(datestr,'(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)')Year,Mnth,Day,Hour,Min,Sec
      !write(*,*) 'Trouve date initiale : '//datestr

      end subroutine
!-----------------------------------------------------------------------------      
      subroutine read_cdf_chimere_lonlat(file_cdf, nx, ny, lon, lat, idebug )

      use netcdf
      use typesizes
      
      implicit none  
      
      character(len=*) :: file_cdf      
      integer :: nx, ny
      real    :: lon(nx,ny)
      real    :: lat(nx,ny)
      integer :: dxmin, dymin 
      
      logical :: idebug     
       
      ! Grille CHIMERE
      if (idebug)write(*,*) 'Lit '//trim(file_cdf)      
      dxmin = 1
      dymin = 1                        
      call read_cdf_var2d(file_cdf ,'lon'	 ,nx ,ny, dxmin, dymin ,lon	 )
      call read_cdf_var2d(file_cdf ,'lat'	 ,nx ,ny, dxmin, dymin ,lat	 )	  	  
  
      end subroutine    
!-----------------------------------------------------------------------------      
      subroutine read_chimere_ground_cdf(file_cdf, nx, ny, lon, lat, landuse, topo, idebug )

      use netcdf
      use typesizes
      
      implicit none  
      
      character(len=*) :: file_cdf      
      integer :: nx, ny
      real    :: lon(nx,ny)
      real    :: lat(nx,ny)
      real    :: landuse(nx,ny)
      real    :: topo(nx,ny)
      integer :: dxmin, dymin 
      
      logical :: idebug     
       
      ! Grille CHIMERE

      if (idebug)write(*,*) 'Lit NetCDF '//trim(file_cdf)
      
      dxmin = 1
      dymin = 1      
                   
      call read_cdf_var2d(file_cdf ,'LON'	 ,nx ,ny, dxmin, dymin ,lon	 )
      call read_cdf_var2d(file_cdf ,'LAT'	 ,nx ,ny, dxmin, dymin ,lat	 )	  
      call read_cdf_var2d(file_cdf ,'GroundClass',nx ,ny, dxmin, dymin ,landuse)
      call read_cdf_var2d(file_cdf ,'topo'	 ,nx ,ny, dxmin, dymin ,topo	 )	  
  
      end subroutine
!---------------------------------------------------------------------------------------------------------------------    
      subroutine create_netcdf
            
      use netcdf
      use typesizes
      
      use params
      use params_cdf      
      
      implicit none
      
      !entree     
      character(len=16) :: version_str      
      
      ! local
      INTEGER :: lonDimID,latDimID,zDimID,frTimeDimID      
      INTEGER :: lonVarID,latVarID,frTimeStrDimID
      integer :: deptDimId, zoneDimId 

      ! local
      character(len=19):: timeStr !YYYY-MM-DD HH:MM:SS
      character(len=4) :: yearStr
      character(len=2) :: mnthStr
      character(len=2) :: dayStr
      character(len=2) :: hourStr
      character(len=2) :: minStr
      character(len=2) :: secStr  
      
      integer :: count4d(4)
      integer :: count3d(3)  
      integer :: count2d(2)           
          
! Grillle de  sorties
      character(len=39):: tunits          
      
      write(version_str,'(F10.2)') version

  !Define timeString  
      Min = 0
      write(yearStr,'(I4)') Year
      if (Mnth.LT.10) then
        write(mnthStr,'(I1)') Mnth
      else
        write(mnthStr,'(I2)') Mnth
      end if
 
      if (Day.LT.10) then
        write(dayStr,'(I1)') Day
      else
        write(dayStr,'(I2)') Day
      end if
 
      if (Hour.LT.10) then
        write(hourStr,'(I1)') Hour
      else
        write(hourStr,'(I2)') Hour
      end if
 
      if (Min.LT.10) then
        write(minStr,'(I1)') Min
      else
        write(minStr,'(I2)') Min
      end if
 
      secStr='00'

    tunits = 'seconds since '//yearStr//'-'//trim(mnthStr)//'-'//trim(dayStr)//' '&
      //trim(hourStr)//':'//trim(minStr)//':'//trim(secStr)
               
    call check(nf90_create(path = fout2, cmode = nf90_clobber, ncid = ncFileID))

    if (idebug)write(*,*) 'Dimensions...'
    call check(nf90_def_dim(ncid = ncFileID, name = 'west_east'  ,len = nx2, dimid = lonDimID))        
    call check(nf90_def_dim(ncid = ncFileID, name = 'south_north',len = ny2, dimid = latDimId  ))     
    call check(nf90_def_dim(ncid = ncFileID, name = 'bottom_top' ,len = 1, dimid = zDimId  ))      
    call check(nf90_def_dim(ncid = ncFileID, name = 'Time',len = nf90_unlimited, dimid = frTimeDimID)) 
    call check(nf90_def_dim(ncid = ncFileID, name = 'DateStrLen',len = 19, dimid = frTimeStrDimID))  
    call check(nf90_def_dim(ncid = ncFileID, name = 'departements',len = 99, dimid = deptDimId  ))  
    call check(nf90_def_dim(ncid = ncFileID, name = 'zones',len = 13, dimid = zoneDimId  ))                 
      
    if (idebug)write(*,*) 'Variables...'
    
    count4d = (/lonDimID,latDimId,zDimId,frTimeDimID/)
    count3d = (/lonDimID,latDimId       ,frTimeDimID/)
    count2d = (/lonDimID,latDimId                   /)       
      
    call check(nf90_def_var(ncFileID, 'Time', nf90_float, frTimeDimID, XVarID) )
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','forecast time'))
    call check(nf90_put_att(ncFileID, XVarID, 'units',tunits))     

    call check(nf90_def_var(ncFileID, 'Times', nf90_char, (/frTimeStrDimID,frTimeDimID/), XVarID) )
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','forecast time'))
  
    call check(nf90_def_var(ncFileID,'lon'    ,nf90_float,dimids=(/lonDimID,latDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Longitude centre mailles'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','degrees'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )    
       
    call check(nf90_def_var(ncFileID,'lat'    ,nf90_float,dimids=(/lonDimID,latDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Latitude centre mailles'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','degrees'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )  
    
    call check(nf90_def_var(ncFileID,'topo_grille_fin'    ,nf90_float,dimids=(/lonDimID,latDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Topographie grille FINE'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','meters'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )    
       
    call check(nf90_def_var(ncFileID,'topo_grille_coa'    ,nf90_float,dimids=(/lonDimID,latDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Topographie grille COARSE'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','meters'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )            
         	         
    do ivar = 1, nvar   
      if (ioutput_var(ivar)) then

        call check(nf90_def_var(ncFileID,trim(var_name(ivar)),nf90_float,dimids=count4d,varID=XVarID))    
        call check(nf90_put_att(ncFileID, XVarID, 'long_name',trim(var_name(ivar))))
        call check(nf90_put_att(ncFileID, XVarID, 'units',trim(var_unit(ivar))))       
        call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )   
        
	if (ioutput_krig_var(ivar)) then
	  call check(nf90_def_var(ncFileID,trim(var_name(ivar))//'_KRIG',nf90_float,dimids=count4d,varID=XVarID))    
          call check(nf90_put_att(ncFileID, XVarID, 'long_name',trim(var_name(ivar))//' assimile'))
          call check(nf90_put_att(ncFileID, XVarID, 'units',trim(var_unit(ivar))))
          call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) ) 
	end if
	
        if (idebug.and.icorrect_z_var(ivar).and.ioutput_dz_var(ivar)) then
          write(*,*)' DEBUG: creer dz_'//trim(var_name(ivar))
          call check(nf90_def_var(ncFileID,'dz_'//trim(var_name(ivar))    ,nf90_float,dimids=count4d,varID=XVarID))    
          call check(nf90_put_att(ncFileID, XVarID, 'long_name','Correction altitude '//trim(var_name(ivar))))
          call check(nf90_put_att(ncFileID, XVarID, 'units',trim(var_unit(ivar))))
          call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) ) 
       end if
       
        if (idebug.and.icorrect_z_var(ivar).and.ioutput_de_var(ivar)) then       
          write(*,*)' DEBUG: creer de_'//trim(var_name(ivar))
          call check(nf90_def_var(ncFileID,'de_'//trim(var_name(ivar))    ,nf90_float,dimids=count4d,varID=XVarID))    
          call check(nf90_put_att(ncFileID, XVarID, 'long_name','Correction emissions '//trim(var_name(ivar))))
          call check(nf90_put_att(ncFileID, XVarID, 'units',trim(var_unit(ivar))))
          call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) ) 	  
        end if 	
	             
       end if
    end do
    
    !Ecrit les emissions ANNUELLES utilises dans GIRAF (02/2011)
    do ivar = 1, nvar_emis 
     
      if (ivar_emis(ivar).and.idebug) then
        write(*,*)' DEBUG: creer emis_'//trim(var_name(ivar))
        call check(nf90_def_var(ncFileID,'emis_'//trim(var_emis_name(ivar)),nf90_float,dimids=count4d,varID=XVarID))    
        call check(nf90_put_att(ncFileID, XVarID, 'long_name',trim(var_emis_name(ivar))))
        call check(nf90_put_att(ncFileID, XVarID, 'units','kg/an'))       
        call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )         		             
       end if
       
      if (ivar_emis(ivar).and.idebug) then
        write(*,*)' DEBUG: creer ratio_emis_'//trim(var_name(ivar))
        call check(nf90_def_var(ncFileID,'ratio_emis_'//trim(var_emis_name(ivar)),nf90_float,dimids=count4d,varID=XVarID))    
        call check(nf90_put_att(ncFileID, XVarID, 'long_name',trim(var_emis_name(ivar))))
        call check(nf90_put_att(ncFileID, XVarID, 'units','kg/an'))       
        call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )         		             
       end if       
       
    end do
    
    
    if (iwin) then    
      call check(nf90_def_var(ncFileID,'fac_w10m',nf90_float,dimids=count3d,varID=XVarID))    
      call check(nf90_put_att(ncFileID, XVarID, 'long_name','Vitesse du vent a 10 metres'))
      call check(nf90_put_att(ncFileID, XVarID, 'units','m/s'))
      call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )        
    end if      
    
    if (ipm) then
      call check(nf90_def_var(ncFileID,'PM10_nonvolat',nf90_float,dimids=count4d,varID=XVarID))    
      call check(nf90_put_att(ncFileID, XVarID, 'long_name','PM10 TEOM 50 degC'))
      call check(nf90_put_att(ncFileID, XVarID, 'units','microg/m3'))
      call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )    
    
      !call check(nf90_def_var(ncFileID,'PM10_ecart_fdms',nf90_float,dimids=count4d,varID=XVarID))    
      !call check(nf90_put_att(ncFileID, XVarID, 'long_name','ecart FDMS (Militon, 2009)'))
      !call check(nf90_put_att(ncFileID, XVarID, 'units','microg/m3'))
      !call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )
    end if          

!    call check(nf90_def_var(ncFileID,'tem2'    ,nf90_float,dimids=count4d,varID=XVarID))    
!    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Temperature a 2 metres'))
!    call check(nf90_put_att(ncFileID, XVarID, 'units','K'))
!    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) ) 
      
!    call check(nf90_def_var(ncFileID,'sreh'    ,nf90_float,dimids=count4d,varID=XVarID))    
!    call check(nf90_put_att(ncFileID, XVarID, 'long_name','HR a 2 metres'))
!    call check(nf90_put_att(ncFileID, XVarID, 'units','%'))
!    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )     


    if (region1km) then
    call check(nf90_def_var(ncFileID,'contour_region'    ,nf90_float,dimids=(/lonDimID,latDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Contour de la regoin Rhone-Alpes'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','0/1'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) ) 
           
    call check(nf90_def_var(ncFileID,'dept'    ,nf90_float,dimids=(/lonDimID,latDimId,deptDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Contour departements francais'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','0/1'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )                   

    call check(nf90_def_var(ncFileID,'zone'    ,nf90_float,dimids=(/lonDimID,latDimId,zoneDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Contour departements francais'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','0/1'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )   

    call check(nf90_def_var(ncFileID,'pop99'   ,nf90_float,dimids=(/lonDimID,latDimId/),varID=XVarID))    
    call check(nf90_put_att(ncFileID, XVarID, 'long_name','Population au km carre'))
    call check(nf90_put_att(ncFileID, XVarID, 'units','hab/km2'))
    call check(nf90_put_att(ncFileID, XVarID,  "_FillValue", -9999.0 ) )    
        
    end if
    
    ! Global attributes
    write(yearStr,'(I4)')annee
    call check(nf90_put_att(ncFileID, nf90_global, 'title'  , 'Raffinement de maillage v'//trim(adjustl(version_str))//' - (c)Air-RA '//yearStr ))      
    call check(nf90_put_att(ncFileID, nf90_global, 'version', trim(adjustl(version_str))))   
    
    call check(nf90_close(ncFileID))     

      end subroutine
!---------------------------------------------------------------------------------------------------------------------   
      subroutine open_cdf
      
      use netcdf
      use typesizes

      use params
      use params_cdf      
        
      implicit none 
      
      call check(nf90_open(fout1, nf90_nowrite, out1FileID))
      if (imeteo)call check(nf90_open(fmet1, nf90_nowrite, met1FileID))

      end subroutine open_cdf      
!--------------------------------------------------------------------
      subroutine read_cdf_1h
      
      use netcdf
      use typesizes

      use params
      use params_cdf    
  
      implicit none 
      
      real    :: val1_tmp_3D(nx0,ny0,nz0)
      real    :: val2_tmp_3D(nx0,ny0,nz0)
            
      real    :: val_tmp_2D(nx0,ny0    )      
      integer :: ivar_ppm
            
      integer :: start4(4)
      integer :: start3(3)
      integer :: start2(2)
            
      !if (it.eq.1) write(*,*) 'Lit la date'
      ! date             	      
      call check(nf90_inq_varid(out1FileID, 'Times', xVarID))    
      call check(nf90_get_var(out1FileID, xVarID, datestr, start = (/1,it/) ))
     	      
      read(datestr,'(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)')&
     	     Year,Mnth,Day,Hour,Min,Sec
      write(*,*) '*******************************************************************'  
      if (idebug)write(*,*) '*******************************************************************' 
      if (idebug)write(*,*) '***                                                             ***'       
      write(*,'(" *** Iteration ",I4,"/",I4)') it,numFrTimes 
      write(*,*) '*** Traite date : '//datestr
      if (idebug)write(*,*) '***                                                             ***' 
      if (idebug)write(*,*) '*******************************************************************'  
      if (idebug)write(*,*) '*******************************************************************'                
      ! out.*.nc
      start4 = (/dx1min,dy1min,   1, it/)      
      start3 = (/dx1min,dy1min,      it/)  
                  
      if (idebug) write(*,*) 'Lit les variables chimie dans '//trim(fout1)      
      
      ippm = .false.
      
      do ivar = 1, nvar
      
        ! SUPPORT DES PPM DANS CHIMERE 2008
        if ((ivar.eq.chim_ippm).and.ipm) then 
	  if (idebug) write(*,*) '... '//trim(var_name_cdf(ivar)) 
	  val_g0_var(ivar,:,:,:) = 0.
	  do ivar_ppm = 1, nvar_ppm
	    if (nf90_inq_varid(out1fileID, trim(var_ppm_name_cdf(ivar_ppm)), xVarID).eq.nf90_noerr) then
	      if (idebug) write(*,*) '...... '//trim(var_ppm_name_cdf(ivar_ppm)) 
	      ippm = .true.            
	      call check(nf90_get_var(out1fileID, xVarID, val1_tmp_3D,  start = start4 ))
	      val_g0_var(ivar,:,:,:) = val_g0_var(ivar,:,:,:) + val1_tmp_3D
	      if (idebug) write(*,*) 'valmax=',maxval(val1_tmp_3D)
	    end if
	  end do
	  if (.not.ippm) write(*,*) 'ECHEC: ne trouve aucune espece PPM'	
	  
        ! SUPPORT VITESSE DU VENT
        else if (ivar.eq.met_iwins) then 
	    if (iwin) then
	      if (idebug)write(*,*) 'info : lit winm et winz'	      
	      if (nf90_inq_varid(out1fileID, 'winz' , xVarID).eq.nf90_noerr) &
	         call check(nf90_get_var(out1fileID, xVarID, val1_tmp_3D,  start = start4 ))
	      if (nf90_inq_varid(out1fileID, 'winm' , xVarID).eq.nf90_noerr) &
	        call check(nf90_get_var(out1fileID, xVarID, val2_tmp_3D,  start = start4 ))	    	      
	      if (idebug)write(*,*) 'info : calcul wins'
	      call calcul_vitesse_vent_3d(nx0,ny0,nz0,val1_tmp_3D,val2_tmp_3D,val_g0_var(ivar,:,:,:))
	      if (idebug)write(*,*) 'info : Lecture vent -> OK'
	    end if 	    
	! VARIABLES AUTRES QUE PPM	
	else  
          if (idebug)write(*,*) '... '//trim(var_name_cdf(ivar))           
	  if ( nf90_inq_varid(out1fileID, trim(var_name_cdf(ivar)), xVarID) .eq.nf90_noerr) then
	    if (idebug)write(*,*) 'info : lit '//trim(var_name_cdf(ivar)) 	    
	    call check(nf90_get_var(out1fileID, xVarID, val_g0_var(ivar,:,:,:),  start = start4 ))
	  else
	    if (idebug)write(*,*) 'info : ignore '//trim(var_name_cdf(ivar)) 
	  end if   
	     	
	end if
      end do
      
      ! SUPPORT VITESSE DU VENT 10 M
      if (iwin) then
	if (idebug)write(*,*) 'info : lit w10m'	    
	if (nf90_inq_varid(out1fileID, 'w10m' , xVarID).eq.nf90_noerr) then
	  call check(nf90_get_var(out1fileID, xVarID, w10m_g0,  start = start3 ))
	else
	  if (idebug)write(*,*) '-> echoue'	    
	end if
      end if		 	         

      ! met.*.nc ! PAS UTILISE     
      if (imeteo) then
        if (it.eq.1) write(*,*) 'Lit les variables meteo'        
        call check(nf90_inq_varid(met1fileID, 'alti', xVarID))         
        call check(nf90_get_var(met1fileID, xVarID, alti_agl_g0,start = start4 ))

        call check(nf90_inq_varid(met1fileID, 'hght', xVarID))         
        call check(nf90_get_var(met1fileID, xVarID, hght_g0,    start = start4 ))

        call check(nf90_inq_varid(met1fileID, 'tem2', xVarID))         
        call check(nf90_get_var(met1fileID, xVarID, tem2_g0,    start = start3 ))

      else ! version 2008
        ! lit l'altitude des mailles dans le fichier NetCDF CHIMERE: variabe hlay ou alti
        if(idebug) write(*,*) 'Lit variable meteo hlay/alti' 	
        if (nf90_inq_varid(out1fileID, 'hlay', xVarID).eq.nf90_noerr) &
          call check(nf90_get_var(out1fileID, xVarID, alti_agl_g0,start = start4 ))
        if (nf90_inq_varid(out1fileID, 'alti', xVarID).eq.nf90_noerr) &
          call check(nf90_get_var(out1fileID, xVarID, alti_agl_g0,start = start4 ))	  
      end if                 

      end subroutine read_cdf_1h
!--------------------------------------------------------------------
      subroutine write_emis_cdf
      
      use netcdf
      use typesizes

      use params
      use params_cdf    
 
      implicit none
      
      integer :: start4(4)

      start4 = (/   1,   1,   1, it/)     
                
      call check(nf90_open(fout2, nf90_write, ncFileID)) 
            
      ! Emissions annuelles sur grille fine
      do ivar = 1, nvar_emis           	
	if (nf90_inq_varid(ncFileID, 'emis_'//trim(var_emis_name(ivar)), xVarID).eq.nf90_noerr) then
           if (idebug)write(*,*) 'emis_'//trim(var_emis_name(ivar))		    
           call check(nf90_put_var(ncFileID, XVarID, emis_g2(ivar,:,:), start = start4 ))		     
	end if  
	
	if (nf90_inq_varid(ncFileID, 'ratio_emis_'//trim(var_emis_name(ivar)), xVarID).eq.nf90_noerr) then
           if (idebug)write(*,*) 'ratio_emis_'//trim(var_emis_name(ivar))		    
           call check(nf90_put_var(ncFileID, XVarID, ratio_emis_g2(ivar,:,:), start = start4 ))		     
	end if  	
	     	
      end do
           	      
      call check(nf90_close(ncFileID))  
            
      end subroutine write_emis_cdf
!--------------------------------------------------------------------
      subroutine write_cdf_1h
      
      use netcdf
      use typesizes

      use params
      use params_cdf    
  
      implicit none            
                
      call check(nf90_open(fout2, nf90_write, ncFileID)) 
      
      if (it.eq.1) then
        !write(*,*) 'Longitude Latitude'
	call check(nf90_inq_varid(ncFileID, 'lon', xVarID))    	
        call check(nf90_put_var(ncFileID, xVarID, lon_g2, start = (/ 1, 1/) )) 
	call check(nf90_inq_varid(ncFileID, 'lat', xVarID))	
        call check(nf90_put_var(ncFileID, xVarID, lat_g2, start = (/ 1, 1/) )) 
	call check(nf90_inq_varid(ncFileID, 'topo_grille_fin', xVarID))	
        call check(nf90_put_var(ncFileID, xVarID, topo_g2, start = (/ 1, 1/) )) 
	call check(nf90_inq_varid(ncFileID, 'topo_grille_coa', xVarID))	
        call check(nf90_put_var(ncFileID, xVarID, topo_g1, start = (/ 1, 1/) )) 		
	
	if ( region1km )call check(nf90_inq_varid(ncFileID, 'contour_region', xVarID))	
        if ( region1km )call check(nf90_put_var(ncFileID, xVarID, rregion, start = (/ 1, 1/) )) 
	if ( region1km )call check(nf90_inq_varid(ncFileID, 'dept', xVarID))	
        if ( region1km )call check(nf90_put_var(ncFileID, xVarID, idept_region, start = (/ 1, 1, 1/) )) 
	if ( region1km )call check(nf90_inq_varid(ncFileID, 'zone', xVarID))	
        if ( region1km )call check(nf90_put_var(ncFileID, xVarID, izone_region, start = (/ 1, 1, 1/) )) 
	if ( region1km )call check(nf90_inq_varid(ncFileID, 'pop99', xVarID))	
        if ( region1km )call check(nf90_put_var(ncFileID, xVarID, pop_region, start = (/ 1, 1/) )) 					
      end if

      if (idebug)write(*,*) 'Time'      
      call check(nf90_inq_varid(ncFileID, 'Time', xVarID))
      call check(nf90_put_var(ncFileID, XVarID, (it-1)*3600., start = (/it/) )) 

      if (idebug)write(*,*) 'Times'      
      call check(nf90_inq_varid(ncFileID, 'Times', xVarID))
      call check(nf90_put_var(ncFileID, XVarID, datestr, start = (/1,it/) )) 

      ! Profils verticaux
      do ivar = 1, nvar
       if ( nf90_inq_varid(ncFileID, 'dz_'//trim(var_name(ivar)), xVarID) .eq. nf90_noerr ) then          
         if (idebug)write(*,*) 'DEBUG : ecrit dz_'//trim(var_name(ivar))
         call check(nf90_put_var(ncFileID, XVarID, dvar_dz_loc(ivar,:,:), start = (/ 1,1,1,it/) ))            
       end if 
       
       if ( nf90_inq_varid(ncFileID, 'de_'//trim(var_name(ivar)), xVarID) .eq. nf90_noerr ) then          
         if (idebug)write(*,*) 'DEBUG : ecrit de_'//trim(var_name(ivar))
         call check(nf90_put_var(ncFileID, XVarID, dvar_de_loc(ivar,:,:), start = (/ 1,1,1,it/) ))            
       end if                    
      end do         
      
      ! Variables sorties
      do ivar = 1, nvar      	
	if (nf90_inq_varid(ncFileID, trim(var_name(ivar)), xVarID).eq.nf90_noerr) then
           if (idebug)write(*,*) trim(var_name(ivar))		    
           call check(nf90_put_var(ncFileID, XVarID, val_g2_var(ivar,:,:,:), start = (/ 1,1,1,it/) ))		     
	end if       	
      end do
      
      ! facteur correctif fin/coarse pour la vitesse 10 m (version 1.11)
      if ( iwin ) then
	if (nf90_inq_varid(ncFileID, 'fac_w10m', xVarID).eq.nf90_noerr) then
          if (idebug)write(*,*) 'fac_w10m'
          call check(nf90_put_var(ncFileID, XVarID, w10m_g2, start = (/ 1,1,it/) ))	     
	end if
      end if          
      
      ! FDMS  
      if (nf90_inq_varid(ncFileID, 'PM10_ecart_fdms', xVarID).eq.nf90_noerr) then
         if (idebug)write(*,*) 'PM10_ecart_fdms'         
        call check(nf90_put_var(ncFileID, XVarID, ecart_fdms_g2(:,:,:), start = (/ 1,1,1,it/) )) 
      end if
      
      if (nf90_inq_varid(ncFileID, 'PM10_nonvolat', xVarID).eq.nf90_noerr) then
          if (idebug)write(*,*) 'PM10_nonvolat'         
        call check(nf90_put_var(ncFileID, XVarID, pm10_nonvolat_g2(:,:,:), start = (/ 1,1,1,it/) )) 
      end if      
	      
      call check(nf90_close(ncFileID))  
      
      
      end subroutine write_cdf_1h
!--------------------------------------------------------------------
