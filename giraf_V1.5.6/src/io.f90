      subroutine header
      
      use params
      
      implicit none
      
      character(len=16)  :: version_str
      character(len=4)   :: annee_str        

      write(version_str,'(F10.2)') version
      write(annee_str,'(I4)') annee      

      write(*,*) '        ##                              #       #'  
      write(*,*) '      #####         #                    #  ####'
      write(*,*) '     ########     ###    #                #     #'
      write(*,*) '##  ########## #######  ###              #          #'
      write(*,*) '##########################################         #' 
      write(*,*) '# GIRAF - Raffinement maillage chimie/meteo #  ####'
      write(*,*) '# (c)Air Rhone-Alpes '//trim(annee_str)//'                ##      #'
      write(*,*) '# version '//trim(version_str)
      write(*,*) '# contact: echaxel@air-rhonealpes.fr'                                                                                
      write(*,*) '###########################################################'
      write(*,*) '# # ## ### #### ##### ####### ######## ######### ##########'
      
      end subroutine

!--------------------------------------------------------------------
      subroutine help_me
      
      use params
      
      implicit none      
      write(*,*) 'POUR AVOIR DE L AIDE: '//trim(exe_str)//' -h'
      write(*,*) ''
      write(*,*) 'BUT     : A partir d un calcul CHIMERE a resolution n km obtenir des resultats a resolution 1 km' 
      write(*,*) ''      
      write(*,*) 'PRINCIPE: A partir des emissions et de la topographie sur une maille de 1 km, le programme'
      write(*,*) '          calcul les concentrations en NO, NO2, O3 et PM10 sur une grille fine a 1 km'
      write(*,*) '          dans la suite COARSE=domaine CHIMERE a dx=n km et FIN=domaine rafine a dx=1 km'
      write(*,*) ''      
      write(*,*) 'SYNTAXE : '//trim(exe_str)//' -i out_chimere_COARSE.nc -t2 topo_FIN -e2 emissions_FIN [ -m meteo_chimere_COARSE.nc -t1 topo_COARSE -e1 emissions_COARSE -o out_chimere_FIN.nc]'  
      write(*,*) ''	   
      write(*,*) 'ENTREE  : -i out_chimere_COARSE.nc = fichier de sortie de CHIMERE au format NetCDF'
      write(*,*) '          -t2 topo_FIN             = fichier contenant [lon, lat, topo] du domaine FIN au format ASCII convention CHIMERE ou NetCDF'
      write(*,*) '          -e2 emissions_FIN        = emissions annuelles (surfaciques au sol) au format ASCII [lon, lat, NOx, PMs] convention CHIMERE ou NetCDF' 
      write(*,*) '                                     du domaine FIN au format ASCII ou NetCDF (unites au choix)'
      write(*,*) ''            
      write(*,*) 'OPTIONS : '  
      write(*,*) '          -r  mailles_region.txt      = fichier avec des donnes sur la region x_utm y_utm population_km2 zone	dept'        
      write(*,*) '          -m  meteo_chimere_COARSE.nc = utilise le fichier meteo de CHIMERE() pour calculer la temperature 2 m reelle'
      write(*,*) '          -t1 topo_COARSE             = fichier contenant [lon, lat, topo] du domaine COARSE (CHIMERE)'        
      write(*,*) '          -e1 emissions_COARSE        = emissions annuelles au format [lon, lat, NOx, COVNM, PPM] du domaine COARSE (CHIMERE)'
      write(*,*) '          -o  out_chimere_FIN.nc      = speficie un nom de fichier de resultat du rafinement (par defaut out.giraf.nc)'
      write(*,*) '          -idebug                     = mode de debogage avec beaucoup affichage'      
      write(*,*) ''
      write(*,*) 'REMARQUE : en cas d absence des fichiers -t2 ou -e2, le rafinement de maillage pour l altitude ou les emissions seront desactivees'       
      write(*,*) ''
      write(*,*) 'SORTIES : out_chimere_FIN.nc'
      write(*,*) '          - format NetCDF de type CHIMERE'
      write(*,*) '          - concentrations en NO, NO2, O3, PM10, PM25, pHNO3, PM10_nonvolat et PM10_ecart_fdms'
      write(*,*) '          - variables meteo : tem2, [u10 et v10 (si fichier meteo.nc disponible) dans version 1.1]'      
      write(*,*) ''
  
                   
      stop    
      
      
      end subroutine help_me
!--------------------------------------------------------------------
      subroutine read_args
      
      use params
      
      implicit none  
      
      character(len=fnstrl) :: argstr
      integer :: iarg
      logical :: ihelp

      ! logical du programme
      idebug    =.false. ! debogage
      iraf_alti =.true.  ! rafinement en altitude
      iraf_chim =.true.  ! rafinement en emissions

      ! logical local
      iout1=.false.
      iout2=.false.      
      iground1=.false.
      iground2=.false.
      imet=.false.
      iemis1=.false.
      iemis2=.false.
      ihelp=.false.
      ifregion=.false.
      
      proj='geo' !defaut: projection geographique
      
      ! lit l'argument
      iarg = 0    
      argstr='null'	  
      do while ( trim(adjustl(argstr)) .ne. '' )
	iarg = iarg + 1	
        !write(*,*) 'Lit argument', iarg	
        call getarg(iarg,argstr)	
        if ( trim(adjustl(argstr)).eq. '-i' )  then
	  call getarg(iarg+1,fout1)
	  iout1 = .true.
        else if ( trim(adjustl(argstr)).eq. '-o' ) then
	  call getarg(iarg+1,fout2)
          iout2=.true.			
        else if ( trim(adjustl(argstr)).eq. '-t1' ) then
	  call getarg(iarg+1,fground1)
          iground1=.true.			
        else if ( trim(adjustl(argstr)).eq. '-t2' ) then
	  call getarg(iarg+1,fground2)	
          iground2=.true.		
        else if ( trim(adjustl(argstr)).eq. '-e1' ) then
	  call getarg(iarg+1,femis1(1))	
          iemis1=.true.		
        else if ( trim(adjustl(argstr)).eq. '-e2' ) then
	  call getarg(iarg+1,femis2(1))
          iemis2=.true.		
        else if ( trim(adjustl(argstr)).eq. '-m' )  then
	  call getarg(iarg+1,fmet1)
          imet=.true.	
        else if ( trim(adjustl(argstr)).eq. '-r' )  then
	  call getarg(iarg+1,fregion)
          ifregion=.true.
        else if ( trim(adjustl(argstr)).eq. '-proj' )  then
	  call getarg(iarg+1,proj)	  		  	
        else if ( trim(adjustl(argstr)).eq. '-debug'.or.trim(adjustl(argstr)).eq. '-d' )  then
          idebug=.true.		
        else if ( trim(adjustl(argstr)).eq. '-h'.or.trim(adjustl(argstr)).eq. '--help'.or.trim(adjustl(argstr)).eq. '-help' )  then
          ihelp=.true.		
	endif	
		
      end do
      
      if (ihelp) then
        call help_me
        stop
      end if
      
      if ( iarg.eq.1) then
        write(*,*) 'warning : Ne trouve pas d argument'
	write(*,*) 'POUR AVOIR DE L AIDE: '//trim(exe_str)//' -h'
	stop
      end if
            
      ! noms par defaut
      if (.not.iout1   )fout1    = 'out1.nc' !defaut      
      if (.not.iout2   )fout2    = 'out.giraf.nc' !defaut
      if (.not.iground1)fground1 = 'topo1.nc' !defaut      
      if (.not.iground2)fground2 = 'topo2.nc' !defaut  
      if (.not.iemis1  )femis1(1)= 'emissions_an_1.nc' !defaut
      if (.not.iemis2  )femis2(1)= 'emissions_an_2.nc' !defaut         
      if (.not.imet    )fmet1    = 'meteo.nc' !defaut
      if (.not.ifregion)fregion  = 'MAILLES_CADASTRE_CADASTRE' !defaut      
      	              
      ! topo au format netCDF ?      
      inc_topo=.false.  
      if ( index(fground2,'.nc')+index(fground2,'.cdf').ne. 0 ) then 
        inc_topo=.true.
      end if
      
      ! emis au format netCDF ?
      inc_emis=.false.
      if ( index(femis2(1),'.nc')+index(femis2(1),'.cdf').ne. 0 ) then 
        inc_emis=.true.
      end if          
                  
      ! DEPOSITION SECHE (obsolete)
      fdep1    ='dep.nc'     ! inutile
      
      write(*,*) '***** FICHIERS DISPONIBLES **************************************************************************'
      
      if (iout1   )write(*,*) 'out.*.nc CHIMERE (fout1)        :'//trim(adjustl(fout1))
      if (iground1)write(*,*) 'sol CHIMERE (fground1)          :'//trim(adjustl(fground1))
      if (iemis1  )write(*,*) 'emissions an.CHIMERE (femis1)   :'//trim(adjustl(femis1(1)))

      if (iout2   )write(*,*) 'sortie domaine FIN (fout2)      :'//trim(adjustl(fout2))
      if (iground2)write(*,*) 'sol domaine FIN (fground2)      :'//trim(adjustl(fground2))
      if (iemis2  )write(*,*) 'emissions an.domaine FIN(femis2):'//trim(adjustl(femis2(1)))
      
                   write(*,*) 'projection proj='//trim(adjustl(proj))
          
      if (ifregion)write(*,*) 'donnees population Rhone-Alpes  :'//trim(adjustl(fregion)) 
      
      if (imet    )write(*,*) 'meteo horaire CHIMERE (fmet1)   :'//trim(adjustl(fmet1))    
      
      if (imet    )write(*,*) 'meteo horaire CHIMERE (fmet1)   :'//trim(adjustl(fmet1))                 
      
      write(*,*) '*****************************************************************************************************'

      !write(*,*) femis1(2)='emi1_2.nc'  
      !write(*,*) femis1(3)='emi1_3.nc'  
      !write(*,*) femis2(2)='emi2_2.nc'  
      !write(*,*) femis2(3)='emi2_3.nc'  
      !write(*,*) fdep1	 ='dep.nc'     
                  
      end subroutine read_args     
