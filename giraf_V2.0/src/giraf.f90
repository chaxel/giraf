      program raffinement_maillage
      
!********************************************************************
!*         PREVALP                               		    *
!*     Raffinement de maillage CHIMERE/PREVAIR                      *
!*          auteur: E. Chaxel (Air Rhone-Alpes)                     *
!*          email:  echaxel@air-rhonealpes.fr                       *
!********************************************************************
      use netcdf    ! parametres NetCDF
      use typesizes ! parametres NetCDF
      use params    ! parametres generaux du programmes
      
      implicit none  
      
      ! local
      integer :: ivar,ific,i,j,k                           
!--------------------------------------------------------------                  


!--------------------------------------------------------------------
      call header

! Lit les entrees
      write(*,*) '>>>>>>>>>>>> Lit les arguments...'
      call read_args
!--------------------------------------------------------------------
!      INITIALISATION
!--------------------------------------------------------------------  
   write(*,*) '>>>>>>>>>>>> Initialise le calcul...'
             
! Fichier out.*.nc COARSE
   if (idebug)write(*,*) 'Lit dimensions de la grille CHIMERE dans '//trim(fout1)
   call read_cdf_ndim(fout1,'west_east'  ,nxc	    )
   call read_cdf_ndim(fout1,'south_north',nyc	    )
   call read_cdf_ndim(fout1,'bottom_top' ,nzc	    )	   
   call read_cdf_ndim(fout1,'Time'	 ,numFrTimes)
   if (idebug)write(*,*) 'Lit la date dans '//trim(fout1)      
   call read_cdf_chimere_date(fout1,'Times',Year,Mnth,Day,Hour,Min,Sec)
   TimeStringLen=60.	      
   if (idebug)write(*,'("Dimensions COARSE nx=",I3," ny=",I3," nz=",I3)') nxc, nyc, nzc     
   if (idebug)write(*,*) 'Lit lon,lat,topo,landuse dans '//trim(fground2)

! Fichier topo.nc FIN (ajoute lecture de ASCII)
   if (.not.iground2) then  
    write(*,*) '*** erreur : vous devez fournir un fichier avec les donnees longitude latitude altitude'
    write(*,*) '*** erreur : consulter l aide avec la commande: '//trim(exe_str)//' -h'
    write(*,*) '*** STOP' 
    stop
   end if
   if (idebug)write(*,*) 'Lit dimensions de la grille FINE dans '//trim(fground2)
   if (inc_topo) then 
      call read_cdf_ndim(fground2,'X',nx2)
      call read_cdf_ndim(fground2,'Y',ny2)
      nz2 = 1
   else
      call read_ascii_ndim(fground2,nx2,ny2)
      nz2 = 1   
   end if       
   if (idebug)write(*,'("Dimensions FINE nx=",I3," ny=",I3," nz=",I3)') nx2, ny2, nz2 

!-------------------------------------------------------------------- 
! Allocations des grilles
!--------------------------------------------------------------------      
   if (idebug)write(*,*) 'Allocations des grilles CHIMERE brutes et FINE'   
   call allocation_grilles
!-------------------------------------------------------------------- 
! Lit coordonnees de CHIMERE
!--------------------------------------------------------------------      
   if (idebug)write(*,*) 'Lit lon,lat CHIMERE dans '//trim(fout1)	
   call read_cdf_chimere_lonlat(fout1, nxc, nyc, lon_gc, lat_gc, idebug )        
!--------------------------------------------------------------------                                   
! Grille FINE (ajoute lecture de ASCII)
   if (idebug)write(*,*) 'Lit lon,lat,topo,landuse FINE dans '//trim(fground2)		
   if (inc_topo)      call read_chimere_ground_cdf(  fground2, nx2, ny2, lon_g2, lat_g2, landuse_g2, topo_g2, idebug )  
   if (.not.inc_topo) call read_chimere_ground_ascii(fground2, nx2, ny2, lon_g2, lat_g2, landuse_g2, topo_g2, idebug )  	  

   !projection vers long/lat au cas ou...
   call proj_xy(proj,nx2,ny2,lon_g2,lat_g2)
   
   write(*,*)'Lon=',lon_g2(1,1),lon_g2(nx2, ny2)
   write(*,*)'Lat=',lat_g2(1,1),lat_g2(nx2, ny2)
!-------------------------------------------------------------------- 
! Grille COARSE (est-ce necessaire ? NON -> remplacé par la moyenne )
!   if (idebug)write(*,*) 'Lit lon,lat,topo,landuse COARSE dans '//trim(fground1)
!   call read_chimere_ground_cdf(fground1, nx0, ny0, dx1min, dy1min, lon_g0, lat_g0, landuse_g0, topo_g0, idebug )
!   if (idebug)write(*,'("Dimensions COARSE nx=",I3," ny=",I3)') nx0, ny0 
!---------------------------------
! Calcule les bornes du domaine reduit nx1, ny1 -> dx1min, dx1max, dy1min, dy1max     
   if (idebug)write(*,*) 'Calcul du domaine reduit'
   if (idebug)write(*,'("Dimensions CHIMERE nxc=",I3," nyc=",I3," nzc=",I3)') nxc, nyc, nzc    
   if (idebug)write(*,*) '-> ATTENTION les grilles des domaines doivent etre les memes si emissions !'
   call calcul_domaine_reduit( &
   nxc,nyc,nzc,lon_gc,lat_gc,&
   nx2,ny2,nz2,lon_g2,lat_g2,&
   nx0,ny0,nz0,dx1min,dx1max,dy1min,dy1max,idebug)
   if (idebug)write(*,'("Dimensions COARSE reduit nx0=",I3," ny0=",I3," nz0=",I3)') nx0, ny0, nz0   
! Domaine COARSE raffine
   nx1=nx2
   ny1=ny2
   nz1=nz0                 
 ! Sort les resultats sur une grille regionale
   if (idebug)write(*,'("Dimensions COARSE raffine nx1=",I3," ny1=",I3," nz1=",I3)') nx1, ny1, nz1
!-------------------------------------------------------------------- 
! Donnees population region Rhone-Alpes [x_utm,y_utm,pop]
!-------------------------------------------------------------------- 
   if (idebug)write(*,*) 'Lit donnees regionales...'   
   call read_mailles_region_ra_1km
!-------------------------------------------------------------------- 
! Allocations pour la suite
!--------------------------------------------------------------------      
   if (idebug)write(*,*) 'Allocations des variables 1D, 2D et 3D...'
   call allocation_var_3d      
!-------------------------------------------------------------------- 
! Defini les coordonnees du domaine COARSE
!--------------------------------------------------------------------                   
  if (idebug)write(*,*) 'Fixe les coordonnees du domaine COARSE...'
  lon_g0 = lon_gc(dx1min:dx1max,dy1min:dy1max)
  lat_g0 = lat_gc(dx1min:dx1max,dy1min:dy1max)    
  distGrid_dg = (lon_g0(nx0,int(ny0/2)) - lon_g0( 1,int(ny0/2))) / ( nx0 - 1 ) 
  distGrid_km = rayon_km * distGrid_dg * d2r
  if (idebug)write(*,*) 'Maille COARSE(deg)=',distGrid_dg     
  if (idebug)write(*,*) 'Maille COARSE(km) =',distGrid_km
!--------------------------------------------------------------------                   
! Emissions grille fine (OBLIGATOIRE sinon iraf_chim=.false.)
!--------------------------------------------------------------------    
   !Indices des emissions NOx et PM10
   emis_inox  = 0
   emis_ipm10 = 0
   do ivar=1, nvar_emis
      !if (trim(var_emis_name(ivar)).eq.'VOC' )emis_icovnm    = ivar
      if (trim(var_emis_name(ivar)).eq.'NO'  )emis_ino	 = ivar
      if (trim(var_emis_name(ivar)).eq.'NO2' )emis_ino2  = ivar   
      if (trim(var_emis_name(ivar)).eq.'NOX' )emis_inox  = ivar
      if (trim(var_emis_name(ivar)).eq.'PM10')emis_ipm10 = ivar
      if (trim(var_emis_name(ivar)).eq.'PM25')emis_ipm25 = ivar   
   end do   
   if (emis_inox*emis_ipm10.eq.0) then
      write(*,*) '***erreur: ne trouve pas NOX et PM10 dans la liste des emissions'
      stop 1
   end if

   if(iemis2) then 
    
    emis_g2 = 0.
    
    !Lit les emissions FINE
    if (inc_emis) then
      if (idebug)write(*,*) 'Lit les emissions NetCDF du domaine FINE...'    
      do ific=1, nfile_emis
        it = 1
        call read_chimere_emis_cdf( femis2(ific), it, nx2, ny2, nsnap1, 1, 1, &
	                            nvar_emis, var_emis_name, ivar_emis, tbl_ratio_emi_snap, emis_g2(:,:,:), &
				    iemis_1h, idebug  )
      end do
      !Calcul les ratios NOx
      call calcul_emissions_nox_pm      
      
    else
      ! en ASCII, plus simple un seul fichier     
      if (idebug)write(*,*) 'Lit les emissions ASCII du domaine FINE...'
      if (idebug)write(*,*) 'Lit fichier 4 colonnes X, Y, NOx, PM10'
      emis_g2 = 0.
      !Lit les emissions dans un fichier sans remplacer lon_g2(:,:) et lat_g2(:,:)
      call read_x_y_ncol_ascii( femis2(1), 2, nx2, ny2, lon_g2(:,:),lat_g2(:,:),emis_g2(1:2,:,:), idebug )               
      ivar_emis(emis_inox)  = .true.
      ivar_emis(emis_ipm10) = .true.
      ivar_emis(emis_ino)   = .true.
      ivar_emis(emis_ino2)  = .true.    
    end if
    
    !Calcul les ratios NOx
    call calcul_emissions_nox_pm
    
   else
      iraf_chim = .false.  
      write(*,'("*** info : emissions non disponibles... pas de raffinement avec emissions")')
      write(*,'("*** info : pour activer cette fonction preparer un fichier avec les emissions (t/an/km2) : longitude latitude ",10A5)' ) var_emis_name
      write(*,'("*** info : continue le calcul...")')
   end if
!--------------------------------------------------------------------                   
! Emissions CHIMERE (FACULTATIF)
!--------------------------------------------------------------------    
! Lit les emissions de la SNAP 2, 4 et 7
   if (iemis1.and.inc_emis) then 
    emis_g0 = 0.
    if (idebug)write(*,*) 'Lit les emissions NetCDF du domaine COARSE...' 
    do ific=1, nfile_emis
    call read_chimere_emis_cdf( femis1(ific), it, nx0, ny0, nsnap1, dx1min, dy1min, &
                                nvar_emis, var_emis_name, ivar_emis, tbl_ratio_emi_snap, emis_g0(:,:,:), &
				iemis_1h, idebug  )
    end do
    do ivar=1, nvar_emis
      if (idebug)write(*,*) 'Domaine COARSE emissions max.(kg/an) '//trim(var_emis_name(ivar)),maxval(emis_g0(ivar,:,:))
    end do      
   end if  
!--------------------------------------------------------------------        
! Ouvre les fichiers de chimie et verifie la presence de PM10
!--------------------------------------------------------------------  
   if (idebug)write(*,*) 'Lit les fichiers NetCDF pour presence PM10/meteo/AEROSOLS'
   call check_cdf
   
   !if (idebug)write(*,*) 'Defini les variables en entrees/sorties'    		  
!-------------------------------------------------------------------- 
! Determination des indices des polluants pour la chimie
!--------------------------------------------------------------------         
   if (idebug)write(*,*) 'Cherche les identifiants des variables dans raffinement CHIMIE+AEROSOLS' 
   
   ! Concentrations 3D
   call check_cdf_var
   
   ! Emissions  
   do ivar=1, nvar_emis
    if (trim(var_emis_name(ivar)).eq.'VOC'     )emis_icovnm    = ivar
    if (trim(var_emis_name(ivar)).eq.'NOX'     )emis_inox      = ivar
    if (trim(var_emis_name(ivar)).eq.'NO'      )emis_ino       = ivar
    if (trim(var_emis_name(ivar)).eq.'NO2'     )emis_ino2      = ivar        
    if (trim(var_emis_name(ivar)).eq.'PM10'    )emis_ipm10     = ivar
   end do
   
   if (idebug) then
   write(*,*)'Conc NO    code',chim_ino,'is ok?',iinput_var(chim_ino)
   write(*,*)'Conc NO2   code',chim_ino2,'is ok?',iinput_var(chim_ino2)
   write(*,*)'Conc O3    code',chim_io3,'is ok?',iinput_var(chim_io3)
   write(*,*)'Conc PM10  code',chim_ipm10,'is ok?',iinput_var(chim_ipm10)
   write(*,*)'Conc PPM   code',chim_ippm,'is ok?',iinput_var(chim_ippm)
   write(*,*)'Conc PHNO3 code',chim_iphno3,'is ok?',iinput_var(chim_iphno3)
   write(*,*)'Temp temp  code',met_itemp,'is ok?',iinput_var(met_itemp)
   write(*,*)'Vent wins  code',met_iwins,'is ok?',iinput_var(met_iwins)
   
   write(*,*)'Emis COVNM code',emis_icovnm
   write(*,*)'Emis NOx   code',emis_inox
   write(*,*)'Emis PM10  code',emis_ipm10
   end if
      
!-------------------------------------------------------------------- 
! Lit les coordonnées latl/long, la topographie, le landuse dans un fichier
!-------------------------------------------------------------------- 
  if (iground1.and.idebug)write(*,*) 'Lit lon,lat,topo,landuse COARSE dans '//trim(fground1)	
  if(iground1)call read_chimere_ground_cdf(fground1, nx0, ny0, lon_g0, lat_g0, landuse_g0, topo_g0, idebug )  
!-------------------------------------------------------------------- 
! Calcul la correspondance des grilles
!-------------------------------------------------------------------- 
   if (idebug)write(*,*) 'Calcul des correspondances de grilles...'   
   call calcul_correspondance    
   if (idebug)write(*,*) 'Calcul des poids de l interpolation...'
   call calcul_poids
!-------------------------------------------------------------------- 
! Moyenne la topologie et les emissions FINE -> COARSE
!--------------------------------------------------------------------    
! calcule les valeurs de topo_g0 a partir de topo_g2
  if (.not.iemis1.and.idebug)write(*,*) 'Somme les emissions grille FINE -> grille COARSE'
  if(.not.iemis1)call somme_var( nvar_emis, nx2, ny2, nx0, ny0, ix2, iy2, emis_g2(:,:,:), emis_g0(:,:,:) )
  
  if (.not.iground1.and.idebug)write(*,*) 'Moyenne la topographie grille FINE -> grille COARSE'
  if(.not.iground1)call moyenne_2d(         nx2, ny2, nx0, ny0, ix2, iy2, topo_g2(:,:)  , topo_g0(:,:)   )                      	                        
!--------------------------------------------------------------------        
! Creer un NetCDF de sorties
!-------------------------------------------------------------------- 
   if (idebug)write(*,*) 'Creer le NetCDF de sortie...'	
   call create_netcdf
   write(*,*) '>>>>>>>>>>>> Creation du NetCDF de sortie OK'
   write(*,*) '-> '//trim(fout2)
   if (idebug)write(*,*) 'Creation OK'			 
!---------------------------------------------------------------
! Calcul des profils d'altitude 
!---------------------------------------------------------------      
   if (idebug)write(*,*) 'Initialisation des profils verticaux...'	 
   call definition_profil_z  

!-------------------------------------------------------------------- 
! Calcul des ratios d'emissions
!-------------------------------------------------------------------- 
!   if(iraf_chim) then
!    if (idebug)write(*,*) 'Calcul du ratio des emissions...'          	 
!    call calcul_ratio_emissions     
!   end if		   
!--------------------------------------------------------------------        
! Ouvre les fichiers en lecture ou lecture/ecriture
!-------------------------------------------------------------------  
   if (idebug)write(*,*) 'Ouvre les fichiers NetCDF chimie et meteo'
   call open_cdf     
!--------------------------------------------------------------------       
!  BOUCLE SUR LES HEURES
!--------------------------------------------------------------------            
   do it=1,numFrTimes
      !-------------------------------------------------------------------- 
      ! Lit les emissions horaires (NetCDF seulement)
      !-------------------------------------------------------------------- 
      if(iraf_chim) then 
	if (inc_emis.and.iemis_1h) then
          emis_g2 = 0.
          if (idebug)write(*,*) 'Lit les emissions NetCDF du domaine FINE...'    
	  ific=1	
	  call read_chimere_emis_cdf( femis2(ific), it, nx2, ny2, nsnap1, 1, 1, nvar_emis, var_emis_name, ivar_emis, tbl_ratio_emi_snap , emis_g2(:,:,:), iemis_1h, idebug  )
          !Calcul les ratios NOx
	  call calcul_emissions_nox_pm
	end if
      end if

      !-------------------------------------------------------------------- 
      ! Ecrit les emissions dans fichier de sortie
      !-------------------------------------------------------------------- 
      !if(iemis_1h.and.idebug)write(*,*) 'Ecrit emissions horaires dans fichier sortie...'	 
      !if(iemis_1h.and.idebug)call write_emis_cdf
      if(idebug)write(*,*) 'Ecrit emissions horaires dans fichier sortie...'	 
      if(idebug)call write_emis_cdf      
   
      !-------------------------------------------------------------------- 
      ! Calcul des ratios d'emissions
      !-------------------------------------------------------------------- 
      if(iraf_chim) then
        if (idebug)write(*,*) 'Calcul du ratio des emissions...'          	 
        call calcul_ratio_emissions     
      end if  
   
      !---------------------------------------------------------------
      ! Recalcul des variables
      !---------------------------------------------------------------   
       call read_cdf_1h
      
      !---------------------------------------------------------------
      ! Recalcul de la temperature 2-m avec hypothese adiabatique (si meteo)
      !---------------------------------------------------------------
      if (imeteo) then   
        if (idebug)write(*,*) 'Calcul de la temperature 2 metres...'
        call calcul_temperature_2m
      end if
            
      if (iraf_chim) then
      !> PROFILS EMIS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      iinterp = .false. ! utilise interpolation ou recopy ?
      !-------------------------------------------------------------------- 
      ! Interpolation des données geographiques de grille 0 -> grille 1
      !--------------------------------------------------------------------
       if (idebug)write(*,*) 'Interpolation des donnees grille 0 -> grille 1... iinterp=',iinterp
       ismooth=.false.
       call interpolation_topo_grille
      !---------------------------------------------------------------
      ! Recalcul des variables
      !---------------------------------------------------------------   
      !if (it.eq.1) write(*,*) 'Recalcul des variables alti et pollution...'  
      if (idebug)write(*,*) 'Recalcul des variables alti et pollution grille 0 -> grille 1... iinterp=',iinterp
      call interpolation_var_1h       		                
      !---------------------------------------------------------------
      ! Calcul des concentrations avec emissions GAZ
      !---------------------------------------------------------------
      if (correct_no.or.correct_o3.or.correct_no2) then      
      write(*,*) '>>>>>>>>>>>> Calcul du raffinement de maillage GAZ...'
      call rafinement_gaz(nx1,ny1,nz1,&
      emis_g2(emis_ino,:,:)      ,emis_g2(emis_ino2,:,:), &
      ratio_emis_g2(emis_ino,:,:),ratio_emis_g2(emis_ino2,:,:), &
      val_g1_var (chim_ino,:,:,:),val_g1_var (chim_ino2,:,:,:),val_g1_var (chim_io3,:,:,:), &
      dvar_de_loc(chim_ino,:,:  ),dvar_de_loc(chim_ino2,:,:  ),dvar_de_loc(chim_io3,:,:  ), &
      correct_no,correct_no2,correct_o3, idebug )
      end if
      !---------------------------------------------------------------
      ! Calcul des concentrations avec emissions AEROSOLS
      !---------------------------------------------------------------
      if (ipm.and.ippm.and.correct_pm10) then
      write(*,*) '>>>>>>>>>>>> Calcul du raffinement de maillage AEROSOLS...'
      call rafinement_aerosols(nx1,ny1,nz1,ratio_emis_g2(emis_ipm10,:,:), &
      val_g1_var (chim_ipm10,:,:,:),val_g1_var (chim_ippm,:,:,:),dvar_de_loc(chim_ipm10,:,:), idebug)      
      end if
      !> PROFILS EMIS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                
      end if  

      if (iraf_alti) then
      !> PROFILS Z>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      iinterp = .true. ! utilise interpolation ou recopy ?
      write(*,*) '>>>>>>>>>>>> Calcul du raffinement de maillage ALTITUDE...'         
      !-------------------------------------------------------------------- 
      ! Interpolation des données geographiques de grille 0 -> grille 1
      !--------------------------------------------------------------------
       if (idebug)write(*,*) 'Interpolation des donnees grille 0 -> grille 1... iinterp=',iinterp
       ismooth=.true. ! smooth_ratio -> params.f90
       call interpolation_topo_grille
      !---------------------------------------------------------------
      ! Recalcul des variables
      !---------------------------------------------------------------   
      !if (it.eq.1) write(*,*) 'Recalcul des variables alti et pollution...'  
      if (idebug)write(*,*) 'Recalcul des variables alti agl et pollution grille 0 -> grille 1... iinterp=',iinterp
      call interpolation_var_1h
      !---------------------------------------------------------------
      ! Recalcul de l'alti niveau mer
      !---------------------------------------------------------------         
      if (idebug) write(*,*) 'Calcul de l altitude asl...'
      if (it.le.2) call calcul_altitude_asl
      !---------------------------------------------------------------
      ! allocation des mailles X, Y, Z au profil de Z (it=1)
      !---------------------------------------------------------------
      if (idebug)write(*,*) 'Allocation des mailles X, Y, Z au profil de Z'
      if (it.eq.1)call allocation_altitude_maille
      !---------------------------------------------------------------
      ! Calcul des profils MOYEN d'ozone, de PM10, de NO2
      !---------------------------------------------------------------
      if (idebug)write(*,*) 'Calcul des profils moyens des variables...'      
      ! calcul du profil d'ozone moyen sur le domaine en ignorant le niveau 1 dans le calcul
      profil_1D_var    = 0.   
      call calcul_profil_moyen      
      !---------------------------------------------------------------
      ! Calcul des profils VOISINAGE d'ozone, de PM10, de NO2
      !---------------------------------------------------------------
      ! calcul du profil d'ozone moyen sur le domaine en ignorant le niveau 1 dans le calcul
      profil_3D_var   = 0.    	 
      if (idebug)write(*,*) 'Calcul des profils sur voisinage local...'   
      call calcul_profil_voisinage                          
      !---------------------------------------------------------------
      ! Calcul des concentrations avec profils Z
      !---------------------------------------------------------------
      if (idebug)write(*,*) 'Calcule les profils Z sur la grille fine...'            
      call calcul_var_z
      !> PROFILS Z>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      end if

      !---------------------------------------------------------------       
      ! Calcul des concentrations finales
      !---------------------------------------------------------------
      val_g2_var(:,:,:,1)  = val_g1_var(:,:,:,1) + dvar_dz_loc(:,:,:) + dvar_de_loc(:,:,:)
      
      if (idebug) then
        do ivar=1, nvar
          write(*,'(A8," COA=",F10.2," DZ=",F10.2" DE=",F10.2)') &
	     trim(var_name(ivar)),val_g1_var(ivar,ixp,iyp,1),dvar_dz_loc(ivar,ixp,iyp),dvar_de_loc(ivar,ixp,iyp)
        end do        
      end if
      
      do i=1,nx2
      do j=1,ny2
      do ivar=1, nvar
        if (val_g2_var(ivar,i,j,1).lt. 0.0)val_g2_var(ivar,i,j,1) = 0.0
      end do
      end do
      end do      
      !---------------------------------------------------------------       
      ! Calcul FDMS
      !---------------------------------------------------------------           
      if (ipm)then      
        write(*,*) '>>>>>>>>>>>> Calcul des ecarts FDMS sur la grille...'
	if (met_itemp.ne.0) then 
	 if (icorrect_z_var(met_itemp)) then
	  tem2_g2(:,:) = val_g2_var(met_itemp,:,:,1)
          if (idebug)write(*,*) 'Trouve temperature 2m interpolee:',maxval(tem2_g2)
	 end if	  
	end if         
        call  calcul_fdms(nx2, ny2, val_g2_var(chim_ipm10,:,:,1), val_g2_var(chim_iphno3,:,:,1),&
                          tem2_g2(:,:), pm10_nonvolat_g2(:,:,1), ecart_fdms_g2(:,:,1), idebug )     
      end if      
      !---------------------------------------------------------------       
      ! Calcul vent 10 mètres (a mettre a la fin)
      !---------------------------------------------------------------           
      if (iwin) then      
        write(*,*) '>>>>>>>>>>>> Calcul du facteur de vent 10 m...'
        call  calcul_vent_10m
        do i=1,nx2
        do j=1,ny2
	  val_g2_var(met_iwins,i,j,1) = w10m_g2(i,j) ! vent wins -> vent 10 m
	  w10m_g2(i,j) = w10m_g2(i,j) / w10m_g1(i,j) ! vent 10 m -> facteur *
        end do
        end do	
        if (idebug)write(*,*) 'Trouve vent interpole:',maxval(w10m_g2)
      end if            
      !---------------------------------------------------------------       
      ! ecrit dans un NetCDF pour le pas de temps it
      !--------------------------------------------------------------- 
      call write_cdf_1h       	   
   end do ! boucle sur les iterations 
!--------------------------------------------------------------------       
!  FIN DE LA BOUCLE SUR LES HEURES
!--------------------------------------------------------------------         
           
!---------------------------------------------------------------
! Ferme les NetCDF
!---------------------------------------------------------------      
    call check(nf90_close(out1fileID))
    if (imeteo) call check(nf90_close(met1fileID))      
      
99  continue      

   
      end program
   
