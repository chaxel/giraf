
      module params_chimie
      
      real, parameter :: taux_no_nox  = 0.92
      real, parameter :: taux_no2_nox = 1.0 - taux_no_nox 
      
      
      end module params_chimie

!--------------------------------------------------------------------        
      subroutine read_chimere_emis_cdf(femis, nx, ny, nsnap1, dxmin, dymin, nvar_emis, var_emis_name, ivar_emis, tbl_ratio_emi_snap,  emis, idebug )

      use netcdf
      use typesizes
      
      implicit none        
      
      character(len=*)   :: femis
      integer            :: nx, ny, nsnap1, dxmin, dymin, nvar_emis
      character(len=*)   :: var_emis_name(nvar_emis)
      logical            :: ivar_emis(nvar_emis)
      real               :: emis(nvar_emis, nx, ny)
      logical            :: idebug

! fichiers d'emissions SNAP 3 
      integer            :: ncode_snap1,ncode_snap2, snap1
      integer            :: ncFileID, xVarID, snapVarID, snapDimID
      integer,allocatable:: tab_code_snap(:)
      integer            :: ns, ivar
      
      ! local
      integer            :: start3(3), count3(3)      
      real               :: emis_tmp(nvar_emis, nx, ny)            
      real               :: tbl_ratio_emi_snap(nsnap1) ! ratios d'emissions au sol pour snap 1-10      
      
                            
      ! Tente d'ouvrir le fichier d'emissions     
      if ( nf90_open(femis, nf90_nowrite, ncFileID).eq.nf90_noerr ) then  
            
        ! affiche le tableau des ratios par SNAP 1
        if (idebug)write(*,'("tbl_ratio_emi_snap=",11F5.2)')tbl_ratio_emi_snap  
	        
        call check(nf90_inq_dimid(ncFileID, 'snap'  , snapDimID))         
        call check(nf90_Inquire_Dimension(ncFileID, snapDimID , len=ncode_snap1))	      
        
	if (idebug)write(*,*) trim(femis)//' - >trouve nombre de SNAP ',ncode_snap1	           
	
	allocate( tab_code_snap(ncode_snap1))    
        !write(*,*) 'Lit les codes SNAP '             
        call check(nf90_inq_varid(ncFileID,'code_snap',snapVarID))
        call check(nf90_get_var(ncFileID,snapVarID,tab_code_snap))	
        do ns = 1, ncode_snap1  
         !initialisation
         emis_tmp(:,:,:) = 0.	 
	 start3 = (/dxmin,dymin, ns/)
         count3 = (/   nx,   ny, 1 /)       
         do ivar= 1, nvar_emis   	    
           if (nf90_inq_varid(ncFileID, trim(var_emis_name(ivar)), xVarID).eq.nf90_noerr) then	 
             ivar_emis(ivar) = .true.
	     call check(nf90_get_var(ncFileID, xVarID, emis_tmp(ivar,:,:),start = start3, count = count3 ))
	   else
	     ivar_emis(ivar) = .false.
	     emis_tmp(ivar,:,:) = 0.
	   end if     
         end do       
         snap1=int(tab_code_snap(ns)/10000)       
         if (idebug.and.tbl_ratio_emi_snap(snap1).ne.0.)write(*,'("Lit les emissions de la SNAP",I7," facteur*",F4.2)') tab_code_snap(ns),tbl_ratio_emi_snap(snap1)
	! aditionne les emissions         
	 emis = emis + emis_tmp * tbl_ratio_emi_snap(snap1)      
        end do 
        call check(nf90_close(ncFileID))
        deallocate(tab_code_snap)	     
      else     
        !Erreur fatale (version 1.3)
	!if (idebug)&
	write(*,*) '***ERREUR : pas d emissions disponibles dans '//trim(femis)//'>>>>>>>>>>' 
	stop 1  	  
      end if
      
      end subroutine 
!----------------------------------------------------------------------------------------------------------------
      subroutine calcul_ratio_emissions
      
      use params
      
      implicit none
      
      !local
      integer :: i2, j2, ix, iy, ivar
      real, parameter :: zerolimit = 1.0E-5
       
      ratio_emis_g0  = 0.
      ratio_emis_g2  = 0.
      nratio_emis_g0 = 0.
           
      do i2=1,nx2
        do j2=1,ny2	
	  ix	= ix2(i2,j2)
	  iy	= iy2(i2,j2)	  	  
          do ivar=1,nvar_emis
	    if ( emis_g0(ivar,ix,iy).ge. 0.01/36.*distgrid_km*distgrid_km ) then ! limite 10 kg/an/mailles (teste avec maille 6 km)
              nratio_emis_g0(ivar,ix,iy) = nratio_emis_g0(ivar,ix,iy) + 1
	      if ( emis_g0(ivar,ix,iy) .ge. zerolimit ) then
	        ratio_emis_g2(ivar,i2,j2) = emis_g2(ivar,i2,j2)  /  emis_g0(ivar,ix,iy)
	      else
	        ratio_emis_g2(ivar,i2,j2) = 1.00
	      end if
	      
	    else
	      nratio_emis_g0(ivar,ix,iy) = nratio_emis_g0(ivar,ix,iy) + 1
              ratio_emis_g2(ivar,i2,j2) = 1.00
	    end if

            ! Normalisation
	    ratio_emis_g0(ivar,ix,iy) = ratio_emis_g0(ivar,ix,iy) + ratio_emis_g2(ivar,i2,j2) 
	     	     
          end do
	end do
      end do
 
      ! Normalisation et ecriture dabns fichier
      open(unit=10,file='ratio_emissions.txt',status='unknown')
      write(10,'(5A10)') 'lon','lat','ratio_nox','ratio_cov','ratio_pm'
      do i2=1,nx2
        do j2=1,ny2	
	  ix = ix2(i2,j2)
	  iy = iy2(i2,j2)      
          do ivar=1,nvar_emis
	    if ( ratio_emis_g0(ivar,ix,iy).ge. zerolimit ) then 	  
	       ratio_emis_g2(ivar,i2,j2) = ( ratio_emis_g2(ivar,i2,j2) / ratio_emis_g0(ivar,ix,iy) ) * nratio_emis_g0(ivar,ix,iy)
	    else
	       ratio_emis_g2(ivar,i2,j2) = 1.00
	    end if
            ! valeur inferieure a la zero limit
	    if ( ratio_emis_g2(ivar,i2,j2).lt. zerolimit ) ratio_emis_g2(ivar,i2,j2) = 0.0
	  end do
	  	  	  
	  write(10,'(2F10.3,3F10.4)') lon_g2(i2,j2),lat_g2(i2,j2),ratio_emis_g2(1:nvar_emis,i2,j2)	
	  !if (idebug) write(*,*)   lon_g2(i2,j2),lat_g2(i2,j2),ratio_emis_g2(1:nvar_emis,i2,j2)	
	end do
      end do
      
      close(10)
      if (idebug)write(*,*) 'Emissions ration OK... ecrit ratio_emissions.txt'
      
      end subroutine

!----------------------------------------------------------------------------------------------------------------
      subroutine rafinement_gaz(nx,ny,nz,&
      ratio_emis_nox, no_g1,no2_g1, o3_g1,&
      dno, dno2, do3,&
      correct_no, correct_no2, correct_o3, idebug)
     
      use params_chimie
      
      implicit none   

      integer :: nx,ny,nz    ! domaine CHIMERE rafin�

! Corrections      
      logical :: correct_no 
      logical :: correct_o3 
      logical :: correct_no2
      logical :: idebug

      integer :: nh
      integer :: emis_inox 
      real    :: ratio_emis_nox(nx,ny)
      real    :: dno(nx,ny),  dno2(nx,ny), do3(nx,ny) ! varaitions sur niveau 1 apport�e par raffinement
      real    :: no_g1(nx,ny,nz) ,  no_g2(nx,ny,nz)
      real    :: no2_g1(nx,ny,nz),  no2_g2(nx,ny,nz)
      real    :: o3_g1(nx,ny,nz) ,  o3_g2(nx,ny,nz)
             
      ! LOCAL
      integer :: i, j
      real    :: tem2, sreh
      real    :: dno_loc, dno2_loc
      character(len=10) :: actifstr
                           
      if(idebug)write(*,*) 'Entre dans la routine'
      if(idebug)write(*,*) 'Dimensions nx,ny,nz=',nx,ny,nz
      if(idebug)write(*,*) 'Taux NO2/NOx=',taux_no2_nox
      !---------------------------------------------------------------
      ! Initialisation
      !---------------------------------------------------------------   
       no_g2(:,:,:)  = no_g1(:,:, :)
       no2_g2(:,:,:) = no2_g1(:,:,:) 
       o3_g2(:,:,:)  = o3_g1(:,:,:) 	

      !---------------------------------------------------------------
      ! ReCalcul du NO
      !---------------------------------------------------------------             
      ! Reaffectation du NO sur la maille de 1 km = le surplus dans la premiere maille	   
      if (correct_no) then
        actifstr='ACTIVE'
        do i = 1, nx
          do j = 1, ny  	 

	    if ( no_g1(i,j,1) .gt. no_g1(i,j,2) ) then ! dans panache ? -> NON car niveau 1 > niveau 2             	     
	     ! Corrige le NO au sol qui est d'origine primaire (niveau 2 - niveau 1)
	     dno_loc = ( no_g1(i,j,1) - no_g1(i,j,2) ) * ratio_emis_nox(i,j)	            
	    else	    
	     dno_loc = 0.	    
	    end if
            ! corrige le NO	    
	    no_g2(i,j,1)  =  no_g1(i,j,1) + dno_loc	     
	    ! Utilise le rapport d'emission en NO/NO2 = 0.08 NO/NOx = 0.08 pour corriger le NO2 primaire
	    no2_g2(i,j,1) = no2_g1(i,j,1) + dno_loc * (taux_no2_nox/taux_no_nox)	    	    	    	    	    
	  end do
         end do	     	     	     
      else
        actifstr='ARRETE' 	     	   	   
      end if !correct_no      
      if (idebug)write(*,*) '*** traitement 1.a-'//trim(actifstr)//': rafinement de NO avec emissions primaires de NOx'          

      !---------------------------------------------------------------
      ! ReCalcul du NO2
      !---------------------------------------------------------------             	  
      ! Reaffectation du NO2 sur la maille de 1 km (ne pas utiliser, conduit a trop de NO2)
      if (correct_no2) then
        actifstr='ACTIVE'          
        do i=1,nx
          do j=1,ny		
	    if ( no2_g1(i,j,1) .gt. no2_g1(i,j,2) ) then ! dans panache ? -> NON car niveau 1 > niveau 2
	     ! meme methode que NO2
	     dno2_loc = ( no2_g1(i,j,1) - no2_g1(i,j,2) ) * ratio_emis_nox(i,j)
	     no2_g2(i,j,1) =  no2_g1(i,j,2) + dno2_loc
	   !else
	   !  no2_g2(i,j, 1) = no2_g1(i,j,1) ! dans panache ? -> OUI NE FAIT RIEN
	   end if
	  end do
        end do  		 
      else
        actifstr='ARRETE'        	 
      end if	   
        if (idebug)write(*,*) '*** traitement 1.b-'//trim(actifstr)//': rafinement de NO2 avec emissions primaires de NOx (deconseille)'    
      !---------------------------------------------------------------
      ! Correction du O3
      !--------------------------------------------------------------- 
      if (correct_o3.and.correct_no2) then
       actifstr='ACTIVE'  
       do i= 1, nx
        do j= 1, ny	   	   
	   if ( o3_g2(i,j,1) - (  no2_g2(i,j,1) - no2_g1(i,j,1) ) .gt.0 ) then
	     o3_g2(i,j,1) =  o3_g2(i,j,1) - (  no2_g2(i,j,1) - no2_g1(i,j,1) )	     
	   else
	     o3_g2(i,j,1) = 0.
	   end if  	   	   
	end do
       end do 
      else
        actifstr='ARRETE'     
      end if 
      
      if (idebug)write(*,*) '*** traitement 1.c-'//trim(actifstr)//': NO2 apparait <-> O3 disparait (deconseille)'                 
      
      !---------------------------------------------------------------
      ! 2nd traitement : relation O3 + NO -> NO2 + O2
      !---------------------------------------------------------------                
      if (correct_no)  then
      actifstr='ACTIVE'           
      do i=1,nx
        do j=1,ny	   
	   ! Reactions lineaires de titration de O3 par NO (actif si on corrige le NO auparavant
	   ! O3 + NO -> NO2

	   ! Conversion de NO en NO2 :  O3 + NO -> NO2 mais pas de conversion NO2 -> NO + O3
	   ! Agit sur la quantite de NO du premier niveau

	   dno_loc = no_g2(i,j,1) - no_g1(i,j,1)

	   if ( dno_loc .gt. 0. ) then

	   if ( o3_g2(i,j,1) - dno_loc .gt.0 ) then

	     ! tout le NO est converti en NO2 :  O3 + NO -> NO2
	     o3_g2 (i,j, 1) = o3_g2 (i,j, 1) - dno_loc
	     no2_g2(i,j, 1) = no2_g2(i,j, 1) + dno_loc
	          
	   else

	     ! une partie dno du NO est convertie en NO2 :  dno*O3 + NO -> dno*NO2 + (1-dno)*NO	     
	     dno_loc = o3_g2(i,j, 1)
	     o3_g2 (i,j, 1) = 0.
	     no2_g2(i,j, 1) = no2_g2(i,j,1) + dno_loc
	     no_g2 (i,j, 1) = no_g2 (i,j,1) - dno_loc
	     
	   end if

	   end if
	   
	   ! remplace le NO par du NO2 en presence d'ozone
	   !if ( ( o3_g2(i,j, 1).gt.0 ) .and. (correct_no2) ) then
	   
	   !if ( no2_g2(i,j,1) + no_g2(i,j,1) - no_g1(i,j,1)  .gt.0 ) then	   	     
	   !  no2_g2(i,j,1) = no2_g2(i,j,1) + no_g2(i,j,1) - no_g1(i,j,1) 
	   !else
	   !  no2_g2(i,j,1) = 0.
	   !end if
	  ! 
	   !end if
	   
	   ! Ajoute une condition sur l'ozone au sol, il est toujours inferieur au niveau 2
	   !if ( o3_g2(i,j,1).gt.o3_g1(i,j,2) ) o3_g2(i,j,1) = o3_g1(i,j,2)
	   
	   !end if ! o3_g2(i,j, 1) .gt. 1.
	   
	   !end if ! si NO[k=2] > NO[k=1] : mode panache, ne fait rien !
	   
	end do
      end do 

      !write(*,*) '... recalcul O3, NO, NO2 avec emissions OK'    
      else
         actifstr='ARRETE'    
      end if

      if (idebug)write(*,*) '*** traitement 2.0-'//trim(actifstr)//': applique relation O3 + NO -> O2 + NO2' 
            
      dno =  no_g2(:,:,1) -  no_g1(:,:,1)
      do3 =  o3_g2(:,:,1) -  o3_g1(:,:,1)      
      dno2= no2_g2(:,:,1) - no2_g1(:,:,1)      
           
      !no_g1 = no_g2
      !o3_g1 = o3_g2
      !no2_g1= no2_g2      

      end subroutine rafinement_gaz
!----------------------------------------------------------------------------------------------------------------
      subroutine rafinement_aerosols(nx,ny,nz,ratio_emis_ppm, pm10_g1 ,ppm_g1, dpm10, idebug )
      
      implicit none  
      
      ! entrees
      integer ::  nx,ny,nz  
      real :: ratio_emis_ppm(nx,ny)
      real :: pm10_g1(nx,ny,nz), pm10_g2(nx,ny,nz) 
      real :: ppm_g1(nx,ny,nz),  ppm_g2(nx,ny,nz)
      logical :: idebug
      
      ! sorties
      real :: dpm10(nx,ny)
       
      ! LOCAL
      integer :: i, j
      real    :: dppm_loc
 
      if(idebug)write(*,*) 'Entre dans la routine'
      if(idebug)write(*,*) 'Dimensions nx,ny,nz=',nx,ny,nz                           

      !---------------------------------------------------------------
      ! Verification de la presence de valeurs dans ppm_g1
      !---------------------------------------------------------------  
      if(idebug) then
      if ( maxval(ppm_g1(:,:, 1)) .lt. 1.E-5)then
        write(*,*) '***erreur: ppm_g1 = 0'
        stop
      end if      
      end if
      !---------------------------------------------------------------
      ! Initialisation
      !---------------------------------------------------------------   
       pm10_g2(:,:,:) = pm10_g1(:,:, :)
       ppm_g2(:,:,:)  = ppm_g1(:,:,:) 
      !---------------------------------------------------------------
      ! Recalcul du PM10
      !---------------------------------------------------------------            
      if (idebug)write(*,*) 'Calcul du PM10'      
     
      do i=1,nx
        do j=1,ny
	   
	   ! base le calcul sur les PPM (PPM_fin + PPM_coa + OCAR_fin + BCAR_fin)
	   dppm_loc = ppm_g1(i,j,1) -  ppm_g1(i,j,2)
	   if (  dppm_loc  .gt. 0 ) then
	     ppm_g2(i,j,1) =  ppm_g1(i,j,2) + dppm_loc * ratio_emis_ppm(i,j) 
           else 
             ppm_g2(i,j,1) =  ppm_g1(i,j,1)  
	   end if

	   ! PM10
           pm10_g2(i,j,1) = pm10_g1(i,j,1) - ppm_g1(i,j,1) + ppm_g2(i,j, 1)	            
	end do
      end do     
   
      dpm10= pm10_g2(:,:,1)  - pm10_g1(:,:,1)              
      pm10_g1 = pm10_g2 
  
      end subroutine rafinement_aerosols

!----------------------------------------------------------------------------------------------------------------
      subroutine calcul_fdms(nx,ny, pm10_g2 ,phno3_g2, tem2_g2, pm10_nonvolat_g2 ,  ecart_fdms_g2, idebug )
            
      implicit none   

      integer ::  nx,ny

      ! entrees            
      real :: pm10_g2(nx,ny)
      real :: phno3_g2(nx,ny)
      real :: tem2_g2(nx,ny)
      logical :: idebug
            
      ! sorties
      real    :: pm10_nonvolat_g2(nx,ny) ,  ecart_fdms_g2(nx,ny)      
       
      ! LOCAL
      integer :: i, j
      real    :: tem2, sreh, phno3, pm10, dmassm, nh4no3, pm10_nonvolat, pm10_ecart_fdms
      
      real    :: fac_teom
         
      !---------------------------------------------------------------
      ! Pour les polluants NO et PM10 on fait la moyenne de la correction et de la conc initiale
      !--------------------------------------------------------------- 
      if (idebug)write(*,*) 'Calcul de ecart FDMS avec methode Militon (2009)'       
      
      do i= 1, nx
        do j= 1, ny
	   tem2 = tem2_g2(i,j)  
	   if ((tem2.lt.200).or.(tem2.ge.350))tem2=288. ! si meteo disponible	   
	   tem2  = tem2 - 273.15                  ! conversion K-> degresC
	   fac_teom = 0.2 + 0.4 * exp(-0.1*tem2)
	   phno3    = phno3_g2(i,j)
	   pm10     = pm10_g2 (i,j)
	   dmassm   = (14.+ 4*1.+ 14.+ 3*16.) / (14.+ 3*16.) ! rapport des masses molaires NH4NO3/NO3   
	   
	   ! evaluation INERIS
	   nh4no3   = dmassm * phno3	   	   
	   
	   ! PM10 non volat 50�C (Chaxel, 2010)
	   pm10_nonvolat = ( pm10 - nh4no3 ) / (1 + fac_teom)
	   if (  pm10_nonvolat .lt. 0. ) pm10_nonvolat = 0.
	   	   
	   ! ecart FDMS calcul� (Militon, 2009)
	   pm10_ecart_fdms = pm10 - pm10_nonvolat	   	   
           if (  pm10_ecart_fdms .lt. 0. ) pm10_ecart_fdms = 0.
	   
	   pm10_nonvolat_g2(i,j) = pm10_nonvolat
	   ecart_fdms_g2(i,j)    = pm10_ecart_fdms
	   
	end do
      end do  
      if (idebug)write(*,'("PM10 max.=",F7.2," ug/m3")') maxval(pm10_g2)	   
      if (idebug)write(*,'("PHNO3 max.=",F7.2," ug/m3")') maxval(phno3_g2)    
      if (idebug)write(*,'("PM10 non volat max.=",F7.2," ug/m3")') maxval(pm10_nonvolat_g2)	   
      if (idebug)write(*,'("Ecart max.=",F7.2," ug/m3")') maxval(ecart_fdms_g2)
      
      	   ! Calcul de la correction du d�bit des tubes pour temperature > 15�C (selon formule du LCSQA, 2001)
	   !tem2 = 273.15
	   !sreh = 80	   
	   !if (imeteo) tem2 = tem2_g2(i,j)		      
	   !tem2 = tem2 - 273.15 ! K -> degresC 	   
	   !if ( tem2 .gt. 15 ) then
	   !  fac = 2.77e-4*sreh**2-3.33e-3*tem2**2+3.38e-4*tem2*sreh-3.6e-2*sreh+9.98e-2*tem2+1.135
	   !else
	   !  fac = 2.77e-4*sreh**2-3.33e-3*15**2+3.38e-4*15*sreh-3.6e-2*sreh+9.98e-2*15+1.135
	   !end if   	   
	   ! version 8.1 
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	   
	   ! regression obtenue a Lyon centre sur 2007/2008 sur moyenne jour
	   !ecart_fdms_g2(i,j) = -5.27835 + 0.05736* sreh2(i,j)+ 0.22519*pm10_g2(i,j,1) + 1.33492 * phno3_g1(i,j,1)	   
	   ! regression obtenue a Lyon centre sur janvier 2009 (version 8)
	   ! Principe : 
	   ! la correction FDMS est la composante d'une composante nitrate et d'uine partie volatile de primaire
	   ! coefficient nitrate 1.647 (Bessagnet, LCSQA, 2007)
	   ! coefficent PM10 = 0.4 (Chaxel, Atmo RA, 01/2009) ATTENTION il semble dependre de la temperature !!!	   
	   ! version 8.2 
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! Pour calcul 2008 : ecart FDMS calcul� avec r�gression sur Lyon
           ! ecart_fdms_g2(i,j) = 0.4 * pm10_g2(i,j,1) + 1.647 * phno3_g1(i,j,1)
	   ! version 8.3
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! EN PREMIERE APPROXIMATION : l'ecart FDMS est juste des nitrates (version 8.3)
	   !ecart_fdms_g2(i,j) = 1.647 * phno3_g1(i,j,1)
	   ! version 9.0 : d'apres travail d'Alice MILITON (juillet 2009)
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   ! corr_fdms = 1.647 * pHNO3 + fac * TEOM + e (e=0 dans cette version)
	   ! resoud l'equation	   
	   ! facteur de la temperature tem2_g2 est en �C
	   
      
      end subroutine calcul_fdms