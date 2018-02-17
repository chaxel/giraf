      module params_physique
      
       ! local a physique
       real    :: topo
       integer :: ix, iy, iz, k1, ivar, i2, j2, k2
       real, parameter :: gradient_adiabatique = -0.006 ! K/m
       
      end module params_physique
!---------------------------------------------------------------------------------------------------------------------     
      subroutine definition_profil_z

      use params
      use params_physique      
      
      implicit none
      
      real :: altitude_min
      
      ! Calcul de l'altitude mini du domaine FINE
      altitude_min = minval(topo_g2)
      if (idebug)write(*,*) 'Altitude mini(m)=' , int(altitude_min)
                   
      ! Anciens profils non orthonorme
      if ( nz1.eq.8 ) then
        alti_1D = (/ 25.,  75.,  150., 500., 1500., 2500., 3500., 4500. /)
		     
      else if ( nz1.eq.12 ) then
        alti_1D = (/ 25.,  75.,  150.,  250., 500. , 750.,  &
                     1250., 2000., 2500.,  3000., 3500., 4500. /)  		     	
!  ancien profil
      else if ( nz1.eq.14 ) then
        alti_1D = (/ 25.,  75.,  125.,  250., 500. , 750.,  1000., &
                     1250., 1500. ,  2000., 2500.,  3000., 3500., 4500. /)  	
! nouveau profil
      else if ( nz1.eq.10 ) then
        alti_1D = (/ 25.,  75.,  150., 250., 500., 750.,  1500., 2500., 3500., 4500. /)	      
      else
        write(*,*) 'Pas de profil d altitude defini pour nz1=',nz1
      end if
            
      ! remonte ce profil de 100 m
      do k1=1,nz1
	  !alti_1D(k1) = alti_1D(k1) + 100.
	  alti_1D(k1) = alti_1D(k1) + altitude_min
      end do
      
      ! PRINT pour affichage, milieu du domaine (Grenoble)
      ixp = int(nx2/2.)
      iyp = int(ny2/2.) 
      
      end subroutine definition_profil_z
!---------------------------------------------------------------------------------------------------------------------     
      subroutine interpolation_donnees_grille

      use params
      use params_physique      
      
      implicit none 
	
      lon_g1 = lon_g2
      lat_g1 = lat_g2
            
      if (iinterp)     call interp2d(nx0,ny0,nx1,ny1,wx1,wy1,ix1,iy1,topo_g0,topo_g1)
      if (.not.iinterp)call recopy2d(nx0,ny0,nx1,ny1,        ix2,iy2,topo_g0,topo_g1)
            
      end subroutine interpolation_donnees_grille      
!---------------------------------------------------------------------------------------------------------------------     
      subroutine interpolation_var_1h
      
      use params
      use params_physique      
      
      implicit none      
                               
! variables 3D

!... altitude	
        if (it.le.2) then	 
	 if (idebug)write(*,*) '...alti_agl_g1' !,nx2, ny2, nz1
	 if (iinterp)     call interp3d(nx0,ny0,nz0,nx1,ny1,nz1,wx1,wy1,ix1,iy1,alti_agl_g0, alti_agl_g1 )
	 if (.not.iinterp)call recopy3d(nx0,ny0,nz0,nx1,ny1,nz1,        ix2,iy2,alti_agl_g0, alti_agl_g1 )
        end if 	     
	 
!... variables chimie
	do ivar = 1, nvar
	  if (icorrect_z_var(ivar)) then
	  if (idebug)write(*,*) '...'//trim(var_name(ivar))
	  if (iinterp)     call interp3d(nx0,ny0,nz0,nx1,ny1,nz1,wx1,wy1,ix1,iy1,val_g0_var(ivar,:,:,:),val_g1_var(ivar,:,:,:))
	  if (.not.iinterp)call recopy3d(nx0,ny0,nz0,nx1,ny1,nz1,        ix2,iy2,val_g0_var(ivar,:,:,:),val_g1_var(ivar,:,:,:))   
	  end if 
	end do  
	     
! variables 2D
      if (idebug)write(*,*) '...hght'
      if (iinterp)     call interp2d(nx0,ny0,nx1,ny1,wx1,wy1,ix1,iy1,hght_g0,hght_g1)
      if (.not.iinterp)call recopy2d(nx0,ny0,nx1,ny1,        ix2,iy2,hght_g0,hght_g1)      
      
      if (imeteo) then
        if (idebug)write(*,*) '...tem2'
        if (iinterp)     call interp2d(nx0,ny0,nx1,ny1,wx1,wy1,ix1,iy1,tem2_g0,tem2_g1)
        if (.not.iinterp)call recopy2d(nx0,ny0,nx1,ny1,        ix2,iy2,tem2_g0,tem2_g1)
      end if
          

      end subroutine interpolation_var_1h
!---------------------------------------------------------------------------------------------------------------------            
      subroutine calcul_temperature_2m
      
      use params
      use params_physique      
      
      implicit none      
      
      do i2=1,nx2
	do j2=1,ny2	  
	  tem2_g2(i2,j2) = tem2_g1(i2,j2) + ( topo_g1(i2,j2) - topo_g2(i2,j2) ) * gradient_adiabatique
       end do
      end do  
      
      end subroutine calcul_temperature_2m
!---------------------------------------------------------------------------------------------------------------------            
      subroutine calcul_vent_10m !version 1.1 (06/10/2010)
      
      use params
      use params_physique      
      
      implicit none   
      
      real :: fac_10m !facteur de vent altitude->10m   

      call interp2d(nx0,ny0,nx1,ny1,wx1,wy1,ix1,iy1,w10m_g0,w10m_g1)    
      
      do i2=1,nx2
	do j2=1,ny2	  
	  if (val_g1_var(met_iwins,i2,j2,1).ge.0.1) then	  
	    ! méthode au 06/10/2010: le facteur d'attenuation se calcule à partir du ratio vent raffiné/vent brut
	    fac_10m = val_g2_var(met_iwins,i2,j2,1) / val_g1_var(met_iwins,i2,j2,1)	    
	    w10m_g2(i2,j2) = w10m_g1(i2,j2) * fac_10m
	  else
	    w10m_g2(i2,j2) = w10m_g1(i2,j2)
	  end if
       end do
      end do  
      
      end subroutine calcul_vent_10m      
!---------------------------------------------------------------------------------------------------------------------     
      subroutine allocation_altitude_maille
 
      use params
      use params_physique      
      
      implicit none   
    
      iz_alti_g1 = 1
      i_alti_g1  = 0. 
       
      do ix=1, nx2
       do iy=1, ny2
           
	   ! calcul des correspondances d'altitude
	   do iz=1,nz1	
	     k1 = 1
	     if ( alti_asl_g1(ix,iy,iz).le.alti_1D(k1) ) then 
	       iz_alti_g1(ix,iy,iz) = k1
	       i_alti_g1(ix,iy,iz)  = 1.	       
	     end if 
	     
	     do k1 = 2, nz1
               if ( (alti_asl_g1(ix,iy,iz).gt.alti_1D(k1-1) ).and.(alti_asl_g1(ix,iy,iz).le.alti_1D(k1) )  ) then
	         iz_alti_g1(ix,iy,iz) = k1
	         i_alti_g1(ix,iy,iz)  = 1.		 
	       end if	
	     end do
	   end do
	   
	 end do
      end do
      
      end subroutine allocation_altitude_maille      
!---------------------------------------------------------------------------------------------------------------------     
      subroutine calcul_altitude_asl

      use params
      use params_physique      
      
      implicit none
      
      ! Si la meteo est disponible
      do i2=1,nx2
        do j2=1,ny2
          k1 = 1
          topo = topo_g1(i2,j2)
          alti_asl_g1(i2,j2,k1) = topo + alti_agl_g1(i2,j2,k1)/2	
          do k1=2,nz1
            alti_asl_g1(i2,j2,k1) = topo + ( alti_agl_g1(i2,j2,k1) + alti_agl_g1(i2,j2,k1-1) )  /2
          end do	
        end do
      end do 
             
      if (idebug) then
        write(*,'("Altitude des niveaux au point ix=",I3," iy=",I3)')ixp,iyp  
        do k1=1,nz1
	  write(*,*)  alti_asl_g1(ixp,iyp,k1)
        end do
      end if  
        
      end subroutine calcul_altitude_asl
!---------------------------------------------------------------------------------------------------------------------  
      subroutine calcul_profil_moyen
      
      use params
      use params_physique      
      
      implicit none
      
      real :: nval_var(nz1)       
      
      iprint=.false.
      if (it.eq.itp)iprint=.true.
      
      do ivar= 1, nvar  
      
      if (icorrect_z_var(ivar)) then  
      
      nval_var = 0.
      
      do i2=1,nx2
        do j2=1,ny2  	          
	  do iz=1,nz1
	      k1 = iz_alti_g1(i2,j2,iz)	 
	      profil_1D_var(ivar,k1)    = profil_1D_var(ivar,k1)    + val_g1_var(ivar,i2,j2,iz)    * i_alti_g1(i2,j2,iz)	      
	      nval_var(k1) = nval_var(k1) + i_alti_g1(i2,j2,iz) 
	   end do    
         end do 	     	  
       end do	 

        do k1=1,nz1	 
	  profil_1D_var(ivar,k1)  = profil_1D_var(ivar,k1) / nval_var(k1)	 	 
        end do    
	  
	! avec l'hypothese de profil d'altitude orthonorme
        k1 = 1	
        profil_1D_var(ivar,k1)    = profil_1D_var(ivar,2)    - ( profil_1D_var(ivar,3)    - profil_1D_var(ivar,2)   )              
        if ( profil_1D_var(ivar,k1)    .lt.0 ) profil_1D_var(ivar,k1)    = 0. 
                   
      ! affichage du profile au pas de temps 2

      if (idebug) then  
      	if (voisinage_xy_var(ivar).eq.0) then ! utilise un voisinage ? NON
        write(*,*) '#######################'        
	write(*,*) 'Profil MOYEN de '//trim(var_name(ivar))//' moyen'
        do k1=1,nz1
	  write(*,'(A10," a ",F6.0," m =",F8.2,A10)')  var_name(ivar),alti_1D(k1),profil_1D_var(ivar,k1),trim(var_unit(ivar))
        end do
	end if
      end if  
      
      end if ! correct_z ?
      
      end do ! ivar ? 

      end subroutine calcul_profil_moyen
!---------------------------------------------------------------------------------------------------------------------     
      subroutine calcul_profil_voisinage
      
      use params
      use params_physique      
      
      implicit none
      
      ! local
      real,parameter  :: pourcentage_limite = 0.05 ! pourcentage limite pour calcul sur les niveaux
      integer :: ixmin_loc,ixmax_loc
      integer :: iymin_loc,iymax_loc
      real :: nval_var(nz1)
      real :: pourcentage_dispo 
      real :: a, b
            
! fonction externe      
      integer :: minint  

      ! remet les valeurs a 0 pour le point (ix, iy)           
      do ivar= 1, nvar 
                     
        ! variable corrige en Z ?
        if (.not.icorrect_z_var(ivar)) goto 89
	
	if (idebug)write(*,*) '...'//trim(var_name(ivar))
      
        do ix=1, nx2
          do iy=1, ny2
	  
	  ! initialisation variable n point ix, iy
	  nval_var   = 0. 
    
          ! affichage TRUE/FALSE
          iprint=.false.
          if ((ix.eq.ixp).and.(iy.eq.iyp))iprint=.true.
    
	 ! OPTIMISATION : si le point n'est pas dans la region, le voisinage n'est pas calculé
	  if ( (region1km).and.( .not.iregion(ix,iy) ) ) then
	    profil_3D_var(:,ix,iy,:) = 0.0
	    goto 88
	  end if
  
	  !!!!! pour polluants sur voisinage    
	  if (voisinage_xy_var(ivar).ne.0) then ! utilise un voisinage ? OUI
	  
            ixmin_loc = max   (   1, ix-voisinage_xy_var(ivar) )
	    ixmax_loc = minint( nx2, ix+voisinage_xy_var(ivar) )
            iymin_loc = max   (   1, iy-voisinage_xy_var(ivar) )
	    iymax_loc = minint( ny2, iy+voisinage_xy_var(ivar) )	  
	    
	    ! balaie toute la grille -> cette boucle prend beaucoup de temps
	    do iz = 1, nz1 ! niveaux verticaux de l'axe d'altitude
	      if (iz.eq.1.and.iprint.and.idebug)&
	         write(*,'("var ",I2," ixp=",I3," iyp=",I3," voisinage ix1..ix2=",2I4," iy1..iy2=",2I4)' ) ivar,ix,iy,ixmin_loc,ixmax_loc, iymin_loc,iymax_loc
	      do i2 = ixmin_loc,ixmax_loc
                do j2 = iymin_loc,iymax_loc
	          k1 = iz_alti_g1(i2,j2,iz)
	          profil_3D_var(ivar,ix,iy,k1) = profil_3D_var(ivar,ix,iy,k1)  + val_g1_var(ivar,i2,j2,iz) * i_alti_g1(i2,j2,iz)
	          nval_var(k1) = nval_var(k1) + i_alti_g1(i2,j2,iz)
                end do     
	      end do	  
	    end do !!!! iz=1,nz1	
	  	  
            ! moyennes par niveaux (le nombre de données doit etre superieur à 5 % des données potentielles)
	    ! pourcentage_limite = 0.05	  	    
            ! sur les autres niveau, applique le profil par voisinage ou profil moyen   
	    
	    do k1 = nz1-1 ,1, -1
	    
	      pourcentage_dispo = nval_var(k1) / voisinage_xy_var(ivar)**2
	      
	      if ( pourcentage_dispo .gt. pourcentage_limite ) then	    
	        profil_3D_var(ivar,ix,iy,k1) = profil_3D_var(ivar,ix,iy,k1) / nval_var(k1)
	      else
	        if (iprint.and.idebug) then 
		  !write(*,*) 'ATTENTION valeurs '//trim(var_name(ivar))//' nulles/insuffisantes au niveau',int(k1)
	          write(*,'("ATTENTION : niveau ",I2," valeurs=",F8.2,"% limite=",F8.2,"%")') k1,pourcentage_dispo*100., pourcentage_limite*100.
		end if
		
		! METHODE 1 : concentrations constantes
	        profil_3D_var(ivar,ix,iy,k1) =  profil_3D_var(ivar,ix,iy,k1+1)
		
		! METHODE 2 : concentration interpolée par régression linéaire
		if (k1.lt.nz1-1) then
		  ! coefficient directeur de la droite
		  a = (profil_3D_var(ivar,ix,iy,k1+2) - profil_3D_var(ivar,ix,iy,k1+1)) / ( alti_1D(k1+2) - alti_1D(k1+1) ) 
		  profil_3D_var(ivar,ix,iy,k1) = profil_3D_var(ivar,ix,iy,k1+1) - a * ( alti_1D(k1+1) - alti_1D(k1) )
		  if ( profil_3D_var(ivar,ix,iy,k1) .lt.0.0 ) profil_3D_var(ivar,ix,iy,k1) = 0.0		
		end if
		
	      end if
	    end do !!!! k1=1,nz1 
	    	    
	  else ! utilise un voisinage ? NON
	  
	    if (iprint.and.idebug) write(*,*) 'ATTENTION : utilise profil '//trim(var_name(ivar))//' moyen'
	    profil_3D_var(ivar,ix,iy,:) = profil_1D_var(ivar,:)
	    
	  end if  
	          
	  ! au sommet, utilise le profil moyen
	  k1=nz1
	  profil_3D_var(ivar,ix,iy,k1) = profil_1D_var(ivar,k1)
	        
          ! avec l'hypothese de profil d'altitude orthonorme
          !k1=1
	  !profil_3D_var(ivar,ix,iy,k1) = profil_3D_var(ivar,ix,iy,2)   -  ( profil_3D_var(ivar,ix,iy,3)    - profil_3D_var(ivar,ix,iy,2)    )		 
	  !if ( profil_3D_var(ivar,ix,iy,k1) .lt.0.0 ) profil_3D_var(ivar,ix,iy,k1) = 0.0
          	  	  
88       continue ! pas de calcul de profils pour les point hors de la region  

       end do !!! iy
      end do !!! ix
      
      
      ! affichage TRUE/FALSE
      iprint=.false.
      if (it.eq.itp)iprint=.true.	
      if (idebug) then
	 if (voisinage_xy_var(ivar).ne.0) then ! utilise un voisinage ? OUI        
         ! pour affichage, milieu du domaine (Grenoble)
	 write(*,*) '##################################'
	 write(*,*) 'Profil LOCAL de '//trim(var_name(ivar))//' (ix,iy)=',ixp,iyp
         do k1=1,nz1
	   write(*,'(A6," a ",F6.0," m =",F7.1,A10)')  var_name(ivar),alti_1D(k1),profil_3D_var(ivar,ixp,iyp,k1),trim(var_unit(ivar))
         end do
	 end if
      end if  

89    continue

      end do ! ivar
              
      end subroutine calcul_profil_voisinage    
!---------------------------------------------------------------------------------------------------------------------     
      subroutine calcul_var_z
      
      use params
      use params_physique      
      
      implicit none
      
      real :: topo1, topo2
      
      ! Initialisation
      ! bien specifie que si la recherche echoue, le gradient est nul   
      dvar_dz_loc =0.

      do ivar= 1, nvar
      
      ! Corrige-t-on la varaible en profil Z
      if (icorrect_z_var(ivar)) then            
   
      do ix=1,nx2
         do iy=1,ny2
	   
	   ! reinitialisation
	   topo1 = topo_g1(ix,iy)	   
	   topo2 = topo_g2(ix,iy)	   	   	    
	   	   
	   do k1 = 2 , nz1-1	      
	      if ( ( ((topo1+topo2)/2.).lt.alti_1D(k1)).and.(  ((topo1+topo2)/2.).ge.alti_1D(k1-1))) then
		dvar_dz_loc(ivar,ix,iy)   = (profil_3D_var(ivar,ix,iy,k1+1)    - profil_3D_var(ivar,ix,iy,k1-1)   ) / (alti_1D(k1+1) - alti_1D(k1-1) )				    
	      end if  
	   end do
	   	   
	   !if (.not.correct_z_var(ivar))  dvar_dz_loc(ivar,i2,j2)   = 0.
	   	   
	   ! calcul de l'O3 sur la grille fine pour les altitudes asl entre 300 et 3000 m
	   if ( topo2.ge.5000.) then
	     dvar_dz_loc(ivar,ix,iy)   = 0.	     
	   end if

	   dvar_dz_loc(ivar,ix,iy) = dvar_dz_loc(ivar,ix,iy)   * ( topo2 - topo1 )
	      
	   ! grille RAFINEE normale		       
	   !val_g2_var(ivar,ix,iy,1)  = val_g1_var(ivar,ix,iy,1)    !+ dvar_dz_loc(ivar,ix,iy)
	   	   			     
	   ! corrige egalement la grille FINE
	   !val_g1_var(ivar,ix,iy,1)  = val_g1_var(ivar,ix,iy,1)    !+ dvar_dz_loc(ivar,ix,iy)
	   	       	            	          
           ! Verifie que les concentrations restent positives	   	   
	   !if ( val_g2_var(ivar,ix,iy, 1) .lt. 0. ) val_g2_var(ivar,ix,iy, 1) = 0.
	   !if ( val_g1_var(ivar,ix,iy, 1) .lt. 0. ) val_g1_var(ivar,ix,iy, 1) = 0.
	   	   	   	   
	 end do
        end do
        if (idebug)write(*,*) '... recalcul '//trim(var_name(ivar))//' avec profils Z OK'
	
      end if ! correct_z_var
      end do ! nvar  
      
      end subroutine calcul_var_z 
!---------------------------------------------------------------------------------------------------------------------
      subroutine calcul_vitesse_vent_3d(nx,ny,nz,uu,vv,ff)
      
      implicit none 
      !entrees
      integer :: nx, ny, nz  
      real :: uu(nx,ny,nz), vv(nx,ny,nz), ff(nx,ny,nz)
      !local
      integer :: i, j, k         
     
      do i= 1, nx
      do j= 1, ny
      do k= 1, nz            
        ff(i,j,k) = ( uu(i,j,k)**2 + vv(i,j,k)**2 )**.5
      end do	
      end do
      end do                  
      end subroutine
