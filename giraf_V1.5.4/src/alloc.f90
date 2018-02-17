    subroutine allocation_grilles

    use params
    
    implicit none

    ! Grille CHIMERE brutes
    allocate( lon_gc    (nxc,nyc) )    
    allocate( lat_gc    (nxc,nyc) ) 
    allocate( landuse_gc(nxc,nyc) )    
    allocate( topo_gc   (nxc,nyc) )    

    ! Grille FINE de sortie
    allocate( lon_g2    (nx2,ny2) )    
    allocate( lat_g2    (nx2,ny2) ) 
    allocate( landuse_g2(nx2,ny2) )    
    allocate( topo_g2   (nx2,ny2) ) 
    
    ! initialise indices de depart de la grille FINE dans la grille COARSE
    dx1min = 1	   
    dy1min = 1 
               
    end subroutine allocation_grilles    
!------------------------------------------------------    
    subroutine allocation_var_3d
    
    use params
    
    implicit none

    ! Grille CHIMERE coarse
            
      allocate( lon_g0    (nx0,ny0) )    
      allocate( lat_g0    (nx0,ny0) ) 
      allocate( landuse_g0(nx0,ny0) )    
      allocate( topo_g0   (nx0,ny0) )    
    
      ! out.*.nc        
      allocate(val_g0_var (nvar,nx0,ny0,nz0))

      ! met.*.nc       
      allocate(hght_g0    (nx0,ny0    ))
      allocate(tem2_g0    (nx0,ny0    ))
      allocate(w10m_g0    (nx0,ny0    ))      
      allocate(alti_agl_g0(nx0,ny0,nz0))      

      ! Grille CHIMERE fine
      ! ground
      allocate( lon_g1    (nx1,ny1) )    
      allocate( lat_g1    (nx1,ny1) ) 
      allocate( landuse_g1(nx1,ny1) )    
      allocate( topo_g1   (nx1,ny1) ) 

      ! out.*.nc        
      allocate(val_g1_var(nvar,nx1,ny1,nz1))          
      
      ! met.*.nc       
      allocate(hght_g1    (nx1,ny1    ))
      allocate(tem2_g1    (nx1,ny1    ))
      allocate(w10m_g1    (nx1,ny1    ))               
      allocate(alti_asl_g1(nx1,ny1,nz1)) 
      allocate(alti_agl_g1(nx1,ny1,nz1))
      allocate(iz_alti_g1 (nx1,ny1,nz1))   
      allocate(i_alti_g1  (nx1,ny1,nz1))              
      
! Grille fine  
      ! met.*.nc       
      allocate(hght_g2 (nx2,ny2))
      allocate(tem2_g2 (nx2,ny2)) 
      allocate(w10m_g2 (nx2,ny2))        
    
      ! out.*.nc            
      allocate(val_g2_var  (nvar,nx2,ny2,nz2))
      
      ! FDMS
      allocate(pm10_nonvolat_g2(nx2,ny2,nz2))
      allocate(ecart_fdms_g2   (nx2,ny2,nz2))
      
! variables speciales PROFILS sur grille FINE
      ! profils 1D
      allocate(profil_1D_var(nvar,nz1)) 
      ! profils altitude 3D
      allocate(profil_3D_var  (nvar,nx1,ny1,nz1))
               
      ! profils altitude 1D
      allocate(alti_1D(nz1))                   

      ! gradient local dVAR/dZ
      allocate(dvar_dz_loc     (nvar,nx1,ny1)) 
      ! gradient local dVAR/demis      
      allocate(dvar_de_loc     (nvar,nx1,ny1)) 
      
      ! initialisations des gradients
      dvar_de_loc = 0.
      dvar_dz_loc = 0.

      ! emissions
      allocate(emis_g0       (nvar_emis,nx0,ny0))  
      allocate(ratio_emis_g0 (nvar_emis,nx0,ny0))
      allocate(nratio_emis_g0(nvar_emis,nx0,ny0))         
      
      allocate(emis_g2       (nvar_emis,nx2,ny2)) 
      allocate(ratio_emis_g2 (nvar_emis,nx2,ny2)) 
            
      ! interpolation
      allocate(ix2(nx2,ny2))    
      allocate(iy2(nx2,ny2)) 
      allocate(ix1(nx2,ny2))
      allocate(iy1(nx2,ny2)) 
      allocate(wx1(nx2,ny2))
      allocate(wy1(nx2,ny2))                
    
    end subroutine allocation_var_3d
