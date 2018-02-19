!*****************************************************************************      
subroutine check(status)
  ! Internal subroutine - checks error status after each netcdf, prints out text message each time
  !   an error code is returned. 

  use netcdf
  use typesizes

  implicit none

  integer, intent ( in) :: status
   
  if(status /= nf90_noerr) then 
    print*, 'Error occured'
    print *, 'status=',status
    print*, trim(nf90_strerror(status))
    stop    
  end if
end subroutine check  

