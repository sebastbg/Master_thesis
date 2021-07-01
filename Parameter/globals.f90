module globals
    
    implicit none
    public
    
    real(8), parameter      ::  pi = 4.d0*atan(1.d0)
    integer                 ::  Nslips                           
    real(8)                 ::  factor								
    real(8), allocatable    ::  P(:,:), Omega(:,:) 

end module globals
    
