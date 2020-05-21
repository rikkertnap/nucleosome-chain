                                                                  
!     .. module file of global variables 

module globals

    use precision_definition
    use mathconst
  
    implicit none
  
    !     .. variables



    integer  :: nsize         ! size lattice, numer of layers
    integer  :: nseg          ! length of polymer 
    integer  :: nsegtypes     ! number of segment types 
    integer  :: cuantas       ! number of polymer configurations
    integer  :: max_confor    ! maximum number of polymer configurations
   
    integer(8)  :: neq        ! number of non-linear equations
    integer(8)  :: neqmax     ! maximum number of non-linear equations
    integer  :: neqint        ! number of non-linear equations, for mpi fnc bindings
    
    character(len=15) :: systype   ! systype selects fcn    
    character(len=15) :: runtype   ! runtype
    character(len=2)  :: bcflag(2) ! bcflag selects bc surface 

    integer, parameter :: LEFT = 1
    integer, parameter :: RIGHT = 2

    logical, parameter :: DEBUG = .TRUE. ! switch for debug information

end module globals

