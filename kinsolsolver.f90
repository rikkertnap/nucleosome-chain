! -----------------------------------------------------------------------------| 
!     kinsolver.f90                                                            |  
!     module containing functions and subroutines for the usage of             |  
!     the numerical kinsol library routine                                     |
!------------------------------------------------------------------------------| 


module kinsolvars

    use precision_definition
    implicit none
    !     .. variables

    real(dp), dimension(:), allocatable ::   pp ! pre-processor variable

end module kinsolvars

!     The routine fkpset is the preconditioner setup routine. It must have
!     that specific name to be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)

    use globals, only : neq
    use kinsolvars
    implicit none

    integer :: ier
    integer :: i
    double precision :: udata(*), uscale(*), fdata(*), fscale(*)
    double precision :: vtemp1(*), vtemp2(*)

    do i = 1, neq
        pp(i) = 0.5_dp / (udata(i) + 5.0_dp)
    enddo
    ier = 0

end subroutine fkpset

!     The routine fkpsol is the preconditioner solve routine. It must have
!     that specific name to be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
  
    use globals, only : neq
    use kinsolvars
    implicit none
  
    integer :: ier
    integer :: i
    double precision :: udata(*), uscale(*), fdata(*), fscale(*)
    double precision :: vv(*), ftem(*)
  
    do  i = 1, neq
        vv(i) = vv(i) * pp(i)
    enddo

    ier = 0
  
end subroutine fkpsol



!     Initialiazaition and run of the kinsolver                                 
!     Scale Preconditioned GMRES solver                                        
!     pre:  input vector x and guess of x (xguess) and n (number of equations) 
!            and fnorm                                                         
!     post: solution of SCMFT stored in x and  fnorm residual error            


subroutine kinsol_gmres_solver(x, xguess, n, error, fnorm)
  
    use globals, only : nsize, neq 
    use kinsolvars
    use parameters, only : iter
    use myutils

    implicit none

    interface 
        subroutine set_contraints(constr)
            use precision_definition
            implicit none
            real(dp) :: constr(:)
        end subroutine
    end interface    


    ! .. neq iout(15) and msbre match C type long int.
    ! .. variables and constant declaractions 

    ! .. array arguments      
    real(dp) :: x(neq)
    real(dp) :: xguess(neq)

    !  .. scalar arguments
    real(dp) :: fnorm 
    real(dp) :: error
    integer(8) :: n                 !  number of equations

    !  .. local arguments 

    integer(8) :: iout(15)         ! Kinsol additional output information
    integer(8) :: msbpre           ! maximum number of iterations without prec. setup 

    real(dp) :: rout(2)           ! Kinsol additional out information
    integer  :: i                 ! dummy index 
    integer  :: ier               ! Kinsol error flag
    integer  ::  maxniter
    real(dp) :: fnormtol, scsteptol
    real(dp) :: fscale(neq)
    real(dp) :: constr(neq)
    integer  ::  globalstrat, maxl, maxlrst
    character(len=lenText) :: text, rstr, istr
  
    !     .. executable statements 

    !     .. init of kinsol variables 
                            
    msbpre  = 5               ! maximum number of iterations without prec. setup 
    fnormtol = error          ! Function-norm stopping tolerance
    scsteptol = error         ! Function-norm stopping tolerance
    maxl = neq                ! maximum Krylov subspace dimension 
    maxlrst = 20              ! maximum number of restarts
    globalstrat = 0           ! inexact Newton  
    maxniter =1000            ! maximum of nonlinear iterations default 200  


    allocate(pp(neq))


    do i = 1, neq             
        constr(i) = 1.0_dp      ! constraint vector  
        fscale(i) = 1.0_dp      ! scaling vector  
        x(i) = xguess(i)        ! initial guess
    enddo
  
    !do i =neq/4+1, neq/2
    !    constr(i)=0.0_dp
    !enddo   
    
    call set_contraints(constr)
  
    call fnvinits(3, neq, ier) ! inits NVECTOR module
    
    if (ier .ne. 0) then      ! 3 for Kinsol, neq number of equations, ier error flag (must be 0 on output)
        print*, 'SUNDIALS_ERROR: FNVINITS returned IER = ', ier
        stop
    endif
  
    call fkinmalloc(iout, rout, ier) ! Allocates memory and output additional information
  
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
        stop
    endif
  
    ! Additional input information

    call fkinsetiin('MAX_SETUPS', msbpre, ier) 

    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSETIIN returned IER = ', ier
        call fkinfree          ! free memory
        stop
    endif

    call fkinsetiin('MAX_NITERS', maxniter, ier)

    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSETIIN returned IER = ', ier
        call fkinfree          ! free memory
        stop
    endif

    call fkinsetrin('FNORM_TOL', fnormtol, ier)
    call fkinsetrin('SSTEP_TOL', scsteptol, ier)
    call fkinsetvin('CONSTR_VEC', constr, ier) 

    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSETVIN returned IER = ', ier
        call fkinfree          ! free memory
        stop
    endif
    !     .. init Scale Preconditioned GMRES solver 
  
    call fkinspgmr(maxl, maxlrst, ier)   
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
        call fkinfree          ! free memory
        stop
    endif
  
    !     .. preconditioner  zero no preconditioner  
    call fkinspilssetprec(0, ier) 

    !     .. call solver
    call fkinsol(x, globalstrat, fscale, fscale, ier) 
    fnorm=rout(1) 
  
    if (ier .lt. 0) then
        print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
        print*, 'Linear Solver returned IER = ', iout(9)
        call fkinfree
        stop
    endif
    
    write(rstr,'(E25.16)')fnorm
    text="Found solution: fnorm = "//trim(rstr)//" "
    call print_to_log(LogUnit,text)

    write(istr,'(I8)')iter
    text="number of iterations  = "//trim(istr)
    call print_to_log(LogUnit,text)
    
    call fkinfree             ! free memory
    deallocate(pp)
    return
  
end subroutine kinsol_gmres_solver

! set constrains on vector x  depending on sysflag value
! 
subroutine set_contraints(constr)
    
    use precision_definition
    use globals, only : sysflag, neq

    implicit none
        
    real(dp):: constr(:)

    integer :: i
 
    if(sysflag=="fcnelectdouble") then
        do i =neq/4+1, neq/2
            constr(i)=0.0_dp
        enddo   
    endif    

end 



!     .. wrapper function 

subroutine fkfun(x,f,ier)
  
    use precision_definition
    use globals,  only  : neq
    use fcnpointer

    implicit none

    real(dp), dimension(neq) :: x  ! expliciet size array   
    real(dp), dimension(neq) :: f
    integer ::  ier

    call fcnptr(x,f,neq)

    ier=0  

end subroutine fkfun























