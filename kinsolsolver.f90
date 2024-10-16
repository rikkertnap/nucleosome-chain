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


subroutine kinsol_gmres_solver(x, xguess, error, fnorm,isSolution)
  
    !use mpivars
    !use ieee_arithmetic, only : ieee_is_nan    ! alternative for function isNaN in myutils
    use globals, only : nsize, neq, systype
    use kinsolvars
    use parameters, only : iter, precondition, maxniter, maxfkfunevals
    use myutils
    use listfcn, only : set_contraints

    implicit none


    ! .. neq iout(15) and msbre match C type long int.
    ! .. variables and constant declaractions 

    ! .. array arguments      
    real(dp) :: x(neq)
    real(dp) :: xguess(neq)

    !  .. scalar arguments
    real(dp), intent(out) :: fnorm 
    real(dp), intent(in)  :: error
    logical,  intent(out) :: isSolution
   
    !  .. local arguments 

    integer(8) :: iout(15)         ! Kinsol additional output information
    integer(8) :: msbpre           ! maximum number of iterations without prec. setup 

    real(dp) :: rout(2)           ! Kinsol additional out information
    integer  :: i                 ! dummy index 
    integer  :: ier               ! Kinsol error flag
    real(dp) :: fnormtol, scsteptol
    real(dp) :: fscale(neq)
    real(dp) :: constr(neq)
    integer  ::  globalstrat, maxl, maxlrst
    character(len=lenText) :: text, rstr, istr
    integer(8) :: nfe,nni 
  
    ! .. executable statements 

    ! .. init of kinsol variables 
                            
    msbpre  = 5               ! maximum number of iterations without prec. setup 
    fnormtol = error          ! Function-norm stopping tolerance
    scsteptol = error         ! Function-norm stopping tolerance
    maxl = 1000               ! maximum Krylov subspace dimension 
    maxlrst = 20              ! maximum number of restarts
    globalstrat = 0           ! inexact Newton  
    !maxniter =1000            ! maximum of nonlinear iterations default 200 , see read_input_file

    allocate(pp(neq))

    do i = 1, neq             
        constr(i) = 1.0_dp      ! constraint vector  
        fscale(i) = 1.0_dp      ! scaling vector  
        x(i) = xguess(i)        ! initial guess
    enddo
  
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
    ! .. init Scale Preconditioned GMRES solver 
  
    call fkinspgmr(maxl, maxlrst, ier)   
    if (ier .ne. 0) then
        print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
        call fkinfree          ! free memory
        stop
    endif
  
    if(precondition) then
        call fkinspilssetprec(1, ier) 
    else     
        !     .. preconditioner  zero no preconditioner  
        call fkinspilssetprec(0, ier) 
    endif    
    ! .. call solver
    call fkinsol(x, globalstrat, fscale, fscale, ier) 
    
    fnorm=rout(1)
    nfe=iout(4) 
    nni=iout(3)
  
    ! .. determine quality of solution
    
    !isSolution=(ier==0).and.(.not.ieee_is_nan(fnorm))
    
    isSolution=(ier==0).and.(.not.isNaN(fnorm))

    if(isSolution) then  
    
        write(rstr,'(E25.16)')fnorm
        text="Found solution: fnorm = "//trim(rstr)
        call print_to_log(LogUnit,text)
        write(istr,'(I8)')iter
        text="number of iterations  = "//trim(istr)
        call print_to_log(LogUnit,text)
        write(istr,'(I8)')ier
        text="kinsol return value  = "//trim(istr)
        call print_to_log(LogUnit,text)
        write(istr,'(I8)')nfe
        text="nfe  = "//trim(istr)
        call print_to_log(LogUnit,text)
        write(istr,'(I8)')nni
        text="nni  = "//trim(istr)
        call print_to_log(LogUnit,text)
    
    else
        write(istr,'(I5)')ier
        text='SUNDIALS_ERROR: FKINSOL returned IER = '//trim(istr)
        call print_to_log(LogUnit,text) 

        if(ier==1) then 

            write(rstr,'(E25.16)')fnorm
            text='Input allready a solution: fnorm = '//trim(rstr)
            call print_to_log(LogUnit,text)  
            write(istr,'(I8)')iter
            text="number of iterations  = "//trim(istr)
            call print_to_log(LogUnit,text)

            isSolution=.true.             ! overrule and accept
        
        elseif(ier==2) then
        
            write(rstr,'(E25.16)')fnorm
            text='Kinsol stalling : fnorm = '//trim(rstr)
            call print_to_log(LogUnit,text) 
            write(istr,'(I8)')iter
            text="number of iterations  = "//trim(istr)
            call print_to_log(LogUnit,text)
        
        elseif(iter>maxfkfunevals) then

            write(rstr,'(E25.16)')fnorm
            text="Maximum of fkfun evals exceeded: fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)
            write(istr,'(I8)')iter
            text="number of iterations  = "//trim(istr)
            call print_to_log(LogUnit,text)
    
        else
    
            write(rstr,'(E25.16)')fnorm
            text="No solution: fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)
            write(istr,'(I8)')iter
            text="number of iterations  = "//trim(istr)
            call print_to_log(LogUnit,text)    

        endif
        
    endif    
    
    call fkinfree             ! free memory
    deallocate(pp)
  
end subroutine kinsol_gmres_solver


!  .. wrapper function 

!  .. since fkun can not be placed in module, because of kinsol
!  .. x and f need to be explicit size array's 


subroutine fkfun(x,f,ier)
  
    use precision_definition
    use globals,  only  : neq
    use fcnpointer
    use parameters, only : iter,maxfkfunevals

    implicit none

    real(dp), dimension(neq) :: x  ! explicit size array   
    real(dp), dimension(neq) :: f
    integer ::  ier

    call fcnptr(x,f,neq)

    if(iter<=maxfkfunevals) then 
        ier=0
    else
        ier=-1 ! stop 
    endif
    
end subroutine fkfun























