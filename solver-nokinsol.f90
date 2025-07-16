subroutine solver(x, xguess, accuracy, residual, isSolution)

    use globals
    use parameters
    use listfcn
    use fcnpointer
    use anderson

    implicit none

    !     .. arguments
    real(dp), dimension(neq) :: x  ! explicit size array
    !real(dp) :: x(neq)
    real(dp), dimension(neq) :: xguess 
    real(dp), intent(in) :: accuracy
    real(dp), intent(out) :: residual
    logical, intent(out)  :: isSolution

    call set_size_neq
    call set_fcn

    if(method.eq."kinsol") then

!        call kinsol_gmres_solver(x, xguess, accuracy, residual, issolution)

    else if(method.eq."anderson") then

        call anderson_min_loop(xguess, x, accuracy, residual, isSolution, maxfkfunevals, neq)

    else if(method.eq."simple") then
    
        call simple_min_loop(xguess, x, accuracy, residual, isSolution, maxfkfunevals, neq)


    else 
        print*,"Solver method incorrect"
        stop
    endif

end subroutine solver
