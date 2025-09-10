subroutine solver(x, xguess, accuracy, residual, isSolution)

    use precision_definition
    use globals ,only : neq
    use parameters, only : method, set_size_neq, maxfkfunevals
    use listfcn, only :  set_fcn
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

        print*,"Warning kinsol cannot be called."
        stop

    else if(method.eq."anderson") then

        call anderson_min_loop(xguess, x, accuracy, residual, isSolution, maxfkfunevals, neq)

    else if(method.eq."simple") then
    
        call simple_min_loop(xguess, x, accuracy, residual, isSolution, maxfkfunevals, neq)


    else 
        print*,"Solver method incorrect."
        stop
    endif

end subroutine solver
