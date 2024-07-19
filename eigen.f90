module eigenvalues

    use precision_definition   
    implicit none

    real(dp), parameter :: eps_eigen=1.0e-8_dp

    private :: eps_eigen

contains


! Computes usign lapack double precision routine dsyev the eigenvalues of a real symmetric matrix 
! input integer :: k : dimension of matrix
!       real(dp) :: A(k,k) : real symmetri matrix
! output real(dp) eigenvalues(k) :: ordered accending list of eigenvalues        
!        integer :: info : retrun value dsyev 

subroutine comp_eigenvalues(k,A,eigenvalues,info)

    integer, intent(in)     :: k
    real(dp), intent(in)    :: A(k,k)
    real(dp), intent(inout) :: eigenvalues(k)
    integer, intent(inout)  :: info

    !local 
    real(dp),allocatable :: work(:)
    integer              :: lwork
    real(dp)    :: B(k,k)

    B=A ! prevent overwritting of input matrix A 

    lwork = max(1,3*k-1)
    allocate(work(lwork))
    
    call dsyev('N','U',k,B,k,eigenvalues,WORK,LWORK,info) 

end subroutine

function Asphericty_parameter(Rgsqr,gyr_tensor)result(Ap)
    
    real(dp), intent(in) :: Rgsqr
    real(dp), intent(in) :: gyr_tensor(3,3)
    real(dp) :: Ap

    ! local 
    real(dp):: eigenvalues(3)
    integer :: info_eigen
    real(dp) :: bsph,ccyl
    real(dp) :: sumeigen,diffRg2

    call comp_eigenvalues(3,gyr_tensor,eigenvalues,info_eigen)
    
    bsph = eigenvalues(3)-(eigenvalues(1)+eigenvalues(2))/2.0_dp
    ccyl = eigenvalues(2)-eigenvalues(1)

    Ap=(bsph*bsph+(3.0_dp/4.0_dp)*ccyl*ccyl)/Rgsqr**2

    sumeigen = sum(eigenvalues)
    diffRg2=abs(sumeigen-Rgsqr)

    if(diffRg2>eps_eigen) then
        print*,"Difference trace gryation tensor and Rg squared larger tolerance:"
        print*,"eigenvalues = ",eigenvalues
        print*,"sum eigenvalues = ",sumeigen," Rg2 = ",Rgsqr, " diff=",diffRg2
    endif 

    if(info_eigen .ne. 0) then
        print*, "Warning in Asphericty_parameter : info = ", info_eigen
    endif

    !print*,"eigenvalues = ",eigenvalues
    !print*,"sum eigenvalues = ",sumeigen," Rg2 = ",Rgsqr, " diff=",diffRg2
    !print*,"Gyration tensor:"
    !call print_matrix(gyr_tensor,3)

end function

function eigenvalue_of_avgyr_tensor(avgyr_tensor)result(eigenvalues)

    real(dp), intent(in) :: avgyr_tensor(3,3) 
    real(dp):: eigenvalues(3)

    integer :: info_eigen

    call comp_eigenvalues(3,avgyr_tensor,eigenvalues,info_eigen)

end function eigenvalue_of_avgyr_tensor

subroutine print_matrix(mat,k)

    integer, intent(in) :: k 
    real(dp), intent(in) :: mat(k,k)

    integer :: i, j

    do i = 1,k
        print *,(mat(i,j),j=1,k)
    enddo

end subroutine 

end module eigenvalues

