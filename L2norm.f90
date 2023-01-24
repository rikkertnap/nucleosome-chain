
module vectornorm

    use precision_definition
    implicit none

contains

    ! computes L2norm of vector f 

    function l2norm(f,n)result(norm)
      
        implicit none

        integer, intent(in)  :: n 
        real(dp), intent(in) :: f(n)
        real(dp)             :: norm ! output

        integer :: i ! dummy index

        norm=0.0_dp
        do i=1,n
            norm = norm + f(i)**2
        enddo
        norm=sqrt(norm)
        
    end function l2norm

    ! computes L2norm of vector f

    function l2norm_f90(f)result(norm)
      
        implicit none

        real(dp), intent(in) :: f(:)
        real(dp)             :: norm ! output
        integer              :: i ,n

        n=size(f)
        norm=0.0_dp
        do i=1,n
            norm = norm + f(i)**2
        enddo
        norm=sqrt(norm)
        
    end function l2norm_f90

    ! computes L2norm of a sub vector f, range nb,ne

    function l2norm_sub(f,nb,ne)result(norm)
      
        implicit none

        real(dp), intent(in) :: f(:)
        integer, intent(in)  :: nb, ne
        real(dp)             :: norm ! output
        integer              :: i ,n

        norm=0.0_dp
        do i=nb,ne
            norm = norm + f(i)**2
        enddo
        norm=sqrt(norm)
        
    end function l2norm_sub

end module vectornorm
      
