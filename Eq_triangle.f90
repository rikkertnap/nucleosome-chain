module Eqtriangle
    
    use mathconst
    use precision_definition  
    
    implicit none

contains    


    ! Equilateralness =  ratio beteen triangle and equilateral triangle with same perimeter

    function equilateralness(a,b,c)result(equilat)
        
        real(dp), intent(in) :: a, b, c
        real(dp) :: equilat 

        real(dp) :: p,s,ee

        p = a+b+c
        s = p/2.0_dp
        ee = p/3.0_dp

        equilat=sqrt( (s-a)*(s-b)*(s-c)/(s-ee)**3)

    end function equilateralness

    function distance(p,q) result(dist)
        
        real(dp), intent(in) :: p(3), q(3)
        real(dp) :: dist 
        real(dp) :: deltapq(3)

        deltapq=p-q
        
        dist = sqrt(sum( deltapq** 2 ))

    end function distance
        
    function equilateralness_vectors(p,q,r) result(equilat)

        real(dp), intent(in) :: p(3), q(3), r(3)

        real(dp) :: equilat
        real(dp) :: a, b, c

        a = distance(p,q)
        b = distance(p,r)
        c = distance(q,r)

        equilat =equilateralness(a,b,c)
       
    end function equilateralness_vectors   


end module