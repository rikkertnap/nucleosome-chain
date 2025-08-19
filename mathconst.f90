!     module of mathconst

module mathconst

    use precision_definition  
    implicit none
  
    real(dp) :: pi
  
    contains 
    
        subroutine make_mathconst()

            pi = acos(-1.0_dp)

        end subroutine make_mathconst

end module mathconst
