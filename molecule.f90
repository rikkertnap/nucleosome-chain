module molecules

    use precision_definition
    implicit none

    type moleclist
        real(dp) :: sol
        real(dp) :: Na
        real(dp) :: Cl
        real(dp) :: K
        real(dp) :: Rb
        real(dp) :: Ca
        real(dp) :: Mg
        real(dp) :: NaCl
        real(dp) :: KCl
        real(dp) :: Hplus
        real(dp) :: OHmin
    end type moleclist


    type bornmoleclist
        real(dp) :: pol
        real(dp) :: polCa
        real(dp) :: polMg
        real(dp) :: Na
        real(dp) :: Cl
        real(dp) :: K
        real(dp) :: Rb 
        real(dp) :: Ca
        real(dp) :: Mg
        real(dp) :: Hplus
        real(dp) :: OHmin
    end type bornmoleclist 

  
end module molecules

