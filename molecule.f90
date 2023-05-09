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

contains

    subroutine init_zero_moleclist(moleclist_record)
    
        type(moleclist), intent(inout) :: moleclist_record

        moleclist_record%sol = 0.0_dp
        moleclist_record%Na = 0.0_dp
        moleclist_record%Cl = 0.0_dp
        moleclist_record%K = 0.0_dp
        moleclist_record%Rb = 0.0_dp
        moleclist_record%Ca = 0.0_dp
        moleclist_record%Mg= 0.0_dp
        moleclist_record%NaCl = 0.0_dp
        moleclist_record%KCl = 0.0_dp
        moleclist_record%Hplus = 0.0_dp
        moleclist_record%OHmin = 0.0_dp


    end subroutine init_zero_moleclist

end module molecules

