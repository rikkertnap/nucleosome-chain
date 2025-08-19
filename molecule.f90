module molecules

    use precision_definition
    implicit none

    type moleclist
        real(dp) :: sol
        real(dp) :: Na
        real(dp) :: Cl
        real(dp) :: K
        real(dp) :: Fe2
        real(dp) :: Fe3
        real(dp) :: Ca
        real(dp) :: Mg
        real(dp) :: NaCl
        real(dp) :: KCl
        real(dp) :: Hplus
        real(dp) :: OHmin
        real(dp) :: O2
    end type moleclist


    type bornmoleclist
        real(dp) :: pol
        real(dp) :: polCa
        real(dp) :: polMg
        real(dp) :: Na
        real(dp) :: Cl
        real(dp) :: K
        real(dp) :: Fe2 
        real(dp) :: Fe3
        real(dp) :: Ca
        real(dp) :: Mg
        real(dp) :: Hplus
        real(dp) :: OHmin
    end type bornmoleclist 

contains

    subroutine init_zero_moleclist(struct)
    
        type(moleclist), intent(inout) :: struct

        struct%sol = 0.0_dp
        struct%Na = 0.0_dp
        struct%Cl = 0.0_dp
        struct%K = 0.0_dp
        struct%Fe2 = 0.0_dp
        struct%Fe3 = 0.0_dp
        struct%Ca = 0.0_dp
        struct%Mg= 0.0_dp
        struct%NaCl = 0.0_dp
        struct%KCl = 0.0_dp
        struct%Hplus = 0.0_dp
        struct%OHmin = 0.0_dp
        struct%O2 = 0.0_dp

    end subroutine init_zero_moleclist


     function get_value_moleclist(struct,member) result(val)

        type(moleclist) , intent(in) :: struct 
        character(len=5),  intent(in) :: member
        real(dp) :: val

        select case (member)
            case ("sol") 
                val=struct%sol
            case ("Na") 
                val=struct%Na
            case ("Cl") 
                val=struct%Cl
            case ("K") 
                val=struct%K
            case ("Hplus") 
                val=struct%Hplus
            case ("OHmin") 
                val=struct%OHmin
            case ("Fe2") 
                val=struct%Fe2
            case ("Fe3") 
                val=struct%Fe3
            case ("Mg") 
                val=struct%Mg
            case default
                print*,"Wrong value member molecule list :  ",member
                stop
        end select      
            
    end function get_value_moleclist


    function sum_value_moleclist(struct) result(val)

        type(moleclist) , intent(in) :: struct
        real(dp) :: val

        
        val = struct%sol + struct%Na +struct%Cl + struct%K + struct%Fe2 + struct%Fe3 + &
        struct%Ca + struct%Mg+struct%NaCl + struct%KCl + struct%Hplus + struct%OHmin + struct%O2 
        
    end function sum_value_moleclist   

end module molecules

