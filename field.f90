module field
  
    !     .. variables
    use precision_definition

    implicit none
    
    real(dp), dimension(:), allocatable   :: xpol     ! volume fraction of polymer 
    real(dp), dimension(:,:), allocatable :: xpol_t   ! volume fraction of polymer in layer i of type t
    real(dp), dimension(:,:), allocatable :: rhopol   ! density monomer of polymer in layer i of type t
    real(dp), dimension(:,:), allocatable :: rhopolin 
    real(dp), dimension(:,:), allocatable :: rhopol_charge ! density chargeable monomer of polymer in layer i of type t
    real(dp), dimension(:), allocatable   :: rhoqpol  ! charge density  monomer of polymer in layer i 

    real(dp), dimension(:), allocatable :: xsol     ! volume fraction solvent
    real(dp), dimension(:), allocatable :: psi      ! electrostatic potential 
    real(dp), dimension(:), allocatable :: xNa      ! volume fraction of positive Na+ ion
    real(dp), dimension(:), allocatable :: xK       ! volume fraction of positive K+ ion
    real(dp), dimension(:), allocatable :: xFe2     ! volume fraction of positive Fe2+ ion 
    real(dp), dimension(:), allocatable :: xFe3     ! volume fraction of positive Fe3+ ion
    real(dp), dimension(:), allocatable :: xCa      ! volume fraction of positive Ca2+ ion
    real(dp), dimension(:), allocatable :: xMg      ! volume fraction of positive Mg2+ ion    
    real(dp), dimension(:), allocatable :: xNaCl    ! volume fraction of NaCl ion pair
    real(dp), dimension(:), allocatable :: xKCl     ! volume fraction of KCl  ion pair
    real(dp), dimension(:), allocatable :: xCl      ! volume fraction of negative ion
    real(dp), dimension(:), allocatable :: xHplus   ! volume fraction of Hplus
    real(dp), dimension(:), allocatable :: xOHmin   ! volume fraction of OHmin 
    real(dp), dimension(:), allocatable :: xO2      ! volume fraction of O2 

    real(dp), dimension(:), allocatable :: rhoq     ! total free charge density in units of vsol  
    real(dp), dimension(:), allocatable :: epsfcn   ! relative dielectric constant 
    real(dp), dimension(:), allocatable :: Depsfcn  ! relative derivative dielectric constant

    real(dp), dimension(:,:), allocatable   :: fdis    ! degree of dissociation of acid monomer and base monomer
                                                       ! acid: AH<=> A^- +H^+ f_A^-=fdis, base : BH^+<=> B+ H^+ f_B=fdis 
    real(dp), dimension(:,:), allocatable   :: fdisA   ! degree of dissociation of acid including condensed states
    real(dp), dimension(:,:), allocatable   :: fdisB   ! degree of dissociation
    real(dp), dimension(:,:,:), allocatable :: gdisA   ! degree of dissociation of acid including condensed states  
    real(dp), dimension(:,:,:), allocatable :: gdisB   ! degree of dissociation of base including condensed states  

    real(dp) :: q          ! normalization partion fnc polymer 
    real(dp) :: lnq        ! exponent of normalization partion fnc polymer 
    real(dp) :: lnproshift ! shift in exponent palpha

    real(dp), dimension(:), allocatable       :: rhoqphos       ! charged density of phosphate needed systype="nucl_ionbin_Mg"
    real(dp), dimension(:,:,:,:), allocatable   :: fdisPP       ! fraction  fdisPP(i,k,J,K)  
    real(dp), dimension(:,:), allocatable       :: fdisP2Mg     ! fraction  fdisP2Mg(i,k)  
    real(dp), dimension(:,:), allocatable       :: fdisP2Fe2    ! fraction  fdisP2Fe2(i,k)  
    real(dp), dimension(:,:), allocatable       :: fdisP2Fe3    ! fraction  fdisP2Fe3(i,k)  

    real(dp), dimension(:,:), allocatable   :: fdisPP_loc, fdisPP_loc_swap     ! fdisPP(J,K) local equivalent of fraction of fdisPP(i,k,J,K)  
    real(dp)                                :: fdisP2Mg_loc, fdisP2Mg_loc_swap ! fdisP2Mg    local equivalent of fraction of fdisP2Mg(i,k)     
    real(dp)                                :: fdisP2Fe2_loc, fdisP2Fe2_loc_swap ! fdisP2Fe2    local equivalent of fraction of fdisP2Fe2(i,k)     
    real(dp)                                :: fdisP2Fe3_loc, fdisP2Fe3_loc_swap ! fdisP2Fe2 
    
    real(dp), dimension(:,:,:), allocatable   :: fdisPPP_loc1 , fdisPPP_loc1_swap    ! fdisPPP(J,K,L)  
    real(dp), dimension(:,:,:), allocatable   :: fdisPPP_loc2 , fdisPPP_loc2_swap    ! fdisPPP(J,K,L)  
    real(dp), dimension(:,:,:), allocatable   :: fdisPPP_loc3 , fdisPPP_loc3_swap    ! fdisPPP(J,K,L)  


    real(dp), parameter  :: eps_val = 1.0e-7_dp 
    
    private :: eps_val 

contains

    subroutine allocate_field(Nx,Ny,Nz,nsegtypes)
 
        integer, intent(in) :: Nx,Ny,Nz,nsegtypes
        
        integer :: N
        integer :: ier(28), i

        N=Nx*Ny*Nz

        allocate(xpol(N),stat=ier(1))
        allocate(xpol_t(N,nsegtypes),stat=ier(2))
        allocate(rhopol(N,nsegtypes),stat=ier(3)) 
        allocate(rhopolin(N,nsegtypes),stat=ier(4)) 
        allocate(rhopol_charge(N,nsegtypes),stat=ier(5)) 
        allocate(rhoqpol(N),stat=ier(6)) 
        allocate(xsol(N),stat=ier(7))
        allocate(psi(N+2*Nx*Ny),stat=ier(8))    !allocate(psi(N),stat=ier(6))
        allocate(xNa(N),stat=ier(9))
        allocate(xK(N),stat=ier(10))
        allocate(xFe2(N),stat=ier(11)) 
        allocate(xFe3(N),stat=ier(12))
        allocate(xCa(N),stat=ier(13))
        allocate(xMg(N),stat=ier(14))
        allocate(xNaCl(N),stat=ier(15)) 
        allocate(xKCl(N),stat=ier(16)) 
        allocate(xCl(N),stat=ier(17)) 
        allocate(xHplus(N),stat=ier(18))
        allocate(xOHmin(N),stat=ier(19))
        allocate(xO2(N),stat=ier(20))
        allocate(rhoq(N),stat=ier(21))
        allocate(epsfcn(N),stat=ier(22))    
        allocate(Depsfcn(N),stat=ier(23))  
        allocate(fdis(N,nsegtypes),stat=ier(24))
        allocate(fdisA(N,8),stat=ier(25))
        allocate(fdisB(N,5),stat=ier(26))
        allocate(gdisA(N,4,nsegtypes),stat=ier(27))
        allocate(gdisB(N,3,nsegtypes),stat=ier(28))
        
        do i=1,28
            if( ier(i)/=0 ) then
                print*, 'Allocation error : stat =', ier(i),' for i= ',i
                stop
            endif
        enddo    
        
    end subroutine allocate_field


    subroutine deallocate_field()
        
        deallocate(xpol)
        deallocate(xpol_t)
        deallocate(rhopol)
        deallocate(rhoqpol)
        deallocate(rhopol_charge)
        deallocate(xsol)
        deallocate(psi)
        deallocate(xNa)
        deallocate(xK)
        deallocate(xFe2)
        deallocate(xFe3)
        deallocate(xCa)
        deallocate(xMg)
        deallocate(xNaCl) 
        deallocate(xKCl) 
        deallocate(xCl) 
        deallocate(xHplus)
        deallocate(xOHmin)
        deallocate(xO2)
        deallocate(rhoq)
        deallocate(epsfcn)
        deallocate(Depsfcn)
        deallocate(fdis)
        deallocate(fdisA)
        deallocate(fdisB)
        deallocate(gdisA)
        deallocate(gdisB)
        
    end subroutine deallocate_field


    ! set all densities to zero
    
    subroutine init_field()

        xpol=0.0_dp
        xpol_t=0.0_dp
        rhopol=0.0_dp
        rhoqpol=0.0_dp
        rhopol_charge=0.0_dp
        xsol=0.0_dp
        xNa=0.0_dp
        xK=0.0_dp
        xFe2=0.0_dp
        xFe3=0.0_dp
        xCa=0.0_dp
        xMg=0.0_dp
        xNaCl=0.0_dp 
        xKCl =0.0_dp
        xCl=0.0_dp
        xHplus=0.0_dp
        xOHmin=0.0_dp
        xO2=0.0_dp
        rhoq=0.0_dp
        psi=0.0_dp
        fdis=0.0_dp
        fdisA=0.0_dp
        fdisB=0.0_dp
        gdisA=0.0_dp
        gdisB=0.0_dp

    end subroutine init_field

    subroutine allocate_field_pairs(Nx,Ny,Nz,maxneigh,maxfdisPP,len_index_phos)

        use globals, only : systype 
        integer, intent(in) :: Nx,Ny,Nz,maxneigh, maxfdisPP,len_index_phos

        integer :: N, Nindex
       
        if(systype=="nucl_ionbin_Mg") then 

            N=Nx*Ny*Nz
            Nindex=len_index_phos
            allocate(rhoqphos(N))    
            allocate(fdisPP(Nindex,maxneigh,maxfdisPP,maxfdisPP)) 
            allocate(fdisP2Mg(Nindex,maxneigh)) 
            allocate(fdisP2Fe2(Nindex,maxneigh)) 
            allocate(fdisP2Fe3(Nindex,maxneigh))

        endif
 
        if(systype=="nucl_ionbin_MgA" .or. systype =="nucl_ionbin_Fe".or. &
            systype=="nucl_ionbin_Fe_ST".or.systype=="nucl_ionbin_Fe_ST_mu") then

            N=Nx*Ny*Nz
            allocate(rhoqphos(N))
            allocate(fdisPP_loc(maxfdisPP,maxfdisPP))
            allocate(fdisPP_loc_swap(maxfdisPP,maxfdisPP))

        endif

    end subroutine allocate_field_pairs

    subroutine allocate_field_triplets(maxfdisPPP)

        use globals, only : systype 
        integer, intent(in) :: maxfdisPPP

        if(systype=="nucl_ionbin_Fe".or.systype=="nucl_ionbin_Fe_ST".or.systype=="nucl_ionbin_Fe_ST_mu") then

            allocate(fdisPPP_loc1(maxfdisPPP,maxfdisPPP,maxfdisPPP))
            allocate(fdisPPP_loc1_swap(maxfdisPPP,maxfdisPPP,maxfdisPPP))
            allocate(fdisPPP_loc2(maxfdisPPP,maxfdisPPP,maxfdisPPP))
            allocate(fdisPPP_loc2_swap(maxfdisPPP,maxfdisPPP,maxfdisPPP))
            allocate(fdisPPP_loc3(maxfdisPPP,maxfdisPPP,maxfdisPPP))
            allocate(fdisPPP_loc3_swap(maxfdisPPP,maxfdisPPP,maxfdisPPP))

        endif

    end subroutine allocate_field_triplets

    subroutine init_field_pairs()
       
        use globals, only : systype
 
        rhoqphos=0.0_dp
        
        if(systype=="nucl_ionbin_Mg") then
            fdisPP=0.0_dp
            fdisP2Mg=0.0_dp
            fdisP2Fe2=0.0_dp
            fdisP2Fe3=0.0_dp
        endif

        if(systype=="nucl_ionbin_MgA") then
            fdisPP_loc=0.0_dp
            fdisP2Mg_loc=0.0_dp
            fdisPP_loc_swap=0.0_dp
            fdisP2Mg_loc_swap=0.0_dp
        endif

    end subroutine init_field_pairs

    subroutine init_field_triplets()
       
        use globals, only : systype
 
        if(systype=="nucl_ionbin_Fe".or.systype=="nucl_ionbin_Fe_ST".or.systype=="nucl_ionbin_Fe_ST_mu") then
            fdisPPP_loc1=0.0_dp
            fdisPPP_loc1_swap=0.0_dp
            fdisPPP_loc2=0.0_dp
            fdisPPP_loc2_swap=0.0_dp
            fdisPPP_loc3=0.0_dp
            fdisPPP_loc3_swap=0.0_dp
        endif

    end subroutine init_field_triplets

    !  compute routines 


    subroutine check_integral_rholpol_multi(sumrhopol, checkintegral)

        use volume, only : volcell
        use globals, only : nsize, nseg, nsegtypes

        real(dp), intent(inout) :: sumrhopol,checkintegral 
        integer :: t,i
        real(dp) :: intrhopol

        sumrhopol=0.0_dp
        do t=1,nsegtypes
            do i=1,nsize
                sumrhopol=sumrhopol+rhopol(i,t)
            enddo  
        enddo      
        sumrhopol=sumrhopol*volcell

        intrhopol=nseg

        checkintegral=sumrhopol-intrhopol

    end subroutine

    subroutine check_integral_rholpolAB(sumrhopol, checkintegral)

        use volume, only : volcell
        use globals, only : nsize, nseg

        real(dp), intent(inout) :: sumrhopol,checkintegral 
        integer :: i
        real(dp) :: intrhopol

        sumrhopol=0.0_dp
        do i=1,nsize
            sumrhopol=sumrhopol+(rhopol(i,1)+rhopol(i,2))
        enddo    
        sumrhopol=sumrhopol*volcell

        intrhopol=nseg  

        checkintegral=sumrhopol-intrhopol

    end subroutine
       

    subroutine charge_polymer()

        use globals, only : systype
        
        select case (systype) 
        case ("brush_mul","brush_mulnoVdW")

            call charge_polymer_multi()

        case ("brushdna","brushborn")

            call charge_polymer_dna()

        case ("nucl_ionbin")

            call charge_nucl_ionbin()

        case ("nucl_ionbin_sv")

            call charge_nucl_ionbin_sv()  

        case ("nucl_ionbin_Mg","nucl_ionbin_MgA","nucl_ionbin_Fe","nucl_ionbin_Fe_ST","nucl_ionbin_Fe_ST_mu")

            call charge_nucl_ionbin_Mg() 

        case ("elect")  

            call charge_polymer_binary()

        case default

            print*,"Error in charge_polymer subroutine"    
            print*,"Wrong value systype : ", systype
            stop

        end select  
       

    end subroutine charge_polymer


    subroutine charge_polymer_dna()

        use globals, only : nsize, nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, qpol, qpol_tot, tA

        integer :: i, t

        qpol_tot=0.0_dp
        do t=1,nsegtypes
            qpol(t)=0.0_dp
            if(t/=tA) then    
                do i=1,nsize
                    qpol(t)=qpol(t)+(fdis(i,t)*zpol(t,2)+(1.0_dp-fdis(i,t))*zpol(t,1))*rhopol(i,t)
                enddo
            else
                do i=1,nsize
                    qpol(t)=qpol(t)+ (-fdisA(i,1)+fdisA(i,4)+fdisA(i,6))*rhopol(i,tA)
                enddo
            endif    

            qpol(t)=qpol(t)*volcell
            qpol_tot=qpol_tot+qpol(t)
        enddo

    end subroutine charge_polymer_dna


    subroutine charge_nucl_ionbin()

        use globals, only : nsize, nsegtypes
        use volume, only : volcell
        use parameters, only : qpol, qpol_tot, tA
        use chains, only : ismonomer_chargeable, type_of_charge 

        integer :: i, t

        qpol_tot=0.0_dp
        
        do t=1,nsegtypes
            qpol(t)=0.0_dp
            if(ismonomer_chargeable(t)) then 
                if(type_of_charge(t)=="A") then   ! acid 
                    if(t/=ta) then
                        do i=1,nsize
                            qpol(t)=qpol(t)-gdisA(i,1,t)*rhopol(i,t)
                        enddo
                    else ! phosphate
                        do i=1,nsize
                            qpol(ta)=qpol(t)+ (-fdisA(i,1)+fdisA(i,4)+fdisA(i,6))*rhopol(i,t)
                        enddo    
                    endif
                else  ! base   
                    do i=1,nsize
                            qpol(t)=qpol(t)+gdisB(i,1,t)*rhopol(i,t)
                    enddo
                endif        
            endif    
            qpol(t)=qpol(t)*volcell
            qpol_tot=qpol_tot+qpol(t)
        enddo

    end subroutine charge_nucl_ionbin


    subroutine charge_nucl_ionbin_sv()

        use globals, only : nsize, nsegtypes
        use volume, only : volcell
        use parameters, only : qpol, qpol_tot, tA
        use chains, only : ismonomer_chargeable, type_of_charge

        integer :: i, t
        real(dp) :: qpoltmp

        qpol_tot=0.0_dp
     

        do t=1,nsegtypes
            qpol(t)=0.0_dp
           
            if(ismonomer_chargeable(t)) then
           
                if(type_of_charge(t)=="A") then   ! acid 
           
                    if(t/=ta) then
                        do i=1,nsize
                            qpol(t)=qpol(t)-gdisA(i,1,t)*rhopol_charge(i,t)
                        enddo
                    else ! phosphate
                        do i=1,nsize
                            qpoltmp=(-fdisA(i,1)+fdisA(i,4)+fdisA(i,6))*rhopol_charge(i,t)
                            qpol(t)=qpol(t)+ qpoltmp
                        enddo    
                    endif
                else  ! base   
                    do i=1,nsize
                        qpoltmp=gdisB(i,1,t)*rhopol_charge(i,t)
                        qpol(t)=qpol(t)+qpoltmp
                    enddo
                endif        
            endif  
            qpol(t)=qpol(t)*volcell
            qpol_tot=qpol_tot+qpol(t)
        enddo
        
    end subroutine charge_nucl_ionbin_sv


    ! computes avarage charge nucl residue for systype==nucl_ionbin_Mg and nucl_ionbin_MgA: 
    ! charge phosphate rhopolqphos(i) seperate computed 

    subroutine charge_nucl_ionbin_Mg()

        use globals, only : nsize, nsegtypes
        use volume, only : volcell
        use parameters, only : qpol, qpol_tot, tA
        use chains, only : ismonomer_chargeable, type_of_charge

        integer :: i, t
        real(dp) :: qpoltmp 

        qpol_tot=0.0_dp
     

        do t=1,nsegtypes
            qpol(t)=0.0_dp
           
            if(ismonomer_chargeable(t)) then
           
                if(type_of_charge(t)=="A") then   ! acid 
           
                    if(t/=ta) then
                        do i=1,nsize
                            qpol(t)=qpol(t)-gdisA(i,1,t)*rhopol_charge(i,t)
                        enddo
                    else ! phosphate
                        do i=1,nsize
                            qpol(t)=qpol(t)+ rhoqphos(i) ! !!!! units m
                        enddo    
                    endif
                else  ! base   
                    do i=1,nsize
                        qpoltmp=gdisB(i,1,t)*rhopol_charge(i,t)
                        qpol(t)=qpol(t)+qpoltmp
                    enddo
                endif        
            endif  
            qpol(t)=qpol(t)*volcell
            qpol_tot=qpol_tot+qpol(t)
        enddo

    end subroutine charge_nucl_ionbin_Mg    

    subroutine charge_polymer_multi()

        use globals, only : nsize, nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, qpol, qpol_tot

        integer :: i, t

        qpol_tot=0.0_dp
        do t=1,nsegtypes
            qpol(t)=0.0_dp
            do i=1,nsize
                qpol(t)=qpol(t)+(fdis(i,t)*zpol(t,2)+(1.0_dp-fdis(i,t))*zpol(t,1))*rhopol(i,t)
            enddo
            qpol(t)=qpol(t)*volcell
            qpol_tot=qpol_tot+qpol(t)
        enddo

    end subroutine charge_polymer_multi

    subroutine charge_polymer_binary()

        use globals, only : nsize
        use volume, only : volcell
        use parameters, only : zpolA, zpolB, qpolA,qpolB, qpol_tot

        integer :: i

        qpolA=0.0_dp
        qpolB=0.0_dp

        do i=1,nsize
            qpolA=qpolA+(zpolA(1)*fdisA(i,1)*rhopol(i,1)+zpolA(4)*fdisA(i,4)*rhopol(i,1))
            qpolB=qpolB+(zpolB(1)*fdisB(i,1)*rhopol(i,2)+zpolB(4)*fdisB(i,4)*rhopol(i,2))
        enddo

        qpolA=qpolA*volcell
        qpolB=qpolB*volcell
        qpol_tot=qpolA+qpolB

    end subroutine charge_polymer_binary

    ! .. post : return average charge of state of polymer

    subroutine average_charge_polymer()

        use globals, only : systype
        
        select case (systype) 
        case ("brush_mul","brush_mulnoVdW")

            call average_charge_polymer_multi()

        case ("brushdna","brushborn")

            call average_charge_polymer_dna()

        case ("nucl_ionbin")

            call average_charge_nucl_ionbin()

        case ("nucl_ionbin_sv")

            call average_charge_nucl_ionbin_sv()

        case ("nucl_ionbin_Mg","nucl_ionbin_MgA")

            call average_charge_nucl_ionbin_Mg()

        case ("nucl_ionbin_Fe","nucl_ionbin_Fe_ST","nucl_ionbin_Fe_ST_mu")
            
            call average_charge_nucl_ionbin_Fe()

        case ("elect","electA","electVdWAB","electdouble") 

            call average_charge_polymer_binary()

        case default

            print*,"Error in average_charge_polymer subroutine"    
            print*,"Wrong value systype : ", systype
            stop

        end select  

    end subroutine average_charge_polymer
        

    subroutine average_charge_polymer_dna()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, avfdis, avfdisA, tA
        use chains, only: type_of_monomer,ismonomer_chargeable

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t,k
        real(dp) :: sumrhopolt ! average density of polymer of type t 

        allocate(npol(nsegtypes))
        
        npol=0

        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   
            
        do t=1,nsegtypes
            avfdis(t)=0.0_dp
            if(ismonomer_chargeable(t)) then 
                sumrhopolt=npol(t)/volcell
                if(npol(t)/=0) then
                    if(t/=tA) then    
                        do i=1,nsize
                            avfdis(t)=avfdis(t)+(fdis(i,t)*zpol(t,2)+(1.0_dp-fdis(i,t))*zpol(t,1))*rhopol(i,t)
                        enddo
                        avfdis(t)=avfdis(t)/sumrhopolt        
                    else
                        do k=1,8
                            avfdisA(k)=0.0_dp
                            do i=1,nsize
                                avfdisA(k)=avfdisA(k)+fdisA(i,k)*rhopol(i,t)
                            enddo
                            avfdisA(k)=avfdisA(k)/sumrhopolt  
                        enddo
                        avfdis(t)=avfdisA(1)
                    endif       
                endif
            endif    
        enddo         

        deallocate(npol)    

    end subroutine average_charge_polymer_dna


    subroutine average_charge_nucl_ionbin()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, avfdis, avfdisA, tA, avgdisA, avgdisB
        use chains, only: type_of_monomer,ismonomer_chargeable

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t,k
        real(dp) :: sumrhopolt ! average density of polymer of type t 

        allocate(npol(nsegtypes))
        
        npol=0

        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   
            
        do t=1,nsegtypes
            ! init 
            avfdis(t)=0.0_dp ! A^-
            do k=1,4               ! A^-, AH, ANa, AK for AA that are acid
                avgdisA(t,k)=0.0_dp 
            enddo
            do k=1,3 !             ! BH^+, B, BHCl for AA that are base
                avgdisB(t,k)=0.0_dp
            enddo
                        
            if(ismonomer_chargeable(t)) then 
                sumrhopolt=npol(t)/volcell
                if(npol(t)/=0) then
                    if(t/=tA) then 
                        if(zpol(t,1)==0) then ! acid
                            do k=1,4
                                avgdisA(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisA(t,k)=avgdisA(t,k)+gdisA(i,k,t)*rhopol(i,t)
                                enddo
                                avgdisA(t,k)=avgdisA(t,k)/sumrhopolt  
                            enddo
                            avfdis(t)=zpol(t,2)*avgdisA(t,1) ! signed charged fraction   
                        else ! base
                            do k=1,3
                                avgdisB(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisB(t,k)=avgdisB(t,k)+gdisB(i,k,t)*rhopol(i,t)
                                enddo
                                avgdisB(t,k)=avgdisB(t,k)/sumrhopolt 
                            enddo
                            avfdis(t)=zpol(t,1)*avgdisB(t,1) 
                        endif            
                    else
                        do k=1,8
                            avfdisA(k)=0.0_dp
                            do i=1,nsize
                                avfdisA(k)=avfdisA(k)+fdisA(i,k)*rhopol(i,t)
                            enddo
                            avfdisA(k)=avfdisA(k)/sumrhopolt  
                        enddo
                        avfdis(t)=avfdisA(1)
                    endif               
                endif
            endif    
        enddo         

        deallocate(npol)    

    end subroutine average_charge_nucl_ionbin

    ! compute average charge fraction of  nucleosome for systype nucl_ionbin_Mg

    subroutine average_charge_nucl_ionbin_Mg()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, tA, avfdis, avfdisA, avgdisA, avgdisB
        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, PhosFe2, PhosFe3
        use parameters, only : avfdisPP, avfdisP2Mg, avfdisP2Fe2, avfdisP2Fe3
        use chains, only: type_of_monomer,ismonomer_chargeable

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t,k,JJ, KK
        real(dp) :: sumrhopolt ! average density of polymer of type t 

        allocate(npol(nsegtypes))
        
        npol=0

        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   
            
        do t=1,nsegtypes
            ! init 
            avfdis(t)=0.0_dp ! A^-
            do k=1,4               ! A^-, AH, ANa, AK for AA that are acid
                avgdisA(t,k)=0.0_dp 
            enddo
            do k=1,3 !             ! BH^+, B, BHCl for AA that are base
                avgdisB(t,k)=0.0_dp
            enddo

            if(ismonomer_chargeable(t)) then 
                sumrhopolt=npol(t)/volcell
                if(npol(t)/=0) then
                    if(t/=tA) then 
                        if(zpol(t,1)==0) then ! acid
                            do k=1,4
                                avgdisA(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisA(t,k)=avgdisA(t,k)+gdisA(i,k,t)*rhopol_charge(i,t)
                                enddo
                                avgdisA(t,k)=avgdisA(t,k)/sumrhopolt  
                            enddo
                            avfdis(t)=zpol(t,2)*avgdisA(t,1) ! signed charged fraction   
                        else ! base
                            do k=1,3
                                avgdisB(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisB(t,k)=avgdisB(t,k)+gdisB(i,k,t)*rhopol_charge(i,t)
                                enddo
                                avgdisB(t,k)=avgdisB(t,k)/sumrhopolt 
                            enddo
                            avfdis(t)=zpol(t,1)*avgdisB(t,1) 
                        endif            
                    else
                        ! t=tA phophates
 
                        do k=1,12
                            avfdisA(k)=0.0_dp
                        enddo   
                            
                        ! charged phosphates
                        do JJ=1,7
                            KK=Phos
                            avfdisA(1)=avfdisA(1) + avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
                        enddo
 
                        ! protonated phosphates
                        do JJ=1,7
                            KK=PhosH    
                            avfdisA(2) = avfdisA(2)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
                        enddo
 
                        ! Na bound  phosphates
                        do JJ=1,7
                            KK=PhosNa
                            avfdisA(3) = avfdisA(3)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
                        enddo
    
                        ! Ca bound phosphates             
                        avfdisA(4) = 0.0_dp

                        ! P2Ca bound phosphates
                        avfdisA(5) = 0.0_dp

                        ! K bound  phosphates
                        do JJ=1,7
                            KK=PhosK
                            avfdisA(8) = avfdisA(8)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
                        enddo

                        ! Mg bound phosphates
                        do JJ=1,7
                            KK=PhosMg
                            avfdisA(6) = avfdisA(6)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
                        enddo

                         ! Fe2 bound phosphates
                        do JJ=1,7
                            KK=PhosFe2
                            avfdisA(9) = avfdisA(9)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
                        enddo

                        ! Fe3 bound phosphates
                        do JJ=1,7
                            KK=PhosFe3
                            avfdisA(11) = avfdisA(11)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
                        enddo

                        ! P2Mg bound phophates 
                        avfdisA(7)=2.0_dp*avfdisP2Mg
                        
                        ! P2Fe2 bound phophates 
                        avfdisA(10)=2.0_dp*avfdisP2Fe2

                        ! P2Fe3 bound phophates 
                        avfdisA(12)=2.0_dp*avfdisP2Fe3

                        do k=1,12
                            avfdisA(k)=avfdisA(k)/2.0_dp
                        enddo  
                        ! divide by 2 because avfdisPP fraction of pairs i.e normed with total number of pairs!
                         
                        avfdis(ta)= - avfdisA(1)+avfdisA(4)+avfdisA(6)  ! signed charged fraction  

                    endif               
                endif
            endif    
        enddo         

        deallocate(npol)    

    end subroutine average_charge_nucl_ionbin_Mg


     ! compute average charge fraction of  nucleosome for systype nucl_ionbin_Fe

    subroutine average_charge_nucl_ionbin_Fe()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, tA, avfdis, avgdisA, avgdisB, avfdisA
        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, PhosFe2, PhosFe3
        use parameters, only : avfdisA_pairs,  avfdisA_triplets
        use chains, only: type_of_monomer,ismonomer_chargeable

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t,k
        real(dp) :: sumrhopolt ! average density of polymer of type t 

        allocate(npol(nsegtypes))
        
        npol=0

        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   
            
        do t=1,nsegtypes
            ! init 
            avfdis(t)=0.0_dp ! A^-
            do k=1,4               ! A^-, AH, ANa, AK for AA that are acid
                avgdisA(t,k)=0.0_dp 
            enddo
            do k=1,3 !             ! BH^+, B, BHCl for AA that are base
                avgdisB(t,k)=0.0_dp
            enddo

            if(ismonomer_chargeable(t)) then 
                sumrhopolt=npol(t)/volcell
                if(npol(t)/=0) then
                    if(t/=tA) then 
                        if(zpol(t,1)==0) then ! acid
                            do k=1,4
                                avgdisA(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisA(t,k)=avgdisA(t,k)+gdisA(i,k,t)*rhopol_charge(i,t)
                                enddo
                                avgdisA(t,k)=avgdisA(t,k)/sumrhopolt  
                            enddo
                            avfdis(t)=zpol(t,2)*avgdisA(t,1) ! signed charged fraction   
                        else ! base
                            do k=1,3
                                avgdisB(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisB(t,k)=avgdisB(t,k)+gdisB(i,k,t)*rhopol_charge(i,t)
                                enddo
                                avgdisB(t,k)=avgdisB(t,k)/sumrhopolt 
                            enddo
                            avfdis(t)=zpol(t,1)*avgdisB(t,1) 
                        endif            
                    else
                        ! t=tA phophates
                        
                        call average_charge_fraction_phos_pairs(avfdisA_pairs)
                        call average_charge_fraction_phos_triplets(avfdisA_triplets)
                        call average_charge_fraction_phos(avfdisA_pairs,avfdisA_triplets, avfdisA)
                         
                        avfdis(ta) = - avfdisA(1)+avfdisA(4)+avfdisA(6)  ! signed charged fraction   

                    endif               
                endif
            endif    
        enddo         

        deallocate(npol)    

    end subroutine average_charge_nucl_ionbin_Fe

    !  Compute average charge fraction orchemical state of phosphate  monomer belong to a phospahtres pairs pairs
    !  avfdisPP fraction of pairs i.e normed with total number of pairs!
    !  index of avfdispairs <=> index avfdisA 
    !  1 == P^- 2 == PH  , 3  == PNa,   4  == PCa,  5  == P2Ca,  6  == PMg, 7 == P2Mg 
    !  8 == PK, 9 == PFe2, 10 == P2Fe2, 11 == PFe3, 12 == P2Fe3, 13 == P3Fe3  
    
    subroutine average_charge_fraction_phos_pairs(avfdisA_pairs)

        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, PhosFe2, PhosFe3
        use parameters, only : avfdisPP, avfdisP2Mg, avfdisP2Fe2, avfdisP2Fe3
        

        real(dp), dimension(:) ,intent(inout) :: avfdisA_pairs
        real(dp) :: check_val

        ! local variables
        integer :: k, JJ, KK    
        integer :: dim_avfdisA_pairs,  dim_avfdisPP

        dim_avfdisA_pairs=size(avfdisA_pairs)
        dim_avfdisPP=size(avfdisPP,dim=1) ! dim1=dim2

        ! init 
        do k=1,dim_avfdisA_pairs
            avfdisA_pairs(k)=0.0_dp
        enddo   
            
        ! charged phosphates
        do JJ=1,dim_avfdisPP
            KK=Phos
            avfdisA_pairs(1)=avfdisA_pairs(1) + avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
        enddo

        ! protonated phosphates
        do JJ=1,dim_avfdisPP
            KK=PhosH    
            avfdisA_pairs(2) = avfdisA_pairs(2)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
        enddo

        ! Na bound  phosphates
        do JJ=1,dim_avfdisPP
            KK=PhosNa
            avfdisA_pairs(3) = avfdisA_pairs(3)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
        enddo

        ! Ca bound phosphates             
        avfdisA_pairs(4) = 0.0_dp

        ! P2Ca bound phosphates
        avfdisA_pairs(5) = 0.0_dp

        ! K bound  phosphates
        do JJ=1,dim_avfdisPP
            KK=PhosK
            avfdisA_pairs(8) = avfdisA_pairs(8)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
        enddo

        ! Mg bound phosphates
        do JJ=1,dim_avfdisPP
            KK=PhosMg
            avfdisA_pairs(6) = avfdisA_pairs(6)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
        enddo

            ! Fe2 bound phosphates
        do JJ=1,dim_avfdisPP
            KK=PhosFe2
            avfdisA_pairs(9) = avfdisA_pairs(9)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
        enddo

        ! Fe3 bound phosphates
        do JJ=1,dim_avfdisPP
            KK=PhosFe3
            avfdisA_pairs(11) = avfdisA_pairs(11)+avfdisPP(JJ,KK)+avfdisPP(KK,JJ)
        enddo

        ! P2Mg bound phophates 
        avfdisA_pairs(7) = 2.0_dp*avfdisP2Mg
        
        ! P2Fe2 bound phophates 
        avfdisA_pairs(10) = 2.0_dp*avfdisP2Fe2

        ! P2Fe3 bound phophates 
        avfdisA_pairs(12) = 2.0_dp*avfdisP2Fe3

        avfdisA_pairs = avfdisA_pairs /2
                        
        ! divide by 2 because avfdisPP fraction of pairs i.e normed with total number of pairs!
        ! Here norming  with total number of monomer that are part of a pair
         
        check_val = sum(avfdisA_pairs) 
         if(abs(check_val-1.0_dp)> eps_val) then 
            print*,"sum avfdisA_pairs not equal to 1"
            print*,"sum avfdisA_pairs = ",check_val       
        endif        

        print*,"sum avfdisA_pairs=",sum(avfdisA_pairs) 
    
    end subroutine average_charge_fraction_phos_pairs

    subroutine average_charge_fraction_phos_triplets(avfdisA_triplets)

        use parameters, only :  Phos, PhosH, PhosK, PhosNa, PhosMg, PhosFe2, PhosFe3
        use parameters, only : Phos2Mg, Phos2Fe2, Phos2Fe3, Phos3Fe3
        use parameters, only : avfdisPPP
        
        real(dp), dimension(:), intent(inout) :: avfdisA_triplets

        real(dp) :: check_val 
        real(dp) :: sum_avfJ
        integer :: dim_avfdisA_triplets,  dim_avfdisPPP
        integer :: J, K, L
        real(dp), dimension(:), allocatable :: avfdisA_triplets_tmp

        dim_avfdisA_triplets=size(avfdisA_triplets)
        dim_avfdisPPP=size(avfdisPPP,dim=1) ! dim1=dim2=dim3

        allocate(avfdisA_triplets_tmp(dim_avfdisPPP))

        ! init  
        avfdisA_triplets = 0.0_dp ! im
        avfdisA_triplets_tmp =0.0_dp
                    
        do J=1, dim_avfdisPPP
            sum_avfJ =0.0_dp
            do K=1,dim_avfdisPPP
                do L=1,dim_avfdisPPP 
                    sum_avfJ = sum_avfJ + avfdisPPP(J,K,L)+avfdisPPP(K,J,L)+avfdisPPP(K,L,J) 
                enddo
            enddo
            avfdisA_triplets_tmp(J)= sum_avfJ/3.0_dp 
        enddo        

        check_val = sum(avfdisA_triplets_tmp)

        if(abs(check_val-1.0_dp)> eps_val) then 
            print*,"sum avfdisA_triplets_tmp not equal to 1"
            print*,"sum avfdisA_triplets_tmp = ",check_val       
        endif        

        ! reshuffle indices
        !  1 == P^-, 2 == PH  , 3  == PNa,   4  == PCa,  5  == P2Ca,  6  == PMg, 7 == P2Mg 
        !  8 == PK,  9 == PFe2, 10 == P2Fe2, 11 == PFe3, 12 == P2Fe3, 13 == P3Fe3  

        avfdisA_triplets(1) = avfdisA_triplets_tmp(Phos)
        avfdisA_triplets(2) = avfdisA_triplets_tmp(PhosH)
        avfdisA_triplets(3) = avfdisA_triplets_tmp(PhosNa)
        avfdisA_triplets(4) = 0.0_dp
        avfdisA_triplets(5) = 0.0_dp
        avfdisA_triplets(6) = avfdisA_triplets_tmp(PhosMg)
        avfdisA_triplets(7) = avfdisA_triplets_tmp(Phos2Mg)
        avfdisA_triplets(8) = avfdisA_triplets_tmp(PhosK)
        avfdisA_triplets(9) = avfdisA_triplets_tmp(PhosFe2)
        avfdisA_triplets(10) = avfdisA_triplets_tmp(Phos2Fe2)
        avfdisA_triplets(11) = avfdisA_triplets_tmp(PhosFe3)
        avfdisA_triplets(12) = avfdisA_triplets_tmp(Phos2Fe3) 
        avfdisA_triplets(13) = avfdisA_triplets_tmp(Phos3Fe3)
        
        check_val = sum(avfdisA_triplets)

        if(abs(check_val-1.0_dp)> eps_val) then 
            print*,"sum avfdisA_triplets not equal to 1"
            print*,"sum avfdisA_triplets = ",check_val       
        endif        

        print*,"sum avfdisA_triplets = ",check_val
        
    end subroutine average_charge_fraction_phos_triplets


    subroutine average_charge_fraction_phos(avfdisA_pairs, avfdisA_triplets, avfdisA)

        real(dp), dimension(:), intent(in) :: avfdisA_pairs
        real(dp), dimension(:), intent(in) :: avfdisA_triplets
        real(dp), dimension(:), intent(inout) :: avfdisA

        real(dp) :: check_val, numpairs, numtriplets

        ! init   
        avfdisA = 0.0_dp            

        numpairs = numbers_pairs()
        numtriplets = numbers_triplets()

        
        avfdisA = 2.0_dp * numpairs * avfdisA_pairs + 3.0_dp * numtriplets * avfdisA_triplets
        avfdisA = avfdisA / ( 2.0_dp * numpairs  + 3.0_dp * numtriplets )
   
        check_val = sum(avfdisA)

        if(abs(check_val-1.0_dp)> eps_val) then 
            print*,"sum avfdisA not equal to 1"
            print*,"sum avfdisA = ",check_val       
        endif        

        print*,"sum avfdisA = ",check_val
        
    end subroutine average_charge_fraction_phos

    function numbers_pairs() result(sumrhopairs)
        
        use globals, only    : nseg, local_conf
        use chains, only      : type_of_monomer, nneigh
        use parameters, only : ta

        real(dp) :: sumrhopairs

        ! local variables
        integer :: c, s, j

        sumrhopairs = 0.0_dp

        do c=local_conf,local_conf       
            do s=1,nseg
                if(type_of_monomer(s)==ta) then
                    do j=1,nneigh(s,c)
                        sumrhopairs = sumrhopairs + 1.0_dp/(2.0_dp*nneigh(s,c)) 
                    enddo      
                endif
            enddo
        enddo

        sumrhopairs = sumrhopairs                 ! * 2.0_dp / 2.0_dp
        
    end function numbers_pairs

    function numbers_triplets() result(sumrhotriplets)
        
        use globals, only    : nseg, local_conf
        use chains, only      : type_of_monomer, ntriplet
        use parameters, only : ta

        real(dp) ::  sumrhotriplets 

        ! local variables
        integer :: c, s, j

        sumrhotriplets = 0.0_dp

        do c=local_conf,local_conf       
            do s=1,nseg
                if(type_of_monomer(s)==ta) then
                    do j=1,ntriplet(s,c) 
                        sumrhotriplets  = sumrhotriplets  + 1.0_dp/(6.0_dp*ntriplet(s,c))
                    enddo   
                endif
            enddo
        enddo

        sumrhotriplets = sumrhotriplets * 2.0_dp  ! *  3.0_dp /3.0_dp

    end function numbers_triplets


    ! compute average charge fractionof nucleosome for systype nucl_ionbin_sv

    subroutine average_charge_nucl_ionbin_sv()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, tA, avfdis, avfdisA, avgdisA, avgdisB
        use chains, only: type_of_monomer,ismonomer_chargeable

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t,k
        real(dp) :: sumrhopolt ! average density of polymer of type t 

        allocate(npol(nsegtypes))
        
        npol=0

        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   
            
        do t=1,nsegtypes
            ! init 
            avfdis(t)=0.0_dp ! A^-
            do k=1,4               ! A^-, AH, ANa, AK for AA that are acid
                avgdisA(t,k)=0.0_dp 
            enddo
            do k=1,3 !             ! BH^+, B, BHCl for AA that are base
                avgdisB(t,k)=0.0_dp
            enddo
                        
            if(ismonomer_chargeable(t)) then 
                sumrhopolt=npol(t)/volcell
                if(npol(t)/=0) then
                    if(t/=tA) then 
                        if(zpol(t,1)==0) then ! acid
                            do k=1,4
                                avgdisA(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisA(t,k)=avgdisA(t,k)+gdisA(i,k,t)*rhopol_charge(i,t)
                                enddo
                                avgdisA(t,k)=avgdisA(t,k)/sumrhopolt  
                            enddo
                            avfdis(t)=zpol(t,2)*avgdisA(t,1) ! signed charged fraction   
                        else ! base
                            do k=1,3
                                avgdisB(t,k)=0.0_dp
                                do i=1,nsize
                                    avgdisB(t,k)=avgdisB(t,k)+gdisB(i,k,t)*rhopol_charge(i,t)
                                enddo
                                avgdisB(t,k)=avgdisB(t,k)/sumrhopolt 
                            enddo
                            avfdis(t)=zpol(t,1)*avgdisB(t,1) 
                        endif            
                    else
                        do k=1,8
                            avfdisA(k)=0.0_dp
                            do i=1,nsize
                                avfdisA(k)=avfdisA(k)+fdisA(i,k)*rhopol_charge(i,t)
                            enddo
                            avfdisA(k)=avfdisA(k)/sumrhopolt  
                        enddo
                        avfdis(t)=avfdisA(1)
                    endif               
                endif
            endif    
        enddo         

        deallocate(npol)    

    end subroutine average_charge_nucl_ionbin_sv


    ! computes local nucleosome charge distrubution fo systype=nucl_ionbin_sv"

    subroutine distribution_charge_nucl_ionbin_sv(qpol_local)

        use globals, only : nsize,nsegtypes
        use parameters, only : zpol, tA
        use chains, only: ismonomer_chargeable

        real(dp), dimension(:), allocatable, intent(inout)  :: qpol_local

       
        integer :: i,t
      
        ! init 

        qpol_local=0.0_dp


        do i=1,nsize
            do t=1,nsegtypes
    
                if(ismonomer_chargeable(t)) then 
                    if(t/=tA) then ! not phosthate 
                        if(zpol(t,1)==0) then ! acid
                            qpol_local(i)=qpol_local(i)- gdisA(i,1,t)*rhopol_charge(i,t)
                        else ! base
                            qpol_local(i)=qpol_local(i)+gdisB(i,1,t)*rhopol_charge(i,t)
                        endif            
                    else
                        qpol_local(i)=qpol_local(i)- fdisA(i,1)*rhopol_charge(i,t)
                    endif               
                endif

            enddo    
        enddo            

    end subroutine distribution_charge_nucl_ionbin_sv


    subroutine average_charge_polymer_multi()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, avfdis
        use chains, only: type_of_monomer,ismonomer_chargeable

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t
        real(dp) :: sumrhopolt ! average density of polymer of type t 

        allocate(npol(nsegtypes))
        
        npol=0
        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   

        do t=1,nsegtypes
            avfdis(t)=0.0_dp
            if(ismonomer_chargeable(t)) then 
                sumrhopolt=npol(t)/volcell
                if(npol(t)/=0) then
                    avfdis(t)=0.0_dp
                    do i=1,nsize
                        avfdis(t)=avfdis(t)+(fdis(i,t)*zpol(t,2)+(1.0_dp-fdis(i,t))*zpol(t,1))*rhopol(i,t)
                    enddo
                    avfdis(t)=avfdis(t)/sumrhopolt       
                else
                    avfdis(t)=0.0_dp
                endif
            endif    
        enddo         

        deallocate(npol)    

    end subroutine average_charge_polymer_multi
        

    subroutine average_charge_polymer_binary()
        
        use globals, only : nseg,nsize
        use volume, only : volcell
        use parameters
        use chains, only : isAmonomer

        integer :: i,s,k
        integer   :: npolA,npolB
        integer, parameter :: A=1, B=2
        real(dp) :: sumrhopolA, sumrhopolB ! average density of polymer of type A and B
        ! .. number of A and B monomors 
        npolA=0
        do s=1,nseg
           if(isAmonomer(s).eqv..true.) then
              npolA=npolA+1
           endif
        enddo
        npolB=nseg-npolA
        sumrhopolA=npolA/volcell
        sumrhopolB=npolB/volcell
          

        if(npolA/=0) then
           do k=1,5
              avfdisA(k)=0.0_dp
              do i=1,nsize
                 avfdisA(k)=avfdisA(k)+fdisA(i,k)*rhopol(i,A)
              enddo
              avfdisA(k)=avfdisA(k)/sumrhopolA
           enddo
        else
           avfdisA=0.0_dp
        endif

        if(npolB/=0) then
           do k=1,5
              avfdisB(k)=0.0_dp
              do i=1,nsize
                 avfdisB(k)=avfdisB(k)+fdisB(i,k)*rhopol(i,B)
              enddo
              avfdisB(k)=avfdisB(k)/sumrhopolB
           enddo
        else
           do k=1,5
              avfdisB(k)=0.0_dp
           enddo
        endif

    end subroutine average_charge_polymer_binary

    !     .. compute average height of denisty provile
    !     .. first moment of density profile 
  
    function average_height_z(rho) result(meanz)

        use volume, only : nz,nx,ny,delta,linearIndexFromCoordinate

        real(dp), intent(in) :: rho(:)
        real(dp) :: meanz    

        integer :: ix, iy, iz, id
        real(dp) :: sumrhoz, sumrho

        sumrhoz = 0.0_dp
        meanz= 0.0_dp

        do iz = 1, nz

            sumrho = 0.0_dp       
            do ix=1, nx
                do iy=1, ny
                    call linearIndexFromCoordinate(ix,iy,iz ,id)
                    sumrho=sumrho+ rho(id)
                enddo
            enddo

            meanz=meanz+sumrho*(iz-0.5_dp)*delta
            sumrhoz=sumrhoz+sumrho
        enddo

        if(sumrhoz>0.0_dp) then 
            meanz=meanz/sumrhoz
        else
            meanz=0.0_dp
        endif

    end function average_height_z

    !     .. compute average of density or volume fraction profile in z-direction
  
    subroutine average_density_z(xvol,xvolz,meanz)

        use volume, only : nz,nx,ny,delta,linearIndexFromCoordinate

        real(dp), intent(in) :: xvol(:)
        real(dp), intent(out) :: xvolz(:)
        real(dp), intent(out), optional :: meanz

        integer :: ix, iy, iz, id
        real(dp) :: sumrhoz, sumxvol

        if(present(meanz)) then 

            sumrhoz = 0.0_dp
            do iz = 1, nz
                sumxvol = 0.0_dp       
                do ix=1, nx
                    do iy=1, ny
                        call linearIndexFromCoordinate(ix,iy,iz ,id)
                        sumxvol=sumxvol+ xvol(id)
                    enddo
                enddo
                xvolz(iz)=sumxvol/(1.0_dp*nx*ny)

                meanz=meanz+sumxvol*(iz-0.5_dp)*delta
                sumrhoz=sumrhoz+sumxvol
            enddo

            if(sumrhoz>0.0_dp) then 
                meanz=meanz/sumrhoz
            else
                meanz=0.0_dp
            endif
        else

            sumrhoz = 0.0_dp
            do iz = 1, nz
                sumxvol = 0.0_dp       
                do ix=1, nx
                    do iy=1, ny
                        call linearIndexFromCoordinate(ix,iy,iz ,id)
                        sumxvol=sumxvol+ xvol(id)
                    enddo
                enddo
                xvolz(iz)=sumxvol/(1.0_dp*nx*ny)
            enddo

        endif
            
    end subroutine average_density_z


    ! Computes ion_exces , gamma_i 
    ! gamma_i = \int dV (\rho_i(r) -\rho_bulk)
    ! input  xion     : array of real(dp) containing volume fraction of ion
    !        xionbulk : real(dp) bulk volume fraction of ion
    !        vol      : real(dp) relative vol ion : true volume vol*vsol
    ! output assigns : ionexcess

    function fcn_ion_excess(xion,xionbulk,vol) result(ionexcess)

        use volume, only : delta
        use globals, only : nsize
        use parameters, only : vsol

        real(dp), intent(in) :: xion(:)
        real(dp), intent(in) :: xionbulk
        real(dp), intent(in) :: vol  ! volume ion divide by vsol 

        real(dp) :: ionexcess

        integer :: i 

        ionexcess=0.0_dp
        do i=1,nsize
            ionexcess=ionexcess+xion(i)
        enddo  
        ionexcess=(ionexcess -xionbulk*nsize)*(delta**3)/(vol*vsol)

    end function fcn_ion_excess

    ! pre : make_ion_excess needed to be called before ion_excess
    ! pre : ion_excess needs to be computed. make_ion_excess called first
    ! input numberelem :  array contain the total elements of all types
    ! output assigns : beta_ion_excess  
    ! definitions beta_i = z_i * \Gamma_i / |qnucl|
  
    subroutine make_beta(numberelem)

        use parameters, only : ion_excess, beta_ion_excess, avgdisA,avgdisB,avfdisA
        use parameters, only : index_Phos=>ta ! index of phosphate 
        use globals, only : nsegtypes
        use molecules, only : moleclist, init_zero_moleclist
    !    use chains, only : mapping_num_to_char

        ! argument list

        real(dp), dimension(:) , allocatable :: numberelem ! array contain the total elements of all types
        
        ! local variable

        type(moleclist) :: ion_excess_tot ,ion_excess_ads
        real(dp) :: qnucl
        integer :: t

        ! calculate ion_excess_tot
        
        ! init 
        call init_zero_moleclist(ion_excess_ads) 
        call init_zero_moleclist(ion_excess_tot)
        call init_zero_moleclist(beta_ion_excess)
        
        ! ion_excess adsorped

        do t=1,nsegtypes
            ion_excess_ads%Na = ion_excess_ads%Na  + (avgdisA(t,3) * numberelem(t)) 
            ion_excess_ads%K  = ion_excess_ads%K   + (avgdisA(t,4) * numberelem(t))
            ion_excess_ads%Cl = ion_excess_ads%Cl  + (avgdisB(t,3) * numberelem(t)) 
        enddo
  
        ion_excess_ads%Na = ion_excess_ads%Na+ avfdisA(3) * numberelem(index_Phos) ! Na-phosphate 
        ion_excess_ads%K  = ion_excess_ads%K + avfdisA(8) * numberelem(index_Phos) ! K-phosphate
        ion_excess_ads%Mg =       (avfdisA(6)+avfdisA(7)) * numberelem(index_Phos) ! Mg-phosphate
        ion_excess_ads%Fe2 =      (avfdisA(9)+avfdisA(10)) * numberelem(index_Phos) ! Fe2-phosphate
        ion_excess_ads%Fe3 =      (avfdisA(11)+avfdisA(12)+avfdisA(13)) * numberelem(index_Phos) ! Fe3-phosphate
    
        
        ! calculate ion_excess_tot = sum of free adsorped ion excess
        
        ion_excess_tot%Na  = ion_excess%Na  + ion_excess_ads%Na   
        ion_excess_tot%K   = ion_excess%K   + ion_excess_ads%K   
        ion_excess_tot%Cl  = ion_excess%Cl  + ion_excess_ads%Cl   
        ion_excess_tot%Mg  = ion_excess%Mg  + ion_excess_ads%Mg
        ion_excess_tot%Fe2 = ion_excess%Fe2 + ion_excess_ads%Fe2
        ion_excess_tot%Fe3 = ion_excess%Fe3 + ion_excess_ads%Fe3
        
        ! Calculate qnucl 
        qnucl = abs(ion_excess_tot%Na + ion_excess_tot%K - ion_excess_tot%Cl + 2.0_dp*ion_excess_tot%Mg +& 
                2.0_dp*ion_excess_tot%Fe2+ 3.0_dp*ion_excess_tot%Fe3 )
        ! Calculate individual betas

        beta_ion_excess%Na  =  ion_excess_tot%Na / qnucl
        beta_ion_excess%K   =  ion_excess_tot%K  / qnucl
        beta_ion_excess%Cl  = -ion_excess_tot%Cl / qnucl
        beta_ion_excess%Mg  = 2.0_dp*ion_excess_tot%Mg / qnucl
        beta_ion_excess%Fe2 = 2.0_dp*ion_excess_tot%Fe2 / qnucl
        beta_ion_excess%Fe3 = 3.0_dp*ion_excess_tot%Fe3 / qnucl
            
    end subroutine make_beta



    subroutine make_ion_excess

        use parameters, only : vNa,vK,vMg,vCl,vCa,vFe2,vFe3
        use parameters, only : xbulk,ion_excess,sum_ion_excess

        ion_excess%Na=fcn_ion_excess(xNa,xbulk%Na,vNa)
        ion_excess%Cl=fcn_ion_excess(xCl,xbulk%Cl,vCl)
        ion_excess%K =fcn_ion_excess(xK,xbulk%K,vK)
        ion_excess%Mg=fcn_ion_excess(xMg,xbulk%Mg,vMg)
        ion_excess%Ca=fcn_ion_excess(xCa,xbulk%Ca,vCa)
        ion_excess%Fe2=fcn_ion_excess(xFe2,xbulk%Fe2,vFe2)
        ion_excess%Fe3=fcn_ion_excess(xFe3,xbulk%Fe3,vFe3)
        ion_excess%Hplus=fcn_ion_excess(xHplus,xbulk%Hplus,1.0_dp)
        ion_excess%OHmin=fcn_ion_excess(xOHmin,xbulk%OHmin,1.0_dp)
       
        ! sum of ion_excess weighted with valence of ion
       
        sum_ion_excess = ion_excess%Na -ion_excess%Cl + ion_excess%K +2.0_dp*ion_excess%Ca+2.0_dp*ion_excess%Mg +&
         ion_excess%Hplus -ion_excess%OHmin + 2.0_dp*ion_excess%Fe2 +3.0_dp*ion_excess%Fe3
        
    end subroutine make_ion_excess


    ! pre : make_ion_excess needed to be called before ion_excess
    ! pre : ion_excess needs to be computed. make_ion_excess called first
    ! input numberelem :  array containing the total elements of all types
    ! output assigns : beta_ion_excess  
    ! definitions beta_i = z_i * \Gamma_i / |qnucl|
  
    subroutine make_beta_old(numberelem)

        use parameters, only : ion_excess, beta_ion_excess, avgdisA, avgdisB, avfdisA
        use parameters, only : index_Phos=>ta ! index of phosphate 
        use globals, only : nsegtypes
        use molecules, only : moleclist

        ! argument list

        real(dp), dimension(:) , allocatable :: numberelem ! array contain the total elements of all types
        
        ! local variable

        type(moleclist) :: ion_excess_tot 
        real(dp) :: qnucl
        integer :: t


        ! run part

        do t=1,nsegtypes

              ion_excess_tot%Na = ion_excess_tot%Na + ion_excess%Na + (avgdisA(t,3) * numberelem(t)) 
              ion_excess_tot%K  = ion_excess_tot%K  + ion_excess%K  + (avgdisA(t,4) * numberelem(t))
              ion_excess_tot%Cl = ion_excess_tot%Cl + ion_excess%Cl + (avgdisB(t,3) * numberelem(t)) 
              print*,"t=", t," ion_excess_tot%Cl=", ion_excess_tot%Cl,"ion_excess Cl bound =",&
              avgdisB(t,3)*numberelem(t) 
        enddo

        ion_excess_tot%Na = ion_excess_tot%Na+ (avfdisA(3)*numberelem(index_Phos)) ! Na-phosphate 
        ion_excess_tot%K  = ion_excess_tot%K + (avfdisA(8)*numberelem(index_Phos)) ! K-phosphate
        
       ! Calculating little q 

       qnucl = ion_excess_tot%Na + ion_excess_tot%K - ion_excess_tot%Cl
       

       ! Calculating individual betas

       beta_ion_excess%Na =  ion_excess_tot%Na / abs(qnucl)
       beta_ion_excess%K =   ion_excess_tot%K  / abs(qnucl)
       beta_ion_excess%Cl = -ion_excess_tot%Cl / abs(qnucl)
    
    end subroutine make_beta_old
    


    ! Calculates the absolute value of the electostatic potential for each face of the lattice 
    ! Output assignment to  max_psi(6)  in mod parameters

    subroutine max_potential()

        use volume, only : indextocoord
        use volume, only : nx, ny, nz
        use globals, only : nsize
        use parameters, only : max_psi

        ! local variable

        real(dp) :: psi_back,psi_front,psi_left,psi_right,psi_up,psi_down
        logical:: is_back, is_front, is_left, is_right,is_up, is_down

        integer :: ind

        ! init
        psi_back = 0.0_dp
        psi_front = 0.0_dp
        psi_left = 0.0_dp
        psi_right = 0.0_dp
        psi_up = 0.0_dp
        psi_down = 0.0_dp
    
        do ind=1,nsize
            is_back =indextocoord(ind,1)==1
            is_front=indextocoord(ind,1)==nx 
            is_left =indextocoord(ind,2)==1
            is_right=indextocoord(ind,2)==ny 
            is_up   =indextocoord(ind,3)==1
            is_down =indextocoord(ind,3)==nz 

            if(is_back)  psi_back =max(psi_back ,abs(psi(ind)))
            if(is_front) psi_front=max(psi_front,abs(psi(ind)))
            if(is_left)  psi_left =max(psi_left ,abs(psi(ind)))
            if(is_right) psi_right=max(psi_right,abs(psi(ind)))
            if(is_up)    psi_up   =max(psi_up   ,abs(psi(ind)))
            if(is_down)  psi_down =max(psi_down ,abs(psi(ind)))
        
        enddo
        
        max_psi = (/psi_back, psi_front, psi_left, psi_right, psi_up, psi_down/) 
        
    end subroutine max_potential




      
end module field

