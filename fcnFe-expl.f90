! --------------------------------------------------------------|
! fcnFe-expl.f90:                                               |
! Constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for                |
! nucleosome with inbinding and distubed volume and             |
! mono, di adn tri valent  binding                              |
! using phosphate pairs.                                        |
! --------------------------------------------------------------|

module modfcnFeexpl
   
    use precision_definition

    implicit none

    real(dp) :: exponent_third
    logical, dimension(:,:,:), allocatable :: allowedstates

    integer :: singlet_list(7), doublet_list(3), triplet_list(1)
    character(len=8) :: char_state(11)
    real(dp), parameter  :: eps_val =1.0e-7_dp 
    
    private 
    public :: compute_fdisPPP, init_var_compute_fdisppp
    public :: fcnnucl_Fe_expl, test_compute_fdisPPP, compute_average_charge_PPP_expl, fcnnonucl_cp
    public :: compute_FEchem_react_PPP_expl

contains
   
    subroutine allocate_allowedstates(dim)

        integer, intent(in) :: dim

        integer :: ier

        allocate(allowedstates(dim,dim,dim),stat=ier)
        if( ier/=0 ) then
            print*, 'Allocation error : stat =', ier,' for allowstates"' 
            stop
        endif

    end subroutine allocate_allowedstates

    ! count number of phoshates in singlet, doublet and triplet state=(state1,state2,state3)
    ! input integer state1, state2, state3
    ! output integer numPhosSinglet,numPhosDoublet,numPhosTriplet)

    subroutine CountPhosphates(state1,state2,state3,numPhosSinglet,numPhosDoublet,numPhosTriplet)

        integer, intent(in) :: state1, state2, state3
        integer, intent(inout) :: numPhosSinglet, numPhosDoublet, numPhosTriplet

        integer :: i
  
        numPhosSinglet = 0
        numPhosDoublet = 0
        numPhosTriplet = 0
    
        do i=1,size(singlet_list)
            if(state1 == singlet_list(i)) numPhosSinglet = numPhosSinglet+1
            if(state2 == singlet_list(i)) numPhosSinglet = numPhosSinglet+1
            if(state3 == singlet_list(i)) numPhosSinglet = numPhosSinglet+1
        enddo

        do i=1,size(doublet_list)
            if(state1 == doublet_list(i)) numPhosDoublet = numPhosDoublet+1
            if(state2 == doublet_list(i)) numPhosDoublet = numPhosDoublet+1
            if(state3 == doublet_list(i)) numPhosDoublet = numPhosDoublet+1
        enddo
        
        do i=1,size(triplet_list)
            if(state1 == triplet_list(i)) numPhosTriplet = numPhosTriplet+1
            if(state2 == triplet_list(i)) numPhosTriplet = numPhosTriplet+1
            if(state3 == triplet_list(i)) numPhosTriplet = numPhosTriplet+1
        enddo

    end subroutine CountPhosphates

    ! Check if the two doublet state of (state1,state2,state3) are is the same doublet state 
    ! input integer :: state1, state2,state3 
    ! output logical :: isSame

    function IsSameDoubletState(state1,state2,state3) result(isSame)
        
        integer, intent(in) :: state1, state2, state3
        logical :: isSame

        integer :: i
        logical :: state1doublet, state2doublet ,state3doublet

        isSame = .false.
        state1doublet = .false.
        state2doublet = .false.
        state3doublet = .false.

        do i=1,size(doublet_list)
            if(state1 == doublet_list(i)) state1doublet = .true.
            if(state2 == doublet_list(i)) state2doublet = .true.
            if(state3 == doublet_list(i)) state3doublet = .true.
        enddo

        if(.not.state1doublet) isSame = (state2 == state3)
        if(.not.state2doublet) isSame = (state1 == state3)
        if(.not.state3doublet) isSame = (state1 == state2)

    end function IsSameDoubletState


    ! Check if (state1,state2,state3) are is physical allowed same chemical state 
    ! input integer :: state1, state2,state3 
    ! output logical :: isAllowed 

    function CheckAllowedState(state1,state2,state3)result(isAllowed)

        integer, intent(in) :: state1, state2, state3
        logical :: isAllowed

        integer :: numPhosSinglet, numPhosDoublet, numPhosTriplet

        call CountPhosphates(state1,state2,state3,numPhosSinglet,numPhosDoublet,numPhosTriplet)

        ! analyze combination numPhosSinglet, numPhosDoublet,numPhosTriplet 
        ! only combination (3,0,0), ( 1,2,0) and (0,0,3) are allowed or possible: conserves the number of phosphates
        ! not sufficient 
        ! if (1,2,0) the doublet states need to be the same chemical state 

        isAllowed=.false.

        if(numPhosSinglet == 3 .and. numPhosDoublet == 0 .and. numPhosTriplet == 0 ) isAllowed=.true.
        if(numPhosSinglet == 1 .and. numPhosDoublet == 2 .and. numPhosTriplet == 0 ) then 
            if(IsSameDoubletState(state1,state2,state3)) isAllowed=.true.
        endif    
        if(numPhosSinglet == 0 .and. numPhosDoublet == 0 .and. numPhosTriplet == 3 ) isAllowed=.true.
        
    end function CheckAllowedState

    subroutine init_allowedStates

        integer :: JJ, KK, LL

        allowedstates = .false.
        
        do JJ=1,11
            do KK=1,11
                do LL=1,11
                    if ( CheckAllowedState(JJ,KK,LL))  allowedstates(JJ,KK,LL) = .true.
                enddo
            enddo
        enddo

    end subroutine init_allowedstates


    subroutine print_allowedStates

        integer :: JJ, KK, LL 
        integer :: numstates

        numstates=0
        
        do JJ=1,11
            do KK=1,11
                do LL=1,11
                    print*, char_state(JJ),' ', char_state(KK),' ',char_state(LL),' ',allowedstates(JJ,KK,LL)
                    if(allowedstates(JJ,KK,LL)) numstates=numstates+1
                enddo
            enddo
        enddo

        print*,"numstates = ",numstates, "totalstates = ",11*11*11

    end subroutine print_allowedstates

    subroutine  init_charstate(char_state) 

        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, PhosFe2, PhosFe3
        use parameters, only : Phos2Mg, Phos2Fe2, Phos2Fe3, Phos3Fe3
    
        character(len=8), intent(inout) :: char_state(11)

        ! associate chemical state give by number with character
    
        char_state(Phos)= "Phos"
        char_state(PhosH)= "PhosH"
        char_state(PhosK)= "PhosK"
        char_state(PhosNa)= "PhosNa"
        char_state(PhosMg)= "PhosMg"
        char_state(PhosFe2)= "PhosFe2"
        char_state(PhosFe3)= "PhosFe3"
        char_state(Phos2Mg)= "Phos2Mg"
        char_state(Phos2Fe2)= "Phos2Fe2"
        char_state(Phos2Fe3)= "Phos2Fe3"
        char_state(Phos3Fe3)= "Phos3Fe3"

    end subroutine  init_charstate

    subroutine init_state_list 
    
        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, PhosFe2, PhosFe3
        use parameters, only : Phos2Mg, Phos2Fe2, Phos2Fe3, Phos3Fe3

        singlet_list = (/Phos, PhosH, PhosK, PhosNa, PhosMg, PhosFe2, PhosFe3/)
        doublet_list = (/Phos2Mg, Phos2Fe2, Phos2Fe3/)
        triplet_list = (/Phos3Fe3/)
    
    end subroutine init_state_list         


    subroutine init_var_compute_fdisPPP

        call allocate_allowedstates(11)
        call init_state_list 
        call init_charstate(char_state)
        call init_allowedStates

        exponent_third = 1.0_dp/3.0_dp
        
        !call print_allowedStates

    end subroutine init_var_compute_fdisPPP
    

    subroutine compute_fdisPPP(fdisPPP,position1,position2,position3)

        use field, only : xHplus, xNa, xK, xMg, xFe2, xFe3, xsol
        use parameters, only : vNa, vK, vMg, vFe2, vFe3, deltavAA
        use parameters, only : K0aAA
        use parameters, only : Phos,PhosH, PhosK, PhosNa, PhosMg, Phos2Mg, PhosFe2, PhosFe3
        use parameters, only : Phos2Mg,Phos2Fe2, Phos2Fe3,PhosFe2, Phos3Fe3

        real(dp), intent(inout), dimension(:,:,:) :: fdisPPP
        integer , intent(in) :: position1, position2,  position3

        real(dp) :: xP(11,3), fPPP, sumxP      
        integer :: i, k, coord(3)     
        integer  :: JJ, KK, LL

        ! .. executable statements 

        coord(1) = position1 ! position in lattice numbers
        coord(2) = position2
        coord(3) = position3


        do k=1,3
            i = coord(k)
            xP(Phos,k)   = 1.0_dp
            ! monovalent binding
            xP(PhosH,k)  = xHplus(i)      /(K0aAA(1)*(xsol(i)**deltavAA(1)))   ! (PH)/P-    : f(PH)P(i,j)/fPP(i,j)
            xP(PhosNa,k) = (xNa(i)/vNa)   /(K0aAA(2)*(xsol(i)**deltavAA(2)))   ! PNa/P-     : f(PNa)P(i,j)/fPP(i,j) 
            xP(PhosK,k)  = (xK(i)/vK)     /(K0aAA(7)*(xsol(i)**deltavAA(7)))   ! PK/P-      : f(PK)P(i,j)/fPP(i,j) 
            xP(PhosMg,k) = (xMg(i)/vMg)   /(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! PMgP+/PP2-
            xP(PhosFe2,k) = (xFe2(i)/vFe2)/(K0aAA(8)*(xsol(i)**deltavAA(8)))   ! PFe(2)+/PP2-
            xP(PhosFe3,k) = (xFe3(i)/vFe3)/(K0aAA(10)*(xsol(i)**deltavAA(10))) ! PFe(3)2+/PP2-
            ! divalent binding
            xP(Phos2Mg,k)  = sqrt( (xMg(i)/vMg)  / (K0aAA(6) *(xsol(i)**deltavAA(6))))    ! P2Mg/PP2- 
            xP(Phos2Fe2,k) = sqrt( (xFe2(i)/vFe2)/ (K0aAA(9) *(xsol(i)**deltavAA(9))))    ! P2Fe2/PP2-
            xP(Phos2Fe3,k) = sqrt( (xFe3(i)/vFe3)/ (K0aAA(11) * (xsol(i)**deltavAA(11)))) ! P2Fe3/PP2-
            ! trivalent 
            xP(Phos3Fe3,k) =( (xFe3(i)/vFe3)/ (K0aAA(12) * (xsol(i)**deltavAA(12))))**(exponent_third) ! P3Fe3/PPP
        enddo

       ! do i=1,11
       !     print*,i,' ',char_state(i),' ',xP(i,:)
       ! enddo    

        sumxP = 0.0_dp
       
        do JJ=1,11
            do KK=1,11
                do LL=1,11
                    if( allowedstates( JJ, KK, LL) ) then 
                        sumxP = sumxP + xP(JJ,1) * xP(KK,2) * xP(LL,3) 
                    endif    
                enddo
            enddo
        enddo

    
        fPPP = 1.0_dp/sumxP    ! fraction of phophate triplets that are all charged
          
        fdisPPP =0.0_dp        ! init 

        do JJ=1,11             ! fraction of phophate triplets that form different chemical states
            do KK=1,11
                do LL=1,11
                    if( allowedstates( JJ, KK, LL) ) then 
                        fdisPPP(JJ,KK,LL) = fPPP * xP(JJ,1) * xP(KK,2) * xP(LL,3)
                    endif    
                enddo
            enddo
        enddo    
           
    end subroutine compute_fdisPPP


    subroutine test_compute_fdisPPP
       
        use globals, only    : nsize
        use field, only      : xsol,xNa,xCl,xK,xHplus,xOHmin,xFe2,xFe3,xMg,psi 
        use parameters, only : expmu 
        use parameters, only : vNa,vK,vCl,vFe2,vFe3,vMg
        use parameters, only : zNa,zK,zCl,zFe2,zFe3,zMg
        use volume , only    : nx,ny,nz, linearIndexFromCoordinate
        use random , only    : seed, rands
        use myutils, only    : print_to_log, LogUnit, error_handler, lenText


        ! local parameter
        
        real(dp), parameter :: epssum = 1.0e-8_dp

        ! local variables 

        real(dp), dimension(:,:,:), allocatable :: fdisPPP
        integer :: i, k, l, ix(3), idx(3) , ndim(3)
        integer :: JJ, KK, LL, seed_old
        real(dp) :: sumfPPP,sumallfPPP, ran 
        integer :: info
        character(len=lenText) :: text

        info = 0

        allocate(fdisPPP(11,11,11))
        
        ndim(1)=nx
        ndim(2)=ny
        ndim(3)=nz

        ! make random xsol and psi
        seed_old=seed
            
        do l=1,5

            do i=1,nsize
                xsol(i)=rands(seed)
                psi(i)= rands(seed) -0.5_dp
            enddo
        
            do i=1,nsize                  ! ion densities
                xNa(i)      = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
                xK(i)       = expmu%K* (xsol(i)**vK) *exp(-psi(i)*zK)  ! K+ volume fraction
                xCl(i)      = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
                xHplus(i)   = expmu%Hplus*(xsol(i))  *exp(-psi(i))     ! H+  volume fraction
                xOHmin(i)   = expmu%OHmin*(xsol(i))  *exp(+psi(i))     ! OH- volume fraction
                xFe2(i)     = expmu%Fe2*(xsol(i)**vFe2)*exp(-psi(i)*zFe2) ! Fe++ volume fraction
                xFe3(i)     = expmu%Fe3*(xsol(i)**vFe3)*exp(-psi(i)*zFe3) ! Fe+++ volume fraction
                xMg(i)      = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg)    ! Mg++ volume fraction
            enddo

            do k=1,3
                ! make random position on lattice 
                do i=1,3
                    ran=rands(seed)
                    ix(i)=int(ran*ndim(i))+1
                enddo            
                call linearIndexFromCoordinate(ix(1),ix(2),ix(3),idx(k))
            !    print*,"ix=",ix
            enddo

            call compute_fdisPPP(fdisPPP,idx(1),idx(2),idx(3))

            sumfPPP=0.0_dp
            sumallfPPP=0.0_dp

            do JJ=1,11             ! fraction of phophate triplets that form different chemical states
                do KK=1,11
                    do LL=1,11 
                        sumallfPPP = sumallfPPP + fdisPPP(JJ,KK,LL) 
                        if( allowedstates( JJ, KK, LL) ) then 
                            sumfPPP = sumfPPP + fdisPPP(JJ,KK,LL) 
                        endif    
                    enddo
                enddo
            enddo  

            if((abs(sumfPPP-1.0_dp) >  epssum ).or.( abs(sumallfPPP-1.0_dp)> epssum)) info=1
        
        enddo

        ! reset of seed
        seed=seed_old    

        deallocate(fdisPPP)


        if(info/=0) then 
            text="warning test_compute_fdisPPP failed" 
            print*,text
            call print_to_log(LogUnit,text) 
            call error_handler(info,"find_phosphate_pairs_exclude_triplets")
        endif

    end subroutine test_compute_fdisPPP


    
    ! nucleosome of AA and dna polymers
    ! with ion charegeable group being on one acid (tA) with counterion binding etc 
    ! distribute volume of neighboring cells

    subroutine fcnnucl_Fe_expl(x,f,nn)

        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, neq, local_conf, LEFT, RIGHT, bctype
        use parameters, only : expmu, vsol
        use parameters, only : vNa,vK,vCl,vFe2,vFe3,vCa,vMg,vnucl,vPP,vO2,vPPP
        use parameters, only : zNa,zK,zCl,zFe2,zFe3,zCa,zMg,qPP,qPPP
        use parameters, only : K0aAA, K0a, K0aion, ta, iter
        use parameters, only : Phos, Phos2Mg, Phos2Fe2, Phos2Fe3 
        use volume, only     : volcell
        use chains, only     : indexconf, type_of_monomer, logweightchain, nelem, ismonomer_chargeable
        use chains, only     : type_of_charge, elem_charge, indexconfpair, nneigh
        use chains, only     : energychainLJ, no_overlapchain
        use chains, only     : indexconftriplet, ntriplet
        use field, only      : xsol,xNa,xCl,xK,xHplus,xOHmin,xFe2,xFe3,xMg,xCa,xO2,rhopol,rhoqpol,rhoq 
        use field, only      : psi,gdisA,gdisB,fdis, rhopol_charge
        use field, only      : fdisPP_loc, fdisPP_loc_swap, fdisP2Mg_loc, fdisP2Mg_loc_swap, rhoqphos
        use field, only      : fdisP2Fe2_loc, fdisP2Fe2_loc_swap, fdisP2Fe3_loc, fdisP2Fe3_loc_swap
        use field, only      : q, lnproshift, xpol=>xpol_t, xpol_tot=>xpol
        use field, only      : fdisPPP_loc1 , fdisPPP_loc1_swap , fdisPPP_loc2 , fdisPPP_loc2_swap
        use field, only      : fdisPPP_loc3 , fdisPPP_loc3_swap 
        use field, only      : numbers_pairs, numbers_triplets
        use vectornorm, only : L2norm, L2norm_sub, L2norm_f90
        use Poisson, only    : Poisson_Equation_bc
        use surface, only    : sigmaqSurfL, sigmaqSurfR, psiSurfL, psiSurfR, surface_charge

        use modfcnMgexpl, only : compute_fdisPP

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)                     ! local density nucleosome
        real(dp) :: local_xpol(nsize,nsegtypes)                       ! local volumer fraction nucleosome
        real(dp) :: local_rhopol_charge(nsize,nsegtypes)              ! local density nucleosome chargeable      
        real(dp) :: local_q                                           ! local normalization q 
        real(dp) :: local_rhoqphos(nsize)                             ! local charge dnisty of phosphates      
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize) 
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,c,s,t,jcharge,m, mm, sw        ! dummy indices
        integer  :: JJ, KK,LL
        real(dp) :: norm, normvol,normPE
        real(dp) :: rhopol0 
        real(dp) :: xA(3),xB(2),sgxA,sgxB                            ! disociation ariables 
        real(dp) :: locallnproshift(2)
        real(dp) :: deltavpolstateCl, deltavpolstateNa, deltavpolstateK, deltaxpol
        real(dp) :: sum_rhoqphos,sum_xphos

        real(dp) :: K0aPP   ! Kdis of P2Mg pair temporarily define 

        real(dp) :: numphos, numphos_comp, numpairs,numtriplets
        logical  :: testnumphos

        ! .. executable statements 

        ! print*,"K0aAA=",K0aAA
        K0aPP=K0aAA(6) ! P2Mg ???

        testnumphos =.true.

        n=nsize

        ! read out x 
        k=n
        do i=1,n                     
            xsol(i) = x(i)        ! volume fraction solvent
            psi(i)  = x(i+k)      ! potential
        enddo  

        !  .. assign global and local polymer density 
        do t=1,nsegtypes
            do i=1,n
                xpol(i,t)  = 0.0_dp 
                rhopol(i,t) = 0.0_dp 
                local_xpol(i,t) = 0.0_dp
                local_rhopol(i,t) = 0.0_dp
                local_rhopol_charge(i,t) = 0.0_dp
                rhopol_charge(i,t) = 0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol_tot(i) = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i)  = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)      = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
            xK(i)       = expmu%K* (xsol(i)**vK) *exp(-psi(i)*zK)  ! K+ volume fraction
            xCl(i)      = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)   = expmu%Hplus*(xsol(i))  *exp(-psi(i))     ! H+  volume fraction
            xOHmin(i)   = expmu%OHmin*(xsol(i))  *exp(+psi(i))     ! OH- volume fraction
            xFe2(i)     = expmu%Fe2*(xsol(i)**vFe2)*exp(-psi(i)*zFe2) ! Fe++ volume fraction
            xFe3(i)     = expmu%Fe3*(xsol(i)**vFe3)*exp(-psi(i)*zFe3) ! Fe+++ volume fraction
            xCa(i)      = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)      = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
            xO2(i)      = expmu%O2*(xsol(i)**vO2)                  ! O2 volume fraction
            lnexppivw(i) = log(xsol(i))/vsol                       ! auxilary variable  divide by vsol  !!
            local_rhoqphos(i) = 0.0_dp 
        enddo

        do t=1,nsegtypes   
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    if(type_of_charge(t)=="A") then  !  acid
                     
                        do i=1,n

                            xA(1) = xHplus(i)/(K0a(t)*xsol(i))           ! AH/A!
                            xA(2) = (xNa(i)/vNa)/(K0aion(t,2))!*xsol(i)) ! ANa/A- :xsol(i)**deltav = xsol(i)**0= 1 
                            xA(3) = (xK(i)/vK)/(K0aion(t,3))!*xsol(i))   ! AK/A-
                            sgxA = 1.0_dp+xA(1)+xA(2)+xA(3)  
                            gdisA(i,1,t) = 1.0_dp/sgxA                    ! A^- 
                            gdisA(i,2,t) = gdisA(i,1,t)*xA(1)             ! AH 
                            gdisA(i,3,t) = gdisA(i,1,t)*xA(2)             ! ANa 
                            gdisA(i,4,t) = gdisA(i,1,t)*xA(3)             ! AK
                       
                            fdis(i,t) = gdisA(i,1,t)
                            lnexppi(i,t) = psi(i) -log(gdisA(i,1,t))      ! auxilary variable palpha log(xsol)*(delta vpol+0) =0 
                        enddo

                    else !  base
                        do i=1,n
                            xB(1) = (K0a(t)*xsol(i))/xHplus(i)            ! B/BH+
                            xB(2) = (xCl(i)/vCl)/(K0aion(t,2))!*xsol(i))  ! BHCl/BH+
                            sgxB =  1.0_dp+xB(1)+xB(2)  
                            gdisB(i,1,t) = 1.0_dp/sgxB                    ! BH^+
                            gdisB(i,2,t) = gdisB(i,1,t)*xB(1)             ! B
                            gdisB(i,3,t) = gdisB(i,1,t)*xB(2)             ! BHCl     
                    
                            lnexppi(i,t) = -log(gdisB(i,2,t))             ! auxilary variable palpha lo 

                            fdis(i,t) = gdisB(i,2,t)  
                        enddo
        
                    endif  
                                
                else
                    ! t=ta : phosphate
                           
                    !do ind=1,len_index_phos ! loop over index of  location of phosphates
                     !i = index_phos(ind)  ! give the lattice location 
                    do i=1,nsize  
                        lnexppi(i,t) =  psi(i)!!   ! auxilary variable palpha
                       ! here used to be computation of fdisPP
                    enddo

                endif
            else  

                fdis(:,t)  = 0.0_dp
                lnexppi(:,t) = 0.0_dp

            endif   
        enddo   
    
        !  .. computation polymer density fraction      
 
        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=local_conf,local_conf                         ! loop over cuantas

            if(no_overlapchain(c)) then 

                lnpro = lnpro+logweightchain(c) - energychainLJ(c)
                
                do s=1,nseg                           ! loop over segments 
                    t=type_of_monomer(s)
                    if(t/=ta) then 
                        do j=1,nelem(s)               ! loop over elements of segment 
                            k = indexconf(s,c)%elem(j)
                            lnpro = lnpro +lnexppivw(k)*vnucl(j,t)   ! excluded-volume contribution      
                        enddo
                        if(ismonomer_chargeable(t)) then
                            jcharge=elem_charge(t)
                            k = indexconf(s,c)%elem(jcharge) 
                            lnpro = lnpro + lnexppi(k,t)  ! electrostatic, VdW and chemical contribution
                        endif
                    else 
                        ! phosphates 
                       
                        k = indexconf(s,c)%elem(1)

                        do jj=1,nneigh(s,c)           ! loop neighbors 

                            m = indexconfpair(s,c)%elem(jj)
                           
                            call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, fdisP2Fe2_loc, fdisP2Fe3_loc, k , m)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                          -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))    
                        enddo

                        do jj=1,ntriplet(s,c)
                        
                            m  = indexconftriplet(1,s,c)%elem(jj) ! ntriplet for monomer s
                            mm = indexconftriplet(2,s,c)%elem(jj)   

                            call  compute_fdisPPP(fdisPPP_loc1, k , m, mm)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta) + lnexppi(mm,ta)+ &
                                            (lnexppivw(k) + lnexppivw(m)  + lnexppivw(mm))*vnucl(1,ta) &
                                          -log(fdisPPP_loc1(Phos,Phos,Phos))  )/(3.0_dp*ntriplet(s,c))   

                            ! divide by 3 not 6 because need to permute m and mm but symmetric                        
                             
                        enddo    

                    endif        
                enddo
            endif         
        enddo

        locallnproshift(1)=lnpro
        locallnproshift(2)=1  ! rank  
    
        ! call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        ! call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        ! lnproshift=globallnproshift(1)
        lnproshift=locallnproshift(1)
         
        do c=local_conf,local_conf                            ! loop over cuantas

            if(no_overlapchain(c)) then 

                lnpro=logweightchain(c) - energychainLJ(c)
           
                do s=1,nseg                           ! loop over segments 
                    t=type_of_monomer(s)
                    if(t/=ta) then 
                        do j=1,nelem(s)               ! loop over elements of segment 
                            k = indexconf(s,c)%elem(j)
                            lnpro = lnpro +lnexppivw(k)*vnucl(j,t)   ! excluded-volume contribution        
                        enddo
                        if(ismonomer_chargeable(t)) then
                            jcharge=elem_charge(t)
                            k = indexconf(s,c)%elem(jcharge) 
                            lnpro = lnpro + lnexppi(k,t)  ! electrostatic, VdW and chemical contribution 
                        endif
                    else 
                        ! phosphates 
                        k = indexconf(s,c)%elem(1)

                        do jj=1,nneigh(s,c)           ! loop neighbors 

                            m = indexconfpair(s,c)%elem(jj)

                            call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc,fdisP2Fe2_loc, fdisP2Fe3_loc,  k , m)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                          -log(fdisPP_loc(Phos,Phos)))/(2.0_dp*nneigh(s,c))
                        enddo

                        do jj=1,ntriplet(s,c)

                            m  = indexconftriplet(1,s,c)%elem(jj) ! ntriplet for monomer s
                            mm = indexconftriplet(2,s,c)%elem(jj)   

                            call  compute_fdisPPP(fdisPPP_loc1, k , m, mm)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta) + lnexppi(mm,ta) + &
                                                    (lnexppivw(k) + lnexppivw(m)  + lnexppivw(mm) )*vnucl(1,ta) &
                                                -log(fdisPPP_loc1(Phos,Phos,Phos))  )/(3.0_dp*ntriplet(s,c)) 
                            ! divide by 3 not 6 because need to permute m and mm but symmetric                        
                             
                        enddo    

                    endif        
                enddo    
            

                pro = exp(lnpro-lnproshift)   
                local_q = local_q+pro
               
                do s=1,nseg
                    t=type_of_monomer(s)
                    if(t/=ta) then  ! not phosphates
                        do j=1,nelem(s)
                            k = indexconf(s,c)%elem(j) 
                            local_xpol(k,t)=local_xpol(k,t)+pro*vnucl(j,t)          ! unnormed polymer volume fraction
                        enddo
                        if(ismonomer_chargeable(t)) then
                            jcharge=elem_charge(t)
                            k = indexconf(s,c)%elem(jcharge) 
                            local_rhopol_charge(k,t)=local_rhopol_charge(k,t)+pro   ! unnormed density of charge center
                        endif
                    else
                        ! pair density of phosphates 
                        k = indexconf(s,c)%elem(1)
                        ! k_ind = inverse_index_phos(k) 

                        do j=1,nneigh(s,c)

                            m = indexconfpair(s,c)%elem(j)

                            call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, fdisP2Fe2_loc, fdisP2Fe3_loc, k , m)
                            call  compute_fdisPP(fdisPP_loc_swap, fdisP2Mg_loc_swap,& 
                                        fdisP2Fe2_loc_swap, fdisP2Fe3_loc_swap, m , k)    

                            ! % first part integral

                            sum_rhoqphos=0.0_dp
                            sum_xphos=0.0_dp 
                        
                            do JJ=1,7
                                do KK=1,7
                                    sum_rhoqphos = sum_rhoqphos+&
                                        (fdisPP_loc(JJ,KK)*qPP(JJ)+fdisPP_loc_swap(JJ,KK)*qPP(KK))/2.0_dp
                                    sum_xphos = sum_xphos   +&
                                        (fdisPP_loc(JJ,KK)*vPP(JJ)+fdisPP_loc_swap(JJ,KK)*vPP(KK))/2.0_dp
                                enddo
                            enddo
        
                            sum_xphos=sum_xphos+(fdisP2Mg_loc+fdisP2Mg_loc_swap)*vPP(Phos2Mg)/4.0_dp 
                            sum_xphos=sum_xphos+(fdisP2Fe2_loc+fdisP2Fe2_loc_swap)*vPP(Phos2Fe2)/4.0_dp 
                            sum_xphos=sum_xphos+(fdisP2Fe3_loc+fdisP2Fe3_loc_swap)*vPP(Phos2Fe3)/4.0_dp 

                            ! division 4.0_dp  because symmetry and  vPP(Phos2Mg)/2 is volume change per phosphate 

                            sum_rhoqphos = sum_rhoqphos+((fdisP2Fe3_loc+fdisP2Fe3_loc_swap)*qPP(Phos2Fe3))/4.0_dp
                           
                            ! division 4.0_dp  because symmetry and  qPP(Phos2Fe3)/2 is charge per phosphate of bridge

                            local_rhoqphos(k) = local_rhoqphos(k) + pro * sum_rhoqphos /(2.0_dp*nneigh(s,c)) ! nneigh could be zero  hence with in loop 
                            local_xpol(k,ta) = local_xpol(k,ta) + pro * sum_xphos /(2.0_dp*nneigh(s,c))

                            local_rhopol_charge(k,ta)=local_rhopol_charge(k,ta)+pro/(2.0_dp*nneigh(s,c))
                    
                            ! second integral contributes to location m of rhoqpos and xphol  xpol  
                         
                            sum_rhoqphos=0.0_dp
                            sum_xphos=0.0_dp 
                        
                            ! contributes to location k of rhoqpos and xol
                               
                            do JJ=1,7
                                do KK=1,7   
                                    sum_rhoqphos = sum_rhoqphos+&
                                        (fdisPP_loc_swap(JJ,KK)*qPP(JJ)+fdisPP_loc(JJ,KK)*qPP(KK))/2.0_dp

                                    sum_xphos = sum_xphos   +&
                                        (fdisPP_loc_swap(JJ,KK)*vPP(JJ)+fdisPP_loc(JJ,KK)*vPP(KK))/2.0_dp
                                enddo
                            enddo
        
                            sum_xphos=sum_xphos+(fdisP2Mg_loc_swap +fdisP2Mg_loc)*vPP(Phos2Mg)/4.0_dp
                            sum_xphos=sum_xphos+(fdisP2Fe2_loc_swap +fdisP2Fe2_loc)*vPP(Phos2Fe2)/4.0_dp
                            sum_xphos=sum_xphos+(fdisP2Fe3_loc_swap +fdisP2Fe3_loc)*vPP(Phos2Fe3)/4.0_dp

                            ! division 4.0_dp  because symmetry and  vPP(Phos2Mg)/2 is volume change per phosphate 
                            
                            sum_rhoqphos = sum_rhoqphos+((fdisP2Fe3_loc_swap+fdisP2Fe3_loc)*qPP(Phos2Fe3))/4.0_dp

                            ! division 4.0_dp  because symmetry and  qPP(Phos2Fe3)/2 is charge per phosphate of bridge                    

                            local_rhoqphos(m) = local_rhoqphos(m) + pro * sum_rhoqphos /(2.0_dp*nneigh(s,c)) ! nneigh could be zero  hence with in loop 
                            local_xpol(m,ta) = local_xpol(m,ta) + pro * sum_xphos /(2.0_dp*nneigh(s,c))

                            local_rhopol_charge(m,ta)=local_rhopol_charge(m,ta)+pro/(2.0_dp*nneigh(s,c))

                       
                        enddo 


                        
                        do j=1,ntriplet(s,c)

                            do sw=1,2 ! sum of potentail other element of s permute m and mm

                                if(sw==1) then 
                                    m  = indexconftriplet(1,s,c)%elem(j) ! ntriplet for monomer s
                                    mm = indexconftriplet(2,s,c)%elem(j)   
                                else 
                                    m  = indexconftriplet(2,s,c)%elem(j) ! ntriplet for monomer s
                                    mm = indexconftriplet(1,s,c)%elem(j)   
                                endif 

                                ! all permutations of (k, m, mm)

                                call  compute_fdisPPP(fdisPPP_loc1,      k , m, mm)
                                call  compute_fdisPPP(fdisPPP_loc1_swap, k , mm, m)
                                call  compute_fdisPPP(fdisPPP_loc2,      m ,  k, mm)
                                call  compute_fdisPPP(fdisPPP_loc2_swap, mm , k, m)
                                call  compute_fdisPPP(fdisPPP_loc3,      m, mm,  k)
                                call  compute_fdisPPP(fdisPPP_loc3_swap, mm ,m,  k)

                                ! loc1      =  (k , m , mm)
                                ! loc1_swap =  (k , mm, m )
                                ! loc2      =  (m , k , mm)
                                ! loc2_swap =  (mm, k , m )
                                ! loc3      =  (m , mm, k )
                                ! loc3_swap =  (mm ,m , k )
                            
                                ! permutations of integral contrubutions of fPPP and second coordinate permutation of triplet
                                
                                ! first integral contribution to position k

                                sum_rhoqphos=0.0_dp
                                sum_xphos=0.0_dp 
                            
                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11
                                            sum_rhoqphos = sum_rhoqphos + &
                                                ((fdisPPP_loc1(JJ,KK,LL) + fdisPPP_loc1_swap(JJ,KK,LL))*qPPP(JJ) + &
                                                (fdisPPP_loc2(JJ,KK,LL) + fdisPPP_loc2_swap(JJ,KK,LL))*qPPP(KK) + &
                                                (fdisPPP_loc3(JJ,KK,LL) + fdisPPP_loc3_swap(JJ,KK,LL))*qPPP(LL) )/6.0_dp
                                        
                                            sum_xphos = sum_xphos  + &
                                                ((fdisPPP_loc1(JJ,KK,LL) + fdisPPP_loc1_swap(JJ,KK,LL))*vPPP(JJ) + &
                                                (fdisPPP_loc2(JJ,KK,LL) + fdisPPP_loc2_swap(JJ,KK,LL))*vPPP(KK) + &
                                                (fdisPPP_loc3(JJ,KK,LL) + fdisPPP_loc3_swap(JJ,KK,LL))*vPPP(LL) )/6.0_dp
                                        enddo
                                    enddo            
                                enddo
                                
                                local_rhoqphos(k) = local_rhoqphos(k) + pro * sum_rhoqphos /(6.0_dp*ntriplet(s,c)) ! ntriplet can be zero  hence with in loop 
                                local_xpol(k,ta) = local_xpol(k,ta) + pro * sum_xphos /(6.0_dp*ntriplet(s,c))
                                local_rhopol_charge(k,ta) = local_rhopol_charge(k,ta)+pro/(6.0_dp*ntriplet(s,c))

                                ! second integral contibution to position m
                                
                                ! loc1      =  (k , m , mm)
                                ! loc1_swap =  (k , mm, m )
                                ! loc2      =  (m , k , mm)
                                ! loc2_swap =  (mm, k , m )
                                ! loc3      =  (m , mm, k )
                                ! loc3_swap =  (mm ,m , k )
                                
                                sum_rhoqphos=0.0_dp
                                sum_xphos=0.0_dp 
                            
                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11
                                            sum_rhoqphos = sum_rhoqphos + &
                                                ((fdisPPP_loc2(JJ,KK,LL)     + fdisPPP_loc3(JJ,KK,LL))*qPPP(JJ) + &
                                                (fdisPPP_loc1(JJ,KK,LL)      + fdisPPP_loc3_swap(JJ,KK,LL))*qPPP(KK) + &
                                                (fdisPPP_loc1_swap(JJ,KK,LL) + fdisPPP_loc2_swap(JJ,KK,LL))*qPPP(LL) )/6.0_dp
                                        
                                            sum_xphos = sum_xphos  + & 
                                                ((fdisPPP_loc2(JJ,KK,LL)     + fdisPPP_loc3(JJ,KK,LL))*vPPP(JJ) + &
                                                (fdisPPP_loc1(JJ,KK,LL)      + fdisPPP_loc3_swap(JJ,KK,LL))*vPPP(KK) + &
                                                (fdisPPP_loc1_swap(JJ,KK,LL) + fdisPPP_loc2_swap(JJ,KK,LL))*vPPP(LL) )/6.0_dp
                                        enddo
                                    enddo            
                                enddo
                                
                                local_rhoqphos(mm) = local_rhoqphos(mm) + pro * sum_rhoqphos /(6.0_dp*ntriplet(s,c)) ! ntriplet can be zero  hence with in loop 
                                local_xpol(mm,ta) = local_xpol(mm,ta) + pro * sum_xphos /(6.0_dp*ntriplet(s,c))
                                local_rhopol_charge(mm,ta) = local_rhopol_charge(mm,ta)+pro/(6.0_dp*ntriplet(s,c))

                                ! third integral contibution to position mm
                                
                                ! loc1      =  (k , m , mm)
                                ! loc1_swap =  (k , mm, m )
                                ! loc2      =  (m , k , mm)
                                ! loc2_swap =  (mm, k , m )
                                ! loc3      =  (m , mm, k )
                                ! loc3_swap =  (mm ,m , k )

                                sum_rhoqphos=0.0_dp
                                sum_xphos=0.0_dp 
                            
                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11
                                            sum_rhoqphos = sum_rhoqphos + &
                                                ((fdisPPP_loc2_swap(JJ,KK,LL) + fdisPPP_loc3_swap(JJ,KK,LL))*qPPP(JJ) + &
                                                (fdisPPP_loc1_swap(JJ,KK,LL) + fdisPPP_loc3(JJ,KK,LL))*qPPP(KK) + &
                                                (fdisPPP_loc1(JJ,KK,LL)      + fdisPPP_loc2(JJ,KK,LL))*qPPP(LL) )/6.0_dp
                                        
                                            sum_xphos = sum_xphos  + & 
                                                ((fdisPPP_loc2_swap(JJ,KK,LL) + fdisPPP_loc3_swap(JJ,KK,LL))*vPPP(JJ) + &
                                                (fdisPPP_loc1_swap(JJ,KK,LL) + fdisPPP_loc3(JJ,KK,LL))*vPPP(KK) + &
                                                (fdisPPP_loc1(JJ,KK,LL)      + fdisPPP_loc2(JJ,KK,LL))*vPPP(LL) )/6.0_dp
                                        enddo
                                    enddo            
                                enddo
                                
                                local_rhoqphos(mm) = local_rhoqphos(mm) + pro * sum_rhoqphos /(6.0_dp*ntriplet(s,c)) ! ntriplet can be zero  hence with in loop 
                                local_xpol(mm,ta) = local_xpol(mm,ta) + pro * sum_xphos /(6.0_dp*ntriplet(s,c))
                                local_rhopol_charge(mm,ta) = local_rhopol_charge(mm,ta)+pro/(6.0_dp*ntriplet(s,c))

                            enddo
        
                        enddo        
        
                    endif

                enddo
            
            endif 
                
        enddo ! cuantas loop
        
        q = 0.0_dp 
        q = local_q
    
        ! first graft point 
        do t=1,nsegtypes
            do i=1,nsize
                xpol(i,t)=local_xpol(i,t) ! polymer volume fraction density 
            enddo
            if(ismonomer_chargeable(t)) then
                do i=1,nsize
                    rhopol_charge(i,t)=local_rhopol_charge(i,t)   ! polymer density of charge center
                enddo    
            endif   
        enddo

        do i=1,nsize
            rhoqphos(i)=local_rhoqphos(i) 
        enddo     

        !  .. construction of fcn and volume fraction polymer 
        !  .. volume polymer segment per volume cell

        rhopol0=(1.0_dp/volcell)/q 

        do t=1, nsegtypes
            if(ismonomer_chargeable(t)) then 

                if(t/=ta) then
                    if(type_of_charge(t)=="A") then ! acid   

                        deltavpolstateNa=vNa*vsol
                        deltavpolstateK=vK*vsol                                     

                        do i=1,n
                            
                            rhopol_charge(i,t) = rhopol0 * rhopol_charge(i,t)                ! density nucleosome of type t  
                            rhoqpol(i) = rhoqpol(i) - gdisA(i,1,t)*rhopol_charge(i,t)*vsol   ! total charge density nucleosome in units of vsol 

                            ! volume fraction only consider Na and K ionpairing
                            deltaxpol = rhopol_charge(i,t)*(gdisA(i,3,t)*deltavpolstateNa+gdisA(i,4,t)*deltavpolstateK)
                            xpol(i,t) = rhopol0 * xpol(i,t) + deltaxpol                      ! scale xpol(i,t) and add delta xspol due to ionbinding
                        enddo

                    else  ! base   

                        deltavpolstateCl=vCl*vsol

                        do i=1,n
                            
                            rhopol_charge(i,t) = rhopol0 * rhopol_charge(i,t)                ! density nucleosome of type t chargeable 
                            rhoqpol(i) = rhoqpol(i) + gdisB(i,1,t)*rhopol_charge(i,t)*vsol   ! total charge density nucleosome in units of vsol 
                            
                            ! volume fraction only consider Cl ionpairing
                            deltaxpol = rhopol_charge(i,t)*gdisB(i,3,t)*deltavpolstateCl
                            xpol(i,t) = rhopol0 * xpol(i,t) + deltaxpol

                        enddo 
                        
                    endif     

                else
                    ! t=tAA phosphate 
                    
                    do i=1,n

                        rhopol_charge(i,ta) = rhopol0 * rhopol_charge(i,ta) 
                        rhoqphos(i) = rhopol0 * rhoqphos(i) 
                        rhoqpol(i) = rhoqpol(i) + rhoqphos(i)* vsol ! total  charge density in units of vsol 
                        xpol(i,ta) = rhopol0 * xpol(i,ta) 

                    enddo           
                        
                endif    
            else  

                ! volume fraction polymer of type t 
                do i=1,n
                    xpol(i,t)  = rhopol0 * xpol(i,t)   
                enddo

            endif 

            do i=1,n
                xpol_tot(i) = xpol_tot(i)+xpol(i,t)  
            enddo

        enddo    

        do i=1,n

            f(i) = xpol_tot(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xFe2(i)+xCa(i)+xMg(i)+&
                xK(i)+xFe3(i)+xO2(i)-1.0_dp
            rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zFe2*xFe2(i)/vFe2 +zFe3*xFe3(i)/vFe3+zK*xK(i)/vK ! total charge density in units of vsol  

        enddo
        
        ! .. end computation polymer density and charge density  

        ! .. electrostatics 
           
        sigmaqSurfR = surface_charge(bctype(RIGHT),psiSurfR,RIGHT)
        sigmaqSurfL = surface_charge(bctype(LEFT),psiSurfL,LEFT)
            
        ! .. Poisson Eq 
        
        !call Poisson_Equation(f,psi,rhoq)
        
        call Poisson_Equation_bc(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
    
        ! .. boundary conditions only if bctype /= cc or cp 
        !   call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bctype)    
        
        norm=l2norm_f90(f)
        iter=iter+1
                    
        normvol = L2norm_f90(f(1:nsize))
        normPE  = L2norm_f90(f(nsize+1:2*nsize))
        
        print*,'iter=', iter ,'norm=',norm, "normvol=",normvol,"normPE=",normPE
                    
        
        ! test
        if(testnumphos) then  
            numpairs= numbers_pairs()
            numtriplets = numbers_triplets()
            numphos = 2 * numpairs + 3 * numtriplets
            numphos_comp = sum(rhopol_charge(:,ta))* volcell ! computed number of phosphates 
            if(abs(numphos-numphos_comp)> eps_val) then 
                print*,"Computed number of phophates not equal to expected"
                print*,"numphos = ",numphos," numphos computed = ",numphos_comp       
            endif
     
        endif        

    end subroutine fcnnucl_Fe_expl    

    ! no nucleosome only  boundary cp 

    subroutine fcnnonucl_cp(x,f,nn)

        use precision_definition
        use globals, only    : nsize, neq, LEFT, RIGHT, bctype
        use parameters, only : expmu
        use parameters, only : vNa,vK,vCl,vFe2,vFe3,vCa,vMg,vO2
        use parameters, only : zNa,zK,zCl,zFe2,zFe3,zCa,zMg 
        use parameters, only : iter
        use field, only      : xsol,xNa,xCl,xK,xHplus,xOHmin,xFe2,xFe3,xMg,xCa,xO2,rhoq 
        use field, only      : psi
        use vectornorm, only : L2norm, L2norm_sub, L2norm_f90
        use Poisson, only    : Poisson_Equation_bc,  Poisson_Equation_Surface
        use surface, only    : sigmaqSurfL, sigmaqSurfR, psiSurfL, psiSurfR, surface_charge

        !     .. scalar arguments
        integer(8), intent(in) :: nn

        !     .. array arguments
        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        integer  :: n,i,k       ! dummy indices
        real(dp) :: norm, normvol,normPE
              
        ! .. executable statements 
        n = nsize
        
        ! read out x 
        k = n
        do i=1,n                     
            xsol(i) = x(i)        ! volume fraction solvent
            psi(i)  = x(i+k)      ! potential
        enddo  

        do i=1,n                  ! init volume fractions
       
            xNa(i)      = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
            xK(i)       = expmu%K* (xsol(i)**vK) *exp(-psi(i)*zK)  ! K+ volume fraction
            xCl(i)      = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)   = expmu%Hplus*(xsol(i))  *exp(-psi(i))     ! H+  volume fraction
            xOHmin(i)   = expmu%OHmin*(xsol(i))  *exp(+psi(i))     ! OH- volume fraction
            xFe2(i)     = expmu%Fe2*(xsol(i)**vFe2)*exp(-psi(i)*zFe2) ! Fe++ volume fraction
            xFe3(i)     = expmu%Fe3*(xsol(i)**vFe3)*exp(-psi(i)*zFe3) ! Fe+++ volume fraction
            xCa(i)      = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)      = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
            xO2(i)      = expmu%O2*(xsol(i)**vO2)                  ! O2 volume fraction
            
        enddo

        
        do i=1,n

            f(i) = xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xFe2(i)+xCa(i)+xMg(i)+&
                xK(i)+xFe3(i)+xO2(i)-1.0_dp
        
            rhoq(i) = zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zFe2*xFe2(i)/vFe2 +zFe3*xFe3(i)/vFe3+zK*xK(i)/vK ! total charge density in units of vsol  

        enddo
        

        ! .. electrostatics 
           
        sigmaqSurfR = surface_charge(bctype(RIGHT),psiSurfR,RIGHT)
        sigmaqSurfL = surface_charge(bctype(LEFT),psiSurfL,LEFT)
            
        ! .. Poisson Eq 
        
        !call Poisson_Equation(f,psi,rhoq)
        
        call Poisson_Equation_bc(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
    
        ! .. boundary conditions only if bctype /= cc or cp 
        call Poisson_Equation_Surface(f,psi,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bctype)    
        
        norm=l2norm_f90(f)
        iter=iter+1
                    
        normvol = L2norm_f90(f(1:nsize))
        normPE  = L2norm_f90(f(nsize+1:2*nsize))
        
        print*,'iter=', iter ,'norm=',norm, "normvol=",normvol,"normPE=",normPE
                    
        
    end subroutine fcnnonucl_cp   




    ! compute the average fraction of charged state of the phosphate pairs 

    subroutine compute_average_charge_PPP_expl(avfdisP2Mg,avfdisP2Fe2,avfdisP2Fe3,avfdisPP, avfdisPPP)

        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, local_conf, DEBUG
        use parameters, only : vsol, vnucl, ta, K0aAA, Phos 
        use volume, only     : nx, ny
        use chains, only     : indexconf, type_of_monomer, logweightchain, nelem, ismonomer_chargeable
        use chains, only     : type_of_charge, elem_charge, indexconfpair, nneigh
        use chains, only     : energychainLJ, no_overlapchain
        use field, only      : xsol, psi, fdis, fdisPP_loc, fdisPP_loc_swap 
        use field, only      : fdisP2Mg_loc, fdisP2Mg_loc_swap, fdisP2Fe2_loc, fdisP2Fe2_loc_swap
        use field, only      : fdisP2Fe3_loc, fdisP2Fe3_loc_swap
        use field, only      : fdisPPP_loc1 , fdisPPP_loc1_swap , fdisPPP_loc2 , fdisPPP_loc2_swap
        use field, only      : fdisPPP_loc3 , fdisPPP_loc3_swap 
        use field, only      : numbers_pairs, numbers_triplets
        use field, only      : q, lnproshift
        use myutils, only    : error_handler
        use modfcnMgexpl, only : compute_fdisPP
        use chains, only     : indexconftriplet, ntriplet

        real(dp), intent(inout) :: avfdisP2Mg, avfdisP2Fe2, avfdisP2Fe3
        real(dp), intent(inout) :: avfdisPP(7,7)
        real(dp), intent(inout) :: avfdisPPP(11,11,11)

        !     .. local variables
        
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,c,s,m,t,jcharge,mm,sw                ! dummy indices
        integer  :: JJ, KK, LL
        real(dp) :: local_avfdisP2Mg,local_avfdisPP(7,7),local_avfdisP2Fe2,local_avfdisP2Fe3
        real(dp) :: local_avfdisPPP(11,11,11)
        real(dp) :: sumrhopairs, sumrhotriplets
        real(dp) :: check
        real(dp) :: K0aPP   ! Kdis of P2Mg pair temporarily define 
        integer  :: nsizepsi
        ! .. executable statements 

        ! .. communication between processors 

        K0aPP=K0aAA(6) ! P2Mg
        nsizepsi = nsize + 2 * nx * ny


        local_avfdisPP = 0.0_dp
        local_avfdisP2Mg = 0.0_dp    
        local_avfdisP2Fe2 = 0.0_dp 
        local_avfdisP2Fe3 = 0.0_dp 

        local_avfdisPPP = 0.0_dp
        
        n=nsize

        do i=1,nsize
            lnexppivw(i)=log(xsol(i))/vsol
        enddo

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    if(type_of_charge(t)=="A") then  !  acid
                     
                        do i=1,n
                            lnexppi(i,t) = psi(i) -log(fdis(i,t))      ! auxilary variable palpha log(xsol)*(delta vpol+0) =0 
                        enddo

                    else !  base
                        do i=1,n
                            lnexppi(i,t) = -log(fdis(i,t))             ! auxilary variable palpha lo 
                        enddo
                    endif  
                                
                else
                    ! t=ta : phosphate
                    do i=1,n  
                        lnexppi(i,t) =  psi(i)!!   ! auxilary variable palpha
                    enddo

                endif
            else  
                lnexppi(:,t) = 0.0_dp
            endif   
        enddo   


        !  .. computation of probability 

        lnpro = 0.0_dp
              
        do c=local_conf,local_conf        ! loop over cuantas

            if( no_overlapchain(c)) then     
            
                lnpro=logweightchain(c) - energychainLJ(c)
               
                do s=1,nseg                       ! loop over segments 
                    t=type_of_monomer(s)
                    if(t/=ta) then 
                        do j=1,nelem(s)               ! loop over elements of segment 
                            k = indexconf(s,c)%elem(j)
                            lnpro = lnpro +lnexppivw(k)*vnucl(j,t)   ! excluded-volume contribution      
                        enddo
                        if(ismonomer_chargeable(t)) then
                            jcharge=elem_charge(t)
                            k = indexconf(s,c)%elem(jcharge) 
                            lnpro = lnpro + lnexppi(k,t)  ! electrostatic, VdW and chemical contribution
                        endif
                    else 
                        ! phosphates 
                        k = indexconf(s,c)%elem(1)

                        do jj=1,nneigh(s,c) ! loop neighbors 
                            m = indexconfpair(s,c)%elem(jj)

                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc,fdisP2Fe2_loc,fdisP2Fe3_loc, k ,m)
     
                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                    -log(fdisPP_loc(Phos,Phos)))/(2.0_dp*nneigh(s,c))
                        enddo    

                         do jj=1,ntriplet(s,c)
                        
                            m  = indexconftriplet(1,s,c)%elem(jj) ! ntriplet for monomer s
                            mm = indexconftriplet(2,s,c)%elem(jj)   

                            call  compute_fdisPPP(fdisPPP_loc1, k , m, mm)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta) + lnexppi(mm,ta)+ &
                                            (lnexppivw(k) + lnexppivw(m)  + lnexppivw(mm))*vnucl(1,ta) &
                                          -log(fdisPPP_loc1(Phos,Phos,Phos))  )/(3.0_dp*ntriplet(s,c))   

                            ! divide by 3 not 6 because need to permute m and mm but symmetric                        
                             
                        enddo    


                    endif        
                enddo    

                pro = exp(lnpro-lnproshift)   
            
                do s=1,nseg
                    
                    t=type_of_monomer(s)

                    if(t==ta) then 
                                    
                        ! pair density of phosphates 
                        k = indexconf(s,c)%elem(1)
       
                        do j=1,nneigh(s,c)

                            m = indexconfpair(s,c)%elem(j)

                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc,fdisP2Fe2_loc,fdisP2Fe3_loc, k ,m)
                            call compute_fdisPP(fdisPP_loc_swap,fdisP2Mg_loc_swap,&
                                        fdisP2Fe2_loc_swap, fdisP2Fe3_loc_swap,  m, k )

                            do JJ=1,7
                                do KK=1,7
                                    !local_avfdisPP(JJ,KK) = local_avfdisPP(JJ,KK)+&
                                    !    (fdisPP(k_ind,mr,JJ,KK)+fdisPP(m_ind,kr,JJ,KK))*pro/(2.0_dp*nneigh(s,c))

                                    local_avfdisPP(JJ,KK) = local_avfdisPP(JJ,KK)+&
                                        (fdisPP_loc(JJ,KK)+fdisPP_loc_swap(JJ,KK))*pro/(2.0_dp*nneigh(s,c))
                            
                                enddo
                            enddo
                            
                            local_avfdisP2Mg  = local_avfdisP2Mg  + fdisP2Mg_loc*pro/nneigh(s,c)
                            local_avfdisP2Fe2 = local_avfdisP2Fe2 + fdisP2Fe2_loc*pro/nneigh(s,c)
                            local_avfdisP2Fe3 = local_avfdisP2Fe3 + fdisP2Fe3_loc*pro/nneigh(s,c)
                    
                        enddo 

                        do j=1,ntriplet(s,c)

                            do sw=1,2 ! sum of potentail other element of s permute m and mm

                                if(sw==1) then 
                                    m  = indexconftriplet(1,s,c)%elem(j) ! ntriplet for monomer s
                                    mm = indexconftriplet(2,s,c)%elem(j)   
                                else 
                                    m  = indexconftriplet(2,s,c)%elem(j) ! ntriplet for monomer s
                                    mm = indexconftriplet(1,s,c)%elem(j)   
                                endif 

                                ! all permutations of (k, m, mm)

                                call  compute_fdisPPP(fdisPPP_loc1,      k , m, mm)
                                call  compute_fdisPPP(fdisPPP_loc1_swap, k , mm, m)
                                call  compute_fdisPPP(fdisPPP_loc2,      m ,  k, mm)
                                call  compute_fdisPPP(fdisPPP_loc2_swap, mm , k, m)
                                call  compute_fdisPPP(fdisPPP_loc3,      m, mm,  k)
                                call  compute_fdisPPP(fdisPPP_loc3_swap, mm ,m,  k)

                                ! loc1      =  (k , m , mm)
                                ! loc1_swap =  (k , mm, m )
                                ! loc2      =  (m , k , mm)
                                ! loc2_swap =  (mm, k , m )
                                ! loc3      =  (m , mm, k )
                                ! loc3_swap =  (mm ,m , k )
                            
                                ! permutations of integral contrubutions of fPPP and second coordinate permutation of triplet
                                
                                ! first integral contibution to position k

                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11
                                            local_avfdisPPP(JJ,KK,LL) = local_avfdisPPP(JJ,KK,LL)+ &
                                               ( fdisPPP_loc1(JJ,KK,LL) + fdisPPP_loc1_swap(JJ,KK,LL)) *&
                                                pro/(6.0_dp*ntriplet(s,c))
                                        enddo
                                    enddo            
                                enddo

                                ! second integral contibution to position m
                            
                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11
                                            local_avfdisPPP(JJ,KK,LL) = local_avfdisPPP(JJ,KK,LL)+ &
                                               (fdisPPP_loc2(JJ,KK,LL) + fdisPPP_loc3(JJ,KK,LL)) * &
                                                pro/(6.0_dp*ntriplet(s,c))
                                            
                                        enddo
                                    enddo            
                                enddo

                                ! third integral contibution to position mm 
                            
                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11
                                            local_avfdisPPP(JJ,KK,LL) = local_avfdisPPP(JJ,KK,LL)+ &
                                               (fdisPPP_loc2_swap(JJ,KK,LL) + fdisPPP_loc3_swap(JJ,KK,LL)) *&
                                                pro/(6.0_dp*ntriplet(s,c))
                                           
                                        enddo
                                    enddo            
                                enddo

                            enddo

                        enddo      
                    
                    endif
                enddo
            endif    
        enddo
   
        avfdisP2Mg = local_avfdisP2Mg/2.0_dp
        avfdisP2Fe2 = local_avfdisP2Fe2/2.0_dp
        avfdisP2Fe3 = local_avfdisP2Fe3/2.0_dp
        avfdisPP = local_avfdisPP/2.0_dp
        avfdisPPP = local_avfdisPPP/6.0_dp
 
        ! .. normalized avfdisPP with number of average number pairs 
        ! .. normalized avfdisPPP with number of average number triplets

        sumrhopairs = numbers_pairs()
        sumrhotriplets = numbers_triplets()

        avfdisPP    = avfdisPP/(sumrhopairs*q) ! also norm with q
        avfdisP2Mg  = avfdisP2Mg/(sumrhopairs*q)
        avfdisP2Fe2 = avfdisP2Fe2/(sumrhopairs*q)
        avfdisP2Fe3 = avfdisP2Fe3/(sumrhopairs*q)

        avfdisPPP  = avfdisPPP/(sumrhotriplets*q) ! also norm with q

    
        ! test 
        check = sum(avfdisPPP)
            
        if(abs(check-1.0_dp)> eps_val) then 
            print*,"sum avfdisPPP not equal to 1"
            print*,"sum avfdis fPPP = ",check       
        endif        

        print*,"sum fPPP =",check
            
    end subroutine compute_average_charge_PPP_expl


    ! compute the chemical free energy contribution for pairs and  triplets

    subroutine compute_FEchem_react_PPP_expl(FEChemPP, FEchemPPP)


        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, local_conf, DEBUG
        use parameters, only : vsol,vnucl
        use parameters, only : qPP, qPPP , vPP, vPPP , Phos, Phos2Mg,  Phos2Fe2, Phos2Fe3        
        use parameters, only : ta 
        use volume, only     : nx, ny
        use chains, only     : indexconf, type_of_monomer, logweightchain, nelem, ismonomer_chargeable
        use chains, only     : type_of_charge, elem_charge, indexconfpair, nneigh
        use chains, only     : energychainLJ, no_overlapchain
        use field, only      : xsol, psi, fdis
        use field, only      : fdisPP_loc, fdisP2Mg_loc, fdisP2Fe2_loc, fdisP2Fe3_loc
        use field, only      : fdisPPP_loc1, fdisPPP_loc1_swap, fdisPPP_loc2 , fdisPPP_loc2_swap
        use field, only      : fdisPPP_loc3, fdisPPP_loc3_swap 
        use field, only      : q, lnproshift
        use myutils, only    : error_handler
        use modfcnMgexpl, only : compute_fdisPP
        use chains, only     : indexconftriplet, ntriplet

        real(dp), intent(inout) :: FEchemPP, FEchemPPP


        !     .. local variables
        
        real(dp) :: lnexppi(nsize,nsegtypes)                         ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro, lnpro
        integer  :: n, i, j, k, c, s, m, t, jcharge, mm, sw                ! dummy indices
        integer  :: JJ, KK, LL
        integer  :: nsizepsi
        real(dp) :: lambda, sum_pi, sum_psi, psi_m, psi_k, psi_mm, betapi_k, betapi_m , betapi_mm
        real(dp) :: local_FEchempair, local_FEchemtriplet, FEchempair, FEchemtriplet
       
        ! .. executable statements 

        ! K0aPP=K0aAA(6) ! P2Mg
        nsizepsi=nsize + 2 * nx * ny

        Local_FEchempair = 0.0_dp
        local_FEchemtriplet = 0.0_dp

        n=nsize

        do i=1,nsize
            lnexppivw(i)=log(xsol(i))/vsol
        enddo

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    if(type_of_charge(t)=="A") then  !  acid
                     
                        do i=1,n
                            lnexppi(i,t) = psi(i) -log(fdis(i,t))      ! auxilary variable palpha log(xsol)*(delta vpol+0) =0 
                        enddo

                    else !  base
                        do i=1,n
                            lnexppi(i,t) = -log(fdis(i,t))             ! auxilary variable palpha lo 
                        enddo
                    endif  
                                
                else
                    ! t=ta : phosphate
                    do i=1,n  
                        lnexppi(i,t) =  psi(i)!!   ! auxilary variable palpha
                    enddo

                endif
            else  
                lnexppi(:,t) = 0.0_dp
            endif   
        enddo   


        !  .. computation of probability 

        lnpro = 0.0_dp
              
        do c=local_conf,local_conf        ! loop over cuantas

            if( no_overlapchain(c)) then     
            
                lnpro=logweightchain(c) - energychainLJ(c)
               
                do s=1,nseg                       ! loop over segments 
                    t=type_of_monomer(s)
                    if(t/=ta) then 
                        do j=1,nelem(s)               ! loop over elements of segment 
                            k = indexconf(s,c)%elem(j)
                            lnpro = lnpro +lnexppivw(k)*vnucl(j,t)   ! excluded-volume contribution      
                        enddo
                        if(ismonomer_chargeable(t)) then
                            jcharge=elem_charge(t)
                            k = indexconf(s,c)%elem(jcharge) 
                            lnpro = lnpro + lnexppi(k,t)  ! electrostatic, VdW and chemical contribution
                        endif
                    else 
                        ! phosphates 
                        k = indexconf(s,c)%elem(1)

                        do jj=1,nneigh(s,c) ! loop neighbors 
                            m = indexconfpair(s,c)%elem(jj)

                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc,fdisP2Fe2_loc,fdisP2Fe3_loc, k ,m)
     
                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                    -log(fdisPP_loc(Phos,Phos)))/(2.0_dp*nneigh(s,c))
                        enddo    

                         do jj=1,ntriplet(s,c)
                        
                            m  = indexconftriplet(1,s,c)%elem(jj) ! ntriplet for monomer s
                            mm = indexconftriplet(2,s,c)%elem(jj)   

                            call  compute_fdisPPP(fdisPPP_loc1, k , m, mm)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta) + lnexppi(mm,ta)+ &
                                            (lnexppivw(k) + lnexppivw(m)  + lnexppivw(mm))*vnucl(1,ta) &
                                          -log(fdisPPP_loc1(Phos,Phos,Phos))  )/(3.0_dp*ntriplet(s,c))   

                            ! divide by 3 not 6 because need to permute m and mm but symmetric                        
                             
                        enddo    


                    endif        
                enddo    

                pro = exp(lnpro-lnproshift)   
            
                do s=1,nseg
                    
                    t=type_of_monomer(s)

                    if(t==ta) then 
                                    
                        ! pair density of phosphates 
                        k = indexconf(s,c)%elem(1)
                      
                        betapi_k=-log(xsol(k))/vsol
                        psi_k = psi(k)
       
                        do j=1,nneigh(s,c)

                            m = indexconfpair(s,c)%elem(j)

                            betapi_m= -log(xsol(m))/vsol
                            psi_m = psi(m)
                        
                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc,fdisP2Fe2_loc, fdisP2Fe3_loc, k, m)

                             ! Lagrange multiplier lambd(r,r') 

                            lambda = -(betapi_k +betapi_m)*vPP(Phos) -(psi_k+psi_m)*qPP(Phos) &
                                -log(fdisPP_loc(Phos,Phos))

                            lambda = lambda*pro/nneigh(s,c)        

                            sum_pi  = 0.0_dp
                            sum_psi = 0.0_dp

                            do JJ=1,7
                                do KK=1,7
                                    sum_pi=sum_pi-(vPP(JJ)*betapi_k+vPP(KK)*betapi_m)*fdisPP_loc(JJ,KK)*pro/nneigh(s,c)
                                    sum_psi=sum_psi-(qPP(JJ)*psi_k+qPP(KK)*psi_m)*fdisPP_loc(JJ,KK)*pro/nneigh(s,c)
                                enddo
                            enddo
        
                            sum_pi=sum_pi-(vPP(Phos2Mg)/2.0_dp)* (betapi_k+betapi_m)*fdisP2Mg_loc*pro/nneigh(s,c)
                            sum_pi=sum_pi-(vPP(Phos2Fe2)/2.0_dp)*(betapi_k+betapi_m)*fdisP2Fe2_loc*pro/nneigh(s,c)
                            sum_pi=sum_pi-(vPP(Phos2Fe3)/2.0_dp)*(betapi_k+betapi_m)*fdisP2Fe3_loc*pro/nneigh(s,c)

                            ! division 2.0_dp  because  vPP(Phos2Mg)/2 is volume change per phosphate 

                            sum_psi=sum_psi-((((psi_k+psi_m)*qPP(Phos2Fe3))/(2.0_dp))*fdisP2Fe3_loc)*pro/nneigh(s,c)
                            ! division 2.0_dp  because  qPP(Phos2Fe3)/2 is charge  per phosphate 

                            local_FEchempair = local_FEchempair+(-lambda +sum_pi+sum_psi)/2.0_dp             
                       
                        enddo 


                        do j=1,ntriplet(s,c)

                            do sw=1,2 ! sum of potentail other element of s permute m and mm

                                if(sw==1) then 
                                    m  = indexconftriplet(1,s,c)%elem(j) ! ntriplet for monomer s
                                    mm = indexconftriplet(2,s,c)%elem(j)   
                                else 
                                    m  = indexconftriplet(2,s,c)%elem(j) ! ntriplet for monomer s
                                    mm = indexconftriplet(1,s,c)%elem(j)   
                                endif 

                                betapi_m = -log(xsol(m))/vsol
                                psi_m = psi(m)

                                betapi_mm = -log(xsol(m))/vsol
                                psi_mm  = psi(mm)

                                ! all permutations of (k, m, mm)

                                call  compute_fdisPPP(fdisPPP_loc1,      k , m, mm)
                                call  compute_fdisPPP(fdisPPP_loc1_swap, k , mm, m)
                                call  compute_fdisPPP(fdisPPP_loc2,      m ,  k, mm)
                                call  compute_fdisPPP(fdisPPP_loc2_swap, mm , k, m)
                                call  compute_fdisPPP(fdisPPP_loc3,      m, mm,  k)
                                call  compute_fdisPPP(fdisPPP_loc3_swap, mm ,m,  k)


                                ! Lagrange multiplier lambd(r,r',r'') 

                                lambda = -(betapi_k + betapi_m + betapi_mm)*vPP(Phos) & 
                                         -(psi_k    + psi_m + psi_mm  )*qPP(Phos) &
                                         -log(fdisPPP_loc1(Phos,Phos,Phos))

                                lambda = lambda*pro/(6.0_dp * ntriplet(s,c))      

                                sum_pi  = 0.0_dp
                                sum_psi = 0.0_dp
        

                                ! loc1      =  (k , m , mm)
                                ! loc1_swap =  (k , mm, m )
                                ! loc2      =  (m , k , mm)
                                ! loc2_swap =  (mm, k , m )
                                ! loc3      =  (m , mm, k )
                                ! loc3_swap =  (mm ,m , k )
                            
                                ! permutations of integral contributions of fPPP and second coordinate permutation of triplet
                                
                                ! first integral contibution to position k

                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11
                                        
                                            sum_pi = sum_pi- &
                                                ((vPPP(JJ) * betapi_k + vPPP(KK) * betapi_m + vPPP(LL) * betapi_mm ) * & 
                                                    fdisPPP_loc1(JJ,KK,LL) + &
                                                 (vPPP(JJ) * betapi_k + vPPP(KK) * betapi_mm + vPPP(LL) * betapi_m ) * & 
                                                    fdisPPP_loc1_swap(JJ,KK,LL) ) * pro/ (6.0_dp*ntriplet(s,c))


                                            sum_psi = sum_psi - &   
                                                ((qPPP(JJ) * psi_k + qPPP(KK) * psi_m + qPPP(LL) * psi_mm ) * & 
                                                    fdisPPP_loc1(JJ,KK,LL) + &
                                                 (qPPP(JJ) * psi_k + qPPP(KK) * psi_mm + qPPP(LL) * psi_m ) * & 
                                                    fdisPPP_loc1_swap(JJ,KK,LL) ) * pro/ (6.0_dp*ntriplet(s,c))
                                          

                                        enddo
                                    enddo            
                                enddo

                                ! second integral contibution to position m
                                ! loc2      =  (m , k , mm)
                                ! loc3      =  (m , mm, k )

                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11

                                            sum_pi = sum_pi- &
                                                ((vPPP(JJ) * betapi_m + vPPP(KK) * betapi_k + vPPP(LL) * betapi_mm ) * & 
                                                    fdisPPP_loc2(JJ,KK,LL) + &
                                                 (vPPP(JJ) * betapi_m + vPPP(KK) * betapi_mm + vPPP(LL) * betapi_k ) * & 
                                                    fdisPPP_loc3(JJ,KK,LL) ) * pro/ (6.0_dp*ntriplet(s,c))


                                            sum_psi = sum_psi - &   
                                                ((qPPP(JJ) * psi_m + qPPP(KK) * psi_k + qPPP(LL) * psi_mm ) * & 
                                                    fdisPPP_loc2(JJ,KK,LL) + &
                                                 (qPPP(JJ) * psi_m + qPPP(KK) * psi_mm + qPPP(LL) * psi_k ) * & 
                                                    fdisPPP_loc3(JJ,KK,LL) ) * pro/ (6.0_dp*ntriplet(s,c))

                                            
                                        enddo
                                    enddo            
                                enddo

                                ! third integral contibution to position mm 
                                ! loc2_swap =  (mm, k , m )
                                ! loc3_swap =  (mm ,m , k )
                            
                                do JJ=1,11
                                    do KK=1,11
                                        do LL=1,11

                                         sum_pi = sum_pi- &
                                                ((vPPP(JJ) * betapi_mm + vPPP(KK) * betapi_k + vPPP(LL) * betapi_m ) * & 
                                                    fdisPPP_loc2_swap(JJ,KK,LL) + &
                                                 (vPPP(JJ) * betapi_mm + vPPP(KK) * betapi_m + vPPP(LL) * betapi_k ) * & 
                                                    fdisPPP_loc3_swap(JJ,KK,LL) ) * pro/ (6.0_dp*ntriplet(s,c))


                                            sum_psi = sum_psi - &   
                                                ((qPPP(JJ) * psi_mm + qPPP(KK) * psi_k + qPPP(LL) * psi_m ) * & 
                                                    fdisPPP_loc2_swap(JJ,KK,LL) + &
                                                 (qPPP(JJ) * psi_mm + qPPP(KK) * psi_m + qPPP(LL) * psi_k ) * & 
                                                    fdisPPP_loc3_swap(JJ,KK,LL) ) * pro/ (6.0_dp*ntriplet(s,c))

                                        enddo
                                    enddo            
                                enddo

                                local_FEchemtriplet = local_FEchemtriplet+(-lambda +sum_pi+sum_psi)/6.0_dp        

                            enddo

                        enddo      
                    
                    endif
                enddo
            endif    
        enddo


        FEchempair = local_FEchempair
        FEchemtriplet = local_FEchemtriplet

        !  .. normalized FEchempair   with q 
       
        FEchempair = FEchempair/q 
        FEchemPP = FEchempair     
      
           
        !  .. normalized FEchemtripelt with q 
       
        FEchemtriplet= FEchemtriplet/q 
        FEchemPPP = FEchemtriplet     
      
    end subroutine compute_FEchem_react_PPP_expl

end module modfcnFeexpl

   
