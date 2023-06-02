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
    real(dp), dimension(:), allocatable :: xRb      ! volume fraction of positive Rb+ ion
    real(dp), dimension(:), allocatable :: xCa      ! volume fraction of positive Ca2+ ion
    real(dp), dimension(:), allocatable :: xMg      ! volume fraction of positive Mg2+ ion    
    real(dp), dimension(:), allocatable :: xNaCl    ! volume fraction of NaCl ion pair
    real(dp), dimension(:), allocatable :: xKCl     ! volume fraction of KCl  ion pair
    real(dp), dimension(:), allocatable :: xCl      ! volume fraction of negative ion
    real(dp), dimension(:), allocatable :: xHplus   ! volume fraction of Hplus
    real(dp), dimension(:), allocatable :: xOHmin   ! volume fraction of OHmin 

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

contains

    subroutine allocate_field(Nx,Ny,Nz,nsegtypes)
 
        integer, intent(in) :: Nx,Ny,Nz,nsegtypes
        
        integer :: N
        integer :: ier(26), i

        N=Nx*Ny*Nz

        allocate(xpol(N),stat=ier(1))
        allocate(xpol_t(N,nsegtypes),stat=ier(25))
        allocate(rhopol(N,nsegtypes),stat=ier(2)) 
        allocate(rhopolin(N,nsegtypes),stat=ier(3)) 
        allocate(rhopol_charge(N,nsegtypes),stat=ier(26)) 
        allocate(rhoqpol(N),stat=ier(4)) 
        allocate(xsol(N),stat=ier(5))
        allocate(psi(N+2*Nx*Ny),stat=ier(6))
        !allocate(psi(N),stat=ier(6))

        allocate(xNa(N),stat=ier(7))
        allocate(xK(N),stat=ier(8))
        allocate(xRb(N),stat=ier(9))
        allocate(xCa(N),stat=ier(10))
        allocate(xMg(N),stat=ier(11))
        allocate(xNaCl(N),stat=ier(12)) 
        allocate(xKCl(N),stat=ier(13)) 
        allocate(xCl(N),stat=ier(14)) 
        allocate(xHplus(N),stat=ier(15))
        allocate(xOHmin(N),stat=ier(16))
        allocate(rhoq(N),stat=ier(17))
        allocate(epsfcn(N),stat=ier(18))    
        allocate(Depsfcn(N),stat=ier(19))  

        allocate(fdis(N,nsegtypes),stat=ier(20))
        allocate(fdisA(N,8),stat=ier(21))
        allocate(fdisB(N,5),stat=ier(22))
        allocate(gdisA(N,4,nsegtypes),stat=ier(23))
        allocate(gdisB(N,3,nsegtypes),stat=ier(24))
       
        do i=1,26
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
        deallocate(xCa)
        deallocate(xNaCl) 
        deallocate(xKCl) 
        deallocate(xCl) 
        deallocate(xHplus)
        deallocate(xOHmin)
        deallocate(rhoq)
        deallocate(epsfcn)
        deallocate(Depsfcn)
        deallocate(fdis)
        deallocate(fdisA)
        deallocate(fdisB)
        deallocate(gdisA)
        deallocate(gdisB)
      
        
    end subroutine deallocate_field


    subroutine allocate_part_fnc(N)
       
        integer, intent(in) :: N

        ! allocate(lnq(N))
        ! allocate(q(N))

    end subroutine allocate_part_fnc

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
        xCa=0.0_dp
        xNaCl=0.0_dp 
        xKCl =0.0_dp
        xCl=0.0_dp
        xHplus=0.0_dp
        xOHmin=0.0_dp
        rhoq=0.0_dp
        psi=0.0_dp
        fdis=0.0_dp
        fdisA=0.0_dp
        fdisB=0.0_dp
        gdisA=0.0_dp
        gdisB=0.0_dp

    end subroutine init_field


    !  debug routine

    subroutine check_integral_rholpol_multi(sumrhopol, checkintegral)

        use volume, only : volcell
        use globals, only : nsize, systype, nseg, nsegtypes

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
        use globals, only : nsize, systype, nseg

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
        use parameters, only : zpol, qpol, qpol_tot, tA
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
        use parameters, only : zpol, qpol, qpol_tot, tA
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
        use chains, only : mapping_num_to_char

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
  
        ion_excess_ads%Na = ion_excess_ads%Na+ (avfdisA(3)*numberelem(index_Phos)) ! Na-phosphate 
        ion_excess_ads%K  = ion_excess_ads%K + (avfdisA(8)*numberelem(index_Phos)) ! K-phosphate
        
        ! calculate ion_excess_tot = sum of free adsorped ion excess
        
        ion_excess_tot%Na = ion_excess%Na  + ion_excess_ads%Na   
        ion_excess_tot%K  = ion_excess%K   + ion_excess_ads%K   
        ion_excess_tot%Cl = ion_excess%Cl  + ion_excess_ads%Cl   
        
        
        ! Calculate qnucl 
        qnucl = abs(ion_excess_tot%Na + ion_excess_tot%K - ion_excess_tot%Cl)
    
        ! Calculate individual betas

        beta_ion_excess%Na =  ion_excess_tot%Na / qnucl
        beta_ion_excess%K  =  ion_excess_tot%K  / qnucl
        beta_ion_excess%Cl = -ion_excess_tot%Cl / qnucl
        
            
    end subroutine make_beta



    subroutine make_ion_excess

        use parameters, only : vNa,vK,vMg,vCl,vCa
        use parameters, only : xbulk,ion_excess,sum_ion_excess

        ion_excess%Na=fcn_ion_excess(xNa,xbulk%Na,vNa)
        ion_excess%Cl=fcn_ion_excess(xCl,xbulk%Cl,vCl)
        ion_excess%K =fcn_ion_excess(xK,xbulk%K,vK)
        ion_excess%Mg=fcn_ion_excess(xMg,xbulk%Mg,vMg)
        ion_excess%Ca=fcn_ion_excess(xCa,xbulk%Ca,vCa)
        ion_excess%Hplus=fcn_ion_excess(xHplus,xbulk%Hplus,1.0_dp)
        ion_excess%OHmin=fcn_ion_excess(xOHmin,xbulk%OHmin,1.0_dp)
       
        ! sum of ion_excess weighted with valence of ion
       
        sum_ion_excess =ion_excess%Na -ion_excess%Cl+ion_excess%K +2.0_dp*ion_excess%Ca+2.0_dp*ion_excess%Mg +&
         ion_excess%Hplus -ion_excess%OHmin
        
    end subroutine make_ion_excess


    ! pre : make_ion_excess needed to be called before ion_excess
    ! pre : ion_excess needs to be computed. make_ion_excess called first
    ! input numberelem :  array containing the total elements of all types
    ! output assigns : beta_ion_excess  
    ! definitions beta_i = z_i * \Gamma_i / |qnucl|
  
    subroutine make_beta_old(numberelem)

        use parameters, only : ion_excess, beta_ion_excess, avgdisA,avgdisB,avfdisA
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
    
  
end module field

