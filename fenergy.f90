
!  .. module file for free energy variables /calculations

module energy 

    use precision_definition
    use molecules
    
    implicit none
  
    !     .. variables
    
    real(dp) :: FE                  ! free energy
    real(dp) :: FEbulk              ! free energybulk
    real(dp) :: deltaFE             ! free energy difference delteFE=FE-FEbulk
  
    !     .. auxiliary variable used in free energy computation  

    real(dp) :: FEq                 ! partition function poly A and B 
    real(dp) :: FEpi                ! sum over pi
    real(dp) :: FErho               ! sum over densities
    real(dp) :: FEel                ! electrostatics energy
    real(dp) :: FEelsurf(2)         ! electrostatics energy from  surface
    real(dp) :: FEelvar             ! electrostatics energy contrbution to total free energy from palpha due to varying dielectric 
    real(dp) :: FEelvarborn         ! electrostatics energy contrbution to total free energy from palpha due to Born self-energy 
    real(dp) :: FEborn              ! Born self-energy  
    real(dp) :: FEchemsurf(2)       ! chemical free energy surface
    real(dp) :: FEchem
    real(dp) :: FEchempair
    real(dp) :: FEbind,FEbindA,FEbindB    ! complexation contribution
    real(dp) :: FEVdW,FEVdWB,FEVdWC       ! Van der Waals contribution
    real(dp) :: FEconf
    real(dp) :: Econf
    real(dp) :: FEalt               ! free energy
    real(dp) :: FEbulkalt           ! free energybulk
    real(dp) :: deltaFEalt          ! free energy difference delteFE=FE-FEbulk
    real(dp) :: Eshift              ! shift in energy for palpha for neutralnoVdW

    real(dp) :: FEBornbulk          ! Born energy free energybulk

    real(dp) :: FEchemsurfalt(2)    ! chemical free energy surface
    real(dp) :: diffFEchemsurf(2)   ! difference chem

    type(moleclist) :: FEtrans,FEchempot,FEtransbulk,FEchempotbulk
    type(moleclist) :: deltaFEtrans,deltaFEchempot

    real(dp), dimension(:), allocatable :: sumphi, sumxpol, sumrhocharge ! integral over phi and xpol

    real(dp) :: sumphiA             ! integral over phiA
    real(dp) :: sumphiB             ! integral over phiB
    real(dp) :: qres                ! residual charge
    real(dp) :: checkphi            ! check integral over phi
    real(dp) :: checkrhocharge      ! check integral over rhocharge
    real(dp) :: checkxpol           ! check integral over xpol
    
    real(dp), parameter :: epsilon_sumxpol = 1.0e-10_dp     ! tolerance of surface coverage below no polymers 

    private :: epsilon_sumxpol

contains

    subroutine fcnenergy()
 
        use globals, only : systype 
        use myutils, only : print_to_log,LogUnit,lenText

      
        character(len=lenText) :: text 

        select case (systype) 
        case ("brush_mul","brush_mulnoVdW","brushdna","nucl_ionbin")
        
            call fcnenergy_electbrush_mul() 
            call fcnenergy_elect_alternative()

        case ("nucl_ionbin_sv")
        
            call fcnenergy_ionbin_sv()
            call fcnenergy_elect_alternative()

        case ("nucl_ionbin_Mg","nucl_ionbin_MgA")
        
            call fcnenergy_ionbin_sv()
            call fcnenergy_elect_alternative()


        case("elect")
            
            call fcnenergy_elect()
            call fcnenergy_elect_alternative()
        
        case("neutral","neutralnoVdW")
            
            call fcnenergy_neutral()
            call fcnenergy_neutral_alternative()  
        
        case ("brushborn")
            print*,"energy born not completed yet "   
            call fcnenergy_electbrush_mul() 
            call fcnenergy_elect_alternative()   

        case ("nucl_neutral_sv")  
            
            call fcnenergy_neutral_sv() 
            call fcnenergy_neutral_sv_alternative() 

        case default  

            text="fcnenergy: wrong systype: "//systype//"stopping program"
            call print_to_log(LogUnit,text)
            stop
        
        end select 
  
    end subroutine fcnenergy

   
    subroutine fcnenergy_elect()

        use globals
        use volume
        use parameters
        use field
        use surface

        !  .. local arguments 
    
        real(dp) :: sigmaq0,psi0
        real(dp) :: qsurf(2)           ! total charge on surface 
        real(dp) :: qsurfg             ! total charge on grafting surface  
        integer  :: i,j,s,g               ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer  :: nzadius
        real(dp) :: sigmaSurf(2),sigmaqSurf(2,ny*nx),sigmaq0Surf(2,nx*ny),psiSurf(2,nx*ny)
        real(dp) :: FEchemSurftmp
        integer, parameter :: A=1, B=2    

        !  .. computation of free energy 
    
        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
        sumphiA = 0.0_dp
        sumphiB = 0.0_dp
        FEq     = 0.0_dp
        FEbindA = 0.0_dp
        FEbindB = 0.0_dp
        FEchem = 0.0_dp
        FEVdWB = 0.0_dp     
        qres   = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - (xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xCl(i)/vCl+xK(i)/vK +&
                xNaCl(i)/vNaCl +xKCl(i)/vKCl)                 ! sum over  rho_i 
            FEel  = FEel  - rhoq(i) * psi(i)/2.0_dp      
            FEbindA = FEbindA + fdisA(i,5)*rhopol(i,A)
            FEbindB = FEbindB + fdisB(i,5)*rhopol(i,B)

            qres = qres + rhoq(i)
            sumphiA = sumphiA +  rhopol(i,A)
            sumphiB = sumphiB +  rhopol(i,B)
        enddo

        ! .. calcualtion of FEVdW
        FEVdW=0.0_dp
        

        FEel  = (volcell/vsol)*FEel
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
        FEbindA = volcell*FEbindA/2.0_dp !  
        FEbindB = volcell*FEbindB/2.0_dp !
        qres = (volcell/vsol)*qres
        sumphiA = volcell*sumphiA
        sumphiB = volcell*sumphiB
        checkphi= nseg-sumphiA-sumphiB

        FEbind = FEbindA+FEbindB
        
        ! .. calcualtion of FEq
        FEq=-log(q)    
        
         
        ! .. Shift in palpha  i.e q 
        Eshift=lnproshift! *ngr  

        !  .. total free energy  

        FE = FEq  + FEpi + FErho + FEel + FEVdW + FEbind -Eshift
        
    
        qres = qres + (qsurf(RIGHT)+qsurf(LEFT))  ! total residual charge 
  
        volumelat= volcell*nsize                  ! volume lattice

        FEbulk   = log(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ xbulk%Na/vNa +&
            xbulk%Ca/vCa +xbulk%Cl/vCl+ xbulk%K/vK + xbulk%NaCl/vNaCl +xbulk%KCl/vKCl )
        FEbulk = volumelat*FEbulk/(vsol)

        deltaFE = FE - FEbulk
    
    end subroutine fcnenergy_elect

    subroutine fcnenergy_elect_alternative()
    
        use globals
        use volume
        use parameters
        use field
        use VdW
        use surface
        use conform_entropy

        !  .. local arguments 
    
        real(dp) :: sigmaq0,psi0
        real(dp) :: qsurf(2)           ! total charge on surface 
        real(dp) :: qsurfg             ! total charge on grafting surface  
        integer :: i,j,s               ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer :: nzadius
      

        ! .. computation ofalternative computation free energy

        ! .. translational entropy 

        FEtrans%sol   = FEtrans_entropy(xsol,xbulk%sol,vsol,"w")   
        FEtrans%Na    = FEtrans_entropy(xNa,xbulk%Na,vNa)
        FEtrans%Cl    = FEtrans_entropy(xCl,xbulk%Cl,vCl)
        FEtrans%Ca    = FEtrans_entropy(xCa,xbulk%Ca,vCa)
        FEtrans%Mg    = FEtrans_entropy(xMg,xbulk%Mg,vMg)
        FEtrans%K     = FEtrans_entropy(xK,xbulk%K,vK)
        FEtrans%KCl   = FEtrans_entropy(xKCl,xbulk%KCl,vKCl)
        FEtrans%NaCl  = FEtrans_entropy(xNaCl,xbulk%NaCl,vNaCl)
        FEtrans%Hplus = FEtrans_entropy(xHplus,xbulk%Hplus,vsol,"w")
        FEtrans%OHmin = FEtrans_entropy(xOHmin,xbulk%OHmin,vsol,"w") 

        ! .. chemical potential + standard chemical potential 

        FEchempot%sol   = 0.0_dp ! by construction  
        FEchempot%Na    = FEchem_pot(xNa,expmu%Na,vNa)
        FEchempot%Cl    = FEchem_pot(xCl,expmu%Cl,vCl)
        FEchempot%Ca    = FEchem_pot(xCa,expmu%Ca,vCa)
        FEchempot%Mg    = FEchem_pot(xMg,expmu%Mg,vMg)
        FEchempot%K     = FEchem_pot(xK,expmu%K,vK) 
        FEchempot%KCl   = FEchem_pot(xKCl,expmu%KCl,vKCl)
        FEchempot%NaCl  = FEchem_pot(xNaCl,expmu%NaCl,vNaCl)
        FEchempot%Hplus = FEchem_pot(xHplus,expmu%Hplus,vsol,"w")
        FEchempot%OHmin = FEchem_pot(xOHmin,expmu%OHmin,vsol,"w")

        ! .. surface chemical contribution
 

        ! .. summing all contrubutions
        
        FEalt = FEtrans%sol +FEtrans%Na+ FEtrans%Cl +FEtrans%NaCl+FEtrans%Ca +FEtrans%Mg
        FEalt = FEalt+FEtrans%OHmin +FEtrans%Hplus +FEtrans%K +FEtrans%KCl 
        FEalt = FEalt+FEchempot%sol +FEchempot%Na+ FEchempot%Cl +FEchempot%NaCl+FEchempot%Ca +FEchempot%Mg
        FEalt = FEalt+FEchempot%OHmin +FEchempot%Hplus+ FEchempot%K +FEchempot%KCl
      

        ! be vary carefull FE = -1/2 \int dz rho_q(z) psi(z)

         ! .. chemical and binding contribution

        select case (systype) 
        case ("brush_mul","brush_mulnoVdW","brushdna","nucl_ionbin","nucl_ionbin_sv",&
            "brushborn","nucl_ionbin_Mg","nucl_ionbin_MgA")
            FEchem = FEchem_react_multi()
        case default
            FEchem = FEchem_react()
        end select 

        !  .. electrostatic Born self energy
        FEborn=FE_selfenergy_brush()

        if(isVdW) then 
            FEalt = FEalt-FEVdW ! add Van der Waals
        endif

        FEalt = FEalt - FEel + FEconf + Econf + FEchem +FEborn     
        !+ FEelSurf(RIGHT)+FEelSurf(LEFT)+FEchemSurfalt(RIGHT)+FEchemSurfalt(LEFT) 

        ! .. delta translational entropy

        FEtransbulk%sol   = FEtrans_entropy_bulk(xbulk%sol,vsol,"w")   
        FEtransbulk%Na    = FEtrans_entropy_bulk(xbulk%Na,vNa)
        FEtransbulk%Cl    = FEtrans_entropy_bulk(xbulk%Cl,vCl)
        FEtransbulk%Ca    = FEtrans_entropy_bulk(xbulk%Ca,vCa)
        FEtransbulk%Mg    = FEtrans_entropy_bulk(xbulk%Mg,vMg)
        FEtransbulk%K     = FEtrans_entropy_bulk(xbulk%K,vK)
        FEtransbulk%KCl   = FEtrans_entropy_bulk(xbulk%KCl,vKCl)
        FEtransbulk%NaCl  = FEtrans_entropy_bulk(xbulk%NaCl,vNaCl)
        FEtransbulk%Hplus = FEtrans_entropy_bulk(xbulk%Hplus,vsol,"w")
        FEtransbulk%OHmin = FEtrans_entropy_bulk(xbulk%OHmin,vsol,"w") 

        ! .. delta chemical potential + standard chemical potential 

        FEchempotbulk%sol   = 0.0_dp ! by construction  
        FEchempotbulk%Na    = FEchem_pot_bulk(xbulk%Na,expmu%Na,vNa)
        FEchempotbulk%Cl    = FEchem_pot_bulk(xbulk%Cl,expmu%Cl,vCl)
        FEchempotbulk%Ca    = FEchem_pot_bulk(xbulk%Ca,expmu%Ca,vCa)
        FEchempotbulk%Mg    = FEchem_pot_bulk(xbulk%Mg,expmu%Mg,vMg)
        FEchempotbulk%K     = FEchem_pot_bulk(xbulk%K,expmu%K,vK) 
        FEchempotbulk%KCl   = FEchem_pot_bulk(xbulk%KCl,expmu%KCl,vKCl)
        FEchempotbulk%NaCl  = FEchem_pot_bulk(xbulk%NaCl,expmu%NaCl,vNaCl)
        FEchempotbulk%Hplus = FEchem_pot_bulk(xbulk%Hplus,expmu%Hplus,vsol,"w")
        FEchempotbulk%OHmin = FEchem_pot_bulk(xbulk%OHmin,expmu%OHmin,vsol,"w")
        
        ! .. bulk free energy

        volumelat = volcell*nsize   ! volume lattice 
        FEbulkalt = FEtransbulk%sol +FEtransbulk%Na+ FEtransbulk%Cl +FEtransbulk%NaCl+FEtransbulk%Ca +FEtransbulk%Mg 
        FEbulkalt = FEbulkalt+FEtransbulk%OHmin +FEtransbulk%Hplus +FEtransbulk%K +FEtransbulk%KCl 
        FEbulkalt = FEbulkalt+FEchempotbulk%sol +FEchempotbulk%Na+FEchempotbulk%Cl +FEchempotbulk%NaCl+FEchempotbulk%Ca 
        FEbulkalt = FEbulkalt+FEchempotbulk%Mg +FEchempotbulk%OHmin +FEchempotbulk%Hplus +FEchempotbulk%K +FEchempotbulk%KCl
       

        FEBornbulk = (  bornbulk%Na*xbulk%Na/vNa    + bornbulk%Cl*xbulk%Cl/vCl + &
                        bornbulk%Ca*xbulk%Ca/vCa    + bornbulk%Mg*xbulk%Mg/vMg + &
                        bornbulk%Rb*xbulk%Rb/vRb    + bornbulk%K*xbulk%K/vK    + &
                        bornbulk%Hplus*xbulk%Hplus  + bornbulk%OHmin*xbulk%OHmin )/vsol  
        
        FEbulkalt = FEbulkalt+FEBornbulk

        FEbulkalt = volumelat*FEbulkalt

        FEBornbulk = volumelat* FEBornbulk

        ! .. delta

        deltaFEtrans%sol   = FEtrans%sol  - FEtransbulk%sol * volumelat
        deltaFEtrans%Na    = FEtrans%Na   - FEtransbulk%Na * volumelat
        deltaFEtrans%Cl    = FEtrans%Cl   - FEtransbulk%Cl * volumelat
        deltaFEtrans%Ca    = FEtrans%Ca   - FEtransbulk%Ca * volumelat
        deltaFEtrans%Mg    = FEtrans%Mg   - FEtransbulk%Mg * volumelat   
        deltaFEtrans%K     = FEtrans%K    - FEtransbulk%K * volumelat
        deltaFEtrans%KCl   = FEtrans%KCl  - FEtransbulk%KCl * volumelat
        deltaFEtrans%NaCl  = FEtrans%NaCl - FEtransbulk%NaCl * volumelat
        deltaFEtrans%Hplus = FEtrans%Hplus- FEtransbulk%Hplus * volumelat
        deltaFEtrans%OHmin = FEtrans%OHmin- FEtransbulk%OHmin * volumelat
         
        deltaFEchempot%sol   = FEchempot%sol  - FEchempotbulk%sol * volumelat
        deltaFEchempot%Na    = FEchempot%Na   - FEchempotbulk%Na * volumelat
        deltaFEchempot%Cl    = FEchempot%Cl   - FEchempotbulk%Cl * volumelat
        deltaFEchempot%Ca    = FEchempot%Ca   - FEchempotbulk%Ca * volumelat
        deltaFEchempot%Mg    = FEchempot%Mg   - FEchempotbulk%Mg * volumelat
        deltaFEchempot%K     = FEchempot%K    - FEchempotbulk%K * volumelat
        deltaFEchempot%KCl   = FEchempot%KCl  - FEchempotbulk%KCl * volumelat
        deltaFEchempot%NaCl  = FEchempot%NaCl - FEchempotbulk%NaCl * volumelat
        deltaFEchempot%Hplus = FEchempot%Hplus- FEchempotbulk%Hplus * volumelat
        deltaFEchempot%OHmin = FEchempot%OHmin- FEchempotbulk%OHmin * volumelat

        ! .. differences

        deltaFEalt = FEalt - FEbulkalt


    end subroutine fcnenergy_elect_alternative
    

    subroutine fcnenergy_electbrush_mul()

        use globals ! , only : systype, nsize,nsegtypes
        use volume
        use parameters
        use field
        use VdW
        use surface
        use dielectric_const, only : born
        use Poisson, only :  grad_pot_sqr_eps_cubic
        use chains, only : ismonomer_chargeable,type_of_monomer,type_of_charge

        !  .. local arguments 
    
        integer  :: i,j,s,g,t          ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        real(dp) :: FEchemSurftmp
        integer  :: ier
        logical  :: alloc_fail
        real(dp) :: sqrgradpsi(nsize)
        real(dp) :: Etotself,lbr
        real(dp) :: vnucltot
        integer  :: ncharge

        
        if (.not. allocated(sumphi))  then 
            allocate(sumphi(nsegtypes),stat=ier)
            if( ier/=0 ) alloc_fail=.true.
        endif
        if (.not. allocated(sumrhocharge))  then 
            allocate(sumrhocharge(nsegtypes),stat=ier)
            if( ier/=0 ) alloc_fail=.true.
        endif

        !  .. computation of free energy

        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
        FEelsurf = 0.0_dp
        sumphi = 0.0_dp
 
        FEq    = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp
        FEVdW  = 0.0_dp
        qres   = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - (xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xMg(i)/vMg+ xCl(i)/vCl+&
                xK(i)/vK +xNaCl(i)/vNaCl +xKCl(i)/vKCl )                 ! sum over  rho_i 
            FEel  = FEel  - rhoq(i) * psi(i)
            qres = qres + rhoq(i)
        enddo

        sumphi=calculate_sumphi()

        checkphi=nseg-sum(sumphi)
        

        if(systype=="nucl_ionbin_sv") then
            ncharge=0
            do s=1,nseg
                t=type_of_monomer(s)
                if(ismonomer_chargeable(t)) ncharge=ncharge+1 
            enddo
            checkrhocharge = ncharge 
            do t=1,nsegtypes
                sumrhocharge(t) = sum(rhopol_charge(:,t))
                sumrhocharge(t) = volcell*sumrhocharge(t)
                checkrhocharge = checkrhocharge -sumrhocharge(t)
            enddo
            checkrhocharge = ncharge -sum(sumrhocharge)
        endif 


        FEel  = (volcell/vsol)*FEel/2.0_dp  ! carefully check this
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
        qres  = (volcell/vsol)*qres


        ! .. calcualtion of FEVdW
        if(isVdW) then 
            FEVdW=-VdW_energy(rhopol)
        else   
            FEVdW=0.0_dp
        endif

        !  .. calcium and magnesium binding contribution to minimized free energy  

        if(systype/="brush_mul".and.systype/="brush_mulnoVdW") then 
            do i=1,nsize
                FEbind = FEbind + (fdisA(i,5)+fdisA(i,7))*rhopol(i,tA)
            enddo                
            FEbind = volcell*FEbind /2.0_dp                                 
        endif      

        if(systype=="brushborn") then

            ! .. scaled gradient potential contribution 
            call grad_pot_sqr_eps_cubic(psi,epsfcn, Depsfcn,sqrgradpsi)

            FEelvar = 0.0_dp           
            do i=1,nsize
                FEelvar=FEelvar + xpol(i)*sqrgradpsi(i)  
            enddo 
            FEelvar=FEelvar*volcell/vsol ! vsol divsion because of constqE definition in grad_pot_sqr_eps_cubic( 
            
            !  needs to  checked vocell  prefactor


            FEelvarborn=0.0_dp

            do i=1,nsize
                lbr = lb/epsfcn(i)     ! local Bjerrum length
                Etotself = &        ! total self energy  
                    born(lbr,bornrad%pol  ,zpolAA(1))*fdisA(i,1)*rhopol(i,tA) + &
                    born(lbr,bornrad%polCa,zpolAA(4))*fdisA(i,4)*rhopol(i,tA) + &
                    born(lbr,bornrad%polMg,zpolAA(6))*fdisA(i,6)*rhopol(i,tA) + &
                    born(lbr,bornrad%Na,zNa)*xNa(i)/(vNa*vsol)     + & 
                    born(lbr,bornrad%K,zK)*xK(i)/(vK*vsol)     + & 
                    born(lbr,bornrad%Cl,zCl)*xCl(i)/(vCl*vsol)     + &
                    born(lbr,bornrad%Rb,zRb)*xRb(i)/(vRb*vsol)     + & 
                    born(lbr,bornrad%Ca,zCa)*xCa(i)/(vCa*vsol)     + &
                    born(lbr,bornrad%Mg,zMg)*xMg(i)/(vMg*vsol)     + & 
                    born(lbr,bornrad%Hplus,1 )*xHplus(i)/vsol      + &
                    born(lbr,bornrad%OHmin,-1)*xOHmin(i)/vsol  

                FEelvarborn=FEelvarborn+Etotself*(Depsfcn(i)/epsfcn(i))*xpol(i)*volcell
            enddo     

        endif    

        ! .. calcualtion of FEq
     
        FEq=-log(q)    
    
     
        ! .. Shift in palpha  i.e q 
        Eshift=lnproshift !*ngr
       
        ! .. total free energy per area of surface 

        FE = FEq + FEpi + FErho + FEel + FEVdW + FEbind + FEelvar + FEelvarborn - Eshift
        
        
        volumelat= volcell*nsize   ! volume lattice

        FEbulk   = log(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ xbulk%Na/vNa +&
            xbulk%Ca/vCa +xbulk%Mg/vMg +xbulk%Cl/vCl+ xbulk%K/vK + xbulk%NaCl/vNaCl +xbulk%KCl/vKCl )
        
        FEbulk = volumelat*FEbulk/vsol

        deltaFE = FE - FEbulk
    
    end subroutine fcnenergy_electbrush_mul



    subroutine fcnenergy_ionbin_sv()

        use globals ! , only : systype, nsize,nsegtypes
        use volume
        use parameters
        use field
        use VdW
        use surface
        use dielectric_const, only : born
        use Poisson, only :  grad_pot_sqr_eps_cubic
        use chains, only : ismonomer_chargeable,type_of_monomer,type_of_charge

        !  .. local arguments 
    
        integer  :: i,j,s,g,t          ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        real(dp) :: FEchemSurftmp
        integer  :: ier
        logical  :: alloc_fail
        real(dp) :: sqrgradpsi(nsize)
        real(dp) :: Etotself,lbr
        real(dp) :: vnucltot
        integer  :: ncharge

        
        if (.not. allocated(sumphi))  then 
            allocate(sumphi(nsegtypes),stat=ier)
            if( ier/=0 ) alloc_fail=.true.
        endif
        if (.not. allocated(sumrhocharge))  then 
            allocate(sumrhocharge(nsegtypes),stat=ier)
            if( ier/=0 ) alloc_fail=.true.
        endif

        !  .. computation of free energy

        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
        FEelsurf = 0.0_dp
        sumphi = 0.0_dp
 
        FEq    = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp
        FEVdW  = 0.0_dp
        qres   = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - (xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xMg(i)/vMg+ xCl(i)/vCl+&
                xK(i)/vK +xNaCl(i)/vNaCl +xKCl(i)/vKCl )                 ! sum over  rho_i 
            FEel  = FEel  - rhoq(i) * psi(i)
            qres = qres + rhoq(i)
        enddo

        sumphi=calculate_sumphi()
        checkphi=nseg-sum(sumphi)

        checkxpol=check_volume_nucl()
   
        ! charge 
        ncharge=0
        do s=1,nseg
            t=type_of_monomer(s)
            if(ismonomer_chargeable(t)) ncharge=ncharge+1 
        enddo
        checkrhocharge = ncharge 
        do t=1,nsegtypes
            sumrhocharge(t) = sum(rhopol_charge(:,t))
            sumrhocharge(t) = volcell*sumrhocharge(t)
            checkrhocharge = checkrhocharge -sumrhocharge(t)
        enddo
        checkrhocharge = ncharge -sum(sumrhocharge)

        if(DEBUG) then
            print*,"     t       sumrhocharge(t) checkrhocharge  ismonomer_chargeable(t) type_of_charge(t)"
            do t=1,nsegtypes
                print*,t,sumrhocharge(t),checkrhocharge, ismonomer_chargeable(t),type_of_charge(t)
            enddo
        endif
            

        FEel  = (volcell/vsol)*FEel/2.0_dp  ! carefully check this
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
        qres  = (volcell/vsol)*qres


        ! .. calcualtion of FEVdW
        if(isVdW) then 
            FEVdW=-VdW_energy(rhopol)
        else   
            FEVdW=0.0_dp
        endif

        !  .. calcium and magnesium binding contribution to minimized free energy  

        if(systype/="brush_mul".and.systype/="brush_mulnoVdW") then 
            do i=1,nsize
                FEbind = FEbind + (fdisA(i,5)+fdisA(i,7))*rhopol(i,tA)
            enddo                
            FEbind = volcell*FEbind /2.0_dp                                 
        endif      

        ! .. calculaltion of FEq
     
        FEq=-log(q)    
    
     
        ! .. Shift in palpha  i.e q 
        Eshift=lnproshift 
       
        ! .. total free energy per area of surface 

        FE = FEq + FEpi + FErho + FEel + FEVdW + FEbind - Eshift
        
        
        volumelat= volcell*nsize   ! volume lattice

        FEbulk   = log(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ xbulk%Na/vNa +&
            xbulk%Ca/vCa +xbulk%Mg/vMg +xbulk%Cl/vCl+ xbulk%K/vK + xbulk%NaCl/vNaCl +xbulk%KCl/vKCl )
        
        FEbulk = volumelat*FEbulk/vsol

        deltaFE = FE - FEbulk
    
    end subroutine fcnenergy_ionbin_sv


    ! fcnenergy_neutral_sv  equal to fcnenergy_neutral !! 
    
    subroutine fcnenergy_neutral_sv()

        !  .. variable and constant declaractions 
    
        use globals, only : nsize, nseg, nsegtypes
        use volume, only : volcell
        use parameters
        use field, only : xsol, xpol_t, rhopol, q, lnproshift
        use VdW, only : VdW_energy
    

        !     .. local arguments 
    
        integer  :: i,j,t,g             ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer  :: ier
        logical  :: alloc_fail
        real(dp) :: vnucltot

        if (.not. allocated(sumphi))  then 
            allocate(sumphi(nsegtypes),stat=ier)
            if( ier/=0 ) alloc_fail=.true.
        endif

        ! .. computation of free energy 
    
        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEq = 0.0_dp
        FEVdW = 0.0_dp 
        FEel  = 0.0_dp
        FEelsurf = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp

        qres = 0.0_dp
        sumphi= 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - xsol(i) 
        enddo

        ! for systype=nucl_neutral_sv no rhopol used !!

        sumphi=calculate_sumphi()
        checkphi=nseg-sum(sumphi)
        print*,"checkphi=",checkphi
        checkxpol =check_volume_nucl()
        print*,"checkxpol=",checkxpol

        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
    
        ! .. calcualtion of FEVdW
        if(isVdW) then 
            FEVdW=-VdW_energy(rhopol)
        else   
            FEVdW=0.0_dp
        endif

        ! .. shift in palpha  i.e q 
        Eshift=lnproshift
    
        FEq=-log(q)    
       
        !  .. total free energy 

        FE = FEq + FEpi + FErho + FEVdW -Eshift
        
       
        volumelat= volcell*nsize  ! volume lattice divide by area surface
        FEbulk = log(xbulk%sol)-xbulk%sol
        FEbulk = volumelat*FEbulk/(vsol)
        deltaFE  = FE -FEbulk
    
    end subroutine fcnenergy_neutral_sv

    subroutine fcnenergy_neutral()

        !  .. variable and constant declaractions 
    
        use globals, only : nsize, nseg, nsegtypes
        use volume, only : volcell
        use parameters
        use field, only : xsol, rhopol, q, lnproshift
        use VdW, only : VdW_energy
    

        !     .. local arguments 
    
        integer  :: i,j,t,g             ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer  :: ier
        logical  :: alloc_fail

        if (.not. allocated(sumphi))  then 
            allocate(sumphi(nsegtypes),stat=ier)
            if( ier/=0 ) alloc_fail=.true.
        endif

        !     .. computation of free energy 
    
        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEq = 0.0_dp
        FEVdW = 0.0_dp 
        FEel  = 0.0_dp
        FEelsurf = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp

        qres = 0.0_dp
        sumphi= 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - xsol(i) 
        enddo

        checkphi=nseg
        do t=1,nsegtypes
            sumphi(t)=0.0_dp
            do i=1,nsize    
                sumphi(t) = sumphi(t) + rhopol(i,t)
            enddo
            sumphi(t) = volcell*sumphi(t)
            checkphi = checkphi-sumphi(t)
        enddo
    
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
    
        ! .. calcualtion of FEVdW
        if(isVdW) then 
            FEVdW=-VdW_energy(rhopol)
        else   
            FEVdW=0.0_dp
        endif

        ! .. shift in palpha  i.e q 
        Eshift=lnproshift
    
        FEq=-log(q)    
       
        !  .. total free energy 

        FE = FEq + FEpi + FErho + FEVdW -Eshift
        
       
        volumelat= volcell*nsize  ! volume lattice divide by area surface
        FEbulk = log(xbulk%sol)-xbulk%sol
        FEbulk = volumelat*FEbulk/(vsol)
        deltaFE  = FE -FEbulk
    
    end subroutine fcnenergy_neutral
      

    subroutine fcnenergy_neutral_alternative()

        !  .. variable and constant declaractions 
    
        use globals, only : nsize, nseg, nsegtypes
        use volume, only : volcell
        use parameters
        use field
        use VdW
        use conform_entropy
 
        real(dp) :: volumelat          ! volume lattice 
        
        ! .. translational entropy
        FEtrans%sol=FEtrans_entropy(xsol,xbulk%sol,vsol,"w") 
       
        ! .. chemical potential + standard chemical potential 
        FEchempot%sol   = 0.0_dp ! by construction 

    
        ! .. summing all contributions
        
        FEalt = FEtrans%sol + FEchempot%sol 
        FEalt = FEalt +FEconf -FEVdW+Econf

    
        ! .. delta translational entropy
        FEtransbulk%sol   = FEtrans_entropy_bulk(xbulk%sol,vsol,"w") 

        ! .. delta chemical potential + standard chemical potential 
        FEchempotbulk%sol   = 0.0_dp ! by construction 
        
        ! .. bulk free energy
        volumelat = volcell*nsize   ! volume lattice 
        FEbulkalt = FEtransbulk%sol +FEchempotbulk%sol 

        FEbulkalt = volumelat*FEbulkalt

        ! .. delta
        deltaFEtrans%sol = FEtrans%sol  - FEtransbulk%sol * volumelat
        deltaFEchempot%sol = FEchempot%sol  - FEchempotbulk%sol * volumelat

        ! .. differences

        deltaFEalt = FEalt - FEbulkalt

    
    end subroutine fcnenergy_neutral_alternative

    subroutine fcnenergy_neutral_sv_alternative()

        !  .. variable and constant declaractions 
    
        use globals, only : nsize, nseg, nsegtypes
        use volume, only : volcell
        use parameters
        use field
        use VdW
        use conform_entropy
 
        real(dp) :: volumelat          ! volume lattice 
        
        ! .. translational entropy
        FEtrans%sol=FEtrans_entropy(xsol,xbulk%sol,vsol,"w") 
       
        ! .. chemical potential + standard chemical potential 
        FEchempot%sol   = 0.0_dp ! by construction 

    
        ! .. summing all contributions
        
        FEalt = FEtrans%sol + FEchempot%sol + FEconf - FEVdW + Econf

    
        ! .. delta translational entropy
        FEtransbulk%sol = FEtrans_entropy_bulk(xbulk%sol,vsol,"w") 

        ! .. delta chemical potential + standard chemical potential 
        FEchempotbulk%sol = 0.0_dp ! by construction 
        
        ! .. bulk free energy
        volumelat = volcell*nsize   ! volume lattice 
        FEbulkalt = FEtransbulk%sol +FEchempotbulk%sol 

        FEbulkalt = volumelat*FEbulkalt

        ! .. delta
        deltaFEtrans%sol = FEtrans%sol  - FEtransbulk%sol * volumelat
        deltaFEchempot%sol = FEchempot%sol  - FEchempotbulk%sol * volumelat

        ! .. differences

        deltaFEalt = FEalt - FEbulkalt

    
    end subroutine fcnenergy_neutral_sv_alternative


    real(dp) function FEtrans_entropy(xvol,xvolbulk,vol,flag)
    
        use globals, only : nsize
        use parameters, only : vsol 
        use volume, only : volcell
        use field, only : xpol,xsol

        real(dp), intent(in) :: xvol(nsize)
        real(dp), intent(in) :: xvolbulk 
        real(dp), intent(in) :: vol    
        character(len=1), optional, intent(in) :: flag    

        integer :: i

        if(xvolbulk==0.0_dp) then 
            FEtrans_entropy = 0.0_dp
        else
            FEtrans_entropy = 0.0_dp
            if(present(flag)) then   ! water special case because vsol treated different then vi  
                do i=1,nsize
                    FEtrans_entropy = FEtrans_entropy + xvol(i)*(log(xvol(i))-1.0_dp)
                enddo 
                FEtrans_entropy = volcell*FEtrans_entropy/vol
            else 
                do i=1,nsize 
                    if(xvol(i)>0.0_dp) FEtrans_entropy = FEtrans_entropy + xvol(i)*(log(xvol(i)/vol)-1.0_dp) 
                enddo
                FEtrans_entropy = volcell*FEtrans_entropy/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy


    real(dp) function FEchem_pot(xvol,expchempot,vol,flag)
    
        use globals, only : nsize
        use volume, only : volcell
        use parameters, only : vsol
        implicit none

        real(dp), intent(in) :: xvol(nsize)
        real(dp), intent(in) :: expchempot 
        real(dp), intent(in) :: vol    
        character(len=1), optional, intent(in) :: flag    

        ! .. local 
        integer :: i
        real(dp) :: chempot ! chemical potential difference 
        real(dp) :: sumdens 

        if(expchempot==0.0_dp) then 
            FEchem_pot = 0.0_dp
        else     
            sumdens=0.0_dp
            if(present(flag)) then  ! water special case because vsol treated diffetent then vi
                chempot = -log(expchempot)    
                do i=1,nsize
                    sumdens = sumdens +xvol(i)
                enddo
                FEchem_pot = volcell*chempot*sumdens/vol            
            else
                chempot = -log(expchempot/vol)    
                do i=1,nsize
                    sumdens = sumdens +xvol(i)
                enddo
                FEchem_pot = volcell*chempot*sumdens/(vol*vsol)               
            endif
        endif

    end function FEchem_pot


    real(dp) function FEtrans_entropy_bulk(xvolbulk,vol,flag)
    
        use globals
        use parameters
        implicit none

        real(dp), intent(in) :: xvolbulk 
        real(dp), intent(in) :: vol    
        character(len=1), optional, intent(in) :: flag    

        integer :: i

        if(xvolbulk==0.0_dp) then 
            FEtrans_entropy_bulk=0.0_dp
        else
            FEtrans_entropy_bulk=0.0_dp
            if(present(flag)) then
                ! water special case because vsol treated diffetent then vi  
                FEtrans_entropy_bulk=xvolbulk*(log(xvolbulk)-1.0_dp)/vol
            else 
                FEtrans_entropy_bulk=xvolbulk*(log(xvolbulk/vol)-1.0_dp)/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy_bulk

    real(dp) function FEchem_pot_bulk(xvolbulk,expchempot,vol,flag)
    
        use globals
        use field
        use parameters
        implicit none

        real(dp), intent(in) :: xvolbulk
        real(dp), intent(in) :: expchempot 
        real(dp), intent(in) :: vol    
        character(len=1), optional, intent(in) :: flag    

        ! .. local 
        integer :: i
        real(dp) :: chempot ! chemical potential difference 
        real(dp) :: sumdens 

        if(expchempot==0.0_dp) then 
            FEchem_pot_bulk = 0.0_dp
        else     
            if(present(flag)) then  ! water special case because vsol treated diffetent then vi
                FEchem_pot_bulk = -log(expchempot)*xvolbulk/vol            
            else
                FEchem_pot_bulk = -log(expchempot/vol)*xvolbulk/(vol*vsol)               
            endif
        endif

    end function FEchem_pot_bulk


    real(dp) function FEchem_react()

        use globals, only : systype, nsize
        use field
        use volume, only : volcell
        use parameters, only : vpolA,vsol,zpolA,vpolB,zpolB


        integer :: i, k
        real(dp) :: lambdaA, lambdaB, rhopolAq, rhopolBq, xpolA, xpolB
        real(dp) :: betapi
        integer, parameter :: A=1, B=2


        select case (systype)
        case ("elect")

            FEchem_react = 0.0_dp

            do i=1,nsize

                betapi=-log(xsol(i))/vsol

                lambdaA=-log(fdisA(i,1)) -psi(i)*zpolA(1)-betapi*vpolA(1)*vsol
                lambdaB=-log(fdisB(i,1)) -psi(i)*zpolB(1)-betapi*vpolB(1)*vsol

                rhopolAq = 0.0_dp
                rhopolBq = 0.0_dp
                xpolA=0.0_dp
                xpolB=0.0_dp

                do k=1,4
                    rhopolAq=rhopolAq+ zpolA(k)*fdisA(i,k)*rhopol(i,A)
                    xpolA   =xpolA + rhopol(i,A)*fdisA(i,k)*vpolA(k)*vsol
                    rhopolBq=rhopolBq+ zpolB(k)*fdisB(i,k)*rhopol(i,B)
                    xpolB   =xpolB + rhopol(i,B)*fdisB(i,k)*vpolB(k)*vsol
                enddo   
                xpolA = xpolA +rhopol(i,A)*fdisA(i,5)*vpolA(5)*vsol/2.0_dp
                xpolB = xpolB +rhopol(i,B)*fdisB(i,5)*vpolB(5)*vsol/2.0_dp
                
                FEchem_react = FEchem_react + (- rhopol(i,A)*lambdaA -psi(i)*rhopolAq -betapi*xpolA &
                    +fdisA(i,5)*rhopol(i,A)/2.0_dp)

                FEchem_react = FEchem_react + (- rhopol(i,B)*lambdaB -psi(i)*rhopolBq -betapi*xpolB &
                    +fdisB(i,5)*rhopol(i,B)/2.0_dp)

            enddo

            FEchem_react=volcell*FEChem_react    

        case("electA") 

            FEchem_react = 0.0_dp

            do i=1,nsize

                betapi=-log(xsol(i))/vsol

                lambdaA=-log(fdisA(i,1)) -psi(i)*zpolA(1)-betapi*vpolA(1)*vsol

                rhopolAq = 0.0_dp
                xpolA=0.0_dp

                do k=1,4
                    rhopolAq=rhopolAq+ zpolA(k)*fdisA(i,k)*rhopol(i,A)
                    xpolA   =xpolA + rhopol(i,A)*fdisA(i,k)*vpolA(k)*vsol
                enddo   
                xpolA = xpolA +rhopol(i,A)*fdisA(i,5)*vpolA(5)*vsol/2.0_dp
                
                FEchem_react = FEchem_react + (- rhopol(i,A)*lambdaA -psi(i)*rhopolAq -betapi*xpolA &
                    +fdisA(i,5)*rhopol(i,A)/2.0_dp)
            enddo

            FEchem_react=volcell*FEChem_react 

        case("neutral","electnopoly") 
            FEchem_react=0.0_dp
        case default
            print*,"systype in FEchem_react wrong"  
        end select

    end function FEchem_react

    function FEchem_react_multi()result(FEchem_react)

        use globals, only : systype,nsize,nsegtypes
        use field, only : xsol,psi,rhopol,fdis,fdisA,gdisA,gdisB,epsfcn,Depsfcn
        use field, only : xNa,xK,xCl,xHplus,xOHmin,xCa,xMg,xRb,rhopol_charge
        use chains,only : elem_charge 
        use volume, only : volcell
        use parameters, only : vsol, vpol,vpolAA, zpol, zpolAA, lb, bornrad 
        use parameters, only : vNa,vCl,vRb,vCa, vK,vMg,zNa,zK,zCl,zRb,zCa,zMg, tA
        use parameters, only : vnucl
        ! use parameters, only : deltavnucl
        use chains, only : ismonomer_chargeable,type_of_charge
        use dielectric_const, only : born
        use Poisson, only :  grad_pot_sqr_eps_cubic
        ! use fcnaux, only : compute_xpol_chargeable, integral_betapi

        real(dp) :: FEchem_react

        integer  :: i, k, t, jcharge
        real(dp) :: lambda,  rhopolq, betapi, xpol, Eself ,bornene, lbr, Ebornself
        real(dp) :: xvol(nsize),sumbetapi(nsize)
        real(dp) :: sqrgradpsi(nsize) ! make allocatable only for sytype=brushborn
        integer  :: state 
        real(dp) :: vpolstate(4)

        FEchem_react = 0.0_dp

        select case (systype)
        case("brush_mul","brush_mulnoVdW") 
        
            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then
                    do i=1,nsize

                        betapi=-log(xsol(i))/vsol
                        lambda=-log(fdis(i,t)) -psi(i)*zpol(t,2)-betapi*vpol(t)*vsol
                        rhopolq=(zpol(t,1)*(1.0_dp-fdis(i,t))+ zpol(t,2)*fdis(i,t))*rhopol(i,t) ! charge from aciod or base 
            
                        FEchem_react = FEchem_react + &
                            (- rhopol(i,t)*lambda -psi(i)*rhopolq -betapi*rhopol(i,t)*vpol(t)*vsol )
                    enddo
                endif        
            enddo
        
        case("brush","brush_neq") 

            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then

                    do i=1,nsize

                        betapi=-log(xsol(i))/vsol
                        lambda=-log(fdisA(i,1)) -psi(i)*zpolAA(1)-betapi*vpolAA(1)*vsol
                       
            
                        rhopolq = 0.0_dp
                        xpol=0.0_dp
                        do k=1,4
                            rhopolq=rhopolq+ zpolAA(k)*fdisA(i,k)*rhopol(i,t)
                            xpol  =xpol + rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol
                        enddo   
                        rhopolq=rhopolq+ zpolAA(6)*fdisA(i,6)*rhopol(i,t)
                        xpol = xpol +rhopol(i,t)*(fdisA(i,5)*vpolAA(5)*vsol/2.0_dp +&
                                                  fdisA(i,6)*vpolAA(6)*vsol +&
                                                  fdisA(i,7)*vpolAA(7)*vsol/2.0_dp+&
                                                  fdisA(i,8)*vpolAA(8)*vsol )

            
                        FEchem_react = FEchem_react + &
                            ( -rhopol(i,t)*lambda -psi(i)*rhopolq -betapi*xpol+(fdisA(i,5)+fdisA(i,7))*rhopol(i,t)/2.0_dp)
                       
                    enddo

                endif        
            enddo


        case("brushdna") 
            
            do t=1,nsegtypes

                if(ismonomer_chargeable(t)) then

                    if(t==ta) then 
                        do i=1,nsize

                            betapi=-log(xsol(i))/vsol
                            lambda=-log(fdisA(i,1)) -psi(i)*zpolAA(1)-betapi*vpolAA(1)*vsol
                           
                
                            rhopolq = 0.0_dp
                            xpol=0.0_dp
                            do k=1,4
                                rhopolq=rhopolq+ zpolAA(k)*fdisA(i,k)*rhopol(i,t)
                                xpol  =xpol + rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol
                            enddo   
                            rhopolq=rhopolq+ zpolAA(6)*fdisA(i,6)*rhopol(i,t)
                            xpol = xpol +rhopol(i,t)*(fdisA(i,5)*vpolAA(5)*vsol/2.0_dp +&
                                                      fdisA(i,6)*vpolAA(6)*vsol +&
                                                      fdisA(i,7)*vpolAA(7)*vsol/2.0_dp+ &
                                                      fdisA(i,8)*vpolAA(8)*vsol)

                
                            FEchem_react = FEchem_react + ( -rhopol(i,t)*lambda &
                                -psi(i)*rhopolq -betapi*xpol+(fdisA(i,5)+fdisA(i,7))*rhopol(i,t)/2.0_dp)
                           
                        enddo

                    else
                        
                        do i=1,nsize

                            betapi = -log(xsol(i))/vsol
                            lambda = -log(fdis(i,t)) -psi(i)*zpol(t,2)-betapi*vpol(t)*vsol
                            rhopolq = (zpol(t,1)*(1.0_dp-fdis(i,t))+ zpol(t,2)*fdis(i,t))*rhopol(i,t)

            
                            FEchem_react = FEchem_react + &
                                (- rhopol(i,t)*lambda -psi(i)*rhopolq -betapi*rhopol(i,t)*vpol(t)*vsol )
                        enddo
                    endif
                        
                endif        
            enddo

        case("nucl_ionbin")
            
            do t=1,nsegtypes

                if(ismonomer_chargeable(t)) then

                    if(t==ta) then 
                        do i=1,nsize

                            betapi=-log(xsol(i))/vsol
                            lambda=-log(fdisA(i,1)) -psi(i)*zpolAA(1)-betapi*vpolAA(1)*vsol
                    
                            rhopolq = 0.0_dp
                            xpol=0.0_dp
                            do k=1,4
                                rhopolq=rhopolq+ zpolAA(k)*fdisA(i,k)*rhopol(i,t)
                                xpol  =xpol + rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol
                            enddo   
                            rhopolq=rhopolq+ zpolAA(6)*fdisA(i,6)*rhopol(i,t)
                            xpol = xpol +rhopol(i,t)*(fdisA(i,5)*vpolAA(5)*vsol/2.0_dp +&
                                                      fdisA(i,6)*vpolAA(6)*vsol +&
                                                      fdisA(i,7)*vpolAA(7)*vsol/2.0_dp+ &
                                                      fdisA(i,8)*vpolAA(8)*vsol)

                
                            FEchem_react = FEchem_react + ( -rhopol(i,t)*lambda &
                                -psi(i)*rhopolq -betapi*xpol+(fdisA(i,5)+fdisA(i,7))*rhopol(i,t)/2.0_dp)
                           
                        enddo

                    else
                        !if(zpol(t,1)==0.and.zpol(t,2)==-1) then ! acid 
                       
                        if(type_of_charge(t)=="A") then ! acid 
                            
                            vpolstate(1)=vpol(t)
                            vpolstate(2)=vpol(t)
                            vpolstate(3)=vpol(t)+vNa
                            vpolstate(4)=vpol(t)+vK
                            
                            do i=1,nsize

                                betapi = -log(xsol(i))/vsol
                                lambda = -log(gdisA(i,1,t)) + psi(i) -betapi*vpol(t)*vsol
                                rhopolq= -gdisA(i,1,t)*rhopol(i,t)
                                
                                xpol=0.0_dp
                                do state=1,4
                                    xpol = xpol + rhopol(i,t)*gdisA(i,state,t)*vpolstate(state)
                                enddo  
                                xpol=xpol*vsol

            
                                FEchem_react = FEchem_react + &
                                (- rhopol(i,t)*lambda -psi(i)*rhopolq -betapi*xpol )
                            enddo
                            
                        else if(type_of_charge(t)=="B") then ! base

                            vpolstate(1)=vpol(t)
                            vpolstate(2)=vpol(t)
                            vpolstate(3)=vpol(t)+vCl
                    
                            do i=1,nsize

                                betapi = -log(xsol(i))/vsol
                                lambda = -log(gdisB(i,2,t))  - betapi*vpol(t)*vsol
                                rhopolq= gdisB(i,1,t)*rhopol(i,t)

                                xpol=0.0_dp
                                do state=1,3
                                    xpol = xpol + rhopol(i,t)*gdisB(i,state,t)*vpolstate(state)
                                enddo  
                                xpol=xpol*vsol
        

                                FEchem_react = FEchem_react + &
                                (- rhopol(i,t)*lambda -psi(i)*rhopolq -betapi*xpol)
                            enddo
                        else 
                            print*,"Error in FEchem_react_multi for systype=nucl_ionbin"
                            print*,"Charged monomer ",t," not acid or base : zpol(1)=",zpol(t,1)," and zpol(2)=",zpol(t,2)
                        endif    

                    endif
                        
                endif        
            enddo

        case("nucl_ionbin_sv")
            
            do t=1,nsegtypes

                if(ismonomer_chargeable(t)) then

                    if(t==ta) then 
                        do i=1,nsize

                            betapi=-log(xsol(i))/vsol
                            lambda=-log(fdisA(i,1)) -psi(i)*zpolAA(1)-betapi*vpolAA(1)*vsol
                    
                            rhopolq = 0.0_dp
                            xpol=0.0_dp
                            do k=1,4
                                rhopolq=rhopolq+ zpolAA(k)*fdisA(i,k)*rhopol_charge(i,t)
                                xpol  =xpol + rhopol_charge(i,t)*fdisA(i,k)*vpolAA(k)*vsol
                            enddo   
                            rhopolq=rhopolq+ zpolAA(6)*fdisA(i,6)*rhopol(i,t)

                            xpol = xpol +rhopol_charge(i,t)*(fdisA(i,5)*vpolAA(5)*vsol/2.0_dp +&
                                                      fdisA(i,6)*vpolAA(6)*vsol +&
                                                      fdisA(i,7)*vpolAA(7)*vsol/2.0_dp+ &
                                                      fdisA(i,8)*vpolAA(8)*vsol)

                
                            FEchem_react = FEchem_react + ( -rhopol_charge(i,t)*lambda &
                                -psi(i)*rhopolq -betapi*xpol+(fdisA(i,5)+fdisA(i,7))*rhopol_charge(i,t)/2.0_dp)
                           
                        enddo

                    else
                        !if(zpol(t,1)==0.and.zpol(t,2)==-1) then ! acid 
                       
                        jcharge=elem_charge(t)

                        if(type_of_charge(t)=="A") then ! acid 
                            
                            vpolstate(1)=vnucl(jcharge,t)
                            vpolstate(2)=vnucl(jcharge,t)
                            vpolstate(3)=vnucl(jcharge,t)+vNa*vsol
                            vpolstate(4)=vnucl(jcharge,t)+vK*vsol

                            
                            do i=1,nsize

                                betapi = -log(xsol(i))/vsol
                             !   lambda = -log(gdisA(i,1,t)) + psi(i) -betapi*vpol(t)*vsol
                                lambda = -log(gdisA(i,1,t)) + psi(i) -betapi*vpolstate(1)
                                rhopolq= -gdisA(i,1,t)*rhopol_charge(i,t)
                                
                                xpol=0.0_dp
                                do state=1,4
                                    xpol = xpol + rhopol_charge(i,t)*gdisA(i,state,t)*vpolstate(state)
                                enddo  
                                !xpol=xpol*vsol

            
                                FEchem_react = FEchem_react + &
                                (- rhopol_charge(i,t)*lambda -psi(i)*rhopolq -betapi*xpol )
                            enddo
                           
                        else if(type_of_charge(t)=="B") then ! base

                            jcharge=elem_charge(t)
                            vpolstate(1)=vnucl(jcharge,t)
                            vpolstate(2)=vnucl(jcharge,t)
                            vpolstate(3)=vnucl(jcharge,t)+vCl*vsol
                    
                            do i=1,nsize

                                betapi = -log(xsol(i))/vsol
                                !lambda = -log(gdisB(i,2,t))  - betapi*vpol(t)*vsol
                                lambda = -log(gdisB(i,2,t))  - betapi*vpolstate(2)
                                rhopolq= gdisB(i,1,t)*rhopol_charge(i,t)

                                xpol=0.0_dp
                                do state=1,3
                                    xpol = xpol + rhopol_charge(i,t)*gdisB(i,state,t)*vpolstate(state)
                                enddo  
                                !xpol=xpol*vsol
        

                                FEchem_react = FEchem_react + &
                                (- rhopol_charge(i,t)*lambda -psi(i)*rhopolq -betapi*xpol)
                            enddo
                            
                        else 
                            print*,"Error in FEchem_react_multi for systype=nucl_ionbin"
                            print*,"Charged monomer ",t," not acid or base : zpol(1)=",zpol(t,1)," and zpol(2)=",zpol(t,2)
                        endif    

                    endif
                        
                endif        
            enddo


        case("nucl_ionbin_Mg","nucl_ionbin_MgA")
            
            do t=1,nsegtypes

                if(ismonomer_chargeable(t)) then

                    if(t==ta) then 

                        FEchem_react = FEchem_react+FEchempair/volcell ! unit FEchempair allready here in E !!

                    else
                       
                        jcharge=elem_charge(t)

                        if(type_of_charge(t)=="A") then ! acid 
                            
                            vpolstate(1)=vnucl(jcharge,t)
                            vpolstate(2)=vnucl(jcharge,t)
                            vpolstate(3)=vnucl(jcharge,t)+vNa*vsol
                            vpolstate(4)=vnucl(jcharge,t)+vK*vsol
    
                            do i=1,nsize

                                betapi = -log(xsol(i))/vsol
                                lambda = -log(gdisA(i,1,t)) + psi(i) -betapi*vpolstate(1)
                                rhopolq= -gdisA(i,1,t)*rhopol_charge(i,t)
                                
                                xpol=0.0_dp
                                do state=1,4
                                    xpol = xpol + rhopol_charge(i,t)*gdisA(i,state,t)*vpolstate(state)
                                enddo  

                                FEchem_react = FEchem_react + &
                                (- rhopol_charge(i,t)*lambda -psi(i)*rhopolq -betapi*xpol )
                            enddo
                           
                        else if(type_of_charge(t)=="B") then ! base

                            jcharge=elem_charge(t)
                            vpolstate(1)=vnucl(jcharge,t)
                            vpolstate(2)=vnucl(jcharge,t)
                            vpolstate(3)=vnucl(jcharge,t)+vCl*vsol
                    
                            do i=1,nsize

                                betapi = -log(xsol(i))/vsol
                                lambda = -log(gdisB(i,2,t))  - betapi*vpolstate(2)
                                rhopolq= gdisB(i,1,t)*rhopol_charge(i,t)

                                xpol=0.0_dp
                                do state=1,3
                                    xpol = xpol + rhopol_charge(i,t)*gdisB(i,state,t)*vpolstate(state)
                                enddo  
        
                                FEchem_react = FEchem_react + &
                                (- rhopol_charge(i,t)*lambda -psi(i)*rhopolq -betapi*xpol)
                            enddo                            
                        else 
                            print*,"Error in FEchem_react_multi for systype=nucl_ionbin_Mg"
                            print*,"Charged monomer ",t," not acid or base : zpol(1)=",zpol(t,1)," and zpol(2)=",zpol(t,2)
                        endif    

                    endif
                        
                endif        
            enddo


        case("brushborn") 

             ! .. scaled gradient potential contribution 
            call grad_pot_sqr_eps_cubic(psi,epsfcn, Depsfcn,sqrgradpsi)

            do t=1,nsegtypes
        
                if(ismonomer_chargeable(t)) then


                    do i=1,nsize

                        Eself = sqrgradpsi(i)/vsol ! vsol division because of constqE definition check this expresion !!!!!!!!   

                        lbr = lb/epsfcn(i)     ! local Bjerrum length
                        Ebornself = &          ! total self energy  
                            born(lbr,bornrad%pol  ,zpolAA(1))*fdisA(i,1)*rhopol(i,t) + &
                            born(lbr,bornrad%polCa,zpolAA(4))*fdisA(i,4)*rhopol(i,t) + &
                            born(lbr,bornrad%polMg,zpolAA(6))*fdisA(i,6)*rhopol(i,t) + &
                            born(lbr,bornrad%Na,zNa)*xNa(i)/(vNa*vsol)     + & 
                            born(lbr,bornrad%K,zK)*xK(i)/(vK*vsol)         + & 
                            born(lbr,bornrad%Cl,zCl)*xCl(i)/(vCl*vsol)     + &
                            born(lbr,bornrad%Rb,zRb)*xRb(i)/(vRb*vsol)     + & 
                            born(lbr,bornrad%Ca,zCa)*xCa(i)/(vCa*vsol)     + &
                            born(lbr,bornrad%Mg,zMg)*xMg(i)/(vMg*vsol)     + & 
                            born(lbr,bornrad%Hplus,1 )*xHplus(i)/vsol      + &
                            born(lbr,bornrad%OHmin,-1)*xOHmin(i)/vsol  

                        Ebornself=Ebornself*Depsfcn(i)/epsfcn(i)
            
                        betapi=-log(xsol(i))/vsol
                        
                        !  Lagrange multiplier
                        lambda=-log(fdisA(i,1)) -psi(i)*zpolAA(1)-betapi*vpolAA(1)*vsol +& 
                            (Eself+Ebornself)*vpolAA(1)*vsol -born(lbr,bornrad%pol,zpolAA(1))                       
        
                       
                        ! Born energy at i 
                    
                        bornene=(born(lbr,bornrad%pol,zpolAA(1)  )*fdisA(i,1) + &
                                 born(lbr,bornrad%polCa,zpolAA(4))*fdisA(i,4) + &
                                 born(lbr,bornrad%polMg,zpolAA(6))*fdisA(i,6) )*rhopol(i,t)
                        
                        rhopolq = 0.0_dp ! total poly charge at i
                        xpol=0.0_dp      ! total poly volume fraction at i   carefully check with definition of xpol in module field
                        do k=1,4 
                            rhopolq=rhopolq+ zpolAA(k)*fdisA(i,k)*rhopol(i,t)
                            xpol   =xpol + rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol
                        enddo   
                        rhopolq=rhopolq+ zpolAA(6)*fdisA(i,6)*rhopol(i,t)
                        xpol = xpol +rhopol(i,t)*vsol*&
                            (fdisA(i,5)*vpolAA(5)/2.0_dp +fdisA(i,6)*vpolAA(6) +fdisA(i,7)*vpolAA(7)/2.0_dp +fdisA(i,8)*vpolAA(8))

                        FEchem_react = FEchem_react + &
                            ( - rhopol(i,t)*lambda -psi(i)*rhopolq +(-betapi+Eself+Ebornself)*xpol +& 
                            (fdisA(i,5)+fdisA(i,7))*rhopol(i,t)/2.0_dp - bornene )
                       
                    enddo


                endif        
            enddo


        case default
            print*,"systype in FEchem_react_multi wrong"  
        end select      

        FEchem_react=volcell*FEChem_react    
        

    end function FEchem_react_multi



    


    function FEelect_surface() result(FEelsurf)

        use globals, only : bcflag,LEFT,RIGHT, pi
        use volume, only : areacell, nx, ny, delta
        use parameters, only : lb
        use surface, only : sigmaSurfL, sigmaSurfR,sigmaqSurfL, sigmaqSurfR, psiSurfL, psiSurfR

        real(dp) :: FEelsurf(2)

        real(dp) :: sigmaSurf(2),sigmaqSurf(2,ny*nx)
        real(dp) :: sigmaq0Surf(2,nx*ny),psiSurf(2,nx*ny)
        integer :: i, s 

        sigmaSurf(RIGHT)  = sigmaSurfR 
        sigmaSurf(LEFT)   = sigmaSurfL

        do s=1,nx*ny
            sigmaqSurf(RIGHT,s) = sigmaqSurfR(s)
            sigmaqSurf(LEFT,s)  = sigmaqSurfL(s)
            psiSurf(RIGHT,s)    = psiSurfR(s)
            psiSurf(LEFT,s)     = psiSurfL(s)
        enddo    
      
        do i = 1,2
            FEelsurf(i)=0.0_dp
            do s = 1, nx*ny    
                sigmaq0Surf(i,s)=  sigmaqSurf(i,s)/(delta*4.0_dp*pi*lb) ! dimensional charge density  
                FEelsurf(i) = FEelsurf(i)+sigmaq0Surf(i,s) * psiSurf(i,s) /2.0_dp 
            enddo   
            FEelsurf(i)=FEelsurf(i)*areacell
        enddo    

    end function
         
    function FEchem_surface(FEelsurf) result(FEchemsurf)

        use globals, only : bcflag,LEFT,RIGHT, pi
        use volume, only : areacell, nx, ny, delta
        use parameters, only : lb
        use surface, only : sigmaSurfL, sigmaSurfR,sigmaqSurfL, sigmaqSurfR, psiSurfL, psiSurfR
        use surface, only : fdisS, fdisTaR, fdisTaL, qS

        real(dp), intent(in) :: FEelsurf(2)
        real(dp) ::  FEchemsurf(2)

        real(dp) :: sigmaSurf(2),sigmaqSurf(2,ny*nx)
        real(dp) :: sigmaq0Surf(2,nx*ny),psiSurf(2,nx*ny)
        real(dp) :: FEchemSurftmp
        integer :: s 


        sigmaSurf(RIGHT)  = sigmaSurfR 
        sigmaSurf(LEFT)   = sigmaSurfL

        do s=1,nx*ny
            sigmaqSurf(RIGHT,s) = sigmaqSurfR(s)
            sigmaqSurf(LEFT,s)  = sigmaqSurfL(s)
            psiSurf(RIGHT,s)    = psiSurfR(s)
            psiSurf(LEFT,s)     = psiSurfL(s)
        enddo    

        if(bcflag(RIGHT)=='qu') then ! quartz

            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisS(2)) ! area surface integration measure 
            enddo    
            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cl" ) then  ! clay
        
            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+(log(fdisS(2))+qS(2)*psiSurfR(s)) 
            enddo    

            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
        
            FEchemSurf(RIGHT) =(log(fdisS(2))+log(fdisS(5)))*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) - &
                2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
        
            FEchemSurf(RIGHT)=(log(fdisTaR(2))*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb)) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cc") then  
        
            FEchemSurf(RIGHT)=0.0_dp
        
        else
            print*,"Error in FEchem_surface"
            print*,"Wrong value bcflag(RIGHT) : ",bcflag(RIGHT)
            stop
        endif 

        if(bcflag(LEFT)=="ta" ) then ! taurine 

            FEchemSurf(LEFT)= log(fdisTaL(2))*sigmaSurf(LEFT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(LEFT)
        
        elseif(bcflag(LEFT)=="cc") then  
        
            FEchemSurf(LEFT)=0.0_dp
        
        else
            print*,"Error in FEchem_surface"
            print*,"Wrong value bcflag(LEFT) : ",bcflag(LEFT)
        endif 

    end function

 
    function SurfaceCharge(sigmaqSurfR,sigmaqSurfL) result(qsurf)
        
        use globals, only : LEFT,RIGHT
        use volume, only : nx,ny,areacell,delta
        use mathconst
        use parameters, only : lb

        real(dp), intent(in) :: sigmaqSurfR(:),sigmaqSurfL(:)
        real(dp) :: qsurf(2)
        ! local
        integer :: i, s 
        real(dp) :: sigmaSurf(2),sigmaqSurf(2,ny*nx)
  
        do s=1,nx*ny
            sigmaqSurf(RIGHT,s) = sigmaqSurfR(s)
            sigmaqSurf(LEFT,s)  = sigmaqSurfL(s)
        enddo    

        do i=LEFT,RIGHT
            qsurf(i)=0.0_dp
            do s=1,nx*ny     
                qsurf(i) = qsurf(i)+sigmaqSurf(i,s)
            enddo
            qsurf(i)=(qsurf(i)/(4.0_dp*pi*lb*delta))*areacell  ! areacell=area size one surface element, 4*pi*lb*delta make correct dimensional full unit    
        enddo
        
    end function


    function FE_selfenergy_brush()result(FEborn)

        use globals, only : systype, nsize
        use field, only : fdisA,rhopol,xNa,xK,xCl,xCa,xMg,xRb,xHplus,xOHmin,epsfcn
        use volume, only : volcell
        use parameters, only : vsol, vNa,vK, vCl,vCa,vMg,vRb, zNa,zK,zCl,zCa,zMg,zRb,zpolAA 
        use parameters, only : bornrad, lb, tA
        use dielectric_const, only : born

        real(dp) :: FEborn

        ! local variables
        integer :: i
        real(dp) :: Eself,Eborn, lbr

        if(systype/="brushborn") then
            FEborn=0.0_dp
        else 
            Eborn=0.0_dp
            do i=1,nsize
                lbr = lb/epsfcn(i)     ! local Bjerrum length
                Eself = &        ! total self energy  
                    born(lbr,bornrad%pol  ,zpolAA(1))*fdisA(i,1)*rhopol(i,tA) + &
                    born(lbr,bornrad%polCa,zpolAA(4))*fdisA(i,4)*rhopol(i,tA) + &
                    born(lbr,bornrad%polMg,zpolAA(6))*fdisA(i,6)*rhopol(i,tA) + &
                    born(lbr,bornrad%Na,zNa)*xNa(i)/(vNa*vsol)       + & 
                    born(lbr,bornrad%K,zK)*xK(i)/(vK*vsol)           + & 
                    born(lbr,bornrad%Cl,zCl)*xCl(i)/(vCl*vsol)       + &
                    born(lbr,bornrad%Rb,zRb)*xRb(i)/(vRb*vsol)       + & 
                    born(lbr,bornrad%Ca,zCa)*xCa(i)/(vCa*vsol)       + &
                    born(lbr,bornrad%Mg,zMg)*xMg(i)/(vMg*vsol)       + &
                    born(lbr,bornrad%Hplus,1)*xHplus(i)/vsol         + &
                    born(lbr,bornrad%OHmin,-1)*xOHmin(i)/vsol  
                Eborn=Eborn*Eself
            enddo     
            FEborn= Eborn*volcell   
        endif  

    end function FE_selfenergy_brush

    ! Calculates the total number of monomers for all (nsegtypes) monomer type 
    ! post real(d) sumphi(nsegtypes)
 
    function calculate_sumphi()result(sumphi)
    
        use globals, only : nsegtypes, nseg
        use chains, only : type_of_monomer

        real(dp) :: sumphi(nsegtypes)

        ! local arguments
        integer  :: s, t         
        
        sumphi=0

        do s=1,nseg
            t=type_of_monomer(s)
            sumphi(t)=sumphi(t)+1
        enddo

    end function calculate_sumphi

    ! Calculates the total volume  of a conformation that has distributed volume
    ! i.e. systype = nucl_neutral_sv, nucl_ionbin_sv, nucl_ionbin_Mg    
    ! sumvolnucl = \sum_t \sum_j(t) vnucl(j,t)
    
    function calculate_sumvolnucl()result(sumvolnucl)


        use globals, only : nseg
        use chains, only : type_of_monomer, nelem 
        use parameters, only : vnucl

        real(dp) :: sumvolnucl

        ! local arguments
        integer  :: s, t, j        
        
        sumvolnucl=0
        do s=1,nseg
            t=type_of_monomer(s)
            do j=1,nelem(s)
                sumvolnucl=sumvolnucl+vnucl(j,t) 
            enddo
        enddo

    end function calculate_sumvolnucl    

    function check_volume_nucl()result(checksumxpoltot)

        use globals, only : systype

        real(dp) :: checksumxpoltot

        select case (systype)
        case ("nucl_ionbin_sv")

            call check_volume_nucl_ionbin_sv(checksumxpoltot)

        case ("nucl_ionbin_Mg","nucl_ionbin_MgA")

            call check_volume_nucl_ionbin_Mg(checksumxpoltot)

        case ("nucl_neutral_sv ")

            call check_volume_nucl_neutral_sv(checksumxpoltot)

        case ("nucl_ionbin") 

            call check_volume_nucl_ionbin(checksumxpoltot)

        case default

            print*,"Error: systype incorrect in check_volume_nucl"

        end select

    end function check_volume_nucl
        

    ! Check volume for systype="nucl_ionbin_sv"
    ! Integral over xpol compared to expected outcome
    ! output : real(dp)  checksumxpoltot : difference

    subroutine check_volume_nucl_ionbin(checksumxpoltot)

        use globals, only : nsegtypes,nsize
        use field, only   : rhopol, fdisA, gdisA, gdisB, xpol
        use parameters, only : vpolAA,vpol,vsol,zpol,vNa,vK,vCl,ta
        use volume, only : volcell 
        use chains, only : ismonomer_chargeable

        real(dp),intent(inout) :: checksumxpoltot


        integer :: t,i,k,state,ier
        real(dp):: vpolstate(4), sumtmp, sumxpoltot

        if (.not. allocated(sumxpol))  then 
            allocate(sumxpol(nsegtypes),stat=ier)
            if( ier/=0 ) print*,"allocation failure in check_volume_nucl_ionbin"
        endif

        do t=1,nsegtypes

            if(ismonomer_chargeable(t)) then

                if(t==ta) then ! phoshate

                    sumtmp=0.0_dp 
                    do i=1,nsize
                        do k=1,4
                            sumtmp =sumtmp + rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol
                        enddo   
                        sumtmp = sumtmp +rhopol(i,t)*(fdisA(i,5)*vpolAA(5)*vsol/2.0_dp +&
                                                  fdisA(i,6)*vpolAA(6)*vsol +&
                                                  fdisA(i,7)*vpolAA(7)*vsol/2.0_dp+ &
                                                  fdisA(i,8)*vpolAA(8)*vsol)
                    enddo
                    sumxpol(t)=sumtmp *volcell
                else
                     
                    if(zpol(t,1)==0.and.zpol(t,2)==-1) then ! acid 
                        
                        vpolstate(1)=vpol(t)
                        vpolstate(2)=vpol(t)
                        vpolstate(3)=vpol(t)+vNa
                        vpolstate(4)=vpol(t)+vK
                        sumtmp=0.0_dp

                        do i=1,nsize
                            do state=1,4
                                sumtmp = sumtmp + rhopol(i,t)*gdisA(i,state,t)*vpolstate(state)
                            enddo  
                        enddo
                        sumxpol(t)=sumtmp*vsol*volcell
                        
                    else if(zpol(t,1)==1.and.zpol(t,2)==0) then ! base

                        vpolstate(1)=vpol(t)
                        vpolstate(2)=vpol(t)
                        vpolstate(3)=vpol(t)+vCl
                        sumtmp=0.0_dp 
                        
                        do i=1,nsize
                            do state=1,3
                                sumtmp = sumtmp + rhopol(i,t)*gdisB(i,state,t)*vpolstate(state)
                            enddo  
                        enddo
                        sumxpol(t)=sumtmp*vsol*volcell

                    else 
                        print*,"Error in subroutine sumxpol"
                        print*,"Charged monomer ",t," not acid or base : zpol(1)=",zpol(t,1)," and zpol(2)=",zpol(t,2)
                    endif    
                endif
            else ! neutral
                
                sumtmp = sum(rhopol(:,t))
                sumxpol(t)=sumtmp*vpol(t)*vsol*volcell

            endif
        enddo   

        sumxpoltot=sum(xpol)*volcell
        checksumxpoltot=sum(sumxpol)- sumxpoltot

        print*,"sumxpol        = ",(sumxpol(t),t=1,nsegtypes)
        print*,"sumxpoltot     = ",checksumxpoltot         
        print*,"checksumxpoltot= ",checksumxpoltot

    end subroutine check_volume_nucl_ionbin



    ! Check volume for systype="nucl_neutral_sv"
    ! Integral over xpol compared to expected outcome
    ! output : real(dp)  checksumxpoltot : difference
    ! pre sumphi need to be computed 

    subroutine check_volume_nucl_neutral_sv(checksumxpoltot)

        use globals, only : nsegtypes, nseg
        use field, only   : xpol_t  
        use volume, only  : volcell 

        real(dp),intent(inout) :: checksumxpoltot

        integer :: t,i,ier
        real(dp) :: sumvolnucl,sumxpoltot

        if (.not. allocated(sumxpol))  then 
            allocate(sumxpol(nsegtypes),stat=ier)
            if( ier/=0 ) print*,"allocation failure in check_volume_nucl_neutral_sv"
        endif

        do t=1,nsegtypes 
            sumxpol(t) = sum(xpol_t(:,t))*volcell       
        enddo

        sumvolnucl=calculate_sumvolnucl()
        
        checksumxpoltot=sum(sumxpol)- sumvolnucl
        
        if(abs(checksumxpoltot)>epsilon_sumxpol) then 
            print*,"sumxpol        = ",(sumxpol(t),t=1,nsegtypes)
            print*,"sumvolnucl     = ",sumvolnucl
            print*,"checksumxpoltot= ",checksumxpoltot
        endif

    end subroutine check_volume_nucl_neutral_sv

    ! Check volume for systype="nucl_ionbin_sv"
    ! Integral over xpol compared to expected outcome
    ! output : real(dp)  checksumxpoltot : difference

    subroutine check_volume_nucl_ionbin_sv(checksumxpoltot)

        use globals, only : nsegtypes, nsize
        use field, only : xpol_t, rhopol_charge, gdisA, gdisB, fdisA
        use volume, only : volcell 
        use parameters, only : ta
        use parameters, only : vsol,vNa,vK,vCl
        use chains, only     : ismonomer_chargeable, type_of_charge
        
        real(dp),intent(inout) :: checksumxpoltot

        integer :: t,i,ier
        real(dp) :: sumvolnucl,sumxpoltot, sumrhophos
        real(dp), dimension(:), allocatable ::  deltaxpol
        real(dp) :: deltavpolstateCl, deltavpolstateNa, deltavpolstateK

        if (.not. allocated(sumxpol))  then 
            allocate(sumxpol(nsegtypes),stat=ier)
            if( ier/=0 ) print*,"allocation failure in check_volume_nucl_ionbin_sv"
        endif

        if (.not. allocated(deltaxpol))  then 
            allocate(deltaxpol(nsegtypes),stat=ier)
            if( ier/=0 ) print*,"allocation failure in check_volume_nucl_ionbin_sv"
        endif

        do t=1,nsegtypes 
            sumxpol(t) = sum(xpol_t(:,t))     
        enddo
        
       ! add volume due to ion binding 

        deltaxpol =0.0_dp

        do t=1, nsegtypes
            if(ismonomer_chargeable(t)) then 

                if(t/=ta) then
                    if(type_of_charge(t)=="A") then ! acid   

                        deltavpolstateNa=vNa*vsol
                        deltavpolstateK=vK*vsol                                     

                        do i=1,nsize       
                            ! volume fraction only consider Na and K ionpairing
                            deltaxpol(t) = deltaxpol(t)+rhopol_charge(i,t)*(gdisA(i,3,t)*deltavpolstateNa+&
                                gdisA(i,4,t)*deltavpolstateK)
                        enddo

                    else  ! base   

                        deltavpolstateCl=vCl*vsol

                        do i=1,nsize
                            ! volume fraction only consider Cl ionpairing
                            deltaxpol(t) = deltaxpol(t)+rhopol_charge(i,t)*gdisB(i,3,t)*deltavpolstateCl
                        enddo 
                        
                    endif     

                else ! t=tAA phosphate 

                    deltavpolstateNa=vNa*vsol
                    deltavpolstateK=vK*vsol

                    do i=1,nsize
                        ! volume fraction only consider Na and K ionpairing
                        deltaxpol(t) = deltaxpol(t)+rhopol_charge(i,t)*(fdisA(i,3)*deltavpolstateNa+&
                            fdisA(i,8)*deltavpolstateK)
                    enddo
                        
                endif    
            endif
        enddo 

        sumvolnucl=calculate_sumvolnucl()
        checksumxpoltot=sum(sumxpol)*volcell -sum(deltaxpol)*volcell - sumvolnucl
        sumrhophos=sum(rhopol_charge(:,ta))*volcell

        if(abs(checksumxpoltot)>epsilon_sumxpol) then 

            print*,"sumxpol(t)     = ",(sumxpol(t),t=1,nsegtypes)
            print*,"sumxpol        = ",sum(sumxpol)*volcell
            print*,"sumdeltaxpol   = ",sum(deltaxpol)*volcell
            print*,"sumvolnucl     = ",sumvolnucl
            print*,"checksumxpoltot= ",checksumxpoltot
            print*,"sumrhophos     = ",sumrhophos
           
        endif

    end subroutine check_volume_nucl_ionbin_sv


    ! Check volume for systype="nucl_ionbin_Mg"
    ! Integral over xpol compared to expected outcome
    ! output : real(dp)  checksumxpoltot : difference
    ! pre avfdis evaluated first 

    subroutine check_volume_nucl_ionbin_Mg(checksumxpoltot)

        use globals, only : nsegtypes, nsize
        use field, only : xpol_t, rhopol_charge, gdisA, gdisB, fdisA
        use volume, only : volcell 
        use parameters, only : tPhos=>ta
        use parameters, only : vsol,vNa,vK,vCl,vMg,vpol,vPP
        use parameters, only : avfdisA, Phos2Mg
        use chains, only     : ismonomer_chargeable, type_of_charge
        
        real(dp),intent(inout) :: checksumxpoltot

        integer :: t,i,ier
        real(dp) :: sumvolnucl,sumxpoltot, sumrhophos, sumxphos
        real(dp), dimension(:), allocatable ::  deltaxpol
        real(dp) :: deltavpolstateCl, deltavpolstateNa, deltavpolstateK
        real(dp) :: deltavpolstatePNa, deltavpolstatePK, deltavpolstatePMg,deltavpolstateP2Mg
        
        if (.not. allocated(sumxpol))  then 
            allocate(sumxpol(nsegtypes),stat=ier)
            if( ier/=0 ) print*,"allocation failure in check_volume_nucl_ionbin_Mg"
        endif

        if (.not. allocated(deltaxpol))  then 
            allocate(deltaxpol(nsegtypes),stat=ier)
            if( ier/=0 ) print*,"allocation failure in check_volume_nucl_ionbin_Mg"
        endif

        do t=1,nsegtypes 
            sumxpol(t) = sum(xpol_t(:,t))*volcell     
        enddo
        
       ! add volume due to ion binding 

        deltaxpol =0.0_dp

        do t=1, nsegtypes
            if(ismonomer_chargeable(t)) then 

                if(t/=tPhos) then
                    if(type_of_charge(t)=="A") then ! acid   

                        deltavpolstateNa=vNa*vsol
                        deltavpolstateK=vK*vsol                                     

                        do i=1,nsize       
                            ! volume fraction only consider Na and K ionpairing
                            deltaxpol(t) = deltaxpol(t)+rhopol_charge(i,t)*(gdisA(i,3,t)*deltavpolstateNa+&
                                gdisA(i,4,t)*deltavpolstateK)
                        enddo

                         deltaxpol(t)=deltaxpol(t)*volcell

                    else  ! base   

                        deltavpolstateCl=vCl*vsol

                        do i=1,nsize
                            ! volume fraction only consider Cl ionpairing
                            deltaxpol(t) = deltaxpol(t)+rhopol_charge(i,t)*gdisB(i,3,t)*deltavpolstateCl
                        enddo 
                        
                        deltaxpol(t)=deltaxpol(t)*volcell                        
    
                    endif     
                    

                else ! t=tPhos = ta  phosphate 
                
                    ! formula in check_volume_nucl_inbin_sv does not work here because  fdisA(i,k) = 0 !!    

                    deltavpolstatePNa=vNa*vsol
                    deltavpolstatePK=vK*vsol
                    deltavpolstatePMg=vMg*vsol

                    ! deltavpolstateP2Mg=(vpol(ta)+vMg)*vsol  
                    deltavpolstateP2Mg=(vPP(Phos2Mg)-2.0_dp*vpol(tPhos)*vsol)/2.0_dp ! volume change per 1 phosphate
                    
                    sumrhophos=sum(rhopol_charge(:,tPhos))*volcell ! total number of phosphates divide by volcell
                    deltaxpol(tPhos)=sumrhophos*(&
                            avfdisA(3)*deltavpolstatePNa+avfdisA(8)*deltavpolstatePK+&
                            avfdisA(6)*deltavpolstatePMg+avfdisA(7)*deltavpolstateP2Mg) 
           
                endif    
            endif
        enddo 
        
        sumvolnucl=calculate_sumvolnucl() 
        checksumxpoltot=sum(sumxpol)-sum(deltaxpol) - sumvolnucl
        sumxphos=sumrhophos*vpol(tPhos)*vsol


        if(abs(checksumxpoltot)>epsilon_sumxpol) then 
            print*,"Warning: checksumxpoltot larger epsilon_sumxpol ",epsilon_sumxpol
            print*,"checksumxpoltot = ",checksumxpoltot  
            print*,"sumxpol         = ",(sumxpol(t),t=1,nsegtypes)
            print*,"sumdeltaxpol    = ",deltaxpol
            print*,"sumvolnucl      = ",sumvolnucl
            print*,"sumrhophos      = ",sumrhophos
            print*,"sumxphos        = ",sumxphos
            print*,"sumxpol(tPhos)     = ",sumxpol(tPhos)
            print*,"deltaxpol(tPhos)   = ",deltaxpol(tPhos)           

        endif

    end subroutine check_volume_nucl_ionbin_Mg

end module energy
