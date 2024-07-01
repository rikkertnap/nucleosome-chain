! --------------------------------------------------------------|
! fcnMg-excpl.f90:                                              |
! Constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for                |
! nucleosome with inbinding and distubed volume and Mg binding  |
! using phosphate pairs.                                        |
! --------------------------------------------------------------|

module modfcnMgexpl
   
    use mpivars
    use precision_definition

    implicit none

contains
   
     subroutine compute_fdisPP(fdisPP,fdisP2Mg,position1,position2)

        use field, only : xHplus, xOHmin, xNa, xK, xMg, xsol, psi
        use parameters, only : vsol,vNa,vK,vCl,vMg,deltavAA
        use parameters, only : K0aAA,K0a,K0aion
        use parameters, only : Phos,PhosH, PhosK, PhosNa, PhosMg, Phos2Mg 

        real(dp), intent(inout), dimension(:,:) :: fdisPP
        real(dp), intent(inout) :: fdisP2Mg
        integer , intent(in) :: position1, position2

        real(dp) :: xA(3),xB(2),sgxA,sgxB, xP(5,2), xP2Mg, fPP, sumxP         ! disociation variables  
        integer :: i, j     
        integer  :: JJ, KK

        ! .. executable statements 
                                                     
        xP(Phos,1)   = 1.0_dp 
        xP(Phos,2)   = 1.0_dp    
      
        i = position1 ! position in lattice numbers
        j = position2

        xP(PhosH,1)  = xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      !  (PH)/P-    : f(PH)P(i,j)/fPP(i,j)
        xP(PhosH,2)  = xHplus(j)/(K0aAA(1)*(xsol(j)**deltavAA(1)))      !  (PH)/P-    : fP(PH)(i,j)/fPP(i,j)
      
        xP(PhosNa,1) = (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   !  PNa/P-     : f(PNa)P(i,j)/fPP(i,j) 
        xP(PhosNa,2) = (xNa(j)/vNa)/(K0aAA(2)*(xsol(j)**deltavAA(2)))   !  PNa/P-     : fPPNa(i,j)/fPP(i,j) 
      
        xP(PhosK,1)  = (xK(i)/vK)/(K0aAA(7)*(xsol(i)**deltavAA(7)))     !  PK/P-      : f(PK)P(i,j)/fPP(i,j) 
        xP(PhosK,2)  = (xK(j)/vK)/(K0aAA(7)*(xsol(j)**deltavAA(7)))     !  PK/P-      : fPPK(i,j)/fPP(i,j) 
     
        xP(PhosMg,1) = (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! PMgP+/PP2-
        xP(PhosMg,2) = (xMg(j)/vMg)/(K0aAA(5)*(xsol(j)**deltavAA(5)))   ! PPAMg+/PP2-

        xP2Mg  = sqrt( (xMg(i)/vMg)*(xMg(j)/vMg)/ ((K0aAA(6)**2) *(xsol(i)**deltavAA(6))*(xsol(j)**deltavAA(6)))) ! P2Mg/PP2- 
       
        sumxP = 0.0_dp
        do JJ=1,5
            do KK=1,5
                sumxP = sumxP + xP(JJ,1) * xP(KK,2)
            enddo
        enddo
        sumxP=sumxP+xP2Mg

        fPP = 1.0_dp/sumxP    ! fraction of phophate pairs that are both charged
          
        do JJ=1,5             ! fraction of phophate pairs that form a bind with H^+,Na^+,K^+
             do KK=1,5
                 fdisPP(JJ,KK) = fPP * xP(JJ,1) * xP(KK,2)
             enddo
        enddo
            
        fdisP2Mg = fPP * xP2Mg  ! fraction of phophate pairs that form a Mg-bridge
   
                
    end subroutine compute_fdisPP



    ! nucleosome of AA and dna polymers
    ! with ion charegeable group being on one acid (tA) with counterion binding etc 
    ! distribute volume of neighboring cells

    subroutine fcnnucl_Mg_expl(x,f,nn)

        use mpivars
        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, neq, neqint, cuantas, cuantas_no_overlap, DEBUG
        use parameters, only : expmu 
        use parameters, only : vsol,vpol,vNa,vK,vCl,vRb,vCa,vMg,vpolAA,deltavAA,vnucl,vPP
        use parameters, only : zpol,zNa,zK,zCl,zRb,zCa,zMg,qPP,K0aAA,K0a,K0aion
        use parameters, only : ta,isVdW,isrhoselfconsistent,iter
        use parameters, only : Phos,PhosH, PhosK, PhosNa, PhosMg, Phos2Mg 
        use volume, only     : volcell, indexneighbor !, inverse_indexneighbor_phos
        use chains, only     : indexconf, type_of_monomer, logweightchain, nelem, ismonomer_chargeable
        use chains, only     : type_of_charge, elem_charge, indexconfpair, nneigh, maxneigh
        use chains, only     : energychainLJ, no_overlapchain
        use field, only      : xsol,xNa,xCl,xK,xHplus,xOHmin,xRb,xMg,xCa,rhopol,rhopolin,rhoqpol,rhoq
        use field, only      : psi,gdisA,gdisB,fdis,fdisA, rhopol_charge
        use field, only      : fdisPP_loc, fdisPP_loc_swap, fdisP2Mg_loc, fdisP2Mg_loc_swap, rhoqphos
        use field, only      : q, lnproshift, xpol=>xpol_t, xpol_tot=>xpol
        use vectornorm, only : L2norm, L2norm_sub, L2norm_f90
        use Poisson, only    : Poisson_Equation

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
        integer  :: n,i,j,k,l,c,s,ln,t,jcharge,kr,m,mr,ii,ind         ! dummy indices
        integer  :: JJ, KK
        integer  :: ix,iy
        real(dp) :: norm, normvol,normPE, normscf
        real(dp) :: rhopol0 
        real(dp) :: xA(3),xB(2),sgxA,sgxB, xP(5,2),xP2Mg, fPP, sumxP         ! disociation variables 
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)
        integer  :: count_scf
        real(dp) :: deltavpolstateCl, deltavpolstateNa, deltavpolstateK, deltaxpol
        real(dp) :: sum_rhoqphos,sum_xphos

        real(dp) :: K0aPP   ! Kdis of P2Mg pair temporarily define 

        ! .. executable statements 

        ! .. communication between processors

        ! print*,"K0aAA=",K0aAA
        K0aPP=K0aAA(6) ! P2Mg

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER,dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x, neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        endif

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
            xRb(i)      = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)      = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)      = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction

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
        
        do c=1,cuantas                            ! loop over cuantas

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

                            call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                          -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))    
                        enddo
                    endif        
                enddo
            endif         
        enddo

        locallnproshift(1)=lnpro/cuantas_no_overlap
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
         
        do c=1,cuantas                            ! loop over cuantas

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

                            call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc,  k , m)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                          -log(fdisPP_loc(Phos,Phos)))/(2.0_dp*nneigh(s,c))
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
                             
                            call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)
                            call  compute_fdisPP(fdisPP_loc_swap, fdisP2Mg_loc_swap,  m , k)    

                            sum_rhoqphos=0.0_dp
                            sum_xphos=0.0_dp 
                        
                            do JJ=1,5
                                do KK=1,5
                                    sum_rhoqphos = sum_rhoqphos+&
                                        (fdisPP_loc(JJ,KK)*qPP(JJ)+fdisPP_loc_swap(JJ,KK)*qPP(KK))/2.0_dp
                                    sum_xphos = sum_xphos   +&
                                        (fdisPP_loc(JJ,KK)*vPP(JJ)+fdisPP_loc_swap(JJ,KK)*vPP(KK))/2.0_dp
                                enddo
                            enddo
        
                            sum_xphos=sum_xphos+(fdisP2Mg_loc+fdisP2Mg_loc_swap)*vPP(Phos2Mg)/4.0_dp 
                                ! division 4.0_dp  because symmetry and  vPP(Phos2Mg)/2 is volume change per phosphate 
                      
                            local_rhoqphos(k) = local_rhoqphos(k) + pro * sum_rhoqphos /(nneigh(s,c)) ! nneigh could be zero  hence with in loop 
                            local_xpol(k,ta) = local_xpol(k,ta) + pro * sum_xphos /(nneigh(s,c))

                            local_rhopol_charge(k,ta)=local_rhopol_charge(k,ta)+pro/(nneigh(s,c))
                            
                        enddo 
         
                    endif

                enddo
            
            endif 
                
        enddo ! cuantas loop
        
        !   .. import results 

        if (rank==0) then 

            q = 0.0_dp 
            q = local_q
            
            do i=1, numproc-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q = q + local_q
            enddo

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

            do i=1, numproc-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_xpol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        xpol(k,t)=xpol(k,t)+local_xpol(k,t)! polymer density 
                    enddo
                enddo
                
                do t=1,nsegtypes
                    if(ismonomer_chargeable(t)) then
                        call MPI_RECV(local_rhopol_charge(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,&
                                MPI_COMM_WORLD,stat,ierr)
                        do k=1,nsize
                            rhopol_charge(k,t)=rhopol_charge(k,t)+local_rhopol_charge(k,t)! polymer density 
                        enddo
                    endif    
                enddo

                call MPI_RECV(local_rhoqphos, nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                do k=1,nsize
                    rhoqphos(k)=rhoqphos(k)+local_rhoqphos(k) ! density  phosphate charge
                enddo
                
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

                f(i) = xpol_tot(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xK(i)-1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK ! total charge density in units of vsol  

            enddo
          
            ! .. end computation polymer density and charge density  

            ! .. electrostatics 

            call Poisson_Equation(f,psi,rhoq)

            norm=l2norm_f90(f)
            iter=iter+1
                        
            normvol = L2norm_f90(f(1:nsize))
            normPE  = L2norm_f90(f(nsize+1:2*nsize))
           
            print*,'iter=', iter ,'norm=',norm, "normvol=",normvol,"normPE=",normPE
                       
           !  if(DEBUG) call locate_xpol_lager_one(xpol_tot)
           ! print*,"vnucl(1,ta)=",vnucl(1,ta)," ta=",ta
           ! print*,"vPP(Phos)=",vPP(phos)
        
            
        else                      ! Export results 
            
            dest = 0 

            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_xpol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then
                    call MPI_SEND(local_rhopol_charge(:,t), nsize, MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,ierr)
                endif
            enddo

            call MPI_SEND(local_rhoqphos, nsize , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            iter=iter+1

        endif

    end subroutine fcnnucl_Mg_expl



    ! compute the average fraction of charged state of the phosphate pairs 

    subroutine compute_average_charge_PP_expl(avfdisP2Mg,avfdisPP)

        use mpivars
        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, cuantas, DEBUG
        use parameters, only : vsol,vnucl
        use parameters, only : qPP,K0aAA,K0a,K0aion,Phos
        use parameters, only : ta,isVdW! isrhoselfconsistent 
        use volume, only     : nx, ny, nz
        use volume, only     : volcell, inverse_indexneighbor_phos, indexneighbor
        use chains, only     : indexconf, type_of_monomer, logweightchain, nelem, ismonomer_chargeable
        use chains, only     : type_of_charge, elem_charge, indexconfpair, nneigh, maxneigh
        use chains, only     : energychainLJ, no_overlapchain
        use field, only      : xsol,psi,fdis, rhopol_charge, fdisPP_loc, fdisPP_loc_swap, fdisP2Mg_loc, fdisP2Mg_loc_swap
        use field, only      : q, lnproshift
        use myutils, only    : error_handler

        real(dp), intent(inout) :: avfdisP2Mg
        real(dp), intent(inout) :: avfdisPP(5,5)

        !     .. local variables
        
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,kr,m,mr,t,jcharge                ! dummy indices
        integer  :: k_ind, m_ind
        integer  :: JJ, KK
        real(dp) :: local_avfdisP2Mg,local_avfdisPP(5,5)
        real(dp) :: sumrhopairs
        real(dp) :: K0aPP   ! Kdis of P2Mg pair temporarily define 
        integer  :: nsizepsi
        ! .. executable statements 

        ! .. communication between processors 

        K0aPP=K0aAA(6) ! P2Mg
        nsizepsi=nsize+2*Nx*Ny


        local_avfdisPP =0.0_dp
        local_avfdisP2Mg =0.0_dp

        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 

        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsizepsi , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(rhopol_charge(:,ta) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    if(ismonomer_chargeable(t)) then 
                        call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    endif
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(psi , nsizepsi, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(rhopol_charge(:,ta) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then 
                    call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                endif
            enddo

            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

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
              
        do c=1,cuantas         ! loop over cuantas

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

                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc, k ,m)
     
                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                    -log(fdisPP_loc(Phos,Phos)))/(2.0_dp*nneigh(s,c))
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

                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc, k ,m)

                            do JJ=1,5
                                do KK=1,5
                                 local_avfdisPP(JJ,KK) = local_avfdisPP(JJ,KK)+&
                                    fdisPP_loc(JJ,KK)*pro/nneigh(s,c)
                            
                                enddo
                            enddo
                            
                            local_avfdisP2Mg=local_avfdisP2Mg+fdisP2Mg_loc*pro/nneigh(s,c)
                    
                        enddo 
                    endif
                enddo
            endif    
        enddo

       
        !   .. import results 

        if (rank==0) then 

            avfdisP2Mg = local_avfdisP2Mg
            do i=1, numproc-1
                source = i
                call MPI_RECV(local_avfdisP2Mg, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                  avfdisP2Mg = avfdisP2Mg + local_avfdisP2Mg
            enddo

            avfdisPP = local_avfdisPP
            do i=1, numproc-1
                source = i
                call MPI_RECV(local_avfdisPP, 25, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                avfdisPP = avfdisPP + local_avfdisPP
            enddo

            ! .. construction of avfdisP2Mg and avfdisPP 
            ! .. normalized avfdisPP with number of average number pairs = integral of rhopol_charge(:,ta)
 
            sumrhopairs=sum(rhopol_charge(:,tA)) 
            sumrhopairs=sumrhopairs*volcell

            avfdisPP=avfdisPP/(sumrhopairs*q) ! also norm with q
            avfdisP2Mg=avfdisP2Mg/(sumrhopairs*q)
            
        else                      ! Export results 
            
            dest = 0 

            call MPI_SEND(local_avfdisP2Mg, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(local_avfdisPP,25, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

        endif

    end subroutine compute_average_charge_PP_expl


    ! compute the average fraction of charged state of the phosphate pairs 

    subroutine compute_FEchem_react_PP_expl(FEchemPP)

        use mpivars
        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, cuantas, DEBUG
        use parameters, only : vsol,vnucl
        use parameters, only : vPP,qPP,K0aAA,K0a,K0aion,Phos,Phos2Mg, ta 
        use volume, only     : nx, ny, nz
        use volume, only     : volcell, inverse_indexneighbor_phos, indexneighbor
        use chains, only     : indexconf, type_of_monomer, logweightchain, nelem, ismonomer_chargeable
        use chains, only     : type_of_charge, elem_charge, indexconfpair, nneigh, maxneigh
        use chains, only     : energychainLJ, no_overlapchain
        use field, only      : xsol,psi,fdis, rhopol_charge 
        use field, only      : fdisPP_loc, fdisPP_loc_swap, fdisP2Mg_loc, fdisP2Mg_loc_swap
        use field, only      : q, lnproshift
        use myutils, only    : error_handler

        real(dp), intent(inout) :: FEchemPP

        !     .. local variables
        
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,kr,m,mr,t,jcharge                ! dummy indices
        integer  :: JJ, KK
        real(dp) :: local_FEchempair, FEchempair
        real(dp) :: K0aPP   ! Kdis of P2Mg pair temporarily define 
        integer  :: nsizepsi
        real(dp) :: betapi_k, betapi_m, psi_k, psi_m, lambda, sum_pi, sum_psi


        ! .. executable statements 

        ! .. communication between processors 

        K0aPP=K0aAA(6) ! P2Mg
        nsizepsi=nsize+2*Nx*Ny


        local_FEchempair =0.0_dp
       

        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 

        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsizepsi , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(rhopol_charge(:,ta) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    if(ismonomer_chargeable(t)) then 
                        call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    endif
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(psi , nsizepsi, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(rhopol_charge(:,ta) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then 
                    call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                endif
            enddo

            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

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
              
        do c=1,cuantas         ! loop over cuantas

            if(no_overlapchain(c)) then

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

                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc, k, m)

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                    -log(fdisPP_loc(Phos,Phos)))/(2.0_dp*nneigh(s,c))
                        enddo    
                    endif        
                enddo    

                pro = exp(lnpro-lnproshift)   
            
                do s=1,nseg
                    
                    t=type_of_monomer(s)

                    if(t==ta) then 
                                    
                        ! pair density of phosphates 
                        k = indexconf(s,c)%elem(1)
                        !k_ind = inverse_index_phos(k)


                        betapi_k=-log(xsol(k))/vsol
                        psi_k = psi(k)
       
                        do j=1,nneigh(s,c)

                            m = indexconfpair(s,c)%elem(j)

                            betapi_m= -log(xsol(m))/vsol
                            psi_m = psi(m)
                        
                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc, k, m)

                             ! Lagrange multiplier lambd(r,r') 

                            lambda = -(betapi_k +betapi_m)*vPP(Phos) -(psi_k+psi_m)*qPP(Phos) &
                                -log(fdisPP_loc(Phos,Phos))

                            lambda = lambda*pro/nneigh(s,c)        

                            sum_pi  = 0.0_dp
                            sum_psi = 0.0_dp

                            do JJ=1,5
                                do KK=1,5
                                    sum_pi=sum_pi-(vPP(JJ)*betapi_k+vPP(KK)*betapi_m)*fdisPP_loc(JJ,KK)*pro/nneigh(s,c)
                                    sum_psi=sum_psi-(qPP(JJ)*psi_k+qPP(KK)*psi_m)*fdisPP_loc(JJ,KK)*pro/nneigh(s,c)
                                enddo
                            enddo
        
                            sum_pi=sum_pi-(vPP(Phos2Mg)/2.0_dp)*(betapi_k+betapi_m)*fdisP2Mg_loc*pro/nneigh(s,c)

                                ! division 2.0_dp  because  vPP(Phos2Mg)/2 is volume change per phosphate 
                            
                            local_FEchempair = local_FEchempair+(-lambda +sum_pi+sum_psi)/2.0_dp             
                        enddo 
                    endif
                enddo
            endif    
        enddo

       
        !   .. import results 

        if (rank==0) then 

            FEchempair= local_FEchempair
            do i=1, numproc-1
                source = i
                call MPI_RECV(local_FEchempair, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                FEchempair = FEchempair + local_FEchempair
            enddo

            ! .. normalized FEchempair  with q 
            FEchempair = FEchempair/q 

            FEchemPP=FEchempair 

            
        else                      ! Export results 
            
            dest = 0 

            call MPI_SEND(local_FEchempair, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            
        endif

    end subroutine compute_FEchem_react_PP_expl

   

end module modfcnMgexpl

   
