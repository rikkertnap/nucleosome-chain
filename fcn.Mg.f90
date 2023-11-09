! --------------------------------------------------------------|
! fcnMgf90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes/nucleosome                  |
! --------------------------------------------------------------|



module modfcnMg
   
    use mpivars
    use precision_definition

    implicit none

contains
    
    ! nucleosome of AA and dna polymers
    ! with ion charegeable group being on one acid (tA) with counterion binding etc 
    ! distribute volume of neighboring cells

    subroutine fcnnucl_Mg(x,f,nn)

        use mpivars
        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, neq, neqint, cuantas, DEBUG
        use parameters, only : expmu 
        use parameters, only : vsol,vpol,vNa,vK,vCl,vRb,vCa,vMg,vpolAA,deltavAA,vnucl,vPP
        use parameters, only : zpol,zNa,zK,zCl,zRb,zCa,zMg,qPP,K0aAA,K0a,K0aion
        use parameters, only : ta,isVdW,isrhoselfconsistent,iter
        use parameters, only : Phos,PhosH, PhosK, PhosNa, PhosMg, Phos2Mg 
        use volume, only     : volcell, inverse_indexneighbor, indexneighbor
        use chains, only     : indexconf, type_of_monomer, logweightchain, nelem, ismonomer_chargeable
        use chains, only     : type_of_charge, elem_charge, indexconfpair, nneigh, maxneigh
        use field, only      : xsol,xNa,xCl,xK,xHplus,xOHmin,xRb,xMg,xCa,rhopol,rhopolin,rhoqpol,rhoq
        use field, only      : psi,gdisA,gdisB,fdis,fdisA, rhopol_charge, fdisPP, fdisP2Mg, rhopairs
        use field, only      : q, lnproshift, xpol=>xpol_t, xpol_tot=>xpol
        use vectornorm, only : L2norm,L2norm_sub,L2norm_f90
        use VdW, only        : VdW_contribution_lnexp
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
        real(dp) :: local_rhopairs(nsize,maxneigh)
        real(dp) :: local_rhoqphos(nsize)
        real(dp) :: rhoqphos(nsize)

        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize) 
        
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t,jcharge,kr,m,mr                ! dummy indices
        integer  :: JJ, KK
        real(dp) :: norm,normvol,normPE
        real(dp) :: rhopol0 
        real(dp) :: xA(3),xB(2),sgxA,sgxB, xP(5,2),xP2Mg, fPP, sumxP         ! disociation variables 
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)
        integer  :: count_scf
        real(dp) :: deltavpolstateCl, deltavpolstateNa, deltavpolstateK, deltaxpol
        real(dp) :: sum_rhoqphos,sum_xphos

        real(dp) :: K0aPP   ! Kdis of P2Mg pair temporarily define 

        ! real(dp) :: g(neq/2)

        ! .. executable statements 

        ! .. communication between processors 

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
        
        count_scf=0    
        do t=1,nsegtypes
            if(isrhoselfconsistent(t)) then
                count_scf=count_scf+1 
                k=(count_scf+1)*n
                do i=1,n                         
                    rhopolin(i,t) = x(i+k)          ! density 
                enddo
            endif        
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

            lnexppivw(i) = log(xsol(i))/vsol                      ! auxilary variable  divide by vsol  !!
 
        enddo

        print*,'nsegtypes=',nsegtypes
        print*,"maxneigh=",maxneigh

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

                    do i=1,n  
                        print*,"i=",i
                        do kr=1,maxneigh
                            !print*,"i=",i,'kr=',kr
                            j=indexneighbor(i,kr)   ! j=index of neighbors number k of index i 
 
                            xP(Phos,1)   = 1.0_dp 
                            xP(Phos,2)   = 1.0_dp    
                          
                            xP(PhosH,1)  = xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      !  (PH)/P-    : f(PH)P(i,j)/fPP(i,j)
                            xP(PhosH,2)  = xHplus(j)/(K0aAA(1)*(xsol(j)**deltavAA(1)))      !  (PH)/P-    : fP(PH)(i,j)/fPP(i,j)
                          
                            xP(PhosNa,1) = (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   !  PNa/P-     : f(PNa)P(i,j)/fPP(i,j) 
                            xP(PhosNa,2) = (xNa(j)/vNa)/(K0aAA(2)*(xsol(j)**deltavAA(2)))   !  PNa/P-     : fPPNa(i,j)/fPP(i,j) 
                          
                            xP(PhosK,1)  = (xK(i)/vK)/(K0aAA(7)*(xsol(i)**deltavAA(7)))     !  PK/P-      : f(PK)P(i,j)/fPP(i,j) 
                            xP(PhosK,2)  = (xK(j)/vK)/(K0aAA(7)*(xsol(j)**deltavAA(7)))     !  PK/P-      : fPPK(i,j)/fPP(i,j) 
                         
                            xP(PhosMg,1) = (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! PMgP+/PP2-
                            xP(PhosMg,2) = (xMg(j)/vMg)/(K0aAA(5)*(xsol(j)**deltavAA(5)))   ! PPAMg+/PP2-

                            xP2Mg  = sqrt( (xMg(i)/vMg)*(xMg(j)/vMg)/ ( (K0aPP**2)*(xsol(i)**deltavAA(6))*(xsol(j)**deltavAA(6))) ) ! P2Mg/PP2- 
                           
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
                                    fdisPP(i,kr,JJ,KK) = fPP * xP(JJ,1) * xP(KK,2) 
                                enddo
                            enddo
                                
                            fdisP2Mg(i,kr) = fPP * xP2Mg  ! fraction of phophate pairs that form a Mg-bridge
                       
                            lnexppi(i,t) = - psi(i)!!   ! auxilary variable palpha
                        
                        enddo
                    enddo

                endif
            else  

                fdis(:,t)  = 0.0_dp
                lnexppi(:,t) = 0.0_dp

            endif   
        enddo      
               
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes  
                if(isrhoselfconsistent(t)) call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer density fraction      
 
        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro = lnpro+logweightchain(c) 
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
                    print*,"c=",c,"s=",s,"k=",k
                    do jj=1,nneigh(s,c) ! loop neighbors 
                        m = indexconfpair(s,c)%elem(jj)
                        print*,"c=",c,"s=",s,"k=",k," jj=",jj,"m=",m
                        mr = inverse_indexneighbor(k,m) ! relative label of index m relative to k

                        lnpro =lnpro + lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                -log(fdisPP(k,mr,Phos,Phos))
                    enddo    
                endif        
            enddo     

        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
              
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c) 
            do s=1,nseg                       ! loop over segments 
                t = type_of_monomer(s) 
                do j=1,nelem(s)               ! loop over elements of segment 
                    k = indexconf(s,c)%elem(j)
                    lnpro = lnpro +lnexppivw(k)*vnucl(j,t)   ! excluded-volume contribution        
                enddo
                if(ismonomer_chargeable(t)) then
                    jcharge=elem_charge(t)
                    k = indexconf(s,c)%elem(jcharge)
                    lnpro = lnpro + lnexppi(k,t)            !  elect and chemical contribution
                endif
            enddo     

            pro = exp(lnpro-lnproshift)   
            local_q = local_q+pro

            do s=1,nseg
                t=type_of_monomer(s)
                if(t/=ta) then 
                    do j=1,nelem(s)
                        k = indexconf(s,c)%elem(j) 
                        local_xpol(k,t)=local_xpol(k,t)+pro*vnucl(j,t)  ! unnormed polymer volume fraction 
                    enddo
                    if(ismonomer_chargeable(t)) then
                        jcharge=elem_charge(t)
                        k = indexconf(s,c)%elem(jcharge) 
                        local_rhopol_charge(k,t)=local_rhopol_charge(k,t)+pro   ! unnormed density of charge center
                    endif
                else
                
                    ! pair density of phosphates 
                    k = indexconf(s,c)%elem(1)
                    
                    do j=1,nneigh(s,c)

                        m = indexconfpair(s,c)%elem(j)

                        mr = inverse_indexneighbor(k,m) ! mr neighbor label of index m relative to origin at index k
                        kr = inverse_indexneighbor(m,k) ! kr neighbor label of index k relative to origin at index m
                         
                        ! local_rhopairs(k,mr) = local_rhopairs(k,mr)+pro/(2.0_dp*nneigh(s,c))

                        sum_rhoqphos=0.0_dp
                        sum_xphos=0.0_dp
                        do JJ=1,5
                            do KK=1,5
                                sum_rhoqphos = sum_rhoqphos+&
                                    (fdisPP(k,mr,JJ,KK)*qPP(JJ)+fdisPP(m,kr,JJ,KK)*qPP(KK))/2.0_dp
                                sum_xphos = sum_xphos   +&
                                    (fdisPP(k,mr,JJ,KK)*vPP(JJ)+fdisPP(m,kr,JJ,KK)*vPP(KK))/2.0_dp
                            enddo
                        enddo
                        sum_xphos=sum_xphos+fdisP2Mg(k,mr)*vPP(Phos2Mg)/2.0_dp ! check 2.0_dp i.e. single vs double 

                        local_rhoqphos(k) = local_rhoqphos(k) + pro * sum_rhoqphos
                        local_xpol(k,ta) = local_xpol(k,t) + pro * sum_xphos

                    enddo
     
                endif

            enddo
        enddo

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

            !do i=1,nsize
            !    do j=1,nsize
            !        rhopairs(i,j)=local_rhopairs(i,j)
            !    enddo     
            !enddo

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

                ! check how to sent and receive a matrix !
                !call MPI_RECV(local_rhopairs, nsize*nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                !do k=1,nsize
                !    do m=1,nsize
                !        rhopairs(k,m)=rhopairs(k,m)+local_rhopairs(k,m)! density  phosphate  pairs
                !    enddo
                !enddo
                
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
                               
                                rhopol_charge(i,t) = rhopol0 * rhopol_charge(i,t)               ! density nucleosome of type t  
                                rhoqpol(i) = rhoqpol(i) - gdisA(i,1,t)*rhopol_charge(i,t)*vsol   ! total charge density nucleosome in units of vsol 

                                ! volume fraction only consider Na and K ionpairing
                                deltaxpol = rhopol_charge(i,t)*(gdisA(i,3,t)*deltavpolstateNa+gdisA(i,4,t)*deltavpolstateK)
                                xpol(i,t) = rhopol0 * xpol(i,t) + deltaxpol

                            enddo

                        else  ! base   

                            deltavpolstateCl=vCl*vsol

                            do i=1,n
                                
                                rhopol_charge(i,t) = rhopol0 * rhopol_charge(i,t)                ! density nucleosome of type t chargeable 
                                rhoqpol(i) = rhoqpol(i) + gdisB(i,1,t)*rhopol_charge(i,t)*vsol  ! total charge density nucleosome in units of vsol 
                                
                                ! volume fraction only consider Cl ionpairing
                                deltaxpol = rhopol_charge(i,t)*gdisB(i,3,t)*deltavpolstateCl
                                xpol(i,t) = rhopol0 * xpol(i,t) + deltaxpol

                            enddo 
                            
                        endif     

                    else
                        ! t=tAA phosphate 
                        
                        do i=1,n

                            rhopol_charge(i,t)  = rhopol0 * rhopol_charge(i,t)  
                            rhoqpol(i) = rhoqpol(i) + rhopol0 * rhoqphos(i)
                            xpol(i,ta) = rhopol0 * xpol(i,t) 

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

            ! self-consistent equation of densities
            count_scf=0    
            do t=1,nsegtypes
                if(isrhoselfconsistent(t)) then
                    count_scf=count_scf+1 
                    k=(count_scf+1)*n
                    do i=1,n   
                        f(i+k) = rhopol(i,t) - rhopolin(i,t) 
                    enddo
                endif        
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
                        
            normvol = L2norm_f90(f(1:neqint/2))
            normPE  = L2norm_f90(f(neqint/2+1:neqint))
           
            print*,'iter=', iter ,'norm=',norm, "normvol=",normvol,"normPE=",normPE

           !  if(DEBUG) call locate_xpol_lager_one(xpol_tot)
            
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

        endif

    end subroutine fcnnucl_Mg

    


end module modfcnMg

   