! --------------------------------------------------------------|
! fcnCa.f90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes onto a tethered planar surface                   |
! --------------------------------------------------------------|


module fcnpointer

    implicit none

    abstract interface
        subroutine fcn(x,f,n)
            use precision_definition
            implicit none
            integer(8), intent(in) :: n
            real(dp), dimension(n), intent(in) :: x
            real(dp), dimension(n), intent(out) :: f
        end subroutine fcn
    end interface

    procedure(fcn), pointer :: fcnptr => null()

end module fcnpointer

module listfcn
   
    use mpivars
    implicit none

contains

    ! brush of multiblock copolymers

    subroutine fcnelectbrushmulti(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field
        use vectornorm
        use VdW, only : VdW_contribution_lnexp
        use surface
        use Poisson

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t  ! dummy indices
        real(dp) :: norm
        real(dp) :: rhopol0 !integra_q
        integer  :: noffset, count_scf
        real(dp) :: locallnproshift(2), globallnproshift(2)


        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
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
                    rhopolin(i,t) = x(i+k)      ! density 
                enddo
            endif        
        enddo

 
        !  .. assign global and local polymer density 

        do t=1,nsegtypes
            do i=1,n
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction 
            xK(i)      = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK) ! Na+ volume fraction 
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
            xNaCl(i)   = expmu%NaCl*(xsol(i)**vNaCl)
        enddo

        !  fdis(i,t) is assocaited with fraction of monomer of type t at i in state 2 
        !  acid : AH  <=> A- +H+
        !  base : BH+ <=> B + H+

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,n                                         
                    fdis(i,t)  = 1.0_dp/(1.0_dp+xHplus(i)/(K0a(t)*xsol(i)))      
                    lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i) -log(fdis(i,t))   ! auxilary variable palpha
                enddo  
            else
                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    lnexppi(i,t) = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      
       
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes  
                if(isrhoselfconsistent(t)) call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)                    
            enddo
        endif 

        !  .. computation polymer volume fraction      

        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal weight

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
             
             
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)     ! initial weight conformation (1 or 0)
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
            pro=exp(lnpro-lnproshift) 
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo
         
        !   .. import results 

        if (rank==0) then 
         
            q=0.0_dp 
            q=local_q

            do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo 
           
            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
           
            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t) ! polymer density 
                    enddo
                enddo
            enddo     
                
            !  .. construction of fcn and volume fraction polymer   
            !  .. volume polymer segment per volume cell 

            rhopol0=(1.0_dp/volcell)/q  


            do t=1, nsegtypes
                do i=1,n

                    rhopol(i,t) = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                    xpol(i)     = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                    if(ismonomer_chargeable(t)) then 
                        rhoqpol(i)  = rhoqpol(i) + (zpol(t,2)*fdis(i,t)+zpol(t,1)*(1.0_dp-fdis(i,t)))*rhopol(i,t)*vsol ! charge density polymer
                    endif    
                    ! ..scf eq for density
                    !f(i+(t+1)*n)    = rhopol(i,t) - rhopolin(i,t) 
                enddo
            enddo      

            ! .. scf eq for density       
            count_scf=0    
            do t=1,nsegtypes
                if(isrhoselfconsistent(t)) then
                    count_scf=count_scf+1 
                    k=(count_scf+1)*n 
                    do i=1,n   
                        f(i+k)  = rhopol(i,t) - rhopolin(i,t) 
                    enddo
                endif        
            enddo



            do i=1,n
                f(i) = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xNaCl(i) +xK(i) -1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK! total charge density in units of vsol  
            enddo
          
            ! .. end computation polymer density and charge density  

           
            ! .. electrostatics 
    
            
            ! .. Poisson Eq 
            call Poisson_Equation(f,psi,rhoq) 
    

            norm=l2norm(f,neqint)
            iter=iter+1

            print*,'iter=', iter ,'norm=',norm

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

        endif


    end subroutine fcnelectbrushmulti


    ! brush of multiblock copolymers

    subroutine fcnelectbrushmultinoVdW(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field
        use vectornorm
        use VdW, only : VdW_contribution_lnexp
        use surface
        use Poisson

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t  ! dummy indices
        real(dp) :: norm
        real(dp) :: rhopol0 !integra_q
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)

        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
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
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction 
            xK(i)      = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! K+ volume fraction 
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
            xNaCl(i)   = expmu%NaCl*(xsol(i)**vNaCl)
        enddo

        !  fdis(i,t) is assocaited with fraction of monomer of type t at i in state 2 
        !  acid : AH  <=> A- +H+
        !  base : BH+ <=> B + H+

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,n                                         
                    fdis(i,t)  = 1.0_dp/(1.0_dp+xHplus(i)/(K0a(t)*xsol(i)))      
                    !exppi(i,t) = (xsol(i)**vpol(t))*exp(-zpol(t,2)*psi(i) )/fdis(i,t)   

                    lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i) -log(fdis(i,t))   ! auxilary variable palpha

                enddo  
            else
                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    lnexppi(i,t) = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      
              
        !  .. computation polymer volume fraction      

        local_q = 0.0_dp    ! init q             
        lnpro=0.0_dp

        do c=1,cuantas         ! loop over cuantas
            lnpro=lnpro +logweightchain(c)        ! internal weight
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank   
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
               
        lnproshift=globallnproshift(1)
          
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)     ! initial weight conformation (1 or 0)
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo
            pro=exp(lnpro-lnproshift)    
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo
         
        !   .. import results 

        if (rank==0) then 

            q=0.0_dp 
            q=local_q

            do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo 

            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
           
            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t) ! polymer density 
                    enddo
                enddo
            enddo     

            !     .. construction of fcn and volume fraction polymer             
            rhopol0=(1.0_dp/volcell)/q  ! volume polymer segment per volume cell

            do t=1, nsegtypes
                do i=1,n
                    rhopol(i,t) = rhopol0 * rhopol(i,t)       ! density polymer of type t  
                    xpol(i)     = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                    if(ismonomer_chargeable(t)) then 
                        rhoqpol(i)  = rhoqpol(i) + (zpol(t,2)*fdis(i,t)+zpol(t,1)*(1.0_dp-fdis(i,t)))*rhopol(i,t)*vsol ! charge density polymer
                    endif    
                enddo
            enddo        

            do i=1,n
                f(i) = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xNaCl(i) +xK(i)-1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK! total charge density in units of vsol  
            enddo
          
            !  .. end computation polymer density and charge density  

            ! .. electrostatics 
            
            ! .. Poisson Eq 
            call Poisson_Equation(f,psi,rhoq)   

          
            norm=l2norm(f,neqint)
            iter=iter+1

            print*,'iter=', iter ,'norm=',norm

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo
            
        endif


    end subroutine fcnelectbrushmultinoVdW

    ! brush of dna polymers
    ! with ion charegeable group being on one acid (tA) with counterion binding etc 

    subroutine fcnbrushdna(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field
        use vectornorm
        use VdW, only : VdW_contribution_lnexp
        use surface
        use Poisson

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t   ! dummy indices

        real(dp) :: norm
        real(dp) :: rhopol0 
        real(dp) :: xA(7),sumxA, sgxA,qAD, constA, constACa, constAMg ! disociation variables 
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)
        integer  :: count_scf

        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
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
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
            xK(i)      = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! K+ volume fraction
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
        enddo

        !  acid in five chargeable state 
        !  AH   <=> A- + H+  
        !  ANa  <=> A- + Na+ 
        !  ACa+ <=> A- + Ca++  
        !  A2Ca <=> 2A- +Ca++ 
        !  AMg+ <=> A- + Mg++  
        !  A2Mg <=> 2A- +Mg++  
        !  AK   <=> A- + K+ 

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    do i=1,n                                         
                        fdis(i,t)  = 1.0_dp/(1.0_dp+xHplus(i)/(K0a(t)*xsol(i)))      
                        lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i) -log(fdis(i,t))   ! auxilary variable palpha
                    enddo  
                else
                
                    do i=1,n  
                        xA(1)= xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      ! AH/A-
                        xA(2)= (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   ! ANa/A-
                        xA(3)= (xCa(i)/vCa)/(K0aAA(3)*(xsol(i)**deltavAA(3)))   ! ACa+/A-
                        xA(5)= (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! AMg+/A-
                        xA(7)= (xK(i)/vK)/(K0aAA(7)*(xsol(i)**deltavAA(7)))     ! AK/A-
           
                        sgxA=1.0_dp+xA(1)+xA(2)+xA(3)+xA(5) +xA(7)                                                           
                        constACa=(2.0_dp*(rhopolin(i,t)*vsol)*(xCa(i)/vCa))/(K0aAA(4)*(xsol(i)**deltavAA(4))) 
                        constAMg=(2.0_dp*(rhopolin(i,t)*vsol)*(xMg(i)/vMg))/(K0aAA(6)*(xsol(i)**deltavAA(6))) 
                        constA=constACa+constAMg

                        qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

                        fdisA(i,1)  = 1.0_dp/qAD                             ! A-  
                        fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
                        fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
                        fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
                        fdisA(i,5)  = (fdisA(i,1)**2)*constACa               ! A2Ca 
                        fdisA(i,6)  = fdisA(i,1)*xA(5)                       ! AMg+ 
                        fdisA(i,7)  = (fdisA(i,1)**2)*constAMg               ! A2Mg 
                        fdisA(i,8)  = fdisA(i,1)*xA(7)                       ! AK

                        lnexppi(i,t)  = log(xsol(i))*vpol(t)+psi(i)-log(fdisA(i,1))   ! auxilary variable palpha
                        fdis(i,t)   = fdisA(i,1) 
                    enddo  
                endif
            else    
                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    lnexppi(i,t)  = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      
               
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes  
                if(isrhoselfconsistent(t)) call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
 
        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal weight

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
             
         
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c) 
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo 
            pro=exp(lnpro-lnproshift)   
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo

        !   .. import results 

        if (rank==0) then 

            q=0.0_dp 
            q=local_q
            
             do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo

            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
           
            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t)! polymer density 
                    enddo
                enddo
            enddo     


            !  .. construction of fcn and volume fraction polymer 

            !  .. volume polymer segment per volume cell

            rhopol0=(1.0_dp/volcell)/q 

            do t=1, nsegtypes
                if(ismonomer_chargeable(t)) then 

                    if(t/=ta) then
                                                              
                        do i=1,n
                            rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                            rhoqpol(i)   = rhoqpol(i) + (zpol(t,2)*fdis(i,t)+zpol(t,1)*(1.0_dp-fdis(i,t)))*rhopol(i,t)*vsol 
                            xpol(i)      = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                        !    f(i+(t+1)*n) = rhopol(i,t) - rhopolin(i,t)         ! scf eq for density
                        enddo  

                    else

                        do i=1,n
                            rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t 
                            rhoqpol(i)   = rhoqpol(i) + (- fdisA(i,1)+fdisA(i,4)+fdisA(i,6) )*rhopol(i,t)*vsol 
                           ! f(i+(t+1)*n) = rhopol(i,t) - rhopolin(i,t) 
                            do k=1,4               ! polymer volume fraction
                                xpol(i) = xpol(i)+rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol   
                            enddo
                            xpol(i)=xpol(i)+rhopol(i,t)*(fdisA(i,5)*vpolAA(5)/2.0_dp + &
                                                     fdisA(i,6)*vpolAA(6) + &
                                                     fdisA(i,7)*vpolAA(7)/2.0_dp +fdisA(i,8)*vpolAA(8) )*vsol 
                        
                        enddo
                    endif    
                else  
                    do i=1,n
                        rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                        xpol(i)      = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                        !f(i+(t+1)*n) = rhopol(i,t) - rhopolin(i,t)         ! scf eq for density
                    enddo
                endif          
            enddo    

            ! self-consistent equation of densities
            count_scf=0    
            do t=1,nsegtypes
                if(isrhoselfconsistent(t)) then
                    count_scf=count_scf+1 
                    k=(count_scf+1)*n
                    do i=1,n   
                        f(i+k)  = rhopol(i,t) - rhopolin(i,t) 
                    enddo
                endif        
            enddo



            do i=1,n
                f(i) = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xK(i) -1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK ! total charge density in units of vsol  
            enddo
          
            ! .. end computation polymer density and charge density  

            ! .. electrostatics 

            call Poisson_Equation(f,psi,rhoq)

         
            norm=l2norm(f,neqint)
            iter=iter+1
           
            print*,'iter=', iter ,'norm=',norm

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

    
        endif


    end subroutine fcnbrushdna

    ! brush of dna polymers
    ! with ion charegeable group being on one acid (tA) with counterion binding etc 

    subroutine fcnnucl_ionbin(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field
        use vectornorm
        use VdW, only : VdW_contribution_lnexp
        use surface
        use Poisson

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t   ! dummy indices

        real(dp) :: norm, normvol, normPE
        real(dp) :: rhopol0 
        real(dp) :: xA(7),xB(3),sgxA,sgxB
        real(dp) :: qAD, constA, constACa, constAMg ! disociation variables 
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)
        integer  :: count_scf

        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
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
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
            xK(i)      = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! K+ volume fraction
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
        enddo

        !  acid in five chargeable state 
        !  AH   <=> A- + H+  
        !  ANa  <=> A- + Na+ 
        !  ACa+ <=> A- + Ca++  
        !  A2Ca <=> 2A- +Ca++ 
        !  AMg+ <=> A- + Mg++  
        !  A2Mg <=> 2A- +Mg++  
        !  AK   <=> A- + K+ 

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    if(type_of_charge(t)=="A") then  !  acid
                     
                        do i=1,n

                            xA(1) = xHplus(i)/(K0a(t)*xsol(i))           ! AH/A!
                            xA(2) = (xNa(i)/vNa)/(K0aion(t,2))!*xsol(i)) ! ANa/A- :xsol(i)**deltav = xsol(i)**0= 1 
                            xA(3) = (xK(i)/vK)/(K0aion(t,3))!*xsol(i))   ! AK/A-
                            sgxA =  1.0_dp+xA(1)+xA(2)+xA(3)  
                            
                            gdisA(i,1,t) = 1.0_dp/sgxA                    ! A^- 
                            gdisA(i,2,t) = gdisA(i,1,t)*xA(1)             ! AH 
                            gdisA(i,3,t) = gdisA(i,1,t)*xA(2)             ! ANa 
                            gdisA(i,4,t) = gdisA(i,1,t)*xA(3)             ! AK

                            lnexppi(i,t) = log(xsol(i))*vpol(t) +psi(i) -log(gdisA(i,1,t))   ! auxilary variable palpha

                            fdis(i,t) = gdisA(i,1,t)

                        enddo 

                    else !  base

                        do i=1,n
                            xB(1) = (K0a(t)*xsol(i))/xHplus(i)            ! B/BH+
                            xB(2) = (xCl(i)/vCl)/(K0aion(t,2))!*xsol(i))  ! BHCl/BH+
                            sgxB =  1.0_dp+xB(1)+xB(2)  
                            gdisB(i,1,t) = 1.0_dp/sgxB                    ! BH^+
                            gdisB(i,2,t) = gdisB(i,1,t)*xB(1)             ! B
                            gdisB(i,3,t) = gdisB(i,1,t)*xB(2)             ! BHCl 

                            lnexppi(i,t) = log(xsol(i))*vpol(t) -log(gdisB(i,2,t))   !    

                            fdis(i,t) = gdisB(i,2,t)  


                        enddo
            
                    endif 

                else
                
                    do i=1,n  
                        xA(1)= xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      ! AH/A-
                        xA(2)= (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   ! ANa/A-
                        xA(3)= (xCa(i)/vCa)/(K0aAA(3)*(xsol(i)**deltavAA(3)))   ! ACa+/A-
                        xA(5)= (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! AMg+/A-
                        xA(7)= (xK(i)/vK)/(K0aAA(7)*(xsol(i)**deltavAA(7)))     ! AK/A-
           
                        sgxA=1.0_dp+xA(1)+xA(2)+xA(3)+xA(5) +xA(7)                                                           
                        constACa=(2.0_dp*(rhopolin(i,t)*vsol)*(xCa(i)/vCa))/(K0aAA(4)*(xsol(i)**deltavAA(4))) 
                        constAMg=(2.0_dp*(rhopolin(i,t)*vsol)*(xMg(i)/vMg))/(K0aAA(6)*(xsol(i)**deltavAA(6))) 
                        constA=constACa+constAMg

                        qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

                        fdisA(i,1)  = 1.0_dp/qAD                             ! A-  
                        fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
                        fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
                        fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
                        fdisA(i,5)  = (fdisA(i,1)**2)*constACa               ! A2Ca 
                        fdisA(i,6)  = fdisA(i,1)*xA(5)                       ! AMg+ 
                        fdisA(i,7)  = (fdisA(i,1)**2)*constAMg               ! A2Mg 
                        fdisA(i,8)  = fdisA(i,1)*xA(7)                       ! AK

                        lnexppi(i,t)  = log(xsol(i))*vpol(t)+psi(i)-log(fdisA(i,1))   ! auxilary variable palpha
                        fdis(i,t)   = fdisA(i,1) 
                    enddo  
                endif

            else

                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    lnexppi(i,t)  = log(xsol(i))*vpol(t)
                enddo  

            endif   
        enddo      
               
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes  
                if(isrhoselfconsistent(t)) call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
 
        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal weight

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
             
         
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c) 
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo 
            pro=exp(lnpro-lnproshift)   
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo

        !   .. import results 

        if (rank==0) then 

            q=0.0_dp 
            q=local_q
            
             do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo

            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
           
            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t)! polymer density 
                    enddo
                enddo
            enddo     


            !  .. construction of fcn and volume fraction polymer 

            !  .. volume polymer segment per volume cell

            rhopol0=(1.0_dp/volcell)/q 

            do t=1, nsegtypes
                if(ismonomer_chargeable(t)) then 

                    if(t/=ta) then
                        if(zpol(t,1)==0) then ! acid                                        
                            do i=1,n
                                rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                                rhoqpol(i)   = rhoqpol(i) -gdisA(i,1,t)*rhopol(i,t)*vsol 
                                xpol(i)      = xpol(i) + rhopol(i,t)* ( (gdisA(i,1,t) + gdisA(i,2,t) )*vpol(t) + &
                                    gdisA(i,3,t)*(vpol(t)+vNa)+ gdisA(i,4,t)*(vpol(t)+vK))*vsol  ! volume fraction polymer
                            enddo  
                        else  ! base                                        
                            do i=1,n
                                rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                                rhoqpol(i)   = rhoqpol(i) + gdisB(i,1,t)*rhopol(i,t)*vsol 
                                xpol(i)      = xpol(i) + rhopol(i,t)* ( (gdisB(i,1,t) + gdisB(i,2,t) )*vpol(t) + &
                                    gdisB(i,3,t)*(vpol(t)+vCl) )*vsol  ! volume fraction polymer
                            enddo 
                        endif     

                    else

                        do i=1,n
                            rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t 
                            rhoqpol(i)   = rhoqpol(i) + (- fdisA(i,1)+fdisA(i,4)+fdisA(i,6) )*rhopol(i,t)*vsol 
                            do k=1,4               ! polymer volume fraction
                                xpol(i) = xpol(i)+rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol   
                            enddo
                            xpol(i)=xpol(i)+rhopol(i,t)*(fdisA(i,5)*vpolAA(5)/2.0_dp + &
                                                     fdisA(i,6)*vpolAA(6) + &
                                                     fdisA(i,7)*vpolAA(7)/2.0_dp +fdisA(i,8)*vpolAA(8) )*vsol 
                        
                        enddo
                    endif    
                else  
                    do i=1,n
                        rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                        xpol(i)      = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                    enddo
                endif          
            enddo    

            ! self-consistent equation of densities
            count_scf=0    
            do t=1,nsegtypes
                if(isrhoselfconsistent(t)) then
                    count_scf=count_scf+1 
                    k=(count_scf+1)*n
                    do i=1,n   
                        f(i+k)  = rhopol(i,t) - rhopolin(i,t) 
                    enddo
                endif        
            enddo



            do i=1,n
                f(i) = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xK(i)-1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK ! total charge density in units of vsol  
            enddo
          
            ! .. end computation polymer density and charge density  

            ! .. electrostatics 

            call Poisson_Equation(f,psi,rhoq)

         
            norm=l2norm(f,neqint)
            iter=iter+1


            normvol=L2norm_sub(f,1,neqint/2)
            normPE=L2norm_sub(f,neqint/2+1,neqint)
           
            print*,'iter=', iter ,'norm=',norm, "normvol=",normvol, "normPE=",normPE
            print*,'iter=', iter ,'norm=',norm

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

    
        endif


    end subroutine fcnnucl_ionbin


    ! nucleosome of AA and dna polymers
    ! with ion charegeable group being on one acid (tA) with counterion binding etc 
    ! distribute voluem of neighborign cell

    subroutine fcnnucl_ionbin_sv(x,f,nn)

        use mpivars
        use precision_definition
        use globals, only : nsize, nsegtypes, nseg, neq, neqint, cuantas, DEBUG
        !use parameters, Tlocal=>Tref
        use parameters, only : expmu 
        use parameters, only : vsol,vpol,vNa,vK,vCl,vRb,vCa,vMg,vpolAA,deltavAA,deltavnucl
        use parameters, only : zpol,zNa,zK,zCl,zRb,zCa,zMg,K0aAA,K0a,K0aion
        use parameters, only : ta,isVdW,isrhoselfconsistent,iter
        use volume, only : volcell
        use chains, only : indexchain, type_of_monomer, logweightchain, ismonomer_chargeable ,type_of_charge    
        use field, only : xsol,xpol,xNa,xCl,xK,xHplus,xOHmin,xRb,xMg,xCa,rhopol,rhopolin,rhoqpol,rhoq
        use field, only : psi,gdisA,gdisB,fdis,fdisA
        use field, only : q, lnproshift
        use vectornorm, only : L2norm,L2norm_sub,L2norm_f90
        use VdW, only : VdW_contribution_lnexp
        use Poisson, only : Poisson_Equation
        use fcnaux, only : compute_lnexppi_neutral, compute_lnexppi_acid, compute_lnexppi_base
        use fcnaux, only : compute_xpol_neutral, compute_xpol_chargeable

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t   ! dummy indices
        real(dp) :: norm,normvol, normPE
        real(dp) :: rhopol0 
        real(dp) :: xA(7),xB(3),sgxA,sgxB
        real(dp) :: qAD, constA, constACa, constAMg ! disociation variables 
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)
        integer  :: count_scf

        real(dp) :: g(neq/2)

        !     .. executable statements 

        !     .. communication between processors 

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
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
            xK(i)      = expmu%K* (xsol(i)**vK) *exp(-psi(i)*zK)  ! K+ volume fraction
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))  *exp(-psi(i))     ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))  *exp(+psi(i))     ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
        enddo

        !  acid in five chargeable state 
        !  AH   <=> A- + H+  
        !  ANa  <=> A- + Na+ 
        !  ACa+ <=> A- + Ca++  
        !  A2Ca <=> 2A- +Ca++ 
        !  AMg+ <=> A- + Mg++  
        !  A2Mg <=> 2A- +Mg++  
        !  AK   <=> A- + K+ 

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

                        enddo 

                        fdis(:,t) = gdisA(:,1,t)
                        call compute_lnexppi_acid(xsol,psi,gdisA(:,1,t),deltavnucl(:,:,:,1,t),lnexppi(:,t)) 

                    else !  base

                        do i=1,n
                            xB(1) = (K0a(t)*xsol(i))/xHplus(i)            ! B/BH+
                            xB(2) = (xCl(i)/vCl)/(K0aion(t,2))!*xsol(i))  ! BHCl/BH+
                            sgxB =  1.0_dp+xB(1)+xB(2)  
                            gdisB(i,1,t) = 1.0_dp/sgxB                    ! BH^+
                            gdisB(i,2,t) = gdisB(i,1,t)*xB(1)             ! B
                            gdisB(i,3,t) = gdisB(i,1,t)*xB(2)             ! BHCl     

                        enddo

                        fdis(:,t) = gdisB(:,2,t)  
                        call compute_lnexppi_base(xsol,gdisB(:,2,t),deltavnucl(:,:,:,2,t),lnexppi(:,t)) 

                    endif  
                                
                else
                    ! t=ta : phosphate

                    do i=1,n  

                        xA(1)= xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      ! AH/A-
                        xA(2)= (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   ! ANa/A-
                        xA(3)= (xCa(i)/vCa)/(K0aAA(3)*(xsol(i)**deltavAA(3)))   ! ACa+/A-
                        xA(5)= (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! AMg+/A-
                        xA(7)= (xK(i)/vK)/(K0aAA(7)*(xsol(i)**deltavAA(7)))     ! AK/A-
           
                        sgxA=1.0_dp+xA(1)+xA(2)+xA(3)+xA(5) +xA(7)                                                           
                        constACa=(2.0_dp*(rhopolin(i,t)*vsol)*(xCa(i)/vCa))/(K0aAA(4)*(xsol(i)**deltavAA(4))) 
                        constAMg=(2.0_dp*(rhopolin(i,t)*vsol)*(xMg(i)/vMg))/(K0aAA(6)*(xsol(i)**deltavAA(6))) 
                        constA=constACa+constAMg

                        qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

                        fdisA(i,1)  = 1.0_dp/qAD                             ! A-  
                        fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
                        fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
                        fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
                        fdisA(i,5)  = (fdisA(i,1)**2)*constACa               ! A2Ca 
                        fdisA(i,6)  = fdisA(i,1)*xA(5)                       ! AMg+ 
                        fdisA(i,7)  = (fdisA(i,1)**2)*constAMg               ! A2Mg 
                        fdisA(i,8)  = fdisA(i,1)*xA(7)                       ! AK
                        
                    enddo

                    fdis(:,ta)    = fdisA(:,1)

                    gdisA(:,1,ta) = fdisA(:,1)  ! A-
                    gdisA(:,2,ta) = fdisA(:,2)  ! AH
                    gdisA(:,3,ta) = fdisA(:,3)  ! ANa
                    gdisA(:,4,ta) = fdisA(:,8)  ! AK 


                    call compute_lnexppi_acid(xsol,psi,fdisA(:,1),deltavnucl(:,:,:,1,ta),lnexppi(:,ta)) 

                endif
            else  

                fdis(:,t)  = 0.0_dp
                call compute_lnexppi_neutral(xsol,deltavnucl(:,:,:,1,t),lnexppi(:,t))
                 
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

            lnpro=lnpro+logweightchain(c)        ! internal weight

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
              
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c) 
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo 
            pro=exp(lnpro-lnproshift)   
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo

        !   .. import results 

        if (rank==0) then 

            q=0.0_dp 
            q=local_q
            
             do i=1, numproc-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo

            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
           
            do i=1, numproc-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t)! polymer density 
                    enddo
                enddo
            enddo     


            !  .. construction of fcn and volume fraction polymer 

            !  .. volume polymer segment per volume cell

            rhopol0=(1.0_dp/volcell)/q 

            do t=1, nsegtypes
                if(ismonomer_chargeable(t)) then 

                    if(t/=ta) then
                        if(type_of_charge(t)=="A") then ! acid                                        

                            do i=1,n
                                rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                                rhoqpol(i)   = rhoqpol(i) -gdisA(i,1,t)*rhopol(i,t)*vsol 
                            enddo

                            call compute_xpol_chargeable(rhopol(:,t),gdisA(:,:,t),deltavnucl(:,:,:,:,t),xpol)

                        else  ! base                                        
                            
                            do i=1,n
                                rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                                rhoqpol(i)   = rhoqpol(i) + gdisB(i,1,t)*rhopol(i,t)*vsol 
                            enddo 
                            
                            call compute_xpol_chargeable(rhopol(:,t),gdisB(:,:,t),deltavnucl(:,:,:,:,t),xpol)

                        endif     

                    else

                        do i=1,n
                            rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                            rhoqpol(i)   = rhoqpol(i) -gdisA(i,1,t)*rhopol(i,t)*vsol 
                        enddo

                        call compute_xpol_chargeable(rhopol(:,t),gdisA(:,:,t),deltavnucl(:,:,:,:,t),xpol)
                            
                            
                    endif    
                else  

                    ! density polymer of type t 
                    do i=1,n
                        rhopol(i,t)  = rhopol0 * rhopol(i,t)    
                    enddo

                    call compute_xpol_neutral(rhopol(:,t),deltavnucl(:,:,:,1,t),xpol)

                endif          
            enddo    

            ! self-consistent equation of densities
            count_scf=0    
            do t=1,nsegtypes
                if(isrhoselfconsistent(t)) then
                    count_scf=count_scf+1 
                    k=(count_scf+1)*n
                    do i=1,n   
                        f(i+k)  = rhopol(i,t) - rhopolin(i,t) 
                    enddo
                endif        
            enddo

            do i=1,n

                f(i) = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xK(i)-1.0_dp
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
            
        else                      ! Export results 
            
            dest = 0 

            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

        endif


    end subroutine fcnnucl_ionbin_sv

    ! Nucleosome of AA and dna polymers distributed volume
    ! neutral system

    subroutine fcnnucl_neutral_sv(x,f,nn)

        use mpivars
        use precision_definition
        use globals, only : nsize, nsegtypes, nseg, neq, neqint, cuantas, DEBUG
        use parameters, only : vsol,vpol, vnucl
        use parameters, only : iter
        use volume, only : volcell
        use chains, only : indexconf, type_of_monomer, logweightchain, nelem  
        use field, only : xsol, xpol=>xpol_t
        use field, only : q, lnproshift
        use vectornorm, only : L2norm_f90
        

        ! .. scalar arguments

        integer(8), intent(in) :: nn

        ! .. array arguments

        real(dp), intent(in)  :: x(neq)
        real(dp), intent(out) :: f(neq)

        ! .. local variables
        
        real(dp) :: local_xpol(nsize,nsegtypes),xpol_tot(nsize)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize)                    ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,j,k,c,s,t                       ! dummy indices
        real(dp) :: norm, normvol
        real(dp) :: xpol0 
        real(dp) :: locallnproshift(2), globallnproshift(2)
      
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
 
        do i=1,nsize                     
            xsol(i) = x(i)        ! volume fraction solvent
        enddo  
    
        !  .. assign global and local polymer volume fraction 
        do t=1,nsegtypes
            do i=1,nsize
                xpol(i,t)=0.0_dp 
                local_xpol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,nsize                              
            lnexppi(i) = log(xsol(i)/vsol)      ! auxilary variable DO not divide by vsol not nessary !!!!!!!!!
        enddo      
               
    
        !  .. computation polymer density fraction      
 
        local_q = 0.0_dp    ! init q
        lnpro   = 0.0_dp
        
        do c=1,cuantas                           ! loop over cuantas
            lnpro=lnpro+logweightchain(c)        ! internal weight
            do s=1,nseg                          ! loop over segments 
                t=type_of_monomer(s)                
                do j=1,nelem(s)                  ! loop over element of segment
                    k=indexconf(s,c)%elem(j)
                    lnpro = lnpro +lnexppi(k)*vnucl(j,t)
                enddo    
            enddo 
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
              
        do c=1,cuantas                        ! loop over cuantas
            lnpro = logweightchain(c) 
            do s=1,nseg                       ! loop over segments 
                do j=1,nelem(s)               ! loop over elements of segment  
                    k = indexconf(s,c)%elem(j)
                    lnpro = lnpro +lnexppi(k)*vnucl(j,t)              
                enddo
            enddo     

            pro = exp(lnpro-lnproshift)   
            local_q = local_q+pro
            
            do s=1,nseg
                t=type_of_monomer(s)
                do j=1,nelem(s)
                    k = indexconf(s,c)%elem(j) 
                    local_xpol(k,t)=local_xpol(k,t)+pro*vnucl(j,t) 
                enddo
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

            ! conformation on rank zero 
            do t=1,nsegtypes
                do i=1,nsize
                    xpol(i,t)=local_xpol(i,t) ! polymer volume fraction 
                enddo
            enddo
           
            do i=1, numproc-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_xpol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        xpol(k,t)=xpol(k,t)+local_xpol(k,t)! polymer volume fraction
                    enddo
                enddo
            enddo     

            !  .. construction of fcn and volume fraction polymer 
            !  .. normalization volume polymer segment per volume cell

            xpol0=(1.0_dp/volcell)/q 

            do t=1, nsegtypes           ! volumer fraction of polymer of type t 
                do i=1,nsize
                    xpol(i,t)  = xpol0 * xpol(i,t)    
                enddo   
            enddo    

            do i=1,nsize
                xpol_tot(i)=0.0_dp
                do t=1,nsegtypes
                    xpol_tot(i)  = xpol_tot(i) + xpol(i,t)    
                enddo  

                f(i) = xpol_tot(i)+xsol(i)-1.0_dp
            enddo
          
            ! .. end computation polymer volume fraction 

            norm = l2norm_f90(f)
            iter = iter + 1 
                         
            print*,'iter=', iter ,'norm=',norm
            
        else                      ! Export results 
            
            dest = 0 

            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_xpol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

        endif

    end subroutine fcnnucl_neutral_sv




    ! brush of multiblock copolymers
    ! with ion charegeable group being an acid with counterion binding etc 

    subroutine fcnbrush(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field
        use vectornorm
        use VdW, only : VdW_contribution_lnexp
        use surface
        use Poisson

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t  ! dummy indices
        real(dp) :: norm
        real(dp) :: rhopol0 
        real(dp) :: xA(7),sumxA, sgxA,qAD, constA, constACa, constAMg ! disociation variables 
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)


        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
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
        do t=1,nsegtypes
            k=(t+1)*n
            do i=1,n 
                rhopolin(i,t) = x(i+k) ! density 
            enddo    
        enddo
             
        !  .. assign global and local polymer density 
        do t=1,nsegtypes
            do i=1,n
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction 
            xK(i)      = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! K+ volume fraction
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
        enddo

        !  AH   <=> A- + H+  
        !  ANa  <=> A- + Na+ 
        !  ACa+ <=> A- + Ca++  
        !  A2Ca <=> 2A- +Ca++ 
        !  AMg+ <=> A- + Mg++  
        !  A2Mg <=> 2A- +Mg++  

        do t=1,nsegtypes

            if(ismonomer_chargeable(t)) then

                do i=1,n  
                    xA(1)= xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      ! AH/A-
                    xA(2)= (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   ! ANa/A-
                    xA(3)= (xCa(i)/vCa)/(K0aAA(3)*(xsol(i)**deltavAA(3)))   ! ACa+/A-
                    xA(5)= (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! AMg+/A-
                    xA(7)= (xK(i)/vK)/(K0aAA(7)*(xsol(i)**deltavAA(7)))     ! AK/A-
       
                    sgxA=1.0_dp+xA(1)+xA(2)+xA(3)+xA(5)+xA(7)                                                        
                    constACa=(2.0_dp*(rhopolin(i,t)*vsol)*(xCa(i)/vCa))/(K0aAA(4)*(xsol(i)**deltavAA(4))) 
                    constAMg=(2.0_dp*(rhopolin(i,t)*vsol)*(xMg(i)/vMg))/(K0aAA(6)*(xsol(i)**deltavAA(6))) 
                    constA=constACa+constAMg

                    qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

                    fdisA(i,1)  = 1.0_dp/qAD                             ! A-  
                    fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
                    fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
                    fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
                    fdisA(i,5)  = (fdisA(i,1)**2)*constACa               ! A2Ca 
                    fdisA(i,6)  = fdisA(i,1)*xA(5)                       ! AMg+ 
                    fdisA(i,7)  = (fdisA(i,1)**2)*constAMg               ! A2Mg 
                    fdisA(i,8)  = fdisA(i,1)*xA(7)                       ! AK 
                    
                    lnexppi(i,t)  = log(xsol(i))*vpol(t)+psi(i)-log(fdisA(i,1))   ! auxilary variable palpha
                    fdis(i,t)   = fdisA(i,1) 
                enddo  
            else
                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    lnexppi(i,t)  = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      

               
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes  
                call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      

        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal weight

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       

        lnproshift=globallnproshift(1)
              
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)  ! initial weight conformation
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro+lnexppi(k,t)
            enddo 
            pro=exp(lnpro-lnproshift)   
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo

        !   .. import results 

        if (rank==0) then 

            q=0.0_dp 
            q=local_q

            do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo 

            ! first graft point 

            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
           
            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t) ! polymer density 
                    enddo
                enddo
            enddo     

            !     .. construction of fcn and volume fraction polymer             
            rhopol0=(1.0_dp/volcell)/q ! volume polymer segment per volume cell

            do t=1, nsegtypes
                if(ismonomer_chargeable(t)) then 

                    do i=1,n
                        rhopol(i,t) = rhopol0 * rhopol(i,t)               ! density polymer of type t 
                        rhoqpol(i)  = rhoqpol(i) + (- fdisA(i,1)+fdisA(i,4)+fdisA(i,6) )*rhopol(i,t)*vsol 
                        f(i+(t+1)*n) = rhopol(i,t) - rhopolin(i,t) 
                        do k=1,4               ! polymer volume fraction
                            xpol(i)=xpol(i)+rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol   
                        enddo
                        xpol(i)=xpol(i)+rhopol(i,t)*(fdisA(i,5)*vpolAA(5)/2.0_dp + &
                                                 fdisA(i,6)*vpolAA(6) + &
                                                 fdisA(i,7)*vpolAA(7)/2.0_dp + fdisA(i,8)*vpolAA(8) )*vsol
                    
                    enddo
                else  
                    do i=1,n
                        rhopol(i,t) = rhopol0 * rhopol(i,t)               ! density polymer of type t 
                        xpol(i) = xpol(i) + rhopol(i,t)*vpol(t)*vsol 
                        f(i+(t+1)*n) = rhopol(i,t) - rhopolin(i,t)         ! scf eq for density
                    enddo
                endif          
            enddo 

            do i=1,n
                f(i) = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xK(i) -1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK! total charge density in units of vsol  
            enddo
          
            ! .. end computation polymer density and charge density  

            ! .. electrostatics 
             
            call Poisson_Equation(f,psi,rhoq)

            !  .. end electrostatics   
         
            norm=l2norm(f,neqint)
            iter=iter+1

            print*,'iter=', iter ,'norm=',norm

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo
         
        endif


    end subroutine fcnbrush

    subroutine fcnbrushborn(x,f,nn)

        use mpivars
        use globals
        use parameters
        use volume
        use chains, only : indexchain, type_of_monomer,logweightchain, ismonomer_chargeable
        use field
        use vectornorm
        use VdW, only : VdW_contribution_lnexp
        use dielectric_const, only : dielectfcn, born
        use surface
        use Poisson, only : Poisson_Equation_Eps, grad_pot_sqr_eps_cubic

        !  .. scalar arguments

        integer(8), intent(in) :: nn

        ! .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

    !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)         ! auxilairy variable for computing P(\alpha) 
        real(dp) :: phi(nsize)                     ! phi=1-xsol-sum_ix_i 
        real(dp) :: pro,lnpro
        real(dp) :: lbr,expborn,avgvol,Etotself,expsqrgrad, Eself
        real(dp) :: expsqrgradpsi(nsize),expdeltaGAA(nsize,6),expEtotself(nsize)
        integer  :: n,i,j,k,l,c,s,ln,t,k1,k2,k3,k4,k5          !indices
        real(dp) :: norm
        real(dp) :: rhopol0 
        real(dp) :: xA(7),sumxA, sgxA,qAD, constA, constACa, constAMg ! disociation variables 
        integer  :: noffset,neq_bc
        integer  :: count_sc
        real(dp) :: rhopolAA(nsize),rhopolACa(nsize),rhopolAMg(nsize)
                    ! AA because  conflict with definition rhopolA in field.f90
        real(dp) :: locallnproshift(2), globallnproshift(2)
        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER,dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x, neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        endif

        ! read out x 
        n=nsize             
        k1=n
        k2=2*n
        k3=3*n
        k4=4*n
        k5=5*n
        
        do i=1,n                     
            xsol(i)      = x(i)      ! volume fraction solvent
            psi(i)       = x(i+k1)   ! potential
            phi(i)       = x(i+k2)   ! phi = 1-\sum_k={ions,sol}x_k
            rhopolAA(i)  = x(i+k3)   ! rhpolAA(i)  = fdisA(i,1)*rhopolin(i,tA)
            rhopolACa(i) = x(i+k4)   ! rhopolACa(i)= fdisA(i,4)*rhopolin(i,tA)
            rhopolAMg(i) = x(i+k5)   ! rhopolAMg(i)= fdisA(i,6)*rhopolin(i,tA)
        enddo           
        
        count_sc=0   
        do t=1,nsegtypes
            if(isrhoselfconsistent(t)) then
                count_sc=count_sc+1 
                k=count_sc*n+k5 
                do i=1,n                         
                    rhopolin(i,t) = x(i+k)      ! density 
                enddo
            endif        
        enddo

        neq_bc=0
        k=count_sc+1
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nx*ny
            do i=1,neq_bc
                psiSurfR(i) =x(k+i)                 ! surface potentail
            enddo
        endif   
        if(bcflag(LEFT)/="cc") then 
            do i=1,nx*ny
                psiSurfL(i) =x(k+neq_bc+i)          ! surface potentail
            enddo
            neq_bc=neq_bc+nx*ny
        endif
        
        !  .. assign global and local polymer density 

        do t=1,nsegtypes
            do i=1,n
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    

        ! value dielectric permittivity
        call dielectfcn(phi,epsfcn,Depsfcn,dielectP,dielectW,nsize) 

        ! surface charge 
        !sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
        !sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)
        
        ! gradient potential contribution to PDF
        call grad_pot_sqr_eps_cubic(psi,epsfcn, Depsfcn,expsqrgradpsi)

       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer

            lbr = lb/epsfcn(i)     ! local Bjerrum length

            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-born(lbr,bornrad%Na,zNa)-psi(i)*zNa) ! Na+ volume fraction 
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-born(lbr,bornrad%Cl,zCl)-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))  *exp(-born(lbr,bornrad%Hplus,1)-psi(i))    ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))  *exp(-born(lbr,bornrad%OHmin,-1)+psi(i))   ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-born(lbr,bornrad%Rb,zRb)-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-born(lbr,bornrad%Ca,zCa)-psi(i)*zCa) ! Ca++ volume fraction 
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-born(lbr,bornrad%Mg,zMg)-psi(i)*zMg) ! Mg++ volume fraction 



            ! xNaCl(i)   = expmu%NaCl*(xsol(i)**vNaCl) ! for this fcn KionNa=0 => xNaCl=0 

            lbr = lb/epsfcn(i)     ! local Bjerrum length
            Etotself = &           ! total self energy    
                born(lbr,bornrad%pol  ,zpolAA(1))*rhopolAA(i)  + & ! rhpolAA(i)  = fdisA(i,1)*rhopolin(i,tA)
                born(lbr,bornrad%polCa,zpolAA(4))*rhopolACa(i) + & ! rhopolACa(i)= fdisA(i,4)*rhopolin(i,tA)
                born(lbr,bornrad%polMg,zpolAA(6))*rhopolAMg(i) + & ! rhopolAMg(i)= fdisA(i,6)*rhopolin(i,tA)
                born(lbr,bornrad%Na,zNa)*xNa(i)/(vNa*vsol)     + & 
                born(lbr,bornrad%Cl,zCl)*xCl(i)/(vCl*vsol)     + &
                born(lbr,bornrad%Rb,zRb)*xRb(i)/(vRb*vsol)     + & 
                born(lbr,bornrad%Ca,zCa)*xCa(i)/(vCa*vsol)     + &
                born(lbr,bornrad%Mg,zMg)*xMg(i)/(vMg*vsol)     + &
                born(lbr,bornrad%Hplus,1 )*xHplus(i)/vsol      + &
                born(lbr,bornrad%OHmin,-1)*xOHmin(i)/vsol

            expEtotself(i) = Etotself*(Depsfcn(i)/epsfcn(i))  
           
            Eself=expsqrgradpsi(i)+expEtotself(i)

            expdeltaGAA(i,1)=exp(-(born(lbr,bornrad%pol,-1)-bornbulk%pol)-(born(lbr,bornrad%Hplus,1)-bornbulk%Hplus) +&      !  AH   <=> A- + H+
                Eself*(vpolAA(1)-vpolAA(2)))
            expdeltaGAA(i,2)=exp(-(born(lbr,bornrad%pol,-1)-bornbulk%pol)-(born(lbr,bornrad%Na,1)-bornbulk%Na) + & 
                Eself*(vpolAA(1)-vpolAA(3)) )                                                                                !  ANa  <=> A- + Na+ 
            expdeltaGAA(i,3)=exp(-(born(lbr,bornrad%pol,-1)-bornbulk%pol)-(born(lbr,bornrad%Ca,2)-bornbulk%Ca) + &        
                (born(lbr,bornrad%polCa,+1)-bornbulk%polCa)+&
                Eself*(vpolAA(1)-vpolAA(4)))                                                                                 !  ACa+  <=> A- + Ca++
            expdeltaGAA(i,4)=exp(-2.0_dp*(born(lbr,bornrad%pol,-1)-bornbulk%pol)-(born(lbr,bornrad%Ca,2)-bornbulk%Ca) +&     !  A2Ca  <=> 2A- +Ca++ 
                Eself*(2.0_dp*vpolAA(1)-vpolAA(5))) 
            expdeltaGAA(i,5)=exp(-(born(lbr,bornrad%pol,-1)-bornbulk%pol)-(born(lbr,bornrad%Mg,2)-bornbulk%Mg) + &        
                (born(lbr,bornrad%polMg,+1)-bornbulk%polMg)+&
                Eself*(vpolAA(1)-vpolAA(6)))                                                                                 !  AMg+  <=> A- + Mg++
            expdeltaGAA(i,6)=exp(-2.0_dp*(born(lbr,bornrad%pol,-1)-bornbulk%pol)-(born(lbr,bornrad%Mg,2)-bornbulk%Mg) +&     !  A2Mg  <=> 2A- +Mg++ 
                Eself*(2.0_dp*vpolAA(1)-vpolAA(7))) 

        enddo
        
        !  AH   <=> A- + H+  
        !  ANa  <=> A- + Na+ 
        !  ACa+ <=> A- + Ca++  
        !  A2Ca <=> 2A- +Ca++ 
        !  AMg+ <=> A- + Mg++  
        !  A2Mg <=> 2A- +Mg++  

        do i=1,n  
            xA(1)= xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0aAA(3)*(xsol(i)**deltavAA(3)))   ! ACa+/A-
            xA(5)= (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! AMg+/A-

            sgxA=1.0_dp+xA(1)+xA(2)+xA(3)+xA(5)   
            
            constACa=(2.0_dp*(rhopolin(i,tA)*vsol)*(xCa(i)/vCa))/(K0aAA(4)*(xsol(i)**deltavAA(4))) 
            constAMg=(2.0_dp*(rhopolin(i,tA)*vsol)*(xMg(i)/vMg))/(K0aAA(6)*(xsol(i)**deltavAA(6))) 
            constA=constACa+constAMg
            
            qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

            fdisA(i,1)  = 1.0_dp/qAD                             ! A-  
            fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
            fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
            fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
            fdisA(i,5)  = (fdisA(i,1)**2)*constACa               ! A2Ca 
            fdisA(i,6)  = fdisA(i,1)*xA(5)                       ! AMg+ 
            fdisA(i,7)  = (fdisA(i,1)**2)*constAMg               ! A2Mg 
        enddo  

        do t=1,nsegtypes
      
            if(ismonomer_chargeable(t)) then    
                do i=1,n  
                    lbr=lb/epsfcn(i)
                    expborn    = -born(lbr,bornrad%pol,-1)+ expEtotself(i)*vpolAA(1)*vsol 
                    expsqrgrad = expsqrgradpsi(i)*vpolAA(1)      ! no mutipilcation with vsol because defintion constqE
                    !exppi(i,t) = (xsol(i)**vpolAA(1))*exp(psi(i)+expsqrgrad+expborn) /fdisA(i,1)   ! auxilary variable palpha
                    lnexppi(i,t) = log(xsol(i))*vpolAA(1)+psi(i)+expsqrgrad+expborn -log(fdisA(i,1))   ! auxilary variable palpha
                    
                    fdis(i,t)  = fdisA(i,1) 
                enddo  
            else

                do i=1,n
                    expborn    = expEtotself(i)*vpol(t)*vsol 
                    expsqrgrad = expsqrgradpsi(i)*vpol(t)
                    !exppi(i,t) = (xsol(i)**vpol(t))*exp(expsqrgrad+expborn)
                    lnexppi(i,t) = log(xsol(i))*vpol(t)+expsqrgrad+expborn
                    
                    fdis(i,t)  = 0.0_dp
                enddo  
            endif   
        enddo      
      
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes  
                if(isrhoselfconsistent(t)) call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal energy

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
             
        do c=1,cuantas         ! loop over cuantas
            !pro=1.0_dp         ! initial weight conformation (1 or 0)
            lnpro=logweightchain(c)
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo    
            pro=exp(lnpro-lnproshift)
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo

         
        !   .. import results 

        if (rank==0) then 

            q=0.0_dp
            q=local_q

            do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo 

            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
 
            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t) ! polymer density 
                    enddo
                enddo
            enddo     

          
            !     .. construction of fcn and volume fraction polymer             
            rhopol0=(1.0_dp/volcell)/q ! volume polymer segment per volume cell

            do t=1, nsegtypes
                if(ismonomer_chargeable(t)) then 

                    do i=1,n
                        rhopol(i,t) = rhopol0 * rhopol(i,t)               ! density polymer of type t 
                        rhoqpol(i)  = rhoqpol(i) + (- fdisA(i,1)+fdisA(i,4)+fdisA(i,6) )*rhopol(i,t)*vsol 
                        do k=1,4               ! polymer volume fraction
                            xpol(i)=xpol(i)+rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol   
                        enddo
                        xpol(i)=xpol(i)+rhopol(i,t)*(fdisA(i,5)*vpolAA(5)/2.0_dp + &
                                                 fdisA(i,6)*vpolAA(6) + &
                                                 fdisA(i,7)*vpolAA(7)/2.0_dp)*vsol
                    
                    enddo
                else  
                    do i=1,n
                        rhopol(i,t) = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                        xpol(i)     = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                    enddo
                endif          
            enddo    


           !  .. self-consistent equation of densities
            do i=1,n 
                f(i+k2) = phi(i)-xpol(i)
                f(i+k3) = rhopolAA(i)  -fdisA(i,1)*rhopol(i,tA)
                f(i+k4) = rhopolACa(i) -fdisA(i,4)*rhopol(i,tA) 
                f(i+k5) = rhopolAMg(i) -fdisA(i,6)*rhopol(i,tA) 
            enddo    
            count_sc=0    
            do t=1,nsegtypes
                if(isrhoselfconsistent(t)) then
                    count_sc=count_sc+1 
                    k=count_sc*n +k5 
                    do i=1,n   
                        f(i+k)  = rhopol(i,t) - rhopolin(i,t) 
                    enddo
                endif        
            enddo
           
            ! .. packing contraint and total charge 
            do i=1,n
                f(i)    = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)-1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb ! total charge density in units of vsol 
            enddo
          
            !  .. end computation polymer and charge density 

            ! .. electrostatics 
            

            call Poisson_Equation_Eps(f,psi,rhoq,epsfcn) 
            

            norm=l2norm(f,neqint)
            iter=iter+1

            print*,'iter=', iter ,'norm=',norm
            

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

        endif

    end subroutine fcnbrushborn



    subroutine fcnelect(x,f,nn)

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface 
        use vectornorm
        use myutils
        use Poisson

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn

        !     .. declare local variables

        real(dp) :: lnexppiA(nsize),lnexppiB(nsize)  ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize,2)   
        real(dp) :: rhopol_local(nsize,2)
        real(dp) :: q_local
        real(dp) :: xA(3),xB(3),sumxA,sumxB, sgxA, sgxB, qAD, qBD 
        real(dp) :: constA,constB
        real(dp) :: pro,rhopol0,lnpro
        integer :: n,ix,iy,iz,neq_bc                  
        integer :: i,j,k,kL,kR,c,s        ! dummy indices
        real(dp) :: norm
        integer :: conf               ! counts number of conformations
        real(dp) :: cn                ! auxilary variable for Poisson Eq
        character(len=lenText) :: text, istr, rstr
        integer,  parameter  :: A=1, B=2
        real(dp) :: locallnproshift(2), globallnproshift(2)
        
        !     .. executable statements 
        !     .. communication between processors
        
        if (rank.eq.0) then
            flag_solver = 1      !  continue program
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1,   MPI_INTEGER, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x,neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo    
        endif
    
        n=nsize                      ! size vector neq=4*nsizeß x=(pi,psi,rhopolA,rhopolB]
    
        do i=1,n                     ! init x 
            xsol(i)= x(i)            ! solvent volume fraction 
            psi(i) = x(i+n)          ! potential
            rhopolin(i,A)=x(i+2*n)
            rhopolin(i,B)=x(i+3*n)
        enddo
   
        if(rank.eq.0) then 
            neq_bc=0
            if(bcflag(RIGHT)/="cc") then
                neq_bc=nx*ny
                do i=1,neq_bc
                    psiSurfR(i) =x(2*n+i)                 ! surface potentail
                enddo
            endif   
            if(bcflag(LEFT)/="cc") then 
                do i=1,nx*ny
                    psiSurfL(i) =x(2*n+neq_bc+i)          ! surface potentail
                enddo
                neq_bc=neq_bc+nx*ny
            endif    
        endif

        do i=1,n                     ! init volume fractions 
            rhopol_local(i,A) = 0.0_dp     ! A polymer density 
            rhopol_local(i,B) = 0.0_dp     ! B polymer density  
           
            xNa(i)   = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))      ! OH-  volume fraction 
       
            xA(1)= xHplus(i)/(K0aA(1)*(xsol(i)**deltavA(1)))      ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0aA(2)*(xsol(i)**deltavA(2)))   ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0aA(3)*(xsol(i)**deltavA(3)))   ! ACa+/A-
       
            sgxA=1.0_dp+xA(1)+xA(2)+xA(3)                                                         
            constA=(2.0_dp*(rhopolin(i,A)*vsol)*(xCa(i)/vCa))/(K0aA(4)*(xsol(i)**deltavA(4))) 
            qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

            fdisA(i,1)  = 1.0_dp/qAD   ! removed double minus sign 
            fdisA(i,5)  = (fdisA(i,1)**2)*constA
            fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
            fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
            fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
       
            xB(1)= xHplus(i)/(K0aB(1)*(xsol(i) **deltavB(1)))     ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0aB(2)*(xsol(i)**deltavB(2)))   ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0aB(3)*(xsol(i)**deltavB(3)))   ! BCa+/B-
       
            sgxB=xB(1)+xB(2)+xB(3)
            constB=(2.0_dp*(rhopolin(i,B)*vsol)*(xCa(i)/vCa))/(K0aB(4)*(xsol(i)**deltavB(4)))
            qBD = (sgxB+sqrt(sgxB*sgxB+4.0_dp*constB))/2.0_dp  ! remove minus

            fdisB(i,1)  = 1.0_dp/qBD
            fdisB(i,5)  = (fdisB(i,1)**2)*constB
            fdisB(i,2)  = fdisB(i,1)*xB(1)                      ! BH 
            fdisB(i,3)  = fdisB(i,1)*xB(2)                      ! BNa 
            fdisB(i,4)  = fdisB(i,1)*xB(3)                      ! BCa+ 
       
            ! use A^- reference state

            lnexppiA(i)=log(xsol(i))*vpolA(1)-zpolA(1)*psi(i)-log(fdisA(i,1)) ! auxiliary variable
            lnexppiB(i)=log(xsol(i))*vpolB(1)-zpolB(1)*psi(i)-log(fdisB(i,1)) ! auxiliary variable

        enddo

        if(rank==0) then            ! global polymer density
            do i=1,n
                xpol(i) = 0.0_dp 
                rhopol(i,A) = 0.0_dp     ! A polymer density 
                rhopol(i,B) = 0.0_dp     ! B polymer density  
            enddo
            
            q=0.0_dp
            
        endif

        !     .. computation polymer volume fraction 
       
        q_local=0.0_dp       ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal energy

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                if(isAmonomer(s)) then ! A segment 
                    lnpro = lnpro+lnexppiA(k)      
                else
                    lnpro = lnpro+lnexppiB(k)
                endif

            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)

        do c=1,cuantas               ! loop over cuantas
            
            lnpro=logweightchain(c)  

            do s=1,nseg              ! loop over segments 
                k=indexchain(s,c)           
                if(isAmonomer(s)) then ! A segment 
                    lnpro = lnpro+lnexppiA(k)
                else
                    lnpro = lnpro+lnexppiB(k)
                endif
            enddo

            pro=exp(lnpro-lnproshift)
            q_local = q_local+pro

            do s=1,nseg
                k=indexchain(s,c)
                if(isAmonomer(s)) then ! A segment 
                    rhopol_local(k,A)=rhopol_local(k,A)+pro
                else
                    rhopol_local(k,B)=rhopol_local(k,B)+pro
                endif
            enddo
      
            
        enddo   ! end cuantas loop 

       

        !     .. import results

        if (rank==0) then

            q=0.0_dp
            q=q_local

            do i=1, size-1
                source = i
                call MPI_RECV(q_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+q_local
            enddo 
        
            ! first graft point 
            do i=1,n
                rhopol(i,A)=rhopol_local(i,A) ! polymer density 
                rhopol(i,B)=rhopol_local(i,B) ! polymer density 
            enddo
           
            do i=1, size-1
                source = i
                
                call MPI_RECV(rhopol_local(:,A), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(rhopol_local(:,B), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                do k=1,nsize
                    rhopol(k,A)=rhopol(k,A)+rhopol_local(k,A) ! polymer density 
                    rhopol(k,B)=rhopol(k,B)+rhopol_local(k,B) ! polymer
                enddo
            enddo
    
            !  .. construction of fcn and volume fraction polymer        
            rhopol0=(1.0_dp/volcell)/q! volume polymer segment per volume cell

            do i=1,n

                rhopol(i,A) = rhopol0*rhopol(i,A)   ! /deltaG(i) not nessesary since deltaG=1
                rhopol(i,A) = rhopol(i,A) 
                rhopol(i,B) = rhopol0*rhopol(i,B)
                rhopol(i,B)  = rhopol(i,B) 
            
                do k=1,4               ! polymer volume fraction
                    xpol(i)=xpol(i)+rhopol(i,A)*fdisA(i,k)*vpolA(k)*vsol  & 
                        +rhopol(i,B)*fdisB(i,k)*vpolB(k)*vsol
                enddo    

                xpol(i)=xpol(i)+rhopol(i,A)*(fdisA(i,5)*vpolA(5)*vsol/2.0_dp)
                xpol(i)=xpol(i)+rhopol(i,B)*(fdisB(i,5)*vpolB(5)*vsol/2.0_dp)
       
                f(i)=xpol(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp

                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    ((zpolA(1)*fdisA(i,1)+ zpolA(4)*fdisA(i,4))*rhopol(i,A) + &
                     (zpolB(1)*fdisB(i,1)+ zpolB(4)*fdisB(i,4))*rhopol(i,B) )*vsol         
       
                   !  .. total charge density in units of vsol
            enddo  !  .. end computation polymer density and charge density  

            ! .. electrostatics 

            
            ! .. Poisson Eq 
            call Poisson_Equation(f,psi,rhoq)
 
            do i=1,n
                f(2*n+i)=rhopol(i,A)-rhopolin(i,A)
                f(3*n+i)=rhopol(i,B)-rhopolin(i,B)
            enddo

            norm=l2norm(f,4*n)
            iter=iter+1
        
            write(rstr,'(E25.16)')norm
            write(istr,'(I6)')iter
            text="iter = "//trim(istr)//" fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)


        else          ! Export results

            dest = 0
            call MPI_SEND(q_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(rhopol_local(:,A), nsize, MPI_DOUBLE_PRECISION,dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(rhopol_local(:,B), nsize, MPI_DOUBLE_PRECISION,dest, tag, MPI_COMM_WORLD, ierr)
        
        endif


    end subroutine fcnelect


    subroutine fcnneutral(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field
        use vectornorm
        use VdW, only : VdW_contribution_lnexp
       
        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t  ! dummy indices
        real(dp) :: norm
        real(dp) :: rhopol0 !integra_q
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)
        integer  :: count_scf

        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER,dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x, neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        endif

        n=nsize
        
        ! read out x 
        do i=1,n                     
            xsol(i) = x(i)        ! volume fraction solvent
            xpol(i) = 0.0_dp      ! volume fraction polymer
        enddo 

        count_scf=0   
        do t=1,nsegtypes
            if(isrhoselfconsistent(t)) then
                count_scf=count_scf+1 
                k=count_scf*n
                do i=1,n                         
                    rhopolin(i,t) = x(i+k)      ! density 
                enddo
            endif        
        enddo

        do t=1,nsegtypes
            do i=1,n 
                rhopol(i,t) = 0.0_dp      ! .. assign global and local polymer density 
                local_rhopol(i,t) = 0.0_dp
                fdis(i,t)  = 0.0_dp
                lnexppi(i,t) = log(xsol(i))*vpol(t)  
            enddo
        enddo      
       
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes 
                if(isrhoselfconsistent(t)) call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)   
            enddo
        endif 

        !  .. computation polymer volume fraction      

        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal energy

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
             
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)        ! initial weight conformation (1 or 0)
            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo    
            pro=exp(lnpro-lnproshift)
            local_q = local_q+pro
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo
         
        !   .. import results 

        if (rank==0) then 
          
            q=0.0_dp
            q=local_q

            do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo 

            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo
           
            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t) ! polymer density 
                    enddo
                enddo
            enddo     
                
            !     .. construction of fcn and volume fraction polymer             
            rhopol0=(1.0_dp/volcell)/q ! volume polymer segment per volume cell

            do t=1, nsegtypes
                do i=1,n
                    rhopol(i,t) = rhopol0 * rhopol(i,t)       ! density polymer of type t  
                    xpol(i)     = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                    ! .. scf eq for density
                    !f(i+t*n)    = rhopol(i,t) - rhopolin(i,t) 
                enddo
            enddo   

            ! .. scf eq for density       
            count_scf=0    
            do t=1,nsegtypes
                if(isrhoselfconsistent(t)) then
                    count_scf=count_scf+1 
                    k=count_scf*n 
                    do i=1,n   
                        f(i+k)  = rhopol(i,t) - rhopolin(i,t) 
                    enddo
                endif        
            enddo


            do i=1,n
                f(i) = xpol(i)+xsol(i)-1.0_dp
            enddo
          
            !     .. end computation polymer density 

            norm=l2norm(f,neqint)
            iter=iter+1

            print*,'iter=', iter ,'norm=',norm

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo


        endif

    end subroutine fcnneutral


    subroutine fcnneutralnoVdW(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field, only : xpol, xsol, rhopol, q, lnproshift
        use vectornorm
       
        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,l,c,s,ln,t   ! dummy indices
        real(dp) :: norm
        real(dp) :: rhopol0 !integra_q
        integer  :: noffset
        real(dp) :: locallnproshift(2), globallnproshift(2)
        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER,dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x, neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        endif

        n=nsize
        
        ! read out x 
        do i=1,n                     
            xsol(i) = x(i)        ! volume fraction solvent
            xpol(i) = 0.0_dp      ! volume fraction polymer  
        enddo 

        do t=1,nsegtypes
            k=t*n
            do i=1,n 
                rhopol(i,t)=0.0_dp      ! .. assign global and local polymer density 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo

        do t=1,nsegtypes
            do i=1,n
                !fdis(i,t)  = 0.0_dp
                lnexppi(i,t) = log(xsol(i))*vpol(t)  
            enddo
        enddo      

        !  .. computation polymer volume fraction      

        local_q = 0.0_dp    ! init q
             
        lnpro=0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal energy

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
        enddo

        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank
    
        !print*,"rank ",rank ," has local lnproshift value of", locallnproshift(1)
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
    
        !if (rank == 0) then  
        !     print*,"Rank ",globallnproshift(2)," has lowest value of", globallnproshift(1)  
        !endif            

        lnproshift=globallnproshift(1)

        do c=1,cuantas         ! loop over cuantas

            lnpro=logweightchain(c)        ! internal energy

            do s=1,nseg        ! loop over segments 
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro +lnexppi(k,t)
            enddo   
            pro=exp(lnpro-lnproshift)
            local_q = local_q+pro
            
            do s=1,nseg
                k=indexchain(s,c) 
                t=type_of_monomer(s)
                local_rhopol(k,t)=local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
            enddo
        enddo
         
        !   .. import results 

        if (rank==0) then 
          
            q=0.0_dp 
            q=local_q
 
            do i=1, size-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                q=q+local_q
            enddo 
           
           ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t) ! polymer density 
                enddo
            enddo

            do i=1, size-1
                source = i
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t) ! polymer density 
                    enddo
                enddo
            enddo     
                
            !     .. construction of fcn and volume fraction polymer             
            rhopol0=(1.0_dp/volcell)/q  ! volume polymer segment per volume cell

            do t=1, nsegtypes
                do i=1,n
                    rhopol(i,t) = rhopol0 * rhopol(i,t)       ! density polymer of type t  
                    xpol(i)     = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                enddo
            enddo        

            do i=1,n
                f(i) = xpol(i)+xsol(i)-1.0_dp
            enddo
          
            !     .. end computation polymer density 

            norm=l2norm(f,neqint)
            iter=iter+1

            print*,'iter=', iter ,'norm=',norm

        else                      ! Export results 
            
            dest = 0 
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

        endif

    end subroutine fcnneutralnoVdW


    !  .. function solves for bulk volume fraction 

    subroutine fcnbulk(x,f,nn)   

        !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use physconst
        use vectornorm

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out):: f(neq)


        !     .. local variables

        real(dp) :: phiNaCl,phiNa,phiCl,phiK,phiKCl
        real(dp) :: deltavolNaCl,deltavolKCl,norm

        !     .. executable statements 


        deltavolNaCl=(vNaCl-vNa-vCl)
        deltavolKCl=(vKCl-vNa-vCl)

        phiNa   =x(1)
        phiCl   =x(2)
        phiNaCl =x(3)
        phiK    =x(4)
        phiKCl  =x(5)

        !     .. condensation equilibrium eq NaCl

        f(1)= phiNaCl-(phiNa*phiCl*K0ionNa*vsol*vNaCl/(vNa*vCl*vsol))*      & 
         ((1.0_dp -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolNaCl)

        !     .. charge neutrality
        f(2)=phiNa/vNa-phiCl/vCl+2.0_dp*xbulk%Ca/vCa+xbulk%Hplus-xbulk%OHmin+phiK/vK

        !     .. conservation of number NaCl

        f(3)=phiNa/vNa+phiCl/vCl+2.0_dp*phiNaCl/vNaCl + phiKCL/vKCL &    
             -2.0_dp*vsol*(Na/1.0e24_dp)*cNaCl-dabs(xbulk%Hplus-xbulk%OHmin) &
            -2.0_dp*xbulk%Ca/vCa-vsol*(Na/1.0e24_dp)*cKCl

        !     .. condensation equilibrium eq KCl

        f(4)= phiKCl-(phiK*phiCl*K0ionK*vsol*vKCl/(vK*vCl*vsol))* &
             ((1.0_dp -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolKCl)

        !     .. conservation of number K
        f(5)= phiK/vK+phiKCl/vKCl-vsol*(Na/1.0e24_dp)*cKCl

    !    norm=l2norm(f,5)
        iter=iter+1
     
    !    print*,'iter=', iter ,'norm=',norm

    end subroutine fcnbulk


    !     set constrains on vector x  depending on systype value

    subroutine set_contraints(constr)
    
        use precision_definition
        use globals, only : systype, neq , nsize,  LEFT, RIGHT, bcflag
        use volume, only : nx,ny 


        implicit none
            
        real(dp), intent(inout):: constr(:)

        integer :: i, neqint 

        neqint=int(neq,kind(neqint))     ! explict conversion from integer(8) to integer
    
        select case (systype)
            case ("brush_mul","brushdna","nucl_ionbin","nucl_ionbin_sv")      ! multi copolymer:
                do i=1,neqint
                    constr(i)=1.0_dp
                enddo
                do i=1,nsize                    
                    constr(i+nsize)=0.0_dp     ! electrostatic potential 
                enddo  
               
            case ("brushborn")                 
                do i=1,neqint
                    constr(i)=1.0_dp
                enddo
                do i=1,nsize                    
                    constr(i+nsize)=0.0_dp     ! electrostatic potential 
                enddo     

            case ("elect")                     ! AB copolymer: weak acid A weak acid B

                do i=1,nsize                    
                    constr(i)=1.0_dp           ! solvent volume fraction   
                    constr(i+nsize)=0.0_dp     ! electrostatic potential
                    constr(i+2*nsize)=1.0_dp   ! number density A
                    constr(i+3*nsize)=1.0_dp   ! number density B 
                enddo  

            case ("neutral","neutralnoVdW","nucl_neutral_sv")    ! neutral polymers
            
                do i=1,nsize
                    constr(i)=1.0_dp
                enddo 
            
            case ("brush_mulnoVdW")            ! multi copolymer:
                do i=1,neqint
                    constr(i)=1.0_dp
                enddo
                do i=1,nsize                    
                    constr(i+nsize)=0.0_dp     ! electrostatic potential 
                enddo  
            
            case default
                print*,"Error in call to set_constraints subroutine"    
                print*,"Wrong value systype : ", systype
                
                do i=1,neqint
                    constr(i)=1.0_dp
                enddo
                 
        end select  

    end subroutine set_contraints


    ! selects correct fcn function 

    subroutine set_fcn

        use globals, only : systype
        use fcnpointer

        select case (systype) 
        case ("brush_mul")
            fcnptr => fcnelectbrushmulti ! acid and base : no counterion binding VdW
        case ("brush_mulnoVdW")
            fcnptr => fcnelectbrushmultinoVdW ! acid and base : no counterion binding VdW
        case ("brushdna")
            fcnptr => fcnbrushdna   
        case ("nucl_ionbin")
            fcnptr => fcnnucl_ionbin  
        case ("nucl_ionbin_sv")
            fcnptr => fcnnucl_ionbin_sv  
         case ("nucl_neutral_sv")
            fcnptr => fcnnucl_neutral_sv    
        case ("brushborn")
            fcnptr => fcnbrushborn    
        case ("elect")                  ! copolymer weak polyacid, no VdW
             fcnptr => fcnelect    
        case ("neutral")                ! copolymer neutral
            fcnptr => fcnneutral
        case ("neutralnoVdW")           ! homopolymer neutral
            fcnptr => fcnneutralnoVdW
        case ("bulk water")             ! determines compositon bulk electrolyte solution
             fcnptr => fcnbulk
        case default
            print*,"Error in call to set_fcn subroutine"    
            print*,"Wrong value systype : ", systype
            stop
        end select  
        
    end subroutine set_fcn


end module listfcn
