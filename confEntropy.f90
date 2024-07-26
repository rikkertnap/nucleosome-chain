! --------------------------------------------------------------|
! ConfEntropy.f90:                                              |
! constructs the free energy Conformational Entropy             |
! beta Fconf = sigma \sum_alpha P(\alpha)lnP(\alpha)            |
! It also computes other scructure averages involving P(\alpha) |
! namely Rg,Rend_to_end distance, bond and dihedral angles      |
! --------------------------------------------------------------|


module conform_entropy

    use mpivars
    use precision_definition
    implicit none

    private                    ! default all routines in this module private 
    public  ::  FEconf_entropy ! only this subroutine  public
   
contains

    subroutine FEconf_entropy(FEconf,Econf)

        use globals, only : systype
        use myutils, only : lenTEXT, print_to_log, LogUnit
    
        real(dp), intent(out) :: FEconf 
        real(dp), intent(out) :: Econf

        character(len=lenText) :: text 

        select case (systype) 
        case ("elect")
            call FEconf_elect(FEconf,Econf)
        case ("neutral")
            call FEconf_neutral(FEconf,Econf)
        case ("neutralnoVdW")
            call FEconf_neutral_noVdW(FEconf,Econf)
        case ("brush_mul","brushdna","nucl_ionbin")
            call FEconf_brush_mul(FEconf,Econf)
        case ("brush_mulnoVdW")
            call FEconf_brush_mulnoVdW(FEconf,Econf)
        case ("brushborn")
            call FEconf_brush_born(FEconf,Econf) 
        case ("nucl_ionbin_sv")
            call FEconf_nucl_ionbin_sv(FEconf,Econf)
        case ("nucl_ionbin_Mg")
            call FEconf_nucl_ionbin_Mg(FEconf,Econf)
        case ("nucl_ionbin_MgA")
            call FEconf_nucl_ionbin_MgA(FEconf,Econf)
        case ("nucl_neutral_sv")
            call FEconf_nucl_neutral_sv(FEconf,Econf)
        case default
            text="FEconf_entropy: wrong systype: "//systype//"stopping program"
            call print_to_log(LogUnit,text)
            stop
        end select   


    end subroutine FEconf_entropy

   ! computes conformational entropy in neutral state 

    subroutine FEconf_neutral(FEconf,Econf)
    
        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, logweightchain 
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, gyr_tensor, avgyr_tensor
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, avnucl_spacing       
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol, rhopol, q, lnproshift
        use parameters, only : vpol, isVdW, isrhoselfconsistent, write_Palpha
        use myutils, only : newunit, lenText

        real(dp), intent(out) :: FEconf,Econf
       
        ! .. declare local variables

        real(dp) :: lnexppi(nsize,nsegtypes)   ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k             ! dummy indices
        real(dp) :: FEconf_local,Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1) 
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname, istr

        !  .. opens file to save Palpha 
        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            do t=1,nsegtypes
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
        endif    

        !     .. executable statements 

        do t=1,nsegtypes
            do i=1,nsize
                lnexppi(i,t) = log(xsol(i))*vpol(t)
            enddo    
        enddo      
       

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init 
        Econf_local=0.0_dp   
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
            
        do c=1,cuantas         ! loop over cuantas
            if(no_overlapchain(c)) then 

                lnpro=logweightchain(c) -energychainLJ(c)      
                do s=1,nseg        ! loop over segments                     
                    k=indexchain(s,c)
                    t=type_of_monomer(s)                
                    lnpro = lnpro+lnexppi(k,t)
                enddo  
                pro=exp(lnpro-lnproshift)  

                FEconf_local = FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+energychainLJ(c)*pro
                Rgsqr_local = Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local = Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local= bond_angle_local +bond_angle(:,c)*pro  
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local +nucl_spacing(:,c)*pro         
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro
                
                if(write_Palpha) write(un,*)pro/q
            endif   

        enddo
        
        Econf_local=Econf_local/q
        Rgsqr_local = Rgsqr_local/q
        Rendsqr_local = Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q 
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q

        ! communicate local quantities

        if(rank==0) then

            ! normalize
            
            FEconf = FEconf_local
            Econf = Econf_local
            avRgsqr = Rgsqr_local
            avRendsqr = Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local 
            avnucl_spacing = nucl_spacing_local 
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local
            
            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
               
                FEconf=FEconf+FEconf_local
                Econf =Econf +Econf_local
                avRgsqr=avRgsqr+Rgsqr_local
                avRendsqr=avRendsqr+Rendsqr_local
                avbond_angle = avbond_angle+bond_angle_local
                avdihedral_angle =  avdihedral_angle + dihedral_angle_local
                avnucl_spacing = avnucl_spacing + nucl_spacing_local 
                avgyr_tensor = avgyr_tensor + gyr_tensor_local
                avAsphparam = avAsphparam + Asphparam_local 

            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(bond_angle_local, nangles , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(nucl_spacing_local,nbonds, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif

        if(write_Palpha) close(un)


    end subroutine FEconf_neutral


    ! computes conformational entropy in neutral state 
     
    subroutine FEconf_neutral_noVdW(FEconf,Econf) 
    
        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, logweightchain
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr,  nucl_spacing, avnucl_spacing
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, gyr_tensor, avgyr_tensor
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol, rhopol, q, lnproshift
        use parameters, only : vpol, isVdW, VdWscale, write_Palpha 
        use myutils, only : newunit, lenText

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3) 
        real(dp) :: nucl_spacing_local(nnucl-1) 
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname, istr
        
        !  .. opens file to save Palpha 
        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif   

       ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
        endif    

        !     .. executable statements 

        do t=1,nsegtypes
            do i=1,nsize
                lnexppi(i,t) = log(xsol(i))*vpol(t)
            enddo    
        enddo      

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init Econf
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp 
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
            
        do c=1,cuantas         ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c) -energychainLJ(c)         ! internal energy  
                do s=1,nseg        ! loop over segments                     
                    k=indexchain(s,c)
                    t=type_of_monomer(s)                
                    lnpro = lnpro+ lnexppi(k,t)
                enddo    
                pro=exp(lnpro-lnproshift)

                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local = Rgsqr_local + Rgsqr(c)*pro
                Rendsqr_local = Rendsqr_local + Rendsqr(c)*pro
                bond_angle_local = bond_angle_local + bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local + dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local + nucl_spacing(:,c)*pro 
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro

               if(write_Palpha) write(un,*)pro/q
            endif   
        enddo

        Econf_local=Econf_local/q
        Rgsqr_local = Rgsqr_local/q
        Rendsqr_local = Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q 

        ! communicate FEconf

        if(rank==0) then
            ! normalize
           
            FEconf=FEconf_local
            Econf=Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local 
            avnucl_spacing = nucl_spacing_local
            avgyr_tensor = gyr_tensor_local 
            avAsphparam = Asphparam_local

            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local,  1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local,   1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local,   1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                if(nangles>=1) then  
                    call MPI_RECV(bond_angle_local, nangles,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif
                if(ndihedrals>=1) then
                    call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif
                if(nbonds>=1) then 
                    call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf+FEconf_local
                Econf =Econf +Econf_local
                avRgsqr=avRgsqr+Rgsqr_local   
                avRendsqr=avRendsqr+Rendsqr_local
                avbond_angle = avbond_angle +bond_angle_local
                avdihedral_angle = avdihedral_angle+ dihedral_angle_local 
                avnucl_spacing = avnucl_spacing+ nucl_spacing_local 
                avgyr_tensor = avgyr_tensor + gyr_tensor_local
                avAsphparam = avAsphparam + Asphparam_local

            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local,  1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local,   1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local,   1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            if(nangles>=1) then 
                call MPI_SEND(bond_angle_local, nangles,       MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(ndihedrals>=1) then 
                call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(nbonds>=1) then
                call MPI_SEND(nucl_spacing_local, nbonds ,     MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif 
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif

        if(write_Palpha) close(un)

    end subroutine FEconf_neutral_noVdW


    subroutine FEconf_brush_mul(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl,nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, avnucl_spacing
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, gyr_tensor, avgyr_tensor
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol,psi, fdis,rhopol,q, lnproshift
        use parameters, only : vpol, zpol, isVdW, isrhoselfconsistent, write_Palpha
        use myutils, only : lenText, newunit

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1) 
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=LenText) :: fname 

        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif   

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

        !     .. executable statements 

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,nsize                                              
                    lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i)-log(fdis(i,t))   ! auxilary variable palpha
                enddo  
            else
                do i=1,nsize
                     lnexppi(i,t) = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp   !init Econf
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
         
        do c=1,cuantas         ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c) -energychainLJ(c)    
                do s=1,nseg        ! loop over segments                     
                    k=indexchain(s,c)
                    t=type_of_monomer(s)                
                    lnpro = lnpro+lnexppi(k,t)
                enddo 
                pro=exp(lnpro-lnproshift)      
                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local =Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local = bond_angle_local +bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local+nucl_spacing(:,c)*pro
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro
           
                if(write_Palpha) write(un,*)pro/q
            endif    

        enddo

        Econf_local=Econf_local/q
        Rgsqr_local=Rgsqr_local/q
        Rendsqr_local=Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q 
        gyr_tensor_local = gyr_tensor_local/q 
        Asphparam_local = Asphparam_local/q

        ! communicate 

        if(rank==0) then

            ! normalize

            FEconf=FEconf_local
            Econf =Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local
            avnucl_spacing = nucl_spacing_local
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local

            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf+FEconf_local
                Econf =Econf +Econf_local             
                avRgsqr=avRgsqr+Rgsqr_local   
                avRendsqr=avRendsqr+Rendsqr_local
                avbond_angle = avbond_angle+bond_angle_local
                avdihedral_angle =  avdihedral_angle +dihedral_angle_local
                avnucl_spacing = avnucl_spacing+ nucl_spacing_local 
                avgyr_tensor = avgyr_tensor + gyr_tensor_local
                avAsphparam = avAsphparam + Asphparam_local
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(bond_angle_local, nangles, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(nucl_spacing_local,nbonds, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(write_Palpha) close(un)

    end subroutine FEconf_brush_mul


    subroutine FEconf_brush_mulnoVdW(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, avnucl_spacing
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, gyr_tensor, avgyr_tensor
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol, psi, fdis, rhopol, q ,lnproshift
        use parameters, only : vpol, zpol, write_Palpha
        use myutils, only : lenText, newunit
        
        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,c,s,k       ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1) 
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local 
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname


        if(write_Palpha) then
             call make_filename_Palpha(fname,rank)
             open(unit=newunit(un),file=fname)
        endif

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

        !     .. executable statements 

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,nsize                                              
                    lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i)-log(fdis(i,t))   ! auxilary variable palpha
                enddo  
            else
                do i=1,nsize
                     lnexppi(i,t) = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      
       
    
        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp 
        Asphparam_local = 0.0_dp
            
        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
        

        do c=1,cuantas         ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c) -energychainLJ(c)      ! internal energy  
                do s=1,nseg        ! loop over segments                     
                    k=indexchain(s,c)
                    t=type_of_monomer(s)                
                    lnpro = lnpro+ lnexppi(k,t)
                enddo    
                pro=exp(lnpro-lnproshift)
                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local =Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local = bond_angle_local +bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local+ nucl_spacing(:,c)*pro
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro
           
                if(write_Palpha) write(un,*)pro/q
            endif    
        enddo
        
        Econf_local=Econf_local/q
        Rgsqr_local=Rgsqr_local/q
        Rendsqr_local=Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local= nucl_spacing_local/q
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q

        ! communicate FEconf

        if(rank==0) then
            ! normalize
           
            FEconf =FEconf_local
            Econf  =Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local
            avnucl_spacing = nucl_spacing_local
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local

            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf +FEconf_local
                Econf =Econf + Econf_local      
                avRgsqr=avRgsqr+Rgsqr_local   
                avRendsqr=avRendsqr+Rendsqr_local
                avbond_angle = avbond_angle+bond_angle_local
                avdihedral_angle =  avdihedral_angle +dihedral_angle_local
                avnucl_spacing = avnucl_spacing+ nucl_spacing_local
                avgyr_tensor = avgyr_tensor + gyr_tensor_local
                avAsphparam = avAsphparam + Asphparam_local

            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(bond_angle_local, nangles , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(nucl_spacing_local, nbonds , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        
        endif

        if(write_Palpha) close(un)

    end subroutine FEconf_brush_mulnoVdW


    subroutine FEconf_elect(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl,nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain, isAmonomer
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, avnucl_spacing
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, gyr_tensor, avgyr_tensor
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field,  only : xsol, psi, fdisA,fdisB, rhopol, q ,lnproshift
        use parameters
        use myutils, only : lenText, newunit

        real(dp), intent(out) :: FEconf
        real(dp), intent(out) :: Econf
        
        !     .. declare local variables
        real(dp) :: lnexppiA(nsize),lnexppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer  :: i,k,c,s        ! dummy indices
        real(dp) :: pro,lnpro
        real(dp) :: FEconf_local, Econf_local
        real(dp) :: q_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1) 
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname

        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif   
 
        ! .. executable statements 

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(fdisA(:,1),nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(fdisB(:,1),nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(fdisA(:,1), nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)    
            call MPI_RECV(fdisB(:,1), nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif
            
        do i=1,nsize
              lnexppiA(i)=log(xsol(i))*vpolA(1)-zpolA(1)*psi(i)-log(fdisA(i,1)) ! auxiliary variable
              lnexppiB(i)=log(xsol(i))*vpolB(1)-zpolB(1)*psi(i)-log(fdisB(i,1)) ! auxiliary variable
        enddo
       
    
        FEconf_local= 0.0_dp
        Econf_local=0.0_dp 
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local= 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3

        do c=1,cuantas             ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c) -energychainLJ(c)              ! initial weight conformation 
                do s=1,nseg              ! loop over segments 
                    k=indexchain(s,c)         
                    if(isAmonomer(s)) then ! A segment 
                        lnpro = lnpro+lnexppiA(k)
                    else
                        lnpro = lnpro+lnexppiB(k)
                    endif
                enddo
                pro=exp(lnpro-lnproshift)
                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local =Rendsqr_local+Rendsqr(c)*pro 
                bond_angle_local = bond_angle_local +bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro  
                nucl_spacing_local = nucl_spacing_local +nucl_spacing(:,c)*pro
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c) * pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro

                if(write_Palpha) write(un,*) pro/q   
            endif      
        enddo  

        Econf_local=Econf_local/q
        Rgsqr_local=Rgsqr_local/q
        Rendsqr_local=Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q 

        ! communicate FEconf

        if(rank==0) then
            ! normalize
            
            FEconf=FEconf_local
            Econf=Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local
            avnucl_spacing =  nucl_spacing_local
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local
            
            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf + FEconf_local 
                Econf= Econf  + Econf_local  
                avRgsqr=avRgsqr+Rgsqr_local   
                avRendsqr=avRendsqr+Rendsqr_local  
                avbond_angle = avbond_angle+bond_angle_local
                avdihedral_angle =  avdihedral_angle +dihedral_angle_local
                avnucl_spacing = avnucl_spacing+ nucl_spacing_local
                avgyr_tensor = avgyr_tensor + gyr_tensor_local 
                avAsphparam = avAsphparam + Asphparam_local

            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(bond_angle_local, nangles , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(nucl_spacing_local, nbonds , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif

        Econf=0.0_dp
         
        if(write_Palpha) close(un)

    end subroutine FEconf_elect



    subroutine FEconf_brush_born(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl,nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, avnucl_spacing
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, gyr_tensor, avgyr_tensor 
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol,psi, fdis,rhopol,q, lnproshift, fdisA, epsfcn, Depsfcn
        use field, only : xOHmin,xHplus,xNa,xCl,xMg,xCa,xRb
        use parameters, only : bornrad, lb, VdWscale, tA, isrhoselfconsistent, isVdW, write_Palpha
        use parameters, only : vpolAA, vsol, vNa, vCl, vRb, vMg, vCa ,vpol
        use parameters, only : zNa, zCl, zRb, zMg, zCa, zpolAA        
        use Poisson, only : Poisson_Equation_Eps, Poisson_Equation_Surface_Eps, grad_pot_sqr_eps_cubic
        use dielectric_const, only : dielectfcn, born
        use myutils, only : lenText, newunit

        
        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,c,s,k,tc,l    ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        integer  :: tcfdis(3)
        real(dp) :: rhopolAA(nsize),rhopolACa(nsize), rhopolAMg(nsize)
        real(dp) :: lbr,expborn,Etotself,expsqrgrad, Eself
        real(dp) :: expsqrgradpsi(nsize),expEtotself(nsize)
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1)
        real(dp) :: gyr_tensor_local(3,3) 
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname

        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif   

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        tcfdis(1)=1
        tcfdis(2)=4
        tcfdis(3)=6
        
        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(epsfcn, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(Depsfcn, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,3
                    tc=tcfdis(t)
                    call MPI_SEND(fdisA(:,tc) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo    
                do t=1,nsegtypes
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
            call MPI_RECV(epsfcn, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(Depsfcn, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            do t=1,3
                tc=tcfdis(t)
                call MPI_RECV(fdisA(:,tc) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)    
            enddo
            do t=1,nsegtypes
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

        !     .. executable statements 

        ! gradient potential contribution to PDF
        call grad_pot_sqr_eps_cubic(psi,epsfcn, Depsfcn,expsqrgradpsi)

       
        do i=1,nsize
           
            rhopolAA(i) =fdisA(i,1)*rhopol(i,tA)
            rhopolACa(i)=fdisA(i,4)*rhopol(i,tA) 
            rhopolAMg(i)=fdisA(i,6)*rhopol(i,tA) 


            lbr = lb/epsfcn(i)     ! local Bjerrum length

            !xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-born(lbr,bornrad%Na,zNa)-psi(i)*zNa) ! Na+ volume fraction 
            !xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-born(lbr,bornrad%Cl,zCl)-psi(i)*zCl) ! Cl- volume fraction
            !xHplus(i)  = expmu%Hplus*(xsol(i))  *exp(-born(lbr,bornrad%Hplus,1)-psi(i))    ! H+  volume fraction
            !xOHmin(i)  = expmu%OHmin*(xsol(i))  *exp(-born(lbr,bornrad%OHmin,-1)+psi(i))   ! OH- volume fraction
            !xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-born(lbr,bornrad%Rb,zRb)-psi(i)*zRb) ! Rb+ volume fraction
            !xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-born(lbr,bornrad%Ca,zCa)-psi(i)*zCa) ! Ca++ volume fraction 
            !xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-born(lbr,bornrad%Mg,zMg)-psi(i)*zMg) ! Mg++ volume fraction 



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
        
        enddo 


        do t=1,nsegtypes
      
            if(ismonomer_chargeable(t)) then    
                do i=1,nsize  
                    lbr=lb/epsfcn(i)
                    expborn    = -born(lbr,bornrad%pol,-1)+ expEtotself(i)*vpolAA(1)*vsol 
                    expsqrgrad = expsqrgradpsi(i)*vpolAA(1)      ! no mutipilcation with vsol because defintion constqE
                    !exppi(i,t) = (xsol(i)**vpolAA(1))*exp(psi(i)+expsqrgrad+expborn) /fdisA(i,1)   ! auxilary variable palpha
                    lnexppi(i,t) = log(xsol(i))*vpolAA(1)+psi(i)+expsqrgrad+expborn -log(fdisA(i,1))   ! auxilary variable palpha  
    
                enddo  
            else

                do i=1,nsize
                    expborn    = expEtotself(i)*vpol(t)*vsol 
                    expsqrgrad = expsqrgradpsi(i)*vpol(t)
                    !exppi(i,t) = (xsol(i)**vpol(t))*exp(expsqrgrad+expborn)
                    lnexppi(i,t) = log(xsol(i))*vpol(t)+expsqrgrad+expborn
                    
                enddo  
            endif   
        enddo      

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp
         
        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3 

        do c=1,cuantas         ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c)  -energychainLJ(c)   
                do s=1,nseg        ! loop over segments                     
                    k=indexchain(s,c)
                    t=type_of_monomer(s)                
                    lnpro = lnpro+lnexppi(k,t)
                enddo 
                pro=exp(lnpro-lnproshift)      
                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local =Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local = bond_angle_local+bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local + dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local + nucl_spacing(:,c)*pro
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro

                if(write_Palpha) write(un,*)pro/q
            endif    
        enddo

        Econf_local=Econf_local/q
        Rgsqr_local=Rgsqr_local/q
        Rendsqr_local=Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q

        ! communicate 

        if(rank==0) then

            ! normalize
             
            FEconf=FEconf_local
            Econf =Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local
            avnucl_spacing = nucl_spacing_local
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local 
            
            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr) 
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr) 

                FEconf=FEconf + FEconf_local 
                Econf= Econf  + Econf_local  
                avRgsqr=avRgsqr+Rgsqr_local   
                avRendsqr=avRendsqr+Rendsqr_local   
                avbond_angle = avbond_angle+bond_angle_local
                avdihedral_angle =  avdihedral_angle + dihedral_angle_local
                avnucl_spacing = avnucl_spacing + nucl_spacing_local
                avgyr_tensor = avgyr_tensor + gyr_tensor_local
                avAsphparam = avAsphparam + Asphparam_local 
                           
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(bond_angle_local, nangles, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(nucl_spacing_local, nbonds , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif
       
        if(write_Palpha) close(un)

    end subroutine FEconf_brush_born


    subroutine FEconf_nucl_ionbin_sv(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl,nsegtypes, nsize, cuantas
        use chains, only : indexconf,  nelem, type_of_monomer, ismonomer_chargeable, logweightchain,elem_charge
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, avnucl_spacing, gyr_tensor, avgyr_tensor
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle 
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol,psi, fdis,rhopol,q, lnproshift
        use parameters, only : vnucl, vsol, zpol, isVdW,  isrhoselfconsistent, write_Palpha
        use myutils, only : lenText, newunit

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: i,j, t,g,gn,c,s,k, jcharge       ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1)
        real(dp) :: gyr_tensor_local(3,3) 
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname

        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif   

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

        !     .. executable statements 
        do i=1,nsize   
            lnexppivw(i) = log(xsol(i))/vsol
        enddo
            
        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,nsize                                              
                    lnexppi(i,t) =  -zpol(t,2)*psi(i)-log(fdis(i,t))   ! auxilary variable palpha
                enddo  
            endif   
        enddo      

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp   !init Econf
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3

         
        do c=1,cuantas         ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c) -energychainLJ(c)    
                do s=1,nseg        ! loop over segments                     
                    t = type_of_monomer(s)                
                    do j=1,nelem(s)               ! loop over elements of segment  
                        k = indexconf(s,c)%elem(j)
                        lnpro = lnpro +lnexppivw(k)*vnucl(j,t)   ! excluded-volume contribution        
                    enddo
                    if(ismonomer_chargeable(t)) then
                        jcharge=elem_charge(t)
                        k = indexconf(s,c)%elem(jcharge) 
                        lnpro = lnpro + lnexppi(k,t)  ! electrostatic, VdW and chemical contribution
                    endif
                enddo 

                pro=exp(lnpro-lnproshift)      
                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local +pro*energychainLJ(c)
                Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local =Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local = bond_angle_local +bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local+nucl_spacing(:,c)*pro
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro 
                Asphparam_local = Asphparam_local + Asphparam(c) * pro
                
                if(write_Palpha) write(un,*)pro/q
            endif    
        enddo
        
        Econf_local=Econf_local/q
        Rgsqr_local=Rgsqr_local/q
        Rendsqr_local=Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q 
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q
        
        ! communicate 

        if(rank==0) then

            ! normalize

            FEconf=FEconf_local
            Econf =Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local
            avnucl_spacing = nucl_spacing_local
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local

            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                if(nangles>=1) then 
                    call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                endif
                if(ndihedrals>=1) then
                    call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,&
                        tag,MPI_COMM_WORLD,stat,ierr)
                endif
                if(nbonds>=1) then
                    call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif 
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf+ FEconf_local
                Econf =Econf + Econf_local             
                avRgsqr=avRgsqr+ Rgsqr_local   
                avRendsqr=avRendsqr+ Rendsqr_local
                avbond_angle = avbond_angle+ bond_angle_local
                avdihedral_angle =  avdihedral_angle + dihedral_angle_local
                avnucl_spacing = avnucl_spacing + nucl_spacing_local 
                avgyr_tensor = avgyr_tensor + gyr_tensor_local 
                avAsphparam = avAsphparam + Asphparam_local 

            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            if(nangles>=1) then
                call MPI_SEND(bond_angle_local, nangles, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(ndihedrals>=1) then
                call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(nbonds>=1) then
                call MPI_SEND(nucl_spacing_local,nbonds, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
         endif
         
         if(write_Palpha) close(un)

    end subroutine FEconf_nucl_ionbin_sv


    subroutine FEconf_nucl_ionbin_Mg(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl,nsegtypes, nsize, cuantas
        use parameters, only : ta, Phos
        use volume, only : inverse_indexneighbor_phos
        use chains, only : indexconf,  nelem, type_of_monomer, ismonomer_chargeable, logweightchain,elem_charge
        use chains, only : type_of_charge, elem_charge, indexconfpair, nneigh, maxneigh
        use chains, only : index_phos, inverse_index_phos, len_index_phos
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, avnucl_spacing
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, gyr_tensor, avgyr_tensor 
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol,psi, q, lnproshift
        use field, only : gdisA,gdisB, fdisPP
        use parameters, only : vnucl, vsol, zpol, isVdW,  isrhoselfconsistent, write_Palpha
        ! use VdW, only : VdW_contribution_lnexp
        use myutils, only : lenText, newunit
        
        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: i,j, t,g,gn,c,s,k, jcharge, ind, jj, m, k_ind, mr      ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1) 
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname

        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    if(ismonomer_chargeable(t)) then
                        if(t/=ta) then
                            if(type_of_charge(t)=="A") then
                                call MPI_SEND(gdisA(:,1,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                            else 
                                call MPI_SEND(gdisB(:,2,t), nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                            endif
                        else 
                            do jj=1,maxneigh
                                call MPI_SEND(fdisPP(:,jj,Phos,Phos), len_index_phos , MPI_DOUBLE_PRECISION,&
                                    dest, tag,MPI_COMM_WORLD,ierr)  
                            enddo
                        endif
                    endif            
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes 
                if(ismonomer_chargeable(t)) then
                    if(t/=ta) then
                        if(type_of_charge(t)=="A") then
                            call MPI_RECV(gdisA(:,1,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                        else
                            call MPI_RECV(gdisB(:,2,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)
                        endif
                    else
                        do jj=1,maxneigh
                            call MPI_RECV(fdisPP(:,jj,Phos,Phos) , len_index_phos, MPI_DOUBLE_PRECISION, &
                                source,tag, MPI_COMM_WORLD,stat, ierr)   
                        enddo
                    endif
                endif              
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

        !     .. executable statements 
        do i=1,nsize   
            lnexppivw(i) = log(xsol(i))/vsol
        enddo
            
         
        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    if(type_of_charge(t)=="A") then  !  acid
                        do i=1,nsize
                            lnexppi(i,t) = psi(i) -log(gdisA(i,1,t))      ! auxilary variable palpha log(xsol)*(delta vpol+0) =0 
                        enddo

                    else !  base
                        do i=1,nsize
                            lnexppi(i,t) = -log(gdisB(i,2,t))             ! auxilary variable palpha lo  
                        enddo
                    endif  
                                
                else
                    ! t=ta : phosphate
                                   
                    do ind=1,len_index_phos ! loop over index of  location of phosphates
                        i = index_phos(ind)  ! give the lattice location 
                        lnexppi(i,t) = psi(i)!!   ! auxilary variable palpha
                    enddo
                endif
            endif      
        enddo   

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp   !init Econf
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
 
 
        do c=1,cuantas                            ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c) -energychainLJ(c)
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
                        k_ind = inverse_index_phos(k)

                        do jj=1,nneigh(s,c)           ! loop neighbors 
     
                            m = indexconfpair(s,c)%elem(jj)
                            !mr = inverse_indexneighbor(k,m) ! relative label of index m relative to k
                            mr = inverse_indexneighbor_phos(k_ind,m) 

                            lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*vnucl(1,ta) &
                                    -log(fdisPP(k_ind,mr,Phos,Phos)))/(2.0_dp*nneigh(s,c))
                        enddo    
                    endif        
                enddo    

                pro = exp(lnpro-lnproshift)  

                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local =Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local = bond_angle_local +bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local+nucl_spacing(:,c)*pro
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro
                
                if(write_Palpha) write(un,*)pro/q
            endif
        enddo
        
        Econf_local=Econf_local/q
        Rgsqr_local=Rgsqr_local/q
        Rendsqr_local=Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q 
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q

        ! communicate 

        if(rank==0) then

            ! normalize

            FEconf=FEconf_local
            Econf =Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local
            avnucl_spacing = nucl_spacing_local
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local


            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                if(nangles>=1) then
                    call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                endif
                if(ndihedrals>=1) then
                    call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,&
                    tag,MPI_COMM_WORLD,stat,ierr)
                endif
                if(nbonds>=1) then
                    call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf + FEconf_local
                Econf =Econf + Econf_local             
                avRgsqr=avRgsqr + Rgsqr_local   
                avRendsqr=avRendsqr + Rendsqr_local
                avbond_angle = avbond_angle + bond_angle_local
                avdihedral_angle =  avdihedral_angle + dihedral_angle_local
                avnucl_spacing = avnucl_spacing + nucl_spacing_local 
                avgyr_tensor = avgyr_tensor + gyr_tensor_local 
                avAsphparam = avAsphparam + Asphparam_local

            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            if(nangles>=1) then 
                 call MPI_SEND(bond_angle_local, nangles, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(ndihedrals>=1) then
                call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif 
            if(nbonds>=1) then 
                call MPI_SEND(nucl_spacing_local,nbonds, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif
        
        if(write_Palpha) close(un)

    end subroutine FEconf_nucl_ionbin_Mg

    subroutine FEconf_nucl_ionbin_MgA(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl,nsegtypes, nsize, cuantas
        use parameters, only : ta, Phos
        use volume, only : inverse_indexneighbor_phos
        use chains, only : indexconf,  nelem, type_of_monomer, ismonomer_chargeable, logweightchain,elem_charge
        use chains, only : type_of_charge, elem_charge, indexconfpair, nneigh, maxneigh
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, avnucl_spacing
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, gyr_tensor, avgyr_tensor
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol,psi, q, lnproshift
        use field, only : gdisA,gdisB, fdisPP_loc, fdisP2Mg_loc
        use parameters, only : vnucl, vsol, zpol, isVdW,  isrhoselfconsistent, write_Palpha
        use myutils, only : lenText, newunit
        use modfcnMgexpl, only : compute_fdisPP
        
        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: i,j, t,g,gn,c,s,k, jcharge, ind, jj, m, k_ind, mr      ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1) 
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname

        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    if(ismonomer_chargeable(t)) then
                        if(t/=ta) then
                            if(type_of_charge(t)=="A") then
                                call MPI_SEND(gdisA(:,1,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                            else 
                                call MPI_SEND(gdisB(:,2,t), nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                            endif
                        endif
                    endif            
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes 
                if(ismonomer_chargeable(t)) then
                    if(t/=ta) then
                        if(type_of_charge(t)=="A") then
                            call MPI_RECV(gdisA(:,1,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                        else
                            call MPI_RECV(gdisB(:,2,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)
                        endif
                    endif
                endif              
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

        !     .. executable statements 
        do i=1,nsize   
            lnexppivw(i) = log(xsol(i))/vsol
        enddo
            
         
        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    if(type_of_charge(t)=="A") then  !  acid
                        do i=1,nsize
                            lnexppi(i,t) = psi(i) -log(gdisA(i,1,t))      ! auxilary variable palpha log(xsol)*(delta vpol+0) =0 
                        enddo

                    else !  base
                        do i=1,nsize
                            lnexppi(i,t) = -log(gdisB(i,2,t))             ! auxilary variable palpha lo  
                        enddo
                    endif  
                                
                else
                    ! t=ta : phosphate              
                    do i=1,nsize ! loop over index of  location of phosphates
                        lnexppi(i,t) = psi(i)!!   ! auxilary variable palpha
                    enddo
                endif
            endif      
        enddo   

       
        !if(isVdW) then 
        !    print*,"VdW SCF not allowed in FEconf_nucl_ionbin_Mg"
        !endif 

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp   !init Econf
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
 
 
        do c=1,cuantas                            ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro=logweightchain(c)-energychainLJ(c) 
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
                            call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc,k ,m)

                            lnpro =lnpro +(lnexppi(k,ta)+lnexppi(m,ta)+ (lnexppivw(k) +lnexppivw(m))*vnucl(1,ta) &
                                    -log(fdisPP_loc(Phos,Phos)))/(2.0_dp*nneigh(s,c))
                        enddo    
                    endif        
                enddo    

                pro = exp(lnpro-lnproshift)  

                FEconf_local=FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local =Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local = bond_angle_local +bond_angle(:,c)*pro
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local+nucl_spacing(:,c)*pro
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c) *pro
                Asphparam_local = Asphparam_local + Asphparam(c) * pro

                if(write_Palpha) write(un,*)pro/q
            endif    
        enddo

        Econf_local=Econf_local/q
        Rgsqr_local=Rgsqr_local/q
        Rendsqr_local=Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q 
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q
        ! communicate 

        if(rank==0) then

            ! normalize

            FEconf=FEconf_local
            Econf =Econf_local
            avRgsqr=Rgsqr_local
            avRendsqr=Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local
            avnucl_spacing = nucl_spacing_local
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local

            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                if(nangles>=1) then
                    call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                endif
                if(ndihedrals>=1) then
                    call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,&
                    tag,MPI_COMM_WORLD,stat,ierr)
                endif
                if(nbonds>=1) then
                    call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf + FEconf_local
                Econf =Econf + Econf_local             
                avRgsqr=avRgsqr + Rgsqr_local   
                avRendsqr=avRendsqr + Rendsqr_local
                avbond_angle = avbond_angle + bond_angle_local
                avdihedral_angle =  avdihedral_angle + dihedral_angle_local
                avnucl_spacing = avnucl_spacing + nucl_spacing_local 
                avgyr_tensor = avgyr_tensor + gyr_tensor_local 
                avAsphparam = avAsphparam + Asphparam_local

            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            if(nangles>=1) then 
                 call MPI_SEND(bond_angle_local, nangles, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(ndihedrals>=1) then
                call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif 
            if(nbonds>=1) then 
                call MPI_SEND(nucl_spacing_local,nbonds, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif
        
        if(write_Palpha) close(un)

    end subroutine FEconf_nucl_ionbin_MgA    

    subroutine FEconf_nucl_neutral_sv(FEconf,Econf)
    
        !  .. variables and constant declaractions 

        use globals, only : nseg, nnucl, nsegtypes, nsize, cuantas
        use chains, only : indexconf, nelem, type_of_monomer, logweightchain
        use chains, only : Rgsqr, Rendsqr, avRgsqr, avRendsqr, nucl_spacing, gyr_tensor, avgyr_tensor
        use chains, only : bond_angle, dihedral_angle,avbond_angle, avdihedral_angle, avnucl_spacing       
        use chains, only : Asphparam, avAsphparam, energychainLJ, no_overlapchain
        use field, only : xsol, rhopol, q, lnproshift
        use parameters, only : vsol, vnucl, isVdW, isrhoselfconsistent, write_Palpha
        use myutils, only : lenText, newunit

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables

        real(dp) :: lnexppi(nsize)   ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,c,s,k ,j           ! dummy indices
        real(dp) :: FEconf_local,Econf_local
        real(dp) :: Rgsqr_local,Rendsqr_local
        real(dp) :: bond_angle_local(nnucl-2)
        real(dp) :: dihedral_angle_local(nnucl-3)
        real(dp) :: nucl_spacing_local(nnucl-1)
        real(dp) :: gyr_tensor_local(3,3)
        real(dp) :: Asphparam_local 
        integer  :: nbonds,ndihedrals,nangles
        integer  :: un
        character(len=lenText) :: fname

        if(write_Palpha) then
            call make_filename_Palpha(fname,rank)
            open(unit=newunit(un),file=fname)
        endif

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
                call MPI_SEND(q , 1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            do t=1,nsegtypes
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
            call MPI_RECV(q , 1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
        endif    

        !     .. executable statements 

        do i=1,nsize
            lnexppi(i) = log(xsol(i))/vsol 
        enddo    

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init 
        Econf_local=0.0_dp   
        Rgsqr_local=0.0_dp
        Rendsqr_local=0.0_dp
        bond_angle_local = 0.0_dp
        dihedral_angle_local = 0.0_dp
        nucl_spacing_local = 0.0_dp
        gyr_tensor_local = 0.0_dp
        Asphparam_local = 0.0_dp

        nbonds=nnucl-1
        nangles=nnucl-2
        ndihedrals=nnucl-3
            
        do c=1,cuantas                        ! loop over cuantas
            if(no_overlapchain(c)) then 
                lnpro = logweightchain(c)-energychainLJ(c)   
                do s=1,nseg                       ! loop over segments 
                    t=type_of_monomer(s)   
                    do j=1,nelem(s)               ! loop over elements of segment  
                        k = indexconf(s,c)%elem(j)
                        lnpro = lnpro +lnexppi(k)*vnucl(j,t)              
                    enddo
                enddo     

                pro = exp(lnpro-lnproshift)  
            
                FEconf_local = FEconf_local+(pro/q)*(log(pro/q)-logweightchain(c))
                Econf_local=Econf_local+pro*energychainLJ(c)
                Rgsqr_local = Rgsqr_local+Rgsqr(c)*pro
                Rendsqr_local = Rendsqr_local+Rendsqr(c)*pro
                bond_angle_local= bond_angle_local +bond_angle(:,c)*pro  
                dihedral_angle_local = dihedral_angle_local +dihedral_angle(:,c)*pro
                nucl_spacing_local = nucl_spacing_local +nucl_spacing(:,c)*pro   
                gyr_tensor_local = gyr_tensor_local + gyr_tensor(:,:,c)*pro    
                Asphparam_local = Asphparam_local + Asphparam(c) * pro   
                
                if(write_Palpha) write(un,*)pro/q 
            endif    
        enddo

        Econf_local=Econf_local/q
        Rgsqr_local = Rgsqr_local/q
        Rendsqr_local = Rendsqr_local/q
        bond_angle_local = bond_angle_local/q
        dihedral_angle_local = dihedral_angle_local/q 
        nucl_spacing_local = nucl_spacing_local/q 
        gyr_tensor_local = gyr_tensor_local/q
        Asphparam_local = Asphparam_local/q

        ! communicate local quantities

        if(rank==0) then

            ! normalize
            
            FEconf = FEconf_local
            Econf = Econf_local
            avRgsqr = Rgsqr_local
            avRendsqr = Rendsqr_local
            avbond_angle = bond_angle_local
            avdihedral_angle = dihedral_angle_local 
            avnucl_spacing = nucl_spacing_local 
            avgyr_tensor = gyr_tensor_local
            avAsphparam = Asphparam_local
            
            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                if(nangles>=1) then
                    call MPI_RECV(bond_angle_local, nangles, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif
                if(ndihedrals>=1) then
                    call MPI_RECV(dihedral_angle_local,ndihedrals,MPI_DOUBLE_PRECISION,source,&
                    tag,MPI_COMM_WORLD,stat,ierr)
                endif
                if(nbonds>=1) then
                    call MPI_RECV(nucl_spacing_local,nbonds,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                endif
                call MPI_RECV(gyr_tensor_local,9,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                call MPI_RECV(Asphparam_local,1,MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)

                FEconf=FEconf+FEconf_local
                Econf =Econf +Econf_local
                avRgsqr=avRgsqr+Rgsqr_local
                avRendsqr=avRendsqr+Rendsqr_local
                avbond_angle = avbond_angle+bond_angle_local
                avdihedral_angle =  avdihedral_angle + dihedral_angle_local
                avnucl_spacing = avnucl_spacing + nucl_spacing_local 
                avgyr_tensor = avgyr_tensor + gyr_tensor_local
                avAsphparam = avAsphparam + Asphparam_local
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            if(nangles>=1) then
                call MPI_SEND(bond_angle_local, nangles , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(ndihedrals>=1) then    
                call MPI_SEND(dihedral_angle_local,ndihedrals, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            if(nbonds>=1) then
                call MPI_SEND(nucl_spacing_local,nbonds, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            endif
            call MPI_SEND(gyr_tensor_local,9, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Asphparam_local,1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif

        if(write_Palpha) close(un)

    end subroutine FEconf_nucl_neutral_sv



    subroutine make_filename_Palpha(fname,rank)

        use myutils, only : lenText

        character(len=*), intent(inout) :: fname
        integer, intent(in) :: rank

        character(len=lenText) :: fnamelabel, istr

        write(istr,'(I4)')rank
        call make_filename_label(fnamelabel)
        fname='Palpha.'//trim(fnamelabel)//'rank'//trim(adjustl(istr))//'.dat'

    end subroutine make_filename_Palpha


    ! Subroutine similar to subroutine  make_filename_label in module module myio in myio.f90
    ! Subroutine used by make_filename_Palpa 
    ! for output Palpha. Reason 'extra' subroutine module myio compile after conform_entropy
    ! and extra variable in name to no exstention .dat.

    subroutine make_filename_label(fnamelabel)

        use globals, only : systype, runtype
        use parameters, only : cNaCl,cKCl,cCaCl2,cMgCl2,pHbulk,VdWepsBB,init_denspol,VdWscale,pKd,dielectscale

        character(len=*), intent(inout) :: fnamelabel

        character(len=20) :: rstr
        real(dp) :: denspol


        denspol=init_denspol()

        !     .. make label filename

        select case(systype)

        case("elect","electnopoly","electA")

            write(rstr,'(F5.3)')denspol
            fnamelabel="phi"//trim(adjustl(rstr))
            write(rstr,'(F5.3)')cNaCl
            fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))

            if(cKCl/=0.0_dp) then
                if(cKCl>=0.001_dp) then
                    write(rstr,'(F5.3)')cKCl
                else
                    write(rstr,'(ES9.2E2)')cKCl
                endif
                fnamelabel=trim(fnamelabel)//"cKCl"//trim(adjustl(rstr))
            endif

            if(cCaCl2>=0.001) then
                write(rstr,'(F5.3)')cCaCl2
            elseif(cCaCl2>0.0) then
                write(rstr,'(ES9.2E2)')cCaCl2
            else
                write(rstr,'(F3.1)')cCaCl2
            endif
            fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
            write(rstr,'(F7.3)')pHbulk
            fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))

        case("neutral","neutralnoVdW","nucl_neutral_sv")

            if(denspol>=0.001) then
                write(rstr,'(F5.3)')denspol
            else
                write(rstr,'(ES9.2E2)')denspol
            endif

            fnamelabel="phi"//trim(adjustl(rstr))
            write(rstr,'(F5.3)')VdWscale%val
            fnamelabel=trim(fnamelabel)//"VdWscale"//trim(adjustl(rstr))

        case("brush_mul","brush_mulnoVdW","brushdna","nucl_ionbin","nucl_ionbin_sv",&
            "brushborn","nucl_ionbin_Mg","nucl_ionbin_MgA")
            
            if(denspol>=0.001) then
                write(rstr,'(F5.3)')denspol
            else
                write(rstr,'(ES9.2E2)')denspol
            endif
            
            fnamelabel="phi"//trim(adjustl(rstr))

            if(cNaCl>=0.001_dp) then
                write(rstr,'(F5.3)')cNaCl
            else
                write(rstr,'(ES9.2E2)')cNaCl
            endif
            
            fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))

            if(cKCl/=0.0_dp) then
                if(cKCl>=0.001_dp) then
                    write(rstr,'(F5.3)')cKCl
                else
                    write(rstr,'(ES9.2E2)')cKCl
                endif
                fnamelabel=trim(fnamelabel)//"cKCl"//trim(adjustl(rstr))
            endif

            if(cCaCl2/=0.0_dp) then
                if(cCaCl2>=0.001) then
                    if(cCaCl2>=0.01) then
                        write(rstr,'(F5.3)')cCaCl2
                    else
                        write(rstr,'(F6.4)')cCaCl2
                    endif  
                elseif(cCaCl2>0.0) then
                    write(rstr,'(ES9.2E2)')cCaCl2
                else
                    write(rstr,'(F3.1)')cCaCl2
                endif
                fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
            endif
            
            if(cMgCl2/=0.0_dp) then
                if(cMgCl2>=0.001) then
                    if(cMgCl2>=0.01) then
                        write(rstr,'(F5.3)')cMgCl2
                    else
                        write(rstr,'(F6.4)')cMgCl2
                    endif  
                elseif(cMgCl2>0.0) then
                    write(rstr,'(ES9.2E2)')cMgCl2
                else
                    write(rstr,'(F3.1)')cMgCl2
                endif

                fnamelabel=trim(fnamelabel)//"cMgCl2"//trim(adjustl(rstr))
            endif

            write(rstr,'(F7.3)')pHbulk
            fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))

            if(runtype=="rangepKd") then
                if(pKd%val<0) then
                    write(rstr,'(F5.2)')pKd%val
                else
                    write(rstr,'(F5.3)')pKd%val
                endif
                fnamelabel=trim(fnamelabel)//"pKd"//trim(adjustl(rstr))

            elseif(runtype=="rangeVdWeps") then
                write(rstr,'(F5.3)')VdWscale%val
                fnamelabel=trim(fnamelabel)//"VdWscale"//trim(adjustl(rstr))

            elseif(runtype=="rangedielect") then
                write(rstr,'(ES9.2E2)')dielectscale%val
                fnamelabel=trim(fnamelabel)//"dielectscale"//trim(adjustl(rstr))

            endif


        case default
            print*,"Error in make_filename_label_Palpha subroutine"
            print*,"Wrong value systype : ", systype
        endselect

    end subroutine make_filename_label

end module conform_entropy
