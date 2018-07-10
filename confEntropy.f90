! --------------------------------------------------------------|
! ConfEntropy.f90:                                              |
! constructs the free energy Conformational Entropy             |
!  beta Fconf = sigma \sum_alpha P(\alpha)lnP{\alpha)           |
! --------------------------------------------------------------|


module conform_entropy

    use mpivars
    implicit none

    private                    ! default all routines in this module private 
    public  ::  FEconf_entropy ! only this subroutine  public
   

contains

    subroutine FEconf_entropy(FEconfAB,FEconfC)

 
        use globals
        use myutils
        
        implicit none
    
        real(dp), intent(out) :: FEconfAB, FEconfC  

        character(len=lenText) :: text 

        select case (sysflag) 
        case ("elect")
            call FEconf_elect(FEconfAB)
            FEconfC=0.0_dp
        case ("electdouble")
            call FEconf_electdouble(FEconfAB)
            FEconfC=0.0_dp
        case ("neutral")
            call FEconf_neutral(FEconfAB,FeconfC)
        case ("electHC")     
            text="fcenergy: sysflag: "//sysflag//" not implemented yet"
            call print_to_log(LogUnit,text)
            print*,text
        case default
            text="FEconf_entropy: wrong sysflag: "//sysflag//"stopping program"
            call print_to_log(LogUnit,text)
            stop
        end select   

    end subroutine FEconf_entropy

  ! computes conformational entropy in neutral state 

    subroutine FEconf_neutral(FEconfAB,FeconfC)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW

        implicit none
        
        real(dp), intent(out) :: FEconfAB, FEconfC  

        !     .. declare local variables

        ! real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 

        ! integer :: i,j,k,c,s         ! dummy indices
        ! real(dp) :: pro,tmp,expVdW 
        ! integer :: conf              ! counts number of conformations
        

        ! real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 


        ! !     .. executable statements 

        ! do i=1,nz    
        !     exppiA(i)=(xsol(i)**vpolA(3)) !*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
        !     exppiB(i)=(xsol(i)**vpolB(3)) !*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   
        !     exppiC(i)=(xsol(i)**vpolC)
       
        !     !     .. VdW interaction   
        !     tmp = 0.0_dp
        !     if((i+VdWcutoffdelta)<=nsize) then 
        !         do j=minrange(i),i+VdWcutoffdelta
        !             tmp = tmp + chis(i,j)*rhopolB(j)*vpolB(3)*vsol
        !         enddo
        !     endif
        !     expVdW=dexp(-VdWepsB*tmp)
        !     exppiB(i)=exppiB(i)*expVdW ! auxiliary variable
        ! enddo

        
        ! FEconfAB=0.0_dp

        ! do c=1,cuantasAB            ! loop over cuantas
        !     pro=1.0_dp                ! initial weight conformation 
        !     do s=1,nsegAB            ! loop over segments 
        !         k=indexchainAB(c,s)
        !         if(isAmonomer(s)) then ! A segment 
        !             pro = pro*exppiA(k)
        !         else
        !             pro = pro*exppiB(k)
        !         endif
        !     enddo
        !     FEconfAB=FEconfAB+pro*log(pro)
        ! enddo

        ! ! normalize needs to changed 
        ! !FEconfAB=(FEconfAB/qAB-log(qAB))*(sigmaAB*delta)    

        ! FEconfC = 0.0_dp                   
        ! do c=1,cuantasC            ! loop over cuantas                                                      
        !     pro=1.0_dp               ! initial weight conformation                                                   
        !     do s=1,nsegC            ! loop over segments                
        !         k=indexchainC(c,s)
        !         pro = pro*exppiC(k)
        !     enddo
        !     FEconfC = FEconfC+pro*log(pro)
        ! enddo
        ! ! normalize
        !FEconfC=(FEConfC/qC -log(qC))*(sigmaC*delta)    

    end subroutine FEconf_neutral



    subroutine FEconf_elect(FEconfAB)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none

        real(dp), intent(out) :: FEconfAB
        
        !     .. declare local variables

        ! real(dp) :: exppiA(nsize),exppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        ! integer :: i,k,c,s         ! dummy indices
        ! real(dp) :: pro


        ! !     .. executable statements 

       
        ! do i=1,nz
        !       exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
        !       exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable
        ! enddo
        

        ! FEconfAB=0.0_dp

        ! do c=1,cuantasAB            ! loop over cuantas
        !     pro=1.0_dp                ! initial weight conformation 
        !     do s=1,nsegAB            ! loop over segments 
        !         k=indexchainAB(c,s)
        !         if(isAmonomer(s)) then ! A segment 
        !             pro = pro*exppiA(k)
        !         else
        !             pro = pro*exppiB(k)
        !         endif
        !     enddo
        !     FEconfAB=FEconfAB+pro*log(pro)
        ! enddo

        ! ! normalize

        ! !FEconfAB=(FEconfAB/qAB-log(qAB))*(sigmaAB*delta)     ! check sigma

        
    end subroutine FEconf_elect


    subroutine FEconf_electdouble(FEconfAB)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none

        real(dp), intent(out) :: FEconfAB
        
        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer  :: i,k,c,s, kL,kR, g,gn         ! dummy indices
        real(dp) :: proL,proR
        real(dp) :: FEconfABL,FEconfABR
        real(dp) :: FEconfABL_local(ngr_node),FEconfABR_local(ngr_node)
        real(dp) :: FEconfABL_array(size*ngr_node),FEconfABR_array(size*ngr_node)
        real(dp) :: qABL_local(ngr_node),qABR_local(ngr_node)
        real(dp) :: qABL_array(size*ngr_node),qABR_array(size*ngr_node)

        ! .. executable statements 

        ! .. communicate xsol,psi and fdsiA(1:) and fdisB(1,:) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(fdisA(1,:),nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(fdisB(1,:),nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(fdisA(1,:), nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)    
            call MPI_RECV(fdisB(1,:), nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
        endif
            
        do i=1,nz
              exppiA(i)=(xsol(i)**vpolA(1))*exp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
              exppiB(i)=(xsol(i)**vpolB(1))*exp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable
        enddo
       
        FEconfABL=0.0_dp
        FEconfABR=0.0_dp

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            FEconfABL_local(gn)= 0.0_dp
            FEconfABR_local(gn)= 0.0_dp

            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                proL=1.0_dp                ! initial weight conformation 
                proR=1.0_dp
            
                if(weightchainAB(gn,c)) then ! initial weight conformation 

                    proL=1.0_dp
                    proR=1.0_dp

                    do s=1,nsegAB              ! loop over segments 
                        kL=indexchainAB(s,gn,c)      
                        kR=mirror_index(kL,nz)    
                        if(isAmonomer(s)) then ! A segment 
                            proL = proL*exppiA(kL)
                            proR = proR*exppiA(kR)
                        else
                            proL = proL*exppiB(kL)
                            proR = proR*exppiB(kR)
                        endif
                    enddo

                   FEconfABL_local(gn)=FEconfABL_local(gn)+proL*log(proL)
                   FEconfABR_local(gn)=FEconfABR_local(gn)+proR*log(proR)
                
                endif
            
            enddo

        enddo   

        ! communicate FEconfABL and FEcopnd_ABR

        if(rank==0) then

            do gn=1,ngr_node
                g = (0)*ngr_node+gn
                FEconfABL_array(g)=FEconfABL_local(gn)
                FEconfABR_array(g)=FEconfABR_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(FEconfABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(FEconfABR_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(qABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(qABR_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    FEconfABL_array(g)=FEconfABL_local(gn)
                    FEconfABR_array(g)=FEconfABR_local(gn)
                    qABL_array(g)=qABL_local(gn)
                    qABR_array(g)=qABR_local(gn)
                enddo
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconfABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(FEconfABR_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABR_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            do g=1,ngr_node*size-1
                FEconfABL = FEconfABL+ (FEconfABL_array(g)/qABL_array(g)-log(qABL_array(g)))     ! check sigma
                FEconfABR = FEconfABR+ (FEconfABR_array(g)/qABR_array(g)-log(qABR_array(g)))    
            enddo    

            FEconfAB=FEconfABL+FEconfABR
        endif


    end subroutine FEconf_electdouble
    

end module conform_entropy