module Poisson  

    use precision_definition

    implicit none

    real(dp), parameter :: epsabsDpsi = 1.0e-8_dp  ! tolerance for absDpsi

    private :: epsabsDpsi
    private :: ipbc 

contains

     subroutine Poisson_Equation(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)

        use volume, only : geometry

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)

        if(geometry=="cubic") then 

            call Poisson_Equation_cubic(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
        
        else if (geometry=="prism") then 
        
            call Poisson_Equation_prism(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
        
        endif
        
    end subroutine
        

    subroutine Poisson_Equation_cubic(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)

        use globals, only : nsize, neq
        use parameters, only : constqW
        use volume, only : nx,ny,nz, linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)
        
        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: idxR, idxL, id2D  
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        noffset=nsize

        do ix=1,nx
            do iy=1,ny
                do iz=2,nz-1
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                    call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)

                    fvec(noffset+id)= -0.5_dp*( psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)
                enddo
            enddo
        enddo    

        ! boundary iz=1 

        do ix=1,nx
            do iy=1,ny
                iz=1
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)

                fvec(noffset+id)= -0.5_dp*( psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls) +sigmaqSurfL(id) &
                    - 5.0_dp*psi(id) +rhoq(id)*constqW)
            
            enddo
        enddo    

        ! boundary iz=nz 

        do ix=1,nx
            do iy=1,ny
                iz=nz
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                
                id2D=id-(nsize-nx*ny)
                
                fvec(noffset+id)= -0.5_dp*( psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzmin) +sigmaqSurfR(id2D) &
                    -5.0_dp*psi(id) +rhoq(id)*constqW)
            enddo
        enddo    


    end subroutine Poisson_Equation_cubic


    subroutine Poisson_Equation_prism(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)

        use globals, only : nsize, neq
        use parameters, only : constqW
        use volume, only : nx,ny,nz, linearIndexFromCoordinate
        use volume, only : cos_two_beta, sin_two_beta
        
        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)
        
        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: idxR, idxL, id2D  
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        integer :: idxypls, idxymin, idxplsymin, idxminypls

        ! .. electrostatics 
       
        noffset=nsize

        do ix=1,nx
            do iy=1,ny
                do iz=2,nz-1
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                    call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy+1,ny),iz ,idxypls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy-1,ny),iz ,idxymin)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy-1,ny),iz ,idxplsymin)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy+1,ny),iz ,idxminypls)

                    
                    fvec(noffset+id)= -0.5_dp*(                                             &
                        (psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id)     &
                        -sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls) +psi(idxymin))/2.0_dp) & 
                        /cos_two_beta+psi(idzpls)-2.0_dp*psi(id)+ psi(idzmin)  +rhoq(id)*constqW)
                enddo
            enddo
        enddo    

        ! boundary iz=1 

        do ix=1,nx
            do iy=1,ny

                iz=1
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy+1,ny),iz ,idxypls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy-1,ny),iz ,idxymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy-1,ny),iz ,idxplsymin)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy+1,ny),iz ,idxminypls)

            
                fvec(noffset+id)= -0.5_dp*(                                                   &
                        (psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id)       &
                        -sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls)+psi(idxymin) )/2.0_dp &
                        )/cos_two_beta+ psi(idzpls)-psi(id)+sigmaqSurfL(id)+rhoq(id)*constqW)
            
            enddo
        enddo    

        ! boundary iz=nz 

        do ix=1,nx
            do iy=1,ny
                iz=nz
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy+1,ny),iz ,idxypls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy-1,ny),iz ,idxymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy-1,ny),iz ,idxplsymin)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy+1,ny),iz ,idxminypls)

                id2D=id-(nsize-nx*ny)
                
                fvec(noffset+id)= -0.5_dp*(                                                 &
                    (psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id)         &
                    -sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls)+psi(idxymin) )/2.0_dp &
                    )/cos_two_beta+ sigmaqSurfR(id2D)-psi(id)+ psi(idzmin)+rhoq(id)*constqW)
            enddo
        enddo    
        
    end subroutine Poisson_Equation_prism


    
    subroutine Poisson_Equation_Surface(fvec,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

        use globals, only : nsize, neq, LEFT, RIGHT, systype
        use parameters, only : constqW
        use volume, only : nx,ny,nz, linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(inout)  :: fvec(:)
        real(dp), intent(in)  :: psi(:)
        real(dp), intent(in)  :: rhoq(:)
        real(dp), intent(inout) :: psisurfR(:)
        real(dp), intent(inout) :: psisurfL(:)
        real(dp), intent(in)  :: sigmaqSurfR(:)
        real(dp), intent(in)  :: sigmaqSurfL(:)
        character(len=2), intent(in)  :: bcflag(2)

        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: idxR, idxL, idxR2D
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        integer :: neq_bc

        ! .. electrostatics: self consistent boundary conditions
        
        if(systype=="elect") then 
            noffset=4*nsize
        else if(systype=="electA") then 
            noffset=3*nsize
        else if(systype=="electdouble") then 
            noffset=4*nsize
        else if(systype=="electnopoly") then 
            noffset=2*nsize
        else if(systype=="dipolarstrong") then 
            noffset=2*nsize
        else if(systype=="dipolarweak") then 
            noffset=2*nsize
        else if(systype=="dipolarnopoly") then 
            noffset=2*nsize
        else 
            print*,"error: systype wrong value for Poisson_equation_surface "    
        endif     

        neq_bc=0

        if(bcflag(RIGHT)/='cc') then
            neq_bc=nx*ny
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,nz,idxR)
                    idxR2D=idxR-(nsize-nx*ny) ! check this 
                    fvec(noffset+idxR2D)=psisurfR(idxR2D)-psi(idxR)-sigmaqSurfR(idxR2D)/2.0_dp
                enddo
            enddo
        else
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,nz,idxR)
                    idxR2D=idxR-(nsize-nx*ny)
                    psisurfR(idxR2D) = psi(idxR)+sigmaqSurfR(idxR2D)/2.0_dp 
                enddo
            enddo
        endif    

        noffset=noffset +neq_bc

        if(bcflag(LEFT)/='cc') then 
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,1,idxL)
                    fvec(noffset+idxL)=psi(idxL)-psisurfL(idxL)+sigmaqSurfL(idxL)/2.0_dp
                enddo
            enddo
        else    
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,1,idxL)
                    psisurfL(idxL) = psi(idxL)+sigmaqSurfL(idxL)/2.0_dp
                enddo
            enddo
        endif   

    end subroutine

        
    subroutine Poisson_Pol_Equation(fvec,D2psi,rhoq,rhob)

        use globals, only : nsize, neq
        use parameters, only : constq0
        use volume, only : nx,ny,nz, linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: D2psi(:) ! nabla of elect potential
        real(dp), intent(in) :: rhoq(:) ! free charge density 
        real(dp), intent(in) :: rhob(:) ! bound charg density 
        
        ! local variables
        integer :: i, noffset

        ! .. electrostatics 
       
        noffset=nsize

        do i=1,nsize    
            fvec(noffset+i)= -0.5_dp*( D2psi(i)+( rhoq(i)+rhob(i) )*constq0 )
        enddo    

    end subroutine


     subroutine grad_and_nabla_pot(psi,Dpsi,D2psi,absDpsi,unitdirDpsi,sigmaqSurfR,sigmaqSurfL)
        
        use volume, only : geometry

        implicit none

        ! input arguments 
        real(dp), intent(in) :: psi(:)
        real(dp), intent(inout) :: Dpsi(:,:),unitdirDpsi(:,:)
        real(dp), intent(inout) :: D2psi(:),absDpsi(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)


        if(geometry=="cubic") then 

            call grad_and_nabla_pot_cubic(psi,Dpsi,D2psi,absDpsi,unitdirDpsi,sigmaqSurfR,sigmaqSurfL)
        
        else if (geometry=="prism") then 
        
            call grad_and_nabla_pot_prism(psi,Dpsi,D2psi,absDpsi,unitdirDpsi,sigmaqSurfR,sigmaqSurfL)
        
        endif

    end subroutine grad_and_nabla_pot

        

    ! computes gradient, nalba (double derivative) , absolute value of gradient and unit direction of gradient of 
    ! electrostaic potential

    subroutine grad_and_nabla_pot_cubic(psi,Dpsi,D2psi,absDpsi,unitdirDpsi,sigmaqSurfR,sigmaqSurfL)
        
        use globals, only : nsize
        use volume, only : nx,ny,nz,delta,linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(in) :: psi(:)
        real(dp), intent(inout) :: Dpsi(:,:),unitdirDpsi(:,:)
        real(dp), intent(inout) :: D2psi(:),absDpsi(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)
        
        ! local variables
        integer :: ix, iy, iz
        integer :: idxR, idxL, id2D  
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        do ix=1,nx
            do iy=1,ny
                do iz=2,nz-1
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                    call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
    
                    ! grad
                    Dpsi(1,id)=(psi(idxpls)-psi(idxmin))/(2.0_dp*delta)
                    Dpsi(2,id)=(psi(idypls)-psi(idymin))/(2.0_dp*delta)
                    Dpsi(3,id)=(psi(idzpls)-psi(idzmin))/(2.0_dp*delta)
    
                    ! nabla
                    D2psi(id)=psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)+psi(idzpls)+&
                        psi(idzmin)-6.0_dp*psi(id) 
                enddo
            enddo
        enddo    

        ! boundary iz=1 

        do ix=1,nx
            do iy=1,ny
                iz=1
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)


                Dpsi(1,id)=(psi(idxpls)-psi(idxmin))/(2.0_dp*delta)
                Dpsi(2,id)=(psi(idypls)-psi(idymin))/(2.0_dp*delta)
                Dpsi(3,id)=(psi(idzpls)-(psi(id)+sigmaqSurfL(id)))/(2.0_dp*delta)
                
                D2psi(id)=psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)+psi(idzpls) +sigmaqSurfL(id) &
                    - 5.0_dp*psi(id)
            enddo
        enddo    

        ! boundary iz=nz 

        do ix=1,nx
            do iy=1,ny
                iz=nz
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                
                id2D=id-(nsize-nx*ny)
                
                Dpsi(1,id)=(psi(idxpls)-psi(idxmin))/(2.0_dp*delta)
                Dpsi(2,id)=(psi(idypls)-psi(idymin))/(2.0_dp*delta)
                Dpsi(3,id)=((sigmaqSurfR(id2D)+psi(id))-psi(idzmin))/(2.0_dp*delta)
            
                D2psi(id)= psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzmin) +sigmaqSurfR(id2D) &
                    -5.0_dp*psi(id) 
            enddo
        enddo    

        do id=1,nsize
           
            ! absolute value grad
            absDpsi(id)=sqrt(Dpsi(1,id)**2+Dpsi(2,id)**2+Dpsi(3,id)**2) 

            if(absDpsi(id)>epsabsDpsi) then   
                ! unit vector grad
                unitdirDpsi(1,id)=Dpsi(1,id)/absDpsi(id)
                unitdirDpsi(2,id)=Dpsi(2,id)/absDpsi(id)
                unitdirDpsi(3,id)=Dpsi(3,id)/absDpsi(id)
            else
                ! lenght should be unimportat if unit is small 
                unitdirDpsi(1,id)=1.0/3.0_dp
                unitdirDpsi(2,id)=1.0/3.0_dp
                unitdirDpsi(3,id)=1.0/3.0_dp
            endif    

        enddo    

    end subroutine grad_and_nabla_pot_cubic
         


    ! computes gradient, nalba ( double derivative), absolute value of gradient and unit direction of gradient of 
    ! electrostaic potential

    subroutine grad_and_nabla_pot_prism(psi,Dpsi,D2psi,absDpsi,unitdirDpsi,sigmaqSurfR,sigmaqSurfL)
        
        use globals, only : nsize
        use volume, only : nx,ny,nz,delta,linearIndexFromCoordinate
        use volume, only : cos_two_beta, sin_two_beta
        
        implicit none

        ! input arguments 
        real(dp), intent(in) :: psi(:)
        real(dp), intent(inout) :: Dpsi(:,:),unitdirDpsi(:,:)
        real(dp), intent(inout) :: D2psi(:),absDpsi(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)
        
        ! local variables
        integer :: ix, iy, iz
        integer :: idxR, idxL, id2D 
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin 
        integer :: idxypls, idxymin, idxplsymin, idxminypls
        real(dp) :: Dpsiu,Dpsiv

        ! .. electrostatics  

        do ix=1,nx
            do iy=1,ny
                do iz=2,nz-1
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                    call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
    
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy+1,ny),iz ,idxypls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy-1,ny),iz ,idxymin)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy-1,ny),iz ,idxplsymin)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy+1,ny),iz ,idxminypls)

                    ! grad

                    Dpsiu=(psi(idxpls)-psi(idxmin))/(2.0_dp*delta)
                    Dpsiv=(psi(idypls)-psi(idymin))/(2.0_dp*delta)

                    Dpsi(1,id)=Dpsiu-sin_two_beta*Dpsiv
                    Dpsi(2,id)=Dpsiv-sin_two_beta*Dpsiu
                    Dpsi(3,id)=(psi(idzpls)-psi(idzmin))/(2.0_dp*delta)
    
                    ! nabla
                    !D2psi(id)=( psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id) + &
                    !   2.0_dp*sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls) +     &
                    !    psi(idxymin) ))/cos_two_beta+psi(idzpls)-2.0_dp*psi(id)+ psi(idzmin)
                    
                    D2psi(id)=(psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id) &
                        -sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls) + psi(idxymin))/2.0_dp &
                        )/cos_two_beta+psi(idzpls)-2.0_dp*psi(id)+ psi(idzmin)

                enddo
            enddo
        enddo    

        ! boundary iz=1 

        do ix=1,nx
            do iy=1,ny
                iz=1
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy+1,ny),iz ,idxypls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy-1,ny),iz ,idxymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy-1,ny),iz ,idxplsymin)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy+1,ny),iz ,idxminypls)

                Dpsiu=(psi(idxpls)-psi(idxmin))/(2.0_dp*delta)
                Dpsiv=(psi(idypls)-psi(idymin))/(2.0_dp*delta)

                Dpsi(1,id)=Dpsiu-sin_two_beta*Dpsiv
                Dpsi(2,id)=Dpsiv-sin_two_beta*Dpsiu
                Dpsi(3,id)=(psi(idzpls)-(psi(id)+sigmaqSurfL(id)))/(2.0_dp*delta)

                ! nabla
                !D2psi(id)=( psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id) + &
                !        2.0_dp*sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls) + &
                !        psi(idxymin) ))/cos_two_beta+psi(idzpls)-psi(id)+sigmaqSurfL(id)
                D2psi(id)=( psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id) &
                    -sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls)+psi(idxymin))/2.0_dp &
                    )/cos_two_beta+psi(idzpls)-psi(id)+sigmaqSurfL(id)
            enddo
        enddo    

        ! bound`ary iz=nz 

        do ix=1,nx
            do iy=1,ny
                iz=nz
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy+1,ny),iz ,idxypls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy-1,ny),iz ,idxymin)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),ipbc(iy-1,ny),iz ,idxplsymin)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),ipbc(iy+1,ny),iz ,idxminypls)

                id2D=id-(nsize-nx*ny)
                Dpsiu=(psi(idxpls)-psi(idxmin))/(2.0_dp*delta)
                Dpsiv=(psi(idypls)-psi(idymin))/(2.0_dp*delta)

                Dpsi(1,id)=Dpsiu-sin_two_beta*Dpsiv
                Dpsi(2,id)=Dpsiv-sin_two_beta*Dpsiu
                Dpsi(3,id)=((sigmaqSurfR(id2D)+psi(id))-psi(idzmin))/(2.0_dp*delta)
            
                ! nabla
                D2psi(id)=( psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id)            &
                    -sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls)+psi(idxymin))/2.0_dp &
                    )/cos_two_beta+sigmaqSurfR(id2D)-psi(id)+ psi(idzmin)
            enddo
        enddo    

        do id=1,nsize
            ! absolute value gradient
            absDpsi(id)=sqrt(Dpsi(1,id)**2+Dpsi(2,id)**2+Dpsi(3,id)**2) 
            
            ! this must be 


            if(absDpsi(id)>epsabsDpsi) then   
                ! unit vector grad
                unitdirDpsi(1,id)=Dpsi(1,id)/absDpsi(id)
                unitdirDpsi(2,id)=Dpsi(2,id)/absDpsi(id)
                unitdirDpsi(3,id)=Dpsi(3,id)/absDpsi(id)
            else
                ! lenght should be unimportant if unit is very small 
                unitdirDpsi(1,id)=1.0/3.0_dp
                unitdirDpsi(2,id)=1.0/3.0_dp
                unitdirDpsi(3,id)=1.0/3.0_dp
            endif    
        enddo    

    end subroutine grad_and_nabla_pot_prism


    ! .. This routine computes the bound charge density
    ! .. rhob= -div.P is the divergence of the polarization density vector


    subroutine charge_density_bound(electPol,rhob)
 
        use volume, only : nx,ny,nz,delta,linearIndexFromCoordinate
        use volume, only : cos_two_beta

        implicit none

        ! input arguments 
        real(dp), intent(in)  :: electPol(:,:)
        real(dp), intent(inout)  :: rhob(:)
        
        ! local variables
        integer :: ix, iy, iz
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        real(dp) :: sqrtcos

        sqrtcos=sqrt(cos_two_beta)

        do ix=1,nx
            do iy=1,ny
                do iz=2,nz-1
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                    call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)

                    rhob(id)=-(sqrtcos*((electPol(1,idxpls)-electPol(1,idxmin))+ &
                        (electPol(2,idypls)-electPol(2,idymin)))+     &
                        (electPol(3,idzpls)-electPol(3,idzmin)))/(2.0_dp*delta)    
                enddo
            enddo
        enddo    

        ! boundary iz=1 

        do ix=1,nx
            do iy=1,ny
                iz=1
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)

                rhob(id)=-(sqrtcos*((electPol(1,idxpls)-electPol(1,idxmin)) + &
                    (electPol(2,idypls)-electPol(2,idymin))) + &
                    (electPol(3,idzpls)+electPol(3,id    )) )/(2.0_dp*delta)    
            enddo
        enddo    

        ! boundary iz=nz 

        do ix=1,nx
            do iy=1,ny
                iz=nz
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                
                rhob(id)=-(sqrtcos*((electPol(1,idxpls)-electPol(1,idxmin)) + &
                    (electPol(2,idypls)-electPol(2,idymin))) + &
                    ( -electPol(3,id)  -electPol(3,idzmin)) )/(2.0_dp*delta)
            enddo
        enddo    

    end subroutine charge_density_bound

            

    function ipbc(ival,imax) result(intpbc)
        implicit none 
        integer, intent(in) :: ival
        integer, intent(in) :: imax
        integer :: intpbc

        if(ival>0) then
            intpbc=ival-int((ival-1)/imax)*imax
        else
            intpbc=ival-(int((ival-1)/imax)-1)*imax
        endif

    end function


end module
