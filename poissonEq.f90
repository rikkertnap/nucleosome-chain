module Poisson  

    use precision_definition

    implicit none

    real(dp), parameter :: epsabsDpsi = 1.0e-8_dp  ! tolerance for absDpsi
    integer :: noffset_bc 

    private :: epsabsDpsi
    private :: ipbc
    private :: noffset_bc 

contains

    subroutine Poisson_Equation(fvec,psi,rhoq)

        use volume, only : geometry

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
    
        if(geometry=="cubic") then 

            call Poisson_Equation_cubic(fvec,psi,rhoq)
        
        else if (geometry=="prism") then 
        
            call Poisson_Equation_prism(fvec,psi,rhoq)

        else 

            print*,"Wrong geometry in Poisson_equation"   
        
        endif
        
    end subroutine Poisson_equation

    subroutine Poisson_Equation_bc(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)

        use volume, only : geometry

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)
    
        if(geometry=="cubic") then 

            call Poisson_Equation_cubic_bc(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
           
        else 

            print*,"Wrong geometry in Poisson_equation_bc"   
        
        endif
        
    end subroutine Poisson_equation_bc


    subroutine Poisson_Equation_ST(fvec,psi,rhoq)

        use volume, only : geometry
        use parameters, only : psiSL, psiSR

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
    
        if(geometry=="cubic") then 

            call Poisson_Equation_cubic_ST(fvec,psi,rhoq,psiSL,psiSR)
    
        else 

            print*,"Wrong geometry in Poisson_equation"   
        
        endif
        
    end subroutine Poisson_equation_ST

    subroutine Poisson_Equation_Eps(fvec,psi,rhoq,eps)

        use volume, only : geometry

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: eps(:)
    
        if(geometry=="cubic") then 

            call Poisson_Equation_Eps_cubic(fvec,psi,rhoq,eps)
        
        else if (geometry=="prism") then 
        
              print*,"Poisson Eq with varying dielectric constant for prism coordinates not implemented yet."
    
        endif
        
    end subroutine

    subroutine Poisson_Equation_cubic(fvec,psi,rhoq)

        use globals, only : nsize
        use parameters, only : constqW
        use volume, only : nx,ny,nz, coordtoindex ! linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        
        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        noffset = nsize

        do ix=1,nx
            do iy=1,ny
                do iz=1,nz      

                    id     = coordtoindex(ix,           iy,iz)
                    idxpls = coordtoindex(ipbc(ix+1,nx),iy,iz)
                    idxmin = coordtoindex(ipbc(ix-1,nx),iy,iz)
                    idzpls = coordtoindex(ix,           iy,ipbc(iz+1,nz))
                    idzmin = coordtoindex(ix,           iy,ipbc(iz-1,nz))
                    idypls = coordtoindex(ix,ipbc(iy+1,ny),iz)
                    idymin = coordtoindex(ix,ipbc(iy-1,ny),iz)

                    fvec(noffset+id)= -0.5_dp*(&
                        psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)
                enddo
            enddo
        enddo    

    end subroutine Poisson_Equation_cubic


    subroutine Poisson_Equation_prism(fvec,psi,rhoq)

        use globals, only : nsize
        use parameters, only : constqW
        use volume, only : nx,ny,nz, coordtoindex
        use volume, only : cos_two_beta, sin_two_beta
        
        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        
        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        integer :: idxypls, idxymin, idxplsymin, idxminypls

        ! .. electrostatics 
       
        noffset=nsize

        do ix=1,nx
            do iy=1,ny
                do iz=1,nz
                        
                    id         = coordtoindex(ix,           iy,iz  )
                    idxpls     = coordtoindex(ipbc(ix+1,nx),iy,iz  )
                    idxmin     = coordtoindex(ipbc(ix-1,nx),iy,iz  )
                    idzpls     = coordtoindex(ix,           iy,ipbc(iz+1,nz))
                    idzmin     = coordtoindex(ix,           iy,ipbc(iz-1,nz))
                    idypls     = coordtoindex(ix,ipbc(iy+1,ny),iz  )
                    idymin     = coordtoindex(ix,ipbc(iy-1,ny),iz  )
                    idxypls    = coordtoindex(ipbc(ix+1,nx),ipbc(iy+1,ny),iz )
                    idxymin    = coordtoindex(ipbc(ix-1,nx),ipbc(iy-1,ny),iz )
                    idxplsymin = coordtoindex(ipbc(ix+1,nx),ipbc(iy-1,ny),iz )
                    idxminypls = coordtoindex(ipbc(ix-1,nx),ipbc(iy+1,ny),iz )
                    
                    fvec(noffset+id)= -0.5_dp*(                                             &
                        (psi(idxpls)+psi(idxmin)+psi(idypls)+psi(idymin)-4.0_dp*psi(id)     &
                        -sin_two_beta*(psi(idxypls)-psi(idxplsymin)-psi(idxminypls) +psi(idxymin))/2.0_dp) & 
                        /cos_two_beta+psi(idzpls)-2.0_dp*psi(id)+ psi(idzmin)  +rhoq(id)*constqW)
                enddo
            enddo
        enddo    
 
        
    end subroutine Poisson_Equation_prism


    subroutine Poisson_Equation_Eps_cubic(fvec,psi,rhoq,eps)

        use globals, only : nsize
        use parameters, only : constqW
        use volume, only : nx,ny,nz, coordtoindex  

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: eps(:)
        
        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        noffset=nsize

        do ix=1,nx
            do iy=1,ny
                do iz=1,nz
               
                    id     = coordtoindex(ix,           iy,iz)
                    idxpls = coordtoindex(ipbc(ix+1,nx),iy,iz)
                    idxmin = coordtoindex(ipbc(ix-1,nx),iy,iz)
                    idzpls = coordtoindex(ix,           iy,ipbc(iz+1,nz))
                    idzmin = coordtoindex(ix,           iy,ipbc(iz-1,nz))
                    idypls = coordtoindex(ix,ipbc(iy+1,ny),iz)
                    idymin = coordtoindex(ix,ipbc(iy-1,ny),iz)


                    fvec(noffset+id)= -0.5_dp*( &
                         (eps(idxpls)+eps(id)    )*psi(idxpls) +& 
                         (eps(id)    +eps(idxmin))*psi(idxmin) +&
                         (eps(idypls)+eps(id)    )*psi(idypls) +&
                         (eps(id)    +eps(idymin))*psi(idymin) +& 
                         (eps(idzpls)+eps(id)    )*psi(idzpls) +&
                         (eps(id)    +eps(idzmin))*psi(idzmin) &
                        -(eps(idxpls)+eps(idxmin)+eps(idypls)+eps(idymin) + &
                          eps(idzpls)+eps(idzmin)+6.0_dp*eps(id)) *psi(id) +2.0_dp*rhoq(id)*constqW)
                enddo
            enddo
        enddo    

    end subroutine Poisson_Equation_Eps_cubic


    subroutine Poisson_Equation_cubic_bc(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)

        use globals, only : nsize
        use parameters, only : constqW
        use volume, only : nx, ny, nz, coordtoindex ! linearIndexFromCoordinate

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)
        
        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: id2D  
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        noffset = nsize

        do ix=1,nx
            do iy=1,ny
                do iz=2,nz-1
                  
                    id      = coordtoindex(ix,           iy,iz)
                    idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)
                    idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                    idzpls  = coordtoindex(ix,           iy,iz+1)
                    idzmin  = coordtoindex(ix,           iy,iz-1)
                    idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                    idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)

                    fvec(noffset+id)= -0.5_dp*(&
                         psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)
                enddo
            enddo
        enddo    

        ! boundary iz=1 

        do ix=1,nx
            do iy=1,ny
                iz=1
                id      = coordtoindex(ix,           iy,iz)
                idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)
                idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                idzpls  = coordtoindex(ix,           iy,iz+1)
                idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)

                fvec(noffset+id)= -0.5_dp*(&
                     psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls) +sigmaqSurfL(id) &
                    - 5.0_dp*psi(id) +rhoq(id)*constqW)
            
            enddo
        enddo    

        ! boundary iz=nz 

        do ix=1,nx
            do iy=1,ny
                iz=nz
                id      = coordtoindex(ix,           iy,iz)
                idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)
                idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                idzmin  = coordtoindex(ix,           iy,iz-1)
                idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)
                
                id2D = id - (nsize-nx*ny)
                
                fvec(noffset+id)= -0.5_dp*(&
                     psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzmin) +sigmaqSurfR(id2D) &
                    -5.0_dp*psi(id) +rhoq(id)*constqW)
            enddo
        enddo    


    end subroutine Poisson_Equation_cubic_bc


    subroutine Poisson_Equation_cubic_ST_pbc(fvec,psi,rhoq,psizmin,psizmax)

        use globals, only : nsize
        use parameters, only : constqW
        use volume, only : nx, ny, nz, coordtoindex ! linearIndexFromCoordinate

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: psizmin
        real(dp), intent(in) :: psizmax
        
        ! local variables
        integer :: ix, iy, iz, nshift
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        nshift = nsize

        ! inside 
        do iz=2,nz-1
            do iy=1,ny
                do ix=1,nx
                  
                    id      = coordtoindex(ix,           iy,iz)
                    idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)
                    idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                    idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                    idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)
                    idzpls  = coordtoindex(ix,           iy,iz+1)
                    idzmin  = coordtoindex(ix,           iy,iz-1)

                    fvec(nshift+id)= -0.5_dp*(&
                         psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)
                enddo
            enddo
        enddo    

        ! faces 
        
        ! bottom z-face
        ! boundary z=0 iz=1 

        iz=1
        do iy=1,ny
            do ix=1,nx
                
                id      = coordtoindex(ix,           iy,iz)
                idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)
                idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)
                idzpls  = coordtoindex(ix,           iy,iz+1)
                !idzmin  = coordtoindex(ix,           iy,iz-1)
                
                fvec(nshift+id)= -0.5_dp*(&
                    psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psizmin &
                    -6.0_dp*psi(id) +rhoq(id)*constqW)
            
            enddo
        enddo    

        ! top z-face
        ! boundary z=nz*delta iz=nz

        iz=nz
        do iy=1,ny
            do ix=1,nx

                id      = coordtoindex(ix,           iy,iz)
                idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)
                idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)
                !idzpls  = coordtoindex(ix,           iy,iz+1)
                idzmin  = coordtoindex(ix,           iy,iz-1) 

                fvec(nshift+id)= -0.5_dp*(&
                         psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psizmax+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)

            enddo
        enddo    


    end subroutine Poisson_Equation_cubic_ST_pbc

    subroutine Poisson_Equation_cubic_ST(fvec,psi,rhoq,psizmin,psizmax)

        use globals, only : nsize
        use parameters, only : constqW
        use volume, only : nx, ny, nz, coordtoindex ! linearIndexFromCoordinate

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: psizmin
        real(dp), intent(in) :: psizmax
        
        ! local variables
        integer :: ix, iy, iz, nshift
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        real(dp) :: factor6, factor5, factor4, factor3

        ! .. electrostatics 
       
        nshift = nsize

        factor6 = -1.0_dp/6.0_dp
        factor5 = -1.0_dp/5.0_dp 
        factor4 = -1.0_dp/4.0_dp
        factor3 = -1.0_dp/3.0_dp

        ! inside 

        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                  
                    id      = coordtoindex(ix  ,iy  ,iz  )
                    idxpls  = coordtoindex(ix+1,iy  ,iz  )
                    idxmin  = coordtoindex(ix-1,iy  ,iz  )
                    idypls  = coordtoindex(ix  ,iy+1,iz  )
                    idymin  = coordtoindex(ix  ,iy-1,iz  )
                    idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                    idzmin  = coordtoindex(ix  ,iy  ,iz-1)

                    fvec(nshift+id)= factor6*(&
                         psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)
                   
                enddo
            enddo
        enddo    

        ! faces 
        
        ! bottom z-face
        ! boundary z=0 iz=1 

        iz=1
        do iy=2,ny-1
            do ix=2,nx-1
                
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
               ! idzmin  = coordtoindex(ix  ,iy  ,iz-1)
                
                fvec(nshift+id)= factor6*(&
                    psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psizmin &
                    -6.0_dp*psi(id) +rhoq(id)*constqW)
            
            enddo
        enddo    

        ! top z-face
        ! boundary z=nz*delta iz=nz

        iz=nz
        do iy=2,ny-1
            do ix=2,nx-1

                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )
                !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1)

                fvec(nshift+id)= factor6*(&
                         psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psizmax+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)

            enddo
        enddo    

        ! fvec(nshift+id)= -0.5_dp*(&
        !            psi(idxpls) + psi(idypls) + psi(idzpls) &
        !                   -6.0_dp*psi(id)  + &
        !           psi(idxmin) + psi(idymin) + psi(idzmin)  +rhoq(id)*constqW)
                

        ! rear x-face 
        ! boundary x= 0 plane ix=1  
        
        ix = 1
        do iz=2,nz-1
            do iy=2,ny-1
                
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )
                !idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  ) 
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1)
               

                fvec(nshift+id)= factor6*(&
                    psi(idxpls) + psi(idypls) + psi(idzpls) -5.0_dp*psi(id)  + &
                    psi(idymin) + psi(idzmin) +rhoq(id)*constqW)
                
            enddo
        enddo    


        ! front x-face
        ! boundary x= nx delta  plane ix=nx  
        
        ix = nx
        do iz=2,nz-1
            do iy=2,ny-1
                
                id      = coordtoindex(ix  ,iy  ,iz  )
                !idxpls  = coordtoindex(ix+1,iy  ,iz  )
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  ) 
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1)
               
                fvec(nshift+id)= factor5*(&
                    psi(idypls) + psi(idzpls) -5.0_dp*psi(id)  + &
                    psi(idxmin) + psi(idymin) + psi(idzmin)  +rhoq(id)*constqW)
                 
            enddo
        enddo    



        ! left y-face
        ! boundary y= 0 plane iy=1  

        iy = 1
        do iz=2,nz-1
            do ix=2,nx-1
                
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                !idymin  = coordtoindex(ix  ,iy-1,iz  ) 
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1)
               
                fvec(nshift+id)= factor5*(&
                    psi(idxpls) + psi(idypls) + psi(idzpls) -5.0_dp*psi(id)  + &
                    psi(idxmin) + psi(idzmin) + rhoq(id)*constqW)+rhoq(id)
                 
            enddo
        enddo    

        ! right y-face
        ! boundary y= ny delta  plane iy=ny  

        iy = ny
        do iz=2,nz-1
            do ix=2,nx-1
                
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                !idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  ) 
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1)
               
                fvec(nshift+id)= factor5*(&
                    psi(idxpls) + psi(idzpls) -5.0_dp*psi(id)  + &
                    psi(idxmin) + psi(idymin) + psi(idzmin) + rhoq(id)*constqW)
                 
            enddo
        enddo    
        

        ! corners 

        ! ix=1  iy=1  iz=1 
        ix = 1
        iy = 1
        iz = 1

        id      = coordtoindex(ix  ,iy  ,iz  )
        idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        !idxmin  = coordtoindex(ix-1,iy  ,iz  )
        idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        !idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        idypls  = coordtoindex(ix  ,iy+1,iz  )
        !idymin  = coordtoindex(ix  ,iy-1,iz  )

        fvec(nshift+id)= factor3*(&
            psi(idxpls) + psi(idypls) + psi(idzpls) -3.0_dp*psi(id) +rhoq(id)*constqW)


        ! ix=nx iy=1  iz=1 :
        ix = nx
        iy = 1
        iz = 1


        id      = coordtoindex(ix  ,iy  ,iz  )
        !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        idxmin  = coordtoindex(ix-1,iy  ,iz  )
        idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        !idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        idypls  = coordtoindex(ix  ,iy+1,iz  )
        !idymin  = coordtoindex(ix  ,iy-1,iz  )

        fvec(nshift+id)= factor3*(&
            psi(idypls) + psi(idzpls) -3.0_dp*psi(id) +psi(idxmin) +rhoq(id)*constqW)


        ! ix=nx iy=ny iz=1 : 
        ix = nx
        iy = ny
        iz = 1

        id      = coordtoindex(ix  ,iy  ,iz  )
        !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        idxmin  = coordtoindex(ix-1,iy  ,iz  )
        idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        !idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        !idypls  = coordtoindex(ix  ,iy+1,iz  )
        idymin  = coordtoindex(ix  ,iy-1,iz  )

        fvec(nshift+id)= factor3*(&
            psi(idzpls) -3.0_dp*psi(id) +psi(idxmin) + psi(idymin)  +rhoq(id)*constqW)



         ! ix=1  iy=ny iz=1 :
        ix = 1
        iy = ny
        iz = 1
        
        id      = coordtoindex(ix  ,iy  ,iz  )
        idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        !idxmin  = coordtoindex(ix-1,iy  ,iz  )
        idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        !idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        !idypls  = coordtoindex(ix  ,iy+1,iz  )
        idymin  = coordtoindex(ix  ,iy-1,iz  )

        fvec(nshift+id)= factor3*(&
            psi(idxpls) + psi(idzpls) -3.0_dp*psi(id)  + psi(idymin) +rhoq(id)*constqW)

        ! ix=1  iy=1  iz=nz :
        ix = 1
        iy = 1
        iz = nz

        id      = coordtoindex(ix  ,iy  ,iz  )
        idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        !idxmin  = coordtoindex(ix-1,iy  ,iz  )
        !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        idypls  = coordtoindex(ix  ,iy+1,iz  )
        !idymin  = coordtoindex(ix  ,iy-1,iz  )

        fvec(nshift+id)= factor3*(&
            psi(idxpls) + psi(idypls) -3.0_dp*psi(id) + psi(idzmin)  +rhoq(id)*constqW)



          ! ix=nx iy=1  iz=nz :
        ix = nx
        iy = 1
        iz = nz 

        id      = coordtoindex(ix  ,iy  ,iz  )
        !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        idxmin  = coordtoindex(ix-1,iy  ,iz  )
        !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        idypls  = coordtoindex(ix  ,iy+1,iz  )
        !idymin  = coordtoindex(ix  ,iy-1,iz  )

        fvec(nshift+id)= factor3*(&
            psi(idypls) -3.0_dp*psi(id) + psi(idxmin) + psi(idzmin)  +rhoq(id)*constqW)

        ! ix=nx iy=ny iz=nz :
        ix = nx
        iy = ny
        iz = nz

        id      = coordtoindex(ix  ,iy  ,iz  )
        !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        idxmin  = coordtoindex(ix-1,iy  ,iz  )
        !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        !idypls  = coordtoindex(ix  ,iy+1,iz  )
        idymin  = coordtoindex(ix  ,iy-1,iz  )

        fvec(nshift+id)= factor3*(&
            psi(idypls) -3.0_dp*psi(id) + psi(idxmin) + psi(idzmin)  +rhoq(id)*constqW)

        ! edges

        ! (1,1,1)  -> (nx,1,1)  line concencting corner  1-2
        ! (nx,1,1) -> (nx,ny,1) line concencting corner  2-3
        ! (nx,ny,1)-> (1,ny,1)  line concencting corner  3-4
        ! (1,ny,1) -> (1,1,1)   line concencting corner  4-1

        ! (1,1,1)  -> (1,1,nz)   line concencting corner  1-5
        ! (nx,1,1) -> (nx,1,nz)  line concencting corner  2-5
        ! (nx,ny,1)-> (nx,ny,nz) line concencting corner  3-5
        ! (1,ny,1) -> (1,ny,nz)  line concencting corner  4-5

        ! (1,1,nz)  -> (nx,1,nz)  line concencting corner  5-6
        ! (nx,1,nz) -> (nx,ny,nz) line concencting corner  6-7
        ! (nx,ny,nz)-> (1,ny,nz)  line concencting corner  7-8
        ! (1,ny,nz) -> (1,1,nz)   line concencting corner  9-5
        

        ! (1,1,1)  -> (nx,1,1)  line concencting corner  1-2

        iy=1
        iz=1
        
        do ix=2,nx-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            !idymin  = coordtoindex(ix  ,iy-1,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            !idzmin  = coordtoindex(ix  ,iy  ,iz-1)

            fvec(nshift+id)= factor4*(&
                psi(idxpls) + psi(idypls) + psi(idzpls) -4.0_dp*psi(id) + psi(idxmin) + rhoq(id)*constqW)

        enddo    
        
        ! (nx,1,1) -> (nx,ny,1) line concencting corner  2-3
         
        ix=nx
        iz=1

         do iy=2,ny-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            !idzmin  = coordtoindex(ix  ,iy  ,iz-1)!
            
            fvec(nshift+id)= factor4*(&
                psi(idypls) + psi(idzpls) -4.0_dp*psi(id) + psi(idxmin) + psi(idymin) + rhoq(id)*constqW)

        enddo   

        ! (nx,ny,1)-> (1,ny,1)  line concencting corner  3-4

        iy=ny
        iz=1
        
        do ix=2,nx-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            !idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            !idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            
            fvec(nshift+id)= factor4*(&
                psi(idxpls) + psi(idzpls) -4.0_dp*psi(id) + psi(idxmin) + psi(idymin) + rhoq(id)*constqW)

        enddo

        ! (1,ny,1) -> (1,1,1)   line concencting corner  4-1

        ix=1
        iz=1
        
        do iy=2,ny-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            !idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            !idzmin  = coordtoindex(ix  ,iy  ,iz-1)

            fvec(nshift+id)= factor4*(&
                psi(idxpls) + psi(idypls) + psi(idzpls) -4.0_dp*psi(id) + psi(idymin) + rhoq(id)*constqW)

        enddo

         ! (1,1,1)  -> (1,1,nz)   line concencting corner  1-5
        ix=1
        iy=1

        do iz=2,nz-1
        
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            !idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            !idymin  = coordtoindex(ix  ,iy-1,iz  )
            
            fvec(nshift+id)= factor4*(&
                psi(idxpls) + psi(idypls) + psi(idzpls) -4.0_dp*psi(id) + psi(idzmin) + rhoq(id)*constqW)
        
        enddo
      
       ! (nx,1,1) -> (nx,1,nz)  line concencting corner  2-5

        ix=nx
        iy=1

        do iz=2,nz-1
        
            id      = coordtoindex(ix  ,iy  ,iz  )
            !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            !idymin  = coordtoindex(ix  ,iy-1,iz  )
            
            fvec(nshift+id)= factor4*(&
                psi(idypls) + psi(idzpls) -4.0_dp*psi(id) + psi(idxmin) + psi(idzmin) + rhoq(id)*constqW)

        enddo   
        
        
        ! (nx,ny,1)-> (nx,ny,nz) line concencting corner  3-5
        ix=nx
        iy=ny

        do iz=2,nz-1
        
            id      = coordtoindex(ix  ,iy  ,iz  )
            !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            !idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            
            fvec(nshift+id)= factor4*(&
                psi(idzpls) - 4.0_dp*psi(id) + psi(idxmin) + psi(idymin) + psi(idzmin) + rhoq(id)*constqW)

        enddo

        ! (1,ny,1) -> (1,ny,nz)  line concencting corner  4-5
        ix=1
        iy=ny

        do iz=2,nz-1
        
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            !idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            !idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            
            fvec(nshift+id)= factor4*(&
                psi(idxpls) + psi(idzpls) -4.0_dp*psi(id) + psi(idymin) + psi(idzmin) + rhoq(id)*constqW)

        enddo

         ! (1,1,nz)  -> (nx,1,nz)  line concencting corner  5-6

        iy=1
        iz=nz
        
        do ix=2,nx-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            !idymin  = coordtoindex(ix  ,iy-1,iz  )
            !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            
            fvec(nshift+id)= factor4*(&
                psi(idxpls) + psi(idypls) -4.0_dp*psi(id) + psi(idxmin) + psi(idzmin) + rhoq(id)*constqW)

        enddo 

        ! (nx,1,nz) -> (nx,ny,nz) line concencting corner  6-7

        ix=nx
        iz=nz
        
        do iy=2,ny-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            !idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            
            fvec(nshift+id)= factor4*(&
                psi(idypls) -4.0_dp*psi(id) + psi(idxmin) + psi(idymin) + psi(idzmin)  +rhoq(id)*constqW)

        enddo     


        ! (nx,ny,nz)-> (1,ny,nz)  line concencting corner  7-8
        iy=ny
        iz=nz
        
        do ix=2,nx-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            idxmin  = coordtoindex(ix-1,iy  ,iz  )
            !idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            
            fvec(nshift+id)= factor4*(&
                psi(idxpls) -4.0_dp*psi(id) +psi(idxmin) + psi(idymin) + psi(idzmin) + rhoq(id)*constqW)

        enddo
  
         ! (1,ny,nz) -> (1,1,nz)   line concencting corner  9-5

        ix=1
        iz=nz
        
        do iy=2,ny-1
                    
            id      = coordtoindex(ix  ,iy  ,iz  )
            idxpls  = coordtoindex(ix+1,iy  ,iz  )     
            !idxmin  = coordtoindex(ix-1,iy  ,iz  )
            idypls  = coordtoindex(ix  ,iy+1,iz  )
            idymin  = coordtoindex(ix  ,iy-1,iz  )
            !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
            idzmin  = coordtoindex(ix  ,iy  ,iz-1)
            
            fvec(nshift+id)= factor4*(&
                psi(idxpls) + psi(idypls) -4.0_dp*psi(id) + psi(idymin) + psi(idzmin) + rhoq(id)*constqW)
        enddo

    end subroutine Poisson_Equation_cubic_ST


    ! Set value offset use to the assignmet of the PE and potential EL bc to the iteration vector fvec 
    ! noffset, noffset_bc and noffset_bc_left 
    ! only noffset_bc set sofar !!!!


    subroutine set_offset_Poisson_Equation_Surface
       
        use globals, only : nsize, systype, nsegtypes
        
        ! .. electrostatics: self consistent boundary conditions
        
        select case(systype)
        case ("elect") 
            noffset_bc=4*nsize
        case("electA")  
            noffset_bc=3*nsize
        case("electdouble") 
            noffset_bc=4*nsize
        case("electnopoly") 
            noffset_bc=2*nsize
        case("brush_mul") 
            noffset_bc=(2+nsegtypes)*nsize    
        case("brushdna") 
            noffset_bc=(2+nsegtypes)*nsize    
        case default    
            print*,"error: systype wrong value for Poisson_equation_surface "  
        end select
    

    end subroutine set_offset_Poisson_Equation_Surface

    subroutine Poisson_Equation_Surface(fvec,psi,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

        use globals, only : nsize, LEFT, RIGHT
        use volume, only : nx,ny,nz, coordtoindex

        implicit none

        ! input arguments 
        real(dp), intent(inout)  :: fvec(:)
        real(dp), intent(in)  :: psi(:)
        real(dp), intent(inout) :: psisurfR(:)
        real(dp), intent(inout) :: psisurfL(:)
        real(dp), intent(in)  :: sigmaqSurfR(:)
        real(dp), intent(in)  :: sigmaqSurfL(:)
        character(len=2), intent(in)  :: bcflag(2)

        ! local variables
        integer :: ix, iy
        integer :: idxR, idxL, idxR2D
        integer :: neq_bc, noffset2_bc

        ! .. neeed to be placed  outside

        call set_offset_Poisson_Equation_Surface
       
        neq_bc=0

        if(bcflag(RIGHT)/='cc'  .and. bcflag(RIGHT)/='cp') then
            neq_bc=nx*ny
            do ix=1,nx
                do iy=1,ny
                    idxR = coordtoindex(ix,iy,nz)
                    idxR2D = idxR-(nsize-nx*ny) ! check this 
                    fvec(noffset_bc+idxR2D)=psisurfR(idxR2D)-psi(idxR)-sigmaqSurfR(idxR2D)/2.0_dp
                enddo
            enddo
        else
            do ix=1,nx
                do iy=1,ny
                    idxR = coordtoindex(ix,iy,nz)
                    idxR2D=idxR-(nsize-nx*ny)
                    psisurfR(idxR2D) = psi(idxR)+sigmaqSurfR(idxR2D)/2.0_dp 
                enddo
            enddo
        endif    

        noffset2_bc=noffset_bc +neq_bc

        if(bcflag(LEFT)/='cc'.and. bcflag(LEFT)/='cp') then 
            do ix=1,nx
                do iy=1,ny
                    idxL = coordtoindex(ix,iy,1)
                    fvec(noffset2_bc+idxL)=psi(idxL)-psisurfL(idxL)+sigmaqSurfL(idxL)/2.0_dp
                enddo
            enddo
        else    
            do ix=1,nx
                do iy=1,ny
                    idxL = coordtoindex(ix,iy,1)
                    psisurfL(idxL) = psi(idxL)+sigmaqSurfL(idxL)/2.0_dp
                enddo
            enddo
        endif   

    end subroutine Poisson_Equation_Surface



    subroutine Poisson_Equation_Surface_Eps(fvec,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag,eps)

        use globals, only : nsize, LEFT, RIGHT, systype, nsegtypes
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
        real(dp), intent(in)  :: eps(:)
        character(len=2), intent(in)  :: bcflag(2)

        ! local variables
        integer :: ix, iy, noffset
        integer :: idxR, idxL, idxR2D
        integer :: idzpls, idzmin
        integer :: neq_bc
        real(dp) :: epsxR, epsxL

        ! .. electrostatics: self consistent boundary conditions
        
        noffset=(2+nsegtypes)*nsize   
        
        select case(systype)
        case("brushborn")   
        case default    
            print*,"error: systype wrong value for Poisson_Equation_Surface_Eps "    
        end select

        neq_bc=0

        if(bcflag(RIGHT)/='cc') then
            neq_bc=nx*ny
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,nz-1,idzmin)
                    call linearIndexFromCoordinate(ix,iy,nz,idxR)
                    idxR2D=idxR-(nsize-nx*ny) 
                    epsxR=3.0_dp*eps(idxR)-2.0_dp*eps(idzmin)
                    fvec(noffset+idxR2D)=psisurfR(idxR2D)-psi(idxR)-sigmaqSurfR(idxR2D)*epsxR/2.0_dp
                enddo
            enddo
        else
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,nz-1,idzmin)
                    call linearIndexFromCoordinate(ix,iy,nz,idxR)
                    idxR2D=idxR-(nsize-nx*ny)
                    epsxR=3.0_dp*eps(idxR)-2.0_dp*eps(idzmin)
                    psisurfR(idxR2D) = psi(idxR)+sigmaqSurfR(idxR2D)/2.0_dp 
                enddo
            enddo
        endif    

        noffset=noffset +neq_bc

        if(bcflag(LEFT)/='cc') then 
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,2,idzpls)
                    call linearIndexFromCoordinate(ix,iy,1,idxL)
                    epsxL=3.0_dp*eps(idxL)-2.0_dp*eps(idzpls)
                    fvec(noffset+idxL)=psi(idxL)-psisurfL(idxL)+sigmaqSurfL(idxL)*epsxL/2.0_dp
                enddo
            enddo
        else    
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,2,idzpls)
                    call linearIndexFromCoordinate(ix,iy,1,idxL)
                    epsxL=3.0_dp*eps(idxL)-2.0_dp*eps(idzpls)
                    psisurfL(idxL) = psi(idxL)+sigmaqSurfL(idxL)*epsxL/2.0_dp
                enddo
            enddo
        endif   

    end subroutine Poisson_Equation_Surface_Eps


        
    subroutine Poisson_Pol_Equation(fvec,D2psi,rhoq,rhob)

        use globals, only : nsize
        use parameters, only : constq0

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




    ! computes gradient electrostatic potential multiplied by (scaled) derivate of dielectric 
    ! with respect of volume fraction of polymer. Terms occur in PDF 

    subroutine grad_pot_sqr_eps_cubic(psi,eps,Deps,sqrDpsi)
        
        use volume, only : nx,ny,nz,linearIndexFromCoordinate
        use parameters, only : constqE

        implicit none

        ! input arguments 
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: eps(:)   
        real(dp), intent(in) :: Deps(:)
        real(dp), intent(inout) :: sqrDpsi(:)
        
        ! local variables
        integer :: ix, iy, iz
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        real(dp):: Dpsi(3)

        ! .. electrostatics 
       
        do ix=1,nx
            do iy=1,ny
                do iz=1,nz
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz+1,nz),idzpls)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz-1,nz),idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
    
                    ! grad
                    Dpsi(1)=(psi(idxpls)-psi(idxmin)) !/2.0_dp*delta
                    Dpsi(2)=(psi(idypls)-psi(idymin))
                    Dpsi(3)=(psi(idzpls)-psi(idzmin))
    
                    ! squared and scaled
                    sqrDpsi(id)= constqE*Deps(id)*(Dpsi(1)**2+Dpsi(2)**2+Dpsi(3)**2)
                enddo
            enddo
        enddo    

    end subroutine grad_pot_sqr_eps_cubic
         

    subroutine grad_and_nabla_pot(psi,Dpsi,D2psi,absDpsi,unitdirDpsi)
        
        use volume, only : geometry

        implicit none

        ! input arguments 
        real(dp), intent(in) :: psi(:)
        real(dp), intent(inout) :: Dpsi(:,:),unitdirDpsi(:,:)
        real(dp), intent(inout) :: D2psi(:),absDpsi(:)

        if(geometry=="cubic") then 

            call grad_and_nabla_pot_cubic(psi,Dpsi,D2psi,absDpsi,unitdirDpsi)
        
        else if (geometry=="prism") then 
        
            call grad_and_nabla_pot_prism(psi,Dpsi,D2psi,absDpsi,unitdirDpsi)
        
        endif

    end subroutine grad_and_nabla_pot

        

    ! computes gradient, nalba (double derivative) , absolute value of gradient and unit direction of gradient of 
    ! electrostaic potential

    subroutine grad_and_nabla_pot_cubic(psi,Dpsi,D2psi,absDpsi,unitdirDpsi) 
        
        use globals, only : nsize
        use volume, only : nx,ny,nz,delta,linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(in) :: psi(:)
        real(dp), intent(inout) :: Dpsi(:,:),unitdirDpsi(:,:)
        real(dp), intent(inout) :: D2psi(:),absDpsi(:)
    
        
        ! local variables
        integer :: ix, iy, iz
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        do ix=1,nx
            do iy=1,ny
                do iz=1,nz
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz+1,nz),idzpls)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz-1,nz),idzmin)
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

    subroutine grad_and_nabla_pot_prism(psi,Dpsi,D2psi,absDpsi,unitdirDpsi) ! sigmaqSurfR,sigmaqSurfL)
        
        use globals, only : nsize
        use volume, only : nx,ny,nz,delta,linearIndexFromCoordinate
        use volume, only : cos_two_beta, sin_two_beta
        
        implicit none

        ! input arguments 
        real(dp), intent(in) :: psi(:)
        real(dp), intent(inout) :: Dpsi(:,:),unitdirDpsi(:,:)
        real(dp), intent(inout) :: D2psi(:),absDpsi(:)
        !real(dp), intent(in) :: sigmaqSurfR(:)
        !real(dp), intent(in) :: sigmaqSurfL(:)
        
        ! local variables
        integer :: ix, iy, iz
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin 
        integer :: idxypls, idxymin, idxplsymin, idxminypls
        real(dp) :: Dpsiu,Dpsiv

        ! .. electrostatics  

        do ix=1,nx
            do iy=1,ny
                do iz=1,nz
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz+1,nz),idzpls)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz-1,nz),idzmin)
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
                do iz=1,nz
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz+1,nz),idzpls)
                    call linearIndexFromCoordinate(ix,           iy,ipbc(iz-1,nz),idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)

                    rhob(id)=-(sqrtcos*((electPol(1,idxpls)-electPol(1,idxmin))+ &
                        (electPol(2,idypls)-electPol(2,idymin)))+     &
                        (electPol(3,idzpls)-electPol(3,idzmin)))/(2.0_dp*delta)    
                enddo
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
