! module computate divergence of diffusive (flux div_flux)
! and related quantities

module flux

    use precision_definition

    implicit none

    real(dp), allocatable, dimension(:,:) :: divJ    ! divergence of flux range divJ(nsize,niontypes)
    real(dp), allocatable, dimension(:)   :: mu      ! chemical potential range mu(nsize)
    real(dp), allocatable, dimension(:,:,:) :: Jvec  ! vector flux in  range J(nsize,niontypes)
    real(dp), allocatable, dimension(:,:) :: mu_ion     ! chemical potential range mu_ion(nsize,niontypes) 

    public :: divJ, mu, Jvec, mu_ion
    public :: div_flux, allocate_divJ, allocate_mu

contains

    subroutine allocate_divJ()

        use globals, only : nsize
        use parameters, only : niontypes

        allocate(divJ(nsize,niontypes))

    end subroutine allocate_divJ

    subroutine allocate_mu()

        use globals, only : nsize
        use parameters, only : niontypes

        allocate(mu(nsize))

    end subroutine allocate_mu

    subroutine allocate_Jvec()

        use globals, only : nsize
        use parameters, only : niontypes

        allocate(Jvec(nsize,3,niontypes))

    end subroutine allocate_Jvec

    subroutine allocate_mu_ion()

        use globals, only : nsize
        use parameters, only : niontypes

        allocate(mu_ion(nsize,niontypes))

    end subroutine allocate_mu_ion


   ! position dependent chem pot  
   ! \beta \mu_i(r) = \ln(\rho_i(r) v_w) + \beta \pi(r) * v_i + \beta * q_i \psi(r)
   ! \beta \mu_i(r) = \ln(x_i(r) v_w/v_i) - \ln(x_w(r)) * v_i/v_w +  z_i *  e  * \beta \psi(r)
   ! input : real(dp) xsol,xvol,psi,iontype)    

    subroutine chem_potential(mu,xsol,xvol,psi,iontype)

        use globals, only : nsize
        use parameters, only :  vol, zval
        use molecules, only : get_value_moleclist

        ! input arguments 
        real(dp), intent(inout) :: mu(:)
        real(dp), intent(in) :: xvol(:),xsol(:)
        real(dp), intent(in) :: psi(:)
        character(len=*), intent(in)  :: iontype

       
        ! local variables

        integer :: idx
        character(len=5):: key
        real(dp) :: volum, valence
        
        key = trim(iontype)
        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        
        if(.true.) then
           print*,"key=",key  
           print*,"volum=",volum
           print*,"valence=",valence
        endif

        do idx=1,nsize
            mu(idx) = log(xvol(idx)/volum) -log(xsol(idx)) * volum + valence * psi(idx)
        enddo        
 
    end subroutine chem_potential


    subroutine volumefraction(xvol,xsol,mu,psi,volum,valence )

        use globals, only :nsize

        ! input arguments 
        real(dp), intent(inout) :: xvol(:)
        real(dp), intent(in) :: mu(:),xsol(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: volum
        integer, intent(in) :: valence
    
        ! local variables
        integer ::  idx
        
        do idx=1,nsize
            xvol(idx)= volum * exp(  mu(idx) - valence * psi(idx)) *  ( xsol(idx) ** volum) 
        enddo        
 
    end subroutine volumefraction

    subroutine div_flux(divJ,xsol,xvol,psi,iontype)

        use volume, only : geometry

        ! input arguments 
        real(dp), intent(inout) :: divJ(:)
        real(dp), intent(in) :: xvol(:), xsol(:)
        real(dp), intent(in) :: psi(:)
        character(len=*), intent(in)  :: iontype

        if(geometry=="cubic") then

            call div_flux_cubic(divJ,xsol,xvol,psi,iontype)
        
        else if (geometry=="prism") then 

            print*,"Prism geometry not defined for div_flux "  
        
        else 
        
            print*,"Wrong geometry in div_flux "    
        
        endif
                
    end subroutine div_flux
        
    ! Computes div. J with numerical scheme using Gauss's theorem 

    subroutine div_flux_cubic(divJ,xsol,xvol,psi,iontype)

        use globals, only : nsize, DEBUG_ST
        use volume, only : nx, ny, nz,  coordtoindex
        use parameters, only : vol, zval
        use parameters, only : mumin, mumax, xvolmin, xvolmax
        use molecules, only : get_value_moleclist

        ! input arguments 

        real(dp), intent(inout) :: divJ(:)
        real(dp), intent(in) :: xvol(:), xsol(:)
        real(dp), intent(in) :: psi(:)
        character(len=*), intent(in)  :: iontype
        
        ! local variables

        integer :: ix, iy, iz
        integer ::  id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        character(len=5):: key
        real(dp) :: volum, valence, coeff_scaled
        real(dp) :: mu_zmin, mu_zpls, xvol_zmin, xvol_zpls
        real(dp) :: mu(nsize)
        real(dp) :: Jdotxpls, Jdotxmin, Jdotypls, Jdotymin, Jdotzpls, Jdotzmin, divJtmp 

        key = trim(iontype)

        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        mu_zmin = get_value_moleclist(mumin,key)
        mu_zpls = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zpls = get_value_moleclist(xvolmax,key)
        
        !coeff_scaled = 1.0_dp/(volum*delta*2.0_dp)
        coeff_scaled =  1.0_dp

        !  volum from conversion of volumefraction to density  
        !  in cylinder coordinates : 2 delta^2 from division by 2 \pi delta^2 {\bar r}_i 
        !  in cubic coordiantes :  2 delta  :  from interpolation of density

        if(DEBUG_ST) then
           print*,"key=",key  
           print*,"volum=",volum
           print*,"valence=",valence
           print*,"mu_zmin=",mu_zmin
           print*,"mu_zpls=",mu_zpls
           print*,"xvol_zpls=",xvol_zpls
           print*,"xvol_zmin=",xvol_zmin
           print*,"size(xvol)=",size(xvol)
        endif

        do id=1,nsize
            mu(id) = log(xvol(id)/volum) -log(xsol(id))*volum + valence * psi(id)
        enddo     
        ! call chem_potential(mu, xsol, xvol, psi, iontype) ! same as above 

        divJ=0.0_dp
        if(DEBUG_ST) divJ= 123435600.0000_dp ! used to detect unassinged values of divJ 
        
        ! inside 
        ! periodic bc  for x and y planes 
        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    
                    id      = coordtoindex(ix  ,iy  ,iz  )
                    idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                    idxmin  = coordtoindex(ix-1,iy  ,iz  )
                    idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                    idzmin  = coordtoindex(ix  ,iy  ,iz-1)
                    idypls  = coordtoindex(ix  ,iy+1,iz  )
                    idymin  = coordtoindex(ix  ,iy-1,iz  )

                    Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                    Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                    Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                    Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
                    Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
                    Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

                    divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                    divJ(id) = - coeff_scaled * divJtmp

                enddo
            enddo
        enddo    

        ! faces :

        ! bottom z-face
        ! boundary z=0 iz=1 
        ! xvol_zmin and mu_zmin resevoir  conditions
        iz = 1
        do iy=2,ny-1
            do ix=2,nx-1

                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
               ! idzmin  = coordtoindex(ix  ,iy  ,iz-1)
                idypls  = coordtoindex(ix  ,iy+1,iz  )                    
                idymin  = coordtoindex(ix  ,iy-1,iz  )
                

                Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))

                Jdotzpls = (xvol(idzpls)+ xvol(id)     )*(mu(idzpls) - mu(id)    )
                Jdotzmin = (xvol(id)    + xvol_zmin    )*(mu(id)     - mu_zmin   )

                divJtmp  = Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                
                divJ(id) = - coeff_scaled * divJtmp

            enddo
        enddo
       
        ! top z-face
        ! boundary z=nz*delta iz=nz
        ! xvol_zmax and mu_zmax resevoir conditions
        
        iz = nz
        do iy=2,ny-1
            do ix=2,nx-1

                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
               !   idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1)
                idypls  = coordtoindex(ix  ,iy+1,iz  )                    
                idymin  = coordtoindex(ix  ,iy-1,iz  )
            
                Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))

                Jdotzpls = (xvol_zpls    + xvol(id)     )*(mu_zpls   - mu(id)    )
                Jdotzmin = (xvol(id)     + xvol(idzmin) )*(mu(id)    - mu(idzmin))

                divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                
                divJ(id) = - coeff_scaled * divJtmp
            
            enddo    
        enddo

        ! boundary x= 0 plane ix=1  
        
        ix=1
        do iy=2,ny-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
               ! idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )

                Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
                Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
                Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

                divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                divJ(id) = - coeff_scaled * divJtmp

            enddo
        enddo    

        ! boundary x= nx delta  plane ix=nx  
        
        ix=nx
        do iy=2,ny-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
               ! idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )

                Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
                Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
                Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

                divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                divJ(id) = - coeff_scaled * divJtmp

            enddo
        enddo    

        ! boundary y= 0 plane iy=1  
        
        iy=1
        do ix=2,nx-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                !idymin  = coordtoindex(ix  ,iy-1,iz  )

                Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                Jdotymin = 0.0_dp !(xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
                Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
                Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

                divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                divJ(id) = - coeff_scaled * divJtmp

            enddo
        enddo    

        ! boundary y= ny delta  plane iy=ny  
        
        iy=ny
        do ix=2,nx-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                !idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )

                Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                Jdotypls = 0.0_dp ! (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
                Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
                Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

                divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                divJ(id) = - coeff_scaled * divJtmp

            enddo
        enddo    


        ! corners 

        ! ix=1  iy=1  iz=1 : c1
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

        Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = 0.0_dp !(xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = 0.0_dp ! (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
        !Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
        Jdotzmin = (xvol(id)    + xvol_zmin    )*(mu(id)     - mu_zmin  )

        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp


        ! ix=nx iy=1  iz=1 : c2
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

        Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = 0.0_dp ! (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
        !Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
        Jdotzmin = (xvol(id)    + xvol_zmin    )*(mu(id)     - mu_zmin  )

        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp

        
        ! ix=nx iy=ny iz=1 : c3
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

        Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = 0.0_dp !(xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
        !Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
        Jdotzmin = (xvol(id)    + xvol_zmin    )*(mu(id)     - mu_zmin  )

        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp
        
        ! ix=1  iy=ny iz=1 : c4
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

        Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = 0.0_dp !(xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )                    
        !Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
        Jdotzmin = (xvol(id)    + xvol_zmin    )*(mu(id)     - mu_zmin  )

        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp
        
        ! ix=1  iy=1  iz=nz : c5
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

        Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = 0.0_dp ! (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls   - mu(id)    )                   
        Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
       
        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp

        
        ! ix=nx iy=1  iz=nz : c6
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

        Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = 0.0_dp ! (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls   - mu(id)    )                   
        Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
       
        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp


        ! ix=nx iy=ny iz=nz : c7
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

        Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = 0.0_dp ! (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls   - mu(id)    )                   
        Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
       
        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp


        
        ! ix=1  iy=ny iz=nz : c8
        ix = 1
        iy = ny
        iz = nz
        id      = coordtoindex(ix  ,iy  ,iz  )
        idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        !idxmin  = coordtoindex(ix-1,iy  ,iz  )
        !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        !idypls  = coordtoindex(ix  ,iy+1,iz  )
        idymin  = coordtoindex(ix  ,iy-1,iz  )

        Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
        Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
        Jdotypls = 0.0_dp ! (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
        Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
        Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls   - mu(id)    )                   
        Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))
       
        divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
        divJ(id) = - coeff_scaled * divJtmp

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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = 0.0_dp !  (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            !Jdotzmin = 0.0_dp ! (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin)) 
            Jdotzmin = (xvol(id)     + xvol_zmin   )*(mu(id)     - mu_zmin)

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

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
            !idzmin  = coordtoindex(ix  ,iy  ,iz-1)

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = 0.0_dp !(xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            Jdotzmin = (xvol(id)     + xvol_zmin   )*(mu(id)     - mu_zmin   )

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = 0.0_dp ! (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            !Jdotzmin = 0.0_dp ! (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin)) 
            Jdotzmin = (xvol(id)     + xvol_zmin   )*(mu(id)     - mu_zmin   )

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            !Jdotzmin = 0.0_dp ! (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin)) 
            Jdotzmin = (xvol(id)     + xvol_zmin   )*(mu(id)     - mu_zmin   )

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
                
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = 0.0_dp ! (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
                
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp
        
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

            Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
                
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = 0.0_dp ! (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
                
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp
        
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

            Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
                
            Jdotypls = 0.0_dp ! (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
                
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp
        
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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
                
            Jdotypls = 0.0_dp ! (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
                
            Jdotzpls = (xvol(idzpls) + xvol(id)    )*(mu(idzpls) - mu(id)    )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin))

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp
        
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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = 0.0_dp ! (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls      - mu(id)  )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin)) 
            

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

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

            Jdotxpls = 0.0_dp ! (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls      - mu(id)  )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin)) 
            

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = 0.0_dp ! (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls    - mu(id)    )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin)) 
        
            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

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

            Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
            Jdotxmin = 0.0_dp ! (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
            
            Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
            Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))
               
            Jdotzpls = (xvol_zpls    + xvol(id)    )*(mu_zpls    - mu(id)    )
            Jdotzmin = (xvol(id)     + xvol(idzmin))*(mu(id)     - mu(idzmin)) 
            

            divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
            divJ(id) = - coeff_scaled * divJtmp

        enddo    
        
        if(DEBUG_ST) then 
            do id=1,nsize    
                if(divJ(id)== 123435600.0000_dp) print*,"divJ unassiged in id=",id
            enddo     
        endif    

    end subroutine div_flux_cubic


    ! Computes flux J= - D rho * grad mu 

    subroutine fluxJ(J,xsol,xvol,psi,iontype)

        use globals, only : nsize, DEBUG_ST
        use parameters, only : vol,  Diffcoeff, vsol
        use molecules, only : get_value_moleclist
        use volume, only : indextocoord
        
        ! input arguments 

        real(dp), intent(inout) :: J(:,:)
        real(dp), intent(in) :: xvol(:), xsol(:)
        real(dp), intent(in) :: psi(:)
        character(len=*), intent(in)  :: iontype

        ! local arguments

        real(dp) :: grad_mu(nsize,3)
        real(dp) ::  volum, Diffconst, J0
        integer :: i, id, ix, iy, iz, idx 
        character(len=5):: key
 
        call grad_chem_pot(grad_mu,xsol,xvol,psi,iontype)

        key = trim(iontype)
        volum = get_value_moleclist(vol,key)
        Diffconst = get_value_moleclist(Diffcoeff,key)

        ! pre factor in flux factor (1.0e-9)^2 arise from vsol in unit of nm and grad_mu in 1/nm
        ! J in 1/( nm^2) s)
        
        J0 = - Diffconst / ( volum * vsol * (1.0e-9_dp)**2) 
        
        do i=1,3
            do id=1,nsize
                J(id,i) = J0 * xvol(id) * grad_mu(id,i) 
            enddo
        enddo    

        if(DEBUG_ST) then
            if(key=="K") then 
                do idx=1,nsize
                    ix=indextocoord(idx,1)
                    iy=indextocoord(idx,2)
                    iz=indextocoord(idx,3)
                    write(200,*)ix,iy,iz,J(idx,1),J(idx,2),J(idx,3),xvol(idx)
                    write(300,*)ix,iy,iz,grad_mu(idx,1),grad_mu(idx,2),grad_mu(idx,3)
                enddo    
            endif
        endif        


    end subroutine fluxJ 

    ! Computes  total current through plane z=nz/2 delta 

    function current_I() result(currI)

        use globals, only : nsize
        use parameters, only : niontypes, isionselfconsistent, iontype
        use field, only      : xsol,psi
        use field, only      : xNa,xCl,xK,xHplus,xOHmin,xFe2,xFe3,xMg,xCa
        
        ! return arguments 

        real(dp) :: currI

        ! local arguments

        real(dp) :: J(nsize,3)
        integer :: t
        real(dp) :: current(niontypes)
        
        current = 0.0_dp

        do t=1,niontypes
            if(isionselfconsistent(t)) then
                if(iontype(t)=="Na") then 
                    call fluxJ(J,xsol,xNa,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif    
                if(iontype(t)=="K") then 
                    call fluxJ(J,xsol,xK,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif 
                if(iontype(t)=="Cl") then    
                    call fluxJ(J,xsol,xCl,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif
                if(iontype(t)=="Hplus") then
                     call fluxJ(J,xsol,xHplus,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif    
                if(iontype(t)=="OHmin") then 
                    call fluxJ(J,xsol,xOHmin,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif
                if(iontype(t)=="Mg") then 
                    call fluxJ(J,xsol,xMg,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif 
                if(iontype(t)=="Fe2") then 
                    call fluxJ(J,xsol,xFe2,psi,iontype(t))
                    current(t) = current_I_ion(J,iontype(t))
                endif
                if(iontype(t)=="Fe3") then 
                    call fluxJ(J,xsol,xFe3,psi,iontype(t)) 
                    current(t) = current_I_ion(J,iontype(t))
                endif   
                if(iontype(t)=="Ca") then  
                    call fluxJ(J,xsol,xCa,psi,iontype(t)) 
                    current(t) = current_I_ion(J,iontype(t))
                endif     
            endif    
        enddo

        currI = sum(current)

    end function current_I 


    ! Computes  current through plane z=nz/2 delta  for given flux density J 
    ! for given charge velence zval 

    function current_I_ion(J,iontype) result(current)

        use globals, only : nsize
        use parameters, only : zval
        use molecules, only : get_value_moleclist
        use volume, only : nz, delta, coordtoindex
        
        ! input arguments 

        real(dp), intent(inout) :: J(:,:)
        character(len=*), intent(in)  :: iontype

        ! return arguments

        real(dp) :: current

        ! local variables

        character(len=5):: key
        real(dp) :: valence, sum_curJ
        integer ix, iy, izmidplane, id

        key = trim(iontype)
        valence = get_value_moleclist(zval,key)
        izmidplane = int(nz/2.0_dp)
        sum_curJ = 0.0_dp

        do ix=1,nsize
            do iy=1,nsize
                id = coordtoindex(ix,iy,izmidplane)
                sum_curJ = sum_curJ + J(id,3) 
            enddo
        enddo    
        
        current = sum_curJ *( delta **2) * valence
        
    end function current_I_ion

    ! computes grad of chemical potential

    subroutine grad_chem_pot(grad_mu,xsol,xvol,psi,iontype)

        use globals, only : nsize, DEBUG_ST
        use volume, only : delta, nx, ny, nz,  coordtoindex
        use parameters, only : vol, zval
        use parameters, only : mumin, mumax, xvolmin, xvolmax
        use molecules, only : get_value_moleclist

        ! input arguments 

        real(dp), intent(inout) :: grad_mu(:,:)
        real(dp), intent(in) :: xvol(:), xsol(:)
        real(dp), intent(in) :: psi(:)
        character(len=*), intent(in)  :: iontype

        ! local variables

        integer :: ix, iy, iz, k
        integer ::  id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        character(len=5):: key
        real(dp) :: volum, valence
        real(dp) :: mu_zmin, mu_zpls, xvol_zmin, xvol_zpls
        real(dp) :: mu(nsize)

        key = trim(iontype)

        volum   = get_value_moleclist(vol,key)
        valence = get_value_moleclist(zval,key)
        mu_zmin = get_value_moleclist(mumin,key)
        mu_zpls = get_value_moleclist(mumax,key)
        xvol_zmin = get_value_moleclist(xvolmin,key)
        xvol_zpls = get_value_moleclist(xvolmax,key)
        
        if(DEBUG_ST) grad_mu = 1234567.89_dp

        do id=1,nsize
            mu(id) = log(xvol(id)/volum) -log(xsol(id))*volum + valence * psi(id)
        enddo     
       
        ! inside 

        do iz=2,nz-1
            do iy=2,ny-1
                do ix=2,nx-1
                    
                    id      = coordtoindex(ix  ,iy  ,iz  )
                    idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                    idxmin  = coordtoindex(ix-1,iy  ,iz  )
                    idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                    idzmin  = coordtoindex(ix  ,iy  ,iz-1)
                    idypls  = coordtoindex(ix  ,iy+1,iz  )
                    idymin  = coordtoindex(ix  ,iy-1,iz  )

                    grad_mu(id,1)  = mu(idxpls) - mu(idxmin)
                    grad_mu(id,2)  = mu(idypls) - mu(idymin)
                    grad_mu(id,3)  = mu(idzpls) - mu(idzmin)
                   
                enddo
            enddo
        enddo    

        ! faces :

        ! bottom z-face
        ! boundary z=0 iz=1 
        ! xvol_zmin and mu_zmin resevoir  conditions

        iz = 1
        do iy=2,ny-1
            do ix=2,nx-1

                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                ! idzmin = coordtoindex(ix  ,iy  ,iz- 1)
                idypls  = coordtoindex(ix  ,iy+1,iz  )                    
                idymin  = coordtoindex(ix  ,iy-1,iz  )

                grad_mu(id,1)  = (mu(idxpls) - mu(idxmin))
                grad_mu(id,2)  = (mu(idypls) - mu(idymin))
                grad_mu(id,3)  = (mu(idzpls) - mu_zmin)
                   
            enddo
        enddo
       
        ! top z-face
        ! boundary z=nz*delta iz=nz
        ! xvol_zmax and mu_zmax resevoir conditions
        
        iz = nz
        do iy=2,ny-1
            do ix=2,nx-1

                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                ! idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1)
                idypls  = coordtoindex(ix  ,iy+1,iz  )                    
                idymin  = coordtoindex(ix  ,iy-1,iz  )
            
                grad_mu(id,1)  = (mu(idxpls) - mu(idxmin))
                grad_mu(id,2)  = (mu(idypls) - mu(idymin))
                grad_mu(id,3)  = (mu_zpls    - mu(idzmin))    
            
            enddo    
        enddo

        ! boundary x= 0 plane ix=1  
        
        ix=1
        do iy=2,ny-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  ) 
                ! idxmin  = coordtoindex(ix-1,iy  ,iz  )    
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )
            
                grad_mu(id,1)  = 0.0_dp
                grad_mu(id,2)  = (mu(idypls) - mu(idymin))
                grad_mu(id,3)  = (mu(idzpls) - mu(idzmin)) 

            enddo
        enddo    

        ! boundary x= nx delta  plane ix=nx  
        
        ix=nx
        do iy=2,ny-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
               ! idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )

                grad_mu(id,1)  = 0.0_dp
                grad_mu(id,2)  = (mu(idypls) - mu(idymin))
                grad_mu(id,3)  = (mu(idzpls) - mu(idzmin))

            enddo
        enddo    

        ! boundary y= 0 plane iy=1  
        
        iy=1
        do ix=2,nx-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                idypls  = coordtoindex(ix  ,iy+1,iz  )
                !idymin  = coordtoindex(ix  ,iy-1,iz  )

                grad_mu(id,1)  = (mu(idxpls) - mu(idxmin))
                grad_mu(id,2)  = 0.0_dp
                grad_mu(id,3)  = (mu(idzpls) - mu(idzmin))

            enddo
        enddo    

        ! boundary y= ny delta  plane iy=ny  
        
        iy=ny
        do ix=2,nx-1
            do iz=2,nz-1
                  
                id      = coordtoindex(ix  ,iy  ,iz  )
                idxpls  = coordtoindex(ix+1,iy  ,iz  )     
                idxmin  = coordtoindex(ix-1,iy  ,iz  )
                idzpls  = coordtoindex(ix  ,iy  ,iz+1)
                idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
                !idypls  = coordtoindex(ix  ,iy+1,iz  )
                idymin  = coordtoindex(ix  ,iy-1,iz  )

                grad_mu(id,1)  = (mu(idxpls) - mu(idxmin))
                grad_mu(id,2)  = 0.0_dp
                grad_mu(id,3)  = (mu(idzpls) - mu(idzmin))

            enddo
        enddo    

        ! corners 

        ! ix=1  iy=1  iz=1 : c1
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

        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu(idzpls) - mu_zmin)


        ! ix=nx iy=1  iz=1 :  c2 
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

        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu(idzpls) - mu_zmin)

        ! ix=nx iy=ny iz=1 : c3
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

        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu(idzpls) - mu_zmin)

        ! ix=1  iy=ny iz=1 : c4
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

        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu(idzpls) - mu_zmin)
        
        ! ix=1  iy=1  iz=nz : c5
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
        
        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu_zpls - mu(idzmin))
        
        ! ix=nx iy=1  iz=nz : c6
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

        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu_zpls - mu(idzmin))
       

        ! ix=nx iy=ny iz=nz : c7
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

        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu_zpls - mu(idzmin))
        
        ! ix=1  iy=ny iz=nz : c8
        ix = 1
        iy = ny
        iz = nz

        id      = coordtoindex(ix  ,iy  ,iz  )
        idxpls  = coordtoindex(ix+1,iy  ,iz  )     
        !idxmin  = coordtoindex(ix-1,iy  ,iz  )
        !idzpls  = coordtoindex(ix  ,iy  ,iz+1)
        idzmin  = coordtoindex(ix  ,iy  ,iz-1) 
        !idypls  = coordtoindex(ix  ,iy+1,iz  )
        idymin  = coordtoindex(ix  ,iy-1,iz  )

        grad_mu(id,1)  = 0.0_dp
        grad_mu(id,2)  = 0.0_dp
        grad_mu(id,3)  = (mu_zpls - mu(idzmin))


        ! edges

        ! (1,1,1)  -> (nx,1,1)  line concencting corner  1-2
        ! (nx,1,1) -> (nx,ny,1) line concencting corner  2-3
        ! (nx,ny,1)-> (1,ny,1)  line concencting corner  3-4
        ! (1,ny,1) -> (1,1,1)   line concencting corner  4-1

        ! (1,1,1)  -> (1,1,nz)   line concencting corner  1-5
        ! (nx,1,1) -> (nx,1,nz)  line concencting corner  2-6
        ! (nx,ny,1)-> (nx,ny,nz) line concencting corner  3-7
        ! (1,ny,1) -> (1,ny,nz)  line concencting corner  4-8

        ! (1,1,nz)  -> (nx,1,nz)  line concencting corner  5-6
        ! (nx,1,nz) -> (nx,ny,nz) line concencting corner  6-7
        ! (nx,ny,nz)-> (1,ny,nz)  line concencting corner  7-8
        ! (1,ny,nz) -> (1,1,nz)   line concencting corner  8-5
      

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

            grad_mu(id,1)  = (mu(idxpls) - mu(idxmin))
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu(idzpls)  - mu_zmin)   

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
            !idzmin  = coordtoindex(ix  ,iy  ,iz-1)

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = (mu(idypls) - mu(idymin))
            grad_mu(id,3)  = (mu(idzpls) - mu_zmin)    

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

            grad_mu(id,1)  = (mu(idxpls) - mu(idxmin))
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu(idzpls) - mu_zmin)   

            

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

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = (mu(idypls) - mu(idymin))
            grad_mu(id,3)  = (mu(idzpls) - mu_zmin)

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

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu(idzpls) - mu(idzmin))
        
        enddo
        
        ! (nx,1,1) -> (nx,1,nz)  line concencting corner  2-6

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

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu(idzpls) - mu(idzmin))
        
        enddo


        ! (nx,ny,1)-> (nx,ny,nz) line concencting corner  3-7
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

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu(idzpls) - mu(idzmin))
        
        enddo


        ! (1,ny,1) -> (1,ny,nz)  line concencting corner  4-8
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

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu(idzpls) - mu(idzmin))
        
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

            grad_mu(id,1)  = (mu(idxpls) - mu(idxmin))
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu_zpls - mu(idzmin))

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

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = (mu(idypls) - mu(idymin))
            grad_mu(id,3)  = (mu_zpls - mu(idzmin))
            
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

            grad_mu(id,1)  = mu(idxpls) - mu(idxmin)
            grad_mu(id,2)  = 0.0_dp
            grad_mu(id,3)  = (mu_zpls - mu(idzmin))

        enddo    

        ! (1,ny,nz) -> (1,1,nz)   line concencting corner  8-5

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

            grad_mu(id,1)  = 0.0_dp
            grad_mu(id,2)  = mu(idypls) - mu(idymin)
            grad_mu(id,3)  = mu_zpls - mu(idzmin)

        enddo    

        ! divided by 2 delta 

        do k=1,3
            do id=1,nsize 
                grad_mu(id,k)= grad_mu(id,k)/(2.0_dp*delta)
            enddo
        enddo    


    end subroutine  grad_chem_pot

    ! computes flux J

    subroutine calculate_fluxJ() 

        use parameters, only : niontypes, isionselfconsistent, iontype
        use field, only  : xsol,psi
        use field, only  : xNa,xCl,xK,xHplus,xOHmin,xFe2,xFe3,xMg,xCa
    
        ! local arguments
        integer :: t

        Jvec = 0.0_dp

        do t=1,niontypes
            if(isionselfconsistent(t)) then
                if(iontype(t)=="Na")    call fluxJ(Jvec(:,:,t),xsol,xNa,  psi,iontype(t))
                if(iontype(t)=="K")     call fluxJ(Jvec(:,:,t),xsol,xK,   psi,iontype(t))
                if(iontype(t)=="Cl")    call fluxJ(Jvec(:,:,t),xsol,xCl,  psi,iontype(t))
                if(iontype(t)=="Hplus") call fluxJ(Jvec(:,:,t),xsol,xHplus,psi,iontype(t))
                if(iontype(t)=="OHmin") call fluxJ(Jvec(:,:,t),xsol,xOHmin,psi,iontype(t))
                if(iontype(t)=="Mg")    call fluxJ(Jvec(:,:,t),xsol,xMg,   psi,iontype(t))
                if(iontype(t)=="Fe2")   call fluxJ(Jvec(:,:,t),xsol,xFe2,  psi,iontype(t))
                if(iontype(t)=="Fe3")   call fluxJ(Jvec(:,:,t),xsol,xFe3,  psi,iontype(t)) 
                if(iontype(t)=="Ca")    call fluxJ(Jvec(:,:,t),xsol,xCa,   psi,iontype(t)) 
            endif    
        enddo

    end subroutine calculate_fluxJ


    subroutine calculate_mu_ion() 

        use parameters, only : niontypes, isionselfconsistent, iontype
        use field, only  : xsol,psi
        use field, only  : xNa,xCl,xK,xHplus,xOHmin,xFe2,xFe3,xMg,xCa
    
        ! local arguments
        integer :: t

        mu_ion = 0.0_dp

        do t=1,niontypes
            if(isionselfconsistent(t)) then
                                
                if(iontype(t)=="Na")    call chem_potential(mu_ion(:,t),xsol,xNa,  psi,iontype(t))
                if(iontype(t)=="K")     call chem_potential(mu_ion(:,t),xsol,xK,   psi,iontype(t))
                if(iontype(t)=="Cl")    call chem_potential(mu_ion(:,t),xsol,xCl,  psi,iontype(t))
                if(iontype(t)=="Hplus") call chem_potential(mu_ion(:,t),xsol,xHplus,psi,iontype(t))
                if(iontype(t)=="OHmin") call chem_potential(mu_ion(:,t),xsol,xOHmin,psi,iontype(t))
                if(iontype(t)=="Mg")    call chem_potential(mu_ion(:,t),xsol,xMg,   psi,iontype(t))
                if(iontype(t)=="Fe2")   call chem_potential(mu_ion(:,t),xsol,xFe2,  psi,iontype(t))
                if(iontype(t)=="Fe3")   call chem_potential(mu_ion(:,t),xsol,xFe3,  psi,iontype(t)) 
                if(iontype(t)=="Ca")    call chem_potential(mu_ion(:,t),xsol,xCa,   psi,iontype(t)) 
                
            endif    
        enddo

    end subroutine calculate_mu_ion

    


end module flux
