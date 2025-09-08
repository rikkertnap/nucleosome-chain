module flux
    use precision_definition

    implicit none

    public :: div_flux

contains
   

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
        

    subroutine div_flux_cubic(divJ,xsol,xvol,psi,iontype)

        use globals, only : nsize, DEBUG_ST
        use volume, only : nx, ny, nz,  coordtoindex, delta
        use parameters, only : vol, zval
        use parameters, only : mumin, mumax, xvolmin, xvolmax
        use molecules, only : get_value_moleclist
        use volume, only : ipbc

        ! input arguments 

        real(dp), intent(inout) :: divJ(:)
        real(dp), intent(in) :: xvol(:), xsol(:)
        real(dp), intent(in) :: psi(:)
        character(len=*), intent(in)  :: iontype
        
        ! local variables

        integer :: ix, iy, iz
        integer :: idx, id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
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
        coeff_scaled = 1.0_dp

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
        endif

        do idx=1,nsize
            mu(id) = log(xvol(id)/volum) -log(xsol(id))*volum +valence * psi(id)
        enddo     

        !  call chem_potential(mu, xsol, xvol, psi, iontype)


        ! inside 
        ! periodic bc  for x and y planes 

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

        ! boundary z=0 iz=1 
        ! xvol_zmin and mu_zmin resevoir  conditions

        do iy=1,ny
            do ix=1,nx

                iz = 1
                
                id      = coordtoindex(ix,           iy,iz)
                idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)    
                idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                idzpls  = coordtoindex(ix,           iy,iz+1)
                idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)

                Jdotxpls = (xvol(idxpls) + xvol(id)    )*(mu(idxpls) - mu(id)    )
                Jdotxmin = (xvol(id)     + xvol(idxmin))*(mu(id)     - mu(idxmin))
               
                Jdotypls = (xvol(idypls) + xvol(id)    )*(mu(idypls) - mu(id))
                Jdotymin = (xvol(id)     + xvol(idymin))*(mu(id)     - mu(idymin))

                Jdotzpls = (xvol(idzpls)+ xvol(id)     )*(mu(idzpls) - mu(id)    )
                Jdotzmin = (xvol(id)    + xvol_zmin    )*(mu(id)     - mu_zmin  )

                divJtmp  =  Jdotxpls - Jdotxmin + Jdotypls - Jdotymin + Jdotzpls - Jdotzmin
                
                divJ(id) = - coeff_scaled * divJtmp

            enddo
        enddo
       

        ! boundary z=nz*delta iz=nz
        ! xvol_zmax and mu_zmax resevoir conditions

        do iy=1,ny
            do ix=1,nx

                iz = nz
                
                id      = coordtoindex(ix,           iy,iz)
                idxpls  = coordtoindex(ipbc(ix+1,nx),iy,iz)    
                idxmin  = coordtoindex(ipbc(ix-1,nx),iy,iz)
                idzmin  = coordtoindex(ix,           iy,iz-1)
                idypls  = coordtoindex(ix,ipbc(iy+1,ny),iz)
                idymin  = coordtoindex(ix,ipbc(iy-1,ny),iz)
            
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
          
       
    end subroutine div_flux_cubic

end module
