! module constains auxilary  function need by fcn function fcnnucl_ionbin_sv

module fcnaux
   
    ! use mpivars
    use precision_definition, only : dp
    implicit none

    private 

    public :: compute_lnexppi_neutral, compute_lnexppi_acid, compute_lnexppi_base
    public :: compute_xpol_neutral, compute_xpol_chargeable
    public :: integral_betapi

contains


    ! Returns restricted exponent lnexppi of P(alpha) for neutral monomer.
    ! Definition
    ! ln(P(alpha))=-\beta \int dr dr' n(alpha,r') v(r',r) \pi(r) = 
    !             = (discreet) \sum_s \sum_i \beta deltav(j(alpha,s),i) \pi(i)) 
    !             := \sum_s exppi(j(alpha,s))
    ! with deltav(j,i):= deltavnucl(i-j)
    
    subroutine compute_lnexppi_neutral(xsol,deltavnucl,lnexppi)

        use globals, only : nsize
        use volume, only : indextocoord, coordtoindex, ipbc, nx, ny, nz

        real(dp), intent(in)    :: xsol(:)
        real(dp), intent(in)    :: deltavnucl(:,:,:)
        real(dp), intent(inout) :: lnexppi(:) 

        real(dp) :: lnexppitmp
        integer :: ix,iy,iz,i,j,k,ip,jp,kp ! integer cartesian coordinate lattice cell
        integer :: deltaix,deltaiy,deltaiz
        integer :: idx, idxnb ! index lattice cells


        do idx=1,nsize  ! loop over lattice cell

            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            lnexppitmp=0.0_dp

            do deltaix=0,1  ! loop over neigbouring cells
                i=ix+deltaix
                do deltaiy=0,1 
                    j=iy+deltaiy
                    do deltaiz=0,1
                        k=iz+deltaiz
                        ip=ipbc(i,nx) ! apply pbc 
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        idxnb=coordtoindex(ip,jp,kp) ! neighbour
                        lnexppitmp=lnexppitmp+deltavnucl(deltaix+1,deltaiy+1,deltaiz+1)*log(xsol(idxnb))
                    enddo
                enddo
            enddo 

            lnexppi(idx)=lnexppitmp

        enddo    

    end subroutine compute_lnexppi_neutral


    
    ! Returns restricted exponent lnexppi of P(alpha) related to  a multi component system.
    ! Definition 
    ! ln(P(alpha))=-\int dr \int dr' n(alpha;r') v(r',r) \beta \pi(r) 
    !               -\int dr n(alpha;r) \beta \psi(r)q_A^- -ln(f_A^-(r))
    !             = (discreet) \sum_s -\sum_i deltav(j(alpha,s),i) \beta \pi(i)) 
    !               \sum_s - beta \psi(j(alpha,s)) q_A^- -ln(f_A^_(j(alpha,s)))
    !            := \sum_s lnexppi(j(alpha,s))
    ! with deltav(j,i):= deltavnucl(i-j)
    ! Computes vector of exponents of Palpha for a given monomer type that is an acid
    ! gdis = fraction of acid in deprotonated state  A^- : state ==1 

    subroutine compute_lnexppi_acid(xsol,psi,gdis,deltavnucl,lnexppi)

        use globals, only : nsize
        use volume, only : indextocoord, coordtoindex, ipbc, nx, ny, nz

        real(dp), intent(in)    :: xsol(:)
        real(dp), intent(in)    :: psi(:) 
        real(dp), intent(in)    :: gdis(:)
        real(dp), intent(in)    :: deltavnucl(:,:,:)
        real(dp), intent(inout) :: lnexppi(:) 

        real(dp) :: lnexppitmp
        integer  :: ix,iy,iz,i,j,k,ip,jp,kp ! integer cartesian coordinate lattice cell
        integer  :: deltaix,deltaiy,deltaiz
        integer  :: idx, idxnb ! index of lattice cells

        do idx=1,nsize

            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            lnexppitmp=0.0_dp

            do deltaix=0,1
                i=ix+deltaix
                do deltaiy=0,1 
                    j=iy+deltaiy
                    do deltaiz=0,1
                        k=iz+deltaiz
                        ! apply pbc 
                        ip=ipbc(i,nx)
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        idxnb=coordtoindex(ip,jp,kp) ! index of neighbour
                        lnexppitmp=lnexppitmp+deltavnucl(deltaix+1,deltaiy+1,deltaiz+1)*log(xsol(idxnb))    
                    enddo
                enddo
            enddo 

            lnexppi(idx)=lnexppitmp + psi(idx) -log(gdis(idx))  

        enddo    

    end subroutine compute_lnexppi_acid

    ! Returns restricted exponent lnexppi of P(alpha) related to  a multi component system.
    ! Definition 
    ! ln(P(alpha))=-\int dr \int dr' n(alpha;r') v(r',r) \beta \pi(r) 
    !               -ln(f_B(r))
    !             = (discreet) \sum_s -\sum_i deltav(j(alpha,s),i) \beta \pi(i)) 
    !               \sum_s -ln(f_B(j(alpha,s)))
    !            := \sum_s lnexppi(j(alpha,s))
    ! with deltav(j,i):= deltavnucl(i-j) 
    ! Computes vector of exponents of Palpha for a given monomer type that is an base
    ! gdis = fraction of base in deprotonated, neutral state B : state ==2 
    ! State implicitely selected on call input in fcnnucl_ionbin_sv: 

    subroutine compute_lnexppi_base(xsol,gdis,deltavnucl,lnexppi)

        use globals, only : nsize
        use volume, only : indextocoord, coordtoindex, ipbc, nx, ny, nz

        real(dp), intent(in)    :: xsol(:)
        real(dp), intent(in)    :: gdis(:)
        real(dp), intent(in)    :: deltavnucl(:,:,:)
        real(dp), intent(inout) :: lnexppi(:) 

        real(dp) :: lnexppitmp
        integer  :: ix,iy,iz,i,j,k,ip,jp,kp 
        integer  :: deltaix,deltaiy,deltaiz
        integer  :: idx, idxnb 

        do idx=1,nsize 

            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            lnexppitmp=0.0_dp

            do deltaix=0,1
                i=ix+deltaix
                do deltaiy=0,1 
                    j=iy+deltaiy
                    do deltaiz=0,1
                        k=iz+deltaiz
                        ! apply pbc 
                        ip=ipbc(i,nx)
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        idxnb=coordtoindex(ip,jp,kp) ! neighbour
                        lnexppitmp=lnexppitmp+deltavnucl(deltaix+1,deltaiy+1,deltaiz+1)*log(xsol(idxnb))  
                    enddo
                enddo
            enddo  

            lnexppi(idx)=lnexppitmp -log(gdis(idx))  

        enddo    

    end subroutine compute_lnexppi_base


    subroutine compute_xpol_neutral(rho,deltavnucl,xpol)

        use globals, only : nsize 
        use volume, only : indextocoord, coordtoindex, ipbc, nx, ny, nz
        use parameters, only : vsol

        real(dp), intent(in)    :: rho(:)
        real(dp), intent(in)    :: deltavnucl(:,:,:)
        real(dp), intent(inout) :: xpol(:) 

        integer  :: ix,iy,iz,i,j,k,ip,jp,kp 
        integer  :: deltaix,deltaiy,deltaiz 
        integer  :: idx,idxnb 
        real(dp) :: xpoltmp

        do idx=1,nsize

            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            xpoltmp=0.0_dp

            do deltaix=0,1 
                i=ix-deltaix
                do deltaiy=0,1 
                    j=iy-deltaiy
                    do deltaiz=0,1
                        k=iz-deltaiz
                        ip=ipbc(i,nx) ! apply pbc 
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        idxnb=coordtoindex(ip,jp,kp) ! neighbour
                        xpoltmp=xpoltmp+deltavnucl(deltaix+1,deltaiy+1,deltaiz+1)*rho(idxnb)
                    enddo
                enddo
            enddo        
            
            xpol(idx)=xpol(idx)+xpoltmp*vsol
        
        enddo

    end subroutine compute_xpol_neutral


    subroutine compute_xpol_chargeable(rho,fdis,deltavnucl,xpol)

        use globals, only : nsize
        use volume, only : indextocoord, coordtoindex, ipbc, nx, ny, nz
        use parameters, only : vsol

        real(dp), intent(in)    :: rho(:)
        real(dp), intent(in)    :: fdis(:,:)
        real(dp), intent(in)    :: deltavnucl(:,:,:,:)
        real(dp), intent(inout) :: xpol(:) 

        integer  :: ix,iy,iz,i,j,k,ip,jp,kp
        integer  :: deltaix,deltaiy,deltaiz
        integer  :: idx,idxnb 
        real(dp) :: xpoltmp
        integer  :: state, nstates

        nstates=size(fdis,2)  ! number of states equal to dim=2 of fdis
        
        do idx=1,nsize

            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            xpoltmp=0.0_dp

            do deltaix=0,1 
                i=ix-deltaix
                do deltaiy=0,1 
                    j=iy-deltaiy
                    do deltaiz=0,1

                        k=iz-deltaiz
                        
                        ip=ipbc(i,nx) ! apply pbc 
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)

                        idxnb=coordtoindex(ip,jp,kp) ! neighbour
                       
                        do state=1,nstates
                            xpoltmp=xpoltmp+deltavnucl(deltaix+1,deltaiy+1,deltaiz+1,state)*rho(idxnb)*fdis(idxnb,state) 
                        enddo
                   
                    enddo
                enddo
            enddo  
            
            xpol(idx)=xpol(idx)+xpoltmp*vsol
        
        enddo

    end subroutine compute_xpol_chargeable

    

    ! auxilary routine used in fenergy module in subroutine FEchem_react_multi
    ! for systype=nucl_ionbin_sv
    ! integral 
    ! sumbetapi(r) =\beta \int dr' \pi(r') v(r,r')= 
    !              = (discreet)  -\sum_i \beta \pi(i) deltav(j,i) 
    !             := sumbetapi(j)
    ! with deltav(j,i):= deltavnucl(i-j)           
    
    subroutine integral_betapi(xsol,deltavnucl,sumbetapi)

        use globals, only : nsize
        use volume, only : indextocoord, coordtoindex, ipbc, nx, ny, nz

        real(dp), intent(in)    :: xsol(:)
        real(dp), intent(in)    :: deltavnucl(:,:,:)
        real(dp), intent(inout) :: sumbetapi(:) 

        real(dp) :: sumtmp
        integer  :: ix,iy,iz,i,j,k,ip,jp,kp ! integer cartesian coordinate lattice cell
        integer  :: deltaix,deltaiy,deltaiz
        integer  :: idx, idxnb ! index of lattice cells


        do idx=1,nsize 

            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            sumtmp=0.0_dp

            do deltaix=0,1
                i=ix+deltaix
                do deltaiy=0,1 
                    j=iy+deltaiy
                    do deltaiz=0,1
                        k=iz+deltaiz
                        ! apply pbc 
                        ip=ipbc(i,nx)
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        idxnb=coordtoindex(ip,jp,kp) ! neighbour
                        sumtmp=sumtmp-deltavnucl(deltaix+1,deltaiy+1,deltaiz+1)*log(xsol(idxnb))  
                    enddo
                enddo
            enddo  

            sumbetapi(idx)=sumtmp  

        enddo    

    end subroutine integral_betapi


end module fcnaux  

