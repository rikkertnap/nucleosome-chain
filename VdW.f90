!    .. module file for VdWcoeff
!    .. module file for computation of Van der Waals interaction 
!       Van der Waals energy :
!       EvdW = - 1/2 \sum{a,b} \esplison_{ab}\int dr \int dr' \rho_a(r) V_{a,b}(|r-r'|) \rho_b(r')
!       with       V_ab(r) =  (l_ab/r)^6  for l_ab<r < alpha l_ab 
!                  V_ab(r) =  0 otherwise !              
!    .. module comnputes V_ab using Monter Carlo simulation


module VdW
    
    use precision_definition
    implicit none

    real(dp), dimension(:,:,:,:,:), allocatable :: VdWcoeff    
    real(dp), dimension(:,:,:,:), allocatable :: rhopoltmp

    integer, parameter :: range = 2 
    integer, parameter :: MCsteps = 100000000
    
    integer, parameter :: VdW_err_allocation = 1
    integer, parameter :: VdW_err_vdwcoeff   = 2
    integer, parameter :: Vdw_err_systype    = 3
    integer, parameter :: VdW_err_inputfile  = 4
    integer, parameter :: Vdw_err_vdwcoeff_exist = 5 
    integer, parameter :: Vdw_err_vdwcoeff_not_exist = 6

    private
    public :: make_VdWcoeff,VdW_contribution_exp_diblock,VdW_contribution_exp_homo,VdW_contribution_exp
    public :: VdW_contribution_lnexp,VdW_energy_homo,VdW_energy_diblock,VdW_energy,make_VdWeps

contains

 

! function determines range VdW coeffcients 
! range = maxlayer=int(VdWcutoff*lseg/delta)+1
! +1 not neccsarry, for safety 

function set_range(lsegAA,VdWcutoff)result(range)

    use globals, only : nsegtypes
    use volume, only : delta
 
        
    real(dp) , intent(in)  :: lsegAA(:)
    real(dp) , intent(in)  :: VdWcutoff
    integer :: range
  

    real(dp) :: lseg
    logical :: flag
    integer :: rangetmp, t 

    flag=.true.

    lseg=lsegAA(1)
    range=int(VdWcutoff*lseg/delta)+1
    ! maxlayer = range 

    do t=2,nsegtypes
        lseg=lsegAA(t)
        rangetmp=int(VdWcutoff*lseg/delta)+1
        if(range<rangetmp) range=rangetmp
    enddo
        
end function

  
subroutine allocate_VdWcoeff(info)
    
    use globals, only : systype, nsegtypes

    integer :: ier
    logical :: alloc_fail

    integer,  intent(out), optional :: info
    
    if (present(info)) info = 0

    alloc_fail=.FALSE.    

    if (.not. allocated(VdWcoeff))  then 
        allocate(VdWcoeff(-range:range, -range:range, -range:range,nsegtypes,nsegtypes),stat=ier)            
        if( ier/=0 ) alloc_fail=.true.
    endif


    if(alloc_fail) then 
        print*,'Allocation warning: allocate_VdWcoeff failed'
        if(present(info)) info= VdW_err_allocation
    endif   

    
end subroutine allocate_VdWcoeff

  
subroutine allocate_auxdensity(info)

    use globals, only : systype,nsegtypes
    use volume, only : nx, ny, nz 
   
    integer,  intent(out), optional :: info

    integer :: ier
    
    if (present(info)) info = 0
    
    if (.not. allocated(rhopoltmp))  then 
        allocate(rhopoltmp(nx,ny,nz,nsegtypes),stat=ier)
               
        if( ier.ne.0 ) then
            print*, 'Allocation error in auxilary polymer density: stat =', ier
            if(present(info)) info= VdW_err_allocation                
            return
        endif
    else
       print*,'Allocation warning: auxilary polymer density already allocated'
    endif
    

end subroutine allocate_auxdensity
 
subroutine allocate_VdWeps
    
    use globals, only : nsegtypes
    use parameters, only : VdWeps, VdWepsin

    integer :: ier
   
    if (.not. allocated(VdWeps))  then 
        allocate(VdWeps(nsegtypes,nsegtypes),stat=ier)
    endif        

    if(ier/=0) then 
        print*,'Allocation error: allocate_VdWeps failed'
        stop
    endif    

    if (.not. allocated(VdWepsin))  then 
        allocate(VdWepsin(nsegtypes,nsegtypes),stat=ier)
    endif        

    if(ier/=0) then 
        print*,'Allocation error: allocate_VdWeps failed'
        stop
    endif    
       
end subroutine allocate_VdWeps


subroutine make_VdWcoeff(info)

    use mpivars, only : rank
    use globals, only : systype, nsegtypes
    use volume, only : geometry
    use parameters, only : lsegPAA,lsegPAMPS, lsegAA
    use myutils

    integer,  intent(out), optional :: info
    
    real(dp) :: lsegAB, lseg
    character(len=lenText) :: text
    integer :: info_allocate_VdW,info_allocate_dens,info_VdWeps
    integer :: s,t

    if (present(info)) info = 0

    !range=set_range(lsegAA,VdWcutoff)

    call allocate_VdWcoeff(info_allocate_VdW)
    call allocate_auxdensity(info_allocate_dens)
   
    do t=1,nsegtypes
        do s=1,nsegtypes
            lseg=(lsegAA(t)+lsegAA(s))/2.0_dp
       
            call read_VdWcoeff(lseg,VdWcoeff(:,:,:,t,s),info)

            if(info/=0) then  
                call MC_VdWcoeff(lseg, VdWcoeff(:,:,:,t,s))
                info=0 ! reset
                if(rank==0) call write_VdWcoeff(lseg,VdWcoeff(:,:,:,t,s),info)
            endif    
        enddo
    enddo    

    call make_VdWeps(info_VdWeps)   
    if (info_VdWeps == VdW_err_inputfile) then
        text="make_VdWcoeff: make_VdWeps failed"
        call print_to_log(LogUnit,text)
        if (present(info)) info = VdW_err_inputfile
        return
    endif

end subroutine




! Monte Carlo simulation to deterimine VdW coefficients


subroutine MC_VdWcoeff(lseg,VdWcoeff)

    use mpivars
    use mathconst
    use volume, only : xt, yt, ut, vt
    use volume, only : delta
    use random 
    use myutils

    real(dp), intent(in)   :: lseg ! size  segment 
    real(dp), intent(out) :: VdWcoeff(-range:range, -range:range, -range:range)
    
    integer :: ix, iy , iz
    real(dp) :: x,y,z, radius, u, v
    real(dp) :: rn, maxradius, limittimesdelta
    integer :: limit
    !parameter (limit = range+1) 
    !real(dp) :: matriz(-limit:limit, -limit:limit, -limit:limit) ! matrix for chi
    
    real(dp), allocatable, dimension(:,:,:) :: matriz 
    
    integer :: i,status
    real(dp) :: sum
    character(len=lenText) :: text, rstr, istr


    if(rank==0) then 
        text="VdWcoefficient calculation"
        call print_to_log(LogUnit,text)
    endif    

    limit = range+1
    limittimesdelta = limit* delta
    maxradius = (limit*1.0_dp/2.0_dp) * delta

    allocate(matriz(-limit:limit, -limit:limit, -limit:limit),stat=status)! matrix for chi
    if(status/=0) then
        print*,"MC_VdWcoeff: allocation error"
    endif

    do ix = -limit, limit
        do iy = -limit, limit
            do iz = -limit, limit
                matriz(ix, iy, iz) = 0.0_dp
            enddo
        enddo
    enddo

    seed = 1010
      

    do i = 1, MCsteps

        !x = 3.0*(rands(seed)-0.5_dp)*delta ! random number between -1.5 * delta and 1.5 * delta
        !y = 3.0*(rands(seed)-0.5_dp)*delta 
        !z = 3.0*(rands(seed)-0.5_dp)*delta 

        x = limittimesdelta*(rands(seed)-0.5_dp)! random number between - limit/2 * delta and limit/2 * delta
        y = limittimesdelta*(rands(seed)-0.5_dp)
        z = limittimesdelta*(rands(seed)-0.5_dp) 

        u = ut(x, y)
        v = vt(x, y)

        radius = sqrt(x**2 + y**2 + z**2) ! real space
 
        !if(radius<=(1.5*delta)) then  ! It is inside the cut-off distance 1.5 delta
        if(radius<=maxradius) then  ! It is inside the cut-off distance 

            ix = anint(v/delta)   ! 
            iy = anint(u/delta) 
            iz = anint(z/delta)

            if(radius>=lseg) then     ! It is not in side sphere of the segment 
            
                matriz(ix, iy, iz) = matriz(ix, iy, iz) + (lseg/radius)**6
            else
                matriz(ix, iy, iz) = matriz(ix, iy, iz) + 1.0_dp
            endif    
        endif
                
    enddo

    ! do ix = -2, 2
    !     do iy = -2, 2
    !         do iz = -2, 2
    !
    !            VdWcoeff(ix, iy, iz) = matriz(ix, iy, iz)/MCsteps*((3.0_dp*delta)**3)
    !            sum = sum +  matriz(ix, iy, iz)/MCsteps*((3.0_dp*delta)**3)
    !    
    !        enddo
    !    enddo
    !enddo

    sum = 0.0_dp

    do ix = -range, range
        do iy = -range, range
            do iz = -range, range
                VdWcoeff(ix, iy, iz) = matriz(ix, iy, iz)/MCsteps*((limittimesdelta)**3)
                sum = sum +  matriz(ix, iy, iz)/MCsteps*((limittimesdelta)**3)
        
            enddo
        enddo
    enddo


    if(rank.eq.0)then 
        write(istr,'(I3)')range
        text="VdWcoefficient calculation : range  = "//trim(adjustl(istr))
        write(rstr,'(F5.3)') sum
        text="VdWcoefficient calculation : Sum (2xrange+1)x(2xrange+1) = "//trim(adjustl(rstr))
        call print_to_log(LogUnit,text)
    endif    

    sum = 0.0_dp

    do ix = -limit, limit
        do iy = -limit, limit
            do iz = -limit, limit
                sum = sum +  matriz(ix, iy, iz)/MCsteps*((3.0_dp*delta)**3)
            enddo
        enddo
    enddo

    if(rank.eq.0) then 
        write(rstr,'(F5.3)') sum
        text="VdWcoefficient calculation : Sum Total = "//trim(adjustl(rstr))
        call print_to_log(LogUnit,text)
    endif    
 
end subroutine



! this compute contribution to Palpha

subroutine VdW_contribution_exp_homo(rhopol,exppi)

    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate
    use parameters, only : VdWepsAA

!    real(dp), intent(in) :: VdWcoeff(:,:,:)
    real(dp), intent(in) :: rhopol(:)
    real(dp), intent(inout) :: exppi(:)

    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: protemp
    integer,  parameter  :: A=1
    
    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopoltmp(ix,iy,iz,A) = rhopol(id)   
    enddo
    
    !VdWstr = VdWeps

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                protemp = 0.0_dp
                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            jz = ipbc(jz,nz)
                            protemp=protemp + VdWcoeff(ax,ay,az,A,A)*rhopoltmp(jx, jy, jz,A)
                                
                        enddo
                    enddo
                enddo
                
                call linearIndexFromCoordinate(ix,iy,iz,id)
                exppi(id) = exppi(id)*exp(VdWepsAA*protemp)

            enddo
        enddo
    enddo

end subroutine



! this compute contribution to Palpha

subroutine VdW_contribution_exp_diblock(rhopolA,rhopolB,exppiA,exppiB)


    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate
    use parameters, only : VdWepsAA, VdWepsAB, VdWepsBB

    real(dp), intent(in) :: rhopolA(:),rhopolB(:)
    real(dp), intent(inout) :: exppiA(:),exppiB(:)

    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: protempAA,protempAB,protempBA,protempBB
    integer,  parameter  :: A=1, B=2

    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopoltmp(ix,iy,iz,A) = rhopolA(id)  
        rhopoltmp(ix,iy,iz,B) = rhopolB(id)    
    enddo
    
    !VdWstr = VdWeps

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                protempAA = 0.0_dp
                protempAB = 0.0_dp
                protempBA = 0.0_dp
                protempBB = 0.0_dp

                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            jz = ipbc(jz,nz)

                            protempAA=protempAA + VdWcoeff(ax,ay,az,A,A)*rhopoltmp(jx, jy, jz,A)
                            protempAB=protempAB + VdWcoeff(ax,ay,az,A,B)*rhopoltmp(jx, jy, jz,B)
                            protempBA=protempBA + VdWcoeff(ax,ay,az,A,B)*rhopoltmp(jx, jy, jz,A)
                            protempBB=protempBB + VdWcoeff(ax,ay,az,B,B)*rhopoltmp(jx, jy, jz,B)
                            
                        enddo
                    enddo
                enddo
                
                call linearIndexFromCoordinate(ix,iy,iz,id)

                exppiA(id) = exppiA(id)*exp(VdWepsAA*protempAA+VdWepsAB*protempAB)
                exppiB(id) = exppiB(id)*exp(VdWepsAB*protempBA+VdWepsBB*protempAB)
            enddo
        enddo
    enddo

end subroutine



! this compute contribution to Palpha

subroutine VdW_contribution_exp(rhopol,exppi,segtype)

    use globals, only : nsize, nsegtypes
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate
    use parameters, only : VdWeps

    real(dp), intent(in) :: rhopol(:,:)
    real(dp), intent(inout) :: exppi(:)
    integer,  intent(in) :: segtype

    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: protemp
    integer :: t
    
    !  assign rhopoltmp 
    
    do t=1, nsegtypes
        do id=1,nsize
            call coordinateFromLinearIndex(id, ix, iy, iz)
            rhopoltmp(ix,iy,iz,t) = rhopol(id,t)     
        enddo
    enddo

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz
                protemp= 0.0_dp
                do ax = -range,range 
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                    do ay = -range,range
                        jy = iy+ay
                        jy = ipbc(jy,ny) 
                        do az = -range,range
                            jz = iz+az
                            jz = ipbc(jz,nz)
                            do t=1,nsegtypes
                                protemp = protemp + VdWeps(segtype,t) * VdWcoeff(ax,ay,az,segtype,t) &
                                    * rhopoltmp(jx,jy,jz,t)
                            enddo
                        enddo
                    enddo
                enddo
                call linearIndexFromCoordinate(ix,iy,iz,id)
                exppi(id) = exppi(id)*exp(protemp)
            enddo
        enddo
    enddo

end subroutine VdW_contribution_exp

! this compute contribution to Palpha

subroutine VdW_contribution_lnexp(rhopol,lnexppi,segtype)

    use globals, only : nsize, nsegtypes
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate
    use parameters, only : VdWeps

    real(dp), intent(in) :: rhopol(:,:)
    real(dp), intent(inout) :: lnexppi(:)
    integer,  intent(in) :: segtype

    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: protemp
    integer :: t
    
    !  assign rhopoltmp 
    
    do t=1, nsegtypes
        do id=1,nsize
            call coordinateFromLinearIndex(id, ix, iy, iz)
            rhopoltmp(ix,iy,iz,t) = rhopol(id,t)     
        enddo
    enddo

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz
                protemp= 0.0_dp
                do ax = -range,range 
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                    do ay = -range,range
                        jy = iy+ay
                        jy = ipbc(jy,ny) 
                        do az = -range,range
                            jz = iz+az
                            jz = ipbc(jz,nz)
                            do t=1,nsegtypes
                                protemp = protemp + VdWeps(segtype,t) * VdWcoeff(ax,ay,az,segtype,t) &
                                    * rhopoltmp(jx,jy,jz,t)
                            enddo
                        enddo
                    enddo
                enddo
                call linearIndexFromCoordinate(ix,iy,iz,id)
                lnexppi(id) = lnexppi(id) + protemp
            enddo
        enddo
    enddo

end subroutine VdW_contribution_lnexp


function VdW_energy_homo(rhopol)result(EVdW)

    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate,volcell
    use parameters, only : VdWepsAA


    real(dp), intent(in) :: rhopol(:)

    real(dp) :: EVdW
    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: Etemp
    integer, parameter :: A=1
    
    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        print*,id,ix,iy,iz 
        rhopoltmp(ix,iy,iz,A) = rhopol(id) 
    enddo
    
    Etemp =0.0_dp

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            jz = ipbc(jz,nz)

                            Etemp=Etemp + VdWcoeff(ax,ay,az,A,A)*rhopoltmp(jx, jy, jz,A)*rhopoltmp(ix, iy, iz,A)
                                
                        enddo
                    enddo
                enddo

            enddo
        enddo
    enddo

    EVdW = -VdWepsAA*Etemp*volcell/2.0_dp

end function



function VdW_energy_diblock(rhopolA,rhopolB)result(EVdW)

    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate,volcell
    use parameters, only : VdWepsAA, VdWepsBB,VdWepsAB


    real(dp), intent(in) :: rhopolA(:),rhopolB(:)

    real(dp) :: EVdW
    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: EtempAA,EtempAB,EtempBB
    integer,  parameter  :: A=1, B=2

    
    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopoltmp(ix,iy,iz,A) = rhopolA(id)
        rhopoltmp(ix,iy,iz,B) = rhopolB(id)   
    enddo
    
    EtempAA =0.0_dp
    EtempAB =0.0_dp
    EtempBB =0.0_dp

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            jz =ipbc(jz,nz)

                            EtempAA=EtempAA + VdWcoeff(ax,ay,az,A,A)*rhopoltmp(jx, jy, jz,A)*rhopoltmp(ix, iy, iz,A)
                            EtempAB=EtempAB + VdWcoeff(ax,ay,az,A,B)*rhopoltmp(jx, jy, jz,A)*rhopoltmp(ix, iy, iz,B)
                            EtempBB=EtempBB + VdWcoeff(ax,ay,az,B,B)*rhopoltmp(jx, jy, jz,B)*rhopoltmp(ix, iy, iz,B)
                                
                        enddo
                    enddo
                enddo

            enddo
        enddo
    enddo

    EVdW = -(VdWepsAA*EtempAA+2.0_dp*VdWepsAB*EtempAB+VdWepsBB*EtempBB)*volcell/2.0_dp

end function



function VdW_energy(rhopol)result(EVdW)

    use globals, only : nsize, nsegtypes
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate,volcell
    use parameters, only : VdWeps

    real(dp), intent(in) :: rhopol(:,:)

    real(dp) :: EVdW
    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: Etemp
    integer :: s, t
    
    !  assign rhopoltmp

    do t=1, nsegtypes
        do id=1,nsize
            call coordinateFromLinearIndex(id, ix, iy, iz)
            rhopoltmp(ix,iy,iz,t) = rhopol(id,t)     
        enddo
    enddo

    EVdW =0.0_dp

    do s=1,nsegtypes
        do t=1,nsegtypes
    
            Etemp =0.0_dp
   
            do ix=1,nx
                do iy=1,ny
                    do iz=1,nz
                        do ax = -range,range 
                
                            jx = ix+ax
                            jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                        
                            do ay = -range,range
                                jy = iy+ay
                                jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                        
                                do az = -range,range 
                                    jz = iz+az
                                    jz = ipbc(jz,nz)

                                    Etemp=Etemp + VdWcoeff(ax,ay,az,s,t)*rhopoltmp(jx, jy, jz,s)*rhopoltmp(ix, iy, iz,t)
                                        
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo

            EVdW = EVdW-VdWeps(s,t)*Etemp*volcell/2.0_dp
        enddo
    enddo            

end function


subroutine read_VdWeps(info)
    
    use mpivars
    use globals, only : nsegtypes , runtype
    use parameters, only : VdWeps, VdWepsin 
    use myutils, only : newunit

    integer,  intent(out), optional :: info
     !     .. local variables 
    character(len=9) :: fname
    integer :: ios, un_input  ! un = unit number
    character(80) :: str
    integer :: s,t, line
   

    ! .. reading in of variables from file 

    fname="VdWeps.in"
   
    open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening input file VdWeps.in : iostat =', ios
        if (present(info)) info = VdW_err_inputfile
        return
    endif    
    
    
    ios=0
    line=0
    do while (line<(nsegtypes**2).and.ios==0)
        line=line+1
        read(un_input,*,iostat=ios)t,s,VdWeps(t,s)
    enddo

    
    if(line/=(nsegtypes)**2) then 
        str="reached end of VdWeps.in before all elements read"
        print*,str
        str="read file "//trim(adjustl(fname))//" failed"
        print*,str
        call MPI_FINALIZE(ierr)
        stop
    endif

    if(ios >0 ) then
        print*, 'Error parsing VdWeps.in : iostat =', ios
        if (present(info)) info = VdW_err_inputfile
        return
    endif
    
    close(un_input)   


end subroutine read_VdWeps




subroutine make_VdWeps(info)

    use parameters, only : set_VdWepsAAandBB, set_VdWepsin
    integer,  intent(out), optional :: info
    

    call allocate_VdWeps()
    call read_VdWeps(info)
    ! special assign 
    call set_VdWepsAAandBB
    call set_VdWepsin

end subroutine


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

!  return values 
!  info = 0                      : write succesfull 
!  info = VdW_err_vdwcoeff       : error opening file to write 
!  info = VdW_err_vdwcoeff_exist : write file existes
!  info otherwise                : write error 

subroutine write_VdWcoeff(lseg,VdWcoeff,info)
    
    use globals, only : nsize 
    use parameters, only : VdWcutoff
    use volume, only : delta, geometry 
    use myutils, only : newunit

    real(dp), intent(in)    :: lseg ! size  segment 
    real(dp), intent(inout) :: VdWcoeff(-range:range, -range:range, -range:range)
    integer,  intent(out), optional :: info

    character(len=40) :: fname
    integer :: ios, ier, un_VdW  ! un = unit number
    character(len=10)  :: rstr
    logical :: exist
    integer :: ix,iy,iz


    if (present(info)) info = 0

    ! open new file 
    write(rstr,'(F5.3)')lseg
    fname="vdwcoeff_lseg"//trim(adjustl(rstr))
    fname=trim(fname)//"_"//trim(adjustl(geometry))//".dat"
    fname=trim(fname)
    
    inquire(file=fname,exist=exist)

    !print*,"Exist=",exist
    if(.not.exist) then
        open(unit=newunit(un_VdW),file=fname,iostat=ios,status='new')
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            if (present(info)) info = VdW_err_vdwcoeff
            return
        endif
    else
        print*,'VdWcoeff file allready exist'
        if (present(info)) info = VdW_err_vdwcoeff_exist
        return
    endif

!"(A9,I1,A5,ES25.16)"

    ! write preamble 
    write(un_VdW,'(A15,I2)',  iostat=ios)'range          ',range
    write(un_VdW,'(A15,A15)', iostat=ios)'geometry       ',geometry
    write(un_VdW,'(A15,F5.3)',iostat=ios)'lseg           ',lseg
    write(un_VdW,'(A15,F5.3)',iostat=ios)'delta          ',delta
    write(un_VdW,'(A15,F5.3)',iostat=ios)'cutoff         ',VdWcutoff
    write(un_VdW,'(A15)',     iostat=ios)'###            ' ! end of preamble

    ! write value VdWcoefficient
    do ix = -range, range
        do iy = -range, range
            do iz = -range, range
                 write(un_VdW,*,iostat=ios)VdWcoeff(ix, iy, iz) 
                 !print*,ix,iy,iz,VdWcoeff(ix, iy, iz) 
            enddo
        enddo
    enddo

    ! close file
    close(un_VdW)

    ! succesfull write 
    if(ios==0 ) then 
        if (present(info)) info = 0
    endif    

end subroutine write_VdWcoeff


!  return values 
!  info = 0                          : read succesfull 
!  info = VdW_err_vdwcoeff           : error opening file to read 
!  info = VdW_err_vdwcoeff_not_exist : read file does not existes
!  info otherwise                    : read error 

subroutine read_VdWcoeff(lseg,VdWcoeff,info)
    
    use globals, only : nsize 
    use parameters, only : VdWcutoff
    use volume, only : delta, geometry 
    use myutils, only : newunit


    real(dp), intent(in)    :: lseg ! size  segment 
    real(dp), intent(inout) :: VdWcoeff(-range:range, -range:range, -range:range)
    integer,  intent(out), optional :: info
    
    character(len=40) :: fname
    integer :: ios, ier, un_input  ! un = unit number
    character(len=10)  :: rstr
    
    !     .. local variables 
    integer  :: range_file
    real(dp) :: lseg_file
    real(dp) :: delta_file
    real(dp) :: VdWcutoff_file
    character(len=15) :: geometry_file
    integer  :: ix,iy,iz
    character(len=100) :: buffer, label
    integer :: pos, line
    logical :: not_end_preamble, exist
    real(dp) :: val
   
    if (present(info)) info = 0

    ! .. reading in of variables from file 
    write(rstr,'(F5.3)')lseg
    fname="vdwcoeff_lseg"//trim(adjustl(rstr))
    fname=trim(fname)//"_"//trim(adjustl(geometry))//".dat"
    fname=trim(fname)
   
    inquire(file=fname,exist=exist)

    if(exist) then 
        open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
        if(ios >0 ) then
            print*, 'Error opening input file ',fname,': iostat =', ios
            if (present(info)) info = VdW_err_inputfile
            return
        endif
    else
        if (present(info)) info = VdW_err_vdwcoeff_not_exist
        return
    endif        
    
    ios=0 
    line = 0
    not_end_preamble=.true.

    ! ios<0 : if an end of record condition is encountered or if an end of file condition was detected.  
    ! ios>0 : if an error occured 
    ! ios=0 : otherwise.

    do while (ios == 0 .and. not_end_preamble) ! scans preamble 

        read(un_input, '(A)', iostat=ios) buffer

        if (ios == 0) then
        
            line = line + 1

            !  Split label and data based on first occurence of a whitespace
            pos = scan(buffer,'     ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            select case (label) !list-directed The  charackter variable is treated as an 'internal file'
            case ('range')
                read(buffer,*,iostat=ios) range_file
            case ('geometry')
                read(buffer,*,iostat=ios) geometry_file
            case ('lseg')
                read(buffer,*,iostat=ios) lseg_file
            case ('delta')
                read(buffer,*,iostat=ios) delta_file
            case ('cutoff')
                read(buffer,*,iostat=ios) VdWcutoff_file
            case ('###') ! end of preamble
                not_end_preamble=.false. 
            case default
            end select
        endif

    enddo     
    
    if(ios >0 ) then
        print*, 'Error parsing file : iostat =', ios
        close(un_input)
        if (present(info)) info = VdW_err_vdwcoeff
        return
    endif

    if(delta_file /= delta) then 
        print*,'delta in VdWcoeff file not equal delta' 
        print*,'delta_file = ',delta_file, ' delta = ',delta
        close(un_input)
        if (present(info)) info = VdW_err_vdwcoeff
        return
    endif

    if(VdWcutoff_file /= VdWcutoff) then 
        print*,'cuttoff in VdWcoeff file not equal imposed cutoff'
        print*,'cutoff_file = ',VdWcutoff_file, ' cutoff = ',VdWcutoff
        close(un_input)
        if (present(info)) info = VdW_err_vdwcoeff
        return
    endif
    
    if(range_file /= range) then
        print*,'range in VdWcoeff file not large enough'
        print*,'range_file=',range_file,'range =',range
        close(un_input)
        if (present(info)) info = VdW_err_vdwcoeff
        return
    endif
    
    if(geometry_file /= geometry) then 
        print*,'geometry in VdWcoeff file not equal delta'
        close(un_input)
        if (present(info)) info = VdW_err_vdwcoeff
        return
    endif        
    
    ! read file

    do ix=-range,range
        do iy=-range,range
            do iz=-range,range
                read(un_input,*)VdWcoeff(ix,iy,iz)
            enddo
        enddo
    enddo    

    close(un_input)
    
    
end subroutine read_VdWcoeff




end module VdW

      
