!     makes volume elements 

module volume      

    use precision_definition
    use mpivars

    implicit none
  
    !     .. variables
   
    real(dp) :: delta               ! delta  spacing of lattice site in x-, y- and z-direction  
    integer  :: nz                  ! nz number of lattice sites in z-direction 
    integer  :: nx                  ! nx number of lattice sites in x-direction 
    integer  :: ny                  ! ny number of lattice sites in y-direction 
    integer  :: nsurf               ! nsurf number of lattice site  at z=0 ore z=nz*delta
    real(dp) :: volcell             ! volcell=delta**3
    real(dp) :: areacell            ! areacell=delta**2
    real(dp) :: areasurf            ! area surfaces spanned in x- and y- direction
    real(dp) :: gamma               ! angle between oblique basis vectors u and v: cubic = gamma=90=pi/2 hexagonl gamma=60=2pi/3                          
    real(dp) :: beta                ! related beta = (pi/2- gamma)/2, angle between basis vector u and x and v and y
    real(dp) :: cos_two_beta        ! sqrt(cos(beta)**2 - sin(beta)**2)=cos(2beta) scales u and v coordinates 
    real(dp) :: sin_two_beta        ! sin(2beta)  
    character(len=11) :: geometry

    ! variable for grafting position
    
    integer :: ngr                  ! total number of graft points    
    integer :: ngr_node             ! number of grafted areas assigned to an individual node
    integer :: ngr_freq             ! frequence spacing in terms of delta 
    integer :: ngrx                 ! total number of graft points along x-axis 
    integer :: ngry                 ! total number of graft points along y-axis                                                       
    logical :: isRandom_pos_graft   ! true random position, false regual pattern
    integer :: seed_graft           ! seed for graft points 
    real(dp) :: scale_ran_step      ! scale random postition,minimum value =1.0 =half ngr_freq
    integer :: sgraft               ! item number of graft point to which loop is attached 
    integer :: nset_per_graft       ! number of confomation set to read in graft point
    real(dp), dimension(:,:), allocatable :: position_graft
   

    ! hash table
    integer, dimension(:,:,:), allocatable :: coordtoindex 
    integer, dimension(:,:),  allocatable  :: indextocoord
    integer, dimension(:,:), allocatable   :: indexneighbor
    !    integer, dimension(:,:), allocatable   :: inverse_indexneighbor
    integer, dimension(:,:), allocatable   :: inverse_indexneighbor_phos

    private
    
    public :: delta,nx,ny,nz,volcell,areacell,geometry, nsurf
    public :: gamma,cos_two_beta, sin_two_beta
    public :: indexneighbor !, inverse_indexneighbor, 
    public :: inverse_indexneighbor_phos
    public :: coordtoindex, indextocoord
    public :: coordinateFromLinearIndex, linearIndexFromCoordinate
    public :: xt, yt, ut, vt, ipbc
    public :: make_geometry
    public :: allocate_index_neighbors_phos, make_table_index_neighbors_phos
    public :: position_graft, isRandom_pos_graft , ngr, ngr_freq, nset_per_graft, seed_graft
    public :: scale_ran_step

    
contains

    subroutine init_lattice

        use globals, only : nsize, systype
        use mathconst
        use chains, only : distphoscutoff ,maxneigh

        implicit none 

        integer :: rangecutoff


        if(geometry=="cubic") then 
            gamma = pi/2.0_dp
        else
            gamma = gamma * pi /180.0_dp ! convert from degree to radians
        endif        

        beta = (pi/2.0_dp - gamma) / 2.0_dp       ! used by ut, vt, xt and yt functions
        cos_two_beta=cos(2.0_dp*beta) ! sqrt(cos(beta)**2 - sin(beta)**2)  ! scaling of u and v coordinates
        sin_two_beta=sin(2.0_dp*beta)

        ! cubic lattice  or prism surface in x-y direction at z=0 and z=nz  
        ! nz=nzmax
        nsize = nx*ny*nz                       ! total number of cells or layers
        volcell = delta*delta*delta*1.0_dp     ! volume of one latice volume 
        areacell = delta*delta
        nsurf = nx*ny  
        areasurf = nsurf*delta*delta
      
        ! setup graft points

        call make_graftpoints

        ! use of hash table 
        call allocate_hashtable(nx,ny,nz)
        call make_hashtable()

        ! call write_hashtable()

        ! this should go somewhere else
 
        if(systype=="nucl_ionbin_Mg".or.systype=="nucl_ionbin_MgA") then
            rangecutoff=int(distphoscutoff/delta)+2 
            maxneigh = (2*rangecutoff+1)**3
            !call allocate_index_neighbors(maxneigh)
            !call make_table_index_neighbors(rangecutoff)
        endif    

    end subroutine init_lattice
         

    subroutine  make_geometry()
    
        use globals
        use myutils

        implicit none
    
        character(len=lenText) :: text, str 

        call init_lattice

        write(str,'(A20)')geometry
        text="geometry="//trim(adjustl(str)) 
        call print_to_log(LogUnit,text) 


    end subroutine make_geometry

    subroutine linearIndexFromCoordinate(x,y,z,idx)
     
      implicit none 
      
        integer, intent(in)  :: x,y,z
        integer, intent(out) :: idx
      
        integer :: a,b,c,d

        a = 1
        b = nx 
        c = nx * ny 
        d = 1-a-b-c
        idx = a*x + b*y + c*z + d

    end subroutine linearIndexFromCoordinate


    subroutine coordinateFromLinearIndex(idx, x,y, z)

        implicit none

        integer, intent(out)  :: x,y,z
        integer, intent(in)   :: idx
        integer :: idxtmp
   
        idxtmp=idx
        x =  mod(idxtmp-1,nx)+1
        idxtmp =int((idxtmp-1)/nx)+1
        y = mod(idxtmp-1,ny)+1
        idxtmp = int((idxtmp-1)/ny)+1
        z = idxtmp
    
    end subroutine coordinateFromLinearIndex

     
    function mirror_index(idx,nz) result(idx_mirror)

        implicit none 
      
        integer, intent(in)  :: idx, nz
        integer :: idx_mirror

        integer :: x,y,z

        call coordinateFromLinearIndex(idx,x,y,z)
        call LinearIndexFromCoordinate(x,y,nz+1-z,idx_mirror)

    end function mirror_index

    ! coordinate transformation from cartesian x,y to prism or oblique u,v coordinates
    ! u = s(cos(gamma)*x -sin(gamma)*y); s=1/srqt(cos(gamma)^2-sin(gamma)^2)
    ! v = s(-sin(gamma)*x +sin(gamma)*y); 
    ! metric ds^2= s^2 du^2 + s^2 dv^2 + ss^2 sin(2 gamma) du dv +dz^2

    function ut(x, y) result(ut_val)
       
        real(dp), intent(in) :: x, y
        real(dp) :: ut_val
        
        ut_val = cos(beta)*x - sin(beta)*y
        ut_val = ut_val / cos_two_beta
    
    end function ut


    ! coordinate transformation from cartesian x,y to prism u,v  coordinates

    function vt(x, y) result(vt_val)
    
        real(dp), intent(in) :: x, y
        real(dp) :: vt_val
       
        vt_val = -sin(beta)*x + cos(beta)*y
        vt_val = vt_val / cos_two_beta
    
    end function vt



    ! inverse coordinate transformation from prism u,v to  cartesian x,y  coordinates

    function xt(u, v) result(xt_val)
        
        real(dp), intent(in) :: u, v
        real(dp) :: xt_val
    
        xt_val = cos(beta)*u + sin(beta)*v
        xt_val = xt_val / cos_two_beta
    
    end function xt

    ! inverse coordinate transformation from prism u,v to  cartesian x,y  coordinates
    
    function yt(u, v) result(yt_val)

        real(dp), intent(in) :: u, v
        real(dp) :: yt_val
    
        yt_val = sin(beta)*u + cos(beta)*v
        yt_val = yt_val / cos_two_beta
    
    end function yt
 
    function isEven(number) result(val) 
        integer, intent(in) :: number
        logical :: val

        if ( mod(number, 2) ==0) then
            val= .True. 
        else
            val= .False.
        endif 
   
    end function

    subroutine allocate_hashtable(nx,ny,nz)

        use globals, only : nsize
        integer, intent(in) :: nx,ny,nz

        allocate(coordtoindex(nx,ny,nz))
        allocate(indextocoord(nsize,3))
        
    end subroutine allocate_hashtable

    subroutine make_hashtable

        use globals, only : nsize

        integer :: idx, ix, iy, iz ,idxtmp
            
        do idx=1,nsize
            
            call coordinateFromLinearIndex(idx, ix, iy,iz)
            coordtoindex(ix,iy,iz)=idx
            
            indextocoord(idx,1)=ix
            indextocoord(idx,2)=iy
            indextocoord(idx,3)=iz

            ! ..check 
            call LinearIndexFromCoordinate(ix,iy,iz,idxtmp)
            if(idx/=idxtmp)then
                print*,"Error in hashtable :idx=",idx," not equal to idxtmp",idxtmp
            endif
            
        enddo
        
    end subroutine make_hashtable

    subroutine write_hashtable

        use globals, only : nsize
        use myutils, only : lenText, newunit, error_handler

        integer :: idx, ix, iy, iz, idxtmp
        integer :: ios, un, info
        character(len=20) :: fname, istr
        character(len=lenText) :: text   
        character(len=100) :: io_msg  

        write(istr,'(I4)')rank
        fname='hashtable.'//trim(adjustl(istr))//'.dat'
        open(unit=newunit(un),file=fname,iostat=ios,status='new',iomsg=io_msg)
        if(ios >0 ) then
            write(istr,'(I4)')ios
            text='Error opening hashtable file : iostat ='//istr
            print*,text
            print*, 'Error message : ',trim(io_msg)
            info=1
            call error_handler(info,text)
        endif

        write(un,*)nsize

        do idx=1,nsize
            
            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            write(un,*)idx,ix,iy,iz
        enddo

        close(un)

    end subroutine write_hashtable

    ! compute periodic boundary condition in integer units 
    ! real(dp) version  : pbbc is located in chaingenerator  

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
         
   
    subroutine allocate_index_neighbors(maxneigh)

!        use globals, only : nsize

        integer, intent(in) :: maxneigh

!        allocate(indexneighbor(nsize,maxneigh))
!        allocate(inverse_indexneighbor(nsize,nsize))
        
    end subroutine allocate_index_neighbors

    
    ! Calculate indexneighbor(idx,k) 
    ! indexneighbor(idx,k) = return lattice cell index of neighbor number k of lattice cell index idx
    ! range : 1<=idx<=nsize 
    !       : 1<=k<=(2*rangecutoff+1)**3 
    !       : k=1 is same lattice cell         

    subroutine make_table_index_neighbors(rangecutoff)

        use globals, only : nsize

        integer, intent(in) :: rangecutoff
        
        integer :: idx, i, j, k, ix, iy, iz, deltaix,deltaiy, deltaiz
        integer :: neighbornumber, idxneigh
        integer :: ip, jp, kp


        do idx=1,nsize

            ! use hash table to get coordiantes 
            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            neighbornumber=0 ! 

            do deltaix=-rangecutoff,rangecutoff
                i=ix+deltaix
                do deltaiy=-rangecutoff,rangecutoff
                    j=iy+deltaiy
                    do deltaiz=-rangecutoff,rangecutoff
                        k=iz+deltaiz
                        ! apply pbc 
                        ip=ipbc(i,nx)
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        neighbornumber=neighbornumber+1
                        idxneigh=coordtoindex(ip,jp,kp) ! index of neighbour
                        indexneighbor(idx,neighbornumber) = idxneigh
                    !    inverse_indexneighbor(idx,idxneigh) = neighbornumber
                    enddo
                enddo
            enddo

        enddo    

    end subroutine make_table_index_neighbors


    ! Stores in a memory efficient way inverse_index_neighbors
    !  
    ! inverse_indexneighbor(ind,k) evaluated for all lattice positions  
    ! inverse_indexneighbor_phos only evaluate for those position that contain a phosphate   
    ! inverse_indexneighbor_phos(ind,k): argument ind is different. 
    ! Here "ind" is ind th element of array that enumarated 
    ! all different phosphate postions 
    ! The assocaited lattice number = index_phos(ind)  ! give the lattice location 

    subroutine allocate_index_neighbors_phos(maxneigh,len_index_phos)

        use globals, only : nsize

        integer, intent(in) :: maxneigh
        integer, intent(in) :: len_index_phos

        allocate(indexneighbor(nsize,maxneigh))
        allocate(inverse_indexneighbor_phos(len_index_phos,nsize))

    end subroutine allocate_index_neighbors_phos



    ! Calculate indexneighbor(idx,k) and inverse_indexneighbor(idx 
    ! indexneighbor(idx,k) = return lattice cell index of neighbor number k of lattice cell index idx
    ! range : 1<=idx<=len_phos 
    !       : 1<=k<=(2*rangecutoff+1)**3 
    !       : k=1 is same lattice cell         

    subroutine make_table_index_neighbors_phos(rangecutoff,len_index_phos,index_phos)

        integer, intent(in) :: rangecutoff
        integer, intent(in) :: len_index_phos
        integer, intent(in) :: index_phos(len_index_phos)

        integer :: idx, i, j, k, ix, iy, iz, deltaix,deltaiy, deltaiz
        integer :: neighbornumber, idxneigh
        integer :: ip, jp, kp, n

        do n=1,len_index_phos    ! loop over element of unique phospate lattice positions 
            
            idx = index_phos(n)  ! give the lattice number of list in unique phosphate 

            ! use hash table to get coordiantes 
            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            neighbornumber=0 ! 

            do deltaix=-rangecutoff,rangecutoff
                i=ix+deltaix
                do deltaiy=-rangecutoff,rangecutoff
                    j=iy+deltaiy
                    do deltaiz=-rangecutoff,rangecutoff
                        k=iz+deltaiz
                        ! apply pbc 
                        ip=ipbc(i,nx)
                        jp=ipbc(j,ny)
                        kp=ipbc(k,nz)
                        neighbornumber=neighbornumber+1
                        idxneigh=coordtoindex(ip,jp,kp) ! index of neighbour
                        indexneighbor(idx,neighbornumber) = idxneigh
                        inverse_indexneighbor_phos(n,idxneigh) = neighbornumber

                    enddo
                enddo
            enddo

        enddo

    end subroutine make_table_index_neighbors_phos


    ! Assins positions to position_graft
    ! real(dp) position_graft(ngr,2) holds position  

    subroutine init_graftpoints()

        use random 

        integer :: i, j, ig
        real(dp) :: rnd
        real(dp) :: u, v

        ! array of regular graft points positions

        ig = 1   
        do i=1,ngrx            ! location graft points 
            do j=1,ngry
                position_graft(ig,1) = (i-0.5_dp)*delta*ngr_freq
                position_graft(ig,2) = (j-0.5_dp)*delta*ngr_freq
                ig = ig + 1
            enddo
        enddo
        
        ! array of regular graft points plus random shift 

        if(isRandom_pos_graft) then
 
            seed=seed_graft         ! randomize

            do ig=1,ngr          
                rnd = (rands(seed)-0.5_dp)*ngr_freq/scale_ran_step
                position_graft(ig,1) = position_graft(ig,1) + delta*rnd
                rnd = (rands(seed)-0.5_dp)*ngr_freq/scale_ran_step
                position_graft(ig,2) = position_graft(ig,2) + delta*rnd
            enddo
        
        endif  
      
        ! position_graft is on a regular square pattern applicable to geometry ="cubic" 
        ! For hexagonal/oblique geometry one needs to tranform (u,v) to (x,y),
        ! then above square pattern becomes a hexagonal/oblique pattern. 

        if(geometry=="prism") then 

            do ig=1,ngr   
                u = position_graft(ig,1)
                v = position_graft(ig,2)  
                position_graft(ig,1) = xt(u, v)  
                position_graft(ig,2) = yt(u, v)
            enddo    
        
        endif   

    end subroutine init_graftpoints

    ! Writes position_graft(ngr,2) to a file only if DEBUG==.true. or rank ==

    subroutine write_graftpoints(info)

        use globals, only : DEBUG
        use myutils, only : newunit,lenText

        integer, intent(out) :: info

        character(len=lenText) :: fname, istr
        integer :: ios, un_pgpt
        character(len=100) :: io_msg
        integer :: ig 
        
        
        !if(DEBUG.or.rank==0) then
        if(.true.) then
            info = 0

            write(fname,'(A18)')'positiongraft-rank'
            write(istr,'(I4)')rank
            fname=trim(fname)//trim(adjustl(istr))//'.dat'
            open(unit=newunit(un_pgpt),file=fname,iostat=ios,iomsg=io_msg)
            if(ios >0 ) then
                print*, 'Error opening positiongraftpt.dat file : iostat =', ios
                print*, 'Error message : ',trim(io_msg)
                info = 1 !myio_err_inputfile
                return
            endif

            do ig=1,ngr          
                write(un_pgpt,*)position_graft(ig,1),position_graft(ig,2) 
            enddo

            close(un_pgpt)

        else

            info = 1
        
        endif         

    end subroutine write_graftpoints


    ! Set up and assigns variable asociated with grafting points 

    subroutine set_var_graftpoints(info)
   
        integer, intent(out) :: info

        info = 0
         
        ngrx = int(nx/ngr_freq)
        ngry = int(ny/ngr_freq)
        ngr  = int(nx/ngr_freq)*int(ny/ngr_freq) ! number of surface elements to be end-grafted with chains

        ! check ngr_freq compatible with nx and ny 
        
        if(.not.(mod(nx,ngr_freq).eq.0))then 
             print*,"ngr test for x-direction failed: exiting"
             print*,"nx= ",nx," ngr_freq = ",ngr_freq
             info = 1
             stop
        endif  

        if(.not.(mod(ny,ngr_freq).eq.0)) then 
             print*,"ngr test for y-direction  failed: exiting"
             print*,"ny= ",ny," ngr_freq = ",ngr_freq
             info = 1
             stop
        endif      
        
        !  check if values of nset_per_graft and numproc are compatible
        
        if((ngr*nset_per_graft)/=numproc) then
            print*,"nset_per_graft test failed: exiting"
            print*,"nset_per_graft=",nset_per_graft,"ngr=",ngr,"numproc=",numproc
            info = 1
            stop
        endif

        allocate(position_graft(ngr,2)) ! only after ngr has been established position_graft can be allocated

    end subroutine set_var_graftpoints


    subroutine make_graftpoints()

        integer :: info

        call set_var_graftpoints(info)
        call init_graftpoints()
        call write_graftpoints(info)

    end subroutine make_graftpoints


end module volume
  
