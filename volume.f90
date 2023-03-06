!     makes volume elements 
!     for spherical coordinates 
module volume      

    use precision_definition
    use mpivars

    implicit none
  
    !     .. variables
   
    real(dp) :: delta               ! delta  spacing of lattice site in x-, y- and z-direction  
    integer :: nz                   ! nz number of lattice sites in z-direction 
    integer :: nx                   ! nx number of lattice sites in x-direction 
    integer :: ny                   ! ny number of lattice sites in y-direction 
    integer :: nsurf                ! nsurf number of lattice site  at z=0 ore z=nz*delta
    real(dp) :: volcell             ! volcell=delta**3
    real(dp) :: areacell            ! areacell=delta**2
    real(dp) :: areasurf            ! area surfaces spanned in x- and y- direction
    real(dp) :: gamma               ! angle between oblique basis vectors u and v: cubic = gamma=90=pi/2 hexagonl gamma=60=2pi/3                          
    real(dp) :: beta                ! related beta = (pi/2- gamma)/2, angle between basis vector u and x and v and y
    real(dp) :: cos_two_beta        ! sqrt(cos(beta)**2 - sin(beta)**2)=cos(2beta) scales u and v coordinates 
    real(dp) :: sin_two_beta        ! sin(2beta)  
    character(len=11) :: geometry

    ! variable for grafting position

    integer :: sgraftpts(3)         ! triplet of units number of histone that is rotated into fixed orientation
    integer, dimension(:), allocatable :: sRg ! unit numbers of AA in histone close to cm of : nmumber of number = nnucl

    ! hash table
    integer, dimension(:,:,:), allocatable :: coordtoindex 
    integer, dimension(:,:), allocatable   :: indextocoord

    private
    
    public :: delta,nx,ny,nz,volcell,areacell,geometry, sgraftpts, nsurf
    public :: gamma,cos_two_beta, sin_two_beta
    public :: coordtoindex,indextocoord
    public :: coordinateFromLinearIndex, linearIndexFromCoordinate
    public :: xt, yt, ut, vt, ipbc
    public :: make_geometry

contains
    
    subroutine init_lattice

        use globals, only : nsize, systype
        use mathconst

        implicit none 

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

        areasurf=nsurf*delta*delta
      
        if(systype=="nucl_ionbin_sv") then     ! use of hash table 
            call allocate_hashtable(nx,ny,nz)
            call make_hashtable
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

        integer :: idx, ix, iy, iz
            
        do idx=1,nsize
            
            call coordinateFromLinearIndex(idx, ix, iy,iz)

            coordtoindex(ix,iy,iz)=idx
            
            indextocoord(idx,1)=ix
            indextocoord(idx,2)=iy
            indextocoord(idx,3)=iz

        enddo
        
    end subroutine make_hashtable

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

end module volume
  
