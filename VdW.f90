!    .. module file for VdW

module VdW
    
    use precision_definition
    implicit none

    integer, parameter :: range = 2 
    integer, parameter :: MCsteps = 100000000
    
    integer, parameter :: VdW_err_allocation = 1
    integer, parameter :: VdW_err_vdwcoeff   = 2
    integer, parameter :: Vdw_err_systype    = 3
    integer, parameter :: VdW_err_inputfile  = 4
    integer, parameter :: Vdw_err_vdwcoeff_exist = 5 
    integer, parameter :: Vdw_err_vdwcoeff_not_exist = 6

    private
 
    public :: VdW_energy,make_VdWeps

contains

 

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





function VdW_energy(rhopol)result(EVdW)

    real(dp), intent(in) :: rhopol(:,:)

    real(dp) :: EVdW

    EVdW =0.0_dp

end function

! reads VdWeps.in file
! info = 0 ; read succes
! info =  VdW_err_inputfile read failure 
 
subroutine read_VdWeps(info)
    
    use globals, only : nsegtypes , runtype
    use parameters, only : VdWeps, VdWepsin 
    use myutils, only : newunit

    integer,  intent(out), optional :: info
     !     .. local variables 
    character(len=9) :: fname
    integer :: ios, un_input  ! un = unit number
    character(80) :: str
    integer :: s,t, line
   
    if (present(info)) info=0

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
    
    if (present(info)) info = 0

    call allocate_VdWeps()
    call read_VdWeps(info)

    ! special assign 
    call set_VdWepsAAandBB
    call set_VdWepsin

end subroutine


end module VdW

      
