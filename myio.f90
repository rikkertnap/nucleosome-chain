module myio

    use precision_definition
    implicit none

    ! return error values

    integer, parameter ::  myio_err_systype   = 1
    integer, parameter ::  myio_err_runtype   = 2
    integer, parameter ::  myio_err_geometry  = 3
    integer, parameter ::  myio_err_method    = 4
    integer, parameter ::  myio_err_chaintype = 5
    integer, parameter ::  myio_err_domain    = 6
    integer, parameter ::  myio_err_inputfile = 7
    integer, parameter ::  myio_err_input     = 8
    integer, parameter ::  myio_err_bcflag    = 9
    integer, parameter ::  myio_err_label     = 10
    integer, parameter ::  myio_err_VdWeps    = 11
    integer, parameter ::  myio_err_nsegtypes = 12
    integer, parameter ::  myio_err_chainmethod = 13
    integer, parameter ::  myio_err_chainsfile = 14
    integer, parameter ::  myio_err_energyfile = 15
    integer, parameter ::  myio_err_dielect   = 16
    integer, parameter ::  myio_err_graft     = 17
    integer, parameter ::  myio_err_index     = 18
    integer, parameter ::  myio_err_conf      = 19
    integer, parameter ::  myio_err_nseg      = 20
    integer, parameter ::  myio_err_inputlabel = 21
    integer, parameter ::  myio_err_readfile  = 22
    integer, parameter ::  myio_err_file_exist = 23


    integer :: num_cNaCl   ! number of salt concentration considered
    integer :: num_cMgCl2
    integer :: num_cKCl

    real(dp), dimension(:), allocatable, target :: cNaCl_array ! salt concentrations
    real(dp), dimension(:), allocatable, target :: cMgCl2_array
    real(dp), dimension(:), allocatable, target :: cKCl_array

    ! varialble for allow steps in between given salt concentrations
    integer :: maxlist_step

    ! unit number
    integer :: un_sys,un_xpolAB,un_xsol,un_xNa,un_xCl,un_xK,un_xCa,un_xMg,un_xNaCl,un_xKCl
    integer :: un_xOHmin,un_xHplus,un_fdisA,un_fdisB,un_psi,un_charge, un_xpair, un_rhopolAB, un_fe
    integer :: un_dip ,un_dielec,un_xpolABz, un_xpol, un_fdis, un_fdisP, un_angle, un_dist, un_fdision

    ! format specifiers
    character(len=80), parameter  :: fmt = "(A9,I1,A5,ES25.16)"
    character(len=80), parameter  :: fmt1reals = "(ES25.16E3)"
    character(len=80), parameter  :: fmt2reals = "(2ES25.16E3)"
    character(len=80), parameter  :: fmt3reals = "(3ES25.16E3)"
    character(len=80), parameter  :: fmt4reals = "(4ES25.16E3)"
    character(len=80), parameter  :: fmt5reals = "(5ES25.16E3)"
    character(len=80), parameter  :: fmt6reals = "(6ES25.16E3)"

    private

    public :: read_inputfile, output_individualcontr_fe, output, compute_vars_and_output, write_chain_config
    public :: myio_err_chainsfile, myio_err_energyfile, myio_err_chainmethod, myio_err_geometry
    public :: myio_err_graft, myio_err_index, myio_err_conf, myio_err_nseg, myio_err_readfile
    public :: num_cNaCl, num_cMgCl2, num_cKCl, cNaCl_array,  cMgCl2_array, cKCl_array
    public :: set_value_NaCl, set_value_MgCl2, set_value_KCl
    public :: set_value_isVdW_on_values
    public :: maxlist_step

contains

subroutine read_inputfile(info)

    use globals
    use parameters
    use surface
    use myutils, only : newunit
    use chains, only : sgraftpts

    integer, intent(out), optional :: info

    ! .. local arguments

    integer :: info_sys, info_bc, info_run, info_geo, info_meth, info_chaintype, info_combi, info_VdWeps
    integer :: info_chainmethod, info_dielect
    character(len=8) :: fname
    integer :: ios,un_input  ! un = unit number
    character(len=100) :: buffer, label
    integer :: pos
    integer :: line
    logical :: isSet_maxnchains, isSet_maxnchainsxy, isSet_precondition, isSet_savePalpha,  isSet_EnergyShift
    logical :: isSet_maxfkfunevals, isSet_maxniter

    if (present(info)) info = 0

    !     .. reading in of variables from file
    write(fname,'(A8)')'input.in'
    open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening input.in file : iostat =', ios
        if (present(info)) info = myio_err_inputfile
        return
    endif

    ! defaults

    isSet_maxnchains  =.false.
    isSet_maxnchainsxy=.false.
    isSet_precondition=.false.
    isSet_savePalpha  =.false.
    isSet_EnergyShift =.false.
    isSet_maxfkfunevals =.false.
    isSet_maxniter     =.false.
    write_mc_chains   =.false.
    write_struct      =.false.
    write_rotations   =.false.
    
    ! default concentrations
    cKCl=0.0_dp
    cCaCl2=0.0_dp
    cMgCl2=0.0_dp
    cRbCl=0.0_dp
    
    ! init surface charge 
    sigmaSurfL = 0.0_dp
    sigmaSurfR = 0.0_dp
    bcflag(LEFT) = "cc"
    bcflag(RIGHT) = "cc"

    ios = 0
    line = 0

    ! ios<0 : if an end of record condition is encountered or if an end of file condition was detected.
    ! ios>0 : if an error occured
    ! ios=0 : otherwise.

    do while (ios == 0)

        read(un_input, '(A)', iostat=ios) buffer

        if (ios == 0) then

            line = line + 1

            !  Split label and data based on first occurence of a whitespace
            pos = scan(buffer, '     ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            select case (label) !list-directed The character variable is treated as an 'internal file'
            case ('method')
                read(buffer, *,iostat=ios) method
            case ('systype')
                read(buffer, *,iostat=ios) systype
            case ('runtype')
                read(buffer, *,iostat=ios) runtype
            case ('chainmethod')
                read(buffer,*,iostat=ios) chainmethod
            case ('chaintype')
                read(buffer,*,iostat=ios) chaintype
            case ('isChainEnergyFile')
                read(buffer,*,iostat=ios) isChainEnergyFile
            case ('isEnergyShift')
                read(buffer,*,iostat=ios) isEnergyShift
                isSet_EnergyShift=.true.
            case ('tolerance')
                read(buffer,*,iostat=ios) tol_conv
            case ('infile')
                read(buffer,*,iostat=ios) infile              ! guess  1==yes
            case ('pH%val')
                read(buffer,*,iostat=ios) pH%val
            case ('pH%min')
                read(buffer,*,iostat=ios) pH%min
            case ('pH%max')
                read(buffer,*,iostat=ios) pH%max
            case ('pH%stepsize')
                read(buffer,*,iostat=ios) pH%stepsize
            case ('pH%delta')
                read(buffer,*,iostat=ios) pH%delta
            case ('KionNa')
                read(buffer,*,iostat=ios) KionNa
            case ('KionK')
                read(buffer,*,iostat=ios) KionK
            case ('cNaCl')
                read(buffer,*,iostat=ios) cNaCl
            case ('cKCl')
                read(buffer,*,iostat=ios) cKCl
            case ('cRbCl')
                read(buffer,*,iostat=ios) cRbCl
            case ('cCaCl2')
                read(buffer,*,iostat=ios) cCaCl2
            case ('cMgCl2')
                read(buffer,*,iostat=ios) cMgCl2
            case ('nsize')
                read(buffer,*,iostat=ios) nsize
            case ('nnucl')
                read(buffer,*,iostat=ios) nnucl   
            case ('nseg')
                read(buffer,*,iostat=ios) nseg
            case ('nsegAA')
                read(buffer,*,iostat=ios) nsegAA
            case ('nsegtypes')
                read(buffer,*,iostat=ios) nsegtypes        ! carefully need to be overwriiten depending on value systype and or chainmethod
            case ('nsegtypesAA')
                read(buffer,*,iostat=ios) nsegtypesAA  
            case ('cuantas')
                read(buffer,*,iostat=ios) max_confor
            case ('chainperiod')
                read(buffer,*,iostat=ios) chainperiod
            case('vpolfname')
                read(buffer,*,iostat=ios) vpolfname
            case('pKafname')
                read(buffer,*,iostat=ios) pKafname
            case('pKaionfname')
                read(buffer,*,iostat=ios) pKaionfname
            case('typesfname')
                read(buffer,*,iostat=ios) typesfname
            case('lsegfname')
                read(buffer,*,iostat=ios) lsegfname
            case('segcmfname')
                read(buffer,*,iostat=ios) segcmfname
            case('mtpdbfname')
                read(buffer,*,iostat=ios) mtpdbfname
            case('vnuclfname')
                read(buffer,*,iostat=ios) vnuclfname
            case('orientfname')
                read(buffer,*,iostat=ios) orientfname
            case ('nx')
                read(buffer,*,iostat=ios) nx
            case ('ny')
                read(buffer,*,iostat=ios) ny
            case ('nz')
                read(buffer,*,iostat=ios) nz
            case ('verboseflag  ')
                read(buffer,*,iostat=ios) verboseflag
            case ('delta')
                read(buffer,*,iostat=ios) delta
            case('unit_conv')
                read(buffer,*,iostat=ios) unit_conv
            case ('geometry')
                read(buffer,*,iostat=ios) geometry
            case ('sgraft')
                read(buffer,*,iostat=ios) sgraftpts(1),sgraftpts(2),sgraftpts(3)  
            case ('gamma')
                read(buffer,*,iostat=ios) gamma
            case ('write_mc_chains')
                read(buffer,*,iostat=ios) write_mc_chains
            case ('precondition')
                read(buffer,*,iostat=ios) precondition
                isSet_precondition=.true.
            case ('maxniter')
                read(buffer,*,iostat=ios) maxniter
                isSet_maxniter=.true.
            case ('maxfkfunevals')
                read(buffer,*,iostat=ios) maxfkfunevals
                isSet_maxfkfunevals=.true.
            case ('dielect_env')
                read(buffer,*,iostat=ios) dielect_env
            case ('VdWscale%val')
                read(buffer,*,iostat=ios) VdWscale%val
            case ('VdWscale%min')
                read(buffer,*,iostat=ios) VdWscale%min
            case ('VdWscale%max')
                read(buffer,*,iostat=ios) VdWscale%max
            case ('VdWscale%stepsize')
                read(buffer,*,iostat=ios) VdWscale%stepsize
            case ('VdWscale%delta')
                read(buffer,*,iostat=ios) VdWscale%delta
            case ('pKd%val')
                read(buffer,*,iostat=ios) pKd%val
            case ('pKd%min')
                read(buffer,*,iostat=ios) pKd%min
            case ('pKd%max')
                read(buffer,*,iostat=ios) pKd%max
            case ('pKd%stepsize')
                read(buffer,*,iostat=ios) pKd%stepsize
            case ('pKd%delta')
                read(buffer,*,iostat=ios) pKd%delta
            case ("write_struct")
                read(buffer,*,iostat=ios) write_struct
            case ('write_rotations')
                read(buffer,*,iostat=ios) write_rotations
            case ('dielectscale%val')
                read(buffer,*,iostat=ios) dielectscale%val
            case ('dielectscale%min')
                read(buffer,*,iostat=ios) dielectscale%min
            case ('dielectscale%max')
                read(buffer,*,iostat=ios) dielectscale%max
            case ('dielectscale%stepsize')
                read(buffer,*,iostat=ios) dielectscale%stepsize
            case ('dielectscale%delta')
                read(buffer,*,iostat=ios) dielectscale%delta
            case default
                if(pos>1) then
                    print *, 'Invalid label at line', line  ! empty lines are skipped
                    if (present(info)) info = myio_err_inputlabel
                    return
                endif
            end select
        endif
    enddo

    if(ios >0 ) then
        print*, 'Error parsing file : iostat =', ios
        print*, 'Read error at line =', line
        if (present(info)) info = myio_err_inputfile
        return
    endif

    close(un_input)


    ! .. check error flag

    call check_value_systype(systype,info_sys)
    if (info_sys == myio_err_systype) then
        if (present(info)) info = info_sys
        return
    endif

    call check_value_runtype(runtype,info_run)
    if (info_sys == myio_err_runtype) then
        if (present(info)) info = info_run
        return
    endif

    call check_value_geometry(geometry,info_geo)
    if (info_geo == myio_err_geometry) then
        if (present(info)) info = info_geo
        return
    endif

    call check_value_bcflag(bcflag,info_bc)
    if (info_bc == myio_err_bcflag) then
        if (present(info)) info = info_bc
        return
    endif

    call check_value_method(method,info_meth)
    if (info_meth == myio_err_method) then
        if (present(info)) info = info_meth
        return
    endif

    call check_value_chaintype(chaintype,info_chaintype)
    if (info_chaintype == myio_err_chaintype) then
        if (present(info)) info = info_chaintype
        return
    endif

    call check_value_chainmethod(chainmethod,info_chainmethod)
    if (info_chainmethod == myio_err_chainmethod) then
        if (present(info)) info = info_chainmethod
        return
    endif

    call check_value_dielect_env(dielect_env,info_dielect)
    if (info_dielect == myio_err_dielect) then
        if (present(info)) info = info_dielect
        return
    endif

    !  .. set input values

    call set_value_isVdW(systype,isVdW)
    call set_value_isVdWintEne(systype, isVdWintEne)
    call set_value_nsegtypes(nsegtypes,chaintype,systype,info)

    !  .. set input values: use of default or based on input.in

    call set_value_int_var(maxnchainsrotations,isSet_maxnchains,1)
    call set_value_int_var(maxnchainsrotationsxy,isSet_maxnchainsxy,12)
    call set_value_logical_var(precondition,isSet_precondition,.false.)
    call set_value_logical_var(isEnergyShift,isSet_EnergyShift,.false.)
    call set_value_int_var(maxfkfunevals,isSet_maxfkfunevals,1000)
    call set_value_int8_var(maxniter,isSet_maxniter,int(1000,8))

!    call set_value_maxnchains(maxnchainsrotations,isSet_maxnchains)
!    call set_value_maxnchainsxy(maxnchainsrotationsxy,isSet_maxnchainsxy)
!    call set_value_precondition(precondition,isSet_precondition)
!    call set_value_isEnergyShift(isEnergyShift,isSet_EnergyShift)
!    call set_value_maxfkfunevals(maxfkfunevals,isSet_maxfkfunevals)
!    call set_value_maxniter(maxniter,isSet_maxniter)

    ! after set_value_isVdW
    call check_value_VdWeps(systype,isVdW,info_VdWeps)
    if (info_VdWeps == myio_err_VdWeps) then
        if (present(info)) info = info_VdWeps
        return
    endif

    ! .. override certain input values
    if(systype=="brushdna".or.systype=="brushborn".or.systype=="brush_mul".or.&
        systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv") then
        KionNa   = 0.0_dp
        KionK    = 0.0_dp
    endif

end subroutine read_inputfile


subroutine check_value_systype(systype,info)

    character(len=15), intent(in) :: systype
    integer, intent(out),optional :: info

    character(len=15) :: systypestr(11)
    integer :: i
    logical :: flag

    ! permissible values of systype

    systypestr(1)="elect"
    systypestr(2)="neutral"
    systypestr(3)="brush_mul"
    systypestr(4)="brush_mulnoVdW"
    systypestr(5)="brushdna"
    systypestr(6)="brushborn"
    systypestr(7)="bulk water"
    systypestr(8)="neutralnoVdW"
    systypestr(9)="nucl_ionbin"
    systypestr(10)="nucl_ionbin_sv"
    systypestr(11)="nucl_neutral_sv"

    flag=.FALSE.

    do i=1,11
        if(systype==systypestr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of systype is not permissible"
        print*,"systype = ",systype
        if (present(info)) info = myio_err_systype
        return
    endif

end subroutine check_value_systype


subroutine check_value_runtype(runtype,info)

    character(len=15), intent(in) :: runtype
    integer, intent(out),optional :: info

    character(len=15) :: runtypestr(6)
    integer :: i
    logical :: flag

    ! permissible values of runtype

    runtypestr(1)="inputcspH" ! used to be rangpH
    runtypestr(2)="inputMgpH"
    runtypestr(3)="inputcsKClpH"
    runtypestr(4)="rangepKd"
    runtypestr(5)="rangeVdWeps"
    runtypestr(6)="rangedielect"

    flag=.FALSE.

    do i=1,6
        if(runtype==runtypestr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of runtype is not permissible"
        print*,"runtype = ",runtype
        if (present(info)) info = myio_err_runtype
        return
    endif

end subroutine check_value_runtype

subroutine check_value_bcflag(bcflag,info)

    use globals, only : LEFT, RIGHT

    character(len=2), intent(in), dimension(2) :: bcflag
    integer, intent(out), optional :: info

    character(len=15) :: bcvalues(2,5)
    integer :: i
    logical :: flag

    if (present(info)) info = 0 ! init 

    ! permissible values of bcflag

    bcvalues(RIGHT,1)="qu"
    bcvalues(RIGHT,2)="cl"
    bcvalues(RIGHT,3)="ca"
    bcvalues(RIGHT,4)="ta"
    bcvalues(RIGHT,5)="cc"
    bcvalues(LEFT,1)="ta"
    bcvalues(LEFT,2)="cc"

    flag=.FALSE.

    do i=1,5
        if(bcflag(RIGHT)==bcvalues(RIGHT,i)) flag=.TRUE.
    enddo
    if (flag.eqv. .FALSE.) then
        print*,"Error value of bcflag is not permissible"
        print*,"bcflag(RIGHT) = ",bcflag(RIGHT)
        if (present(info)) info = myio_err_systype
        !if (present(fcnname)) print*,"Error in ",fcnname
        stop
    endif

    flag=.FALSE.
    do i=1,2
        if(bcflag(LEFT)==bcvalues(LEFT,i)) flag=.TRUE.
    enddo
    if (flag.eqv. .FALSE.) then
        print*,"Error value of bcflag is not permissible"
        print*,"bcflag(LEFT) = ",bcflag(LEFT)
        if (present(info)) info = myio_err_bcflag
        !if(present(fcnname)) print*,"Error in ",fcnname
        ! stop
    endif

end subroutine check_value_bcflag


subroutine set_value_NaCl(runtype,info)

    use myutils, only : newunit, lenText, LogUnit, print_to_log

    character(len=12), intent(in) :: runtype
    integer, intent(out) :: info

    ! local variables
    character(len=7) :: fname
    integer :: ios, i, un_cs
    logical :: exist
    character(len=lenText) :: text, istr

    info = 0 ! init 

    if(runtype=="inputcspH".or.runtype=="rangedielect") then

        !     .. read salt concentrations from file
        write(fname,'(A7)')'salt.in'
        inquire(file=fname,exist=exist)

        if(exist) then
            open(unit=newunit(un_cs),file=fname,iostat=ios,status='old')
        else
            text='File :'//trim(fname)//' does not exit'
            call print_to_log(LogUnit,text)
            print*,text
            info = myio_err_file_exist
            return
        endif

        if(ios > 0 ) then
            write(istr,'(I4)')ios
            text='Error opening file salt.in : iostat ='//trim(istr)
            call print_to_log(LogUnit,text)
            print*,text
            info = myio_err_inputfile
            return
        endif

        read(un_cs,*)num_cNaCl ! read number of salt concentration form file
        allocate(cNaCl_array(num_cNaCl))

        do i=1,num_cNaCl     ! read value salt concentration
             read(un_cs,*)cNaCl_array(i)
        enddo
        close(un_cs)
    
    else
        text='Error wrong runtype in set_value_NaCl'
        call print_to_log(LogUnit,text)
        print*, text 
        info = myio_err_runtype
    
    endif

end subroutine  set_value_NaCl


subroutine set_value_KCl(runtype,info)

    use myutils, only : newunit, lenText, LogUnit, print_to_log

    character(len=12), intent(in) :: runtype
    integer, intent(out) :: info

    ! local variables
    character(len=7) :: fname
    integer :: ios, i, un_cs
    logical :: exist
    character(len=lenText) :: text, istr

    info = 0

    if(runtype=="inputcsKClpH") then

        !     .. read salt concentrations from file
        write(fname,'(A7)')'salt.in'
        inquire(file=fname,exist=exist)

        if(exist) then
            open(unit=newunit(un_cs),file=fname,iostat=ios,status='old')
        else
            text='File :'//trim(fname)//' does not exit'
            call print_to_log(LogUnit,text)
            print*, text 
            info = myio_err_file_exist
            return
        endif

        if(ios > 0 ) then
            write(istr,'(I4)')ios
            text='Error opening file salt.in : iostat ='//trim(istr)
            call print_to_log(LogUnit,text)
            print*, text 
            info = myio_err_inputfile
            return
        endif

        read(un_cs,*)num_cKCl ! read number of salt concentration form file
        allocate(cKCl_array(num_cKCl))

        do i=1,num_cKCl     ! read value salt concentration
            read(un_cs,*)cKCl_array(i)
        enddo
        close(un_cs)

    else
        text='Error wrong runtype in set_value_KCl'
        call print_to_log(LogUnit,text)
        print*, text 
        info = myio_err_inputfile
    endif

end subroutine set_value_KCl

subroutine set_value_MgCl2(runtype,info)

    use myutils, only : newunit

    character(len=12), intent(in) :: runtype
    integer, intent(out) :: info

    ! local variables
    character(len=9) :: fname
    integer :: ios, i, un_cs
    logical :: exist

    info=0

    if(runtype=="inputMgpH".or.runtype=="rangepKd".or.runtype=="rangeVdWeps") then

        !     .. read salt concentrations from file
        write(fname,'(A9)')'saltMg.in'
        inquire(file=fname,exist=exist)

        if(exist) then
            open(unit=newunit(un_cs),file=fname,iostat=ios,status='old')
        else
            print*,' File :',fname,' does not exit'
            info = myio_err_file_exist
            return
        endif

        if(ios > 0 ) then
            print*, 'Error opening file saltMg.in : iostat =', ios
            info = myio_err_inputfile
            return
        endif

        read(un_cs,*)num_cMgCl2 ! read number of salt concentration form file
        allocate(cMgCl2_array(num_cMgCl2))

        do i=1,num_cMgCl2     ! read value salt concentration
            read(un_cs,*)cMgCl2_array(i)
        enddo
        close(un_cs)

    else
        
        print*,'Error wrong runtype in set_value_MgCl2'
        info = myio_err_inputfile
    
    endif

end subroutine  set_value_MgCl2


subroutine check_value_geometry(geometry,info)

    character(len=11), intent(in) :: geometry
    integer, intent(out),optional :: info

    logical :: flag
    character(len=11) :: geometrystr(2)
    integer :: i

    ! permissible values of geometry

    geometrystr(1)="cubic"
    geometrystr(2)="prism"

    flag=.FALSE.

    do i=1,2
        if(geometry==geometrystr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of geometry is not permissible"
        print*,"geometry = ",geometry
        if (present(info)) info = myio_err_geometry
        return
    endif

end subroutine


subroutine check_value_chaintype(chaintype,info)

    character(len=8), intent(in) :: chaintype
    integer, intent(out),optional :: info

    logical :: flag
    character(len=8) :: chaintypestr(6)
    integer :: i

    ! permissible values of chaintype

    chaintypestr(1)="diblockA"
    chaintypestr(2)="diblockB"
    chaintypestr(3)="altA"
    chaintypestr(4)="altB"
    chaintypestr(5)="copolyAB"
    chaintypestr(6)="multi"

    flag=.FALSE.

    do i=1,6
        if(chaintype==chaintypestr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of chaintype is not permissible"
        print*,"chaintype = ",chaintype
        if (present(info)) info = myio_err_chaintype
        return
    endif

end subroutine check_value_chaintype

subroutine check_value_chainmethod(chainmethod,info)

    character(len=15), intent(in) :: chainmethod
    integer, intent(out),optional :: info

    logical :: flag
    character(len=15) :: chainmethodstr(2)
    integer :: i

    ! permissible values of chainmethod

    chainmethodstr(1)="MC"
    chainmethodstr(2)="FILE_XYZ"

    flag=.FALSE.

    do i=1,2
        if(chainmethod==chainmethodstr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of chainmethod is not permissible"
        print*,"chainmethod = ",chainmethod
        if (present(info)) info = myio_err_chainmethod
        return
    endif

end subroutine check_value_chainmethod


subroutine check_value_method(method,info)

    character(len=8), intent(in) :: method
    integer, intent(out),optional :: info

    character(len=8) :: methodstr
    integer :: i
    logical :: flag

    ! permissible values of runtype

    methodstr="kinsol"

    flag=.FALSE.

    if (method==methodstr) flag=.TRUE.

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of method is not permissible"
        print*,"method = ",method
        if (present(info)) info = myio_err_method
        return
    endif

end subroutine check_value_method

subroutine check_value_dielect_env(dielect_env,info)

    character(len=15), intent(in) :: dielect_env
    integer, intent(out),optional :: info

    logical :: flag
    character(len=15) :: dielect_env_str(3)
    integer :: i

    ! permissible values of dielect_env

    dielect_env_str(1)="constant"
    dielect_env_str(2)="linear"
    dielect_env_str(3)="MaxwellGarnett"

    flag=.FALSE.

    do i=1,3
        if(dielect_env==dielect_env_str(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of dielect_env is not permissible"
        print*,"dielect_env = ",dielect_env
        if (present(info)) info = myio_err_dielect
        return
    endif

end subroutine check_value_dielect_env



subroutine check_value_VdWeps(systype,isVdW,info)

    logical, intent(in) :: isVdW
    character(len=15), intent(in) :: systype
    integer, intent(out), optional :: info

    character(len=15) :: systypestr(7)
    integer :: i
    logical :: flag

    flag=.false.

    if(isVdW) then ! check for correct combination systype

        systypestr(1)="neutral"
        systypestr(2)="brush_mul"
        systypestr(3)="brushdna"
        systypestr(4)="brushborn"
        systypestr(5)="nucl_ionbin"
        systypestr(6)="nucl_ionbin_sv"
        systypestr(7)="nucl_neutral_sv"

        do i=1,7
            if(systype==systypestr(i)) flag=.true.
        enddo

    else  ! isVdW=.false. so oke
        flag=.true.
    endif

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error:  combination systype and isVdW is not permissible"
        print*,"systype = ",systype, " isVdW =",isVdW
        if (present(info)) info = myio_err_VdWeps
        return
    end if

end subroutine check_value_VdWeps



! override input value nzmin
! Sets nzmin=nzmax for which runtype that do not loop over nz
! ensuring that the output files are properly close once

subroutine set_value_nzmin(runtype,nzmin,nzmax)

    character(len=15), intent(in) :: runtype
    integer, intent(inout) :: nzmin
    integer, intent(in) :: nzmax

    if (runtype/="rangedist") nzmin=nzmax

end subroutine

! isVdW = .true. then uses position dependent Van der Waals interaction

subroutine set_value_isVdW(systype, isVdW)

    character(len=15), intent(in) :: systype
    logical, intent(inout)  :: isVdW

    character(len=15) :: systypestr(3)
    integer :: i

    isVdW=.True.

    ! systype that NOT involve VdW interactions 

    systypestr(1)="elect"
    systypestr(2)="neutralnoVdW"
    systypestr(3)="brush_mulnoVdW"

    do i=1,3
        if(systype==systypestr(i)) isVdW=.FALSE.
    enddo

end subroutine

! isVdWintEne = .true. compute internal VdW energy chain
! allways true

subroutine set_value_isVdWintEne(systype, isVdWintEne)

    character(len=15), intent(in) :: systype
    logical, intent(inout)  :: isVdWintEne

    character(len=15) :: systypestr(2)
    integer :: i

    ! all systype that involve internal VdW chain energy

    isVdWintEne=.true.

end subroutine

! changes value isVdW to .false. if and only if all values of VdWeps  
! are smalller than threshold 
! pre : VdWeps need to be initialized see module VdW

subroutine set_value_isVdW_on_values(nsegtypes, VdWeps, isVdW)

    use parameters, only : VdWepsilon

    integer, intent(in) :: nsegtypes
    real(dp), intent(in) :: VdWeps(:,:)
    logical, intent(inout)  :: isVdW

    integer :: s,t
    logical :: isallVdWzero

    isallVdWzero=.true.

    do s=1,nsegtypes
        do t=1,nsegtypes
            if (abs(VdWeps(s,t))>VdWepsilon) isallVdWzero=.false.
        enddo
    enddo

    if(isallVdWzero) isVdW=.false. 

end subroutine set_value_isVdW_on_values


subroutine set_value_nsegtypes(nsegtypes,chaintype,systype,info)


    integer, intent(inout) :: nsegtypes
    integer, intent(out),optional :: info
    character(len=8), intent(in) :: chaintype
    character(len=15), intent(in) :: systype

    logical :: flag
    character(len=8) :: chaintypestr(5)
    integer :: i

    
    ! two components
    ! permissible values of chaintype

    chaintypestr(1)="diblockA"
    chaintypestr(2)="diblockB"
    chaintypestr(3)="altA"
    chaintypestr(4)="altB"
    chaintypestr(5)="copolyAB"

    flag=.false.
    do i=1,5
        if(chaintype==chaintypestr(i)) flag=.true.
    enddo
    if(flag) nsegtypes=2

    if(chaintype=="multi") flag=.true.

    if(systype=="neutral") flag=.true.

    if (present(info)) info = 0

    if (flag.eqv. .false.) then
        print*,"Error: value of nsegtype could not be set for given "
        print*,"chaintype = ",chaintype
        if (present(info)) info = myio_err_nsegtypes
        return
    endif

end subroutine set_value_nsegtypes


subroutine set_value_maxnchains(maxnchainsrotations,isSet_maxnchains)

    integer, intent(inout) :: maxnchainsrotations
    logical, intent(in)  :: isSet_maxnchains

    if(.not.isSet_maxnchains) then
        maxnchainsrotations=12 ! default value
    else
        maxnchainsrotations=abs(maxnchainsrotations) !  make positive
    endif

end subroutine set_value_maxnchains



subroutine set_value_maxnchainsxy(maxnchainsrotationsxy,isSet_maxnchainsxy)

    integer, intent(inout) :: maxnchainsrotationsxy
    logical, intent(in)  :: isSet_maxnchainsxy

    if(.not.isSet_maxnchainsxy) then
        maxnchainsrotationsxy=1 ! default value
    else
        maxnchainsrotationsxy=abs(maxnchainsrotationsxy) !  make positive
    endif

end subroutine set_value_maxnchainsxy


subroutine set_value_precondition(precondition,isSet_precondition)

    logical, intent(inout) :: precondition
    logical, intent(in)  :: isSet_precondition

    if(.not.isSet_precondition) precondition=.false. ! default value


end subroutine set_value_precondition


subroutine set_value_maxniter(maxniter,isSet_maxniter)

    integer(8), intent(inout) :: maxniter
    logical, intent(in)  :: isSet_maxniter

    if(.not.isSet_maxnIter) maxniter=1000 ! default value


end subroutine set_value_maxniter


subroutine set_value_maxfkfunevals(maxfkfunevals,isSet_maxfkfunevals)

    integer, intent(inout) :: maxfkfunevals
    logical, intent(in)  :: isSet_maxfkfunevals

    if(.not.isSet_maxfkfunevals) maxfkfunevals=1000 ! default value

end subroutine set_value_maxfkfunevals



subroutine set_value_isEnergyShift(isEnergyShift,isSet_EnergyShift)

    logical, intent(inout) :: isEnergyShift
    logical, intent(in)  :: isSet_EnergyShift

    if(.not.isSet_EnergyShift) isEnergyShift=.false. ! default value

end subroutine set_value_isEnergyShift


subroutine set_value_logical_var(var,isSet_var,default_value_var)

    logical, intent(inout) :: var
    logical, intent(in) :: isSet_var
    logical, intent(in) :: default_value_var

    if(.not.isSet_var) var=default_value_var! default value

end subroutine set_value_logical_var

subroutine set_value_int_var(var,isSet_var,default_value_var)

    integer, intent(inout) :: var
    logical, intent(in) :: isSet_var
    integer, intent(in) :: default_value_var

    if(.not.isSet_var) var=default_value_var! default value

end subroutine set_value_int_var

subroutine set_value_int8_var(var,isSet_var,default_value_var)

    integer(8), intent(inout) :: var
    logical, intent(in) :: isSet_var
    integer(8), intent(in) :: default_value_var

    if(.not.isSet_var) var=default_value_var! default value

end subroutine set_value_int8_var

subroutine set_value_double_var(var,isSet_var,default_value_var)

    real(dp), intent(inout) :: var
    logical, intent(in) :: isSet_var
    real(dp), intent(in) :: default_value_var

    if(.not.isSet_var) var=default_value_var ! default value

end subroutine set_value_double_var


subroutine output()

    use globals, only : systype

    select case (systype)
    case ("elect") 

        call output_elect

    case("neutral","neutralnoVdW") 

        call output_neutral
        call output_individualcontr_fe

    case("brush_mul","brush_mulnoVdW","brushdna","brushborn",&
        "nucl_ionbin","nucl_ionbin_sv")

        call output_nucl_mul
        call output_individualcontr_fe

    case("nucl_neutral_sv")

        call output_neutral()   

    case default

        print*,"Error in output subroutine"
        print*,"Wrong value systype : ", systype

    end select

end subroutine output


subroutine output_nucl_mul

    !     .. variables and constant declaractions
    use globals
    use volume
    use parameters
    use field
    use energy
    use surface
    use myutils, only : newunit
    use chains, only : isHomopolymer, avRgsqr, avRendsqr, avbond_angle,avdihedral_angle,avnucl_spacing 
    use chains, only : type_of_charge, mapping_num_to_char, sgraftpts
    
    !     .. local arguments

    integer :: t

    !     .. output file names

    character(len=90) :: sysfilename
    character(len=90) :: xsolfilename
    character(len=90) :: xpolfilename
    character(len=90) :: xpolendfilename
    character(len=90) :: xNafilename
    character(len=90) :: xKfilename
    character(len=90) :: xCafilename
    character(len=90) :: xMgfilename
    character(len=90) :: xNaClfilename
    character(len=90) :: xKClfilename
    character(len=90) :: xClfilename
    character(len=90) :: potentialfilename
    character(len=90) :: chargefilename
    character(len=90) :: xHplusfilename
    character(len=90) :: xOHminfilename
    character(len=90) :: densfracfilename
    character(len=90) :: densfracPfilename
    character(len=90) :: densfracionfilename
    character(len=90) :: densfracionpairfilename
    character(len=90) :: anglesfilename
    character(len=90) :: spacingfilename
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

    logical :: isopen
    integer :: i,j,k          ! dummy indexes
    real(dp) :: denspol


    ! .. executable statements

    denspol=init_denspol()

   
    !     .. make label filenames f
    call make_filename_label(fnamelabel)

    sysfilename    = 'system.'//trim(fnamelabel)
    xpolfilename   = 'xpol.'//trim(fnamelabel)
    xsolfilename   = 'xsol.'//trim(fnamelabel)
    xNafilename    = 'xNaions.'//trim(fnamelabel)
    xKfilename     = 'xKions.'//trim(fnamelabel)
    xCafilename    = 'xCaions.'//trim(fnamelabel)
    xMgfilename    = 'xMgions.'//trim(fnamelabel)
    xNaClfilename  = 'xNaClionpair.'//trim(fnamelabel)
    xKClfilename   = 'xKClionpair.'//trim(fnamelabel)
    xClfilename    = 'xClions.'//trim(fnamelabel)
    potentialfilename = 'potential.'//trim(fnamelabel)
    chargefilename = 'charge.'//trim(fnamelabel)
    xHplusfilename = 'xHplus.'//trim(fnamelabel)
    xOHminfilename = 'xOHmin.'//trim(fnamelabel)
    densfracfilename = 'densityfrac.'//trim(fnamelabel)
    densfracPfilename = 'densityfracP.'//trim(fnamelabel)
    densfracionfilename = 'densityfracion.'//trim(fnamelabel)
    densfracionpairfilename = 'densityfracionpair.'//trim(fnamelabel)
    anglesfilename = 'angles.'//trim(fnamelabel)
    spacingfilename = 'spacing.'//trim(fnamelabel)

    !     .. opening files

    open(unit=newunit(un_sys),file=sysfilename)
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_psi),file=potentialfilename)
    open(unit=newunit(un_xpol),file=xpolfilename)
    open(unit=newunit(un_fdis),file=densfracfilename)
    if(systype=="brushdna".or.systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv") then
        open(unit=newunit(un_fdisP),file=densfracPfilename)
    endif
    if(systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv") open(unit=newunit(un_fdision),file=densfracionfilename)
    if(nnucl>1) open(unit=newunit(un_dist),file=spacingfilename)
    if(nnucl>2) open(unit=newunit(un_angle),file=anglesfilename)
        

    if(verboseflag=="yes") then
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xMg),file=xMgfilename)
        open(unit=newunit(un_xNaCl),file=xNaClfilename)
        open(unit=newunit(un_xKCl),file=xKClfilename)
        open(unit=newunit(un_xpair),file=densfracionpairfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
    endif


    !  .. output of bond and dihedral angles
    do i=1,nnucl-3
        write(un_angle,*)avbond_angle(i),avdihedral_angle(i)
    enddo
    if(nnucl>=3)write(un_angle,*)avbond_angle(nnucl-2)

    do i=1,nnucl-1
        write(un_dist,*)avnucl_spacing(i)
    enddo

    do i=1,nsize
        write(un_xsol,*)xsol(i)
        write(un_psi,*)psi(i)
    enddo

    do i=1,nsize
        write(un_xpol,*)xpol(i),(rhopol(i,t),t=1,nsegtypes)
    enddo

    if(systype/="nucl_ionbin".and.systype/="nucl_ionbin_sv") then
        do i=1,nsize
            write(un_fdis,*)(fdis(i,t),t=1,nsegtypes)
        enddo
    endif     

    if(systype=="brushdna".or.systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv") then
        do i=1,nsize
            write(un_fdisP,'(8ES25.16)')(fdisA(i,k),k=1,8)
        enddo
    endif

    if(systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv") then
        do t=1,nsegtypes
            if(type_of_charge(t)=="A") then 
                write(un_fdision,*)"monomer type t=",t,mapping_num_to_char(t)
                do i=1,nsize
                    write(un_fdision,*)(gdisA(i,k,t),k=1,4)
                enddo
            endif 
        enddo    
        do t=1,nsegtypes
            if(type_of_charge(t)=="B") then 
                do i=1,nsize
                    write(un_fdision,*)(gdisB(i,k,t),k=1,3)
                enddo
            endif       
        enddo
    endif

    if(verboseflag=="yes") then
        do i=1,nsize
            write(un_xNa,*)xNa(i)
            write(un_xK,*)xK(i)
            write(un_xCa,*)xCa(i)
            write(un_xMg,*)xMg(i)
            write(un_xNaCl,*)xNaCl(i)
            write(un_xKCl,*)xKCl(i)
            write(un_xpair,*)(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)xCl(i)
            write(un_charge,*)rhoq(i)
            write(un_xHplus,*)xHplus(i)
            write(un_xOHmin,*)xOHmin(i)
        enddo
    endif

    ! .. writing system information

    write(un_sys,*)'system      = nucleosome chain'
    write(un_sys,*)'version     = ',VERSION

    ! chain description
    write(un_sys,*)'chainmethod = ',chainmethod
    write(un_sys,*)'chaintype   = ',chaintype
    write(un_sys,*)"isHomopolymer= ",isHomopolymer
    write(un_sys,*)'nseg        = ',nseg
    write(un_sys,*)'nsegtypes   = ',nsegtypes
    write(un_sys,*)'nnucl       = ',nnucl
    write(un_sys,*)'sgraftpts   = ',(sgraftpts(t),t=1,3)
    write(un_sys,*)'cuantas     = ',cuantas
    write(un_sys,*)'denspol     = ',denspol

    ! system description
    write(un_sys,*)'systype     = ',systype
    write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
    write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
    write(un_sys,*)'delta       = ',delta
    write(un_sys,*)'nx          = ',nx
    write(un_sys,*)'ny          = ',ny
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'tol_conv    = ',tol_conv

    ! concentration
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    write(un_sys,*)'cMgCl2      = ',cMgCl2
    write(un_sys,*)'xsolbulk    = ',xbulk%sol
    write(un_sys,*)'xNabulk     = ',xbulk%Na
    write(un_sys,*)'xClbulk     = ',xbulk%Cl
    write(un_sys,*)'xKbulk      = ',xbulk%K
    write(un_sys,*)'xRbbulk     = ',xbulk%Rb
    write(un_sys,*)'xNaClbulk   = ',xbulk%NaCl
    write(un_sys,*)'xKClbulk    = ',xbulk%KCl
    write(un_sys,*)'xCabulk     = ',xbulk%Ca
    write(un_sys,*)'xMgbulk     = ',xbulk%Mg
    write(un_sys,*)'xHplusbulk  = ',xbulk%Hplus
    write(un_sys,*)'xOHminbulk  = ',xbulk%OHmin
    write(un_sys,*)'pHbulk      = ',pHbulk

    ! disociation constants
    do t=1,nsegtypes
        write(un_sys,*)'pKa(',t,')  = ',pKa(t)
    enddo    
    !
    if(systype=="brushdna".or.systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv"&
        .or.systype=="brushborn") then
        do k=1,7
            write(un_sys,*)'pKaAA(',k,') = ',pKaAA(k)
        enddo
    endif

    write(un_sys,*)'KionNa      = ',KionNa
    write(un_sys,*)'KionK       = ',KionK
    write(un_sys,*)'K0ionNa     = ',K0ionNa
    write(un_sys,*)'K0ionK      = ',K0ionK
    write(un_sys,*)'dielectW    = ',dielectW
    if(runtype=="rangedielect") then 
        write(un_sys,*)'dielectscale = ',dielectscale%val
    endif    
    write(un_sys,*)'lb          = ',lb
    write(un_sys,*)'T           = ',Tref

    ! charge components
    do t=1,nsegtypes
        write(un_sys,*)'zpol(',t,')     = ',zpol(t,1),zpol(t,2)
    enddo
    write(un_sys,*)'zNa         = ',zNa
    write(un_sys,*)'zCa         = ',zCa
    write(un_sys,*)'zMg         = ',zMg
    write(un_sys,*)'zK          = ',zK
    write(un_sys,*)'zCl         = ',zCl
    
    ! volume
    do t=1,nsegtypes
        write(un_sys,*)'vpol(',t,')     = ',vpol(t)
    enddo 

    write(un_sys,*)'vsol        = ',vsol
    write(un_sys,*)'vNa         = ',vNa*vsol
    write(un_sys,*)'vCl         = ',vCl*vsol
    write(un_sys,*)'vCa         = ',vCa*vsol
    write(un_sys,*)'vMg         = ',vMg*vsol
    write(un_sys,*)'vK          = ',vK*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol

    do t=1,nsegtypes
        write(un_sys,*)'lseg(',t,')     = ',lsegAA(t)
    enddo    

    ! structure and thermo 
    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'tol_conv    = ',tol_conv
    write(un_sys,*)'denspol     = ',denspol
    do t=1,nsegtypes
        write(un_sys,*)'sumphi(',t,')    = ',sumphi(t)
    enddo    
    write(un_sys,*)'check phi   = ',checkphi
    write(un_sys,*)'FEq         = ',FEq
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEel        = ',FEel
    write(un_sys,*)'FEelsurf(LEFT)  = ',FEelsurf(LEFT)
    write(un_sys,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
    write(un_sys,*)'FEbind      = ',FEbind
    write(un_sys,*)'FEVdW       = ',FEVdW
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'q           = ',q
    write(un_sys,*)'avRgsqr     = ',avRgsqr 
    write(un_sys,*)'avRendsqr   = ',avRendsqr 
    do t=1,nsegtypes
        write(un_sys,*)'qpol(',t,')      = ',qpol(t)
    enddo
    write(un_sys,*)'qpoltot     = ',qpol_tot

    if(systype=="brushdna".or.systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv"&
        .or.systype=="brushborn")then
        do k=1,8
            write(un_sys,*)'avfdisA(',k,')   = ',avfdisA(k)
        enddo
        do t=1,nsegtypes
            write(un_sys,*)'avfdis(',t,')    = ',avfdis(t)
        enddo    
    else
        do t=1,nsegtypes
            write(un_sys,*)'avfdis(',t,')    = ',avfdis(t)
        enddo
    endif
    if(systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv")then
        do k=1,4
            do t=1,nsegtypes
                write(un_sys,'(A8,I5,A,I5,A5,ES25.16)')'avgdisA(',t,',',k,')  = ',avgdisA(t,k)
            enddo
        enddo 
        do k=1,3 
            do t=1,nsegtypes
                write(un_sys,'(A8,I5,A,I5,A5,ES25.16)')'avgdisB(',t,',',k,')  = ',avgdisB(t,k)
            enddo
        enddo    
    endif
    

    if(bcflag(LEFT)=='ta') then
      do i=1,4
        write(un_sys,fmt)'fdisTaL(',i,')  = ',fdisTaL(i)
      enddo
    endif
    if(bcflag(RIGHT)=='ta') then
      do i=1,4
        write(un_sys,fmt)' fdisTaR(',i,')  = ',fdisTaR(i)
      enddo
    endif
    if((bcflag(RIGHT)/='ta').and.(bcflag(RIGHT)/='cc') ) then
      do i=1,6
        write(un_sys,fmt)' fdisSuR(',i,')  = ',fdisS(i)
      enddo
    endif
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'cuantas     = ',cuantas
    write(un_sys,*)'iterations  = ',iter
    write(un_sys,*)'maxniter    = ',maxniter
    write(un_sys,*)'maxfkfunevals = ', maxfkfunevals
    write(un_sys,*)'pH%val      = ',pH%val
    write(un_sys,*)'VdWscale%val= ',VdWscale%val

    ! output ion_excces

    write(un_sys,*)'gamma%Na        = ',ion_excess%Na
    write(un_sys,*)'gamma%Cl        = ',ion_excess%Cl
    write(un_sys,*)'gamma%K         = ',ion_excess%K
    write(un_sys,*)'gamma%Ca        = ',ion_excess%Ca
    write(un_sys,*)'gamma%Mg        = ',ion_excess%Mg
    write(un_sys,*)'gamma%Hplus     = ',ion_excess%Hplus
    write(un_sys,*)'gamma%OHmin     = ',ion_excess%OHmin
    write(un_sys,*)'sumgamma        = ',sum_ion_excess


    ! .. closing files

    close(un_sys)
    close(un_xsol)
    close(un_psi)
    close(un_xpol)

    if(systype/="nucl_ionbin")close(un_fdis)
    if(systype=="brushdna".or.systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv") close(un_fdisP)
    if(systype=="nucl_ionbin".or.systype=="nucl_ionbin_sv")close(un_fdision)
    if(nnucl>=3) close(un_angle)
    if(nnucl>=2) close(un_dist)

    if(verboseflag=="yes") then
        close(un_xNa)
        close(un_xK)
        close(un_xCa)
        close(un_xMg)
        close(un_xNaCl)
        close(un_xKCl)
        close(un_xpair)
        close(un_xCl)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
    endif

end subroutine output_nucl_mul


subroutine output_elect

    !     .. variables and constant declaractions
    use globals
    use volume
    use parameters
    use field
    use energy
    use surface
    use myutils, only : newunit
    use chains, only : isHomopolymer, avRgsqr, avRendsqr,avbond_angle,avdihedral_angle,avnucl_spacing 

    !     .. local arguments

    !     .. output file names

    character(len=90) :: sysfilename
    character(len=90) :: xsolfilename
    character(len=90) :: xpolfilename
    character(len=90) :: xNafilename
    character(len=90) :: xKfilename
    character(len=90) :: xCafilename
    character(len=90) :: xMgfilename
    character(len=90) :: xNaClfilename
    character(len=90) :: xKClfilename
    character(len=90) :: xClfilename
    character(len=90) :: potentialfilename
    character(len=90) :: chargefilename
    character(len=90) :: xHplusfilename
    character(len=90) :: xOHminfilename
    character(len=90) :: densfracAfilename
    character(len=90) :: densfracBfilename
    character(len=90) :: densfracionpairfilename
    character(len=90) :: anglesfilename
    character(len=90) :: spacingfilename
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    integer :: i,j,k,t          ! dummy indexes
    real(dp) :: denspol

    ! .. executable statements

    denspol=init_denspol()

    !     .. make label filename
    call make_filename_label(fnamelabel)

    fnamelabel=trim(fnamelabel)//".dat"

    sysfilename='system.'//trim(fnamelabel)
    xpolfilename='xpol.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    xNafilename='xNaions.'//trim(fnamelabel)
    xKfilename='xKions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xMgfilename='xMgions.'//trim(fnamelabel)
    xNaClfilename='xNaClionpair.'//trim(fnamelabel)
    xKClfilename='xKClionpair.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    potentialfilename='potential.'//trim(fnamelabel)
    chargefilename='charge.'//trim(fnamelabel)
    xHplusfilename='xHplus.'//trim(fnamelabel)
    xOHminfilename='xOHmin.'//trim(fnamelabel)
    densfracAfilename='densityAfrac.'//trim(fnamelabel)
    densfracBfilename='densityBfrac.'//trim(fnamelabel)
    densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)
    anglesfilename='angles.'//trim(fnamelabel)
    spacingfilename='spacing.'//trim(fnamelabel)

    !     .. opening files

    open(unit=newunit(un_sys),file=sysfilename)
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_psi),file=potentialfilename)
    open(unit=newunit(un_xpol),file=xpolfilename)
    open(unit=newunit(un_fdisA),file=densfracAfilename)
    open(unit=newunit(un_fdisB),file=densfracBfilename)
    open(unit=newunit(un_angle),file=anglesfilename)
    open(unit=newunit(un_dist),file=spacingfilename)

    if(verboseflag=="yes") then
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xMg),file=xMgfilename)
        open(unit=newunit(un_xNaCl),file=xNaClfilename)
        open(unit=newunit(un_xKCl),file=xKClfilename)
        open(unit=newunit(un_xpair),file=densfracionpairfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
    endif

    !   .. writting files


    !  .. output of bond and dihedral 
    
    do i=1,nnucl-3
        write(un_angle,*)avbond_angle(i),avdihedral_angle(i)
    enddo
    if(nnucl>=3)write(un_angle,*)avbond_angle(nnucl-2)

        
    do i=1,nnucl-1
        write(un_dist,*)avnucl_spacing(i)
    enddo

    do i=1,nsize
        write(un_xsol,*)xsol(i)
        write(un_psi,*)psi(i)
    enddo

    do i=1,nsize
        write(un_xpol,fmt3reals)xpol(i),rhopol(i,1),rhopol(i,2)
        write(un_fdisA,fmt5reals)(fdisA(i,k),k=1,5)
        write(un_fdisB,fmt5reals)(fdisB(i,k),k=1,5)
    enddo

    if(verboseflag=="yes") then
        do i=1,nsize
            write(un_xNa,*)xNa(i)
            write(un_xK,*)xK(i)
            write(un_xCa,*)xCa(i)
            write(un_xMg,*)xMg(i)
            write(un_xNaCl,*)xNaCl(i)
            write(un_xKCl,*)xKCl(i)
            write(un_xpair,*)(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)xCl(i)
            write(un_charge,*)rhoq(i)
            write(un_xHplus,*)xHplus(i)
            write(un_xOHmin,*)xOHmin(i)
        enddo
    endif

    ! .. writing system information
   
    write(un_sys,*)'system      =  nucleosome chain'
    write(un_sys,*)'version     = ',VERSION
    ! chain description
    write(un_sys,*)'chainmethod = ',chainmethod
    write(un_sys,*)'chaintype   = ',chaintype
    write(un_sys,*)"isHomopolymer= ",isHomopolymer
    if(chainmethod.eq."FILE") then
       write(un_sys,*)'readinchains = ',readinchains
    endif
    write(un_sys,*)'nseg        = ',nseg
    write(un_sys,*)'nnucl       = ',nnucl
    write(un_sys,*)'chainperiod = ',chainperiod
    write(un_sys,*)'cuantas     = ',cuantas
    write(un_sys,*)'denspol     = ',denspol

    ! system description
    write(un_sys,*)'systype     = ',systype
    write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
    write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
    write(un_sys,*)'delta       = ',delta
    write(un_sys,*)'nx          = ',nx
    write(un_sys,*)'ny          = ',ny
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'tol_conv    = ',tol_conv

    ! concentration
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    write(un_sys,*)'cMgCl2      = ',cMgCl2
    write(un_sys,*)'xNabulk     = ',xbulk%Na
    write(un_sys,*)'xClbulk     = ',xbulk%Cl
    write(un_sys,*)'xKbulk      = ',xbulk%K
    write(un_sys,*)'xNaClbulk   = ',xbulk%NaCl
    write(un_sys,*)'xKClbulk    = ',xbulk%KCl
    write(un_sys,*)'xCabulk     = ',xbulk%Ca
    write(un_sys,*)'xHplusbulk  = ',xbulk%Hplus
    write(un_sys,*)'xOHminbulk  = ',xbulk%OHmin
    write(un_sys,*)'pHbulk      = ',pHbulk
    ! disociation constants
    write(un_sys,*)'pKa         = ',pKaA(1)
    write(un_sys,*)'pKaNa       = ',pKaA(2)
    write(un_sys,*)'pKaACa      = ',pKaA(3)
    write(un_sys,*)'pKaA2Ca     = ',pKaA(4)
    write(un_sys,*)'pKb         = ',pKaB(1)
    write(un_sys,*)'pKbNa       = ',pKaB(2)
    write(un_sys,*)'pKbBCa      = ',pKaB(3)
    write(un_sys,*)'pKbB2Ca     = ',pKaB(4)
    write(un_sys,*)'KionNa      = ',KionNa
    write(un_sys,*)'KionK       = ',KionK
    write(un_sys,*)'K0ionNa     = ',K0ionNa
    write(un_sys,*)'K0ionK      = ',K0ionK
    write(un_sys,*)'dielectW    = ',dielectW 
    if(runtype=="rangedielect") then 
        write(un_sys,*)'dielectscale = ',dielectscale%val
    endif  
    write(un_sys,*)'lb          = ',lb
    write(un_sys,*)'T           = ',Tref
    if(isVdW) then
        write(un_sys,*)'VdWeps  = ',VdWeps
        write(un_sys,*)'VdWepsAA = ',VdWepsAA
        write(un_sys,*)'VdWepsBB = ',VdWepsBB
        write(un_sys,*)'VdWepsAB = ',VdWepsAB
    endif
    ! charge components
    write(un_sys,*)'zpolA(1)    = ',zpolA(1)
    write(un_sys,*)'zpolA(2)    = ',zpolA(2)
    write(un_sys,*)'zpolA(3)    = ',zpolA(3)
    write(un_sys,*)'zpolA(4)    = ',zpolA(4)
    write(un_sys,*)'zpolB(1)    = ',zpolB(1)
    write(un_sys,*)'zpolB(2)    = ',zpolB(2)
    write(un_sys,*)'zpolB(3)    = ',zpolB(3)
    write(un_sys,*)'zpolB(4)    = ',zpolB(4)
    write(un_sys,*)'zpolB(5)    = ',zpolB(5)
    write(un_sys,*)'zNa         = ',zNa
    write(un_sys,*)'zCa         = ',zCa
    write(un_sys,*)'zK          = ',zK
    write(un_sys,*)'zCl         = ',zCl

     ! volume
    write(un_sys,*)'vsol        = ',vsol
    write(un_sys,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(un_sys,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(un_sys,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(un_sys,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(un_sys,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(un_sys,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(un_sys,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(un_sys,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(un_sys,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(un_sys,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(un_sys,*)'vNa         = ',vNa*vsol
    write(un_sys,*)'vCl         = ',vCl*vsol
    write(un_sys,*)'vCa         = ',vCa*vsol
    write(un_sys,*)'vK          = ',vK*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol
    do t=1,nsegtypes
        write(un_sys,*)'vpol(',t,')       = ',vpol(t)*vsol
    enddo   
    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'tol_conv    = ',tol_conv
    write(un_sys,*)'denspol     = ',denspol
    write(un_sys,*)'sumphiA     = ',sumphiA
    write(un_sys,*)'sumphiB     = ',sumphiB
    write(un_sys,*)'check phi   = ',checkphi
    write(un_sys,*)'FEq         = ',FEq
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEel        = ',FEel
    write(un_sys,*)'FEelsurf(LEFT)  = ',FEelsurf(LEFT)
    write(un_sys,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
    write(un_sys,*)'FEbind      = ',FEbind
    write(un_sys,*)'FEVdW       = ',FEVdW
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'q           = ',q
    write(un_sys,*)'avRgsqr     = ',avRgsqr 
    write(un_sys,*)'avRendsqr   = ',avRendsqr 
    write(un_sys,*)'qpolA       = ',qpolA
    write(un_sys,*)'qpolB       = ',qpolB
    write(un_sys,*)'qpoltot     = ',qpol_tot
    write(un_sys,*)'avfdisA(1)  = ',avfdisA(1)
    write(un_sys,*)'avfdisA(2)  = ',avfdisA(2)
    write(un_sys,*)'avfdisA(3)  = ',avfdisA(3)
    write(un_sys,*)'avfdisA(4)  = ',avfdisA(4)
    write(un_sys,*)'avfdisA(5)  = ',avfdisA(5)
    write(un_sys,*)'avfdisB(1)  = ',avfdisB(1)
    write(un_sys,*)'avfdisB(2)  = ',avfdisB(2)
    write(un_sys,*)'avfdisB(3)  = ',avfdisB(3)
    write(un_sys,*)'avfdisB(4)  = ',avfdisB(4)
    write(un_sys,*)'avfdisB(5)  = ',avfdisB(5)
    write(un_sys,*)'sigmaSurfL  = ',sigmaSurfL/((4.0_dp*pi*lb)*delta)
    write(un_sys,*)'sigmaSurfR  = ',sigmaSurfR/((4.0_dp*pi*lb)*delta)

    if(bcflag(LEFT)=='ta') then
      do i=1,4
        write(un_sys,fmt)'fdisTaL(',i,')  = ',fdisTaL(i)
      enddo
    endif
    if(bcflag(RIGHT)=='ta') then
      do i=1,4
        write(un_sys,fmt)' fdisTaR(',i,')  = ',fdisTaR(i)
      enddo
    endif
    if((bcflag(RIGHT)/='ta').and.(bcflag(RIGHT)/='cc') ) then
      do i=1,6
        write(un_sys,fmt)' fdisSuR(',i,')  = ',fdisS(i)
      enddo
    endif
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'cuantas     = ',cuantas
    write(un_sys,*)'iterations  = ',iter
    write(un_sys,*)'maxniter    = ',maxniter
    write(un_sys,*)'maxfkfunevals = ', maxfkfunevals
    write(un_sys,*)'pH%val      = ',pH%val
    write(un_sys,*)'VdWscale%val= ',VdWscale%val

    ! output ion_exces

    write(un_sys,*)'gamma%Na        = ',ion_excess%Na
    write(un_sys,*)'gamma%Cl        = ',ion_excess%Cl
    write(un_sys,*)'gamma%K         = ',ion_excess%K
    write(un_sys,*)'gamma%Ca        = ',ion_excess%Ca
    write(un_sys,*)'gamma%Mg        = ',ion_excess%Mg
    write(un_sys,*)'gamma%Hplus     = ',ion_excess%Hplus
    write(un_sys,*)'gamma%OHmin     = ',ion_excess%OHmin

    ! .. closing files
    
    close(un_sys)
    close(un_xsol)
    close(un_psi)
    close(un_xpol)
    close(un_fdisA)
    close(un_fdisB)
    close(un_angle)
    close(un_dist)

    if(verboseflag=="yes") then
        close(un_xNa)
        close(un_xK)
        close(un_xCa)
        close(un_xMg)
        close(un_xNaCl)
        close(un_xKCl)
        close(un_xpair)
        close(un_xCl)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
    endif

end subroutine output_elect


subroutine output_neutral

    !     .. variables and constant declaractions
    use globals
    use volume
    use parameters
    use field
    use energy
    use myutils, only : newunit
    use chains, only : isHomopolymer, avRgsqr, avRendsqr, avbond_angle,avdihedral_angle,avnucl_spacing 

    !     .. output file names
    character(len=90) :: sysfilename
    character(len=90) :: xsolfilename
    character(len=90) :: xpolfilename
    character(len=90) :: anglesfilename
    character(len=90) :: spacingfilename
    
    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals,fmtNplus1reals

    !     .. local arguments
    integer :: i, t
    character(len=100) :: fnamelabel
    character(len=20) :: rstr,istr
    logical :: isopen
    real(dp) :: denspol

    !     .. executable statements

    denspol=init_denspol()

    fmt2reals = "(2ES25.16)"
    fmt3reals = "(3ES25.16)"
    fmt4reals = "(4ES25.16)"
    fmt5reals = "(5ES25.16)"
    fmt6reals = "(6ES25.16)"


    write(istr,'(I2)')nsegtypes+1
    fmtNplus1reals='('//trim(adjustl(istr))//'ES25.16)' ! format nsegtypes+! reals

    !     .. make label filename
    call make_filename_label(fnamelabel)

    !     .. make filenames
    sysfilename='system.'//trim(fnamelabel)
    xpolfilename='xpol.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    anglesfilename='angles.'//trim(fnamelabel) 
    spacingfilename='spacing.'//trim(fnamelabel) 


    !      .. opening files
    open(unit=newunit(un_sys),file=sysfilename)
    open(unit=newunit(un_xpol),file=xpolfilename)
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_angle),file=anglesfilename)
    open(unit=newunit(un_dist),file=spacingfilename)

    !  .. writting files

    do i=1,nsize
       write(un_xpol,fmtNplus1reals)xpol(i),(rhopol(i,t),t=1,nsegtypes)
       write(un_xsol,fmt1reals)xsol(i)
    enddo

    !  .. output of bond and dihedral angles
   
    do i=1,nnucl-3
        write(un_angle,*)avbond_angle(i),avdihedral_angle(i)
    enddo
    if(nnucl>=3)write(un_angle,*)avbond_angle(nnucl-2)
    
    !  .. output of nuclesome distance 
    do i=1,nnucl-1
        write(un_dist,*)avnucl_spacing(i)
    enddo

    !     .. system information

    write(un_sys,*)'system      = nuclesome chain'
    write(un_sys,*)'version     = ',VERSION
    ! chain description
    write(un_sys,*)'chainmethod = ',chainmethod
    write(un_sys,*)'chaintype   = ',chaintype
    write(un_sys,*)"isHomopolymer= ",isHomopolymer
    if(chainmethod.eq."FILE") then
       write(un_sys,*)'readinchains = ',readinchains
    endif
    write(un_sys,*)'nseg        = ',nseg
    write(un_sys,*)'nnucl       = ',nnucl
    write(un_sys,*)'chainperiod = ',chainperiod
    write(un_sys,*)'cuantas     = ',cuantas
    write(un_sys,*)'maxnchainsrot      = ',maxnchainsrotations
    write(un_sys,*)'maxnchainsrotxy    = ',maxnchainsrotationsxy
    ! system description
    write(un_sys,*)'systype     = ',systype
    write(un_sys,*)'delta       = ',delta
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'nx          = ',nx
    write(un_sys,*)'ny          = ',ny
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'tol_conv    = ',tol_conv
    ! other physical parameters
    write(un_sys,*)'T           = ',Tref
     do t=1,nsegtypes
        write(un_sys,*)'lseg(',t,')        = ',lseg
    enddo 
    ! volume
    write(un_sys,*)'vsol        = ',vsol
    do t=1,nsegtypes
        write(un_sys,*)'vpol(',t,')        = ',vpol(t)*vsol
    enddo 
    write(un_sys,*)'isVdW       = ',isVdW

    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'FEconf      = ',FEconf
    write(un_sys,*)'Econf       = ',Econf
    write(un_sys,*)'FEtrans%sol = ',FEtrans%sol
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'sumphi      = ',(sumphi(t),t=1,nsegtypes)
    write(un_sys,*)'check phi   = ',checkphi
    write(un_sys,*)'denspol     = ',denspol
    write(un_sys,*)'FEq         = ',FEq
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEVdW       = ',FEVdW
    write(un_sys,*)'q           = ',q
    write(un_sys,*)'mu          = ',-log(q)
    write(un_sys,*)'avRgsqr     = ',avRgsqr 
    write(un_sys,*)'avRendsqr   = ',avRendsqr 
    write(un_sys,*)'iterations  = ',iter
    write(un_sys,*)'VdWscale%val = ',VdWscale%val

    ! .. closing files
    
    close(un_xsol)
    close(un_xpol)
    close(un_sys)
    close(un_angle)
    close(un_dist)
    
end subroutine output_neutral


subroutine output_individualcontr_fe

    use globals, only : LEFT,RIGHT, systype
    use energy
    use myutils, only : newunit
    use volume, only : delta,nz
    use parameters, only : isEnergyShift
    use chains, only : energychain_min

    ! local arguments

    character(len=100) :: fenergyfilename
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

   
    !     .. make label filename
    call make_filename_label(fnamelabel)
    fenergyfilename='energy.'//trim(fnamelabel)
    !     .. opening files
    open(unit=newunit(un_fe),file=fenergyfilename)
    
    write(un_fe,*)'FE              = ',FE
    write(un_fe,*)'FEbulk          = ',FEbulk
    write(un_fe,*)'deltaFE         = ',deltaFE
    write(un_fe,*)'FEalt           = ',FEalt
    write(un_fe,*)'FEbulkalt       = ',FEbulkalt
    write(un_fe,*)'deltaFEalt      = ',deltaFEalt
    write(un_fe,*)'deltadeltaFE    = ',deltaFE-deltaFEalt

    write(un_fe,*)'FEq             = ',FEq
    write(un_fe,*)'FEconf          = ',FEconf
    write(un_fe,*)'Econf           = ',Econf
    write(un_fe,*)'Eshift          = ',Eshift
    write(un_fe,*)'isEnergyShift   = ',isEnergyShift
    write(un_fe,*)'Emin            = ',energychain_min

    write(un_fe,*)'FEpi            = ',FEpi
    write(un_fe,*)'FErho           = ',FErho
    write(un_fe,*)'FEel            = ',FEel
    write(un_fe,*)'FEelsurf(LEFT)  = ',FEelsurf(LEFT)
    write(un_fe,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
    write(un_fe,*)'FEelvar         = ',FEelvar
    write(un_fe,*)'FEelvarborn     = ',FEelvarborn
    write(un_fe,*)'FEborn          = ',FEborn
    write(un_fe,*)'FEbornbulk      = ',FEbornbulk
    write(un_fe,*)'deltaFEborn     = ',FEborn - FEbornbulk

    write(un_fe,*)'FEbind          = ',FEbind
    write(un_fe,*)'FEchem          = ',FEchem
    write(un_fe,*)'FEVdW           = ',FEVdW

    write(un_fe,*)"FEtrans%sol     = ",FEtrans%sol
    write(un_fe,*)"FEtrans%Na      = ",FEtrans%Na
    write(un_fe,*)"FEtrans%Cl      = ",FEtrans%Cl
    write(un_fe,*)"FEtrans%Ca      = ",FEtrans%Ca
    write(un_fe,*)"FEtrans%Mg      = ",FEtrans%Mg
    write(un_fe,*)"FEtrans%K       = ",FEtrans%K
    write(un_fe,*)"FEtrans%KCl     = ",FEtrans%KCl
    write(un_fe,*)"FEtrans%NaCl    = ",FEtrans%NaCl
    write(un_fe,*)"FEtrans%Hplus   = ",FEtrans%Hplus
    write(un_fe,*)"FEtrans%OHmin   = ",FEtrans%OHmin

    write(un_fe,*)"FEchempot%Na    = ",FEchempot%Na
    write(un_fe,*)"FEchempot%Cl    = ",FEchempot%Cl
    write(un_fe,*)"FEchempot%Ca    = ",FEchempot%Ca
    write(un_fe,*)"FEchempot%Mg    = ",FEchempot%Mg
    write(un_fe,*)"FEchempot%K     = ",FEchempot%K
    write(un_fe,*)"FEchempot%KCl   = ",FEchempot%KCl
    write(un_fe,*)"FEchempot%NaCl  = ",FEchempot%NaCl
    write(un_fe,*)"FEchempot%Hplus = ",FEchempot%Hplus
    write(un_fe,*)"FEchempot%OHmin = ",FEchempot%OHmin

    write(un_fe,*)"FEchemsurf(LEFT)     = ",FEchemsurf(LEFT)
    write(un_fe,*)"FEchemsurf(RIGHT)    = ",FEchemsurf(RIGHT)
    write(un_fe,*)"FEchemsurfalt(LEFT)  = ",FEchemsurfalt(LEFT)
    write(un_fe,*)"FEchemsurfalt(RIGHT) = ",FEchemsurfalt(RIGHT)

    write(un_fe,*)"delta FEchemsurfalt(LEFT) = ",FEchemsurfalt(LEFT)-FEchemsurf(LEFT)-diffFEchemsurf(LEFT)
    write(un_fe,*)"delta FEchemsurfalt(RIGHT)= ",FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)-diffFEchemsurf(RIGHT)


    close(un_fe)

end subroutine output_individualcontr_fe


subroutine make_filename_label(fnamelabel)

    use globals, only : LEFT,RIGHT, systype, runtype
    use parameters, only : cNaCl,cKCl,cCaCl2,cMgCl2,pHbulk,VdWepsBB,init_denspol,VdWscale,pKd,dielectscale

    character(len=*), intent(inout) :: fnamelabel

    character(len=20) :: rstr
    real(dp) :: denspol


    denspol=init_denspol()

    !     .. make label filename

    select case(systype)

    case("elect","electnopoly","electA")

        write(rstr,'(F5.3)')denspol
        fnamelabel="phi"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))

        if(cKCl/=0.0_dp) then
            if(cKCl>=0.001_dp) then
                write(rstr,'(F5.3)')cKCl
            else
                write(rstr,'(ES9.2E2)')cKCl
            endif
            fnamelabel=trim(fnamelabel)//"cKCl"//trim(adjustl(rstr))
        endif

        if(cCaCl2>=0.001) then
            write(rstr,'(F5.3)')cCaCl2
        elseif(cCaCl2>0.0) then
            write(rstr,'(ES9.2E2)')cCaCl2
        else
            write(rstr,'(F3.1)')cCaCl2
        endif
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"

    case("neutral","neutralnoVdW","nucl_neutral_sv")

        if(denspol>=0.001) then
            write(rstr,'(F5.3)')denspol
        else
            write(rstr,'(ES9.2E2)')denspol
        endif

        fnamelabel="phi"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')VdWscale%val
        fnamelabel=trim(fnamelabel)//"VdWscale"//trim(adjustl(rstr))//".dat"

    case("brush_mul","brush_mulnoVdW","brushdna","nucl_ionbin","nucl_ionbin_sv",&
        "brushborn")

        write(rstr,'(F5.3)')denspol
        fnamelabel="phi"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))

        if(cKCl/=0.0_dp) then
            if(cKCl>=0.001_dp) then
                write(rstr,'(F5.3)')cKCl
            else
                write(rstr,'(ES9.2E2)')cKCl
            endif
            fnamelabel=trim(fnamelabel)//"cKCl"//trim(adjustl(rstr))
        endif

        if(cCaCl2/=0.0_dp) then
            if(cCaCl2>=0.001) then
                if(cCaCl2>=0.01) then
                    write(rstr,'(F5.3)')cCaCl2
                else
                    write(rstr,'(F6.4)')cCaCl2
                endif  
            elseif(cCaCl2>0.0) then
                write(rstr,'(ES9.2E2)')cCaCl2
            else
                write(rstr,'(F3.1)')cCaCl2
            endif
            fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        endif
        
        if(cMgCl2/=0.0_dp) then
            if(cMgCl2>=0.001) then
                if(cMgCl2>=0.01) then
                    write(rstr,'(F5.3)')cMgCl2
                else
                    write(rstr,'(F6.4)')cMgCl2
                endif  
            elseif(cMgCl2>0.0) then
                write(rstr,'(ES9.2E2)')cMgCl2
            else
                write(rstr,'(F3.1)')cMgCl2
            endif

            fnamelabel=trim(fnamelabel)//"cMgCl2"//trim(adjustl(rstr))
        endif

        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))

        if(runtype=="rangepKd") then
            if(pKd%val<0) then
                write(rstr,'(F5.2)')pKd%val
            else
                write(rstr,'(F5.3)')pKd%val
            endif
            fnamelabel=trim(fnamelabel)//"pKd"//trim(adjustl(rstr))//".dat"

        elseif(runtype=="rangeVdWeps") then
            write(rstr,'(F5.3)')VdWscale%val
            fnamelabel=trim(fnamelabel)//"VdWscale"//trim(adjustl(rstr))//".dat"

        elseif(runtype=="rangedielect") then
            write(rstr,'(ES9.2E2)')dielectscale%val
            fnamelabel=trim(fnamelabel)//"dielectscale"//trim(adjustl(rstr))//".dat"

        else
            fnamelabel=trim(fnamelabel)//".dat"
        endif


    case default
        print*,"Error in make_filename_label subroutine"
        print*,"Wrong value systype : ", systype
    endselect

end subroutine make_filename_label

subroutine copy_solution(x)

    use globals, only : systype, neq, nsize, bcflag, LEFT, RIGHT
    use volume, only  : nx,ny
    use surface, only : psiSurfL, psiSurfR
    use field

    real(dp), dimension(neq) :: x  ! expliciet size array

    ! local variable
    integer :: i, neq_bc
    integer, parameter :: A=1, B=2

    select case (systype)
    case ("elect")

        do i=1,nsize
            xsol(i)= x(i)
            psi(i) = x(i+nsize)
            rhopol(i,A)=x(i+2*nsize)
            rhopol(i,B)=x(i+3*nsize)
        enddo

        neq_bc=0 ! surface potential
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nx*ny
            do i=1,neq_bc
                psiSurfR(i) =x(4*nsize+i)
            enddo
        endif
        if(bcflag(LEFT)/="cc") then
            do i=1,nx*ny
                psiSurfL(i) =x(4*nsize+neq_bc+i)
            enddo
            neq_bc=neq_bc+nx*ny
        endif

    case ("neutral")

        do i=1,nsize
            xsol(i)= x(i)
        enddo

    case ("neutralnoVdW")

        do i=1,nsize
            xsol(i)= x(i)
        enddo

    case ("brush_mulnoVdW")

        do i=1,nsize
            xsol(i)= x(i)
            psi(i) = x(i+nsize)
        enddo

    case default

        print*,"Error: systype incorrect in copy_solution"
        print*,"stopping program"
        stop

    end select

end subroutine copy_solution

! output routine

subroutine compute_vars_and_output()

    use globals, only : systype,DEBUG
    use energy, only : fcnenergy,check_volume_xpol
    use field, only : charge_polymer,average_charge_polymer,make_ion_excess

    select case (systype)
    case ("elect")

        call fcnenergy()
        call charge_polymer()
        call average_charge_polymer()
        call make_ion_excess()
        call output()

    case ("neutral","neutralnoVdW")

        call fcnenergy()
        call output()           

    case ("brush_mul","brush_mulnoVdW","brushdna",&
        "nucl_ionbin","brushborn")

        call fcnenergy()
        call charge_polymer()
        call average_charge_polymer()
        call make_ion_excess()
        call output()  
        if(DEBUG) call check_volume_xpol()        

    case ("nucl_ionbin_sv")

        call fcnenergy()
        call charge_polymer()
        call average_charge_polymer()
        call make_ion_excess()
        call output()           
        if(DEBUG) call check_volume_xpol() 

     case ("nucl_neutral_sv")

        call fcnenergy()
        call output()           

    case default

        print*,"Error: systype incorrect in compute_vars_and_output"
        print*,"stopping program"
        stop

    end select

end subroutine compute_vars_and_output


subroutine write_chain_config()

    use mpivars, only : rank
    use globals, only: nseg, nsegtypes
    use chains, only : type_of_monomer,type_of_monomer_char,isAmonomer,type_of_charge
    use chains, only : mapping_num_to_char
    use parameters, only : lseg,lsegAA,vpol, vsol
    use myutils, only : newunit

    character(len=100) :: fname
    integer :: i, un_cc
    character(len=10) ::istr
   

    if(rank==0) then
        
        write(istr,'(I4)')rank
        fname='chain_config.'//trim(adjustl(istr))//'.log'
        !     .. opening file
        open(unit=newunit(un_cc),file=fname)

        write(un_cc,*)"chain configuration summary"
        write(un_cc,*)"#nsegtypes char charge type "
        do i=1,nsegtypes
            write(un_cc,*)i,mapping_num_to_char(i),' ',type_of_charge(i)
        enddo
        write(un_cc,*)"###"
        write(un_cc,*)"lseg=",lseg
        write(un_cc,*)"#nsegtypes   lsegAA    vpol"
        do i=1,nsegtypes
            write(un_cc,*)i,lsegAA(i),vpol(i)*vsol
        enddo
        write(un_cc,*)"#segment type_num  type_char  isAmonomer"
        do i=1,nseg
            write(un_cc,*)i,type_of_monomer(i),type_of_monomer_char(i),isAmonomer(i)
        enddo
        write(un_cc,*)"###"
        close(un_cc)

    endif


end subroutine write_chain_config


end module
