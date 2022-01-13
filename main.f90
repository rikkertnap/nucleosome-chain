! ---------------------------------------------------------------|
! Solves the SCMFT eqs for WEAK polyelectrolytes polymers        |
! coated onto a planar surface,                                  |
! input/output: see myio.f90                                     |
! ---------------------------------------------------------------|


program main

    !     .. variable and constant declaractions
    use mpivars
    use globals       ! parameters definitions
    use physconst
    use mathconst
    use volume
    use random
    use field
    use parameters
    use matrices
    use energy
    use conform_entropy
    use chains
    use chaingenerator
    use VdW
    use listfcn
    use fcnpointer
    use initxvector
    use surface
    use myio
    use myutils
    use dielectric_const

    implicit none

    real(dp),  dimension(:), allocatable :: x         ! iteration vector
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector
    real(dp),  dimension(:), allocatable :: fvec

    integer :: i,c,num
    logical :: use_xstored
    logical :: isfirstguess
    logical :: issolution
    integer :: info
    character(len=lenText) :: text, istr, rstr
    character(len=20) :: fname, conffilename
    integer :: iend , un_conf
    type (looplist), pointer :: loop
    real(dp) :: loopbegin,  loopstepsizebegin
    real(dp), parameter :: loopeps = 1.0e-10_dp   
    real(dp), dimension(:),  pointer :: list
    real(dp), pointer :: list_val
    real(dp) :: list_first


    ! .. executable statements


    ! .. mpi

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

    ! .. logfile
    
    write(istr,'(I4)')rank
    if( size>9999) then 
        text="Error: size to large for status file number"
        call print_to_log(LogUnit,text)
        print*,text
        call MPI_FINALIZE(ierr)
        stop
    endif
    fname='status.'//trim(adjustl(istr))//'.log'
    call open_logfile(logUnit,fname)
    write(istr,'(I4)')rank
    
    text='program begins : rank '//istr
    call print_to_log(LogUnit,text)
    write(istr,'(A40)')VERSION
    text='program version     = '//istr
    call print_to_log(LogUnit,text)
    if(rank==0) print*,text

    ! .. init

    call read_inputfile(info)
    if(info/=0) then
        write(istr,'(I3)')info
        text="Error in input file: info = "//istr//" : end program."
        call print_to_log(LogUnit,text)
        print*,text
        call MPI_FINALIZE(ierr)
        stop
    endif

    call init_constants()
    call make_geometry()            ! generate volume elements lattice
    call allocate_chain_parameters()  
    call init_matrices()            ! init matrices for chain generation
    call allocate_chains(cuantas,nnucl,nseg,nsegtypes,maxnchainsrotations,maxnchainsrotationsxy)
    call init_chain_parameters      ! chain volume, charges, pKa etc
    call make_sequence_chain(chainperiod,chaintype)
    call make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
    call make_segcom(segcm,nnucl,segcmfname)
    call set_properties_chain(chainperiod,chaintype) 

    if(isVdW) then 
        call make_VdWcoeff(info)
        if(info/=0) then
            write(istr,'(I3)')info
            text="Error in make_VdWcoeff: info = "//trim(adjustl(istr))//" : end program."
            call print_to_log(LogUnit,text)
            print*,text
            stop
        endif
    else
        call make_VdWeps(info)    
    endif  

    call make_chains(chainmethod)   
    call chain_filter()
    call allocate_field(nx,ny,nz,nsegtypes)
    call init_field()
    call init_surface(bcflag,nsurf)

    call make_isrhoselfconsistent(isVdW)
    call set_size_neq()             ! number of non-linear equation neq
    call set_fcn()
    call set_dielect_fcn(dielect_env)

    call write_chain_config()
    call write_chain_struct(write_struct,info)

    !  .. computation starts

    allocate(xstored(neq))
    allocate(x(neq))
    allocate(xguess(neq))
    allocate(fvec(neq))

        
    if(runtype=="rangedist") then ! loop over distances

        nz = nzmax
        neqmax = neq

        isfirstguess = .true.
        use_xstored = .false.         ! with both flags set false make_guess will set xguess equal to x
        
        pH%val=pH%min
       
        iter = 0


        do while (nz>=nzmin)        ! loop distances

            call set_size_neq()

            if(.not.allocated(x)) allocate(x(neq))
            if(.not.allocated(xguess)) allocate(xguess(neq))
            if(.not.allocated(fvec)) allocate(fvec(neq))

            call init_vars_input()          ! sets up chem potenitals
            call chain_filter()
            call set_fcn()           

            flag_solver = 0

            if(rank.eq.0) then     ! node rank=0
                call make_guess(x, xguess, isfirstguess, use_xstored, xstored)
                call solver(x, xguess, tol_conv, fnorm, issolution)
                call fcnptr(x, fvec, neq)
                flag_solver = 0   ! stop nodes
                do i = 1, size-1
                    dest =i
                    call MPI_SEND(flag_solver, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
                enddo
            else
                flag_solver = 1
                do while(flag_solver.eq.1)
                    flag_solver = 0
                    source = 0
                    call MPI_RECV(flag_solver, 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
                    if(flag_solver.eq.1) then
                        call MPI_RECV(x, neqint, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)
                        call fcnptr(x, fvec, neq)
                    endif
                enddo
            endif

            call FEconf_entropy(FEconf,Econf) ! parrallel computation of conf entropy

            if(rank==0) then

                call compute_vars_and_output()

                isfirstguess =.false.
                use_xstored = .true.
                iter = 0                ! reset of iteration counter
                nz = nz-nzstep          ! reduce distance
                do i=1,neq
                    xstored(i)=x(i)
                enddo
                ! communicate new values of nz from master to  compute  nodes to advance while loop on compute nodes
                do i = 1, size-1
                    dest = i
                    call MPI_SEND(nz, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
                enddo
            else ! receive values
                source = 0
                call MPI_RECV(nz , 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
            endif

            iter  = 0              ! reset of iteration counter
            deallocate(x)
            deallocate(xguess)
            deallocate(fvec)

        enddo ! end while loop


    else  ! loop over pH, or pKd etc  values
 

        if(runtype=="inputcspH".or.runtype=="inputMgpH".or.runtype=="inputcsKClpH") then 
            loop => pH
        else if(runtype=="rangecpro") then 
            loop => cpro
        else if (runtype=="rangepKd") then
            loop => pKd   
        else if(runtype=="rangeVdWeps") then 
            loop => VdWscale    
        else
            if(associated(loop)) nullify(loop) ! make explict that no association is made
        endif  

         ! .. select variable with which list_array to associate

        if (runtype=="inputMgpH".or.runtype=="rangepKd".or.runtype=="rangeVdWeps") then
            call set_value_MgCl2(runtype,info)
            if(info/=0) then
                print*,"Error in input file: info = ",info," : end program."
                stop
            endif
            num=num_cMgCl2
            list=>cMgCl2_array
            list_val => cMgCl2

        else if(runtype=="inputcsKClpH") then
            call set_value_KCl(runtype,info)
            if(info/=0) then
                print*,"Error in input file: info = ",info," : end program."
                stop
            endif
            num=num_cKCl
            list=>cKCl_array
            list_val => cKCl

        else 
            call set_value_NaCl(runtype,info)
            if(info/=0) then
                print*,"Error in input file: info = ",info," : end program."
                stop
            endif
            num=num_cNaCl
            list=>cNaCl_array
            list_val => cNaCl ! pointer points to target  cNaCl

        endif

        nz = nzmax
        neqmax = neq

        isfirstguess = .true.
        use_xstored = .false.         ! with both flags set false make_guess will set xguess equal to x
        iter = 0                    ! iteration counter

        if(loop%stepsize>0) then
            loop%val=loop%min
        else
            loop%val=loop%max
        endif

        call set_fcn()
        call chain_filter() 
         
        ! free unused variables 
        deallocate(energychain)
        deallocate(energychain_init)
        deallocate(indexchain_init) 

        list_first= list(1)
        loopstepsizebegin=loop%stepsize


        do c=1,num                          ! loop over list items       
   
            list_val=list(c)                ! get value from array
            iter = 0                        ! iteration counter 
            loop%stepsize=loopstepsizebegin ! reset of loopstep size 
                                  
            if(loop%stepsize>0) then        ! sign loops%stepsize determines direction loop
                loop%val=loop%min
                loopbegin=loop%min 
            else
                loop%val=loop%max
                loopbegin=loop%max
            endif    

            do while (loop%min<=loop%val.and.loop%val<=loop%max.and.&
                    (abs(loop%stepsize)>=loop%delta))
                
                isfirstguess=(loop%val==loopbegin) !abs(loop%val-loopbegin)<loopeps)
               
                !print*,"loop%val+",loop%val,"isfirstguess=",isfirstguess," list_val=", list_val
                !print*,"use_xstored=",use_xstored
                call init_vars_input()  ! sets up chem potentials
                
                flag_solver = 0

                if(rank==0) then     ! node rank=0
                    call make_guess(x, xguess, isfirstguess,use_xstored,xstored)
                    call solver(x, xguess, tol_conv, fnorm, issolution)
                    call fcnptr(x, fvec, neq)
                    flag_solver = 0   ! stop nodes
                    do i = 1, size-1
                        dest =i
                        call MPI_SEND(flag_solver, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
                    enddo
                else
                    flag_solver = 1
                    do while(flag_solver.eq.1)
                        flag_solver = 0
                        source = 0
                        call MPI_RECV(flag_solver, 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
                        if(flag_solver.eq.1) then
                            call MPI_RECV(x, neqint, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)
                            call fcnptr(x, fvec, neq)
                        endif
                    enddo
                endif

                call FEconf_entropy(FEconf,Econf) ! parrallel computation of conf FEconf_entropy

                if(rank==0) then

                    if(.not.isSolution) then
                        
                        text="no solution: backstep" 
                        call print_to_log(LogUnit,text)
                        loop%stepsize=loop%stepsize/2.0d0   ! decrease increment
                        loop%val=loop%val-loop%stepsize     ! step back
                        do i=1,neq
                            x(i)=xguess(i)
                         enddo
                    
                    else 

                        call compute_vars_and_output()
                        if(isfirstguess) then
                           do i=1,neqint
                              xstored(i)=x(i)
                           enddo
                           use_xstored=.false.
                        endif
                        
                        isfirstguess =.false.
                        loop%val=loop%val+loop%stepsize

                     endif
                     
                    
                    ! communicate new values of loop from head to compute nodes to advance while loop on compute nodes
                    do i = 1, size-1
                        dest = i
                        call MPI_SEND(loop%val, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,ierr)
                    enddo
                else ! receive values
                    source = 0
                    call MPI_RECV(loop%val, 1, MPI_DOUBLE_PRECISION, source, tag,MPI_COMM_WORLD,stat, ierr)
                endif

                iter  = 0              ! reset of iteration counter

            enddo ! end while loop

            use_xstored=.true.

        enddo ! end loop list item

        deallocate(x)
        deallocate(xguess)
        deallocate(fvec)

    endif

    call MPI_FINALIZE(ierr)

    deallocate(xstored)
      
    call deallocate_field()

    text="program end"

    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program main
