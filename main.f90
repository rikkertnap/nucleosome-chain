
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
    real(dp), parameter :: listeps = 1.0e-7_dp   
    real(dp), dimension(:),  pointer :: list
    real(dp), pointer :: list_val
    real(dp) :: list_first, list_step
    integer  :: nlist_elem, maxlist_elem, nlist_step


    ! .. executable statements

    ! .. mpi

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numproc, ierr)
    size=numproc !! temporary needs to be removed

    ! .. logfile
    
    write(istr,'(I4)')rank
    if( size>9999) then 
        text="Error: size to large for status file number"
        ! call print_to_log(LogUnit,text)
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
    call error_handler(info,"read_inputfile")

    call init_constants()
    call make_geometry()                        ! generate volume elements lattice
    call allocate_chain_parameters()            !  
    call init_matrices()                        ! init matrices for MC chain generation
    call allocate_chains(cuantas,nnucl,nseg,nsegAA,nsegtypes,nsegtypesAA,maxnchainsrotations,maxnchainsrotationsxy)
    call init_chain_parameters()                                   ! reads pKa, pKaions zpol and vpol from file
    call make_sequence_chain(chainperiod,chaintype)                ! set isAmonomer, type_of_monomer and  type_of_monomer_char
    call make_charge_table(ismonomer_chargeable,zpol,nsegtypes)    ! set ismonomer_chargeable,
    call make_type_of_charge_table(type_of_charge,zpol,nsegtypes)  ! set type_of_charge
    call make_segcom(segcm,nnucl,segcmfname)
    call set_properties_chain(chainperiod,chaintype) 
    call set_mapping_num_to_char(mapping_num_to_char)
       
    ! init distributed volume  
    if(systype=="nucl_ionbin_sv".or.systype=="nucl_neutral_sv") call init_vnucl_type() ! ismonomer_chargable etc needs to be set 
        

    call make_VdWeps(info) 
    call error_handler(info,"make_VdWeps")
     
    call set_value_isVdW_on_values(nsegtypes, VdWeps, isVdW) 
    
    if(isVdW) then 
        call make_VdWcoeff(info)
        call error_handler(info,"make_VdWepscoeff")
    endif  

    call make_chains(chainmethod,systype)   
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
        

    ! temp override
    print*,"elemcharge"
    do i=1,nsegtypes
       print*,"elemcharge(",i,")=",elem_charge(i),"type=",mapping_num_to_char(i)
    enddo   


    ! .. loop over pH, or pKd etc  values

    if(runtype=="inputcspH".or.runtype=="inputMgpH".or.runtype=="inputcsKClpH") then 
        loop => pH
    else if (runtype=="rangepKd") then
        loop => pKd   
    else if (runtype=="rangedeltaGd") then
        loop => deltaGd  
    else if(runtype=="rangeVdWeps") then 
        loop => VdWscale    
    else if(runtype=="rangedielect") then 
        loop => dielectscale   
    else
        if(associated(loop)) nullify(loop) ! make explicit that no association is made
    endif  

    ! .. select variable with which list_array to associate

    if (runtype=="inputMgpH".or.runtype=="rangepKd".or.runtype=="rangeVdWeps".or.runtype=="rangedeltaGd") then
        call set_value_MgCl2(runtype,info)
        call error_handler(info,"set_value_MgCl2")

        num = num_cMgCl2
        list => cMgCl2_array
        list_val => cMgCl2

    else if(runtype=="inputcsKClpH") then
        call set_value_KCl(runtype,info)
        call error_handler(info,"set_value_KCl")
        
        num = num_cKCl
        list => cKCl_array
        list_val => cKCl

    else 

        call set_value_NaCl(runtype,info)
        call error_handler(info,"set_value_NaCl")
        num = num_cNaCl
        list => cNaCl_array
        list_val => cNaCl ! pointer points to target  cNaCl

    endif


    isfirstguess = .true.
    use_xstored = .false.       ! with both flags set false make_guess will set xguess equal to x
    iter = 0                    ! iteration counter

    if(loop%stepsize>0) then
        loop%val = loop%min
    else
        loop%val = loop%max
    endif

    call set_fcn()
     
    loopstepsizebegin = loop%stepsize
    list_val = list(1)                ! get value from array
    nlist_elem = 1
    maxlist_elem = num
    nlist_step = 0
    maxlist_step = 99

    do while (nlist_elem<=maxlist_elem .and. nlist_step<=maxlist_step )   ! loop over list items

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
            
            isfirstguess=(loop%val==loopbegin) 

            call init_vars_input()  ! sets up chem potentials
            
            flag_solver = 0

            if(rank==0) then     ! node rank=0
                call make_guess(x, xguess, isfirstguess,use_xstored,xstored)
                call solver(x, xguess, tol_conv, fnorm, isSolution)
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

                if(isSolution) then

                    call compute_vars_and_output()
                    if(isfirstguess) then
                       do i=1,neqint
                          xstored(i)=x(i)
                       enddo
                       use_xstored=.false.
                    endif
                
                    isfirstguess =.false.
                    loop%val=loop%val+loop%stepsize
                
                else if(abs(loop%val-loopbegin)>loopeps) then ! no solution and not first loop value

                    text="no solution: backstep" 
                    call print_to_log(LogUnit,text)
                    loop%stepsize=loop%stepsize/2.0d0   ! decrease increment
                    loop%val=loop%val-loop%stepsize     ! step back
                    
                    do i=1,neq
                        x(i)=xguess(i)
                    enddo
            
                else 
                    ! break while loop over loop%val by making loop%stepsize smaller loop%delta
                    loop%stepsize = loop%delta/2.0_dp 
                endif
                 
                
                ! communicate new values of loop from head to compute nodes to advance while loop on compute nodes
                do i = 1, size-1 !numproc-1
                    dest = i
                    call MPI_SEND(loop%val, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,ierr)
                    call MPI_SEND(loop%stepsize, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,ierr)
                enddo
            else ! receive values
                source = 0
                call MPI_RECV(loop%val, 1, MPI_DOUBLE_PRECISION, source, tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(loop%stepsize, 1, MPI_DOUBLE_PRECISION, source, tag,MPI_COMM_WORLD,stat, ierr)
            endif

            iter  = 0              ! reset of iteration counter

        enddo ! end while loop


        if(rank==0) then

            if(isSolution.or.(abs(loop%val-loopbegin)>loopeps) )then 
                if(abs(list_val-list(nlist_elem))<listeps) then 
                    if(nlist_elem<maxlist_elem) then ! prevents list exceed upper bond
                        list_step = list(nlist_elem+1) - list(nlist_elem)
                    endif    
                    nlist_step = nlist_step + 1
                    nlist_elem = nlist_elem + 1 ! advance element in input list values 
                else  
                    if(nlist_elem<maxlist_elem) then ! prevents list exceed upper bond
                        list_step = list(nlist_elem) - list_val
                    endif
                    nlist_step = nlist_step + 1
                endif        
                list_val = list_val+list_step                
            else 
                if(nlist_elem>1) then     ! not first element list
                    list_step = list_step/2.0_dp      
                    nlist_step = nlist_step +1
                    list_val = list_val-list_step  
                else if(nlist_elem==1) then 
                    nlist_step = maxlist_step+1       ! break out while over list
                    text="no solution first call of double while loop"
                    call print_to_log(LogUnit,text)
                    print*,text
                endif   
            endif       

        endif

        call MPI_Barrier(MPI_COMM_WORLD, ierr) ! synchronize 

        call MPI_Bcast(list_val,   1, MPI_DOUBLE_PRECISION, 0 ,MPI_COMM_WORLD, ierr)
        call MPI_Bcast(list_step,  1, MPI_DOUBLE_PRECISION, 0 ,MPI_COMM_WORLD, ierr)
        call MPI_Bcast(nlist_elem, 1, MPI_INTEGER, 0 ,MPI_COMM_WORLD, ierr)
        call MPI_Bcast(nlist_step, 1, MPI_INTEGER, 0 ,MPI_COMM_WORLD, ierr)

        use_xstored=.true.

    enddo ! end loop list item

    deallocate(x)
    deallocate(xguess)
    deallocate(fvec)

    call MPI_FINALIZE(ierr)

    deallocate(xstored)
      
    call deallocate_field()

    text="program end"

    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program main
