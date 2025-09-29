module initxvector
    
    use precision_definition, only : dp
    
    implicit none

    private 
    public  :: make_guess

contains


! Makes an inital guess vector xguess
! Order of assignment:
! if selectfirstguess true             : xguess => init_guess(x,xguess) 
! else if flagstored present and true  : xguess => xstored
! if flagstored false and isfirstguess true                 
!                                      : xguess => init_guess(x,xguess)  
! else                                 : xguess => x

subroutine make_guess(x, xguess, isfirstguess, selectfirstguess, flagstored, xstored)
  
    use myutils, only : error_handler, lenText

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential 
    logical, intent(in) :: isfirstguess        ! first guess  
    logical, intent(in) :: selectfirstguess    ! use first guess 
    logical, optional, intent(in) :: flagstored
    real(dp), optional, intent(in) :: xstored(:)
    
    !  ..local variables 
    character(len=lenText) :: text

    if(present(flagstored)) then
        if(present(xstored)) then
            if(selectfirstguess) then 
                call init_guess(x,xguess)
            else if(flagstored) then  
                call make_guess_from_xstored(xguess,xstored)
            else if(isfirstguess) then       ! first guess
                call init_guess(x,xguess)
            else
                xguess = x
            endif
        else
            text="Error: argument xstored not present, while flagstored present"
            call error_handler(1,text)
        endif 
    else if(isfirstguess) then       ! first guess
        call init_guess(x,xguess)
    else if(selectfirstguess) then 
        call init_guess(x,xguess)
    else     
        xguess = x     
    endif

end subroutine make_guess

! initial or first guess x and xguess 
! infile = 0  guess xsol = xbulk%sol and psi=0 for equilibrium and linear interpol for psi for steady state
! infile = 1  guess from input files for varialbe xsol.in psi.in etc
! infile = 2  guess from file x.out containing complete x vector
! infile = 3  guess is infile 0 option in use select first guess == true in make_guess

subroutine init_guess(x, xguess)
    
    use globals, only : systype

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential   

   ! print*," init guess systype=",systype," len=",len(systype)

    select case (systype)
        case ("elect")   
            call init_guess_elect(x,xguess)    
        case ("neutral")  
            call init_guess_neutral(x,xguess)
        case ("neutralnoVdW")  
            call init_guess_neutralnoVdW(x,xguess)
        case ("brush_mul")  
            call init_guess_multi(x,xguess)
        case ("brush_mulnoVdW")  
            call init_guess_multinoVdW(x,xguess)
        case ("brushdna","nucl_ionbin","nucl_ionbin_sv","nucl_ionbin_Mg","nucl_ionbin_MgA","nucl_ionbin_Fe") 
            call init_guess_multi(x,xguess)
        case ("nucl_neutral_sv")  
            call init_guess_nucl_neutral_sv(x,xguess) 
        case ("nonucl_ST","nucl_ionbin_Fe_ST","nucl_ionbin_MgA_ST") 
            call init_guess_nonucl_ST(x,xguess)
        case ("nonucl_ST_mu","nucl_ionbin_Fe_ST_mu","nucl_ionbin_MgA_ST_mu") 
            call init_guess_nonucl_ST_mu(x,xguess)
        case ("brushborn") 
            call init_guess_multi_born(x,xguess)
        case default   
            print*,"Init_guess: Wrong value systype : ", systype
    end select 

end subroutine init_guess


!     purpose: initalize x and xguess

subroutine init_guess_elect(x, xguess)

    use globals, only : bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
    use field, only : xsol,psi,rhopol
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname(4)
    integer :: ios,un_file(4)
    integer, parameter :: A=1, B=2   
    character(len=lenText) :: text, istr
  
   
    if(infile==0.or.infile==3) then  
        ! .. init guess all xbulk  
        x = 0.0_dp
        x(1:nsize) = xbulk%sol

    else if(infile==1) then   
        ! .. infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A7)')'rhoA.in'
        write(fname(4),'(A7)')'rhoB.in'
     
        do i=1,4 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                write(istr,'(I5)')un_file(i)
                text='init_guess_elect: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(i)))
                call error_handler(ios,text)
            endif
        enddo

        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)
            enddo
        endif            
        do i=1,nsize
            read(un_file(1),*)xsol(i)     ! solvent
            read(un_file(2),*)psi(i)      ! degree of complexation A
            read(un_file(3),*)rhopol(i,A) ! degree of complexation A
            read(un_file(4),*)rhopol(i,B) ! degree of complexation A
            x(i)         = xsol(i)        ! placing xsol  in vector x
            x(i+nsize)   = psi(i)         ! placing xsol  in vector x
            x(i+2*nsize) = rhopol(i,A)    ! placing xsol  in vector x
            x(i+3*nsize) = rhopol(i,B)    ! placing xsol  in vector x
        enddo
    
        if(bcflag(RIGHT)/="cc") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif            
       
        do i=1,4
            close(un_file(i))
        enddo

    !     .. end init from file 
   
    else if(infile==2) then 

        call read_xout(x)
    
    endif

    xguess = x
    
end subroutine init_guess_elect


subroutine init_guess_nucl_neutral_sv(x, xguess)
    
    use globals, only : nsize
    use field, only : xsol
    use parameters, only : xbulk, infile
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname
    character(len=lenText) :: text, istr
    integer :: ios,un_file
  
     if (infile==0.or.infile==3) then
        ! .. init guess all xbulk   
        x = 0.0_dp
        x(1:nsize) = xbulk%sol
    
    else if (infile==1) then   
        ! .. infile is read in from file/stdio  
        write(fname,'(A7)')'xsol.in'   
        open(unit=newunit(un_file),file=fname,iostat=ios,status='old')
        if(ios >0 ) then
            write(istr,'(I5)')un_file
            text='init_guess_neutral_sv: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname))
            call error_handler(ios,text)
        endif

        do i=1,nsize
            read(un_file,*)xsol(i) ! solvent
            x(i) = xsol(i)            ! placing xsol  in vector x
        enddo
        
        close(un_file)
   
        ! .. end init from file 
  
    else if(infile==2) then 

        call read_xout(x)
    
    endif

    xguess = x 
    
end subroutine init_guess_nucl_neutral_sv

subroutine init_guess_neutral(x, xguess)
    
    use globals, only : nsize,nsegtypes
    use field, only : xsol,rhopol,xpol
    use parameters, only : xbulk, infile, isrhoselfconsistent
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i, k, t
    character(len=8) :: fname(2)
    character(len=lenText) :: text, istr
    integer :: ios,un_file(2),count_scf
  

    if (infile==0.or.infile==3) then
        ! .. init guess all xbulk      
        x =0.0_dp    
        x(1:nsize)=xbulk%sol    
    else if(infile==1) then   
        ! ..infile is read in from file/stdio  
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A7)')'xpol.in'
        do i=1,2 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                write(istr,'(I5)')un_file(i)
                text='init_guess_neutral: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(i)))
                call error_handler(ios,text)
            endif
        enddo
        do i=1,nsize
            read(un_file(1),*)xsol(i) ! solvent
            read(un_file(2),*)xpol(i),(rhopol(i,t),t=1,nsegtypes)
            x(i) = xsol(i)            ! placing xsol  in vector x
        enddo
        
        count_scf=0                     ! placing density in vector x
        do t=1,nsegtypes
            if(isrhoselfconsistent(t)) then
                count_scf=count_scf+1 
                k=count_scf*nsize
                do i=1,nsize
                    x(i+k) = rhopol(i,t)                            
                enddo
            endif        
        enddo    

        close(un_file(1))
        close(un_file(2))

        ! .. end init from file 
  
    else if(infile==2) then 

        call read_xout(x)
    
    endif

    xguess = x
        
end subroutine init_guess_neutral

subroutine init_guess_neutralnoVdW(x, xguess)
    
    use globals, only : nsize
    use field, only : xsol
    use parameters, only : xbulk, infile
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname
    character(len=lenText) :: text, istr
    integer :: ios,un_file
  
    if (infile==0.or.infile==3) then
        !  .. init guess all xbulk      
        x =0.0_dp    
        x(1:nsize)=xbulk%sol
    
    else if(infile==1) then   
        ! .. infile is read in from file/stdio  
   
        write(fname,'(A7)')'xsol.in'
        open(unit=newunit(un_file),file=fname,iostat=ios,status='old')
        if(ios >0 ) then
            write(istr,'(I5)')un_file
            text='init_guess_neutralnoVdW: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname))
            call error_handler(ios,text)
        endif
        do i=1,nsize
            read(un_file,*)xsol(i) ! solvent
            x(i) = xsol(i)            ! placing xsol  in vector x   
        enddo     
        close(un_file)

        !  .. end init from file 

    else if(infile==2) then 

        call read_xout(x)
    
    endif

    xguess = x
    
end subroutine init_guess_neutralnoVdW


subroutine init_guess_multi(x, xguess)

    use globals, only : bcflag,LEFT,RIGHT,nsize,nsegtypes,systype
    use volume, only : nsurf
    use field, only : xsol,psi,rhopol,xpol,xpol_t
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile, isrhoselfconsistent
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i, k, t
    character(len=8) :: fname(4)
    integer :: ios,un_file(4),count_scf
    character(len=lenText) :: text, istr

    if (infile==0.or.infile==3) then
        ! .. init guess all xbulk     
        x=0.0_dp    
        x(1:nsize)=xbulk%sol

    else if (infile==1) then   
        ! i.. nfile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A7)')'xpol.in'
     
        do i=1,3 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                write(istr,'(I5)')un_file(i)
                text='init_guess_multi: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(i)))
                call error_handler(ios,text)
            endif    
        enddo

        if(bcflag(LEFT)/="cc" .and. bcflag(LEFT)/="cp") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)
            enddo
        endif    

        if(systype/="nucl_ionbin_sv".and.systype/="nucl_ionbin_Mg" .and. systype/="nucl_ionbin_MgA" &
            .and. systype/="nucl_ionbin_Fe") then
            do i=1,nsize
                read(un_file(1),*)xsol(i)    ! solvent
                read(un_file(2),*)psi(i)     ! potential
                read(un_file(3),*)xpol(i),(rhopol(i,t),t=1,nsegtypes)

                x(i)         = xsol(i)    ! placing xsol in vector x
                x(i+nsize)   = psi(i)     ! placing psi in vector x
            enddo 
        else
        
            do i=1,nsize
                read(un_file(1),*)xsol(i)    ! solvent
                read(un_file(2),*)psi(i)     ! potential
                read(un_file(3),*)xpol(i),(xpol_t(i,t),t=1,nsegtypes)  

                x(i)         = xsol(i)    ! placing xsol in vector x
                x(i+nsize)   = psi(i)     ! placing psi in vector x
            enddo 

        endif  

    
        count_scf=0                     ! placing density in vector x
        do t=1,nsegtypes
            if(isrhoselfconsistent(t)) then
                if(systype=="nucl_ionbin_sv".or.systype=="nucl_ionbin_Mg".or.systype=="nucl_ionbin_MgA") then
                    text='init_guess_multi: combination '//systype//' with VdW not working !'
                    print*,text
                    !call error_handler(-1,text)
                endif

                count_scf=count_scf+1 
                k=(count_scf+1)*nsize
                do i=1,nsize
                    x(i+k) = rhopol(i,t)                            
                enddo
            endif        
        enddo

        if(bcflag(RIGHT)/="cc" .and. bcflag(LEFT)/="cp") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif            
       
         do i=1,3
            close(un_file(i))
        enddo

    else if(infile==2) then 

        call read_xout(x)
    
    endif

    xguess = x
    
end subroutine init_guess_multi


subroutine init_guess_multinoVdW(x, xguess)

    use globals, only : bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
    use field, only : xsol,psi
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: i
    character(len=8) :: fname(2)
    character(len=lenText) :: text, istr
    integer :: ios,un_file(2)
  
        
    if(infile==0.or.infile==3) then ! .. init guess all xbulk 

        x =0.0_dp    
        x(1:nsize)=xbulk%sol
    
    else if(infile==1) then    ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
     
        do i=1,2 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                write(istr,'(I5)')un_file(i)
                text='init_guess_multinoVdW: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(i)))
                call error_handler(ios,text)
            endif
        enddo
        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)
            enddo
        endif            
        do i=1,nsize
            read(un_file(1),*)xsol(i)    ! solvent
            read(un_file(2),*)psi(i)     ! potential
        
            x(i)         = xsol(i)    ! placing xsol in vector x
            x(i+nsize)   = psi(i)     ! placing psi in vector x
                  
        enddo
    
        if(bcflag(RIGHT)/="cc") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif            
       
         do i=1,2
            close(un_file(i))
        enddo
        !     .. end init from file 

    else if(infile==2) then 

        call read_xout(x)
    
    endif

    xguess = x
  
end subroutine init_guess_multinoVdW


! Makes an inital guess vector xguess for systype nonucl_ST

subroutine init_guess_nonucl_ST(x, xguess)

    use globals, only : LEFT, RIGHT, nsize, DEBUG_ST
    use field, only : xsol, psi, xNa, xCl, xK, xHplus, xOHmin, xMg, xFe2, xFe3, xCa 
    use parameters, only : xbulk, infile, xvolmin, psiSL, psiSR
    use parameters, only : niontypes, iontype, isionselfconsistent
    use myutils, only : newunit, lenText, error_handler
    use molecules, only : get_value_moleclist, sum_value_moleclist

  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !  ..local variables 
    integer :: i, t, k
    character(len=9) :: fname(11)
    character(len=lenText) :: text, istr
    integer :: ios,un_file(11)
    character(len=5) :: key
    real(dp) :: xvol, xtest

    if(infile==0.or.infile==3) then 

        ! .. init guess all xbulk     

        call linear_interpolation(psi,psiSL,psiSR)
    
        do i=1,nsize
            x(i)       = xbulk%sol
            x(i+nsize) = psi(i)
        enddo 
        
        xtest=0.0_dp
        k = 2*nsize
        do t=1,size(iontype)
            if(isionselfconsistent(t)) then 
                key = trim(iontype(t))
                xvol = get_value_moleclist(xvolmin,key)
                xtest = xtest +xvol
                if(DEBUG_ST) print*,"t= ",t," key= ",key," xvol=",xvol
        
                x(k+1:k+nsize)=xvol
                k = k + nsize
            endif           
        enddo

        xtest=xtest+xbulk%sol
        if(DEBUG_ST)  print*,"xtest=",xtest

        xtest = sum_value_moleclist(xbulk)
        if(DEBUG_ST) print*,"xtest=",xtest


    else if(infile==1) then   ! infile is read in from file/stdio  
    
        write(fname(10),'(A7)')'xsol.in'
        write(fname(11),'(A6)')'psi.in'
   
        write(fname(1),'(A6)')'xNa.in'
        write(fname(2),'(A6)')'xCl.in'
        write(fname(3),'(A5)')'xK.in'
        write(fname(4),'(A9)')'xHplus.in'
        write(fname(5),'(A9)')'xOHmin.in'
        write(fname(6),'(A6)')'xMg.in'
        write(fname(7),'(A9)')'xFe2.in'
        write(fname(8),'(A9)')'xFe3.in'
        write(fname(9),'(A6)')'xCa.in'

        do t=1,size(iontype) !  open ion files 
            if(isionselfconsistent(t)) then 
                open(unit=newunit(un_file(t)),file=fname(t),iostat=ios,status='old')
                if(ios >0 ) then
                    write(istr,'(I5)')un_file(t)
                    text='init_guess_nonucl_ST: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(t)))
                    call error_handler(ios,text)
                endif
            endif    
        enddo
        do t=10,11 ! open xsol and psi files 
            open(unit=newunit(un_file(t)),file=fname(t),iostat=ios,status='old')
            if(ios >0 ) then
                write(istr,'(I5)')un_file(t)
                text='init_guess_nonucl_ST: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(t)))
                call error_handler(ios,text)
            endif
        enddo
           
        do i=1,nsize
            read(un_file(10),*)xsol(i)    ! solvent
            read(un_file(11),*)psi(i)     ! potential
        
            x(i)         = xsol(i)    ! placing xsol in vector x
            x(i+nsize)   = psi(i)     ! placing psi in vector x      
        enddo
                
        k=2*nsize
        do t=1,niontypes 
            if(isionselfconsistent(t)) then 
                select case (iontype(t))
                case ("Hplus")
                    do i=1,nsize     
                        read(un_file(t),*)xHplus(i)  
                    enddo
                    x(k+1:k+nsize)=xHplus
                case( "OHmin") 
                    do i=1,nsize     
                        read(un_file(t),*)xOHmin(i)  
                    enddo  
                     x(k+1:k+nsize)=xOHmin
                case("Na")
                    do i=1,nsize     
                        read(un_file(t),*)xNa(i)
                    enddo  
                     x(k+1:k+nsize)=xNa
                case("K")
                    do i=1,nsize     
                        read(un_file(t),*)xK(i)
                    enddo  
                     x(k+1:k+nsize)=xK
                case("Cl")
                    do i=1,nsize   
                        read(un_file(t),*)xCl(i) 
                    enddo  
                     x(k+1:k+nsize)=xCl
                case("Mg")
                    do i=1,nsize     
                        read(un_file(t),*)xMg(i)
                    enddo   
                     x(k+1:k+nsize)=xMg
                case("Fe2")
                    do i=1,nsize     
                        read(un_file(t),*)xFe2(i) 
                    enddo 
                     x(k+1:k+nsize)=xFe2
                case("Fe3")
                    do i=1,nsize     
                        read(un_file(t),*)xFe3(i)
                    enddo  
                    x(k+1:k+nsize)=xFe3  
                case("Ca")
                    do i=1,nsize     
                        read(un_file(t),*)xCa(i)
                    enddo  
                    x(k+1:k+nsize)=xCa     
                case default
                end select
                k = k + nsize
            endif
        enddo


        do i=10,11
            close(un_file(i))
        enddo

        do t=1,size(iontype) !  open ion files 
            if(isionselfconsistent(t)) close(un_file(t))   
        enddo

    else if(infile==2) then   ! x  read  directly from x.out
       
        call read_xout(x)

    endif

    ! assign xguess to x 
  
    xguess =x 
  
end subroutine init_guess_nonucl_ST


! Makes an inital guess vector xguess for systype nonucl_ST_mu

subroutine init_guess_nonucl_ST_mu(x, xguess)

    use globals, only : LEFT, RIGHT, nsize, DEBUG_ST
    use field, only : xsol, psi, xNa, xCl, xK, xHplus, xOHmin, xMg, xFe2, xFe3, xCa 
    use parameters, only : xbulk, infile, mumin,mumax, psiSL, psiSR
    use parameters, only : niontypes, iontype, isionselfconsistent
    use myutils, only : newunit, lenText, error_handler
    use molecules, only : get_value_moleclist, sum_value_moleclist
    use flux, only : chem_potential
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !  ..local variables 
    integer :: i, t, k
    character(len=9) :: fname(11)
    character(len=lenText) :: text, istr
    integer :: ios,un_file(11)
    character(len=5) :: key
    real(dp) :: mu(nsize), muL , muR


    if(infile==0.or.infile==3) then 
        ! .. init guess xsol and psi 

        call linear_interpolation(psi,psiSL,psiSR)
        
        do i=1,nsize
            x(i)       = xbulk%sol
            x(i+nsize) = psi(i)
        enddo 
        
        ! iontype(9)=(/"Na   ","Cl   ","K    ","Hplus","OHmin","Mg   ","Fe2  ","Fe3  ","Ca   "/) 

        k = 2*nsize
        do t=1,size(iontype)
            if(isionselfconsistent(t)) then 
                key = trim(iontype(t))
                muL = get_value_moleclist(mumin,key)
                muR = get_value_moleclist(mumax,key)
                call linear_interpolation(mu,muL,muR)
        
                x(k+1:k+nsize) = mu
                k = k + nsize
            endif           
        enddo

    else if(infile==1) then   
        ! .. infile is read in from file/stdio  
        
        write(fname(10),'(A7)')'xsol.in'
        write(fname(11),'(A6)')'psi.in'
   
        write(fname(1),'(A6)')'xNa.in'
        write(fname(2),'(A6)')'xCl.in'
        write(fname(3),'(A5)')'xK.in'
        write(fname(4),'(A9)')'xHplus.in'
        write(fname(5),'(A9)')'xOHmin.in'
        write(fname(6),'(A6)')'xMg.in'
        write(fname(7),'(A9)')'xFe2.in'
        write(fname(8),'(A9)')'xFe3.in'
        write(fname(9),'(A6)')'xCa.in'

        do t=1,size(iontype) !  open ion files 
            if(isionselfconsistent(t)) then 
                open(unit=newunit(un_file(t)),file=fname(t),iostat=ios,status='old')
                if(ios >0 ) then
                    write(istr,'(I5)')un_file(t)
                    text='init_guess_nonucl_ST_mu: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(t)))
                    call error_handler(ios,text)
                endif
            endif    
        enddo
        do t=10,11 ! open xsol and psi files 
            open(unit=newunit(un_file(t)),file=fname(t),iostat=ios,status='old')
            if(ios >0 ) then
                write(istr,'(I5)')un_file(t)
                text='init_guess_nonucl_ST_mu: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(t)))
                call error_handler(ios,text)
            endif
        enddo
           
        do i=1,nsize
            read(un_file(10),*)xsol(i)    ! solvent
            read(un_file(11),*)psi(i)     ! potential
        
            x(i)         = xsol(i)    ! placing xsol in vector x
            x(i+nsize)   = psi(i)     ! placing psi in vector x      
        enddo
                
        k=2*nsize
        do t=1,niontypes 
            if(isionselfconsistent(t)) then 
                select case (iontype(t))
                case ("Hplus")
                    do i=1,nsize     
                        read(un_file(t),*)xHplus(i)  
                    enddo
                    call chem_potential(mu,xsol,xHplus,psi,"Hplus")
                    x(k+1:k+nsize) = mu
                case( "OHmin") 
                    do i=1,nsize     
                        read(un_file(t),*)xOHmin(i)  
                    enddo  
                    call chem_potential(mu,xsol,xOHmin,psi,"OHmin")
                    x(k+1:k+nsize) = mu
                case("Na")
                    do i=1,nsize     
                        read(un_file(t),*)xNa(i)
                    enddo  
                    call chem_potential(mu,xsol,xNa,psi,"Na")
                    x(k+1:k+nsize) = mu
                case("K")
                    do i=1,nsize     
                        read(un_file(t),*)xK(i)
                    enddo  
                    call chem_potential(mu,xsol,xK,psi,"K") 
                    x(k+1:k+nsize) = mu
                case("Cl")
                    do i=1,nsize   
                        read(un_file(t),*)xCl(i) 
                    enddo  
                    call chem_potential(mu,xsol,xCl,psi,"Cl")  
                    x(k+1:k+nsize) = mu
                case("Mg")
                    do i=1,nsize     
                        read(un_file(t),*)xMg(i)
                    enddo   
                    call chem_potential(mu,xsol,xMg,psi,"Mg") 
                    x(k+1:k+nsize) = mu
                case("Fe2")
                    do i=1,nsize     
                        read(un_file(t),*)xFe2(i) 
                    enddo 
                    call chem_potential(mu,xsol,xFe2,psi,"Fe2") 
                    x(k+1:k+nsize) = mu
                case("Fe3")
                    do i=1,nsize     
                        read(un_file(t),*)xFe3(i)
                    enddo  
                    call chem_potential(mu,xsol,xFe3,psi,"Fe3")
                    x(k+1:k+nsize) = mu
                case("Ca")
                    do i=1,nsize     
                        read(un_file(t),*)xCa(i)
                    enddo 
                    call chem_potential(mu,xsol,xCa,psi,"Ca") 
                    x(k+1:k+nsize) = mu     
                case default
                end select
                k = k + nsize
            endif
        enddo


        do i=10,11
            close(un_file(i))
        enddo

        do t=1,size(iontype) !  open ion files 
            if(isionselfconsistent(t)) close(un_file(t))   
        enddo

    else if(infile==2) then   ! x  read  directly from x.out

        call read_xout(x)

    endif
    
    ! assign xguess to x 
  
    xguess =x 
  
end subroutine init_guess_nonucl_ST_mu

subroutine init_guess_multi_born(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize,nsegtypes
    use volume, only : nsurf
    use field, only : xsol,psi,rhopol,xpol,rhopol,fdisA
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile, isrhoselfconsistent, tA
    use myutils, only : newunit, lenText, error_handler
    real(dp) ::  x(neq)       ! volume fraction solvent iteration vector 
    real(dp) ::  xguess(neq)  ! guess fraction  solvent 

    !     ..local variables 
    integer :: i, t, k, k1, k2, k3, k4, k5,neq_bc
    character(len=8) :: fname(4)
    character(len=lenText) :: text, istr
    integer :: ios, un_file(4), count_sc

    if(infile ==0.or.infile == 3) then

        x=0.0_dp    
        x(1:nsize)=xbulk%sol
   
    else if(infile==1) then   ! infile is read in from file/stdio  
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A7)')'xpol.in'
        write(fname(4),'(A8)')'fdisA.in'   

        do i=1,4
            un_file(i)=newunit()
            open(unit=un_file(i),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                write(istr,'(I5)')un_file(i)
                text='init_guess_multi_born: file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname(i)))
                call error_handler(ios,text)
            endif
        enddo

        k1=nsize
        k2=2*nsize
        k3=3*nsize
        k4=4*nsize
        k5=5*nsize

        if(bcflag(RIGHT)/="cc") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif        

        do i=1,nsize
            read(un_file(1),*)xsol(i) ! solvent
            read(un_file(2),*)psi(i)
            read(un_file(3),*)xpol(i),(rhopol(i,t),t=1,nsegtypes)
            read(un_file(4),*)(fdisA(i,k),k=1,8)
            
            x(i)    = xsol(i)   ! placing xsol in vector x
            x(i+k1) = psi(i)
            x(i+k2) = xpol(i)
            x(i+k3) = fdisA(i,4)*rhopol(i,tA)
            x(i+k4) = fdisA(i,1)*rhopol(i,tA)
            x(i+k5) = fdisA(i,6)*rhopol(i,tA)
        enddo         

        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)
            enddo
        endif  

        count_sc=0    
        do t=1,nsegtypes
            if(isrhoselfconsistent(t)) then
                count_sc=count_sc+1 
                k=count_sc*nsize+k5
                do i=1,nsize  
                    x(i+k) = rhopol(i,t) 
                enddo   
            endif        
        enddo

        
        neq_bc=0
        k=count_sc+1
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nsurf
            do i=1,neq_bc
                x(k+i)   =psiSurfR(i)                  ! surface potentail
            enddo
        endif   
        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf
                x(k+neq_bc+i) = psiSurfL(i)           ! surface potentail
            enddo
        endif

        do i=1,4
            close(un_file(i))
        enddo

    else if (infile==2) then   ! x  read  directly from x.out

        call read_xout(x)

    endif
    
    ! assign xguess to x 
  
    xguess = x 

end subroutine init_guess_multi_born


! .. copy solution of previous solution to create new guess


subroutine make_guess_from_xstored(xguess,xstored)

    use globals, only : neqint

    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potentia
    real(dp), intent(in) :: xstored(:)

    !   .. local variables
    integer :: i
   
    do i=1,neqint
        xguess(i)=xstored(i)     
    enddo 

end subroutine make_guess_from_xstored


! linear function in z-direction of fcn bewteen 
! fcn_begin at z= -delta /2 and fc 
! fcn_end.  at z= nz delta -delta /2 
! value in x and y direction same  

subroutine linear_interpolation(fcn_interp,fcn_begin,fcn_end)

    use globals, only : DEBUG_ST
    use volume, only : delta, nx, ny, nz, coordtoindex


    real(dp), intent(inout) :: fcn_interp(:)
    real(dp), intent(in) :: fcn_begin, fcn_end

    ! local variable 
    real(dp) :: slope, intercept,  fcn_val 
    integer :: i, j, k, indx
     
    slope = (fcn_end-fcn_begin)/((nz+1)*delta)
    intercept = fcn_begin+slope * delta/2.0_dp
    
    do k=1,nz 
        fcn_val = slope * (k - 0.5_dp) * delta  +  intercept  ! middle of lattice cell in z-direction
        do j=1,ny
            do i=1,nx 
                indx=coordtoindex(i,j,k)
                fcn_interp(indx)=fcn_val
            enddo
        enddo
    enddo         
    
end subroutine linear_interpolation

! Open and read file x.out ( dump file from solver) 
! assign vector to input vector x 

subroutine read_xout(x)

    use myutils, only : newunit, lenText, error_handler
    use globals, only : neqint

    real(dp) :: x(:)       
    
    !  ..local variables 
    integer :: i
    character(len=5) :: fname
    character(len=lenText) :: text, istr
    integer :: ios, un_file

    write(fname,'(A5)')'x.out'
    open(unit=newunit(un_file),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        write(istr,'(I5)')un_file
        text='Open file failed : file number = '//trim(adjustl(istr))//' file name = '//trim(adjustl(fname))                
        call error_handler(ios,text)
    endif
    
    do i=1,neqint     
        read(un_file,*)x(i)
    enddo 

    close(un_file)

end subroutine read_xout       



end module initxvector
