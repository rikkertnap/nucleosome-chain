module initxvector
    
    use precision_definition, only : dp
    
    implicit none

    private 
    public  :: make_guess

contains

subroutine make_guess(x, xguess, isfirstguess, flagstored, xstored)
  
    use globals, only : neq,neqmax,systype,bcflag,LEFT,RIGHT
    use volume, only : nsurf 

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential 
    logical, intent(in) :: isfirstguess     ! first guess   
    logical, optional, intent(in) :: flagstored
    real(dp), optional, intent(in) :: xstored(:)
    
    !  ..local variables 
    integer :: i ,neq_bc

    if(present(flagstored)) then
        if(present(xstored)) then
            if(flagstored) then  
                call make_guess_from_xstored(xguess,xstored)
            else if(isfirstguess) then       ! first guess
                call init_guess(x,xguess)
            else  
                do i=1,neq
                    xguess(i)=x(i)      
                enddo
            endif
        else
            print*,"Error: argument xstored not present, while flagstored present"
            stop 
        endif 
    else if(isfirstguess) then       ! first guess
        call init_guess(x,xguess)
    else     
        do i=1,neq
            xguess(i)=x(i)     
        enddo
    endif

end subroutine make_guess



subroutine init_guess(x, xguess)
    
    use globals, only : systype

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential   

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
        case ("brushdna","nucl_ionbin","nucl_ionbin_sv","nucl_ionbin_Mg")  
            call init_guess_multi(x,xguess)
        case ("nucl_neutral_sv")  
            call init_guess_nucl_neutral_sv(x,xguess)
        case ("brushborn") 
            call init_guess_multi_born(x,xguess)
        case default   
            print*,"Init_guess: Wrong value systype : ", systype
    end select 

end subroutine init_guess


!     purpose: initalize x and xguess

subroutine init_guess_elect(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
    use field, only : xsol,psi,rhopol
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(4)
    integer :: ios,un_file(4)
    integer, parameter :: A=1, B=2   
    character(len=lenText) :: text, istr, str
  
    ! .. init guess all xbulk     

    do i=1,nsize
        x(i)=xbulk%sol
        x(i+nsize)=0.0_dp
        x(i+2*nsize)=0.0_dp
        x(i+3*nsize)=0.0_dp
    enddo
   
    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
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
            read(un_file(1),*)xsol(i)    ! solvent
            read(un_file(2),*)psi(i)     ! degree of complexation A
            read(un_file(3),*)rhopol(i,A) ! degree of complexation A
            read(un_file(4),*)rhopol(i,B)   ! degree of complexation A
            x(i)         = xsol(i)    ! placing xsol  in vector x
            x(i+nsize)   = psi(i)     ! placing xsol  in vector x
            x(i+2*nsize) = rhopol(i,A) ! placing xsol  in vector x
            x(i+3*nsize) = rhopol(i,B)   ! placing xsol  in vector x
        enddo
    
        if(bcflag(RIGHT)/="cc") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif            
       
        do i=1,4
            close(un_file(i))
        enddo

    endif

    !     .. end init from file 
    do i=1,neq
        xguess(i)=x(i)
    enddo

    
end subroutine init_guess_elect


subroutine init_guess_nucl_neutral_sv(x, xguess)
    
    use globals, only : neqint,nsize,nsegtypes
    use volume, only : nz
    use field, only : xsol,rhopol,xpol
    use parameters, only : xbulk, infile, isrhoselfconsistent
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, k, t
    character(len=8) :: fname
    character(len=lenText) :: text, istr
    integer :: ios,un_file,count_scf
  
    !     .. init guess all xbulk      

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo

    if (infile.eq.1) then   ! infile is read in from file/stdio  
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
    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo
    
end subroutine init_guess_nucl_neutral_sv



subroutine init_guess_neutral(x, xguess)
    
    use globals, only : neqint,nsize,nsegtypes
    use volume, only : nz
    use field, only : xsol,rhopol,xpol
    use parameters, only : xbulk, infile, isrhoselfconsistent
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, k, t
    character(len=8) :: fname(2)
    character(len=lenText) :: text, istr, str
    integer :: ios,un_file(2),count_scf
  
    !     .. init guess all xbulk      

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo


    if (infile.eq.1) then   ! infile is read in from file/stdio  
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
    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo
    
end subroutine init_guess_neutral

subroutine init_guess_neutralnoVdW(x, xguess)
    
    use globals, only : neqint,nsize,nsegtypes
    use volume, only : nz
    use field, only : xsol,rhopol
    use parameters, only : xbulk, infile
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, t
    character(len=8) :: fname
    character(len=lenText) :: text, istr
    integer :: ios,un_file
  
    !     .. init guess all xbulk      

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo


    if (infile.eq.1) then   ! infile is read in from file/stdio  
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
    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo
    
end subroutine init_guess_neutralnoVdW


subroutine init_guess_multi(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize,neqint,nsegtypes,systype
    use volume, only : nsurf
    use field, only : xsol,psi,rhopol,xpol,xpol_t
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile, isrhoselfconsistent
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, k, t
    character(len=8) :: fname(4)
    integer :: ios,un_file(4),count_scf
    character(len=lenText) :: text, istr

    ! .. init guess all xbulk     

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo


    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
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

        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)
            enddo
        endif    

        if(systype/="nucl_ionbin_sv".and.systype/="nucl_ionbin_Mg") then
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
                if(systype=="nucl_ionbin_sv".or.systype=="nucl_ionbin_Mg") then
                    text='init_guess_multi: combination '//systype//' with VdW not working !'
                    !call error_handler(-1,text)
                endif

                count_scf=count_scf+1 
                k=(count_scf+1)*nsize
                do i=1,nsize
                    x(i+k) = rhopol(i,t)                            
                enddo
            endif        
        enddo

        if(bcflag(RIGHT)/="cc") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif            
       
         do i=1,3
            close(un_file(i))
        enddo

    endif
    !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo

end subroutine init_guess_multi


subroutine init_guess_multinoVdW(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize,neqint
    use volume, only : nsurf
    use field, only : xsol,psi
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit, lenText, error_handler
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, t
    character(len=8) :: fname(2)
    character(len=lenText) :: text, istr
    integer :: ios,un_file(2)
  
    ! .. init guess all xbulk     

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo

    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
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

    endif
    !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo

end subroutine init_guess_multinoVdW


subroutine init_guess_multi_born(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize,neqint,nsegtypes
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
    real(dp) :: val ! dummy variable for reading in files


    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo

    if (infile.eq.1) then   ! infile is read in from file/stdio  
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

    endif

    !  .. end init from file 

    do i=1,neqint
        xguess(i)=x(i)
    enddo

    

end subroutine init_guess_multi_born



! .. copy solution of previous solution to create new guess


subroutine make_guess_from_xstored(xguess,xstored)

    use globals, only : neq

    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potentia
    real(dp), intent(in) :: xstored(:)

    !   .. local variables
    integer :: i
   
    do i=1,neq
        xguess(i)=xstored(i)     
    enddo 

end subroutine make_guess_from_xstored



end module initxvector
