module myutils

    use precision_definition
    implicit none

    character(len=10) :: LogName        ! filename of status log file
    integer, parameter :: LogUnit=321   ! associated unit number with
    integer, parameter :: lenText=80    ! number of characters of one line of status log file
    
contains    
    
    ! printing to a log file 

    subroutine print_to_log(UnitNum,text)

        implicit none

        !Formal argument
        character (len=lenText), intent(in) :: text
        integer, intent (in) :: UnitNum
        !Local variables
        character (len=8)  :: date
        character (len=10) :: time
        character (len=26) :: date_time

        !Executable part
        call date_and_time(date,time)
        date_time='['//date(7:8)//'-'//date(5:6)//'-'//date(1:4)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)//'] '
        write(UnitNum,*) date_time, text

    end subroutine print_to_log


     ! printing to a log file 

    subroutine print_to_log_new(UnitNum,text)

        implicit none

        !Formal argument
        character (len=*), intent(in) :: text
        integer, intent (in) :: UnitNum
        !Local variables
        character (len=8)  :: date
        character (len=10) :: time
        character (len=26) :: date_time

        !Executable part
        call date_and_time(date,time)
        date_time='['//date(7:8)//'-'//date(5:6)//'-'//date(1:4)//' '//time(1:2)//':'//time(3:4)//':'//time(5:10)//'] '
        write(UnitNum,*) date_time//trim(text)

    end subroutine print_to_log_new


    subroutine open_logfile(UnitNum,FileName)  

        implicit none
        integer, intent (in) :: UnitNum
        character (len=*), intent (in) :: FileName
        logical :: exist
        integer :: ios

        inquire(file=FileName, exist=exist)
        if (exist) then
            open(UnitNum,file=FileName, iostat=ios,status="old", position="append", action="write")
        else
            open(UnitNum,file=FileName, iostat=ios,status="new", action="write")
        endif
      
        if(ios > 0 ) then
            print*, 'Error opening file : iostat =', ios
            stop
        endif

    end subroutine open_logfile



    subroutine close_logfile(UnitNum)

        implicit none

        integer, intent (in) :: UnitNum

        close(UnitNum)

    end subroutine close_logfile
        

    ! in fortran 2008 newunit is provided 
 
    integer function newunit(unit)
        implicit none

        integer, intent(out), optional :: unit
        ! local
        integer, parameter :: LUN_MIN=10, LUN_MAX=1000
        logical :: opened
        integer :: lun
        ! begin
        newunit=-1
        do lun=LUN_MIN,LUN_MAX
            inquire(unit=lun,opened=opened)
            if (.not. opened) then
                newunit=lun
                exit
            endif
        enddo
        if (present(unit)) unit=newunit
    end function newunit


    logical function isNaN(x)
        implicit none
        real(dp) :: x
         
        if (x /= x) then
            isNaN=.true.
        else
            isNaN=.false.
        endif 

    end function isNaN


end module myutils
