module error_prt
! Error termination of program
  use globalmod_array, only: idwrite
  implicit none

  private
  public perr

contains

  subroutine perr (msg, code)
! Print error message and stop    
    character(LEN=*), intent(in)      :: msg
    integer, optional, intent(in)     :: code

    if (PRESENT(code)) then
       write (idwrite,'(a,i6)') TRIM(msg), code
    else 
       write (idwrite,'(a)') TRIM(msg)
    end if
    stop
  end subroutine perr
end module error_prt
