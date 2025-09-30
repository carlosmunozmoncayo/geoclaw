module debug_module
    contains
    subroutine debug_log(boundary_flag,t,x,y,q,aux)
    implicit none
    integer, intent(in) :: boundary_flag
    real(kind=8), intent(in) :: t, x,y
    real(kind=8), intent(in) :: q(:)
    real(kind=8), intent(in) :: aux(:)

    character(len=140) :: filepath
    integer :: unit

    character(len=256) :: claw_env
    integer :: length, status

    

    ! Get CLAW environment variable
    call get_environment_variable("CLAW", claw_env, length=length, status=status)

    if (status /= 0) then
        print *, "Error: environment variable CLAW not set!"
        stop
    endif


    filepath = trim(claw_env(1:length)) // "/geoclaw/examples/tsunami/simple_setup_ghost_cells/debug_log.txt"

    ! Open file in append mode
    open(newunit=unit, file=trim(filepath), status='unknown', position='append', action='write')

    ! Write x, y, q, aux to file in a single line
    write(unit,*) boundary_flag, t, x, y, q, aux

    close(unit)
end subroutine debug_log
end module