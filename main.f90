module config
    implicit none
    public

    !Segment boundaries
    real(8), parameter :: fl = 0.d0, cl = 1.d0
    !Integration precision
    integer, parameter :: dots = 100
end module config

program output
    use config
    implicit none

    procedure(func), pointer :: fptr, gptr
    procedure(rkfunc), pointer :: rptr, kptr

    fptr => func
    gptr => gunc

    call cowell(fptr, gptr, fl, cl)

    rptr => rkfunc
    kptr => rkgunc

    call runge(rptr, kptr, fl, cl)


contains

    ! y'' = func(x)*y + gunc(x)

    pure real(8) function func(x)
        implicit none
        real(8), intent(in) :: x
        func = 1
    end function func

    pure real(8) function gunc(x)
        implicit none
        real(8), intent(in) :: x
        gunc = 2*(x**2-1)*cos(x)
    end function gunc

    subroutine cowell(fptr, gptr, st, ed)
        use config
        implicit none

        procedure(func), pointer, intent(in) :: fptr, gptr
        real(8), intent(in) :: st, ed
        real(8) :: x0, x1, x2, y0, y1, y2, step
        integer :: i

        step = (ed - st) / dots

        !Initial conditions
        y0 = 2.d0
        y1 = (2-step**2)*cos(step)+2*step*sin(step) !Main source of error

        open(unit = 1, file = "cowell.dat")
        do i = 1, dots+1
            x0 = st + (i-1)*step
            x1 = x0 + step
            x2 = x1 + step
            y2 = 2*y1*(1+5*(step**2)*fptr(x1)/12) - y0*(1-(step**2)*fptr(x0)/12) + (step**2)*(gptr(x2)+10*gptr(x1)+gptr(x0))/12
            y2 = y2/(1-(step**2)*fptr(x2)/12)
            write(1,*) x2, y2
            y0 = y1
            y1 = y2
        end do
        close(1)

    end subroutine cowell

    pure real(8) function sltn(x)
        !Red dashed line
        implicit none
        real(8), intent(in) :: x
        sltn = (2 - x**2) * cos(x) + 2*x*sin(x)
    end function sltn

    !y'=rkfunc(x, y, z)
    !z'=rkgunc(x, y, z)

    pure real(8) function rkfunc(x, y, z)
        implicit none
        real(8), intent(in) :: x, y, z
        rkfunc = x**2*z
    end function rkfunc

    pure real(8) function rkgunc(x, y, z)
        implicit none
        real(8), intent(in) :: x, y, z
        rkgunc = cos(x)
    end function rkgunc

    subroutine runge(fptr, gptr, st, ed)
        use config
        implicit none

        procedure(rkfunc), pointer, intent(in) :: fptr, gptr
        real(8), intent(in) :: st, ed
        real(8), dimension(4) :: k, l
        real(8) :: x0, x, y, z, step
        integer :: i, j

        step = (ed - st) / dots

        !Initial conditions
        x0 = 0.d0
        y = 2.d0
        z = 0.d0

        open(unit = 1, file = "runge.dat")
        do i = 1, dots+1
            x = x0 + (i-1)*step
            write(1,*) x, y, sltn(x)

            k(1) = fptr(x, y, z)
            l(1) = gptr(x, y, z)

            k(2) = fptr(x+step/2, y+k(1)*step/2, z+k(1)*step/2)
            l(2) = gptr(x+step/2, y+l(1)*step/2, z+l(1)*step/2)

            k(3) = fptr(x+step/2, y+k(2)*step/2, z+k(2)*step/2)
            l(3) = gptr(x+step/2, y+l(2)*step/2, z+l(2)*step/2)

            k(4) = fptr(x+step, y+k(3)*step, z+k(3)*step)
            l(4) = gptr(x+step, y+l(3)*step, z+l(3)*step)

            y = y + step*(k(1)+2*k(2)+2*k(3)+k(4))/6
            z = z + step*(l(1)+2*l(2)+2*l(3)+l(4))/6


        end do
        close(1)

    end subroutine runge

end program output
