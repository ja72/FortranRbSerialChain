    module mod_ode
    use mod_lu
    implicit none

    interface
    function f_rate(x,y) result(yp)
    import
    real(wp), intent(in) :: x, y(:)
    real(wp) :: yp(size(y))
    end function
    end interface
    
    type ode_sys
        procedure(f_rate), pointer, nopass :: f    
        real(wp) :: x, x_end, e, h
        real(wp), allocatable :: y(:), yp(:)
    contains
        procedure :: slope => ode_slope_ext
        procedure :: rk4 => ode_rk4
        procedure :: rk45 => ode_rk45
        procedure :: done => ode_done
    end type
    
    interface ode_sys
        procedure :: ode_sys_new
    end interface

    contains

    function ode_done(u) result(ok)
    class(ode_sys), intent(in) :: u
    logical :: ok
        ok = u%x >= u%x_end
    end function
    
    function ode_sys_new(f,y,x_end,n) result(u)
        procedure(f_rate), intent(in), pointer :: f
        real(wp),intent(in) :: y(:), x_end
        type(ode_sys) :: u
        integer, optional :: n
        real(wp), allocatable :: yt(:), yp(:)
            u%x = 0d0
            u%x_end = x_end
            if( present(n) )then
                u%h = x_end/n
            else
                u%h = x_end/360
            end if
            u%f => f
            u%y = y        
            u%yp = f(u%x, u%y)
            yp = u%yp
            yt = rk45(u%f, u%x, u%y, u%h, u%e, yp)
    end function
    
    function ode_slope_ext(u,x,y) result(yp)
    class(ode_sys), intent(in) :: u
    real(wp), intent(in) :: x, y(:)    
    real(wp), allocatable :: yp(:)
        yp = u%f(x, y)    
    end function
    
    subroutine ode_rk4(u)
    class(ode_sys), intent(inout) :: u
    real(wp) :: y(size(u%y))
    real(wp) :: h, k0(size(u%y)), k1(size(u%y)), k2(size(u%y)), k3(size(u%y))
        h = min(u%h, u%x_end - u%x)        
        k0 = u%yp
        k1 = u%f(u%x+h/2, u%y+h/2*k0)
        k2 = u%f(u%x+h/2, u%y+h/2*k1)
        k3 = u%f(u%x+h, u%y+h*k2)    
        y = u%y + h/6 * (k0 + 2*k1 + 2*k2 + k3)
        u%h = h
        u%x = u%x + h
        u%y = y        
        u%yp = u%f(u%x, y)
    end subroutine
    
    subroutine ode_rk45(u)
    class(ode_sys), intent(inout) :: u
    real(wp), allocatable :: y(:), yp(:)
    real(wp) :: e, h
    logical :: step
        step = .false.
        do while(.not. step)
            h = min(u%h, u%x_end - u%x)
            yp = u%yp
            y = rk45(u%f, u%x, u%y, h, e, yp)
        
            if( e>2d0*u%e .and. u%e>0 ) then
                u%h = h/2
                cycle
            elseif(e<0.5d0*u%e .and. u%e>0) then
                u%x = u%x + h
                u%y = y
                u%yp = yp
                h = h * 2d0
            else
                u%x = u%x + h            
                u%y = y
                u%yp = yp
            end if
            u%h = h
            step = .true.
        end do        
        u%yp = u%f(u%x, u%y)
    end subroutine
    
    function rk45(f,t,x,h,e,xp) result(x1)
    ! https://web.ma.utexas.edu/CNA/cheney-kincaid/f90code/CHP10/rk45ad.f90    
    !
    ! Numerical Mathematics and Computing, Fifth Edition
    ! Ward Cheney & David Kincaid
    ! Brooks/Cole Publ. Co.
    ! Copyright (c) 2003.  All rights reserved.
    ! For educational use with the Cheney-Kincaid textbook.
    ! Absolutely no warranty implied or expressed.
    !
    ! Section 10.3
    !
    ! File: rk45ad.f90
    !
    ! Adaptive scheme based on Runge-Kutta-Fehlberg method (rk45ad,rk45,f)
    procedure(f_rate), intent(in), pointer :: f
    real(wp), intent(in):: t, x(:), h
    real(wp), intent(out):: e
    real(wp), intent(inout):: xp(size(x))
    real(wp) :: x1(size(x))
    real(wp), parameter :: c20=0.25, c21=0.25, c30=0.375, c31=0.09375, c32=0.28125
    real(wp), parameter :: c40=0.92307692307692
    real(wp), parameter :: c41=0.87938097405553, c42=-3.2771961766045, c43=3.3208921256258
    real(wp), parameter :: c51=2.0324074074074, c52=-8.0, c53=7.1734892787524
    real(wp), parameter :: c54=-0.20589668615984
    real(wp), parameter :: c60=0.5, c61=-0.2962962962963, c62=2.0
    real(wp), parameter :: c63=-1.3816764132554, c64=0.45297270955166, c65=-0.275
    real(wp), parameter :: a1=0.11574074074074, a2=0, a3=0.54892787524366
    real(wp), parameter :: a4=0.5353313840156, a5=-0.2
    real(wp), parameter :: b1=0.11851851851852, b2=0.0, b3=0.51898635477583
    real(wp), parameter :: b4=0.50613149034201, b5=-0.18
    real(wp), parameter :: b6=0.036363636363636
    real(wp), dimension(size(x)) :: f1,f2,f3,f4,f5,f6,x5
        f1 = h*xp
        f2 = h*f(t+ c20*h,x + c21*f1)
        f3 = h*f(t+ c30*h,x + c31*f1 + c32*f2)
        f4 = h*f(t+ c40*h,x + c41*f1 + c42*f2 + c43*f3)
        f5 = h*f(t+h,x + c51*f1 + c52*f2 + c53*f3 + c54*f4)
        f6 = h*f(t+ c60*h,x + c61*f1 + c62*f2 + c63*f3 + c64*f4 + c65*f5)
        x5 = x + b1*f1 + b3*f3 + b4*f4 + b5*f5 + b6*f6
        x1  = x + a1*f1 + a3*f3 + a4*f4 + a5*f5
        e = maxval( abs(x1 - x5) )
        xp = f(t+h, x1)
    end function

    end module