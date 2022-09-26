    module mod_screws
    use mod_vectors
    use mod_show
    implicit none


    real(wp), dimension(3), parameter :: gravity = -10*j_

    type :: rigidbody
        real(wp) :: mass
        real(wp) :: mmoi(3)
        real(wp) :: cg(3)
        real(wp), allocatable :: dims(:)
    end type


    interface twist
    module procedure :: s_twist_at, s_twist_pure, s_twist_make
    end interface

    interface wrench
    module procedure :: s_wrench_at, s_wrench_pure, s_wrench_make
    end interface

    contains

    pure function mmoi_from_local(mmoi, q, inverse) result(I)
    real(wp), intent(in) :: mmoi(3), q(4)
    real(wp) :: I(3,3)
    logical, intent(in), optional :: inverse
    real(wp) :: I_body(3), R(3,3)
    if( present(inverse) .and. inverse) then
        I_body = 1.0_wp / mmoi
    else
        I_body = mmoi
    end if
    R = q_to_rotation(q)
    I(1,:) = I_body(1) * R(:,1)
    I(2,:) = I_body(2) * R(:,2)
    I(3,:) = I_body(3) * R(:,3)
    I = matmul(R, I)
    end function

    pure function s_twist_at(value, position, pitch) result(s)
    real(wp), intent(in) :: value(3), position(3), pitch
    real(wp) :: s(6)
    s = [ cross(position, value) + pitch * value, value]
    end function
    pure function s_twist_pure(value) result(s)
    real(wp), intent(in) :: value(3)
    real(wp) :: s(6)
    s = [ value, o_]
    end function
    pure function s_twist_make(value, moment) result(s)
    real(wp), intent(in) :: value(3), moment(3)
    real(wp) :: s(6)
    s = [ moment, value]
    end function
    pure function s_wrench_at(value, position, pitch) result(s)
    real(wp), intent(in) :: value(3), position(3), pitch
    real(wp) :: s(6)
    s = [ value, cross(position, value) + pitch * value]
    end function
    pure function s_wrench_pure(value) result(s)
    real(wp), intent(in) :: value(3)
    real(wp) :: s(6)
    s = [ o_, value]
    end function
    pure function s_wrench_make(value, moment) result(s)
    real(wp), intent(in) :: value(3), moment(3)
    real(wp) :: s(6)
    s = [ value, moment]
    end function

    pure function s_block(a,b,c,d) result(s)
    real(wp), intent(in), dimension(3,3) :: a,b,c,d
    real(wp), dimension(6,6) :: s
    integer :: i
    do i=1,3
        s(:,i) =[ a(:,i), c(:,i) ]
        s(:,i+3) =[ b(:,i), d(:,i) ]
    end do
    end function

    pure function s_inertia_matrix(mass, mmoi, q, cg, inverse) result(I)
    real(wp), intent(in) :: mass, mmoi(3), cg(3), q(4)
    logical, intent(in), optional :: inverse
    real(wp) :: cx(3,3), cxcx(3,3), IC(3,3), I(6,6)
    cx = v_cross_op(cg)
    if( present(inverse) .and. inverse) then
        IC = mmoi_from_local(mmoi, q, inverse)
        I = s_block(eye_/mass - matmul( cx, matmul(IC, cx)), matmul(cx, IC), -matmul(IC, cx), IC)
    else
        IC = mmoi_from_local(mmoi, q)
        cxcx = matmul(cx,cx)
        I = s_block( mass*eye_, -mass*cx, mass*cx, IC - mass*cxcx)
    end if
    end function

    pure function s_twist_cross_twist(t, s) result(u)
    real(wp), intent(in) :: t(6), s(6)
    real(wp) :: u(6), v(3), w(3), a(3), b(3)
    v = t(1:3)
    w = t(4:6)
    a = s(1:3)
    b = s(4:6)
    u = [ v_cross(w,a) + v_cross(v,b), v_cross(w, b) ]
    end function

    pure function s_twist_cross_wrench(t, h) result(u)
    real(wp), intent(in) :: t(6), h(6)
    real(wp) :: u(6), v(3), w(3), a(3), b(3)
    v = t(1:3)
    w = t(4:6)
    a = h(1:3)
    b = h(4:6)
    u = [ v_cross(w,a), v_cross(v,a) + v_cross(w, b) ]
    end function


    end module