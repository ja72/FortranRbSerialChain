    module mod_vectors
    use, intrinsic :: iso_fortran_env
    implicit none

    ! Variables
    integer, parameter :: sp = real32, wp = real64
    real(wp), parameter :: pi = acos(-1.0_wp), tiny = 1.0_wp/8796093022208

    real(wp), dimension(3), parameter :: &
        o_ = [0.0_wp, 0.0_wp, 0.0_wp], &
        i_ = [1.0_wp, 0.0_wp, 0.0_wp], &
        j_ = [0.0_wp, 1.0_wp, 0.0_wp], &
        k_ = [0.0_wp, 0.0_wp, 1.0_wp]
    real(wp), dimension(3, 3), parameter :: &
        zero_= reshape([0.0_wp,0.0_wp,0.0_wp, &
                        0.0_wp,0.0_wp,0.0_wp, &
                        0.0_wp,0.0_wp,0.0_wp], [3,3]), &
        eye_ = reshape([1.0_wp,0.0_wp,0.0_wp, &
                        0.0_wp,1.0_wp,0.0_wp, &
                        0.0_wp,0.0_wp,1.0_wp], [3,3])
    real(wp), dimension(4), parameter :: &
        q_zero = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], &
        q_eye = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]

    interface operator (.x.)
    module procedure :: v_cross
    end interface

    interface operator (.o.)
    module procedure :: q_product
    end interface

    interface cross
    module procedure :: v_cross, v_cross_op
    end interface

    interface rot
    module procedure :: q_to_rotation, q_axis_angle, q_rotate_vector
    end interface

    contains
    pure function identity_matrix(n) result(I)
    integer, intent(in) :: n
    real(wp) :: I(n,n)
    integer :: k
    I = 0.0_wp
    forall(k=1:n)
        I(k,k) = 1.0_wp
    end forall
    end function
    
    pure function v_outer(a,b) result(s)
    real(wp), intent(in) :: a(:),b(:)
    real(wp) :: s( size(a), size(b) )
    integer :: i,j
        forall (i=1:size(a))
          forall(j=1:size(b)) s(i,j) = a(i)*b(j)
        end forall
    end function

    pure function v_cross(v, u) result(s)
    real(wp), intent(in) :: v(3), u(3)
    real(wp) :: s(3)
    s = [ v(2)*u(3) - v(3)*u(2), v(3)*u(1) - v(1)*u(3), v(1)*u(2) - v(2)*u(1) ]
    end function

    pure function v_cross_op(v) result(s)
    real(wp), intent(in) :: v(3)
    real(wp) :: s(3,3)
    s = reshape( [0.0_wp, v(3), -v(2), -v(3), 0.0_wp, v(1), v(2), -v(1), 0.0_wp ], [3,3])
    end function

    pure function q_conjugate(q) result(p)
    real(wp), intent(in) :: q(4)
    real(wp) :: p(4)
    p = [ -q(1:3), q(4) ]
    end function

    pure function q_inverse(q) result(p)
    real(wp), intent(in) :: q(4)
    real(wp) :: p(4), qq
    qq = dot_product(q,q)
    p = [ -q(1:3)/qq, q(4)/qq ]
    end function

    pure function r_axis_angle(axis,theta) result(R)
    real(wp), intent(in) :: axis(3), theta
    real(wp) :: R(3,3), s, v, k(3), kx(3,3), kxkx(3,3)
    s = sin(theta)
    v = 1-cos(theta)
    k = axis / norm2(axis)
    kx = v_cross_op(k)
    kxkx = matmul(kx,kx)
    R = eye_ + s * kx + v * kxkx
    end function

    pure function q_axis_angle(axis,theta) result(q)
    real(wp), intent(in) :: axis(3), theta
    real(wp) :: q(4), k(3)
    k = axis / norm2(axis)
    q = [ sin(theta/2) * k, cos(theta/2) ]
    end function

    pure function q_to_rotation(q, inverse) result(R)
    real(wp), intent(in) :: q(4)
    real(wp) :: R(3,3), s, v(3), kx(3,3), kxkx(3,3)
    logical, intent(in), optional :: inverse
    s = q(4)
    v = q(1:3)
    kx = v_cross_op(v)
    kxkx = matmul(kx,kx)
    if( present(inverse) .and. inverse) then
        s = -s
    end if
    R = eye_ + 2*s*kx + 2*kxkx
    end function
    
    pure function q_rotate_vector(q, x) result(y)
    real(wp), intent(in) :: q(4), x(3)
    real(wp) :: y(3), kx(3), kxkx(3), s, v(3)
        s = q(4)
        v = q(1:3)
        kx = v .x. x
        kxkx = v .x. kx
        y = x + 2*s*kx + 2*kxkx
    end function

    pure function q_product(a,b) result(q)
    real(wp), intent(in) :: a(4), b(4)
    real(wp) :: q(4)
    q =  [a(4)*b(1) - a(3)* b(2) + a(2)* b(3) + a(1)* b(4),a(4)* b(2) + a(3)* b(1) + a(2)* b(4) - a(1)* b(3),a(4)* b(3) + a(3)* b(4) - a(2)* b(1) + a(1)* b(2),a(4)* b(4) - a(3)* b(3) - a(2)* b(2) - a(1)* b(1)]
    end function
    
    end module