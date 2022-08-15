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
    contains
        procedure :: kinematics => rb_kinematics
    end type
    
    type :: kinematics
        real(wp) :: q, qp, tau
        real(wp) :: pos(3), ori(4), cg(3)
        real(wp) :: spi(6,6)
        real(wp) :: s(6), v(6), k(6), p(6)
    end type
    
    type :: articulated
        real(wp) :: spi(6,6), RU(6,6)
        real(wp) :: p(6), T(6)
    end type
    
    type :: dynamics
        real(wp) :: qpp
        real(wp) :: a(6), f(6)
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
    
    pure function rb_kinematics(rb,q,qp,tau) result(kin)
    class(rigidbody), intent(in) :: rb
    real(wp), intent(in) :: q(:), qp(size(q)),tau(size(q))
    type(kinematics) :: kin(size(q))
    real(wp) :: length, r_prev(3), q_prev(4), v_prev(6)
    real(wp) :: w_i(6)    
    integer :: i, n
        n = size(q)
        
        !tex: Recursive Kinematics
        ! $$\begin{aligned}
        !   \vec{r}_i & = \vec{r}_{i-1} + R_{i-1} \vec{\ell} \\
        !   \vec{R}_i & = \vec{R}_{i-1} \mathrm{rot}(\hat{k},q_i) \\
        !   \vec{cg}_i & = \vec{r}_i + R_i \vec{cg}_{\rm body} \\
        !   \bf{I}_i & = \mathrm{spi}(m, I_{\rm body}, R_i, \vec{cg}_i) \\
        !   \boldsymbol{s}_i &= \mathrm{twist}(\hat{k}, \vec{r}_i, 0) \\
        !   \boldsymbol{v}_i &= \boldsymbol{v}_{i-1} + \boldsymbol{s}_i \dot{q}_i \\
        !   \boldsymbol{\kappa}_i &= \boldsymbol{v}_i \times \boldsymbol{s}_i \dot{q}_i \\
        !   \boldsymbol{w}_i &= \mathrm{wrench}(m\,\vec{g}, \vec{cg}, 0) \\
        !   \boldsymbol{p}_i &= \boldsymbol{v}_i \times \bf{I}_i \boldsymbol{v}_i - \boldsymbol{w}_i
        ! \end{aligned}$$
        
        length = rb%dims(1)
        r_prev = 0.0_wp                                                        
        q_prev = q_eye
        v_prev = 0.0_wp
        do i=1, n
            kin(i)%q = q(i)
            kin(i)%qp = qp(i)
            kin(i)%tau = tau(i)            
            kin(i)%pos = r_prev + matmul(rot(q_prev), i_ * length)
            kin(i)%ori = q_prev .o. rot(k_, q(i))
            kin(i)%cg = kin(i)%pos + matmul(rot(kin(i)%ori), rb%cg )
            kin(i)%spi = s_inertia_matrix(rb%mass,rb%mmoi,kin(i)%ori,kin(i)%cg)
            kin(i)%s = twist( k_, kin(i)%pos, 0.0_wp)
            kin(i)%v = v_prev + kin(i)%s*qp(i)
            kin(i)%k = s_twist_cross_twist( v_prev, kin(i)%s*qp(i))
            w_i = wrench( rb%mass*gravity, kin(i)%cg, 0.0_wp)
            kin(i)%p = s_twist_cross_wrench(kin(i)%v, matmul( kin(i)%spi, kin(i)%v)) - w_i
            
            r_prev = kin(i)%pos
            q_prev = kin(i)%ori
            v_prev = kin(i)%v            
        end do
    end function
    
    pure function kin_articulated(kin) result(art)
    type(kinematics), intent(in) :: kin(:)
    type(articulated) :: art(size(kin))
    integer :: i,n
        n = size(kin)
        
        !tex: Recursive Articulated Inertia
        !$$\begin{aligned}\boldsymbol{d}_{i} & =\boldsymbol{p}_{i}+\boldsymbol{T}_{i+1}Q_{i+1}+\boldsymbol{\Phi}_{i+1}\left({\bf A}_{i+1}\boldsymbol{\kappa}_{i+1}+\boldsymbol{d}_{i+1}\right)\\
        !{\bf A}_{i} & ={\bf I}_{i}+\boldsymbol{\Phi}_{i+1}{\bf A}_{i+1}\\
        !\boldsymbol{T}_{i} & ={\bf A}_{i}\boldsymbol{s}_{i}\left(\boldsymbol{s}_{i}^{\intercal}{\bf A}_{i}\boldsymbol{s}_{i}\right)^{-1}\\
        !\boldsymbol{\Phi}_{i} & =1-\boldsymbol{T}_{i}\boldsymbol{s}_{i}^{\intercal}
        !\end{aligned}$$

        
        art(n)%spi = kin(n)%spi
        art(n)%p = kin(n)%p
        art(n)%T = matmul(art(n)%spi,kin(n)%s)/dot_product(kin(n)%s, matmul(art(n)%spi,kin(n)%s))
        art(n)%RU = identity_matrix(6) - v_outer(art(n)%T, kin(n)%s)
                
        do i=n-1,1,-1
            art(i)%p = kin(i)%p + art(i+1)%T*kin(i+1)%tau + matmul( art(i+1)%RU, matmul(art(i+1)%spi,kin(i+1)%k) + art(i+1)%p)
            art(i)%spi = kin(i)%spi + matmul( art(i+1)%RU, art(i+1)%spi )
            art(i)%T = matmul(art(i)%spi,kin(i)%s)/dot_product(kin(i)%s, matmul(art(i)%spi,kin(i)%s))
            art(i)%RU = identity_matrix(6) - v_outer(art(i)%T, kin(i)%s)            
        end do
    end function
    
    pure function art_dynamics(kin, art) result(dyn)
    type(kinematics), intent(in) :: kin(:)
    type(articulated), intent(in) :: art(:)
    type(dynamics) :: dyn(size(kin))
    real(wp) :: a_prev(6), jti, tst(6), f_next(6), f_acc(6), f_err
    integer :: i,n
        n = size(kin)
        
        !tex: Recursive Dynamics
        !$$\begin{aligned}\ddot{q}_{i} & =\left(\boldsymbol{s}_{i}^{\intercal}{\bf A}_{i}\boldsymbol{s}_{i}\right)^{-1}\left(Q_{i}-\boldsymbol{s}_{i}^{\intercal}\left({\bf A}_{i}\left(\boldsymbol{a}_{i-1}+\boldsymbol{\kappa}_{i}\right)+\boldsymbol{d}_{i}\right)\right)\\
        !\boldsymbol{a}_{i} & =\boldsymbol{a}_{i-1}+\boldsymbol{s}_{i}\ddot{q}_{i}+\boldsymbol{\kappa}_{i}\\
        !\boldsymbol{f}_{i} & ={\bf A}_{i}\boldsymbol{a}_{i}+\boldsymbol{d}_{i}
        !\end{aligned}$$
        
        a_prev = 0.0_wp
        do i=1,n
            jti = dot_product(kin(i)%s, matmul( art(i)%spi, kin(i)%s))
            dyn(i)%qpp = ( kin(i)%tau - dot_product(kin(i)%s, matmul(art(i)%spi, a_prev + kin(i)%k) + art(i)%p))/jti
            dyn(i)%a = a_prev + kin(i)%s*dyn(i)%qpp + kin(i)%k
            dyn(i)%f = matmul(art(i)%spi, dyn(i)%a) + art(i)%p            
            a_prev = dyn(i)%a
        end do
        
        !tex: Check Force Balance
        ! $$\boldsymbol{f}_{i}={\bf I}_{i}\boldsymbol{a}_{i}+\boldsymbol{p}_{i}+\boldsymbol{f}_{i+1}$$
        
        f_next = 0.0_wp
        f_err = 0.0_wp
        do i=n,1,-1
            f_acc = matmul( kin(i)%spi, dyn(i)%a) + kin(i)%p
            tst = dyn(i)%f - f_acc - f_next
            f_err = max(f_err, maxval(abs(tst)))
            f_next = dyn(i)%f
        end do
        if( f_err > tiny ) then
            error stop "Force Balance Error."
        end if
    end function
    
    pure function dyn_chain_rate(rb,q,qp,tau) result(qpp)
    ! About: 
    !   Solves a dynamics problem of a serial chain of connected
    !   identical bodies. Uses the recursive articulated inertia method
    !   to solve for the joint accelerations.
    !
    ! Inputs:
    !   rb:     the rigid body properties for all bodies
    !   q:      vector of joint angles. # of bodies = size(q)
    !   qp:     vector of joint speeds.
    !   tau:    vector of joint torques
    !
    ! Outputs
    !   qpp:    vector of joint accelerations
    !
    ! References:
    !   URL:        https://homes.cs.washington.edu/~todorov/courses/amath533/FeatherstoneOrin00.pdf
    !   Title:      Robot Dynamics: Equations and Algorithms
    !   Authors:    Roy Featherstone, David Orin
    type(rigidbody), intent(in) :: rb
    real(wp), intent(in) :: q(:), qp(size(q)),tau(size(q))
    real(wp) :: qpp(size(q))
    type(kinematics), allocatable :: kin(:)
    type(articulated), allocatable :: art(:)
    type(dynamics), allocatable :: dyn(:)
    
        kin = rb_kinematics(rb, q, qp, tau)
        art = kin_articulated(kin)
        dyn = art_dynamics(kin, art)
        
        qpp = dyn(:)%qpp
    
    end function

    end module