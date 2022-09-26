    module mod_joints
    use mod_screws
    implicit none
    
    enum, bind(c)
        enumerator :: revolute_joint
        enumerator :: prismatic_joint
    end enum
    
    enum, bind(c)
        enumerator :: kinematic_joint
        enumerator :: dynamic_joint
    end enum
    
    type :: joint
        integer :: type
        real(wp) :: pos(3), axis(3), pitch        
        type(rigidbody) :: body
    end type
    
    interface joint
        module procedure :: jt_new_joint
    end interface
    
    type :: state
        integer :: status
        real(wp):: q, qp, qpp, tau
    end type
    
    interface state
        module procedure :: st_new_state
    end interface
    
    type :: kinematics
        real(wp) :: pos(3), ori(4), cg(3)
        real(wp) :: spi(6,6)
        real(wp) :: s(6), v(6), k(6), p(6)
        real(wp) :: a(6), f(6)
    end type

    type :: articulated
        real(wp) :: spi(6,6), RU(6,6)
        real(wp) :: p(6), T(6)
    end type
    
    
    contains

    function jt_new_joint(rb, type, pos, axis, pitch) result(jt)
    type(rigidbody), intent(in), target :: rb
    type(joint) :: jt
    integer, intent(in) :: type
    real(wp), intent(in) :: pos(3), axis(3)
    real(wp), intent(in), optional :: pitch
        jt%type = type
        jt%pos = pos
        jt%axis = axis
        if( present(pitch) ) then
            jt%pitch = pitch
        else
            jt%pitch = 0.0_wp
        end if
        jt%body = rb
    end function
    
    pure function st_new_state(status, q, qp, known) result(st)
    type(state) :: st
    integer, intent(in) :: status
    real(wp), intent(in), optional :: q, qp, known
        
        st%status = status
        st%q = 0.0_wp
        st%qp = 0.0_wp
        st%qpp = 0.0_wp
        st%tau = 0.0_wp
        if(present(q)) then
            st%q = q
        end if
        if(present(qp)) then
            st%qp = qp
        end if
        if(present(known)) then
            select case(status)
            case(kinematic_joint)
                st%qpp = known
            case(dynamic_joint)                
                st%tau = known                
            end select            
        end if
    end function
        
    pure function solve_chain(jt, st) result(dyn)
    type(joint), intent(in) :: jt(:)
    type(state), intent(in) :: st(size(jt))
    type(state) :: dyn(size(jt))
    type(kinematics) :: kin(size(jt))
    type(articulated) :: art(size(jt))
    real(wp) :: a_prev(6), jti
    integer :: i, n
    n = size(jt)
    
        kin = calc_kinematics(st)
        art = calc_articulated(kin,st)
        
        a_prev = 0.0_wp
        do i=1,n
            dyn(i) = st(i)
            select case(st(i)%status)
            case(kinematic_joint)                
                kin(i)%a = a_prev + kin(i)%s*st(i)%qpp + kin(i)%k
                kin(i)%f = matmul(art(i)%spi, kin(i)%a) + art(i)%p
                dyn(i)%tau = dot_product(kin(i)%s, kin(i)%f)
            case(dynamic_joint)                        
                jti = dot_product(kin(i)%s, matmul( art(i)%spi, kin(i)%s))
                dyn(i)%qpp = ( st(i)%tau - dot_product(kin(i)%s, matmul(art(i)%spi, a_prev + kin(i)%k) + art(i)%p))/jti
                kin(i)%a = a_prev + kin(i)%s*dyn(i)%qpp + kin(i)%k
                kin(i)%f = matmul(art(i)%spi, kin(i)%a) + art(i)%p
            end select
            a_prev = kin(i)%a
        end do
        
    contains
    
        pure function calc_kinematics(st) result(kin)
        type(state), intent(in) :: st(size(jt))
        type(kinematics) :: kin(size(jt))
        type(rigidbody) :: rb        
        real(wp) :: r_prev(3), q_prev(4), v_prev(6)
        real(wp) :: w_i(6)
        integer :: i
        
        r_prev = 0.0_wp                                                        
        q_prev = q_eye
        v_prev = 0.0_wp
        do i=1, n
            rb = jt(i)%body            
            select case(jt(i)%type)
            case(revolute_joint)
                kin(i)%pos = r_prev + rot(q_prev, jt(i)%pos + jt%pitch * jt(i)%axis * st(i)%q)
                kin(i)%ori = q_prev .o. rot(jt(i)%axis, st(i)%q)
            case(prismatic_joint)
                kin(i)%pos = r_prev + rot(q_prev, jt(i)%pos + jt(i)%axis * st(i)%q)
                kin(i)%ori = q_prev 
            end select
            kin(i)%cg = kin(i)%pos + rot(kin(i)%ori, rb%cg )
            kin(i)%spi = s_inertia_matrix(rb%mass,rb%mmoi,kin(i)%ori,kin(i)%cg)
            select case(jt(i)%type)
            case(revolute_joint)
                kin(i)%s = twist( rot(q_prev, jt(i)%axis), kin(i)%pos, jt(i)%pitch)
            case(prismatic_joint)
                kin(i)%s = twist( rot(q_prev, jt(i)%axis) )
            end select
            kin(i)%v = v_prev + kin(i)%s*st(i)%qp
            kin(i)%k = s_twist_cross_twist( v_prev, kin(i)%s*st(i)%qp)
            w_i = wrench( rb%mass*gravity, kin(i)%cg, 0.0_wp)
            kin(i)%p = s_twist_cross_wrench(kin(i)%v, matmul( kin(i)%spi, kin(i)%v)) - w_i
                        
            r_prev = kin(i)%pos
            q_prev = kin(i)%ori
            v_prev = kin(i)%v            
        end do
        
        end function
        
        pure function calc_articulated(kin,st) result(art)
        type(articulated) :: art(size(jt))
        type(kinematics), intent(in) :: kin(size(jt))
        type(state), intent(in) :: st(size(jt))
        integer :: i
            art(n)%spi = kin(n)%spi
            art(n)%p = kin(n)%p
            art(n)%T = matmul(art(n)%spi,kin(n)%s)/dot_product(kin(n)%s, matmul(art(n)%spi,kin(n)%s))
            art(n)%RU = identity_matrix(6) - v_outer(art(n)%T, kin(n)%s)
                
            do i=n-1,1,-1
                art(i)%p = kin(i)%p + art(i+1)%T*st(i+1)%tau + matmul( art(i+1)%RU, matmul(art(i+1)%spi,kin(i+1)%k) + art(i+1)%p)
                art(i)%spi = kin(i)%spi + matmul( art(i+1)%RU, art(i+1)%spi )
                art(i)%T = matmul(art(i)%spi,kin(i)%s)/dot_product(kin(i)%s, matmul(art(i)%spi,kin(i)%s))
                art(i)%RU = identity_matrix(6) - v_outer(art(i)%T, kin(i)%s)            
            end do
        
        end function
            
    end function
    
        
    
    pure function jt_find_children(parents, index) result(children)
    integer, intent(in) :: parents(:), index
    integer , allocatable :: children(:)
    integer :: i, j, n
        n = size(parents)
        allocate(children(n))
        j = 0
        do i=1, n
            if( parents(i)==index) then
                j = j + 1
                children(j) = i
            end if
        end do
        children = children(1:j)
    end function

    
    end module