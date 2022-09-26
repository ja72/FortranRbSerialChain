!****************************************************************************
!
!  PROGRAM: FortranRbSerialChain
!
!  PURPOSE:  Serial RB chain solver using articulated method with arrays
!
!****************************************************************************

    program FortranRbSerialChain
    use mod_show
    use mod_joints
    implicit none
        
    character(len=*), parameter :: fmt = '(a,*(1x,g0.7))'
    integer, parameter :: n_bodies = 6
    ! Variables
    integer :: i
    type(state) :: st(n_bodies)
    real(wp) :: length
    type(rigidbody),target :: rb
    type(joint) :: joints(size(st))
    !integer, allocatable :: parents(:), children(:)
        
    ! Body of FortranRbSerialChain
    print *, 'Serial Chain with n=', n_bodies, 'bodies.'
    
    length = 0.25_wp
    rb%mass = 0.1_wp
    rb%mmoi = [0.0_wp, rb%mass/12*length**2, rb%mass/12*length**2]
    rb%cg = i_ * (length/2)
    rb%dims = [ length ]

    print *, "Body Properties"
    call show("length = ", length)
    call show("mass = ", rb%mass)
    call show("mmoi = ", rb%mmoi)
    call show("cg = ", rb%cg)
    
    do i=1,n_bodies
        if( i==1 ) then
            joints(i) = joint(rb, revolute_joint, o_, k_)
        else
            joints(i) = joint(rb, revolute_joint, length*i_, k_)
        end if
    end do
    st = [ (state(dynamic_joint, known=0.0_wp), i=1, n_bodies) ]
    st(1)%status = kinematic_joint
    st(1)%qpp = 50.0_wp
    
    st = solve_chain(joints, st)
    
    print *, ''
    print *, " == Serial Chain Dynamics == "
        
    call show("q=", st(:)%q)
    call show("qp=", st(:)%qp)
    call show("tau=", st(:)%tau)
    call show("qpp=", st(:)%qpp)
    
    contains
    

    end program FortranRbSerialChain

