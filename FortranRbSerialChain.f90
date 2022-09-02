!****************************************************************************
!
!  PROGRAM: FortranRbSerialChain
!
!  PURPOSE:  Serial RB chain solver using articulated method with arrays
!
!****************************************************************************

    program FortranRbSerialChain
    use mod_show
    use mod_screws
    implicit none
        
    character(len=*), parameter :: fmt = '(a,*(1x,g0.7))'
    integer, parameter :: n_bodies = 6
    ! Variables
    integer :: i
    real(wp), dimension(n_bodies) :: q,qp,qpp,tau
    real(wp) :: length
    type(rigidbody),target :: rb
    type(joint) :: joints(n_bodies)
    !integer, allocatable :: parents(:), children(:)
        
    ! Body of FortranConsole2
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
            joints(i) = rb%joint(revolute_joint, o_, k_)
        else
            joints(i) = rb%joint(revolute_joint, length*i_, k_)
        end if
    end do
    
    q = 0.0_wp
    qp = 0.0_wp
    tau = 0.0_wp
    qpp = jt_acceleration(joints, q, qp, tau)
    
    print *, ''
    print *, " == Serial Chain Dynamics == "
        
    call show("q=", q)
    call show("qp=", qp)
    call show("tau=", tau)
    call show("qpp=", qpp)
    
    contains
    

    end program FortranRbSerialChain

