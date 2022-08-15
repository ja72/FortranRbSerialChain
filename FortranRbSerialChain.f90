!****************************************************************************
!
!  PROGRAM: FortranRbSerialChain
!
!  PURPOSE:  Serial RB chain solver using articulated method with arrays
!
!****************************************************************************

    program FortranConsole2
    use mod_show
    use mod_screws
    implicit none
    
    character(len=*), parameter :: fmt = '(a,*(1x,g0.7))'
    integer, parameter :: n_bodies = 6
    ! Variables
    real(wp), dimension(n_bodies) :: q,qp,qpp,tau
    real(wp) :: length
    type(rigidbody) :: rb

    ! Body of FortranConsole2
    print *, 'Serial Chain with n=', n_bodies, 'bodies.'
    
    length = 0.25_wp
    rb%mass = 0.1_wp
    rb%mmoi = [0.0_wp, rb%mass/12*length**2, rb%mass/12*length**2]
    rb%cg = i_ * (length/2)
    rb%dims = [ length ]

    print *, "Body Properties"
    print *, "length = ", length
    print *, "mass = ", rb%mass
    print *, "mmoi = ", rb%mmoi
    print *, "cg = ", rb%cg
    
    q = 0.0_wp
    qp = 0.0_wp
    tau = 0.0_wp
    qpp = dyn_chain_rate(rb, q, qp, tau)
    
    print *, ''
    print *, " == Serial Chain Dynamics == "
        
    print *, 'Joint Angles'
    call show(q)
    print *, 'Joint Speeds'
    call show(qp)
    print *, 'Joint Torques'
    call show(tau)
    print *, 'Joint Accelerations'
    call show(qpp)
    
    contains
    

    end program FortranConsole2

