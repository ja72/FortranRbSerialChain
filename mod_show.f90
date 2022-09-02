    module mod_show
    use, intrinsic :: iso_fortran_env
    implicit none
        
    interface show
        procedure :: show_scalar_i, show_vector_i, show_matrix_i
        procedure :: show_scalar_r, show_vector_r, show_matrix_r
        procedure :: show_scalar_d, show_vector_d, show_matrix_d
    end interface
    

    contains

    subroutine show_matrix_r(t,a)
    character(len=*) :: t
    real(real32), intent(in) :: a(:,:)
    integer i,n
        n = size(a,1)
        print *, t
        do i=1,n
            print '(*(1x,g12.4))', a(i,:)
        end do
    end subroutine
    
    subroutine show_matrix_d(t, a)
    character(len=*) :: t    
    real(real64), intent(in) :: a(:,:)
    call show_matrix_r(t, real(a))
    end subroutine

    subroutine show_vector_r(t,a)
    character(len=*) :: t
    real(real32), intent(in) :: a(:)
    integer i,n
        n = size(a,1)
        print *, t
        do i=1,n
            print '(*(1x,g12.4))', a(i)
        end do
    end subroutine
    
    subroutine show_vector_d(t, a)
    character(len=*) :: t
    real(real64), intent(in) :: a(:)
    call show_vector_r(t, real(a))
    end subroutine
    
    subroutine show_vector_i(t,a)
    character(len=*) :: t
    integer, intent(in) :: a(:)
    integer i,n
        n = size(a,1)
        print *, t
        do i=1,n
            print '(*(1x,i12))', a(i)
        end do
    end subroutine
    subroutine show_matrix_i(t,a)
    character(len=*) :: t
    integer, intent(in) :: a(:,:)
    integer i,n
        n = size(a,1)
        print *, t
        do i=1,n
            print '(*(1x,i12))', a(i,:)
        end do
    end subroutine

    subroutine show_scalar_i(t,a)
    character(len=*) :: t
    integer, intent(in) :: a
        print *, t
        print '(*(1x,i12))', a
    end subroutine
    subroutine show_scalar_r(t,a)
    character(len=*) :: t
    real(real32), intent(in) :: a
        print *, t
        print '(*(1x,g12.4))', a
    end subroutine
    subroutine show_scalar_d(t,a)
    character(len=*) :: t
    real(real64), intent(in) :: a
    call show_scalar_r(t, real(a))
    end subroutine

    end module