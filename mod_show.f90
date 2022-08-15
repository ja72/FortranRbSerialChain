    module mod_show
    use, intrinsic :: iso_fortran_env
    implicit none
        
    interface show
        procedure :: show_vector_i, show_matrix_i
        procedure :: show_vector_r, show_vector_d
        procedure :: show_matrix_r, show_matrix_d
    end interface
    

    contains

    subroutine show_matrix_r(a)
    real(real32), intent(in) :: a(:,:)
    integer i,n
        n = size(a,1)
        do i=1,n
            print '(*(1x,g12.4))', a(i,:)
        end do
    end subroutine
    
    subroutine show_matrix_d(a)
    real(real64), intent(in) :: a(:,:)
        call show_matrix_r(real(a))
    end subroutine

    subroutine show_vector_r(a)
    real(real32), intent(in) :: a(:)
    integer i,n
        n = size(a,1)
        do i=1,n
            print '(*(1x,g12.4))', a(i)
        end do
    end subroutine
    
    subroutine show_vector_d(a)
    real(real64), intent(in) :: a(:)
        call show_vector_r(real(a))
    end subroutine
    
    subroutine show_vector_i(a)
    integer, intent(in) :: a(:)
    integer i,n
        n = size(a,1)
        do i=1,n
            print '(*(1x,i12))', a(i)
        end do
    end subroutine
    subroutine show_matrix_i(a)
    integer, intent(in) :: a(:,:)
    integer i,n
        n = size(a,1)
        do i=1,n
            print '(*(1x,i12))', a(i,:)
        end do
    end subroutine


    end module