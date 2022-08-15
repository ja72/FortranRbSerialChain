    !*******************************************************
    !*    LU decomposition routines used by test_lu.f90    *
    !*                                                     *
    !*                 F90 version by J-P Moreau, Paris    *
    !*                        (www.jpmoreau.fr)            *
    !* --------------------------------------------------- *
    !* Reference:                                          *
    !*                                                     *
    !* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
    !*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
    !*  University Press, 1986" [BIBLI 08].                *
    !*                                                     *
    !*******************************************************
    MODULE mod_lu
    use mod_vectors
    implicit none
    
    private LUDCMP, LUBKSB
        
    type lu_matrix
        real(wp), allocatable :: data(:,:)
        integer, allocatable :: indx(:)
        integer :: d, n
    contains
        procedure :: solve => lu_solve
        procedure :: inverse => lu_inverse
        procedure :: determinant => lu_determinant
    end type
    
    interface lu_matrix
        procedure :: lu_from_matrix
    end interface

    CONTAINS
    
    function eye(n) result(u)
    integer, intent(in) :: n
    real(wp) :: u(n,n)
    integer :: i
        u = 0d0
        forall(i=1:n)
            u(i,i)=1d0
        end forall    
    end function    
    
    function lu_from_matrix(A) result(lu)
    real(wp), intent(in) :: A(:,:)
    type(lu_matrix) :: lu
    integer :: ierr
        lu%data = A
        lu%n = size(A,1)
        call LUDCMP(lu%data, lu%indx, lu%d, ierr)
        if( ierr==1) then
            return
        end if        
    end function
    
    function lu_solve(lu,b) result(x)
    class(lu_matrix), intent(in) :: lu
    real(wp), intent(in) :: b(:)
    real(wp), allocatable :: x(:)
        x =  b
        call LUBKSB(lu%data, lu%indx, x)
    end function
    
    function lu_inverse(lu) result(A_inv)
    class(lu_matrix), intent(in) :: lu
    real(wp), allocatable :: A_inv(:,:)
    integer :: i
        A_inv = eye(lu%n)
        do i=1, lu%n
            call LUBKSB(lu%data, lu%indx, A_inv(:,i))
        end do
    end function
    
    function lu_determinant(lu) result(d)
    class(lu_matrix), intent(in) :: lu
    real(wp) :: d
    integer :: i
        d = lu%d        
        do i=1, lu%n
            d = d * lu%data(i,i)
        end do
    end function

    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the LU *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. INDX is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
    Subroutine LUDCMP(A,INDX,D,ierr)
    real(wp), intent(inout), allocatable :: A(:,:)
    integer, intent(out), allocatable :: indx(:)    
    integer, intent(out) :: d, ierr
    integer :: i, j, k, imax, n
    real(wp) ::  AMAX,DUM, s
    real(wp), allocatable ::  VV(:)
    n = size(A,1)
    if( n /= size(A,2) ) then
        return
    end if
    allocate(indx(n))
    allocate(vv(n))
    D=1; ierr=0    
    DO I=1,N
        AMAX=0.d0
        DO J=1,N
            IF (ABS(A(I,J)) > AMAX) AMAX=ABS(A(I,J))
        END DO ! j loop
        IF(AMAX < TINY) THEN
            ierr = 1
            RETURN
        END IF
        VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
        DO I=1,J-1
            s = A(I,J)
            DO K=1,I-1
                s = s - A(I,K)*A(K,J)
            END DO ! k loop
            A(I,J) = s
        END DO ! i loop
        AMAX = 0.d0
        DO I=J,N
            s = A(I,J)
            DO K=1,J-1
                s = s - A(I,K)*A(K,J)
            END DO ! k loop
            A(I,J) = s
            DUM = VV(I)*ABS(s)
            IF(DUM >= AMAX) THEN
                IMAX = I
                AMAX = DUM
            END IF
        END DO ! i loop

        IF(J /= IMAX) THEN
            DO K=1,N
                DUM = A(IMAX,K)
                A(IMAX,K) = A(J,K)
                A(J,K) = DUM
            END DO ! k loop
            D = -D
            VV(IMAX) = VV(J)
        END IF

        INDX(J) = IMAX
        IF(ABS(A(J,J)) < TINY) A(J,J) = TINY

        IF(J.NE.N) THEN
            DUM = 1.d0 / A(J,J)
            DO I=J+1,N
                A(I,J) = A(I,J)*DUM
            END DO ! i loop
        END IF
    END DO ! j loop

    RETURN
    END subroutine LUDCMP


    !  ******************************************************************
    !  * Solves the set of N linear equations A . X = B.  Here A is     *
    !  * input, not as the matrix A but rather as its LU decomposition, *
    !  * determined by the routine LUDCMP. INDX is input as the permuta-*
    !  * tion vector returned by LUDCMP. B is input as the right-hand   *
    !  * side vector B, and returns with the solution vector X. A, N and*
    !  * INDX are not modified by this routine and can be used for suc- *
    !  * cessive calls with different right-hand sides. This routine is *
    !  * also efficient for plain matrix inversion.                     *
    !  ******************************************************************
    Subroutine LUBKSB(A,INDX,B)
    integer, intent(in) :: indx(:)
    real(wp), intent(in) :: A(:,:)
    real(wp), intent(inout) :: B(:)
    real(wp)  :: s
    INTEGER :: i,j,ll,ii,n

    n = size(A,1)
    if( n /= size(A,2) .or. n/=size(indx) .or. n/=size(b) ) then
        return
    end if
    II = 0

    DO I=1,N
        LL = INDX(I)
        s = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
            DO J=II,I-1
                s = s - A(I,J)*B(J)
            END DO ! j loop
        ELSE IF(s.NE.0.d0) THEN
            II = I
        END IF
        B(I) = s
    END DO ! i loop

    DO I=N,1,-1
        s = B(I)
        IF(I < N) THEN
            DO J=I+1,N
                s = s - A(I,J)*B(J)
            END DO ! j loop
        END IF
        B(I) = s / A(I,I)
    END DO ! i loop

    RETURN
    END subroutine LUBKSB

    END MODULE 
