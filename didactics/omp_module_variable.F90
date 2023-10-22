
MODULE var_mod
    IMPLICIT NONE(type, external)
    integer :: alpha
    integer :: beta
END MODULE var_mod


PROGRAM main

    use omp_lib
    use var_mod
    implicit none

    integer, parameter :: len = 100

    integer :: a(len), a_copy(len)
    integer :: b(len)
    integer :: i, n

    a = 1
    b = 2

    CALL omp_set_num_threads(2)

!$omp parallel private( alpha )

    DO n = 1, 10

!$omp single
        IF ( n > 2 ) THEN
            IF ( alpha /= omp_get_thread_num() ) WRITE(*,*) "ERROR alpha"
            IF ( beta /= omp_get_thread_num()+1 ) WRITE(*,*) "ERROR beta"
        END IF
        DO i = 1, len
            a(i) = a(i) + b(len+1-i)
        END DO
        DO i = len, 1, (-1)
            a(i) = a(i) - b(len+1-i)
        END DO
        alpha = omp_get_thread_num()
        beta = omp_get_thread_num() + 1
!$omp end single

!$omp single
        IF ( n > 2 ) THEN
            IF ( alpha /= omp_get_thread_num() ) WRITE(*,*) "ERROR alpha"
            IF ( beta /= omp_get_thread_num()+1 ) WRITE(*,*) "ERROR beta"
        END IF
        DO i = 1, len
            a(i) = a(i) + b(len+1-i)
        END DO
        DO i = len, 1, (-1)
            a(i) = a(i) - b(len+1-i)
        END DO
        alpha = omp_get_thread_num()
        beta = omp_get_thread_num() + 1
!$omp end single
!$omp barrier

    END DO

!$omp end parallel

    WRITE(*,*) a
    WRITE(*,*) b

END PROGRAM main
