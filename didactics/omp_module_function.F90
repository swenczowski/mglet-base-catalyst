

MODULE dummy_mod
    IMPLICIT NONE(type, external)
    integer, allocatable :: buffer(:)
    public :: add
CONTAINS
    SUBROUTINE add( inout, len, diff )
        integer, intent(in) :: len
        integer, intent(in) :: diff
        integer, intent(inout) :: inout(len)
        integer :: i,j
!$omp critical
        allocate( buffer(len) )
            WRITE(*,*) "critial entered"
            DO i = 1, len
                buffer(i) = inout(i) + diff
            END DO
            DO j = len, 1, (-1)
                buffer(j) = buffer(j) - diff
            END DO
            inout = buffer
            WRITE(*,*) "critial left"
        deallocate( buffer )
!$omp end critical
    END SUBROUTINE
END MODULE dummy_mod



PROGRAM main

    use omp_lib
    use dummy_mod
    implicit none

    integer, parameter :: len = 100

    integer :: a(len), a_copy(len)
    integer :: b(len)
    integer :: n, alpha

    a = 1
    b = 2

    CALL omp_set_num_threads(5)

!$omp parallel do num_threads(5)
    DO n = 1, 10
        alpha = omp_get_thread_num()
        CALL add( b, len, alpha )
    END DO
!$omp end parallel do

    WRITE(*,*) a
    WRITE(*,*) b

END PROGRAM main
