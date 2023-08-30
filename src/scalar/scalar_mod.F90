MODULE scalar_mod
    ! Do not 'USE' core_mod or ib_mod here, because they will get exported
    ! (this module had no 'PRIVATE' default)
    USE timeintegrate_scalar_mod
    USE scacore_mod

    IMPLICIT NONE(type, external)

    ! PRIVATE ::

CONTAINS
    SUBROUTINE init_scalar()
        USE core_mod, ONLY: dread, dcont
        USE itinfo_scalar_mod, ONLY: init_itinfo_scalar
        USE scastat_mod, ONLY: init_scastat

        CALL init_scacore()
        IF (.NOT. has_scalar) RETURN

        CALL init_itinfo_scalar(dcont)

        IF (.NOT. dread) THEN
            CALL init_tfield()
        END IF

        CALL init_scastat()
    END SUBROUTINE init_scalar


    SUBROUTINE finish_scalar
        USE itinfo_scalar_mod, ONLY: finish_itinfo_scalar
        USE scastat_mod, ONLY: finish_scastat

        IF (.NOT. has_scalar) RETURN

        CALL finish_scastat()
        CALL finish_itinfo_scalar()
        CALL finish_scacore()
    END SUBROUTINE finish_scalar


    SUBROUTINE init_tfield()
        USE core_mod
        USE ib_mod
        USE setboundarybuffers_scalar_mod, ONLY: setboundarybuffers_scalar

        TYPE(field_t), POINTER :: t
        INTEGER(intk) :: l, ilevel

        DO l = 1, nsca
            CALL get_field(t, scalar(l)%name)

            ! Write boundary conditions into buffers
            DO ilevel = minlevel, maxlevel
                CALL setboundarybuffers_scalar%bound(ilevel, t)
            END DO

            ! TODO: set initial condition

            CALL zero_ghostlayers(t)

            DO ilevel = minlevel, maxlevel
                CALL connect(ilevel, 2, s1=t%arr, corners=.TRUE.)
            END DO

            DO ilevel = maxlevel, minlevel+1, -1
                CALL ftoc(ilevel, t%arr, t%arr, 'T')
            END DO
        END DO
    END SUBROUTINE init_tfield
END MODULE scalar_mod
