MODULE scacore_mod
    USE core_mod
    USE ib_mod, ONLY: ftoc, parent, ib
    USE flow_mod, ONLY: has_flow

    IMPLICIT NONE(type, external)
    PRIVATE

    TYPE :: scalar_t
        CHARACTER(len=nchar_name) :: name
        REAL(realk) :: prmol
        INTEGER(intk) :: kayscrawford
    CONTAINS
        PROCEDURE :: prt
    END TYPE scalar_t

    ! Control parameters
    LOGICAL, PROTECTED :: has_scalar = .FALSE.
    LOGICAL, PROTECTED :: solve_scalar = .FALSE.
    REAL(realk) :: prturb

    ! Scalar/physical paramters
    INTEGER(intk), PROTECTED :: nsca
    TYPE(scalar_t), ALLOCATABLE, PROTECTED :: scalar(:)

    PUBLIC :: init_scacore, finish_scacore, scalar_t, scalar, nsca, prturb, &
      has_scalar, solve_scalar, maskbt

CONTAINS
    SUBROUTINE init_scacore()
        USE blockbt_mod
        ! Subroutine arguments
        ! None...

        ! Local variables
        TYPE(config_t) :: scaconf, sc
        TYPE(field_t), POINTER :: t, bt
        INTEGER(intk) :: l
        CHARACTER(len=nchar_name) :: fieldname
        CHARACTER(len=64) :: jsonptr
        LOGICAL :: kayscrawford

        ! Read configuration values - if not exists no timeintegration is
        ! performed
        has_scalar = .FALSE.
        IF (.NOT. fort7%exists("/scalar")) THEN
            IF (myid == 0) THEN
                WRITE(*, '("NO SCALAR")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF
        has_scalar = .TRUE.

        ! Scalar is dependent on flow to de defined (but not neccesarily solved)
        ! to have variables like gmol and rho defined
        IF (.NOT. has_flow) THEN
            WRITE(*, *) "Scalar needs flow!"
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Required values
        CALL fort7%get(scaconf, "/scalar")
        CALL scaconf%get_size("/scalars", nsca)
        ALLOCATE(scalar(nsca))

        ! Read parameters per scalar
        DO l = 1, nsca
            WRITE(jsonptr, '("/scalars/", I0)') l-1
            CALL scaconf%get(sc, jsonptr)

            CALL sc%get_value("/name", scalar(l)%name)
            ! To allow for field "{name}_AVG" to fit in nchar_name characters
            IF (LEN_TRIM(scalar(l)%name) + 4 > nchar_name) THEN
                CALL errr(__FILE__, __LINE__)
            END IF

            CALL sc%get_value("/prmol", scalar(l)%prmol)

            ! The parameters.json should have a logical, but it is easier
            ! to integrate computations with an integer
            CALL sc%get_value("/kayscrawford", kayscrawford, .FALSE.)
            scalar(l)%kayscrawford = l_to_i(kayscrawford)

            CALL sc%finish()
        END DO

        ! Optional values
        CALL scaconf%get_value("/solve", solve_scalar, .TRUE.)
        CALL scaconf%get_value("/prturb", prturb, 1.0)

        ! Declare fields
        DO l = 1, nsca
            CALL set_field(scalar(l)%name, dread=dread, required=dread, &
                dwrite=dwrite, buffers=.TRUE.)

            ! Get field, set PRMOL and scalar index as attribute
            CALL get_field(t, scalar(l)%name)
            CALL t%set_attr(scalar(l)%prmol, "PRMOL")
            CALL t%set_attr(l, "SCAIDX")

            ! For RK time integration
            WRITE(fieldname, '("D", A)') TRIM(scalar(l)%name)
            CALL set_field(fieldname)
        END DO

        ! Compute BT field from BP
        CALL set_field("BT", dwrite=.TRUE.)
        CALL get_field(bt, "BT")
        CALL blockbt(bt)
    END SUBROUTINE init_scacore


    SUBROUTINE finish_scacore
        IF (ALLOCATED(scalar)) DEALLOCATE(scalar)
    END SUBROUTINE finish_scacore


    PURE ELEMENTAL REAL(realk) FUNCTION prt(this, gtgmol)
        ! Calculation of local turbulent Prandtl number
        !
        ! limitations of Kays/Crawford: 0.5 < prmol < 7,
        !                               any Re
        !                               any dp/dx

        ! Kays/Crawford provides a relatively high prturb
        ! near the wall (in the sublayer) but approaches
        ! 0.85 as y+ increases into the log layer

        ! function arguments
        CLASS(scalar_t), INTENT(in) :: this
        REAL(realk), INTENT(IN) :: gtgmol

        ! Local variables
        REAL(realk) :: kayscrawford
        INTEGER(intk) :: switch

        ! Kays/Crawford solution
        IF (gtgmol > 0.0) THEN
            kayscrawford = 0.5882 + 0.228*gtgmol &
                - 0.0441*gtgmol**2*(1.0 - exp(-5.165/gtgmol))
        ELSE
            kayscrawford = this%prmol
        ENDIF

        ! Switch between Kays/Crawford solution and prturb
        switch = this%kayscrawford
        prt = switch*kayscrawford + (1-switch)*prturb
    END FUNCTION prt


    SUBROUTINE maskbt(t_f)
        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: t_f

        ! Local variables
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        TYPE(field_t), POINTER :: bt_f
        REAL(realk), POINTER, CONTIGUOUS :: t(:, :, :), bt(:, :, :)

        CALL get_field(bt_f, "BT")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL bt_f%get_ptr(bt, igrid)
            CALL t_f%get_ptr(t, igrid)

            CALL maskbt_grid(kk, jj, ii, t, bt)
        END DO
    END SUBROUTINE maskbt


    PURE SUBROUTINE maskbt_grid(kk, jj, ii, t, bt)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: t(kk, jj, ii)
        REAL(realk), INTENT(in) :: bt(kk, jj, ii)

        ! Local variables
        INTEGER :: k, j, i

        ! TODO: Indices?
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    t(k, j, i) = t(k, j, i)*bt(k, j, i)
                END DO
            END DO
        END DO
    END SUBROUTINE maskbt_grid
END MODULE scacore_mod
