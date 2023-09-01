#include "mg41.def"

MODULE render_mod

    USE iso_c_binding, ONLY: C_CHAR, C_INT

    USE precision_mod, ONLY: intk, realk, real64
    USE dimensions_mod, ONLY: idim3d, idim2d, idim1d, nscal

    USE fields_mod, ONLY: x, y, z, dx, dy, dz, ddx, ddy, ddz
    USE fields_mod, ONLY: u, v, w, p, t, hilf, bp
    USE statfields_mod, ONLY: au, av, aw, ap, at
    USE statfields_mod, ONLY: auum, avvm, awwm, auvm, auwm, avwm

    USE comms_mod, ONLY: myid, numprocs

    IMPLICIT NONE

    ! everything is private, unless explicitly declared public
    PRIVATE

    INTERFACE
        SUBROUTINE proto( p_mgdims, p_get_ip1, p_get_ip3, p_list_grids_lvl, p_compute_q, &
                          p_compute_lamda2, p_compute_p_extended, p_compute_max_finecell, &
                          p_mgbasb, p_compute_event_field, &
                          u, v, w, p, t, au, av, aw, ap, at, &
                          auu, avv, aww, auv, auw, avw, &
                          x, y, z, dx, dy, dz, ddx, ddy, ddz, &
                          hilf, myid, numprocs, istep, dim3d, dim2d, dim1d, &
                          nscal, lvlmin, lvlmax ) BIND(C)

            USE ISO_C_BINDING, ONLY: c_funptr, c_int
            USE precision_mod, ONLY: c_realk

            ! pointers must be passed with VALUE (it is already a pointer)
            TYPE(c_funptr), INTENT(in), VALUE :: p_mgdims, p_get_ip1, p_get_ip3, p_list_grids_lvl, &
                                                 p_compute_q, p_compute_lamda2, p_compute_p_extended, &
                                                 p_compute_max_finecell, p_mgbasb, p_compute_event_field

            ! rest is handled as pointers
            ! WARNING: Implementation in C++ currently only supports fields in double precision (TODO)
            REAL(kind=c_realk), INTENT(in), DIMENSION(*) :: u, v, w, p, t, au, av, aw, ap, at, hilf
            REAL(kind=c_realk), INTENT(in), DIMENSION(*) :: x, y, z, dx, dy, dz, ddx, ddy, ddz
            REAL(kind=c_realk), INTENT(in), DIMENSION(*) :: auu, avv, aww, auv, auw, avw
            INTEGER(kind=c_int), INTENT(in) :: myid, numprocs, istep, dim3d, dim2d, dim1d, &
                                               nscal, lvlmin, lvlmax

        END SUBROUTINE proto
    END INTERFACE

    ! function pointer Fortran
    PROCEDURE(proto), POINTER :: render_vtk_shlib

    ! private variables
    CHARACTER(len=*), PARAMETER :: vtkdir = "VTK"
    LOGICAL :: isInit = .false.
    LOGICAL :: exists = .false.

    PUBLIC :: render_checkpoint

CONTAINS

    SUBROUTINE render_checkpoint( istep, lvlmin, lvlmax )

        USE ISO_C_BINDING, ONLY: C_FUNPTR, C_NULL_FUNPTR, &
        C_F_PROCPOINTER, C_ASSOCIATED, C_FUNLOC, C_INT

        USE shlib_mod, ONLY: shlib_get_fun
        USE pointer_mod, ONLY: get_ip1, get_ip3
        USE create_directory_mod, ONLY: create_directory

        ! interface variables
        INTEGER, INTENT(in) :: istep, lvlmin, lvlmax

        ! local variables
        TYPE(C_FUNPTR) :: funptr

        TYPE(C_FUNPTR) :: cp_get_ip1, cp_get_ip3, cp_mgdims, cp_iterate_grids_lvl, &
                          cp_compute_q, cp_compute_lambda2, cp_compute_p_extended, &
                          cp_compute_max_finecell, cp_mgbasb, cp_compute_event_field

        INTEGER(kind=C_INT) :: c_myid, c_numprocs, c_istep, c_lvlmin, c_lvlmax, &
                               c_idim3d, c_idim2d, c_idim1d, c_nscal

        ! setting the function pointers
        ! (necessary to transfer MGLET functions to the library)

        ! get C procedure pointer of the fuctions witk INT(kind=4) arguments
        cp_mgdims = C_FUNLOC( c_mgdims )
        cp_get_ip1 = C_FUNLOC( c_get_ip1 )
        cp_get_ip3 = C_FUNLOC( c_get_ip3 )
        cp_iterate_grids_lvl = C_FUNLOC( c_iterate_grids_lvl )
        cp_mgbasb = C_FUNLOC( c_mgbasb )
        cp_compute_q = C_FUNLOC( c_compute_q )
        cp_compute_lambda2 = C_FUNLOC( c_compute_lambda2 )
        cp_compute_p_extended = C_FUNLOC( c_compute_p_extended )
        cp_compute_max_finecell = C_FUNLOC( c_compute_max_finecell )
        cp_compute_event_field = C_FUNLOC( c_compute_event_field )

        ! convert the remaining integer to ensure consistency
        c_myid = INT( myid, kind=C_INT )
        c_numprocs = INT( numprocs, kind=C_INT )
        c_istep = INT( istep, kind=C_INT )
        c_idim3d = INT( idim3d, kind=C_INT )
        c_idim2d = INT( idim2d, kind=C_INT )
        c_idim1d = INT( idim1d, kind=C_INT )
        c_nscal = INT( nscal, kind=C_INT )
        c_lvlmin = INT( lvlmin, kind=C_INT )
        c_lvlmax = INT( lvlmax, kind=C_INT )

        ! -- Mechanism to include a dynamically linked library --

        ! getting the pointer to the C/C-binding function
        funptr = C_NULL_FUNPTR
        CALL shlib_get_fun("render_vtk_function_c", funptr, required=.FALSE.)

        ! check if the function was found and the pointer is set
        IF ( C_ASSOCIATED(funptr) ) THEN

            ! creating the folder if necessary
            IF ( myid == 0 .and. (.not. isInit) ) THEN
                CALL create_directory( vtkdir )
                isInit = .true.
            END IF

            ! convert C pointer to Fortran function pointer
            CALL C_F_PROCPOINTER( funptr, render_vtk_shlib )

            ! call the function via the pointer
            CALL render_vtk_shlib( cp_mgdims, cp_get_ip1, &
            cp_get_ip3, cp_iterate_grids_lvl, cp_compute_q, &
            cp_compute_lambda2, cp_compute_p_extended, &
            cp_compute_max_finecell, cp_mgbasb, cp_compute_event_field, &
            u, v, w, p, t, au, av, aw, ap, at, &
            auum, avvm, awwm, auvm, auwm, avwm, &
            x, y, z, dx, dy, dz, ddx, ddy, ddz, &
            hilf, c_myid, c_numprocs, c_istep, &
            c_idim3d, c_idim2d, c_idim1d, &
            c_nscal, c_lvlmin, c_lvlmax )

        END IF

    END SUBROUTINE


    ! Return grid dimensions with INT(kind=4) arguments
    SUBROUTINE c_mgdims(kk, jj, ii, igrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(OUT) :: kk, jj, ii
        INTEGER(kind=c_int), INTENT(in) :: igrid
        INTEGER :: kk_tmp, jj_tmp, ii_tmp, igrid_tmp
        ! Function body
        igrid_tmp = INT( igrid )
        CALL mgdims( kk_tmp, jj_tmp, ii_tmp, igrid_tmp )
        kk = INT( kk_tmp, kind=c_int )
        jj = INT( jj_tmp, kind=c_int )
        ii = INT( ii_tmp, kind=c_int )
    END SUBROUTINE c_mgdims


    ! Get 3D pointer with INT(kind=4) arguments
    SUBROUTINE c_get_ip3(ip3, igrid) bind(C)
        USE pointer_mod, only: get_ip3
        USE precision_mod, only: intk
        USE ISO_C_BINDING, ONLY: c_int
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: igrid
        INTEGER(kind=c_int), INTENT(out) :: ip3
        INTEGER(kind=intk) :: igrid_tmp, ip3_tmp
        ! Function body
        igrid_tmp = INT( igrid, kind=intk )
        CALL get_ip3( ip3_tmp, igrid_tmp )
        ip3 = INT( ip3_tmp, kind=c_int )
    END SUBROUTINE c_get_ip3


    ! Get 1D pointer with INT(kind=4) arguments
    SUBROUTINE c_get_ip1(ip1, igrid) bind(C)
        USE pointer_mod, only: get_ip1
        USE precision_mod, only: intk
        USE ISO_C_BINDING, ONLY: c_int
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: igrid
        INTEGER(kind=c_int), INTENT(out) :: ip1
        INTEGER(kind=intk) :: igrid_tmp, ip1_tmp
        ! Function body
        igrid_tmp = INT( igrid, kind=intk )
        CALL get_ip1( ip1_tmp, igrid_tmp )
        ip1 = INT( ip1_tmp, kind=c_int )
    END SUBROUTINE c_get_ip1


    ! Returns the boundary type at front, back, right, left, bottom, top for igrid
    SUBROUTINE c_mgbasb(cnfro, cnbac, cnrgt, cnlft, cnbot, cntop, cigrid) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: cigrid
        INTEGER(kind=c_int), INTENT(out) :: cnfro, cnbac, cnrgt, cnlft, cnbot, cntop
        INTEGER :: nfro, nbac, nrgt, nlft, nbot, ntop, igrid
        igrid = INT(cigrid)
        ! calling the MGLET routine
        CALL mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
        ! transforming back for C interoperability
        cnfro = INT(nfro, kind=c_int)
        cnbac = INT(nbac, kind=c_int)
        cnrgt = INT(nrgt, kind=c_int)
        cnlft = INT(nlft, kind=c_int)
        cnbot = INT(nbot, kind=c_int)
        cntop = INT(ntop, kind=c_int)
    END SUBROUTINE c_mgbasb


    ! Facilitates iteration over local grids on certain level with INT(kind=4) arguments
    SUBROUTINE c_iterate_grids_lvl(ilevel, num, igrid) bind(C)
        USE setmpi_mod, ONLY: nmygridslvl, mygridslvl
        USE ISO_C_BINDING, ONLY: c_int
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: ilevel
        INTEGER(kind=c_int), INTENT(in) :: num
        INTEGER(kind=c_int), INTENT(out) :: igrid
        INTEGER :: igrid_tmp, max
        ! Function body
        max = nmygridslvl(ilevel)
        IF ( num > max ) THEN
            igrid = INT( -1, kind=c_int ) ! indicates end
        ELSE
            igrid_tmp = mygridslvl(num, ilevel)
            igrid = INT( igrid_tmp, kind=c_int )
        END IF
    END SUBROUTINE c_iterate_grids_lvl


    ! Facilitates iteration over local grids on certain level with INT(kind=4) arguments
    SUBROUTINE c_compute_q(ilevel) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE setmpi_mod, ONLY: nmygridslvl, mygridslvl
        USE pointer_mod, ONLY: get_ip1, get_ip3
        USE fields_mod, ONLY: gsaw, gsae, gsan, gsas, gsat, gsab, gsap => geovp, bp
        USE connect2_mod, ONLY: Connect2 => Connect
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: ilevel
        INTEGER :: igrid, i, ii, jj, kk, ip1, ip3
        ! Function body
        CALL Connect2( ilevel, 2, s1=p, corners=.TRUE. ); hilf = 0.0;
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL mgdims(kk, jj, ii, igrid)
            CALL get_ip1(ip1, igrid)
            CALL get_ip3(ip3, igrid)
            CALL laplacephi(kk, jj, ii, p(ip3), hilf(ip3), &
            gsaw(ip1), gsae(ip1), gsan(ip1), gsas(ip1), gsat(ip1), gsab(ip1), gsap(ip3) &
#ifdef _IB_
            , bp(ip3) &
#endif
            )
            CALL maskbp(kk, jj, ii, hilf(ip3), bp(ip3))
        END DO
        CALL Connect2( ilevel, 2, s1=hilf, corners=.TRUE. )
    END SUBROUTINE c_compute_q


    ! Masking the output with a possibly extended field BP
    SUBROUTINE maskbp(kk, jj, ii, field, bp)
        USE precision_mod, only: realk
        IMPLICIT NONE
        INTEGER, INTENT(in) :: kk, jj, ii
        REAL(kind=realk), INTENT(inout) :: field(kk,jj,ii)
        REAL(kind=realk), INTENT(in) :: bp(kk,jj,ii)
        INTEGER :: i, j, k, oi, oj, ok
        INTEGER, PARAMETER :: d = 1
        DO i = 1+d, ii-d, 1
            DO j = 1+d, jj-d, 1
                DO k = 1+d, kk-d, 1
                ! checking the surrounding
                  DO oi = -d, d
                    DO oj = -d, d
                      DO ok = -d, d
                        field(k,j,i) = field(k,j,i) * &
                        bp(k+ok,j+oj,i+oi)
                      END DO
                    END DO
                  END DO
                ! checking the surrounding
                END DO
            END DO
        END DO
    END SUBROUTINE maskbp


    ! Computation of the lambda2 field
    SUBROUTINE c_compute_lambda2(ilevel) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE setmpi_mod, ONLY: nmygridslvl, mygridslvl
        USE pointer_mod, ONLY: get_ip1, get_ip3
        USE fields_mod, ONLY: u, v, w, dx, dy, dz, ddx, ddy, ddz, hilf
        USE connect2_mod, ONLY: connect2 => connect
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: ilevel
        INTEGER :: igrid, i, ii, jj, kk, ip1, ip3
        ! Function body
        CALL connect2( ilevel, 2, v1=u, v2=v, v3=w, corners=.TRUE. );
        hilf = 0.0;
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL mgdims(kk, jj, ii, igrid)
            CALL get_ip1(ip1, igrid)
            CALL get_ip3(ip3, igrid)
            CALL lambda2_grid(kk, jj, ii, u(ip3), v(ip3), w(ip3), &
            dx(ip1), dy(ip1), dz(ip1), ddx(ip1), ddy(ip1), ddz(ip1), hilf(ip3))
            CALL maskbp(kk, jj, ii, hilf(ip3), bp(ip3))
        END DO
        CALL connect2( ilevel, 2, s1=hilf, corners=.TRUE. )
    END SUBROUTINE c_compute_lambda2


    SUBROUTINE lambda2_grid( kk, jj, ii, u, v, w, dx, dy, dz, ddx, ddy, ddz, lambda2 )
        USE precision_mod, only: realk
        IMPLICIT NONE
        INTEGER, INTENT(in) :: kk, jj, ii
        REAL(kind=realk), INTENT(in) :: u(kk,jj,ii), v(kk,jj,ii), w(kk,jj,ii)
        REAL(kind=realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(kind=realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(kind=realk), INTENT(out) :: lambda2(kk,jj,ii)
        ! local variables
        INTEGER :: i, j, k, it, jt
        INTEGER, PARAMETER :: d = 1
        REAL(kind=realk) :: up, uc, um, vp, vc, vm, wp, wc, wm ! p = plus, c = center, m = minus (in direction of derivative)
        REAL(kind=realk) :: e1, e2, e3, eig1, eig2, eig3, p1, p2, p3, q, p, r, phi
        REAL(kind=realk), PARAMETER :: PI = ACOS(-1.0)
        ! tensors in C indexing (this allows to translate routines from MGTOOLS)
        REAL(kind=realk) :: gradient(0:8)
        REAL(kind=realk) :: tensorS(0:2,0:2), tensorW(0:2,0:2), A(0:2,0:2), B(0:2,0:2)

        lambda2 = 0.0;
        DO i = 1+d, ii-d, 1
            DO j = 1+d, jj-d, 1
                DO k = 1+d, kk-d, 1
                    ! compute the velocity gradient tensor
                    ! [dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz]
                    gradient(0) = ( u(k,j,i) - u(k,j,i-1) ) / ddx(i)
                    um = 0.5 * (u(k,j-1,i) + u(k,j-1,i-1))
                    uc = 0.5 * (u(k,j,  i) + u(k,j,  i-1))
                    up = 0.5 * (u(k,j+1,i) + u(k,j+1,i-1))
                    gradient(1) = 0.5 * (uc-um)/dy(j-1) + 0.5 * (up-uc)/dy(j)
                    um = 0.5 * (u(k-1,j,i) + u(k-1,j,i-1))
                    uc = 0.5 * (u(k  ,j,i) + u(k,  j,i-1))
                    up = 0.5 * (u(k+1,j,i) + u(k+1,j,i-1))
                    gradient(2) = 0.5 * (uc-um)/dz(k-1) + 0.5 * (up-uc)/dz(k)

                    vm = 0.5 * (v(k,j,i-1) + v(k,j-1,i-1))
                    vc = 0.5 * (v(k,j,i  ) + v(k,j-1,i  ))
                    vp = 0.5 * (v(k,j,i+1) + v(k,j-1,i+1))
                    gradient(3) = 0.5 * (vc-vm)/dx(i-1) + 0.5 * (vp-vc)/dx(i)
                    gradient(4) = ( v(k,j,i) - v(k,j-1,i) ) / ddy(j)
                    vm = 0.5 * (v(k-1,j,i) + v(k-1,j-1,i))
                    vc = 0.5 * (v(k  ,j,i) + v(k  ,j-1,i))
                    vp = 0.5 * (v(k+1,j,i) + v(k+1,j-1,i))
                    gradient(5) = 0.5 * (vc-vm)/dz(k-1) + 0.5 * (vp-vc)/dz(k)

                    wm = 0.5 * (w(k,j,i-1) + w(k-1,j,i-1))
                    wc = 0.5 * (w(k,j,i  ) + w(k-1,j,i  ))
                    wp = 0.5 * (w(k,j,i+1) + w(k-1,j,i+1))
                    gradient(6) = 0.5 * (wc-wm)/dx(i-1) + 0.5 * (wp-wc)/dx(i)
                    wm = 0.5 * (w(k,j-1,i) + w(k-1,j-1,i))
                    wc = 0.5 * (w(k,j,  i) + w(k-1,j,  i))
                    wp = 0.5 * (w(k,j+1,i) + w(k-1,j+1,i))
                    gradient(7) = 0.5 * (wc-wm)/dy(j-1) + 0.5 * (wp-wc)/dy(j)
                    gradient(8) = ( w(k,j,i) - w(k-1,j,i) ) / ddz(k)

                    ! Computing strain and rotation rate tensor
                    DO jt = 0, 2
                        DO it = 0, 2
                            tensorS(it,jt) = 0.5 * ( gradient(it*3+jt) + gradient(jt*3+it) )
                            tensorW(it,jt) = 0.5 * ( gradient(it*3+jt) - gradient(jt*3+it) )
                        END DO
                    END DO
                    ! Computing matrix A = S2 + Om2
                    DO jt = 0, 2
                        DO it = 0, 2
                            A(it,jt) = &
                            tensorS(0,jt)*tensorS(it,0) + tensorS(1,jt)*tensorS(it,1) + tensorS(2,jt)*tensorS(it,2) + &
                            tensorW(0,jt)*tensorW(it,0) + tensorW(1,jt)*tensorW(it,1) + tensorW(2,jt)*tensorW(it,2)
                        END DO
                    END DO
                    p1 = A(0,1)*A(0,1) + A(0,2)*A(0,2) + A(1,2)*A(1,2);
                    IF ( p1 == 0.0 ) THEN
                        e1 = A(0,0); e2 = A(1,1); e3 = A(2,2)
                        eig2 = e1 + e2 + e3 - MAX(MAX(e1, e2), e3) - MIN(MIN(e1, e2), e3)
                    ELSE
                        q = ( A(0,0) + A(1,1) + A(2,2) ) / 3.0
                        p2 = (A(0,0) - q) * (A(0,0) - q) + &
                             (A(1,1) - q) * (A(1,1) - q) + &
                             (A(2,2) - q) * (A(2,2) - q) + 2.0 * p1;
                        p = SQRT( p2 / 6.0 );
                        ! Filling matrix B
                        DO jt = 0, 2
                            DO it = 0, 2
                                IF ( it == jt ) THEN
                                    B(it,jt) = ( 1.0 / p ) * ( A(it,jt) - q )
                                ELSE
                                    B(it,jt) = ( 1.0 / p ) * A(it,jt)
                                END IF
                            END DO
                        END DO
                        r =  ( B(0,0) * ( B(1,1)*B(2,2) - B(1,2)*B(2,1) ) &
                             - B(0,1) * ( B(1,0)*B(2,2) - B(1,2)*B(2,0) ) &
                             + B(0,2) * ( B(1,0)*B(2,1) - B(2,2)*B(2,0) ) ) / 2.0
                        ! Correction of possible small precision errors
                        IF ( r <= -1.0 ) THEN
                            phi = pi / 3.0
                        ELSE IF ( r >= 1.0 ) THEN
                            phi = 0.0
                        ELSE
                            phi = ACOS(r) / 3.0
                        END IF
                        ! Eigenvalues eig3 <= eig2 <= eig1
                        eig1 = q + 2.0 * p * COS(phi)
                        eig3 = q + 2.0 * p * COS(phi + pi*(2.0/3.0))
                        eig2 = 3.0 * q - eig1 - eig3
                    END IF
                    lambda2(k,j,i) = eig2;
                END DO
            END DO
        END DO
    END SUBROUTINE lambda2_grid


    ! Computation of the pressure field extended into the first row with BP=0
    SUBROUTINE c_compute_p_extended(ilevel) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE setmpi_mod, ONLY: nmygridslvl, mygridslvl
        USE pointer_mod, ONLY: get_ip1, get_ip3
        USE fields_mod, ONLY: hilf, p, bp
        USE connect2_mod, ONLY: Connect2 => Connect
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: ilevel
        INTEGER :: igrid, i, ii, jj, kk, ip3
        ! Function body
        hilf = p
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            ! boundp( igrid, kk, jj, ii, p )
            CALL boundp(igrid, kk, jj, ii, hilf(ip3))
        END DO
        CALL Connect2( ilevel, 2, s1=hilf, corners=.TRUE. )
    END SUBROUTINE c_compute_p_extended


    ! Computation of the maximum value of finecell inside the grid
    SUBROUTINE c_compute_max_finecell(igrid, imax) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE pointer_mod, ONLY: get_ip3
        USE fields_mod, ONLY: finecell
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: igrid
        INTEGER(kind=c_int), INTENT(out) :: imax
        INTEGER :: ii, jj, kk, ip3, result
        ! Function body
        CALL mgdims(kk, jj, ii, INT(igrid))
        CALL get_ip3(ip3, INT(igrid))
        CALL max_finecell_grid(finecell(ip3), kk, jj, ii, result)
        imax = INT( result, kind=c_int )
    END SUBROUTINE c_compute_max_finecell

    SUBROUTINE max_finecell_grid(finecell, kk, jj, ii, result)
        IMPLICIT NONE
        ! Variables
        INTEGER, INTENT(in) :: kk, jj, ii
        REAL(kind=realk), INTENT(in) :: finecell(kk,jj,ii)
        INTEGER, INTENT(out) :: result
        REAL(kind=realk) :: val
        ! Function body
        val = MAXVAL( finecell( 3:(kk-2), 3:(jj-2), 3:(ii-2) ) )
        result = 0
        IF ( val > 0.5 ) THEN
            result = 1
        END IF
    END SUBROUTINE max_finecell_grid


    ! Computation of the underlaying field for sweep / eject event detection
    SUBROUTINE c_compute_event_field(ilevel) bind(C)
        USE ISO_C_BINDING, ONLY: c_int
        USE setmpi_mod, ONLY: nmygridslvl, mygridslvl
        USE pointer_mod, ONLY: get_ip1, get_ip3
        USE fields_mod, ONLY: u, v, w, hilf
        USE statfields_mod, ONLY: au, av, aw
        USE connect2_mod, ONLY: connect2 => connect
        IMPLICIT NONE
        ! Variables
        INTEGER(kind=c_int), INTENT(in) :: ilevel
        INTEGER :: igrid, i, ii, jj, kk, ip3
        ! Function body
        CALL connect2( ilevel, 2, v1=u, v2=v, v3=w, corners=.TRUE. );
        CALL connect2( ilevel, 2, v1=au, v2=av, v3=aw, corners=.TRUE. );
        hilf = 0.0;
        DO i = 1, nmygridslvl(ilevel)
            igrid = mygridslvl(i, ilevel)
            CALL mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)
            CALL event_grid( kk, jj, ii, u(ip3), w(ip3), &
                             au(ip3), aw(ip3), hilf(ip3) )
            CALL maskbp(kk, jj, ii, hilf(ip3), bp(ip3))
        END DO
        CALL connect2( ilevel, 2, s1=hilf, corners=.TRUE. )
    END SUBROUTINE c_compute_event_field

    SUBROUTINE event_grid( kk, jj, ii, u, w, au, aw, hilf )
        USE precision_mod, only: realk, great
        IMPLICIT NONE
        INTEGER, INTENT(in) :: kk, jj, ii
        REAL(kind=realk), INTENT(in) :: u(kk,jj,ii), w(kk,jj,ii)
        REAL(kind=realk), INTENT(in) :: au(kk,jj,ii), aw(kk,jj,ii)
        REAL(kind=realk), INTENT(out) :: hilf(kk,jj,ii)
        ! local variables
        INTEGER :: i, j, k, it, jt
        INTEGER, PARAMETER :: d = 1
        REAL :: inst_uwm, abs_auwm, val1, val2, uc, wc, fac
        ! all computations for position stag = [1 0 1]
        hilf = 0.0
        DO i = 1+d, ii-d, 1
            DO j = 1+d, jj-d, 1
                DO k = 1+d, kk-d, 1
                    ! computing instantaneous fluctuations
                    uc = 0.5 * ( (u(k,j,i-1)-au(k,j,i-1)) + (u(k,j,i)-au(k,j,i)) )
                    wc = 0.5 * ( (w(k-1,j,i)-aw(k-1,j,i)) + (w(k,j,i)-aw(k,j,i)) )
                    inst_uwm = uc * wc
                    ! only negative values indicate sweep / eject events
                    ! implementation via fac = 0 if inst_uwm > 0, inst_uwm = 1 if uc < 0
                    fac = 0.5 * ( 1.0 + SIGN(1.0,-inst_uwm) )
                    ! eject = positive
                    ! sweep = negative
                    ! else = 0
                    hilf(k,j,i) = fac * inst_uwm * SIGN(1.0,uc)
                END DO
            END DO
        END DO
    END SUBROUTINE event_grid

END MODULE render_mod
