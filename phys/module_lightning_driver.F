!WRF:MEDIATION_LAYER:PHYSICS
!
! Contains initialization subroutine lightning_init and driver subroutine
! lightning_driver.
!
! History:
!   3.4.1 - rewritten and added init, separate out flash rate
!           parameterization from emission
!   3.3.x - lightning_driver written by M. Barth called in
!           emission_driver in chem
!
! Contact: J. Wong <johnwong@ucar.edu>
!
!**********************************************************************

 MODULE module_lightning_driver
 CONTAINS


!**********************************************************************
!
! SUBROUTINE lightning_init
!
! Performs compatibility checks and zero out flash arrays at first timestep.
!
!**********************************************************************

 SUBROUTINE lightning_init ( &
                              itimestep, restart, dt,                   &
                            ! Namelist control options
                              cu_physics,                               &
                              lightning_option, lightning_dt,           &
                              lightning_start_seconds,                  &
                              iccg_prescribed_num, iccg_prescribed_den, &
                            ! Order dependent args for domain, mem, and tile dims
                              ids, ide, jds, jde, kds, kde,             &
                              ims, ime, jms, jme, kms, kme,             &
                              its, ite, jts, jte, kts, kte,             &
                            ! IC and CG flash rates and accumulated flash count
                              ic_flashcount, ic_flashrate,              &
                              cg_flashcount, cg_flashrate               &
                            )
!-----------------------------------------------------------------
 USE module_state_description
 USE module_wrf_error
 IMPLICIT NONE
!-----------------------------------------------------------------

 INTEGER,  INTENT(IN)        :: itimestep
 LOGICAL,  INTENT(IN)        :: restart
 REAL,     INTENT(IN)        :: dt
 INTEGER,  INTENT(IN)        :: cu_physics,lightning_option
 REAL,     INTENT(IN)        :: lightning_dt, lightning_start_seconds
 REAL,     INTENT(IN)        :: iccg_prescribed_num, iccg_prescribed_num
 INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte

! Making these optional just in case qualitative lightning indices get implemented
 REAL, OPTIONAL, DIMENSION( ims:ime,kms:kme,jms:jme ), &
                 INTENT(OUT) :: ic_flashcount, ic_flashrate, &
                                cg_flashcount, cg_flashrate

 CHARACTER (LEN=80) :: message

!-----------------------------------------------------------------

!-- do not reset unless it is the first timestep or lightning_option is on
 IF (itimestep .gt. 0 .or. lightning_option .eq. 0) return

!-- check to see if lightning_dt is a proper multiple of dt
 IF ( MOD(lightning_dt,dt) .ne. 0. ) THEN
    CALL wrf_error_fatal (' lightning_init: lightning_dt needs to be a multiple of model time step dt')
 ENDIF

!--  check to see if the prescribed IC:CG ratio is valid (0/0 and -1 are not allowed)
 IF (iccg_prescribed_den .eq. 0. .and. iccg_prescribed_num .eq. 0.) THEN
    CALL wrf_error_fatal (' lightning_init: iccg_prescribed cannot be 0.0/0.0')
 ENDIF
 IF (iccg_prescribed_den .ne. 0.) THEN
    IF (iccg_prescribed_num/iccg_prescribed_den .eq. -1.) THEN
        CALL wrf_error_fatal (' lightning_init: iccg_prescribed cannot be -1')
    ENDIF
 ENDIF

!
!-- check to see if lightning_option is valid
!
!   Add new schemes here so it is recognized and proper checks are performed
!
 ltng_select: SELECT CASE(lightning_option)

    ! Convective resolved/permitted
    CASE (ltng_crm_PR92w)
        WRITE(message, * ) ' lightning_init: CRM lightning option used: ', lightning_option
        CALL wrf_debug ( 200 , message )

    ! Convective parameterized
    CASE (ltng_cpm_PR92z)
        IF ( cu_physics .ne. GDSCHEME .and. cu_physics .ne. G3SCHEME  ) THEN
          CALL wrf_error_fatal( ' lightning_init: Selected lightning option requires GD or G3 convective parameterization' )
        ENDIF

        WRITE(message, * ) ' lightning_init: CPM lightning option selected: ', lightning_option
        CALL wrf_debug ( 200 , message )

    ! Non-existing options
    CASE DEFAULT
        CALL wrf_error_fatal ( ' lightning_init: invalid lightning_option')
 END SELECT ltng_select

!-- do not re-initialize for restarts
 IF (restart) return

!-- zero out arrays
 IF ( PRESENT( ic_flashcount ) .and. PRESENT( ic_flashrate ) .and. &
      PRESENT( cg_flashcount ) .and. PRESENT( cg_flashrate ) ) THEN
    CALL wrf_debug ( 200 , ' lightning_init: flash initializing lightning flash arrays' )

    ic_flashrate(:,:,:)  = 0.
    ic_flashcount(:,:,:) = 0.
    cg_flashrate(:,:,:)  = 0.
    cg_flashcount(:,:,:) = 0.
 ELSE
    CALL wrf_error_fatal ( ' lightning_init: flash arrays not present' )
 ENDIF

 CALL wrf_debug( 200, ' lightning_init: finishing')

 END SUBROUTINE lightning_init


!**********************************************************************
!
! SUBROUTINE lightning_driver
!
! Redirect to the appropriate lightning subroutine.
!
!**********************************************************************

 SUBROUTINE lightning_driver ( &
                          ! Frequently used prognostics
                            itimestep, dt, dx, dy,                &
                            xlat, xlon, xland, ht,                &
                            t_phy, p_phy, rho, u, v, w,           &
                            z, moist,                             &
                          ! Scheme specific prognostics
                            ktop_deep,                            &
                          ! Mandatory namelist inputs
                            lightning_option,                     &
                            lightning_dt,                         &
                            lightning_start_seconds,              &
                            flashrate_factor,                     &
                          ! IC:CG namelist settings
                            iccg_method,                          &
                            iccg_prescribed_num,                  &
                            iccg_prescribed_den,                  &
                          ! Scheme specific namelist inputs
                            cellcount_method,                     &
                            cldtop_adjustment,                    &
                            ltng_temp_upper,                      &
                            ltng_temp_lower,                      &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte,         &
                          ! Mandatory outputs for all quantitative schemes
                            ic_flashcount, ic_flashrate,          &
                            cg_flashcount, cg_flashrate           &
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

! Parameterization options
! USE module_ltng_pickering    ! lightning_option == 1,   ltng_crm_PR92w
 USE module_ltng_cpmpr92z      ! lightning_option == 11,  ltng_cpm_PR92z

! IC:CG methods
 USE module_ltng_iccg

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 INTEGER, INTENT(IN   )    ::       itimestep
 REAL,    INTENT(IN   )    ::       dt, dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: xlat, xlon, xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: t_phy, t_phy, rho
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: u, v, w, z
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme, num_moist), INTENT(IN   ) :: moist

! Scheme specific prognostics
 INTEGER, DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: ktop_deep     ! model LNB from cu_physics

! Mandatory namelist inputs
 INTEGER, INTENT(IN   )    ::       lightning_option
 REAL,    INTENT(IN   )    ::       lightning_dt, lightning_start_seconds, flashrate_factor

! IC:CG namelist settings
 INTEGER, INTENT(IN   )    ::       iccg_method
 REAL,    INTENT(IN   )    ::       iccg_prescribed_num, iccg_prescribed_den

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method                    ! used in CRM
 REAL,    INTENT(IN   )    ::       cldtop_adjustment                   ! used in CPM
 REAL,    INTENT(IN   )    ::       ltng_temp_upper,ltng_temp_lower     ! used in CRM

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       its,ite, jts,jte, kts,kte

! Mandatory outputs for all quantitative schemes
 REAL, OPTIONAL, DIMENSION( ims:ime,          jms:jme ), INTENT(INOUT) :: ic_flashcount , cg_flashcount
 REAL, OPTIONAL, DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: ic_flashrate  , cg_flashrate

! Local variables
 REAL, DIMENSION( ims:ime, jms:jme ) :: total_flashrate
 CHARACTER (LEN=80) :: message

!-----------------------------------------------------------------

 IF ( lightning_option .eq. 0 ) RETURN

 IF ( itimestep * dt .lt. lightning_start_seconds ) RETURN

 IF ( MOD((itimestep * dt - lightning_start_seconds), lightning_dt ) .ne. 0 ) RETURN

!-----------------------------------------------------------------

 CALL wrf_debug ( 100, ' lightning_driver: calculating flash rate' )
 flashrate_select: SELECT CASE(lightning_option)

    ! CRM lightning options
    CASE( ltng_crm_PR92w )
        CALL wrf_debug ( 100, ' lightning_driver: calling Price and Rind 1992 (w_max, CRM)' )

!   CASE ( another_crm_option)
!       CALL ...

    ! CPM lightning options
    CASE( ltng_cpm_PR92z )
        CALL wrf_debug ( 100, ' lightning_driver: calling Price and Rind 1992 (z_top, CPM)' )

        CALL ltng_cpmpr92z ( &
                  ! Frequently used prognostics
                    dx, dy, xland, ht, z, t_phy,      &
                    ktop_deep, cldtop_adjustment,     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,     &
                    ims, ime, jms, jme, kms, kme,     &
                    its, ite, jts, jte, kts, kte,     &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate
                  )

!   CASE ( another_cpm_option)
!       CALL ...

    ! Invalid lightning options
    CASE DEFAULT
        WRITE(wrf_err_message, * ) ' lightning_driver: The lightning option does not exist: lightning = ', lightning_option
        CALL wrf_error_fatal ( wrf_err_message )

 END SELECT flashrate_select

!-----------------------------------------------------------------

 CALL wrf_debug ( 100, ' lightning_driver: partitioning IC:CG')
 iccg_select: SELECT CASE(iccg_method)
    ! Flash rate option defaults
    CASE( 0 )
        CALL wrf_debug( 100, ' lightning_driver: using option-default IC:CG method'
        iccg_method_default: SELECT CASE(lightning_option)

            CASE( ltng_crm_PR92w,ltng_cpm_PR92z )
                CALL iccg_boccippio( &
                            xlat, xlon                                 &
                            iccg_prescribed_num, iccg_prescribed_den,  &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,              &
                            ims, ime, jms, jme, kms, kme,              &
                            its, ite, jts, jte, kts, kte,              &
                          ! Input
                            total_flashrate,                           &
                          ! Output
                            ic_flashrate, cg_flashrate                 &
                          )

            CASE DEFAULT
                CALL wrf_debug ( 100, ' lightning_driver: no method-default IC:CG implemented, using user-prescribed constant')
                CALL iccg_prescribed( &
                            iccg_prescribed_num,                  &
                            iccg_prescribed_den,                  &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte,         &
                          ! Input
                            total_flashrate,                      &
                          ! Output
                            ic_flashrate, cg_flashrate            &
                          )

        END SELECT iccg_method_default

    ! Used-prescribed constant
    CASE( 1 )
        WRITE(message, * ) ' lightning_driver: using user-prescribed IC:CG ratio = ', iccg_prescribed_num, iccg_prescribed_den
        CALL wrf_debug ( 100, message )
        CALL iccg_prescribed( &
                    iccg_prescribed_num,                  &
                    iccg_prescribed_den,                  &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Input
                    total_flashrate,                      &
                  ! Output
                    ic_flashrate, cg_flashrate            &
                  )

    ! Boccippio et al, 2001
    CASE( 2 )
        CALL wrf_debug ( 100, ' lightning_driver: using Boccippio 2001 IC:CG climatology')
        CALL iccg_boccippio( &
                    xlat, xlon                                 &
                    iccg_prescribed_num, iccg_prescribed_den,  &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                  )

    ! Price and Rind, 1993
    CASE( 3 )
        CALL wrf_debug ( 100, ' lightning_driver: using Price and Rind 1993 IC:CG ratio')
        CALL iccg_pr93( &
                    ktop_deep, cldtop_adjustment, t_phy, z,    &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                )

    ! Invalid IC:CG method
    CASE DEFAULT
        WRITE(wrf_err_message, * ) ' lightning_driver: Invalid IC:CG method (iccg_method) = ', lightning_option
        CALL wrf_error_fatal ( wrf_err_message )

 END SELECT iccg_select

!-----------------------------------------------------------------

 CALL wrf_debug( 200, ' lightning_driver: converting flash rate to total')

 ic_flashcount(its:ite,jts:jte) = ic_flashcount(its:ite,jts:jte) + ic_flashrate(its:ite,jts:jte) * flashrate_factor * lightning_dt
 cg_flashcount(its:ite,jts:jte) = cg_flashcount(its:ite,jts:jte) + cg_flashrate(its:ite,jts:jte) * flashrate_factor * lightning_dt

!-----------------------------------------------------------------

 CALL wrf_debug ( 100, ' lightning_driver: returning from')

 END SUBROUTINE lightning_driver