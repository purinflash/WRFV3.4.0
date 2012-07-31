! WRF:MODEL_LAYER:PHYSICS
!
! Lightning flash rate prediction based on cloud-top height. Implemented
! for models using convective parameterization. Assume use of sub-grid LNB.
!
! Price, C., and D. Rind (1992), A Simple Lightning Parameterization for Calculating
!   Global Lightning Distributions, J. Geophys. Res., 97(D9), 9919–9933, doi:10.1029/92JD00719.
!
! Wong, J., M. Barth, and D. Noone (2012), Evaluating a Lightning Parameterization
!   at Resolutions with Partially-Resolved Convection, GMDD, in preparation.
!
! Contact: J. Wong <johnwong@ucar.edu>
!
!**********************************************************************

 MODULE module_ltng_cpmpr92z
 CONTAINS

 SUBROUTINE ltng_cpmpr92z ( &
                          ! Frequently used prognostics
                            dx, dy, xland, ht, z, t,              &
                          ! Scheme specific prognostics
                            kLNB,                                 &
                          ! Scheme specific namelist inputs
                            cldtop_adjustment,                    &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t

! Scheme specific prognostics
 INTEGER, DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: kLNB     ! model LNB from cu_physics

! Scheme specific namelist inputs
 REAL,    INTENT(IN   )    ::       cldtop_adjustment

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,        , jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL                                :: dA              ! grid area dx*dy in km2
 REAL                                :: zkm             ! AGL z in km
 REAL, DIMENSION( ips:ipe, jps:jpe ) :: ztop            ! cloud top height AGL in km

 REAL, PARAMETER:: baseArea=1296. ! base-case area, dx = 36 km

 INTEGER :: i,k,j

 CHARACTER (LEN=80) :: message

!-----------------------------------------------------------------

 dA = dx*dy/1E6

 ztop( ips:ipe,jps:jpe ) = 0.
 total_flashrate( ips:ipe,jps:jpe ) = 0.

! Compute AGL heights in km
 jloop: DO j=jps,jpe
    iloop: DO i=ips,ipe
        IF ( t(i,ktop,j) .gt. 273.15 ) CONTINUE                    ! Cloud top must freeze
        zkm = ( z(i,ktop,j) - ht(i,j) )/1E3 + cldtop_adjustment    ! Compute AGL height in km
        IF ( zkm .lt. 0. ) CONTINUE                                ! Must be above ground
        ztop(i,j) = zkm                                            
    ENDDO iloop
 ENDDO jloop

! Convert to flash rate
 WHERE ( (xland(ips:ipe,jps:jpe) .lt. 1.5) .and. (ztop(ips:ipe,jps:jpe) .gt. 0.) )
    total_flashrate(ips:ipe,jps:jpe) = 3.44E-5 * (ztop(ips:ipe,jps:jpe)**4.9) /60.
 ELSEWHERE
    total_flashrate(ips:ipe,jps:jpe) = 6.57E-6 * (ztop(ips:ipe,jps:jpe)**4.9) /60.
 ENDWHERE

! Scale by areal ratio
 total_flashrate(ips:ipe,jps:jpe) = total_flashrate(ips:ipe,jps:jpe) * dA/baseArea

 END SUBROUTINE ltng_cpmpr92z

 END MODULE module_ltng_cpmpr92z
