!WRF:MODEL_LAYER:CHEMISTRY
!
! Contains subroutine for converting flash rates into NO emissions
! based on Decaria 2000 vertical distirbutions.
!
! Input: flashes (#/s)
! Output: tendency (ppmv/s)
!
! The output will be muliplied by timestep and used to incremeent NO
! concentration and the respective passive tracer in lightning_nox_driver.
!
! See module_lightning_nox_driver for more info.
!
! Contact: M. Barth <barthm@ucar.edu>, J. Wong <johnwong@ucar.edu>
!
!**********************************************************************
 MODULE module_lightning_nox_decaria

 IMPLICIT NONE

 CONTAINS

!**********************************************************************
!
! DeCaria et al, 2000
!
! DeCaria, A. J., K. E. Pickering, G. L. Stenchikov, and L. E. Ott (2005),
! Lightning-generated NOX and its impact on tropospheric ozone production:
! A three-dimensional modeling study of a Stratosphere-Troposphere Experiment:
! Radiation, Aerosols and Ozone (STERAO-A) thunderstorm, J. Geophys. Res.,
! 110, D14303, doi:10.1029/2004JD005556.
!
!**********************************************************************
 SUBROUTINE lightning_nox_decaria ( &
                          ! Frequently used prognostics
                            dx, dy, xland, ht, t, rho, z, p,      &
                            ic_flashrate, cg_flashrate,           & ! flashes (#/s)
                          ! Scheme specific prognostics
                            refl,                                 &
                          ! Namelist inputs
                            N_IC, N_CG,                           &
                            ltng_temp_upper,ltng_temp_lower,      &
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte,         &
                          ! outputs
                            lnox_ic_tend, lnox_cg_tend            & ! tendency (ppmv/s)
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_max_real, wrf_dm_min_real, wrf_dm_sum_real

! Lightning method
 USE module_lightning_driver, only: countCells

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: t, rho, z, p

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl

! Scheme specific namelist inputs
 REAL,    INTENT(IN   )    ::       N_IC, N_CG
 REAL,    INTENT(IN   )    ::       ltng_temp_lower, ltng_temp_upper
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(  OUT) :: lnox_ic_tend,lnox_cg_tend

! Local variables
 INTEGER :: i,k,j
 REAL :: ic_fr, cg_fr ! reconsolidated flashrates
 REAL :: term2, B
 CHARACTER (LEN=250) :: message
 REAL, DIMENSION( ims:ime, jms:jme ) :: cellcount
 REAL, DIMENSION( kps:kpe ) :: z_average, t_average, p_average, rho_average, molesofair
 REAL, DIMENSION( kps:kpe ) :: fd, fd2, dz ! fd = distribution

 REAL, PARAMETER :: ref_threshold = 20.

!-----------------------------------------------------------------

 lnox_ic_tend (ips:ipe,jps:jpe ) = 0.
 lnox_cg_tend (ips:ipe,jps:jpe ) = 0.

! Determine cloud extents. Previously also calculated in physics but
! cellcount was not persistent
 CALL countCells( &
      ! Inputs
        refl, ref_threshold, cellcount_method,     &
      ! Order dependent args for domain, mem, and tile dims
        ids, ide, jds, jde, kds, kde,              &
        ims, ime, jms, jme, kms, kme,              &
        ips, ipe, jps, jpe, kps, kpe,              &
      ! Output
        cellcount )

! Reconsolidate flash counts
 ic_fr = sum(ic_flashrate(ips:ipe,jps:jpe))
 cg_fr = sum(cg_flashrate(ips:ipe,jps:jpe))
 if ( cellcount_method .eq. 2 ) then
    ic_fr = wrf_dm_sum_real(ic_fr)
    cg_fr = wrf_dm_sum_real(cg_fr)
 ENDIF

!-----------------------------------------------------------------

! Average z, t, p, rho
 CALL horizontalAverage( z, ips, ipe, kps, kpe, jps, jpe, z_average )
 CALL horizontalAverage( t, ips, ipe, kps, kpe, jps, jpe, t_average )
! CALL horizontalAverage( p, ips, ipe, kps, kpe, jps, jpe, p_average )
 CALL horizontalAverage( rho, ips, ipe, kps, kpe, jps, jpe, rho_average )

 molesofair(kps:kpe) = rho_average(kps:kpe) * 1E3 * dx * dy / .02897 ! # moles per km in z

! term2 = 30 * 8.3145E6/dx/dy/28.96/100./100.

! Calculates IC distribution
 IF (( ic_fr .gt. 0 ) .and. (( ktop .gt. klower ) .and. (kbtm .lt. ktop) ) )THEN
   call bellcurve(kbot,ktop,klower,z_average, kps,kpe, fd, dz)
   if (ktop .gt. kupper) then
     call bellcurve(kbot,ktop,kupper,z_average, kps,kpe, fd2, dz)
     fd(kbot:ktop) = 0.5*( fd(kbot:ktop) + fd2(kbot:ktop) )
   endif
!   B = N_IC/sum(f(kbot:ktop)*p_average(kbot:ktop))

   DO k=kbot,ktop
     if ( cellcount(k) .gt. 0. ) THEN
!       delta = term2*B*fd(k)*t_average(k)*ic_fr/cellcount(k)/dz(k)/100.
        !*  implementation note: 1) ic_fr * N_IC/cellcount gives moles of NO in the column
        !*                       2) Multiplying by fd gives the # moles per level
        !*                       3) Divide by moles of air gives mixing ratio
        !*                       4) Multiply by 1E6 gives ppmv
        delta = (ic_fr * N_IC / cellcount(k)) * fd(k) / (molesofair(k)*dz(k)) * 1E6
        where(refl(ips:ipe,k,jps:jpe) .gt. ref_threshold )
          lnox_ic_tend(ips:ipe,k,jps:jpe) = delta
        endwhere
     ENDIF
   ENDDO

 ENDIF ! IC lightning

!-----------------------------------------------------------------
! Calculates CG distribution
 IF ((cg_fr .gt. 0 ) .and. (( ktop .gt. k_lower ) .and. (kbtm .lt. ktop) ) ) THEN
   call bellcurve(kts,ktop,klower,z_average, kps,kpe, fd, dz)
!   B = N_CG/(sum(fd(kts:ktop)*p_average(kts:ktop)))

!  delta = term2*B*fd(k)*t_average(k)*cg_fr/cellcount(kbtm)/dz(k)/100.
   delta = (cg_fr * N_CG / cellcount(kbtm)) * fd(k) / (molesofair(k)*dz(k)) * 1E6
   ! below cloud
   !*  implementation note: previously, between kts & kbtm is practically excluded because of
   !*                       the test for refl > 20
   !*  implementation note: using "where" here would end up reevaluating refl(ips:ipe,kbtm,jps:jpe)
   !*                       multiple times because of the kts,kbtm-1 loop.
   do i=ips,ipe
     do j=jps,jpe
       if ( refl(i,kbtm,j) .gt. ref_threshold ) THEN
         lnox_cg_tend(i,kts:(kbtm-1),j) = delta
       ENDIF
     enddo
   enddo

   ! within cloud
   do k = kbtm,ktop
    !delta = term2*B*fd(k)*t_average(k)*cg_fr/cellcount(k)/dz(k)/100.
     delta = (cg_fr * N_CG / cellcount(k)) * fd(k) / (molesofair(k)*dz(k)) * 1E6
     where( refl(ips:ipe,k,jps:jpe) ) .gt. ref_threshold )
       lnox_cg_tend(ips:ipe,k,jps:jpe) = delta
     endwhere
   enddo
 ENDIF

 END SUBROUTINE lightning_nox_decaria



!************************************************************************
! This subroutine prepares a normal distribution between k_min and
! k_max centered at k_mu. Distribution for each level is
! normalized to \int^{z_at_w(k_max}_{z_at_w(k_min-1)}f(z)dz
!
! Unit of f is fraction of total column
!
! Modified from v3.4.1 module_ltng_crm.F, kept the math but changed
! the implementation for better clarity. Removed patch-wide averaging
! of z.
!************************************************************************

 SUBROUTINE bellcurve ( k_min, k_max, k_mu, z, kps,kpe, f, dz )
!-----------------------------------------------------------------

 INTEGER,                      INTENT(IN   ) :: k_min, k_max, k_mu
 REAL,   DIMENSION( kps:kpe ), INTENT(IN   ) :: z       ! at phy
 INTEGER,                      INTENT(IN   ) :: kps,kpe

 REAL,   DIMENSION( kps:kpe ), INTENT(  OUT) :: f, dz

 INTEGER :: i,j,k
 REAL, DIMENSION( kps:kpe ) :: ex
 REAL :: sigma, z_mu, cuml_f_dist
 REAL, PARAMETER :: two_pi = 6.2831854

!-----------------------------------------------------------------

 f(kps:kpe) = 0.
 z_mu = z(k_mu)
 sigma = AMIN1(z(k_max)-z_mu,z_mu-z(k_min))/3.0

 ! distance from mean
 ex(k_min:k_max) = (z(k_min:k_max)-z_mu)/sigma
 
 ! Truncated Gaussian at 3 sigma
 f(k_min:k_max) = (1.0/(sqrt(two_pi)*sigma))*exp(-ex(k_min:k_max)*ex(k_min:k_max)/2.0)

 dz(kps) = 0. ! safe as long as k_min != kps
 dz(kpe) = 0. ! safe as long as k_max != kpe
 DO k=kps+1,kpe-1
!  dz(k) = (z(k+1)+z(k))/2. - (z(k)+z(k-1))/2.
   dz(k) = (z(k+1) - z(k-1))*.5
 ENDDO

 ! Normalize
 cuml_f_dist = sum(f(k_min:k_max) * dz(k_min:k_max))
 f(k_min:k_max) = f(k_min:k_max)/cuml_f_dist

 END SUBROUTINE bellcurve

!************************************************************************
! This subroutine finds the vertical indices (phy grid) of the follow
! within a column:
! 1) ktop - reflectivity cloud top
! 2) kbtm - reflectivity cloud bottom
! 3) kupper - upper mode isotherm
! 3) klower - lower mode isotherm
!************************************************************************
 SUBROUTINE kfind ( &
              ! Prognostics
                refl, t,                              &
              ! Settings
                refl_threshold,                       &
              ! Namelist settings
                ltng_temp_upper,ltng_temp_lower,      &
                cellcount_method,                     &
              ! Order dependent args for domain, mem, and tile dims
                ids, ide, jds, jde, kds, kde,              &
                ims, ime, jms, jme, kms, kme,              &
                ips, ipe, jps, jpe, kps, kpe,              &
              ! Outputs
                ktop,kbtm,kupper,klower
            )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
!-----------------------------------------------------------------

! Prognostics
 REAL, DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: t, refl

! Namelist settings
 REAL,    INTENT(IN   )    ::       ltng_temp_lower, ltng_temp_upper
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Outputs
 INTEGER, INTENT(  OUT)    ::       ktop,kbtm,kupper,klower

! Local vars
 REAL    :: ktop_r, kbtm_r, kupper_r, klower_r
 INTEGER :: k, k_tmin
 REAL    :: tmin

!-----------------------------------------------------------------
 ktop = kps
 kbtm = kps
 kupper = kps
 klower = kps

 ! look for ktop
 k = kte
 DO WHILE (maxval(refl(ips:ipe,k,jps:jpe)) .lt. refl_threshold .and. k .gt. kts )
  k = k-1
 ENDDO
 ktop = k

 ! Look for kbtm
 k = kts
 DO WHILE( maxval(refl(ips:ipe,k,jps:jpe)) .lt. refl_threshold .and. k .lt. ktop )
  k = k+1
 ENDDO
 kbtm = k
 if (kbtm .eq. kts) kbtm = kte

 ! Look for tropopause
 k_tmin = kts
 tmin = t(kts)
 DO k=kts+1,kte
   IF ( minval(t(ips:ipe,k,jps:jpe)) .lt. tmin) THEN
     tmin = t(k)
     k_tmin = k
   ENDIF
 ENDDO

 ! Look for kupper
 k = k_tmin
 DO WHILE ( minval(t(ips:ipe,k,jps:jpe)) .lt. ltng_temp_upper + 273.15 )
   k = k - 1
 ENDDO
 kupper = k

 ! Look for klower
 DO WHILE( minval(t(ips:ipe,k,jps:jpe)) .lt. ltng_temp_lower + 273.15 )
   k = k - 1
 ENDDO
 klower = k
 
 IF ( cellcount_method .eq. 2 ) THEN
   kbot_r = real(kbot)
   ktop_r = real(ktop)
   klower_r = real(klower)
   kupper_r = real(kupper)
   kbot = nint(wrf_dm_min_real(kbot_r))
   ktop = nint(wrf_dm_max_real(ktop_r))
   klower = nint(wrf_dm_max_real(k_ower_r))
   kuppwer = nint(wrf_dm_max_real(kupper_r))
 endif

 END SUBROUTINE kfind

!************************************************************************
! This function averages out a 3D array into 1D in the 2nd dimension
!************************************************************************
 SUBROUTINE horizontalAverage( array3D, ips, ipe, kps, kpe, jps, jpe, array1D )
!-----------------------------------------------------------------
 REAL, DIMENSION(ips:ipe,kps:kpe,jps:jpe), INTENT(IN) :: array3D
 INTEGER, INTENT(IN) :: ips,ipe,kps,kpe,jps,jpe

 REAL, DIMENSION(kps:kpe), INTENT(OUT) :: array1D
!-----------------------------------------------------------------
 DO k=kps,kpe
   array1D(k) = sum(array3D(ips:ipe,k,jps:jpe))/((ipe-ips+1)*(jpe-jps+1))
 ENDDO
    
 END SUBROUTINE horizontalAverage

 END MODULE module_lightning_nox_decaria