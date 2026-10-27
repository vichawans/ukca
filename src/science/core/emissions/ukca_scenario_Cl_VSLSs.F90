! *****************************COPYRIGHT*******************************
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
! Description:
! Copied and adapted from ukca_scenario_wmo1a and ukca_scenario_rcp modules,
! but adapted for four Cl-VSLS, with latitude-dependant LBCs.
! 08/2018 E.M. Bednarz
!------------------------------------------------------------------
!
MODULE ukca_scenario_Cl_VSLSs_mod

USE umPrintMgr, ONLY: umPrint, umMessage, PrStatus_Diag, PrintStatus, &
                      PrStatus_Normal, PrStatus_Oper, newline
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_SCENARIO_Cl_VSLSs'

PUBLIC :: ukca_scenario_CL_VSLSs
PUBLIC :: test_scenario_cl_vsls_ctl

! logical to turn on testing
LOGICAL, PUBLIC, PARAMETER :: L_ukca_test_scenario_cl_vsls=.FALSE.

CONTAINS

SUBROUTINE ukca_scenario_Cl_VSLSs(n_lbc_specs_2D, lbc_specs_2D,   &
                      i_year, i_day_number,                       &
                      lbc_mmr_2D)

USE ukca_constants,   ONLY: c_ch2cl2, c_chcl3, c_c2cl4, c_c2h4cl2
!USE yomhook,          ONLY: lhook, dr_hook
USE ereport_mod,      ONLY: ereport
USE nlsizes_namelist_mod, ONLY: rows
!USE missing_data_mod, ONLY: rmdi
USE trignometric_mod, ONLY: true_latitude
USE conversions_mod, ONLY: recip_pi_over_180
USE ukca_option_mod, ONLY: l_ukca_fix_clvsls, clvsls_fix_year, l_ukca_dest_clvsls_lbc,       &
                           ukca_clvsls_specify_method,            &
                           userval_ch2cl2, userval_chcl3, userval_c2cl4,                    &
                           userval_c2h4cl2

IMPLICIT NONE


INTEGER, INTENT(IN) :: n_lbc_specs_2D            ! Number of lat-dependent lower BC species
INTEGER, INTENT(IN) :: i_year                    ! year
INTEGER, INTENT(IN) :: i_day_number                ! number of days since start of year
REAL, INTENT(INOUT) :: lbc_mmr_2D(rows,n_lbc_specs_2D)  ! Lower lat-dependent BC mass mixing ratios
CHARACTER(LEN=10), INTENT(IN) :: lbc_specs_2D(n_lbc_specs_2D) ! LBC species

LOGICAL, SAVE :: first = .TRUE.

INTEGER,SAVE :: ich2cl2=0, ichcl3=0, ic2cl4=0, ic2h4cl2=0

INTEGER, PARAMETER :: vsls_nyear = 20
REAL, SAVE :: vsls_year(vsls_nyear) =                             &
 (/2000.5,2001.5,2002.5,2003.5,2004.5,2005.5,2006.5,2007.5,2008.5,&
   2009.5,2010.5,2011.5,2012.5,2013.5,2014.5,2015.5,2016.5,2017.5,&
   2018.5,2019.5/)
REAL :: this_latitude_deg

INTEGER, PARAMETER :: n_full_spec_2D=4

REAL, SAVE :: full_lbc_2D_lat1(n_full_spec_2D)
REAL, SAVE :: full_lbc_2D_lat2(n_full_spec_2D)
REAL, SAVE :: full_lbc_2D_lat3(n_full_spec_2D)
REAL, SAVE :: full_lbc_2D_lat4(n_full_spec_2D)
REAL, SAVE :: full_lbc_2D_lat5(n_full_spec_2D)


REAL, SAVE :: ch2cl2_lat1(vsls_nyear)
REAL, SAVE :: ch2cl2_lat2(vsls_nyear)
REAL, SAVE :: ch2cl2_lat3(vsls_nyear)
REAL, SAVE :: ch2cl2_lat4(vsls_nyear)
REAL, SAVE :: ch2cl2_lat5(vsls_nyear)

REAL, SAVE :: chcl3_lat1(vsls_nyear)
REAL, SAVE :: chcl3_lat2(vsls_nyear)
REAL, SAVE :: chcl3_lat3(vsls_nyear)
REAL, SAVE :: chcl3_lat4(vsls_nyear)
REAL, SAVE :: chcl3_lat5(vsls_nyear)

REAL, SAVE :: c2cl4_lat1(vsls_nyear)
REAL, SAVE :: c2cl4_lat2(vsls_nyear)
REAL, SAVE :: c2cl4_lat3(vsls_nyear)
REAL, SAVE :: c2cl4_lat4(vsls_nyear)
REAL, SAVE :: c2cl4_lat5(vsls_nyear)

REAL, SAVE :: c2h4cl2_lat1
REAL, SAVE :: c2h4cl2_lat2
REAL, SAVE :: c2h4cl2_lat3
REAL, SAVE :: c2h4cl2_lat4
REAL, SAVE :: c2h4cl2_lat5


! local variables
REAL :: ftime
INTEGER :: yindex
REAL :: frac
INTEGER  :: myyear
INTEGER :: i
INTEGER :: mylat
INTEGER :: ispecies

! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage

!INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
!INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
!REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_SCENARIO_Cl_VSLSs'

!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!!!!!!!!! --------------------------------------------------------------------
!!! EMB: input data as vmr
IF (first) THEN

!lat1=90-30S
!lat2=30-0S
!lat3=0-30N
!lat4=30-60N
!lat5=60-90N


  ch2cl2_lat1(:) = (/9.5E-12,  9.4E-12,  9.2E-12,  9.1E-12,  9.5E-12,  9.6E-12,  9.9E-12, &
                    11.1E-12, 11.9E-12, 11.7E-12, 12.6E-12, 13.5E-12, 13.7E-12, 15.6E-12, &
                    17.0E-12, 16.5E-12, 16.7E-12, 17.2E-12, 18.8E-12, 18.4E-12/)

  ch2cl2_lat2(:) = (/10.7E-12, 10.5E-12, 10.2E-12, 10.3E-12, 10.7E-12, 11.6E-12, 13.4E-12, &
                     13.6E-12, 14.2E-12, 14.0E-12, 14.9E-12, 15.6E-12, 16.4E-12, 19.3E-12, &
                     20.1E-12, 20.3E-12, 20.1E-12, 22.3E-12, 24.2E-12, 22.3E-12/)

  ch2cl2_lat3(:) = (/23.5E-12, 24.7E-12, 24.4E-12, 26.6E-12, 25.9E-12, 27.1E-12, 28.4E-12, &
                     32.5E-12, 34.2E-12, 35.3E-12, 41.1E-12, 39.8E-12, 43.4E-12, 53.0E-12, &
                     56.5E-12, 55.3E-12, 57.8E-12, 60.0E-12, 59.7E-12, 66.1E-12/)

  ch2cl2_lat4(:) = (/36.3E-12, 34.0E-12, 35.7E-12, 36.0E-12, 33.8E-12, 36.0E-12, 38.1E-12, &
                     41.6E-12, 43.0E-12, 42.8E-12, 49.0E-12, 48.7E-12, 51.1E-12, 62.6E-12, &
                     63.8E-12, 62.5E-12, 64.6E-12, 70.2E-12, 68.9E-12, 72.5E-12/)

  ch2cl2_lat5(:) = (/35.2E-12, 33.0E-12, 32.9E-12, 33.9E-12, 32.6E-12, 33.1E-12, 35.3E-12, &
                     37.7E-12, 41.0E-12, 41.4E-12, 46.4E-12, 46.4E-12, 48.0E-12, 59.3E-12, &
                     61.0E-12, 60.1E-12, 62.4E-12, 68.2E-12, 67.0E-12, 70.4E-12/) 
!!!!! ---------------------------------------------

  chcl3_lat1(:) = (/5.8E-12, 5.6E-12, 5.6E-12, 5.5E-12, 5.6E-12, 5.5E-12, 5.3E-12, &
                    5.6E-12, 5.5E-12, 5.4E-12, 5.4E-12, 5.6E-12, 5.7E-12, 5.4E-12, &
                    5.9E-12, 6.0E-12, 6.0E-12, 5.9E-12, 6.0E-12, 6.0E-12/)

  chcl3_lat2(:) = (/5.2E-12, 5.0E-12, 4.9E-12, 5.0E-12, 4.7E-12, 4.9E-12, 4.8E-12, &
                    5.1E-12, 4.8E-12, 4.8E-12, 4.7E-12, 5.0E-12, 5.4E-12, 5.2E-12, &
                    5.5E-12, 5.6E-12, 5.4E-12, 6.8E-12, 6.1E-12, 6.1E-12/)

  chcl3_lat3(:) = (/7.6E-12, 7.0E-12, 7.1E-12, 7.2E-12, 6.9E-12, 6.9E-12, 7.0E-12, &
                    7.1E-12, 7.3E-12, 6.7E-12, 7.0E-12, 7.2E-12, 7.5E-12, 8.0E-12, &
                    9.1E-12, 8.0E-12, 8.9E-12, 10.6E-12, 9.9E-12, 9.9E-12/)

  chcl3_lat4(:) = (/11.6E-12, 10.6E-12, 10.4E-12, 10.3E-12, 10.4E-12, 10.6E-12, 10.6E-12, &
                    10.6E-12, 10.8E-12, 10.7E-12, 11.3E-12, 11.0E-12, 11.2E-12, 12.1E-12, &
                    13.4E-12, 13.9E-12, 14.2E-12, 14.9E-12, 14.4E-12, 14.4E-12/)

  chcl3_lat5(:) = (/10.8E-12, 10.8E-12, 10.6E-12, 10.9E-12, 10.4E-12, 10.3E-12, 10.3E-12, &
                    10.4E-12, 10.2E-12, 10.6E-12, 11.2E-12, 11.6E-12, 11.9E-12, 12.4E-12, &
                    13.3E-12, 13.2E-12, 14.0E-12, 15.1E-12, 14.6E-12, 14.6E-12/) 

!!!!! ---------------------------------------------
! EMB values of C2Cl4 for 2008 added after missing - 05/09/18

  c2cl4_lat1(:) = (/0.8E-12, 0.7E-12, 0.7E-12, 0.6E-12, 0.6E-12, 0.6E-12, 0.5E-12, 0.5E-12, 0.6E-12,&
                    0.5E-12, 0.5E-12, 0.4E-12, 0.4E-12, 0.6E-12, 0.5E-12, 0.4E-12, 0.4E-12, 0.4E-12,&
                    0.4E-12, 0.4E-12/)

  c2cl4_lat2(:) = (/0.9E-12, 0.8E-12, 0.8E-12, 0.7E-12, 0.7E-12, 0.6E-12, 0.7E-12, 0.6E-12, 0.8E-12,&
                    0.7E-12, 0.6E-12, 0.5E-12, 0.5E-12, 0.6E-12, 0.6E-12, 0.6E-12, 0.5E-12, 0.6E-12,&
                    0.6E-12, 0.5E-12/)

  c2cl4_lat3(:) = (/2.8E-12, 2.8E-12, 2.6E-12, 2.6E-12, 2.3E-12, 2.2E-12, 2.2E-12, 2.0E-12, 2.0E-12,&
                    1.9E-12, 1.9E-12, 1.6E-12, 1.5E-12, 1.6E-12, 1.7E-12, 1.6E-12, 1.6E-12, 1.6E-12,&
                    1.5E-12, 1.6E-12/)

  c2cl4_lat4(:) = (/7.3E-12, 7.1E-12 ,6.4E-12, 5.8E-12, 5.2E-12, 5.6E-12, 5.9E-12, 5.5E-12, 4.6E-12,&
                    4.3E-12, 4.3E-12, 3.7E-12, 3.3E-12, 3.6E-12, 3.4E-12, 3.4E-12, 3.3E-12, 3.4E-12,&
                    3.0E-12, 3.2E-12/)

  c2cl4_lat5(:) = (/5.0E-12, 4.5E-12, 4.6E-12, 4.4E-12, 3.9E-12, 3.7E-12, 3.6E-12, 3.3E-12, 3.4E-12,&
                    2.8E-12, 2.9E-12, 2.5E-12, 2.3E-12, 2.3E-12, 2.4E-12, 2.5E-12, 2.4E-12, 2.4E-12,&
                    2.1E-12, 2.2E-12/)

!!!!! ---------------------------------------------

  c2h4cl2_lat1 =  2.0E-12
  c2h4cl2_lat2 =  3.0E-12
  c2h4cl2_lat3 = 10.0E-12
  c2h4cl2_lat4 = 15.0E-12
  c2h4cl2_lat5 = 15.0E-12

  first = .FALSE.
END IF

  ich2cl2=0
  ichcl3=0
  ic2cl4=0
  ic2h4cl2=0

  DO i=1,n_lbc_specs_2D
    SELECT CASE (lbc_specs_2D(i))
    CASE ('CH2Cl2    ')
      ich2cl2=i
    CASE ('CHCl3     ')
      ichcl3=i
    CASE ('C2Cl4     ')
      ic2cl4=i
    CASE ('C2H4Cl2   ')
      ic2h4cl2=i
    END SELECT
  END DO


!EMB: this paragraph was copied and adapted from ukca_scenario_rcp
!
! work out values to interpolate to in time
! currently the routine will give values according to the model
! year and day - however, if the model year is outside the range
! in the file then the output from the routine will be as the
! first/last point would be, i.e.
! if the file goes from 1750-2500 then years <1750 would have the
! 1750 value, and years >2500 would have the 2500 value
ierr = 0
IF ((vsls_nyear  >   1) .AND. (l_ukca_dest_clvsls_lbc)) THEN
!!! Constant CL-VSLSs values
  IF (l_ukca_fix_clvsls) THEN
!!! SET YEAR
    IF (ukca_clvsls_specify_method == 1) THEN
      myyear=clvsls_fix_year
      yindex=myyear-1999
      frac=0.0
!!
!!! USER-DEFINED Cl-VSLSs 
    ELSE IF (ukca_clvsls_specify_method == 2) THEN
      full_lbc_2D_lat1(1) = userval_ch2cl2(1)
      full_lbc_2D_lat2(1) = userval_ch2cl2(2)
      full_lbc_2D_lat3(1) = userval_ch2cl2(3)
      full_lbc_2D_lat4(1) = userval_ch2cl2(4)
      full_lbc_2D_lat5(1) = userval_ch2cl2(5)

      full_lbc_2D_lat1(2) = userval_chcl3(1)
      full_lbc_2D_lat2(2) = userval_chcl3(2)
      full_lbc_2D_lat3(2) = userval_chcl3(3)
      full_lbc_2D_lat4(2) = userval_chcl3(4)
      full_lbc_2D_lat5(2) = userval_chcl3(5)

      full_lbc_2D_lat1(3) = userval_c2cl4(1)
      full_lbc_2D_lat2(3) = userval_c2cl4(2)
      full_lbc_2D_lat3(3) = userval_c2cl4(3)
      full_lbc_2D_lat4(3) = userval_c2cl4(4)
      full_lbc_2D_lat5(3) = userval_c2cl4(5)

      full_lbc_2D_lat1(4) = userval_c2h4cl2(1)
      full_lbc_2D_lat2(4) = userval_c2h4cl2(2)
      full_lbc_2D_lat3(4) = userval_c2h4cl2(3)
      full_lbc_2D_lat4(4) = userval_c2h4cl2(4)
      full_lbc_2D_lat5(4) = userval_c2h4cl2(5)
    END IF
!!! Time-varying Cl-VSLSs  
  ELSE
    ftime = REAL(i_year) + REAL(i_day_number - 1)/360.0
    IF (ftime < vsls_year(1)) THEN
      yindex = 1
      frac = 0.0
      IF (PrintStatus >= PrStatus_Diag) THEN
        ! print warning message to this effect
        ierr = -1*INT(vsls_year(1))
        WRITE(cmessage,FMT='(A,F12.6)')                          &
             'Warning, model year is less than ', vsls_year(1)
        CALL ereport('UKCA_SCENARIO_Cl_VSLSs',ierr,cmessage)
      END IF
    ELSE IF (ftime >= vsls_year(vsls_nyear)) THEN
      yindex = vsls_nyear - 1
      frac = 1.0
      IF (PrintStatus >= PrStatus_Diag) THEN
        ! print warning message to this effect
        ierr = -1*INT(vsls_year(vsls_nyear))
        WRITE(cmessage,FMT='(A,F12.6)')                          &
             'Warning, model year is greater than ', vsls_year(vsls_nyear)
        CALL ereport('UKCA_SCENARIO_Cl_VSLSs',ierr,cmessage)
      END IF
    ELSE
      yindex = 1
      DO WHILE (ftime > vsls_year(yindex + 1))
        yindex = yindex + 1
      END DO
      frac = (ftime-vsls_year(yindex))/(vsls_year(yindex+1)-vsls_year(yindex))
    END IF
  ENDIF

  IF ((.NOT. l_ukca_fix_clvsls) .OR. (ukca_clvsls_specify_method == 1)) THEN
    ! CH2Cl2
    full_lbc_2D_lat1(1) = (1.0 - frac)* ch2cl2_lat1(yindex)+frac* ch2cl2_lat1(yindex+1)
    full_lbc_2D_lat2(1) = (1.0 - frac)* ch2cl2_lat2(yindex)+frac* ch2cl2_lat2(yindex+1)
    full_lbc_2D_lat3(1) = (1.0 - frac)* ch2cl2_lat3(yindex)+frac* ch2cl2_lat3(yindex+1)
    full_lbc_2D_lat4(1) = (1.0 - frac)* ch2cl2_lat4(yindex)+frac* ch2cl2_lat4(yindex+1)
    full_lbc_2D_lat5(1) = (1.0 - frac)* ch2cl2_lat5(yindex)+frac* ch2cl2_lat5(yindex+1)
    ! CHCl3
    full_lbc_2D_lat1(2) = (1.0 - frac)* chcl3_lat1(yindex)+frac* chcl3_lat1(yindex+1)
    full_lbc_2D_lat2(2) = (1.0 - frac)* chcl3_lat2(yindex)+frac* chcl3_lat2(yindex+1)
    full_lbc_2D_lat3(2) = (1.0 - frac)* chcl3_lat3(yindex)+frac* chcl3_lat3(yindex+1)
    full_lbc_2D_lat4(2) = (1.0 - frac)* chcl3_lat4(yindex)+frac* chcl3_lat4(yindex+1)
    full_lbc_2D_lat5(2) = (1.0 - frac)* chcl3_lat5(yindex)+frac* chcl3_lat5(yindex+1)
    ! C2Cl4
    full_lbc_2D_lat1(3) = (1.0 - frac)* c2cl4_lat1(yindex)+frac* c2cl4_lat1(yindex+1)
    full_lbc_2D_lat2(3) = (1.0 - frac)* c2cl4_lat2(yindex)+frac* c2cl4_lat2(yindex+1)
    full_lbc_2D_lat3(3) = (1.0 - frac)* c2cl4_lat3(yindex)+frac* c2cl4_lat3(yindex+1)
    full_lbc_2D_lat4(3) = (1.0 - frac)* c2cl4_lat4(yindex)+frac* c2cl4_lat4(yindex+1)
    full_lbc_2D_lat5(3) = (1.0 - frac)* c2cl4_lat5(yindex)+frac* c2cl4_lat5(yindex+1)
    ! C2H4Cl2
    full_lbc_2D_lat1(4)=c2h4cl2_lat1
    full_lbc_2D_lat2(4)=c2h4cl2_lat2
    full_lbc_2D_lat3(4)=c2h4cl2_lat3
    full_lbc_2D_lat4(4)=c2h4cl2_lat4
    full_lbc_2D_lat5(4)=c2h4cl2_lat5
  END IF

  ! Assign values to LBC array, and convert from vmr to mmr
  ! Default values to 0
  lbc_mmr_2D(:,:) = 0.0

  DO mylat=1,rows
    this_latitude_deg=true_latitude(1,mylat)*recip_pi_over_180

    IF (ich2cl2     > 0) THEN
      IF ((this_latitude_deg .GE. -90) .AND. (-30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ich2cl2    ) = full_lbc_2D_lat1(1)*c_ch2cl2
 
      ELSE IF ((this_latitude_deg .GE. -30) .AND. (0 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ich2cl2    ) = full_lbc_2D_lat2(1)*c_ch2cl2

      ELSE IF ((this_latitude_deg .GE. 0) .AND. (30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ich2cl2    ) = full_lbc_2D_lat3(1)*c_ch2cl2

      ELSE IF ((this_latitude_deg .GE. 30) .AND. (60 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ich2cl2    ) = full_lbc_2D_lat4(1)*c_ch2cl2

      ELSE IF ((this_latitude_deg .GE. 60) .AND. (90 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ich2cl2    ) = full_lbc_2D_lat5(1)*c_ch2cl2

      ELSE
        cmessage='EMB INPUTTING LBCs FOR VSLSs: unsuported grid!!!'
        CALL ereport('UKCA_SCENARIO_Cl_VSLSs',ierr,cmessage)
      END IF
    END IF

    IF (ichcl3      > 0) THEN
      IF ((this_latitude_deg .GE. -90) .AND. (-30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ichcl3     ) = full_lbc_2D_lat1(2)*c_chcl3

      ELSE IF ((this_latitude_deg .GE. -30) .AND. (0 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ichcl3     ) = full_lbc_2D_lat2(2)*c_chcl3

      ELSE IF ((this_latitude_deg .GE. 0) .AND. (30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ichcl3     ) = full_lbc_2D_lat3(2)*c_chcl3

      ELSE IF ((this_latitude_deg .GE. 30) .AND. (60 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ichcl3     ) = full_lbc_2D_lat4(2)*c_chcl3

      ELSE IF ((this_latitude_deg .GE. 60) .AND. (90 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ichcl3     ) = full_lbc_2D_lat5(2)*c_chcl3

      ELSE
        cmessage='EMB INPUTTING LBCs FOR VSLSs: unsuported grid!!!'
        CALL ereport('UKCA_SCENARIO_Cl_VSLSs',ierr,cmessage)
      END IF
    END IF


    IF (ic2cl4      > 0) THEN
      IF ((this_latitude_deg .GE. -90) .AND. (-30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2cl4     ) = full_lbc_2D_lat1(3)*c_c2cl4

      ELSE IF ((this_latitude_deg .GE. -30) .AND. (0 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2cl4     ) = full_lbc_2D_lat2(3)*c_c2cl4

      ELSE IF ((this_latitude_deg .GE. 0) .AND. (30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2cl4     ) = full_lbc_2D_lat3(3)*c_c2cl4

      ELSE IF ((this_latitude_deg .GE. 30) .AND. (60 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2cl4     ) = full_lbc_2D_lat4(3)*c_c2cl4

      ELSE IF ((this_latitude_deg .GE. 60) .AND. (90 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2cl4     ) = full_lbc_2D_lat5(3)*c_c2cl4

      ELSE
        cmessage='EMB INPUTTING LBCs FOR VSLSs: unsuported grid!!!'
        CALL ereport('UKCA_SCENARIO_Cl_VSLSs',ierr,cmessage)
      END IF
    END IF



    IF (ic2h4cl2    > 0) THEN
      IF ((this_latitude_deg .GE. -90) .AND. (-30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2h4cl2   ) = full_lbc_2D_lat1(4)*c_c2h4cl2

      ELSE IF ((this_latitude_deg .GE. -30) .AND. (0 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2h4cl2   ) = full_lbc_2D_lat2(4)*c_c2h4cl2

      ELSE IF ((this_latitude_deg .GE. 0) .AND. (30 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2h4cl2   ) = full_lbc_2D_lat3(4)*c_c2h4cl2

      ELSE IF ((this_latitude_deg .GE. 30) .AND. (60 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2h4cl2   ) = full_lbc_2D_lat4(4)*c_c2h4cl2

      ELSE IF ((this_latitude_deg .GE. 60) .AND. (90 .GT. this_latitude_deg)) THEN
        lbc_mmr_2D(mylat,ic2h4cl2   ) = full_lbc_2D_lat5(4)*c_c2h4cl2

      ELSE
        cmessage='EMB INPUTTING LBCs FOR VSLSs: unsuported grid!!!'
        CALL ereport('UKCA_SCENARIO_Cl_VSLSs',ierr,cmessage)
      END IF
    END IF

  END DO

!!! debuging 

  DO mylat=1,rows !! format (A,F4.2,A,I1,E9.4) ?
    WRITE(umMessage,*) 'lbc_mmr_2D at',(true_latitude(1,mylat)*recip_pi_over_180),'deg;',ich2cl2,'=',lbc_mmr_2D(mylat,ich2cl2    )
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

    WRITE(umMessage,*) 'lbc_mmr_2D at',(true_latitude(1,mylat)*recip_pi_over_180),'deg;',ichcl3,'=',lbc_mmr_2D(mylat,ichcl3     )
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

    WRITE(umMessage,*) 'lbc_mmr_2D at',(true_latitude(1,mylat)*recip_pi_over_180),'deg;',ic2cl4,'=',lbc_mmr_2D(mylat,ic2cl4     )
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

    WRITE(umMessage,*) 'lbc_mmr_2D at',(true_latitude(1,mylat)*recip_pi_over_180),'deg;',ic2h4cl2,'=',lbc_mmr_2D(mylat,ic2h4cl2   )
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')
  END DO
!
  DO ispecies=1,4
    WRITE(umMessage,*) 'full_lbc_2D_lat1(',ispecies,')=',full_lbc_2D_lat1(ispecies)
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

    WRITE(umMessage,*) 'full_lbc_2D_lat2(',ispecies,')=',full_lbc_2D_lat2(ispecies)
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

    WRITE(umMessage,*) 'full_lbc_2D_lat3(',ispecies,')=',full_lbc_2D_lat3(ispecies)
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

    WRITE(umMessage,*) 'full_lbc_2D_lat4(',ispecies,')=',full_lbc_2D_lat4(ispecies)
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

    WRITE(umMessage,*) 'full_lbc_2D_lat5(',ispecies,')=',full_lbc_2D_lat5(ispecies)
    CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')
  END DO

ELSE IF (.NOT. l_ukca_dest_clvsls_lbc) THEN
  DO i=1,n_lbc_specs_2D
    lbc_mmr_2D(:,i) = 0.0e+00
  END DO
END IF

WRITE(umMessage,*) 'n_lbc_specs_2D:',n_lbc_specs_2D,'; lbc_specs_2D:',lbc_specs_2D
CALL umPrint(umMessage,src='ukca_scenario_Cl_VSLSs')

!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_scenario_Cl_VSLSs

!!! -------------------------------------------------------------------

SUBROUTINE test_scenario_cl_vsls_ctl(n_lbc_specs_2D, lbc_specs_2D)
!USE yomhook,          ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN)      :: n_lbc_specs_2D          ! Number of LBC species
CHARACTER(LEN=10), INTENT(IN) :: lbc_specs_2D(n_lbc_specs_2D)! LBC species

! loop counters
INTEGER :: i
! to simulate many years to calculate the lbc values
INTEGER :: jday
INTEGER :: jyear

CHARACTER(LEN=256), PARAMETER :: ukca_test_file_vsls1=                 &
                                 'Test_Cl_VSLS_UKCA_2D.dat'
CHARACTER(LEN=256), PARAMETER :: ukca_test_file_vsls2=                 &
                                 'Test_Cl_VSLS_UKCA_SHmidhigh_mmr.dat'
CHARACTER(LEN=256), PARAMETER :: ukca_test_file_vsls3=                 &
                                 'Test_Cl_VSLS_UKCA_SHmidhigh_vmr.dat'
CHARACTER(LEN=256) :: file_header

!INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
!INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
!REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TEST_SCENARIO_Cl_VSLS_CTL'

! Body of subroutine
!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! this routine calls the routine which tests the
! the scenario routine by
! cycling through all 1960-2100 and writing the
! lbc_mmr values out to a file so that they can be
! plotted and compared to the ones in the original
! specification file(s)
! this has been separated off, as we want to call
! twice - once for the UKCA lbc_specs and once
! for the full listing of lbc_specs

! check UKCA values
file_header = 'LBC_Cl_VSLS MMR VALUES REQUESTED FOR USE IN UKCA'
CALL test_scenario_cl_vsls(n_lbc_specs_2D, lbc_specs_2D,          &
                       ukca_test_file_vsls1, ukca_test_file_vsls2,&
                       ukca_test_file_vsls3, file_header)


!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE test_scenario_cl_vsls_ctl

!!! -------------------------------------------------------------------------------------

SUBROUTINE test_scenario_cl_vsls(n_lbc_specs_2D, lbc_specs_2D,      &
                             test_file_vsls1, test_file_vsls2,     &
                             test_file_vsls3, file_header)
USE ukca_constants,   ONLY:  c_ch2cl2,               &
    c_chcl3, c_c2cl4, c_c2h4cl2
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE ereport_mod,      ONLY: ereport
!USE yomhook,          ONLY: lhook, dr_hook
USE nlsizes_namelist_mod, ONLY: rows
USE trignometric_mod, ONLY: true_latitude
USE conversions_mod, ONLY: recip_pi_over_180
IMPLICIT NONE

INTEGER, INTENT(IN)      :: n_lbc_specs_2D          ! Number of LBC species
CHARACTER(LEN=10), INTENT(IN) :: lbc_specs_2D(n_lbc_specs_2D)! LBC species
CHARACTER(LEN=256), INTENT(IN) :: test_file_vsls1      ! file to write to
CHARACTER(LEN=256), INTENT(IN) :: test_file_vsls2      ! file to write to
CHARACTER(LEN=256), INTENT(IN) :: test_file_vsls3      ! file to write to
CHARACTER(LEN=256), INTENT(IN) :: file_header     ! header info for file

REAL, ALLOCATABLE :: lbc_mmr_2D_test(:,:)  ! Lower BC mass mixing ratios

INTEGER :: file_unit
! loop counters
INTEGER :: i
INTEGER :: j
! to simulate many years to calculate the lbc values
INTEGER :: jday
INTEGER :: jyear
INTEGER :: mylat

REAL, ALLOCATABLE :: lbc_c_species_2D_test(:)

! FORMAT STRINGS FOR OUTPUT STATEMENTS - NEED TO BE VARIABLE LENGTH
! format of the name of the strings
CHARACTER(LEN=80)  :: tr_name_fmt_str
! format of the conversion factors
CHARACTER(LEN=80)  :: tr_conv_fmt_str
! format of the mmr values
CHARACTER(LEN=80)  :: tr_float_fmt_str
! format to number the lines for ease of reading
CHARACTER(LEN=80)  :: tr_num_fmt_str


! error handling
INTEGER :: ierr
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=errormessagelength) :: iomessage

!INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
!INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
!REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='TEST_SCENARIO_CL_VSLS'

! Body of subroutine
!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! allocate arrays
!ALLOCATE(lbc_mmr_2D_test(1:rows,1:n_lbc_specs_2D))
ALLOCATE(lbc_c_species_2D_test(1:n_lbc_specs_2D))

! populate lbc_c_species for printing purposes
lbc_c_species_2D_test(:) = 0.0
WHERE (lbc_specs_2D == 'CH2Cl2    ') lbc_c_species_2D_test = c_ch2cl2
WHERE (lbc_specs_2D == 'CHCl3     ') lbc_c_species_2D_test = c_chcl3
WHERE (lbc_specs_2D == 'C2Cl4     ') lbc_c_species_2D_test = c_c2cl4
WHERE (lbc_specs_2D == 'C2H4Cl2   ') lbc_c_species_2D_test = c_c2h4cl2

!! START OF FILE 1 !!!!
ALLOCATE(lbc_mmr_2D_test(1:rows,1:n_lbc_specs_2D))

! find a file unit to write to
CALL assign_file_unit(TRIM(ADJUSTL(test_file_vsls1)), &
                      file_unit, handler="fortran")

! open file to write test output to
ierr = 0
OPEN(UNIT=file_unit,FILE=TRIM(ADJUSTL(test_file_vsls1)),                    &
     POSITION='rewind',STATUS='unknown',IOSTAT=ierr, IOMSG=iomessage)
IF (ierr /= 0) THEN
  cmessage='Error opening file '//TRIM(ADJUSTL(test_file_vsls1)) // ' : '// &
           TRIM(iomessage)
  CALL ereport('TEST_SCENARIO_CL_VSLS',ierr,cmessage)
END IF

! write some header to the file
! start header lines with #
! generate what the format of the header lines should be first
! e.g. "#",X,A10,X,7(A20,X) etc.
! it is necessary to generate these 'on the fly' as there will be a
! variable number of columns dependent on scheme used, or if the output
! is 'full' or just the values UKCA will see (i.e. with lumping).
WRITE(tr_name_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(A20,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_conv_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(E20.10,1X))'
! e.g. 2X,F12.6,7(E20.10,X) etc.
WRITE(tr_float_fmt_str,FMT='(A,I2,A)')                            &
                     '(2X,F12.6,1X,',n_lbc_specs_2D,'(E20.10,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_num_fmt_str,FMT='(A,I2,A)')                              &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(I20,1X))'

! now write header information to the file
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(test_file_vsls1))
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(file_header))
WRITE(file_unit,FMT=tr_name_fmt_str) '#','SPECIES',               &
                    (TRIM(ADJUSTL(lbc_specs_2D(i))),i=1,n_lbc_specs_2D)
WRITE(file_unit,FMT=tr_conv_fmt_str) '#','CONV FAC',              &
                    (lbc_c_species_2D_test(i),i=1,n_lbc_specs_2D)
WRITE(file_unit,FMT=tr_num_fmt_str) '#','LAT (deg)',(i,i=1,n_lbc_specs_2D)

IF (PrintStatus >= Prstatus_Oper) THEN
  WRITE(umMessage,'(A,A)') 'Outputting Cl VSLS LBC data to ',   &
       TRIM(ADJUSTL(test_file_vsls1))
  CALL umPrint(umMessage,src='test_scenario_cl-vsls')
END IF

jyear=1981
jday=241
CALL  ukca_scenario_Cl_VSLSs(n_lbc_specs_2D, lbc_specs_2D,  &
                    jyear, jday,                              &
                    lbc_mmr_2D_test)
DO mylat=1,rows
  WRITE(file_unit,FMT=tr_float_fmt_str)                       &
      (true_latitude(1,mylat)*recip_pi_over_180),         &
      (lbc_mmr_2D_test(mylat,i),i=1,n_lbc_specs_2D)
END DO


! all finished now - close the file
CLOSE(file_unit)
CALL release_file_unit(file_unit, handler="fortran")

! deallocate arrays
IF (ALLOCATED(lbc_mmr_2D_test)) DEALLOCATE(lbc_mmr_2D_test)

!!! END OF FILE 1
!!!
!!! REPEAT FOR FILE 2 !!!!

ALLOCATE(lbc_mmr_2D_test(1:rows,1:n_lbc_specs_2D))

! find a file unit to write to
CALL assign_file_unit(TRIM(ADJUSTL(test_file_vsls2)), &
                      file_unit, handler="fortran")

! open file to write test output to
ierr = 0
OPEN(UNIT=file_unit,FILE=TRIM(ADJUSTL(test_file_vsls2)),                    &
     POSITION='rewind',STATUS='unknown',IOSTAT=ierr, IOMSG=iomessage)
IF (ierr /= 0) THEN
  cmessage='Error opening file '//TRIM(ADJUSTL(test_file_vsls2)) // ' : '// &
           TRIM(iomessage)
  CALL ereport('TEST_SCENARIO_CL_VSLS',ierr,cmessage)
END IF

! write some header to the file
! start header lines with #
! generate what the format of the header lines should be first
! e.g. "#",X,A10,X,7(A20,X) etc.
! it is necessary to generate these 'on the fly' as there will be a
! variable number of columns dependent on scheme used, or if the output
! is 'full' or just the values UKCA will see (i.e. with lumping).
WRITE(tr_name_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(A20,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_conv_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(E20.10,1X))'
! e.g. 2X,F12.6,7(E20.10,X) etc.
WRITE(tr_float_fmt_str,FMT='(A,I2,A)')                            &
                     '(2X,F12.6,1X,',n_lbc_specs_2D,'(E20.10,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_num_fmt_str,FMT='(A,I2,A)')                              &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(I20,1X))'

! now write header information to the file
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(test_file_vsls2))
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(file_header))
WRITE(file_unit,FMT=tr_name_fmt_str) '#','SPECIES',               &
                    (TRIM(ADJUSTL(lbc_specs_2D(i))),i=1,n_lbc_specs_2D)
WRITE(file_unit,FMT=tr_conv_fmt_str) '#','CONV FAC',              &
                    (lbc_c_species_2D_test(i),i=1,n_lbc_specs_2D)
WRITE(file_unit,FMT=tr_num_fmt_str) '#','YEAR',(i,i=1,n_lbc_specs_2D)

IF (PrintStatus >= Prstatus_Oper) THEN
  WRITE(umMessage,'(A,A)') 'Outputting Cl VSLS LBC data to ',   &
       TRIM(ADJUSTL(test_file_vsls2))
  CALL umPrint(umMessage,src='test_scenario_cl-vsls')
END IF

! call the scenario_rcp routine, as if cycling through many centuries
! do up to 2020 (i.e. 2019+360 days)
DO jyear=1980,2019
  DO jday=1,360
    CALL  ukca_scenario_Cl_VSLSs(n_lbc_specs_2D, lbc_specs_2D,  &
                    jyear, jday,                                &
                    lbc_mmr_2D_test)

    ! output mmr values to our output file
      WRITE(file_unit,FMT=tr_float_fmt_str)                       &
                        (REAL(jyear)+REAL(jday)/360.0),         &
                        (lbc_mmr_2D_test(1,i),i=1,n_lbc_specs_2D)
      WRITE(file_unit,*) '-------'
  END DO
END DO

! all finished now - close the file
CLOSE(file_unit)
CALL release_file_unit(file_unit, handler="fortran")

! deallocate arrays
IF (ALLOCATED(lbc_mmr_2D_test)) DEALLOCATE(lbc_mmr_2D_test)

!!! END OF FILE 2
!!!
!!! REPEAT FOR FILE 3 !!!!

ALLOCATE(lbc_mmr_2D_test(1:rows,1:n_lbc_specs_2D))

! find a file unit to write to
CALL assign_file_unit(TRIM(ADJUSTL(test_file_vsls3)), &
                      file_unit, handler="fortran")

! open file to write test output to
ierr = 0
OPEN(UNIT=file_unit,FILE=TRIM(ADJUSTL(test_file_vsls3)),                    &
     POSITION='rewind',STATUS='unknown',IOSTAT=ierr, IOMSG=iomessage)
IF (ierr /= 0) THEN
  cmessage='Error opening file '//TRIM(ADJUSTL(test_file_vsls3)) // ' : '// &
           TRIM(iomessage)
  CALL ereport('TEST_SCENARIO_CL_VSLS',ierr,cmessage)
END IF

! write some header to the file
! start header lines with #
! generate what the format of the header lines should be first
! e.g. "#",X,A10,X,7(A20,X) etc.
! it is necessary to generate these 'on the fly' as there will be a
! variable number of columns dependent on scheme used, or if the output
! is 'full' or just the values UKCA will see (i.e. with lumping).
WRITE(tr_name_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(A20,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_conv_fmt_str,FMT='(A,I2,A)')                             &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(E20.10,1X))'
! e.g. 2X,F12.6,7(E20.10,X) etc.
WRITE(tr_float_fmt_str,FMT='(A,I2,A)')                            &
                     '(2X,F12.6,1X,',n_lbc_specs_2D,'(E20.10,1X))'
! e.g. "#",X,A10,X,7(E20.10,X) etc.
WRITE(tr_num_fmt_str,FMT='(A,I2,A)')                              &
                     '(A1,1X,A12,1X,',n_lbc_specs_2D,'(I20,1X))'

! now write header information to the file
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(test_file_vsls3))
WRITE(file_unit,FMT='(A1,1X,A)') '#',TRIM(ADJUSTL(file_header))
WRITE(file_unit,FMT=tr_name_fmt_str) '#','SPECIES',               &
                    (TRIM(ADJUSTL(lbc_specs_2D(i))),i=1,n_lbc_specs_2D)
WRITE(file_unit,FMT=tr_conv_fmt_str) '#','CONV FAC',              &
                    (lbc_c_species_2D_test(i),i=1,n_lbc_specs_2D)
WRITE(file_unit,FMT=tr_num_fmt_str) '#','YEAR',(i,i=1,n_lbc_specs_2D)

IF (PrintStatus >= Prstatus_Oper) THEN
  WRITE(umMessage,'(A,A)') 'Outputting Cl VSLS LBC data to ',   &
       TRIM(ADJUSTL(test_file_vsls3))
  CALL umPrint(umMessage,src='test_scenario_cl-vsls')
END IF

! call the scenario_rcp routine, as if cycling through many centuries
! do up to 2020 (i.e. 2019+360 days)
DO jyear=1980,2019
  DO jday=1,360
    CALL  ukca_scenario_Cl_VSLSs(n_lbc_specs_2D, lbc_specs_2D,  &
                    jyear, jday,                                &
                    lbc_mmr_2D_test)

    ! output mmr values to our output file
      WRITE(file_unit,FMT=tr_float_fmt_str)                       &
                        (REAL(jyear)+REAL(jday)/360.0),         &
                        ((lbc_mmr_2D_test(1,i))/lbc_c_species_2D_test(i),i=1,n_lbc_specs_2D)
      WRITE(file_unit,*) '-------'
  END DO
END DO

! all finished now - close the file
CLOSE(file_unit)
CALL release_file_unit(file_unit, handler="fortran")

! deallocate arrays
IF (ALLOCATED(lbc_mmr_2D_test)) DEALLOCATE(lbc_mmr_2D_test)

!!! END OF FILE 3

IF (ALLOCATED(lbc_c_species_2D_test)) DEALLOCATE(lbc_c_species_2D_test)
!IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE test_scenario_cl_vsls

END MODULE ukca_scenario_Cl_VSLSs_mod
