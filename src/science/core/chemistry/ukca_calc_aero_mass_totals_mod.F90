! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!  Module containing subroutine to calculate aerosol mass-mixing-ratio totals
!  from MODE aerosol tracers:
!   - total sulfate (cp_su) over all active modes
!   - total OM1 (cp_oc) over all active modes
!   - total OM2 (cp_so) over all active modes
!   - total aerosol mass over all active components and modes
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!
! ######################################################################

MODULE ukca_calc_aero_mass_totals_mod

IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_CALC_AERO_MASS_TOTALS_MOD'

PUBLIC :: ukca_calc_aero_mass_totals

CONTAINS

SUBROUTINE ukca_calc_aero_mass_totals(tracer, n_chemistry_tracers,            &
                                      glomap_variables_local,                 &
                                      row_length, rows, model_levels,         &
                                      so4_total_mmr, om1_total_mmr,            &
                                      om2_total_mmr, aero_total_mmr)

USE ukca_mode_setup,               ONLY: nmodes, cp_su, cp_oc, cp_so,         &
                                        glomap_variables_type
USE ukca_mode_tracer_maps_mod,     ONLY: mmr_index
USE parkind1,                      ONLY: jprb, jpim
USE yomhook,                       ONLY: lhook, dr_hook

IMPLICIT NONE

! variables from outside
INTEGER, INTENT(IN) :: n_chemistry_tracers
TYPE(glomap_variables_type), TARGET, INTENT(IN) :: glomap_variables_local
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

! tracer MMR (row_length,rows,model_levels,ntracers)
REAL, INTENT(IN)    :: tracer(:,:,:,:)

! output totals in kg/kg(dry air)
REAL, INTENT(OUT) :: so4_total_mmr(row_length, rows, model_levels)
REAL, INTENT(OUT) :: om1_total_mmr(row_length, rows, model_levels)
REAL, INTENT(OUT) :: om2_total_mmr(row_length, rows, model_levels)
REAL, INTENT(OUT) :: aero_total_mmr(row_length, rows, model_levels)

! local variables
INTEGER :: imode
INTEGER :: icp
INTEGER :: iaer
INTEGER :: ifirst
INTEGER :: ntr_all
INTEGER :: ncp

LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode(:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_AERO_MASS_TOTALS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

component => glomap_variables_local%component
mode      => glomap_variables_local%mode
ncp       => glomap_variables_local%ncp

so4_total_mmr(:,:,:) = 0.0
om1_total_mmr(:,:,:) = 0.0
om2_total_mmr(:,:,:) = 0.0
aero_total_mmr(:,:,:) = 0.0

IF (.NOT. ALLOCATED(mmr_index)) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

ifirst = n_chemistry_tracers + 1
ntr_all = SIZE(tracer,4)

DO imode = 1, nmodes
  IF (.NOT. mode(imode)) CYCLE

  ! Sulfate total across all active modes.
  IF (cp_su <= ncp) THEN
    IF (component(imode, cp_su) .AND. mmr_index(imode, cp_su) > 0) THEN
      iaer = n_chemistry_tracers + mmr_index(imode, cp_su) - ifirst + 1
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        so4_total_mmr(:,:,:) = so4_total_mmr(:,:,:) + tracer(:,:,:,iaer)
      END IF
    END IF
  END IF

  ! OM1 (cp_oc) total across all active modes.
  IF (cp_oc <= ncp) THEN
    IF (component(imode, cp_oc) .AND. mmr_index(imode, cp_oc) > 0) THEN
      iaer = n_chemistry_tracers + mmr_index(imode, cp_oc) - ifirst + 1
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        om1_total_mmr(:,:,:) = om1_total_mmr(:,:,:) + tracer(:,:,:,iaer)
      END IF
    END IF
  END IF

  ! OM2 (cp_so) may be unavailable in reduced-component setups.
  IF (cp_so <= ncp) THEN
    IF (component(imode, cp_so) .AND. mmr_index(imode, cp_so) > 0) THEN
      iaer = n_chemistry_tracers + mmr_index(imode, cp_so) - ifirst + 1
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        om2_total_mmr(:,:,:) = om2_total_mmr(:,:,:) + tracer(:,:,:,iaer)
      END IF
    END IF
  END IF

  ! Sum all active aerosol components for total aerosol mass.
  DO icp = 1, ncp
    IF (component(imode, icp) .AND. mmr_index(imode, icp) > 0) THEN
      iaer = n_chemistry_tracers + mmr_index(imode, icp) - ifirst + 1
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        aero_total_mmr(:,:,:) = aero_total_mmr(:,:,:) + tracer(:,:,:,iaer)
      END IF
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_calc_aero_mass_totals

END MODULE ukca_calc_aero_mass_totals_mod