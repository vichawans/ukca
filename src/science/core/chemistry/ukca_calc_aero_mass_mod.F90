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
!   - total sulfate (cp_su) over all soluble modes
!   - total OM1 (cp_oc) over all soluble modes
!   - total OM2 (cp_so) over all soluble modes
!   - total aerosol mass over all soluble components and modes
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
! Code Description:
!   Language:  FORTRAN 90
!
! ######################################################################

MODULE ukca_calc_aero_mass_mod
IMPLICIT NONE
PRIVATE

CHARACTER(LEN=*), PARAMETER :: ModuleName='UKCA_CALC_AERO_MASS_MOD'

PUBLIC :: ukca_calc_aero_mass

CONTAINS

SUBROUTINE ukca_calc_aero_mass(tracer, n_chemistry_tracers,            &
                                row_length, rows, model_levels,         &
                                so4_total_mmr, om1_total_mmr,            &
                                om2_total_mmr, aero_total_mmr)

USE ukca_mode_setup,               ONLY: cp_su, cp_oc, cp_so,                  &
                                        glomap_variables_type,                 & 
                                        mode_nuc_sol, mode_cor_sol
USE ukca_config_specification_mod, ONLY: glomap_variables
USE ukca_mode_tracer_maps_mod,     ONLY: mmr_index
USE parkind1,                      ONLY: jprb, jpim
USE yomhook,                       ONLY: lhook, dr_hook

IMPLICIT NONE


! tracer MMR (row_length,rows,model_levels,ntracers)
REAL, INTENT(IN)    :: tracer(:,:,:,:)

! variables from outside
INTEGER, INTENT(IN) :: n_chemistry_tracers
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels

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

INTEGER, POINTER :: ncp
LOGICAL, POINTER :: component(:,:)
LOGICAL, POINTER :: mode(:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_CALC_AERO_MASS_TOTALS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

component => glomap_variables%component
mode      => glomap_variables%mode
ncp       => glomap_variables%ncp

! WRITE (*,*) '!--VS component=', component
! WRITE (*,*) '!--VS mode=', mode
! WRITE (*,*) '!--VS ncp=', ncp

so4_total_mmr(:,:,:) = 0.0
om1_total_mmr(:,:,:) = 0.0
om2_total_mmr(:,:,:) = 0.0
aero_total_mmr(:,:,:) = 0.0

ifirst = n_chemistry_tracers + 1
ntr_all = SIZE(tracer,4)
WRITE (*,*) '!--VS ntr_all=', ntr_all
WRITE (*,*) '!--VS n_chemistry_tracers=', n_chemistry_tracers


DO imode = mode_nuc_sol, mode_cor_sol ! only take the soluble modes
  IF (.NOT. mode(imode)) CYCLE
  ! Sulfate total across all soluble modes.
  IF (cp_su <= ncp) THEN
    IF (component(imode, cp_su) .AND. mmr_index(imode, cp_su) > 0) THEN
      iaer = mmr_index(imode, cp_su) - ifirst + 1
      WRITE (*,*) '!--VS sulfate iaer=', iaer
      WRITE (*,*) '!--VS sulfate imode=', imode
      WRITE (*,*) '!--VS sulfate mmr_index(imode, cp_su)=', mmr_index(imode, cp_su)
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        WRITE (*,*) '!--VS adding sulfate to mixing ratio'
        WRITE (*,*) '!--VS sum(tracer(:,:,:,mmr_index(imode, cp_su)))=',       &
                    sum(tracer(:,:,:,mmr_index(imode, cp_su)))
        so4_total_mmr(:,:,:) =                                                 & 
          so4_total_mmr(:,:,:) + tracer(:,:,:,mmr_index(imode, cp_su))
      END IF
    END IF
  END IF

  ! OM1 (cp_oc) total across all soluble modes.
  IF (cp_oc <= ncp) THEN
    IF (component(imode, cp_oc) .AND. mmr_index(imode, cp_oc) > 0) THEN
      iaer = mmr_index(imode, cp_oc) - ifirst + 1
      WRITE (*,*) '!--VS oc iaer=', iaer
      WRITE (*,*) '!--VS oc imode=', imode
      WRITE (*,*) '!--VS oc mmr_index(imode, cp_oc)=', mmr_index(imode, cp_oc)
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        WRITE (*,*) '!--VS adding om1 to mixing ratio'
        WRITE (*,*) '!--VS um(tracer(:,:,:,mmr_index(imode, cp_oc)))=',        &
                    sum(tracer(:,:,:,mmr_index(imode, cp_oc)))
        om1_total_mmr(:,:,:) =                                                 &
          om1_total_mmr(:,:,:) + tracer(:,:,:,mmr_index(imode, cp_oc))
      END IF
    END IF
  END IF

  ! OM2 (cp_so) may be unavailable in reduced-component setups.
  IF (cp_so <= ncp) THEN

    IF (component(imode, cp_so) .AND. mmr_index(imode, cp_so) > 0) THEN
      iaer = mmr_index(imode, cp_so) - ifirst + 1
      WRITE (*,*) '!--VS so iaer=', iaer
      WRITE (*,*) '!--VS so imode=', imode
      WRITE (*,*) '!--VS so mmr_index(imode, cp_so)=', mmr_index(imode, cp_so)
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        WRITE (*,*) '!--VS adding om2 to mixing ratio'
        WRITE (*,*) '!--VS ssum(tracer(:,:,:,mmr_index(imode, cp_so)))=',      &
                    sum(tracer(:,:,:,mmr_index(imode, cp_so)))
        om2_total_mmr(:,:,:) =                                                 &
          om2_total_mmr(:,:,:) + tracer(:,:,:,mmr_index(imode, cp_so))
      END IF
    END IF
  END IF

  ! Sum all soluble aerosol components for total aerosol mass.
  DO icp = 1, ncp
    IF (component(imode, icp) .AND. mmr_index(imode, icp) > 0) THEN
      iaer = mmr_index(imode, icp) - ifirst + 1
      WRITE (*,*) '!--VS aer iaer=', iaer
      WRITE (*,*) '!--VS aer imode=', imode
      WRITE (*,*) '!--VS aer mmr_index(imode, icp)=', mmr_index(imode, icp)
      IF (iaer >= 1 .AND. iaer <= ntr_all) THEN
        WRITE (*,*) '!--VS adding aer to mixing ratio'
        WRITE (*,*) '!--VS sum(tracer(:,:,:,mmr_index(imode, icp)))=',         &
                    sum(tracer(:,:,:,mmr_index(imode, icp)))
        aero_total_mmr(:,:,:) =                                                &
          aero_total_mmr(:,:,:) + tracer(:,:,:,mmr_index(imode, icp))
      END IF
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ukca_calc_aero_mass

END MODULE ukca_calc_aero_mass_mod