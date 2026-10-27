! copied from acsf22_mod, but updated for CFC12 using 
! SPARC2013 recommendations 

MODULE acsf12_updt_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF12_UPDT_MOD'

CONTAINS
SUBROUTINE acsf12_updt(t,jpwav,wavenm,af12)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of F12 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: af12(jpwav)

! Local variables.

!     Polynomial coefficients.
! DEST - updated to SPARC2013
REAL, PARAMETER :: a0=-654.48602
REAL, PARAMETER :: a1=11.3316594
REAL, PARAMETER :: a2=-0.07413773
REAL, PARAMETER :: a3=0.00021176449
REAL, PARAMETER :: a4=-2.2472927e-7

REAL, PARAMETER :: b0=7.2196167
REAL, PARAMETER :: b1=-0.142178066
REAL, PARAMETER :: b2=0.00104589285
REAL, PARAMETER :: b3=-3.4086704e-6
REAL, PARAMETER :: b4=4.1577462e-9

INTEGER :: jw
REAL :: tc

REAL :: lam
REAL :: lam2
REAL :: lam3
REAL :: lam4
REAL :: arg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF12_UPDT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 295.0 ) tc = 295.0
IF ( tc < 225.0 ) tc = 225.0

!     Wavelength.
DO jw = 45,72

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3
  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +                              &
    (tc-273.0)*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)
  af12(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf12_updt
END MODULE acsf12_updt_mod
