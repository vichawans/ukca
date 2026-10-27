! copied from acsf22_mod, updated for CFC114
!

MODULE acsf114_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF114_MOD'

CONTAINS
SUBROUTINE acsf114(t,jpwav,wavenm,af114)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: af114(jpwav)

! Local variables.

!     Polynomial coefficients.
! SPARC2013
REAL, PARAMETER :: a0=-184.86
REAL, PARAMETER :: a1=2.98117
REAL, PARAMETER :: a2=-0.0190482
REAL, PARAMETER :: a3=5.152e-5
REAL, PARAMETER :: a4=-5.1086e-8

REAL, PARAMETER :: b0=-1.9446
REAL, PARAMETER :: b1=0.0437746
REAL, PARAMETER :: b2=-0.000365055
REAL, PARAMETER :: b3=1.33636e-6
REAL, PARAMETER :: b4=-1.81075e-9

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF114'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 300.0 ) tc = 300.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength.
! 172-220nm
DO jw = 44,68

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3
  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +                              &
    (tc-273.0)*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)
  af114(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf114
END MODULE acsf114_mod
