! copied from acsf22_mod, updated for CFC113
!

MODULE acsf113_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF113_MOD'

CONTAINS
SUBROUTINE acsf113(t,jpwav,wavenm,af113)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: af113(jpwav)

! Local variables.

!     Polynomial coefficients.
! SPARC2013
REAL, PARAMETER :: a0=-1112.111832
REAL, PARAMETER :: a1=20.47564047
REAL, PARAMETER :: a2=-0.1426332745
REAL, PARAMETER :: a3=0.0004393747238
REAL, PARAMETER :: a4=-5.072508872e-7

REAL, PARAMETER :: b0=11.62056201
REAL, PARAMETER :: b1=-0.2223995208
REAL, PARAMETER :: b2=0.001590706707
REAL, PARAMETER :: b3=-5.040602321e-6
REAL, PARAMETER :: b4=5.973112715e-9

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF113'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 300.0 ) tc = 300.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength.
! 182-230nm
DO jw = 50,72

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3
  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +                              &
    (tc-273.0)*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)
  af113(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf113
END MODULE acsf113_mod
