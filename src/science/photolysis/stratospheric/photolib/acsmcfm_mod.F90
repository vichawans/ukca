! copied from acsf22_mod, updated for MeCCl3
!

MODULE acsmcfm_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSMCFM_MOD'

CONTAINS
SUBROUTINE acsmcfm(t,jpwav,wavenm,amcfm)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: amcfm(jpwav)

! Local variables.

!     Polynomial coefficients.
! SPARC2013
REAL, PARAMETER :: a0=164.15
REAL, PARAMETER :: a1=-3.866
REAL, PARAMETER :: a2=0.030461
REAL, PARAMETER :: a3=-0.000104583
REAL, PARAMETER :: a4=1.30734e-7

REAL, PARAMETER :: b0=0.65591
REAL, PARAMETER :: b1=-0.0131953
REAL, PARAMETER :: b2=0.000101047
REAL, PARAMETER :: b3=-3.4977e-7
REAL, PARAMETER :: b4=4.6223e-10

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSMCFM'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 300.0 ) tc = 300.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength.
! 182-240nm
DO jw = 50,76

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3
  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +                              &
    (tc-273.0)*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)
  amcfm(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsmcfm
END MODULE acsmcfm_mod
