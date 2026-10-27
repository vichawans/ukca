! copied from acsf22_mod, but updated for CFC11 using 
! SPARC2013 recommendations 

MODULE acsf11_updt_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF11_UPDT_MOD'

CONTAINS
SUBROUTINE acsf11_updt(t,jpwav,wavenm,af11)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of F11 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: af11(jpwav)

! Local variables.

!     Polynomial coefficients.
! DEST - updated to SPARC2013
REAL, PARAMETER :: a0=-82.55
REAL, PARAMETER :: a1=0.75448
REAL, PARAMETER :: a2=-0.001749
REAL, PARAMETER :: a3=-5.493e-6
REAL, PARAMETER :: a4=1.7089e-8

REAL, PARAMETER :: b0=-5.7463
REAL, PARAMETER :: b1=0.115997
REAL, PARAMETER :: b2=-0.000874014
REAL, PARAMETER :: b3=2.91147e-6
REAL, PARAMETER :: b4=-3.6148e-9

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF11_UPDT'


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
  af11(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf11_updt
END MODULE acsf11_updt_mod
