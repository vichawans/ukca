! *****************************COPYRIGHT*******************************
! copied from acsf22_mod.F90 but adapted for halon 1301 using sparc 2013 parametrization 

MODULE acsf1301_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF1301_MOD'

CONTAINS
SUBROUTINE acsf1301(t,jpwav,wavenm,af13b1)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of F22 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: af13b1(jpwav)

! Local variables.

!     Polynomial coefficients.
! DEST - updated to SPARC2013
REAL, PARAMETER :: a0=-46.70542d0
REAL, PARAMETER :: a1=-1.55047d-1
REAL, PARAMETER :: a2=-1.020187d-3
REAL, PARAMETER :: a3=2.246169d-5
REAL, PARAMETER :: a4=-1.300982d-7

REAL, PARAMETER :: b0=1.694026d-4
REAL, PARAMETER :: b1=8.723247d-6
REAL, PARAMETER :: b2=5.953165d-9
REAL, PARAMETER :: b3=-3.872168d-9
REAL, PARAMETER :: b4=-1.803325d-11

INTEGER :: k
REAL :: tc

REAL :: dlam
REAL :: dlam2
REAL :: dlam3
REAL :: dlam4
REAL :: arg_a
REAL :: arg_b

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF1301'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 296.0 ) tc = 296.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength.
DO k = 1,jpwav
  IF ( (wavenm(k) .ge. 214.d0) .and. (wavenm(k) .le. 285.d0) ) THEN
    dlam = wavenm(k)-242.2466d0
    dlam2 = dlam*dlam
    dlam3 = dlam*dlam2
    dlam4 = dlam*dlam3
    arg_a = (a0 + a1*dlam + a2*dlam2 + a3*dlam3 + a4*dlam4)
    arg_b = (b0 + b1*dlam + b2*dlam2 + b3*dlam3 + b4*dlam4)
    af13b1(k) = DEXP(arg_a + (296.d0-tc)*arg_b*arg_a)
  ELSE
   af13b1(k) = af13b1(k)
  ENDIF

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf1301
END MODULE acsf1301_mod
