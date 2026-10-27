! *****************************COPYRIGHT*******************************
! copied from acsf22_mod.F90 but adapted for halon 1211 using sparc 2013 parametrization 

MODULE acsf1211_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF1211_MOD'

CONTAINS
SUBROUTINE acsf1211(t,jpwav,wavenm,af12b1)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of F22 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: af12b1(jpwav)

! Local variables.

!     Polynomial coefficients.
! DEST - updated to SPARC2013
REAL, PARAMETER :: a0=-45.4087d0
REAL, PARAMETER :: a1=-1.304811d-1
REAL, PARAMETER :: a2=-6.995443d-4
REAL, PARAMETER :: a3=6.159709d-6
REAL, PARAMETER :: a4=-9.384074d-9

REAL, PARAMETER :: b0=1.528905d-4
REAL, PARAMETER :: b1=6.024833d-6
REAL, PARAMETER :: b2=1.030995d-7
REAL, PARAMETER :: b3=-6.387931d-11
REAL, PARAMETER :: b4=-3.718503d-11

REAL, PARAMETER :: aa0=-48.3578d0
REAL, PARAMETER :: aa1=-0.1547325d0
REAL, PARAMETER :: aa2=-4.966942d-4
REAL, PARAMETER :: aa3=1.56338d-6
REAL, PARAMETER :: aa4=3.664034d-8

REAL, PARAMETER :: bb0=0.0002989d0
REAL, PARAMETER :: bb1=8.5306d-6
REAL, PARAMETER :: bb2=4.26d-8
REAL, PARAMETER :: bb3=-1.84d-9
REAL, PARAMETER :: bb4=1.284d-11

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF1211'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 296.0 ) tc = 296.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength.
DO k = 1,jpwav
  IF ( (wavenm(k) .ge. 220.d0) .and. (wavenm(k) .lt. 260.d0) ) THEN
    dlam = wavenm(k)-259.8989d0
    dlam2 = dlam*dlam
    dlam3 = dlam*dlam2
    dlam4 = dlam*dlam3
    arg_a = (a0 + a1*dlam + a2*dlam2 + a3*dlam3 + a4*dlam4)
    arg_b = (b0 + b1*dlam + b2*dlam2 + b3*dlam3 + b4*dlam4)
    af12b1(k) = DEXP(arg_a + (296.d0-tc)*arg_b*arg_a)
  ELSE IF ( (wavenm(k) .ge. 260.d0) .and. (wavenm(k) .le. 340.d0) ) THEN
    dlam = wavenm(k)-280.376d0
    dlam2 = dlam*dlam
    dlam3 = dlam*dlam2
    dlam4 = dlam*dlam3
    arg_a = (aa0 + aa1*dlam + aa2*dlam2 + aa3*dlam3 + aa4*dlam4)
    arg_b = (bb0 + bb1*dlam + bb2*dlam2 + bb3*dlam3 + bb4*dlam4)
    af12b1(k) = DEXP(arg_a + (296.d0-tc)*arg_b*arg_a)
  ELSE
   af12b1(k) = af12b1(k)
  ENDIF

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf1211
END MODULE acsf1211_mod
