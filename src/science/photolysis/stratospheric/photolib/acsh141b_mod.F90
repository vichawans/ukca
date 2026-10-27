! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!
! UM-UKCA DEST
! EMB: module based on acsf22_mod.f but edited for H141b
!
MODULE acsh141b_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSH141B_MOD'

CONTAINS
SUBROUTINE acsh141b(t,jpwav,wavenm,ah141b)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of H141b in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(INOUT) :: ah141b(jpwav)

! Local variables.

!     Polynomial coefficients.
REAL, PARAMETER :: a0=-682.913042 
REAL, PARAMETER :: a1=12.122290
REAL, PARAMETER :: a2=-8.187699e-2
REAL, PARAMETER :: a3=2.437244e-4
REAL, PARAMETER :: a4=-2.719103e-7

REAL, PARAMETER :: b0=4.074747
REAL, PARAMETER :: b1=-8.053899e-2
REAL, PARAMETER :: b2=5.946552e-4
REAL, PARAMETER :: b3=-1.945048e-6
REAL, PARAMETER :: b4=2.380143e-9

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSH141B'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 300.0 ) tc = 300.0
IF ( tc < 210.0 ) tc = 210.0

! H141b: T-dep over 172-240nm
!     Wavelength.
DO jw = 44,76

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3
  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +   &
    (tc-273.0)*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)
  ah141b(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsh141b
END MODULE acsh141b_mod
