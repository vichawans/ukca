! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsf22_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     A subroutine which calculates the F22 (CHF2Cl) absorption
!     cross section based on P. C. Simon et al. (1988).

! Method:
!     Journal of atmospheric chemistry, Vol. 7, pp. 107-135, 1988.

!     This is done via a polynomial expression of the form;
!     log10(sigma)=A(lamda) + T B(lamda)

!     The expression is valid for the wavelength range; 174-204 nm
!                            and the temperature range; 210-300 K.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSF22_MOD'

CONTAINS
SUBROUTINE acsf22(t,jpwav,wavenm,af22)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: t              ! Temperature in kelvin
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Wavelength of each interval in nm
! Absorption cross-section of F22 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(IN OUT) :: af22(jpwav)

! Local variables.

!     Polynomial coefficients.
! DEST - updated to SPARC2013
REAL, PARAMETER :: a0=-369.398
REAL, PARAMETER :: a1=7.09844
REAL, PARAMETER :: a2=-0.0527699
REAL, PARAMETER :: a3=0.0001715165
REAL, PARAMETER :: a4=-2.0822e-7

REAL, PARAMETER :: b0=-6.5412
REAL, PARAMETER :: b1=0.138288
REAL, PARAMETER :: b2=-0.00109207
REAL, PARAMETER :: b3=3.8164e-6
REAL, PARAMETER :: b4=-4.9772e-9

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSF22'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check that temperature is in range.
tc = t
IF ( tc > 300.0 ) tc = 300.0
IF ( tc < 210.0 ) tc = 210.0

!     Wavelength.
DO jw = 45 , 61

  !        Wavelength in nm.
  lam = wavenm(jw)

  lam2 = lam*lam
  lam3 = lam*lam2
  lam4 = lam*lam3
  arg = a0 + a1*lam + a2*lam2 + a3*lam3 + a4*lam4 +                              &
    (tc-273.0)*(b0+b1*lam+b2*lam2+b3*lam3+b4*lam4)
  af22(jw) = 10.0**arg

END DO
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsf22
END MODULE acsf22_mod
