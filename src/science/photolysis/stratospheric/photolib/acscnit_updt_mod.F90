! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acscnit_updt_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate ClONO2 temperature dependent cross sections,

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSCNIT_UPDT_MOD'

CONTAINS
SUBROUTINE acscnit_updt(temp,wavenm,acnita, acnitb)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: temp
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Contains the wavelength intervals
! ACNITA, ACNITB Contain the cross sections times quantum yields
! at model wavelengths.Contains the result on exit.
REAL, INTENT(INOUT) :: acnita(jpwav)
REAL, INTENT(OUT) :: acnitb(jpwav)

! Local variables

INTEGER :: i

REAL :: tc
LOGICAL, SAVE :: first = .TRUE.
!     ACNIT at 296K; JPL2010/TOMCAT
REAL,PARAMETER :: acnitt(jpwav) = [(0.0,i=1,54),                              &
   4.32E-19, 1.98E-18, 2.91E-18, 3.01E-18, 2.87E-18, 2.78E-18, &
   2.78E-18, 2.84E-18, 2.96E-18, 3.10E-18, 3.26E-18, 3.39E-18, &
   3.45E-18, 3.39E-18, 3.24E-18, 2.97E-18, 2.64E-18, 2.27E-18, &
   1.92E-18, 1.59E-18, 1.31E-18, 1.09E-18, 8.97E-19, 7.44E-19, &
   6.07E-19, 5.13E-19, 4.35E-19, 3.70E-19, 3.15E-19, 2.66E-19, &
   2.21E-19, 1.84E-19, 1.50E-19, 1.21E-19, 9.52E-20, 7.33E-20, &
   5.50E-20, 4.01E-20, 2.97E-20, 2.19E-20, 1.60E-20, 1.14E-20, &
   8.31E-21, 6.11E-21, 4.66E-21, 3.66E-21, 3.02E-21, 2.58E-21, &
   2.29E-21, 2.08E-21, 2.00E-21, 1.79E-21, 1.59E-21, 1.41E-21, &
   1.21E-21, 1.06E-21, 9.09E-22, 7.59E-22, 6.38E-22, 5.38E-22, &
   4.44E-22, 3.67E-22, 3.16E-22, 2.31E-22, 1.89E-22, 5.26E-23, &
  (0.0,i=1,83)]
!     Coeffs A1, A2 from JPL2010/TOMCAT
REAL,PARAMETER :: a1(jpwav) = [(0.0,i=1,54),                                  &
   1.73E-05, 7.50E-05, 1.00E-04, 8.82E-05, 3.61E-05,-5.88E-05, &
  -1.95E-04,-3.44E-04,-5.11E-04,-6.59E-04,-7.85E-04,-8.71E-04, &
  -9.03E-04,-8.73E-04,-7.83E-04,-6.38E-04,-4.53E-04,-2.35E-04, &
  -6.97E-06, 2.19E-04, 4.16E-04, 5.64E-04, 6.78E-04, 7.81E-04, &
   9.08E-04, 1.08E-03, 1.26E-03, 1.44E-03, 1.59E-03, 1.74E-03, &
   1.88E-03, 2.02E-03, 2.20E-03, 2.37E-03, 2.55E-03, 2.74E-03, &
   2.95E-03, 3.28E-03, 3.72E-03, 4.12E-03, 4.53E-03, 4.98E-03, &
   5.40E-03, 5.75E-03, 5.92E-03, 5.85E-03, 5.51E-03, 4.92E-03, &
   4.02E-03, 3.27E-03, 2.70E-03, 2.08E-03, 1.33E-03, 7.65E-04, &
   3.53E-04, 2.39E-04, 4.10E-04, 7.77E-04, 1.38E-03, 2.15E-03, &
   3.38E-03, 4.88E-03, 6.70E-03, 8.10E-03, 9.72E-03, 0.00E+00, &
  (0.0,i=1,83)]
REAL,PARAMETER :: a2(jpwav) = [(0.0,i=1,54),                                  &
  -1.16E-06,-5.32E-06,-7.85E-06,-8.21E-06,-7.81E-06,-7.52E-06, &
  -7.46E-06,-7.62E-06,-7.93E-06,-8.30E-06,-8.69E-06,-9.03E-06, &
  -9.25E-06,-9.37E-06,-9.37E-06,-9.27E-06,-9.06E-06,-8.65E-06, &
  -7.98E-06,-7.13E-06,-6.36E-06,-6.00E-06,-6.09E-06,-6.45E-06, &
  -6.57E-06,-6.03E-06,-5.08E-06,-4.20E-06,-3.50E-06,-2.74E-06, &
  -1.87E-06,-9.85E-07, 6.13E-08, 1.13E-06, 2.14E-06, 3.05E-06, &
   3.74E-06, 5.29E-06, 7.78E-06, 9.89E-06, 1.20E-05, 1.49E-05, &
   1.84E-05, 2.27E-05, 2.70E-05, 3.01E-05, 3.11E-05, 2.86E-05, &
   2.07E-05, 1.38E-05, 8.59E-06, 2.01E-06,-7.40E-06,-1.44E-05, &
  -1.91E-05,-2.11E-05,-2.05E-05,-1.87E-05,-1.42E-05,-7.14E-06, &
   4.47E-06, 1.93E-05, 3.87E-05, 5.57E-05, 7.52E-05, 0.00E+00, &
  (0.0,i=1,83)]
REAL, SAVE :: quanty(jpwav) ! Quantum yield for Cl + NO3 channel

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSCNIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (first) THEN
  quanty = MIN(MAX(7.143e-3 * wavenm - 1.6, 0.6), 1.0)
  first = .FALSE.
END IF

tc=temp
tc=MAX(tc, 220.0)
tc=MIN(tc, 296.0)

! Initialise acnitb here as acnita already initialised in fill_spectra
acnitb(:) = 0.0

! ClONO2 cross section, following JPL (2002)
acnita(55:120) = acnitt(55:120)*(1.0 + a1(55:120)*(tc-296.0)                   &
  + a2(55:120)*((tc-296.0)**2.0))
! ClONO2 cross section, times quantum yield for ClO + NO2 channel
acnitb(55:120) = (1.0 - quanty(55:120)) * acnita(55:120)
! ClONO2 cross section, times quantum yield for Cl  + NO3 channel
acnita(55:120) = quanty(55:120) * acnita(55:120)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acscnit_updt
END MODULE acscnit_updt_mod
