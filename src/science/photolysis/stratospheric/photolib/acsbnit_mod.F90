! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
! copied from acscnit_mod but updated for BrONO2 using TOMCAT library
!
MODULE acsbnit_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ACSBNIT_MOD'

CONTAINS
SUBROUTINE acsbnit(temp,wavenm,abrno3)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)    :: temp
REAL, INTENT(IN)    :: wavenm(jpwav)  ! Contains the wavelength intervals
REAL, INTENT(INOUT) :: abrno3(jpwav)

! Local variables

INTEGER :: i

REAL :: tc
LOGICAL, SAVE :: first = .TRUE.
!     ABNIT at 296K; JPL2006/TOMCAT
REAL,PARAMETER :: abnit(jpwav) = (/(0.0,i=1,48),                              &
  0.000E+00,0.000E+00,0.000E+00,7.816E-18,1.410E-17,1.320E-17, &
  1.216E-17,1.106E-17,9.936E-18,8.853E-18,7.754E-18,6.468E-18, &
  5.819E-18,5.155E-18,4.482E-18,3.803E-18,3.385E-18,3.085E-18, &
  2.846E-18,2.673E-18,2.519E-18,2.384E-18,2.253E-18,2.116E-18, &
  1.971E-18,1.800E-18,1.616E-18,1.419E-18,1.228E-18,1.043E-18, &
  8.730E-19,7.301E-19,6.124E-19,5.186E-19,4.450E-19,3.918E-19, &
  3.518E-19,3.228E-19,2.997E-19,2.794E-19,2.612E-19,2.420E-19, &
  2.216E-19,2.005E-19,1.803E-19,1.603E-19,1.418E-19,1.246E-19, &
  1.102E-19,9.861E-20,8.944E-20,8.224E-20,7.640E-20,7.149E-20, &
  6.655E-20,6.188E-20,5.679E-20,5.158E-20,4.648E-20,4.154E-20, &
  3.681E-20,3.235E-20,2.868E-20,2.527E-20,2.244E-20,1.992E-20, &
  1.807E-20,1.654E-20,1.499E-20,1.380E-20,1.295E-20,1.200E-20, &
  1.112E-20,1.033E-20,9.308E-21,8.330E-21,7.409E-21,6.510E-21, &
  5.609E-21,4.627E-21,3.804E-21,2.830E-21,2.385E-21,2.011E-21, &
  7.476E-22,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, &
  (0.0,i=1,65)/)
!     Coeffs A1, A2 from JPL2006/TOMCAT
REAL,PARAMETER :: a1(jpwav) = (/(0.0,i=1,54),                                  &
  0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00,7.300E-04, &
  4.430E-04,2.145E-04,1.510E-04,1.762E-04,3.576E-04,6.442E-04, &
  8.988E-04,1.107E-03,1.215E-03,1.138E-03,1.140E-03,1.302E-03, &
  1.400E-03,1.409E-03,1.512E-03,1.695E-03,1.632E-03,1.602E-03, &
  1.965E-03,2.236E-03,2.282E-03,2.322E-03,2.374E-03,2.391E-03, &
  2.288E-03,2.188E-03,2.103E-03,1.925E-03,1.791E-03,1.480E-03, &
  1.010E-03,6.405E-04,7.318E-04,8.713E-04,6.600E-04,6.475E-04, &
  7.190E-04,7.315E-04,8.430E-04,8.030E-04,8.180E-04,9.865E-04, &
  8.420E-04,1.065E-03,9.690E-04,1.105E-03,1.240E-03,1.395E-03, &
  1.660E-03,1.580E-03,2.110E-03,2.215E-03,1.980E-03,1.885E-03, &
  1.700E-03,1.445E-03,1.290E-03,6.465E-04,8.770E-04,5.640E-04, &
  4.470E-04,3.265E-04,5.670E-04,1.008E-03,1.390E-03,1.315E-03, &
  2.320E-03,3.465E-03,5.240E-03,7.355E-03,1.000E-02,1.400E-02, &
  1.800E-02,0.000E+00,0.000E+00,0.000E+00,0.000E+00,0.000E+00, &
  (0.0,i=1,65)/)
REAL,PARAMETER :: a2(jpwav) = (/(0.0,i=1,54),                                  &
   0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
   -2.605E-05, &
   -2.615E-05,-2.483E-05,-1.967E-05,-1.161E-05,-3.906E-06, &
   4.648E-06, &
   1.161E-05, 1.662E-05, 1.970E-05, 1.962E-05, 2.030E-05, &
   2.228E-05, &
   2.322E-05, 2.281E-05, 2.301E-05, 2.385E-05, 2.016E-05, &
   1.601E-05, &
   1.722E-05, 1.762E-05, 1.556E-05, 1.457E-05, 1.431E-05, &
   1.466E-05, &
   1.385E-05, 1.361E-05, 1.371E-05, 1.278E-05, 1.233E-05, &
   9.229E-06, &
   4.110E-06,-1.894E-07, 1.894E-07, 1.711E-06,-2.280E-06, &
   -1.500E-06, &
   -1.560E-06,-2.105E-07, 1.410E-06, 1.275E-06,1.700E-06, &
   3.430E-06, &
   1.980E-06, 3.920E-06, 3.620E-06, 4.865E-06, 6.630E-06, &
   7.495E-06, &
   9.150E-06, 7.030E-06, 9.670E-06, 9.440E-06, 5.340E-06, &
   3.135E-06, &
   3.160E-07,-2.485E-06,-3.810E-07,-4.370E-06, 8.800E-07, &
   7.965E-08, &
   1.200E-06,-1.940E-06, 3.180E-06, 6.340E-06, 9.160E-06, &
   2.165E-06, &
   7.950E-06, 1.615E-05, 2.540E-05, 5.620E-05, 6.520E-05, &
   8.075E-05, &
   1.190E-04, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, &
   0.000E+00, &
  (0.0,i=1,65)/)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ACSBNIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

tc=temp
tc=MAX(tc, 200.0)
tc=MIN(tc, 300.0)

abrno3(60:133) = abnit(60:133)*(1.0 + a1(60:133)*(tc-296.0)                   &
  + a2(60:133)*((tc-296.0)**2.0))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsbnit
END MODULE acsbnit_mod
