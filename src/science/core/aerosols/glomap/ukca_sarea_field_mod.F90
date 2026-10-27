! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Field to hold partial surface area to be used in ukca_hetero_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: UKCA
!
! Code Description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.

MODULE ukca_sarea_field_mod

IMPLICIT NONE
! surface area of sulfate aerosols [cm^2/cm^3]
REAL, ALLOCATABLE, SAVE :: sarea_su(:,:,:)
! surface area of black carbon aerosols [cm^2/cm^3]
REAL, ALLOCATABLE, SAVE :: sarea_bc(:,:,:)
! surface area of organic matter 1 aerosols [cm^2/cm^3]
REAL, ALLOCATABLE, SAVE :: sarea_om(:,:,:)

END MODULE ukca_sarea_field_mod
