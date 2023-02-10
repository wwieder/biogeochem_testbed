! ==============================================================================
!
! File:    mimics_variable_CN.f90
!
! Purpose: defines/allocates variables for MIMICS model
!   MODULE mimicsdimension
!   MODULE mimicsparam
!   MODULE mimicsvariable 
!     with subroutine alloc_mimicsvariable
!
! Contact: Melannie Hartman
!          melannie@ucar.edu
!
! History: 
!   12/22/2014 - Created 
!   1/5/2015   - Compiled
!   6/10/2019  - Added MIMICS-CN variables
! ==============================================================================

MODULE mimicsdimension

  IMPLICIT NONE

  INTEGER, PARAMETER :: nmimicslit =  2    ! Number of litter pools (METBC, STRUC)
  INTEGER, PARAMETER :: nmimicsmcrb = 2    ! Number of microbial pools (RSEL, KSEL)
  INTEGER, PARAMETER :: nmimicssom =  3    ! Number of SOM pools (AVAL, CHEM, PHYS)
  INTEGER, PARAMETER :: nRK   = 6       ! r1, r2, r3, k1, k2, k3

END MODULE mimicsdimension
!-------------------------------------------------------------------------------- 

MODULE mimicsparam
  USE mimicsdimension

  !Litter Pools
  INTEGER, PARAMETER :: METBC   = 1     ! Metabolic
  INTEGER, PARAMETER :: STRUC   = 2     ! Structural

  !Microbial Pools
  INTEGER, PARAMETER :: RSEL    = 1     ! r-selected
  INTEGER, PARAMETER :: KSEL    = 2     ! K-selected

  !Soil Organic Matter Pools
  INTEGER, PARAMETER :: AVAL    = 1     ! Available
  INTEGER, PARAMETER :: CHEM    = 2     ! Chemically Protected
  INTEGER, PARAMETER :: PHYS    = 3     ! Physically Protected

  !Decomposition of litter and active SOM by microbial pools:
  INTEGER, PARAMETER :: R1 = 1  ! LITm to MICr
  INTEGER, PARAMETER :: R2 = 2  ! LITs to MICr
  INTEGER, PARAMETER :: R3 = 3  ! SOMa to MICr
  INTEGER, PARAMETER :: K1 = 4  ! LITm to MICk
  INTEGER, PARAMETER :: K2 = 5  ! LITs to MICk
  INTEGER, PARAMETER :: K3 = 6  ! SOMa to MICk

END MODULE mimicsparam
!-------------------------------------------------------------------------------- 

MODULE mimicsvariable
  USE casadimension
  USE mimicsdimension

!-------------------------
  TYPE mimics_biome

    REAL(r_2), POINTER :: tauModDenom, tauMod_MIN, tauMod_MAX

    REAL(r_2), POINTER :: densDep, &
                          fracDINavailMIC, &
                          cnModNum, CNr, CNk,  &
                          desorpQ10, desorpTref

    REAL(r_2), DIMENSION(:),    POINTER :: fmet,    &
                                           ligninNratioAvg, &
                                           CN_r,    &
                                           CN_k,    &
                                           Vslope,  &
                                           Kslope,  &
                                           Vmod,    &
                                           Vint,    &
                                           av,      &
                                           Kint,    &
                                           ak,    &
                                           Fi,      &
                                           MGE,     &
                                           NUE,     &
                                           KO,      &
                                           Pscalar, &
                                           desorp,  &
                                           tauR,    &
                                           tau_r,   &
                                           tauK,    &
                                           tau_k,   &
                                           tauMod,  &
                                           depth,   &
                                           fPHYS_r, &
                                           fPHYS_K, &
                                           fCHEM_r, &
                                           fCHEM_K, &
                                           fSOM_p, &
                                           fmet_p, &
                                           phys_scalar


    REAL(r_2), DIMENSION(:,:),  POINTER :: ligninNratio, &
                                           Kmod, &
                                           Vmax, &
                                           Km,   &
                                           fPHYS,  &
                                           fCHEM,  &
                                           fAVAL  
    INTEGER:: fWFunction
 
  END TYPE mimics_biome

!-------------------------
  TYPE mimics_flux
    REAL(r_2), DIMENSION(:,:),POINTER :: ClitInput
    REAL(r_2), DIMENSION(:)  ,POINTER :: Chresp
    REAL(r_2), DIMENSION(:)  ,POINTER :: CSOMpInput
    REAL(r_2), DIMENSION(:)  ,POINTER :: Overflow_r
    REAL(r_2), DIMENSION(:)  ,POINTER :: Overflow_k
    REAL(r_2), DIMENSION(:,:),POINTER :: NlitInput
  END TYPE mimics_flux

!-------------------------
  TYPE mimics_pool
    REAL(r_2), DIMENSION(:),POINTER :: LITm, &
                                       LITs, &
                                       MICr, &
                                       MICk, &
                                       SOMa, &
                                       SOMc, &
                                       SOMp, &
                                       thetaLiq, &
                                       thetaFrzn, &
                                       fT,   &
                                       fW

    REAL(r_2), DIMENSION(:),POINTER :: LITmN, &
                                       LITsN, &
                                       MICrN, &
                                       MICkN, &
                                       SOMaN, &
                                       SOMcN, &
                                       SOMpN, &
                                       DIN
  END TYPE mimics_pool

!-------------------------
  TYPE mimics_annual_pool
      
    REAL(r_2), DIMENSION(:,:), POINTER :: ClitterAn,    &
                                          CsoilAn,      &
                                          CmicrobeAn

    REAL(r_2), DIMENSION(:,:), POINTER :: NlitterAn,    &
                                          NsoilAn,      &
                                          NmicrobeAn

    REAL(r_2), DIMENSION(:), POINTER :: thetaLiqAn,  &
                                        thetaFrznAn, &
                                        fTAn,        &
                                        fWAn,        &
                                        DINAn

  END TYPE mimics_annual_pool

!-------------------------
  TYPE mimics_annual_flux
      
    REAL(r_2), DIMENSION(:,:), POINTER :: ClitInputAn
    REAL(r_2), DIMENSION(:),   POINTER :: ChrespAn
    REAL(r_2), DIMENSION(:),   POINTER :: CSOMpInputAn
    REAL(r_2), DIMENSION(:),   POINTER :: Overflow_rAn
    REAL(r_2), DIMENSION(:),   POINTER :: Overflow_kAn
    REAL(r_2), DIMENSION(:,:), POINTER :: NlitInputAn

  END TYPE mimics_annual_flux

!-------------------------
  TYPE (mimics_biome)              :: mimicsbiome
  TYPE (mimics_flux)               :: mimicsflux
  TYPE (mimics_pool)               :: mimicspool
  TYPE (mimics_annual_pool)        :: mimicspoolAn
  TYPE (mimics_annual_flux)        :: mimicsfluxAn

  character(len=100) :: filename_mimicsbiome  !! file for biome-specific MIMICS parameters
  character(len=100) :: sPtFileNameMIMICS
  integer :: iptToSave_mimics

CONTAINS

!-------------------------
SUBROUTINE alloc_mimicsvariable(mp,mvtype,mplant)
  USE mimicsdimension
  implicit none
  integer, INTENT(IN) :: mp     ! number of grid points
  integer, INTENT(IN) :: mvtype ! number of vegetation (biome) types
  integer, INTENT(IN) :: mplant ! number of plant pools
  integer :: arraysize
  arraysize = mp

  ALLOCATE(mimicsbiome%tauModDenom, &
           mimicsbiome%tauMod_MIN, &
           mimicsbiome%tauMod_MAX, &
           mimicsbiome%densDep, &
           mimicsbiome%CNr, &
           mimicsbiome%CNk, &
           mimicsbiome%cnModNum, &
           mimicsbiome%fracDINavailMIC, &
           mimicsbiome%desorpQ10, &
           mimicsbiome%desorpTref)

  ALLOCATE(mimicsbiome%Vint(nRK), &
           mimicsbiome%av(nRK),   &
           mimicsbiome%Kint(nRK), & 
           mimicsbiome%ak(nRK),   &
           mimicsbiome%Vslope(nRK),    &
           mimicsbiome%Kslope(nRK),    &
           mimicsbiome%Vmod(nRK),      &
           mimicsbiome%Kmod(arraysize,nRK),      &
           mimicsbiome%Vmax(arraysize,nRK),      &
           mimicsbiome%Km(arraysize,nRK),        &
           mimicsbiome%Fi(nmimicslit), &
           mimicsbiome%MGE(4),         &
           mimicsbiome%NUE(4),         &
           mimicsbiome%KO(2),          &
           mimicsbiome%depth(mvtype),      &
           mimicsbiome%Pscalar(arraysize), &
           mimicsbiome%desorp(arraysize),  &
           mimicsbiome%fmet(arraysize),    &
           mimicsbiome%ligninNratioAvg(arraysize), &
           mimicsbiome%CN_r(arraysize),    &
           mimicsbiome%CN_k(arraysize),    &
           mimicsbiome%tauR(arraysize),    &      
           mimicsbiome%tauK(arraysize),    &      
           mimicsbiome%tauMod(arraysize),  & 
           mimicsbiome%fPHYS(arraysize,2), &     
           mimicsbiome%fCHEM(arraysize,2), &     
           mimicsbiome%fAVAL(arraysize,2), &
           mimicsbiome%tau_r(2),   &
           mimicsbiome%tau_k(2),   &
           mimicsbiome%fPHYS_r(2), &
           mimicsbiome%fPHYS_K(2), &
           mimicsbiome%fCHEM_r(3), &
           mimicsbiome%fCHEM_K(3), &
           mimicsbiome%fSOM_p(2), &
           mimicsbiome%phys_scalar(2), &
           mimicsbiome%fmet_p(3), &
           mimicsbiome%ligninNratio(arraysize,mplant))

  ALLOCATE(mimicsflux%ClitInput(arraysize,nmimicslit))
  ALLOCATE(mimicsflux%Chresp(arraysize))
  ALLOCATE(mimicsflux%CSOMpInput(arraysize))
  ALLOCATE(mimicsflux%Overflow_r(arraysize))
  ALLOCATE(mimicsflux%Overflow_k(arraysize))
  ALLOCATE(mimicsflux%NlitInput(arraysize,nmimicslit))

  ALLOCATE(mimicspool%LITm(arraysize),        &
           mimicspool%LITs(arraysize),        &
           mimicspool%MICr(arraysize),        &
           mimicspool%MICk(arraysize),        &
           mimicspool%SOMa(arraysize),        &
           mimicspool%SOMc(arraysize),        &
           mimicspool%SOMp(arraysize),        &
           mimicspool%thetaLiq(arraysize),    &
           mimicspool%thetaFrzn(arraysize),   &
           mimicspool%fT(arraysize),          &
           mimicspool%fW(arraysize),          &

           mimicspool%LITmN(arraysize),       &
           mimicspool%LITsN(arraysize),       &
           mimicspool%MICrN(arraysize),       &
           mimicspool%MICkN(arraysize),       &
           mimicspool%SOMaN(arraysize),       &
           mimicspool%SOMcN(arraysize),       &
           mimicspool%SOMpN(arraysize),       &
           mimicspool%DIN(arraysize))

  ALLOCATE(mimicspoolAn%ClitterAn(arraysize,nmimicslit),   &
           mimicspoolAn%CsoilAn(arraysize,nmimicssom),     &
           mimicspoolAn%CmicrobeAn(arraysize,nmimicsmcrb), &
           mimicspoolAn%thetaLiqAn(arraysize),             &
           mimicspoolAn%thetaFrznAn(arraysize),            &
           mimicspoolAn%fTAn(arraysize),                   &
           mimicspoolAn%fWAn(arraysize),                   &
           mimicspoolAn%DINAn(arraysize))

  ALLOCATE(mimicspoolAn%NlitterAn(arraysize,nmimicslit),   &
           mimicspoolAn%NsoilAn(arraysize,nmimicssom),     &
           mimicspoolAn%NmicrobeAn(arraysize,nmimicsmcrb))

  ALLOCATE(mimicsfluxAn%ClitInputAn(arraysize,nmimicslit))
  ALLOCATE(mimicsfluxAn%ChrespAn(arraysize))
  ALLOCATE(mimicsfluxAn%CSOMpInputAn(arraysize))
  ALLOCATE(mimicsfluxAn%Overflow_rAn(arraysize))
  ALLOCATE(mimicsfluxAn%Overflow_kAn(arraysize))
  ALLOCATE(mimicsfluxAn%NlitInputAn(arraysize,nmimicslit))

END SUBROUTINE alloc_mimicsvariable

END MODULE mimicsvariable
!-------------------------------------------------------------------------------- 
