!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: defines/allocates variables for CASA-CNP
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed for offline CASA-CNP, code revision likely to better
!          suit ACCESS and to merge more consistently with CABLE code
!
! ==============================================================================
! casa_variable.f90
!
! This file contains:
!   MODULE casadimension
!   MODULE casaparm
!   MODULE casavariable
!     with subroutine alloc_casavariable
!   MODULE phenvariable
!     with subroutine alloc_phenvariable
!   MODULE clmgridvariable
!     with subroutine alloc_CLMgridVariable
!
! the following modules are used when "casacnp" is coupled to "cable"
!   casadimension
!   casaparm
!   casavariable with subroutine alloc_casavariable
!   phenvariable with subroutine alloc_phenvariable

!================================================================================
MODULE casadimension
   
   IMPLICIT NONE
  
  !------------------------------------------------------------------------------------------------
  ! Added this section (-MDH 6/9/2014)
  !
  ! i_d is default kind for representing integer values.
  INTEGER, PARAMETER :: i_d = KIND(9)
  ! r_1 is default kind for representing REAL values (typically 32 or 64 bits).
  INTEGER, PARAMETER :: r_1  = KIND(1.0)
  ! r_2 is kind for representing REAL values with at least 10-digit precision
  ! (typically 64 bits).
  INTEGER, PARAMETER :: r_2  = SELECTED_REAL_KIND(12, 50)
  !------------------------------------------------------------------------------------------------
  
  INTEGER, PARAMETER :: mdyear=365         ! days per year
  INTEGER, PARAMETER :: mdmonth=30         ! days per month
  INTEGER, PARAMETER :: mdweek=7           ! days per week
  INTEGER, PARAMETER :: mmyear=12          ! month per year
  INTEGER, PARAMETER :: mt=36500           ! integration time step
  INTEGER, PARAMETER :: mpftmax=2          ! max. PFT/cell
  INTEGER, PARAMETER :: mplant = 3         ! plant pools
  INTEGER, PARAMETER :: mlitter= 3         ! litter pools
  INTEGER, PARAMETER :: msoil  = 3         ! soil pools
  INTEGER, PARAMETER :: mso    = 12        ! soil order number
! Put icycle into namelist file (icycle=1 is default) -MDH 9/21/2014
  INTEGER            :: icycle=1           ! =1 for C, =2 for C+N; =3 for C+N+P
  INTEGER, PARAMETER :: mstart=1           ! starting time step
  INTEGER, PARAMETER :: mphase=4           ! phen. phases
  REAL(r_2),    PARAMETER :: deltcasa=1.0/365.0 ! year
  REAL(r_2),    PARAMETER :: deltpool=1.0       ! pool delt(1day)

  INTEGER, PARAMETER :: CASACNP=1 
  INTEGER, PARAMETER :: MIMICS=2 
  INTEGER, PARAMETER :: CORPSE=3
  !isomModel is new option in namelist file, SOM model: 1=CASACNP, 2=MIMICS, 3=CORPSE) -MDH 2/8/2016 
  INTEGER :: isomModel = CASACNP

END MODULE casadimension

!================================================================================
MODULE casaparm
  USE casadimension

  !! Comment out the next 2 lines.  initcasa is read from startup file (-MDH 6/9/2014)
  integer :: initcasa = 0               ! =0 spin; 1 restart file; 2 transient simulation
  !! Updated cropland, iceland values (-MDH 7/7/2014)
  INTEGER, PARAMETER :: iceland  = 15   ! =13 for casa vegtype =15 for IGBP vegtype
  !! Added "barren" and "tundra" (-MDH 6/9/2014)
  integer, parameter :: barren   = 16   ! IGBP barren or sparsely vegetated
  integer, parameter :: tundra  = 18    ! IGBP+tundra classification
  INTEGER, PARAMETER :: cropland = 12   ! IGBP cropland
  INTEGER, PARAMETER :: croplnd2 = 14   ! IGBP cropland mosaic

  !! Vegetation types
  INTEGER, PARAMETER :: forest  = 3
  INTEGER, PARAMETER :: shrub   = 2
  INTEGER, PARAMETER :: grass   = 1
  INTEGER, PARAMETER :: icewater= 0

  INTEGER, PARAMETER :: LEAF    = 1
  INTEGER, PARAMETER :: WOOD    = 2
  INTEGER, PARAMETER :: FROOT   = 3
! INTEGER, PARAMETER :: LABILE  = 4
  INTEGER, PARAMETER :: METB    = 1
  INTEGER, PARAMETER :: STR     = 2
  INTEGER, PARAMETER :: CWD     = 3
  INTEGER, PARAMETER :: MIC     = 1
  INTEGER, PARAMETER :: SLOW    = 2
  INTEGER, PARAMETER :: PASS    = 3
  INTEGER, PARAMETER :: PLAB    = 1
  INTEGER, PARAMETER :: PSORB   = 2
  INTEGER, PARAMETER :: POCC    = 3
  INTEGER, PARAMETER :: LALLOC  = 0      !C allocation: fixed=0; dynamic=1; wolf=2
  REAL(r_2), PARAMETER :: z30=0.3
  REAL(r_2), PARAMETER :: R0=0.3
  REAL(r_2), PARAMETER :: S0=0.3
  REAL(r_2), PARAMETER :: fixed_stem=1.0/3.0
  REAL(r_2), PARAMETER :: Q10alloc=2.0
  REAL(r_2), PARAMETER :: ratioNCstrfix = 1.0/150.0
  REAL(r_2), PARAMETER :: ratioNPstrfix = 25.0                  
  REAL(r_2), PARAMETER :: fracCbiomass = 0.50
  REAL(r_2), PARAMETER :: tsoilrefc=25.0
  REAL(r_2), PARAMETER :: tkzeroc=273.15
  REAL(r_2), PARAMETER :: frootparma = 0.3192
  REAL(r_2), PARAMETER :: frootparmb =-0.0485
  REAL(r_2), PARAMETER :: frootparmc = 0.1755
  REAL(r_2), PARAMETER :: xweightalloc = 0.2
!  REAL(r_2), PARAMETER :: xkplab=0.5*deltcasa
!  REAL(r_2), PARAMETER :: xkpsorb=0.01*deltcasa
!  REAL(r_2), PARAMETER :: xkpocc =0.01*deltcasa
END MODULE casaparm

!================================================================================
MODULE casavariable
  USE casadimension
  IMPLICIT NONE

  TYPE casa_biome
    INTEGER,   DIMENSION(:),POINTER :: ivt2
    REAL(r_2), DIMENSION(:),POINTER :: xkleafcoldmax,  &
                                       xkleafcoldexp,  &
                                       xkleafdrymax,   &
                                       xkleafdryexp,   &
                                       glaimax,        &
                                       glaimin,        &
                                       sla,            &
                                       ratiofrootleaf, &
                                       kroot,          &
                                       krootlen,       &
                                       rootdepth,      &
                                       kuptake,        &
                                       kminN,          &
                                       kuplabP,        &
                                       kclabrate,      &
                                       xnpmax,         &
                                       q10soil,        &
                                       xkoptlitter,    &
                                       xkoptsoil,      &
                                       xkplab,         &
                                       xkpsorb,        &
                                       xkpocc,         &
                                       prodptase,      &
                                       costnpup,       &
                                       maxfinelitter,  &
                                       maxcwd,         &             
                                       nintercept,     &  
                                       nslope,         &

                                       xkNlimitmin,    &             
                                       xkNlimitmax,    &             
                                       fracRootExud,   &             
                                       CNRootExud,     &    
                                       CUEmetbmic,     &         
                                       CUEstrmic,      &         
                                       CUEstrslow,     &         
                                       CUEcwdmic,      &         
                                       CUEcwdslow,     &         
                                       CUEmicslow,     &         
                                       CUEmicpass,     &         
                                       CUEpassslow

    REAL(r_2), DIMENSION(:,:),POINTER :: plantrate,     &
                                       rmplant,         &
                                       fracnpptoP,      &
                                       fraclignin,      &
                                       fraclabile,      &
                                       ratioNCplantmin, &
                                       ratioNCplantmax, &
                                       ratioNPplantmin, &
                                       ratioNPplantmax, &
                                       fracLigninplant, &
                                       ftransNPtoL,     &
                                       ftransPPtoL,     &
                                       litterrate
    REAL(r_2), DIMENSION(:,:),POINTER :: soilrate
  END TYPE casa_biome

  TYPE casa_pool
    REAL(r_2), DIMENSION(:),POINTER :: Clabile,       &
                                       dClabiledt               
    REAL(r_2), DIMENSION(:,:),POINTER :: Cplant,      &
                                       Nplant,        &
                                       Pplant,        &
                                       dCplantdt,     &
                                       dNplantdt,     &
                                       dPplantdt,     &
                                       ratioNCplant,  &
                                       ratioNPplant
    REAL(r_2), DIMENSION(:),POINTER :: Nsoilmin,      &
                                       Psoillab,      &
                                       Psoilsorb,     &
                                       Psoilocc,      & 
                                       dNsoilmindt,   &
                                       dPsoillabdt,   &
                                       dPsoilsorbdt,  &
                                       dPsoiloccdt    
    REAL(r_2), DIMENSION(:,:), POINTER :: Clitter,    &
                                       Nlitter,       &
                                       Plitter,       &
                                       dClitterdt,    &
                                       dNlitterdt,    &
                                       dPlitterdt,    &
                                       ratioNClitter, &
                                       ratioNPlitter
    REAL(r_2), DIMENSION(:,:),POINTER :: Csoil,       &
                                       Nsoil,         &
                                       Psoil,         &
                                       dCsoildt,      &
                                       dNsoildt,      &
                                       dPsoildt,      &
                                       ratioNCsoil,   &
                                       ratioNCsoilnew,&
                                       ratioNPsoil,   &
                                       ratioNCsoilmin,&
                                       ratioNCsoilmax

    REAL(r_2), DIMENSION(:),POINTER :: fT, fW, thetaLiq

  END TYPE casa_pool

  !! casa_annual_pool is used to store average annual pools
  !! for output to casa_pool_flux.nc file only. -MDH 9/29/2014
  !! Added average annual air and soil temperature. -MDH 08/17/2015
  TYPE casa_annual_pool
      
    REAL(r_2), DIMENSION(:,:),POINTER :: CplantAn,    &
                                       NplantAn,      &
                                       PplantAn

    
    REAL(r_2), DIMENSION(:,:), POINTER :: ClitterAn,  &
                                       NlitterAn,     &
                                       PlitterAn

    REAL(r_2), DIMENSION(:,:),POINTER :: CsoilAn,     &
                                       NsoilAn,       &
                                       PsoilAn

    REAL(r_2), DIMENSION(:),POINTER :: tairAn,        &
                                       tsoilAn,       &
                                       fTAn,          &
                                       fWAn,          &
                                       thetaLiqAn

    REAL(r_2), DIMENSION(:), POINTER :: NsoilminAn


  END TYPE casa_annual_pool

  TYPE casa_flux
    REAL(r_2), DIMENSION(:),POINTER :: Cgpp,          &
                                       Cnpp,          &
                                       CnppAn,        &
                                       Crp,           &
                                       Crgplant,      &
                                       ClitInptMet,        &
                                       ClitInptMetAn,      &
                                       ClitInptStruc,      &
                                       ClitInptStrucAn,    &
                                       NlitInptMet,        &
                                       NlitInptMetAn,      &
                                       NlitInptStruc,      &
                                       NlitInptStrucAn,    &
                                       Nminuptake,    &
                                       NminuptakeAn,  &
                                       Plabuptake,    &
                                       Clabloss,      &
                                       fracClabile,   &
                                       xkNlimiting,   &  
                                       CpassInpt,     &  
                                       CpassInptAn   
    REAL(r_2), DIMENSION(:,:),POINTER :: fracCalloc,  &
                                       fracNalloc,    &
                                       fracPalloc,    &
                                       Crmplant,      &
                                       kplant
    REAL(r_2), DIMENSION(:,:,:),POINTER :: fromPtoL
    REAL(r_2), DIMENSION(:),POINTER :: Cnep,        &
                                       Crsoil,      &
                                       Nmindep,     &
                                       NmindepAn,   &
                                       Nminfix,       &
                                       NminfixAn,     &
                                       Nminloss,    &
                                       NminlossAn,  &
                                       Nminleach,   &
                                       NminleachAn, &
                                       Nlittermin,  &
                                       NlitterminAn,&
                                       Nsmin,       &
                                       NsminAn,     &
                                       Nsimm,       &
                                       NsimmAn,     &
                                       Nsnet,       &
                                       NsnetAn,     &
                                       Nupland,     &
                                       fNminloss,   &
                                       fNminleach,  &
                                       Pdep,        &
                                       Pwea,        &
                                       Pleach,      &
                                       Ploss,       &
                                       Pupland,     &
                                       Plittermin,  &
                                       Psmin,       &
                                       Psimm,       &
                                       Psnet,       &
                                       fPleach,     &
                                       kplab,       &
                                       kpsorb,      &
                                       kpocc,       &
                                       kmlabp,      &
                                       Psorbmax
    REAL(r_2), DIMENSION(:,:),POINTER    :: klitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: ksoil
    REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromLtoS
    REAL(r_2), DIMENSION(:,:,:),POINTER  :: fromStoS
    REAL(r_2), DIMENSION(:,:),POINTER    :: fromLtoCO2
    REAL(r_2), DIMENSION(:,:),POINTER    :: fromStoCO2
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtolitter
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxCtosoil
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxNtosoil
    REAL(r_2), DIMENSION(:,:),POINTER    :: FluxPtosoil
    REAL(r_2), DIMENSION(:),POINTER      :: FluxCtoCO2
  END TYPE casa_flux

  TYPE casa_met
    REAL(r_2), DIMENSION(:),POINTER    :: glai,     &
                                          Tairk,    &
                                          precip,   &
                                          tsoilavg, &
                                          moistavg, &
                                          frznmoistavg, &
                                          btran
    INTEGER, DIMENSION(:), POINTER     :: lnonwood
    REAL(r_2), DIMENSION(:,:), POINTER :: Tsoil,     &
                                          frznmoist, &
                                          moist
! Added ilat(*) and ilon(*) to casamet so that casacnp output can be
! remapped to a 2D grid (-Melannie Hartman 3/15/2014).
    INTEGER, DIMENSION(:), POINTER     :: iveg2,    &  
                                          ijgcm,    &
                                          isorder,  &
                                          ilat, ilon
    REAL(r_2), DIMENSION(:), POINTER   :: lat,      &
                                          lon,      &
                                          areacell
    ! added yp wang 5/nov/2012
!!  ! Commented out. Not needed. -mdh 3/13/2017
!!  REAL(r_2), DIMENSION(:,:), POINTER :: Tairkspin,&
!!                                        cgppspin,&
!!                                        crmplantspin_1,&
!!                                        crmplantspin_2,&
!!                                        crmplantspin_3,&
!!                                        Tsoilspin_1,&
!!                                        Tsoilspin_2,&
!!                                        Tsoilspin_3,&
!!                                        Tsoilspin_4,&
!!                                        Tsoilspin_5,&
!!                                        Tsoilspin_6,&
!!                                        moistspin_1,&
!!                                        moistspin_2,&
!!                                        moistspin_3,&
!!                                        moistspin_4,&
!!                                        moistspin_5,&
!!                                        moistspin_6

  END TYPE casa_met

  TYPE casa_balance
    REAL(r_2), DIMENSION(:),POINTER   :: FCgppyear,FCnppyear,                 &
            FCrmleafyear,FCrmwoodyear,FCrmrootyear,FCrgrowyear,               &
            FCrpyear, FCrsyear,FCneeyear,                                     &
            FNdepyear,FNfixyear, FNsnetyear,FNupyear, FNleachyear,FNlossyear, &
            FPweayear,FPdustyear,FPsnetyear,FPupyear, FPleachyear,FPlossyear
    REAL(r_2), DIMENSION(:,:),POINTER :: glaimon,glaimonx
    REAL(r_2), DIMENSION(:,:),POINTER :: cplantlast,nplantlast,pplantlast
    REAL(r_2), DIMENSION(:,:),POINTER :: clitterlast,nlitterlast,plitterlast
    REAL(r_2), DIMENSION(:,:),POINTER :: csoillast,nsoillast,psoillast
    REAL(r_2), DIMENSION(:),  POINTER :: nsoilminlast,psoillablast,  &
                                         psoilsorblast,psoilocclast, &
                                         cbalance,nbalance,pbalance, &
                                         sumcbal,sumnbal,sumpbal
    REAL(r_2), DIMENSION(:),POINTER   :: clabilelast
  END TYPE casa_balance

! Added filename type for casaCNP (BP apr2010)
  TYPE casafiles_type
!   CHARACTER(LEN=100) :: cnpbiome    ! file for biome-specific BGC parameters
!   CHARACTER(LEN=100) :: cnppoint    ! file for point-specific BGC inputs
!   CHARACTER(LEN=100) :: cnpepool    ! file for end-of-run pool sizes
!   CHARACTER(LEN=100) :: cnpipool    ! file for inital pool sizes
!   CHARACTER(LEN=100) :: cnpmetin      ! met file for spin up 
!   CHARACTER(LEN=100) :: cnpmetout     ! met file for spin up 
!   CHARACTER(LEN=100) :: ndep          ! N deposition input file   
! added yp wang
!   CHARACTER(LEN=100) :: cnpspin       ! input file for spin up
!   CHARACTER(LEN=100) :: dump_cnpspin  ! name of dump file for spinning casa-cnp
!
!   CHARACTER(LEN=100) :: phen        ! leaf phenology datafile
!   CHARACTER(LEN=100) :: cnpflux     ! modelled mean yearly CNP fluxes

    CHARACTER(LEN=100) :: sPtFileNameCASA
    INTEGER :: iptToSave = 0
    INTEGER :: iptToSaveIndx = 0
!   LOGICAL           :: l_ndep
  END TYPE casafiles_type
  TYPE(casafiles_type) :: casafile

! Inserted this section (-MDH 6/9/2014)
  TYPE (casa_biome)              :: casabiome
  TYPE (casa_pool)               :: casapool
  TYPE (casa_flux)               :: casaflux
  TYPE (casa_met)                :: casamet
  TYPE (casa_balance)            :: casabal
! Inserted new type (-MDH 9/29/2014)
  TYPE (casa_annual_pool)        :: casapoolAn

  character(len=100) :: filename_cnpbiome !! file for biome-specific CASACCNP parameters

Contains


! Inserted the previous subroutine declaration (-MDH 6/9/2014)
SUBROUTINE alloc_casavariable(mp,mvtype,ms)
  USE casadimension
  implicit none
  integer, INTENT(IN)  :: mp, mvtype, ms
  integer :: arraysize
  arraysize = mp

  ALLOCATE(casabiome%ivt2(mvtype),                   &
           casabiome%xkleafcoldmax(mvtype),          &
           casabiome%xkleafcoldexp(mvtype),          &
           casabiome%xkleafdrymax(mvtype),           &
           casabiome%xkleafdryexp(mvtype),           &
           casabiome%glaimax(mvtype),                &
           casabiome%glaimin(mvtype),                &
           casabiome%sla(mvtype),                    &
           casabiome%ratiofrootleaf(mvtype),         &
           casabiome%kroot(mvtype),                  &
           casabiome%krootlen(mvtype),               &
           casabiome%rootdepth(mvtype),              & 
           casabiome%kuptake(mvtype),                &
           casabiome%kminN(mvtype),                  &
           casabiome%KuplabP(mvtype),                &
           casabiome%kclabrate(mvtype),              &
           casabiome%xnpmax(mvtype),                 &
           casabiome%q10soil(mvtype),                &
           casabiome%xkoptlitter(mvtype),            &
           casabiome%xkoptsoil(mvtype),              &
           casabiome%xkplab(mso),                    &
           casabiome%xkpsorb(mso),                   &
           casabiome%xkpocc(mso),                    &
           casabiome%prodptase(mvtype),              &
           casabiome%costnpup(mvtype),               &
           casabiome%maxfinelitter(mvtype),          &
           casabiome%maxcwd(mvtype),                 &
           casabiome%nintercept(mvtype),             &
           casabiome%nslope(mvtype),                 &

           casabiome%xkNlimitmin(mvtype),            &             
           casabiome%xkNlimitmax(mvtype),            &             
           casabiome%fracRootExud(mvtype),           &             
           casabiome%CNRootExud(mvtype),             &    
           casabiome%CUEmetbmic(mvtype),             &         
           casabiome%CUEstrmic(mvtype),              &         
           casabiome%CUEstrslow(mvtype),             &         
           casabiome%CUEcwdmic(mvtype),              &         
           casabiome%CUEcwdslow(mvtype),             &         
           casabiome%CUEmicslow(mvtype),             &         
           casabiome%CUEmicpass(mvtype),             &         
           casabiome%CUEpassslow(mvtype),            &

           casabiome%plantrate(mvtype,mplant),       &
           casabiome%rmplant(mvtype,mplant),         &
           casabiome%fracnpptoP(mvtype,mplant),      &
           casabiome%fraclignin(mvtype,mplant),      &
           casabiome%fraclabile(mvtype,mplant),      &
           casabiome%ratioNCplantmin(mvtype,mplant), &
           casabiome%ratioNCplantmax(mvtype,mplant), &
           casabiome%ratioNPplantmin(mvtype,mplant), &
           casabiome%ratioNPplantmax(mvtype,mplant), &
           casabiome%fracLigninplant(mvtype,mplant), &
           casabiome%ftransNPtoL(mvtype,mplant),     &
           casabiome%ftransPPtoL(mvtype,mplant),     &
           casabiome%litterrate(mvtype,mlitter),     &
           casabiome%soilrate(mvtype,msoil))

  ALLOCATE(casapool%Clabile(arraysize),               &
           casapool%dClabiledt(arraysize),            &
           casapool%Cplant(arraysize,mplant),         &
           casapool%Nplant(arraysize,mplant),         &
           casapool%Pplant(arraysize,mplant),         &
           casapool%dCplantdt(arraysize,mplant),      &
           casapool%dNplantdt(arraysize,mplant),      &
           casapool%dPplantdt(arraysize,mplant),      &
           casapool%ratioNCplant(arraysize,mplant),   &
           casapool%ratioNPplant(arraysize,mplant),   &
           casapool%Nsoilmin(arraysize),              &
           casapool%Psoillab(arraysize),              &
           casapool%Psoilsorb(arraysize),             &
           casapool%Psoilocc(arraysize),              &
           casapool%dNsoilmindt(arraysize),           &
           casapool%dPsoillabdt(arraysize),           &
           casapool%dPsoilsorbdt(arraysize),          &
           casapool%dPsoiloccdt(arraysize),           &
           casapool%Clitter(arraysize,mlitter),       &
           casapool%Nlitter(arraysize,mlitter),       &
           casapool%Plitter(arraysize,mlitter),       &
           casapool%dClitterdt(arraysize,mlitter),    &
           casapool%dNlitterdt(arraysize,mlitter),    &
           casapool%dPlitterdt(arraysize,mlitter),    &
           casapool%ratioNClitter(arraysize,mlitter), &
           casapool%ratioNPlitter(arraysize,mlitter), &
           casapool%Csoil(arraysize,msoil),           &
           casapool%Nsoil(arraysize,msoil),           &
           casapool%Psoil(arraysize,msoil),           &
           casapool%dCsoildt(arraysize,msoil),        &
           casapool%dNsoildt(arraysize,msoil),        &
           casapool%dPsoildt(arraysize,msoil),        &
           casapool%ratioNCsoil(arraysize,msoil),     &
           casapool%ratioNPsoil(arraysize,msoil),     &
           casapool%ratioNCsoilnew(arraysize,msoil),  &
           casapool%ratioNCsoilmin(arraysize,msoil),  &
           casapool%ratioNCsoilmax(arraysize,msoil),  &
           casapool%fT(arraysize),                    &
           casapool%fW(arraysize),                    &
           casapool%thetaLiq(arraysize))

  !! Added new variable type. -MDH 9/29/2014
  ALLOCATE(casapoolAn%CplantAn(arraysize,mplant),        &
           casapoolAn%NplantAn(arraysize,mplant),        &
           casapoolAn%PplantAn(arraysize,mplant),        &
           casapoolAn%CsoilAn(arraysize,msoil),          &
           casapoolAn%NsoilAn(arraysize,msoil),          &
           casapoolAn%PsoilAn(arraysize,msoil),          &
           casapoolAn%ClitterAn(arraysize,mlitter),      &
           casapoolAn%NlitterAn(arraysize,mlitter),      &
           casapoolAn%PlitterAn(arraysize,mlitter),      &
           casapoolAn%tairAn(arraysize),                 &
           casapoolAn%tsoilAn(arraysize),                &
           casapoolAn%fTAn(arraysize),                   &
           casapoolAn%fWAn(arraysize),                   &
           casapoolAn%thetaLiqAn(arraysize),             &
           casapoolAn%NsoilminAn(arraysize))


  ALLOCATE(casaflux%Cgpp(arraysize),                     &
           casaflux%Cnpp(arraysize),                     &
           casaflux%CnppAn(arraysize),                   &
           casaflux%Crp(arraysize),                      &
           casaflux%Crgplant(arraysize),                 &
           casaflux%ClitInptMet(arraysize),              &
           casaflux%ClitInptMetAn(arraysize),            &
           casaflux%ClitInptStruc(arraysize),            &
           casaflux%ClitInptStrucAn(arraysize),          &
           casaflux%NlitInptMet(arraysize),              &
           casaflux%NlitInptMetAn(arraysize),            &
           casaflux%NlitInptStruc(arraysize),            &
           casaflux%NlitInptStrucAn(arraysize),          &

           casaflux%xkNlimiting(arraysize),              &
           casaflux%CpassInpt(arraysize),                &
           casaflux%CpassInptAn(arraysize),              &
           casaflux%Nminfix(arraysize),                  &
           casaflux%NminfixAn(arraysize),                &
           casaflux%Nminuptake(arraysize),               &
           casaflux%NminuptakeAn(arraysize),             &
           casaflux%Plabuptake(arraysize),               &
           casaflux%Clabloss(arraysize),                 &
           casaflux%fracClabile(arraysize),              &
           casaflux%fracCalloc(arraysize,mplant),        &
           casaflux%fracNalloc(arraysize,mplant),        &
           casaflux%fracPalloc(arraysize,mplant),        &
           casaflux%kplant(arraysize,mplant),            &
           casaflux%Crmplant(arraysize,mplant),          &
           casaflux%fromPtoL(arraysize,mlitter,mplant),  &
           casaflux%Cnep(arraysize),                     &
           casaflux%Crsoil(arraysize),                   &
           casaflux%Nmindep(arraysize),                  &
           casaflux%NmindepAn(arraysize),                &
           casaflux%Nminloss(arraysize),                 &
           casaflux%NminlossAn(arraysize),               &
           casaflux%Nminleach(arraysize),                &
           casaflux%NminleachAn(arraysize),              &
           casaflux%Nupland(arraysize),                  &
           casaflux%Nlittermin(arraysize),               &
           casaflux%NlitterminAn(arraysize),             &
           casaflux%Nsmin(arraysize),                    &
           casaflux%NsminAn(arraysize),                  &
           casaflux%Nsimm(arraysize),                    &
           casaflux%NsimmAn(arraysize),                  &
           casaflux%Nsnet(arraysize),                    &
           casaflux%NsnetAn(arraysize),                  &
           casaflux%fNminloss(arraysize),                &
           casaflux%fNminleach(arraysize),               &
           casaflux%Pdep(arraysize),                     &
           casaflux%Pwea(arraysize),                     &
           casaflux%Pleach(arraysize),                   &
           casaflux%Ploss(arraysize),                    &
           casaflux%Pupland(arraysize),                  &
           casaflux%Plittermin(arraysize),               &
           casaflux%Psmin(arraysize),                    &
           casaflux%Psimm(arraysize),                    &
           casaflux%Psnet(arraysize),                    &
           casaflux%fPleach(arraysize),                  &
           casaflux%kplab(arraysize),                    &
           casaflux%kpsorb(arraysize),                   &
           casaflux%kpocc(arraysize),                    &
           casaflux%kmlabP(arraysize),                   &
           casaflux%Psorbmax(arraysize),                 &
           casaflux%klitter(arraysize,mlitter),          &
           casaflux%ksoil(arraysize,msoil),              &
           casaflux%fromLtoS(arraysize,msoil,mlitter),   &
           casaflux%fromStoS(arraysize,msoil,msoil),     &
           casaflux%fromLtoCO2(arraysize,mlitter),       &
           casaflux%fromStoCO2(arraysize,msoil))

  ALLOCATE(casaflux%FluxCtolitter(arraysize,mlitter),    &
           casaflux%FluxNtolitter(arraysize,mlitter),    &
           casaflux%FluxPtolitter(arraysize,mlitter))

  ALLOCATE(casaflux%FluxCtosoil(arraysize,msoil),        &
           casaflux%FluxNtosoil(arraysize,msoil),        &
           casaflux%FluxPtosoil(arraysize,msoil))

  ALLOCATE(casaflux%FluxCtoco2(arraysize))

! Added ilat(*) and ilon(*) to casamet so that casacnp output can be
! remapped to a 2D grid (-Melannie Hartman 3/15/2014).
! frznmoist and frznmostavg added 3/13/2017. -mdh
  ALLOCATE(casamet%glai(arraysize),                &
           casamet%lnonwood(arraysize),            &
           casamet%Tairk(arraysize),               &
           casamet%precip(arraysize),              &
           casamet%tsoilavg(arraysize),            &
           casamet%moistavg(arraysize),            &
           casamet%frznmoistavg(arraysize),        &
           casamet%btran(arraysize),               &
           casamet%Tsoil(arraysize,ms),            &
           casamet%moist(arraysize,ms),            &
           casamet%frznmoist(arraysize,ms),        &
           casamet%iveg2(arraysize),               &  
           casamet%ijgcm(arraysize),               &
           casamet%isorder(arraysize),             &
           casamet%lat(arraysize),                 &
           casamet%lon(arraysize),                 &
           casamet%areacell(arraysize),            &
           casamet%ilat(arraysize),                &
           casamet%ilon(arraysize))                 

!! These were commented out. Not needed. -mdh 3/13/2017
!!         casamet%Tairkspin(arraysize,mdyear),     &
!!         casamet%cgppspin(arraysize,mdyear),      &
!!         casamet%crmplantspin_1(arraysize,mdyear),&
!!         casamet%crmplantspin_2(arraysize,mdyear),&
!!         casamet%crmplantspin_3(arraysize,mdyear),&
!!         casamet%Tsoilspin_1(arraysize,mdyear),   &
!!         casamet%Tsoilspin_2(arraysize,mdyear),   &
!!         casamet%Tsoilspin_3(arraysize,mdyear),   &
!!         casamet%Tsoilspin_4(arraysize,mdyear),   &
!!         casamet%Tsoilspin_5(arraysize,mdyear),   &
!!         casamet%Tsoilspin_6(arraysize,mdyear),   &
!!         casamet%moistspin_1(arraysize,mdyear),   &
!!         casamet%moistspin_2(arraysize,mdyear),   &
!!         casamet%moistspin_3(arraysize,mdyear),   &
!!         casamet%moistspin_4(arraysize,mdyear),   &
!!         casamet%moistspin_5(arraysize,mdyear),   &
!!         casamet%moistspin_6(arraysize,mdyear))


  ALLOCATE(casabal%FCgppyear(arraysize),           &
           casabal%FCnppyear(arraysize),           &
           casabal%FCrpyear(arraysize),            &
           casabal%FCrmleafyear(arraysize),        &
           casabal%FCrmwoodyear(arraysize),        &
           casabal%FCrmrootyear(arraysize),        &
           casabal%FCrgrowyear(arraysize),         &
           casabal%FCrsyear(arraysize),            &
           casabal%FCneeyear(arraysize),           &
           casabal%FNdepyear(arraysize),           &
           casabal%FNfixyear(arraysize),           &
           casabal%FNsnetyear(arraysize),          &
           casabal%FNupyear(arraysize),            &
           casabal%FNleachyear(arraysize),         &
           casabal%FNlossyear(arraysize),          &
           casabal%FPweayear(arraysize),           &
           casabal%FPdustyear(arraysize),          &
           casabal%FPsnetyear(arraysize),          &
           casabal%FPupyear(arraysize),            &
           casabal%FPleachyear(arraysize),         &
           casabal%FPlossyear(arraysize))
            
  ALLOCATE(casabal%glaimon(arraysize,12),          &
           casabal%glaimonx(arraysize,12))

  ALLOCATE(casabal%cplantlast(arraysize,mplant),   &
           casabal%nplantlast(arraysize,mplant),   &
           casabal%pplantlast(arraysize,mplant))

  ALLOCATE(casabal%clitterlast(arraysize,mlitter), &
           casabal%nlitterlast(arraysize,mlitter), &
           casabal%plitterlast(arraysize,mlitter))

  ALLOCATE(casabal%csoillast(arraysize,msoil),     &
           casabal%nsoillast(arraysize,msoil),     &
           casabal%psoillast(arraysize,msoil))

  ALLOCATE(casabal%nsoilminlast(arraysize),        &
           casabal%psoillablast(arraysize),        &
           casabal%psoilsorblast(arraysize),       &
           casabal%psoilocclast(arraysize),        &
           casabal%cbalance(arraysize),            &
           casabal%nbalance(arraysize),            &
           casabal%pbalance(arraysize),            &
           casabal%sumcbal(arraysize),             &
           casabal%sumnbal(arraysize),             &
           casabal%sumpbal(arraysize),             &
           casabal%clabilelast(arraysize))
END SUBROUTINE alloc_casavariable

END MODULE casavariable

!================================================================================
MODULE phenvariable
  USE casadimension
  IMPLICIT NONE
  TYPE phen_variable
    INTEGER,   DIMENSION(:),  POINTER :: phase        
    REAL(r_2), DIMENSION(:),  POINTER :: TKshed
    INTEGER,   DIMENSION(:,:),POINTER :: doyphase
  END type phen_variable

! Inserted this statement (-MDH 6/9/2014)
  TYPE(phen_variable)        :: phen  

CONTAINS

! Updated this subroutine declaration (-MDH 6/9/2014)
SUBROUTINE alloc_phenvariable(mp,mvtype)
  INTEGER,             INTENT(IN) :: mp, mvtype
  INTEGER                         :: arraysize
  arraysize = mp

  ALLOCATE(phen%Tkshed(mvtype))
  ALLOCATE(phen%phase(arraysize),         &
           phen%doyphase(arraysize,mphase))
END SUBROUTINE alloc_phenvariable

End MODULE phenvariable

!================================================================================
! MODULE clmgridvariable has been added in order to run the
! CASACNP model with surface and meteorological data derived from CLM
! history files (Melannie D. Hartman, March 2014)
! 
MODULE clmgridvariable
   use casadimension
   implicit none

   TYPE clm_grid_variable
      integer :: nlat, nlon
      real(4), dimension(:), pointer :: lat1d          ! lon1d(nlon) coordinate longitude (degrees east)
      real(4), dimension(:), pointer :: lon1d          ! lat1d(nlat) coordinate latitude (degrees north)
      integer, dimension(:,:), pointer :: cellMissing  ! cellMissing(nlon,nlat) 0=no missing data, 1=missing data
      integer, dimension(:,:), pointer :: cellid       ! cellid(nlon,nlat) grid cell ids (1..nlat*nlon)
!     real(4), dimension(:,:), pointer  :: landarea    ! landarea(nlon,nlat) area of land in grid cell (m^2)?
   END TYPE clm_grid_variable

   TYPE (clm_grid_variable) :: clmgrid

   CONTAINS

   SUBROUTINE alloc_CLMgridVariable(nlon, nlat)
      integer, intent(in) :: nlon, nlat
      clmgrid%nlat = nlat
      clmgrid%nlon = nlon
      allocate(clmgrid%lat1d(nlat))
      allocate(clmgrid%lon1d(nlon))
      allocate(clmgrid%cellMissing(nlon,nlat))
      allocate(clmgrid%cellid(nlon,nlat))
!     allocate(clmgrid%landarea(nlon,nlat))
   END SUBROUTINE alloc_clmgridvariable

END MODULE clmgridvariable

