!--------------------------------------------------------------------------------
! FILE: mimics_cycle.f90
!
! Purpose: 
!   Subroutines for C cycling in the MIMICS model
!
!     SUBROUTINE mimics_soil_forwardMM 
!         - calculate daily changes in litter, microbial, and SOM pools using forward
!           Michaelis-Menten dynamics
!         - replaces casa_coeffsoil and casa_delsoil when MIMICS or CORPSE models are used
!     SUBROUTINE mimics_soil_reverseMM
!         - calculate daily changes in litter, microbial, and SOM pools
!         - replaces casa_coeffsoil and casa_delsoil when MIMICS or CORPSE models are used
!         - alternative to mimics_soil_forwardMM
!         - Implements the reverse Michaelis-Menten kinetics
!     SUBROUTINE mimics_coeffplant
!         - calculate litter inputs to MIMICS model
!         - replaces casa_coeffplant when MIMICS model is used
!     SUBROUTINE mimics_delplant
!         - replaces casa_delplant when MIMICS model is used
!     SUBROUTINE mimics_xratesoil (-mdh 4/26/2015) 
!         - replaces casa_xratesoil and casa_coeffsoil when MIMICS or CORPSE models are used
!         - used to compute the dcomposition rate for CWD, casaflux%klitter(npt,cwd)
!     SUBROUTINE mimics_ccycle
!         - called daily
!         - update all C pools that were not updated in mimics_soil_*MM 
!         - replaces casa_cnpcyle when mimics or corpse models are used
!     SUBROUTINE mimics_caccum
!         - called daily
!         - accumulate annual C fluxes and pool values
!
! Contact: Melannie Hartman
!          melannie@ucar.edu
!
! History:
!   1/12/2015 - Created
!   8/24/2015 - Added reverse Michaelis-Menton option (subroutine mimics_soil_reverseMM)
!--------------------------------------------------------------------------------

MODULE mimics_cycle_module

USE define_dimensions
USE define_types
USE casadimension
USE casaparm
USE casavariable
USE phenvariable
USE mimicsdimension
USE mimicsparam
USE mimicsvariable
IMPLICIT NONE
CONTAINS

SUBROUTINE mimics_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                             casaflux,casamet)

! Calculate the plant litter fall rate, and litter inputs to the MIMICS model
!   * Modify casaflux%fromPtoL(npt,nL,nP) calculations
!
! inputs:
!     xkleafcold(mp):  cold stress induced leaf senescence rate (1/day)
!     xkleafdry(mp):   drought-induced leaf senescence rate (1/day)
!     xkleaf(mp):      set to 0.0 during maximal leaf growth phase
!
! outputs:
!     kplant(mp,mplant):           senescence rate of plant pool (1/day)
!     fromPtoL(mp,mlitter,mplant): fraction of senesced plant biomass to litter pool (fraction)

  REAL(r_2), DIMENSION(mp),    INTENT(IN)    :: xkleafcold,xkleafdry,xkleaf
  TYPE (veg_parameter_type),   INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),           INTENT(INOUT) :: casabiome
  TYPE (casa_pool),            INTENT(INOUT) :: casapool
  TYPE (casa_flux),            INTENT(INOUT) :: casaflux
  TYPE (casa_met),             INTENT(INOUT) :: casamet

  !Local variables
! REAL(r_2), DIMENSION(mp)  :: xk
  INTEGER :: npt    

  casaflux%fromPtoL(:,:,:)      = 0.0
  casaflux%kplant(:,:)          = 0.0   

  WHERE(casamet%iveg2 /= icewater)
      !Fraction of senesced plant biomass transferred to litter pools
      ! Use the fmet parameters also used in mimics_delplant (w.wieder 11/07/2016).
      !!casaflux%fromPtoL(:,metb,leaf)   = max(0.001, 0.85 - 0.013 * mimicsbiome%ligninNratio(:,leaf)) 
      !!casaflux%fromPtoL(:,metb,froot)  = max(0.001, 0.85 - 0.013 * mimicsbiome%ligninNratio(:,froot))
      casaflux%fromPtoL(:,metb,leaf)   = max(0.001, mimicsbiome%fmet_p(1)*(mimicsbiome%fmet_p(2)-mimicsbiome%fmet_p(3) &
                                                    *mimicsbiome%ligninNratio(:,leaf))) 
      casaflux%fromPtoL(:,metb,froot)  = max(0.001, mimicsbiome%fmet_p(1)*(mimicsbiome%fmet_p(2)-mimicsbiome%fmet_p(3) &
                                                    *mimicsbiome%ligninNratio(:,froot)))
      casaflux%fromPtoL(:,str,leaf)    = 1.0 - casaflux%fromPtoL(:,metb,leaf)
      casaflux%fromPtoL(:,str,froot)   = 1.0 - casaflux%fromPtoL(:,metb,froot)
      casaflux%fromPtoL(:,cwd,wood)    = 1.0

      !Live plant turnover (death) rates
      casaflux%kplant(:,leaf)  = casabiome%plantrate(veg%iveg(:),leaf)*xkleaf(:) &
                                 + xkleafcold(:) + xkleafdry(:)
      casaflux%kplant(:,wood)  = casabiome%plantrate(veg%iveg(:),wood) 
      casaflux%kplant(:,froot) = casabiome%plantrate(veg%iveg(:),froot) 
  ENDWHERE

  ! When glai<glaimin,leaf biomass will not decrease anymore. (Q.Zhang 10/03/2011)
  DO npt = 1,mp 
      if(casamet%glai(npt).le.casabiome%glaimin(veg%iveg(npt))) casaflux%kplant(npt,leaf) = 0.0
  ENDDO


END SUBROUTINE mimics_coeffplant

!--------------------------------------------------------------------------------
SUBROUTINE mimics_xratesoil(veg,soil,casamet,casabiome)
!  to account for effects of T and W on litter decomposition: xklitter
!  inputs:
!     ivt(mp)  :       biome type
!     tsoilavg(mp):    soil temperature in K
!     moistavg(mp):    volumetric soil moisture

  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome

  ! local variables
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps 
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)
! Kirschbaum function parameters
  REAL(r_2), parameter :: xkalpha=-3.764   ! Kirschbaum (1995, SBB)
  REAL(r_2), parameter :: xkbeta=0.204
  REAL(r_2), parameter :: xktoptc=36.9
  REAL(r_2), DIMENSION(mp)       :: xkwater,xktemp,xklitter
  REAL(r_2), DIMENSION(mp)       :: fwps,tsavg
  INTEGER :: npt

  xklitter(:) = 1.0
  casaflux%klitter(:,:) = 0.0        ! initialize klitter for all 3 casa litter pool.  Only cwd will be reset.
  casaflux%fromLtoCO2(:,:) = 0.0     ! flow from L to CO2

  fwps(:)     =  min(1.0, casamet%moistavg(:)/soil%ssat(:))
  tsavg(:)    =  casamet%tsoilavg(:) 

  DO npt=1,mp
    IF(casamet%iveg2(npt) /= icewater) THEN
      xktemp(npt)  = casabiome%q10soil(veg%iveg(npt))**(0.1*(tsavg(npt)-TKzeroC-35.0))
      xkwater(npt) = ((fwps(npt)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
                 * ((fwps(npt)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
      IF (veg%iveg(npt) == cropland .OR. veg%iveg(npt) == croplnd2) &
                 xkwater(npt)=1.0
      xklitter(npt) = casabiome%xkoptlitter(veg%iveg(npt)) * xktemp(npt) * xkwater(npt)
      casaflux%klitter(npt,cwd) = xklitter(npt) * casabiome%litterrate(veg%iveg(npt),cwd)   

      !! from casa_coeffsoil for reference
      !! casaflux%fromLtoS(:,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(:),wood)) ! CWD -> fmic
      !! casaflux%fromLtoS(:,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(:),wood)        ! CWD -> slow

      ! Fraction of cwd decomposition that goes to heterotrophic respiration
      casaflux%fromLtoCO2(npt,cwd) = 1.0 - 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(npt),wood)) &
                                         - 0.7 * casabiome%fracLigninplant(veg%iveg(npt),wood) 
      casaflux%fromLtoCO2(npt,cwd) = MAX(0.0, casaflux%fromLtoCO2(npt,cwd))

    END IF
  END DO

END SUBROUTINE mimics_xratesoil

!--------------------------------------------------------------------------------

SUBROUTINE mimics_delplant(veg,casabiome,casapool,casaflux,casamet,            &
                           cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                           nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                           pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                           cwd2co2,cwd2str)

! Calculate the change in plant C pools when the MIMICS model is used

  TYPE (veg_parameter_type), INTENT(INOUT) :: veg        ! vegetation parameters
  TYPE (casa_biome),         INTENT(INOUT) :: casabiome
  TYPE (casa_pool),          INTENT(INOUT) :: casapool
  TYPE (casa_flux),          INTENT(INOUT) :: casaflux
  TYPE (casa_met),           INTENT(INOUT) :: casamet

  real, dimension(mp),INTENT(OUT) :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                                     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                                     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                                     cwd2co2,cwd2str
  ! Local Variables
  integer ::  npt,nL,nP,nland
  real(r_2):: fracLeafMetabolic

  casaflux%FluxCtolitter = 0.0
  casaflux%FluxNtolitter = 0.0
  casaflux%FluxPtolitter = 0.0

  mimicsbiome%ligninNratioAvg(:) = 0.0
  casapool%dClitterdt(:,:) = 0.0

  cleaf2met = 0.0
  cleaf2str = 0.0
  croot2met = 0.0
  croot2str = 0.0
  cwood2cwd = 0.0
  cwd2co2 = 0.0
  cwd2str = 0.0

  nleaf2met = 0.0
  nleaf2str = 0.0
  nroot2met = 0.0
  nroot2str = 0.0
  nwood2cwd = 0.0

  pleaf2met = 0.0
  pleaf2str = 0.0
  proot2met = 0.0
  proot2str = 0.0
  pwood2cwd = 0.0

  DO npt=1,mp
      IF (casamet%iveg2(npt) /= icewater) THEN

          casapool%dCplantdt(npt,:)= casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)   &
                                     - casaflux%kplant(npt,:) * casapool%cplant(npt,:)
   
          ! Calculate fraction C to labile pool as a fraction of GPP, not NPP
          casapool%dClabiledt(npt) = casaflux%Cgpp(npt) * casaflux%fracClabile(npt) - casaflux%clabloss(npt)
   
          !Compute litter inputs from casa (gC/m2/day)
          cleaf2met(npt) = casaflux%fromPtoL(npt,metb,leaf)  * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          cleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf)   * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
          croot2str(npt) = casaflux%fromPtoL(npt,str,froot)  * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
          cwood2cwd(npt) = casaflux%fromPtoL(npt,cwd,wood)   * casaflux%kplant(npt,wood)  * casapool%cplant(npt,wood)

          if ((cleaf2met(npt) + cleaf2str(npt)) .gt. 0.0 ) then
              fracLeafMetabolic = cleaf2met(npt) / (cleaf2met(npt) + cleaf2str(npt))
          else
              fracLeafMetabolic = 0.0
          endif

          !For MIMICS, litter inputs are converted from gC/m2/day to mgC/cm3/day  

          !! Will suggested that CWD does not go into structural material all at once (4/20/2015).
          !! The decay rate of CWD and fraction lost to Rh are computed in mimics_xratesoil (-mdh 4/26/2015)
          !!   casaflux%klitter(npt,cwd) - decay rate of CWD
          !!   fromLtoCO2(mp,mlitter)    - fraction of decomposed litter emitted as CO2 (fraction)

          cwd2co2(npt) =  casaflux%fromLtoCO2(npt,cwd) * casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)
          cwd2str(npt) = (1.0-casaflux%fromLtoCO2(npt,cwd)) * casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)
          ! Calculate change in CWD pool.  CWD pool is updated in mimics_ccyle (-mdh 5/4/2015)
          casapool%dClitterdt(npt,cwd) = cwood2cwd(npt) - casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)

          ! mimicsflux%ClitInput(*,*) is both a flux to the MIMICS litter poosl and an output variable
          mimicsflux%ClitInput(npt,metbc) = ( cleaf2met(npt) &
                                            + croot2met(npt)) * 0.1 / mimicsbiome%depth(veg%iveg(npt))

          mimicsflux%ClitInput(npt,struc) = ( cleaf2str(npt) &
                                            + croot2str(npt) &
                                            + cwd2str(npt)) * 0.1 / mimicsbiome%depth(veg%iveg(npt))

          ! For consistency, compute the Litter Inputs that are written to casa's netcdf file. -mdh 1/23/2017
          casaflux%ClitInptMet(npt) = cleaf2met(npt) + croot2met(npt)
          casaflux%ClitInptStruc(npt) = cleaf2str(npt) + croot2str(npt)+ cwd2str(npt)


          ! Weighted average of lignin:N ratio for all litter inputs. Used for mimicsbiome%fmet(npt) calculation below
          mimicsbiome%ligninNratioAvg(npt) =  &
                                ( mimicsbiome%ligninNratio(npt,leaf)  * (cleaf2met(npt) + cleaf2str(npt)) &
                                + mimicsbiome%ligninNratio(npt,froot) * (croot2met(npt) + croot2str(npt)) &
                                + mimicsbiome%ligninNratio(npt,wood)  * (cwd2str(npt)) ) &
                                / max(0.001,cleaf2met(npt)+cleaf2str(npt)+croot2met(npt)+croot2str(npt)+cwd2str(npt))

      ENDIF
  ENDDO
 

  DO npt = 1,mp 

      ! fmet = fraction of plant residue transferred to metabolic litter (0.0 - 1.0)
      ! It is computed from the weighted average lignin:N ratio of all plant pools 
      !   (g lignin / g N) =  (g lignin / g C) / (g N / g C) 

      ! fmet is too high, causing MICk to go to zero for many sites.
      ! Will suggested we create a global value that is 85% of the local value. -mdh 5/11/2015
      ! mimicsbiome%fmet(npt) =  0.85 - 0.013 * mimicsbiome%ligninNratioAvg(npt) 
      ! mimicsbiome%fmet(npt) =  0.85 *(0.85 - 0.013 * mimicsbiome%ligninNratioAvg(npt)) 
      mimicsbiome%fmet(npt) = mimicsbiome%fmet_p(1) &
                              * (mimicsbiome%fmet_p(2) - mimicsbiome%fmet_p(3)*mimicsbiome%ligninNratioAvg(npt))
      if(mimicsbiome%fmet(npt) < 0.0) then
          write(*,'(a28,i5,a2,f10.6)') '  WARNING: mimicsbiome%fmet(', npt, ')=', mimicsbiome%fmet(npt)
          write(*,*) '  Resetting fmet to zero...'
          mimicsbiome%fmet(npt) = 0.0
      endif

      ! fPHYS = fraction of tauR and tauK partitioned into SOMp (0.0 - 1.0)
      ! mimicsbiome%fPHYS(npt,1) = 0.3 * exp(1.3 * soil%clay(npt))
      ! mimicsbiome%fPHYS(npt,2) = 0.2 * exp(0.8 * soil%clay(npt))
      mimicsbiome%fPHYS(npt,1) = mimicsbiome%fPHYS_r(1) * exp(mimicsbiome%fPHYS_r(2) * soil%clay(npt))
      mimicsbiome%fPHYS(npt,2) = mimicsbiome%fPHYS_K(1) * exp(mimicsbiome%fPHYS_K(2) * soil%clay(npt))

      ! fCHEM = fraction of tauR and tauK partitioned into SOMc (0.0 - 1.0)
      ! mimicsbiome%fCHEM(npt,1) = 0.1 * exp(-3.0 * mimicsbiome%fmet(npt))
      ! mimicsbiome%fCHEM(npt,2) = 0.3 * exp(-3.0 * mimicsbiome%fmet(npt))
      ! Updated based on Will's (-mdh 6/1/2015)
      ! mimicsbiome%fCHEM(npt,1) = mimicsbiome%fCHEM(npt,1) * 4.0
      ! mimicsbiome%fCHEM(npt,2) = mimicsbiome%fCHEM(npt,2) * 4.0

      mimicsbiome%fCHEM(npt,1) = mimicsbiome%fCHEM_r(1) * exp(mimicsbiome%fCHEM_r(2) * mimicsbiome%fmet(npt)) &
                                 * mimicsbiome%fCHEM_r(3)
      mimicsbiome%fCHEM(npt,2) = mimicsbiome%fCHEM_K(1) * exp(mimicsbiome%fCHEM_K(2) * mimicsbiome%fmet(npt)) &
                                 * mimicsbiome%fCHEM_K(3)

      ! fAVAL = fraction of tauR and tauK partitioned into SOMa (0.0 - 1.0)
      mimicsbiome%fAVAL(npt,1) = 1.0 - mimicsbiome%fCHEM(npt,1) - mimicsbiome%fPHYS(npt,1)
      mimicsbiome%fAVAL(npt,2) = 1.0 - mimicsbiome%fCHEM(npt,2) - mimicsbiome%fPHYS(npt,2)
 
!     if (casamet%ijgcm(npt) .eq. 9401) then
!         write(*,*) 'mimicsbiome%ligninNratioAvg(',npt,') = ', mimicsbiome%ligninNratioAvg(npt)
!         write(*,*) 'mimicsbiome%fmet(',npt,') = ', mimicsbiome%fmet(npt)
!         write(*,*) 'fAVAL(1), fAVAL(2) = ', mimicsbiome%fAVAL(npt,1), mimicsbiome%fAVAL(npt,2)
!         write(*,*) 'fCHEM(1), fCHEM(2) = ', mimicsbiome%fCHEM(npt,1), mimicsbiome%fCHEM(npt,2)
!         write(*,*) 'fPHYS(1), fPHYS(2) = ', mimicsbiome%fPHYS(npt,1), mimicsbiome%fPHYS(npt,2)
!     endif
  ENDDO


END SUBROUTINE mimics_delplant


!--------------------------------------------------------------------------------
!  SUBROUTINE mimics_soil_forward_MM
!    - calculate daily changes in litter, microbial, and SOM pools

SUBROUTINE mimics_soil_forwardMM(mp,iYrCnt,idoy,mdaily,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)

  ! Function Arguments
  integer, INTENT(IN) :: mp,iYrCnt,idoy  ! number of grid points, simulation year count, day of year
  integer, INTENT(IN) :: mdaily
  real, dimension(mp),INTENT(IN) :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd

  ! Local Variables
  integer :: npt, ihr
  integer, parameter :: NHOURS = 24
  real(r_2) :: LITmin(4), MICtrn(6), SOMmin(2), DEsorb, OXIDAT
  real(r_2) :: dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk
  real(r_2) :: NHOURSf
  real(r_2) :: Tsoil	!! average soil temperature for the day (degrees C)
  real(r_2) :: theta_liq       ! WW average liquid soil water
  real(r_2) :: theta_frzn      ! WW average frozen soil water
  real(r_2) :: air_filled_porosity !Fraction of 1.0.  Different from 1.0-theta_liq since it includes ice
  real(r_2) :: fW              ! CORPSE moisture function theta_liq^3*(1-air_filled_porosity)^2.5,
                               ! WW adjusted to give max values of 1, min = 0.01
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps 
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)

! if (iptToSave_mimics > 0) then
!     open(214,file=sPtFileNameMIMICS, access='APPEND')
! endif

  NHOURSf = real(NHOURS)

  do npt=1,mp

  IF(casamet%iveg2(npt) /= icewater) THEN

      mimicsflux%Chresp(npt) = 0.0

      ! casaflux%CnppAn(npt) =  average annual NPP (gC/m2/yr)
      ! Restrict tauMod to values 0.8 - 1.2 (Will Wieder, 5/4/2015)
      ! mimicsbiome%tauMod(npt) = SQRT(casaflux%CnppAn(npt)/100)
      mimicsbiome%tauMod(npt) = SQRT(casaflux%CnppAn(npt)/mimicsbiome%tauModDenom)
      mimicsbiome%tauMod(npt) = MAX(mimicsbiome%tauMod_MIN, mimicsbiome%tauMod(npt))
      mimicsbiome%tauMod(npt) = MIN(mimicsbiome%tauMod_MAX, mimicsbiome%tauMod(npt))

      ! Use site-level value for tauR (-mdh 4/20/2015)
      !! mimicsbiome%tauR(npt) = 5.2 * 0.0001 * exp(0.3 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      !! Global value for tauR (-mdh 4/6/2015)
      !! mimicsbiome%tauR(npt) = 5.2 * 0.0001 * exp(0.4 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      ! Use site-level value for tauK (-mdh 4/20/2015)
      !! mimicsbiome%tauK(npt) = 2.4 * 0.0001 * exp(0.1 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
    
      ! Updated based on Will's (-mdh 6/1/2015).  I had been using site-level values, so only tauR changed.
      ! mimicsbiome%tauR(npt) = 5.2 * 0.0001 * exp(0.4 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      ! mimicsbiome%tauK(npt) = 2.4 * 0.0001 * exp(0.1 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      mimicsbiome%tauR(npt) = mimicsbiome%tau_r(1) * &
                                exp(mimicsbiome%tau_r(2) * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      mimicsbiome%tauK(npt) = mimicsbiome%tau_k(1) * &
                                exp(mimicsbiome%tau_k(2) * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)

      ! Vmax - temperature sensitive maximum reaction velocities (mg C (mg MIC)-1 h-1) 
      Tsoil = casamet%tsoilavg(npt) - tkzeroc

     ! Read in soil moisture data as in CORPSE
      theta_liq  = min(1.0, casamet%moistavg(npt)/soil%ssat(npt))     ! fraction of liquid water-filled pore space (0.0 - 1.0)
      theta_frzn = min(1.0, casamet%frznmoistavg(npt)/soil%ssat(npt)) ! fraction of frozen water-filled pore space (0.0 - 1.0)
      air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)

      if (mimicsbiome%fWFunction .eq. CORPSE) then
        ! CORPSE water scalar, adjusted to give maximum values of 1
        fW = (theta_liq**3 * air_filled_porosity**2.5)/0.022600567942709
        fW = max(0.05, fW) 
      elseif (mimicsbiome%fWFunction .eq. CASACNP) then
        ! CASA water scalar, does not use frozen water in the calculation!
        ! local variables
        fW = ((theta_liq-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe &
           * ((theta_liq-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd      
        fW = min(fW, 1.0)
        fW = max(0.01, fW)
      else
        fW = 1.0
      endif

      mimicspool%fW(npt) =  fW
      mimicspool%thetaLiq(npt)  =  theta_liq
      mimicspool%thetaFrzn(npt) =  theta_frzn

      mimicsbiome%Vmax(npt,R1) = exp(mimicsbiome%Vslope(R1) * Tsoil + mimicsbiome%Vint(R1)) &
                                 * mimicsbiome%av(R1) * mimicsbiome%Vmod(R1) *fW
      mimicsbiome%Vmax(npt,R2) = exp(mimicsbiome%Vslope(R2) * Tsoil + mimicsbiome%Vint(R2)) &
                                 * mimicsbiome%av(R2) * mimicsbiome%Vmod(R2) *fW
      mimicsbiome%Vmax(npt,R3) = exp(mimicsbiome%Vslope(R3) * Tsoil + mimicsbiome%Vint(R3)) &
                                 * mimicsbiome%av(R3) * mimicsbiome%Vmod(R3) *fW
      mimicsbiome%Vmax(npt,K1) = exp(mimicsbiome%Vslope(K1) * Tsoil + mimicsbiome%Vint(K1)) &
                                 * mimicsbiome%av(K1) * mimicsbiome%Vmod(K1) *fW
      mimicsbiome%Vmax(npt,K2) = exp(mimicsbiome%Vslope(K2) * Tsoil + mimicsbiome%Vint(K2)) &
                                 * mimicsbiome%av(K2) * mimicsbiome%Vmod(K2) *fW
      mimicsbiome%Vmax(npt,K3) = exp(mimicsbiome%Vslope(K3) * Tsoil + mimicsbiome%Vint(K3)) &
                                 * mimicsbiome%av(K3) * mimicsbiome%Vmod(K3) *fW
   
      ! WW also modify TAU as a function of soil moisture, so things don't
      ! colapse in frozen soils...
      mimicsbiome%tauR(npt) = mimicsbiome%tauR(npt) * fW
      mimicsbiome%tauK(npt) = mimicsbiome%tauK(npt) * fW
 
      ! Km - half saturation constants (temperature sensitive) (mg C cm-3)
      mimicsbiome%Km(npt,R1) = exp(mimicsbiome%Kslope(R1) * Tsoil + mimicsbiome%Kint(R1)) &
                               * mimicsbiome%ak(R1) / mimicsbiome%Kmod(npt,R1)
      mimicsbiome%Km(npt,R2) = exp(mimicsbiome%Kslope(R2) * Tsoil + mimicsbiome%Kint(R2)) &
                               * mimicsbiome%ak(R2) / mimicsbiome%Kmod(npt,R2)
      mimicsbiome%Km(npt,R3) = exp(mimicsbiome%Kslope(R3) * Tsoil + mimicsbiome%Kint(R3)) &
                               * mimicsbiome%ak(R3) / mimicsbiome%Kmod(npt,R3)
      mimicsbiome%Km(npt,K1) = exp(mimicsbiome%Kslope(K1) * Tsoil + mimicsbiome%Kint(K1)) &
                               * mimicsbiome%ak(K1) / mimicsbiome%Kmod(npt,K1)
      mimicsbiome%Km(npt,K2) = exp(mimicsbiome%Kslope(K2) * Tsoil + mimicsbiome%Kint(K2)) &
                               * mimicsbiome%ak(K2) / mimicsbiome%Kmod(npt,K2)
      mimicsbiome%Km(npt,K3) = exp(mimicsbiome%Kslope(K3) * Tsoil + mimicsbiome%Kint(K3)) &
                               * mimicsbiome%ak(K3) / mimicsbiome%Kmod(npt,K3)
    
      do ihr = 1, NHOURS

          ! Flows to and from MICr
  
          !MICr decomp of METABOLIC litter (f1. LITm-->MICr)
          LITmin(1) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R1) * mimicspool%LITm(npt) &
                      / (mimicsbiome%Km(npt,R1) + mimicspool%LITm(npt))
  
          !MICr decomp of STRUCTURAL litter (f2. LITs-->MICr)   
          LITmin(2) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%LITs(npt) &
                      / (mimicsbiome%Km(npt,R2) + mimicspool%LITs(npt))  
      
          !MICr turnover to SOMp (f41. MICr-->SOMp) 
          MICtrn(1) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fPHYS(npt,1)
      
          !MICr turnover to SOMc (f42. MICr-->SOMc)                  
          MICtrn(2) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fCHEM(npt,1)   
      
          !MICr turnover to SOMa (f43. MICr-->SOMa)               
          MICtrn(3) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fAVAL(npt,1)    
      
          !decomp of SOMa by MICr (f3. SOMa-->MICr)
          SOMmin(1) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R3) * mimicspool%SOMa(npt) &
                      / (mimicsbiome%Km(npt,R3) + mimicspool%SOMa(npt))   
      
          !Flows to and from MICk
      
          !decomp of METABOLIC litter (f5. LITm-->MICk)
          LITmin(3) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K1) * mimicspool%LITm(npt) &
                      / (mimicsbiome%Km(npt,K1) + mimicspool%LITm(npt)) 
      
          !decomp of STRUCTURAL litter (f6. LITs-->MICk)  
          LITmin(4) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%LITs(npt) &
                      / (mimicsbiome%Km(npt,K2) + mimicspool%LITs(npt))  
      
          !MICk turnover to SOMp (f81. MICk-->SOMp) 
          MICtrn(4) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fPHYS(npt,2)     
      
          !MICk turnover to SOMc (f82. MICk-->SOMc)             
          MICtrn(5) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fCHEM(npt,2) 
      
          !MICk turnover to SOMa (f83. MICk-->SOMa)
          MICtrn(6) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fAVAL(npt,2)       
      
          !decomp of SOMa by MICk (f7. SOMa-->MICk)         
          SOMmin(2) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K3) * mimicspool%SOMa(npt) &
                      / (mimicsbiome%Km(npt,K3) + mimicspool%SOMa(npt))   
      
          ! Desorbtion of SOMp to SOMa (function of fCLAY) (f9. SOMp-->SOMa)
          DEsorb = mimicspool%SOMp(npt) * mimicsbiome%desorb(npt)  
			    	       
          ! Oxidation of SOMc to SOMa (f10. SOMc-->SOMa)
          OXIDAT = (mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%SOMc(npt)    &
                    / (mimicsbiome%KO(2)*mimicsbiome%Km(npt,K2) + mimicspool%SOMc(npt)))      &
                    + (mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%SOMc(npt) &
                    / (mimicsbiome%KO(1)*mimicsbiome%Km(npt,R2) + mimicspool%SOMc(npt)))
    
          !Divide total litter inputs (mgC/cm3/day) by number of hours in the daily timestep (NHOURSf)
          dLITm = mimicsflux%ClitInput(npt,metbc)/NHOURSf * (1.0-mimicsbiome%Fi(metbc)) - LITmin(1) - LITmin(3)
          dMICr = mimicsbiome%MGE(1)*(LITmin(1)+ SOMmin(1)) + mimicsbiome%MGE(2)*(LITmin(2)) - sum(MICtrn(1:3))
          dSOMp = mimicsflux%ClitInput(npt,metbc)/NHOURSf * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4)- DEsorb 
          dLITs = mimicsflux%ClitInput(npt,struc)/NHOURSf * (1.0-mimicsbiome%Fi(struc)) - LITmin(2) - LITmin(4)
          dMICk = mimicsbiome%MGE(3)*(LITmin(3)+ SOMmin(2)) + mimicsbiome%MGE(4)*(LITmin(4)) - sum(MICtrn(4:6)) 
          dSOMc = mimicsflux%ClitInput(npt,struc)/NHOURSf * mimicsbiome%Fi(struc) + MICtrn(2) + MICtrn(5) - OXIDAT
          dSOMa = MICtrn(3) + MICtrn(6) + DEsorb + OXIDAT - SOMmin(1) - SOMmin(2)
      
          !Sum daily heterotrphic respiration flux (mgC/cm3) 
          mimicsflux%Chresp(npt) = mimicsflux%Chresp(npt)  &
                                   + (1.0 - mimicsbiome%MGE(1)) * (LITmin(1) + SOMmin(1)) &
                                   + (1.0 - mimicsbiome%MGE(2)) * (LITmin(2))             &
                                   + (1.0 - mimicsbiome%MGE(3)) * (LITmin(3) + SOMmin(2)) &
                                   + (1.0 - mimicsbiome%MGE(4)) * (LITmin(4))
    
          !Update pools (mgC/cm3)
          mimicspool%LITm(npt) = mimicspool%LITm(npt) + dLITm
          mimicspool%LITs(npt) = mimicspool%LITs(npt) + dLITs
          mimicspool%MICr(npt) = mimicspool%MICr(npt) + dMICr
          mimicspool%MICk(npt) = mimicspool%MICk(npt) + dMICk
          mimicspool%SOMa(npt) = mimicspool%SOMa(npt) + dSOMa
          mimicspool%SOMc(npt) = mimicspool%SOMc(npt) + dSOMc
          mimicspool%SOMp(npt) = mimicspool%SOMp(npt) + dSOMp

      end do

!     if (casamet%ijgcm(npt) .eq. 2945) then     ! Crop (12)
!     if (casamet%ijgcm(npt) .eq. 10354) then    ! Tundra (18)
!     if (casamet%ijgcm(npt) .eq. 10052) then    ! Grassland (10)
!     if (casamet%ijgcm(npt) .eq. 10919) then    ! Deciduous broadleaf (4)
!     if (casamet%ijgcm(npt) .eq. 11018) then    ! Evergreen needleleaf (1)
!     if (casamet%ijgcm(npt) .eq. 11569) then    ! Evergreen broadleaf (2)

      if (casamet%ijgcm(npt) .eq. iptToSave_mimics) then 
          ! Write to point file daily if mdaily==1, else write once a year on day 365. -mdh 5/14/2018
          if ((mdaily == 1) .or. (idoy==365)) then
              call WritePointMIMICS(214, sPtFileNameMIMICS, npt, mp, iYrCnt, idoy, &
                  cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
                  LITmin, MICtrn, SOMmin, DEsorb, OXIDAT, &
                  dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil)
          endif
      endif 

  ENDIF

  end do

END SUBROUTINE mimics_soil_forwardMM

!--------------------------------------------------------------------------------
!  SUBROUTINE mimics_soil_reverseMM
!    - calculate daily changes in litter, microbial, and SOM pools
!    - alternative to mimics_soil_forwardMM; implements the reverse Michaelis-Mention kinetics

SUBROUTINE mimics_soil_reverseMM(mp,iYrCnt,idoy,mdaily,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)

  ! Function Arguments
  integer, INTENT(IN) :: mp,iYrCnt,idoy  ! number of grid points, simulation year count, day of year
  integer, INTENT(IN) :: mdaily
  real, dimension(mp),INTENT(IN) :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd

  ! Local Variables
  integer :: npt, ihr
  integer, parameter :: NHOURS = 24
  real(r_2) :: LITmin(4), MICtrn(6), SOMmin(2), DEsorb, OXIDAT
  real(r_2) :: dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk
  real(r_2) :: NHOURSf
  real(r_2) :: Tsoil	!! average soil temperature for the day (degrees C)
  real(r_2) :: theta_liq       ! WW average liquid soil water
  real(r_2) :: theta_frzn      ! WW average frozen soil water
  real(r_2) :: air_filled_porosity !Fraction of 1.0.  Different from 1.0-theta_liq since it includes ice
  real(r_2) :: fW              ! CORPSE moisture function theta_liq^3*(1-air_filled_porosity)^2.5,
                               ! WW adjusted to give max values of 1, min =
                               ! 0.01
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps 
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)

! if (iptToSave_mimics > 0) then
!     open(214,file=sPtFileNameMIMICS, access='APPEND')
! endif

  NHOURSf = real(NHOURS)

  do npt=1,mp

  IF(casamet%iveg2(npt) /= icewater) THEN

      mimicsflux%Chresp(npt) = 0.0

      ! casaflux%CnppAn(npt) =  average annual NPP (gC/m2/yr)
      ! Restrict tauMod to values 0.8 - 1.2 (Will Wieder, 5/4/2015)
      ! mimicsbiome%tauMod(npt) = SQRT(casaflux%CnppAn(npt)/100)
      mimicsbiome%tauMod(npt) = SQRT(casaflux%CnppAn(npt)/mimicsbiome%tauModDenom)
      mimicsbiome%tauMod(npt) = MAX(mimicsbiome%tauMod_MIN, mimicsbiome%tauMod(npt))
      mimicsbiome%tauMod(npt) = MIN(mimicsbiome%tauMod_MAX, mimicsbiome%tauMod(npt))

      ! Use site-level value for tauR (-mdh 4/20/2015)
      !! mimicsbiome%tauR(npt) = 5.2 * 0.0001 * exp(0.3 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      !! Global value for tauR (-mdh 4/6/2015)
      !! mimicsbiome%tauR(npt) = 5.2 * 0.0001 * exp(0.4 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      ! Use site-level value for tauK (-mdh 4/20/2015)
      !! mimicsbiome%tauK(npt) = 2.4 * 0.0001 * exp(0.1 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
    
      ! Updated based on Will's (-mdh 6/1/2015).  I had been using site-level values, so only tauR changed.
      ! mimicsbiome%tauR(npt) = 5.2 * 0.0001 * exp(0.4 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      ! mimicsbiome%tauK(npt) = 2.4 * 0.0001 * exp(0.1 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      mimicsbiome%tauR(npt) = mimicsbiome%tau_r(1) * &
                                exp(mimicsbiome%tau_r(2) * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      mimicsbiome%tauK(npt) = mimicsbiome%tau_k(1) * &
                                exp(mimicsbiome%tau_k(2) * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)

      ! Vmax - temperature sensitive maximum reaction velocities (mg C (mg MIC)-1 h-1) 
      Tsoil = casamet%tsoilavg(npt) - tkzeroc
      ! Read in soil moisture data as in CORPSE
      theta_liq  = min(1.0, casamet%moistavg(npt)/soil%ssat(npt))     ! fraction of liquid water-filled pore space (0.0 - 1.0)
      theta_frzn = min(1.0, casamet%frznmoistavg(npt)/soil%ssat(npt)) ! fraction of frozen water-filled pore space (0.0 - 1.0)
      air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)

      if (mimicsbiome%fWFunction .eq. CORPSE) then
        ! CORPSE water scalar, adjusted to give maximum values of 1
        fW = (theta_liq**3 * air_filled_porosity**2.5)/0.022600567942709
        fW = max(0.05, fW) 
      elseif (mimicsbiome%fWFunction .eq. CASACNP) then
        ! CASA water scalar, does not use frozen water in the calculation!
        ! local variables
        fW = ((theta_liq-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe &
           * ((theta_liq-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd      
        fW = min(fW, 1.0)
        fW = max(0.01, fW)
      else
        fW = 1.0
      endif

      mimicspool%fW(npt) =  fW
      mimicspool%thetaLiq(npt)  =  theta_liq
      mimicspool%thetaFrzn(npt) =  theta_frzn

      mimicsbiome%Vmax(npt,R1) = exp(mimicsbiome%Vslope(R1) * Tsoil + mimicsbiome%Vint(R1)) &
                                 * mimicsbiome%av(R1) * mimicsbiome%Vmod(R1) * fW
      mimicsbiome%Vmax(npt,R2) = exp(mimicsbiome%Vslope(R2) * Tsoil + mimicsbiome%Vint(R2)) &
                                 * mimicsbiome%av(R2) * mimicsbiome%Vmod(R2) * fW
      mimicsbiome%Vmax(npt,R3) = exp(mimicsbiome%Vslope(R3) * Tsoil + mimicsbiome%Vint(R3)) &
                                 * mimicsbiome%av(R3) * mimicsbiome%Vmod(R3) * fW
      mimicsbiome%Vmax(npt,K1) = exp(mimicsbiome%Vslope(K1) * Tsoil + mimicsbiome%Vint(K1)) &
                                 * mimicsbiome%av(K1) * mimicsbiome%Vmod(K1) * fW
      mimicsbiome%Vmax(npt,K2) = exp(mimicsbiome%Vslope(K2) * Tsoil + mimicsbiome%Vint(K2)) &
                                 * mimicsbiome%av(K2) * mimicsbiome%Vmod(K2) * fW
      mimicsbiome%Vmax(npt,K3) = exp(mimicsbiome%Vslope(K3) * Tsoil + mimicsbiome%Vint(K3)) &
                                 * mimicsbiome%av(K3) * mimicsbiome%Vmod(K3) * fW
    
      ! WW also modify TAU as a function of soil moisture, so things don't
      ! colapse in frozen soils...
      mimicsbiome%tauR(npt) = mimicsbiome%tauR(npt) * fW
      mimicsbiome%tauK(npt) = mimicsbiome%tauK(npt) * fW
 
      ! Km - half saturation constants (temperature sensitive) (mg C cm-3)
      mimicsbiome%Km(npt,R1) = exp(mimicsbiome%Kslope(R1) * Tsoil + mimicsbiome%Kint(R1)) &
                               * mimicsbiome%ak(R1) / mimicsbiome%Kmod(npt,R1)
      mimicsbiome%Km(npt,R2) = exp(mimicsbiome%Kslope(R2) * Tsoil + mimicsbiome%Kint(R2)) &
                               * mimicsbiome%ak(R2) / mimicsbiome%Kmod(npt,R2)
      mimicsbiome%Km(npt,R3) = exp(mimicsbiome%Kslope(R3) * Tsoil + mimicsbiome%Kint(R3)) &
                               * mimicsbiome%ak(R3) / mimicsbiome%Kmod(npt,R3)
      mimicsbiome%Km(npt,K1) = exp(mimicsbiome%Kslope(K1) * Tsoil + mimicsbiome%Kint(K1)) &
                               * mimicsbiome%ak(K1) / mimicsbiome%Kmod(npt,K1)
      mimicsbiome%Km(npt,K2) = exp(mimicsbiome%Kslope(K2) * Tsoil + mimicsbiome%Kint(K2)) &
                               * mimicsbiome%ak(K2) / mimicsbiome%Kmod(npt,K2)
      mimicsbiome%Km(npt,K3) = exp(mimicsbiome%Kslope(K3) * Tsoil + mimicsbiome%Kint(K3)) &
                               * mimicsbiome%ak(K3) / mimicsbiome%Kmod(npt,K3)
    
      do ihr = 1, NHOURS

          ! Flows to and from MICr
  
          !MICr decomp of METABOLIC litter (f1. LITm-->MICr), reverse Michaelis-Menton Kinetics
          LITmin(1) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R1) * mimicspool%LITm(npt) &
                      / (mimicsbiome%Km(npt,R1) + mimicspool%MICr(npt))
  
          !MICr decomp of STRUCTURAL litter (f2. LITs-->MICr), reverse Michaelis-Menton Kinetics   
          LITmin(2) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%LITs(npt) &
                      / (mimicsbiome%Km(npt,R2) + mimicspool%MICr(npt))  
      
          !MICr turnover to SOMp (f41. MICr-->SOMp) 
          MICtrn(1) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fPHYS(npt,1)
      
          !MICr turnover to SOMc (f42. MICr-->SOMc)                  
          MICtrn(2) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fCHEM(npt,1)   
      
          !MICr turnover to SOMa (f43. MICr-->SOMa)               
          MICtrn(3) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fAVAL(npt,1)    
      
          !decomp of SOMa by MICr (f3. SOMa-->MICr), reverse Michaelis-Menton Kinetics
          SOMmin(1) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R3) * mimicspool%SOMa(npt) &
                      / (mimicsbiome%Km(npt,R3) + mimicspool%MICr(npt))   
      
          !Flows to and from MICk
      
          !decomp of METABOLIC litter (f5. LITm-->MICk), reverse Michaelis-Menton Kinetics
          LITmin(3) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K1) * mimicspool%LITm(npt) &
                      / (mimicsbiome%Km(npt,K1) + mimicspool%MICk(npt)) 
      
          !decomp of STRUCTURAL litter (f6. LITs-->MICk), reverse Michaelis-Menton Kinetics 
          LITmin(4) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%LITs(npt) &
                      / (mimicsbiome%Km(npt,K2) + mimicspool%MICk(npt))  
      
          !MICk turnover to SOMp (f81. MICk-->SOMp) 
          MICtrn(4) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fPHYS(npt,2)     
      
          !MICk turnover to SOMc (f82. MICk-->SOMc)             
          MICtrn(5) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fCHEM(npt,2) 
      
          !MICk turnover to SOMa (f83. MICk-->SOMa)
          MICtrn(6) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fAVAL(npt,2)       
      
          !decomp of SOMa by MICk (f7. SOMa-->MICk, reverse Michaelis-Menton Kinetics         
          SOMmin(2) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K3) * mimicspool%SOMa(npt) &
                      / (mimicsbiome%Km(npt,K3) + mimicspool%MICk(npt))   
      
          ! Desorbtion of SOMp to SOMa (function of fCLAY) (f9. SOMp-->SOMa)
          DEsorb = mimicspool%SOMp(npt) * mimicsbiome%desorb(npt)  
			    	       
          ! Oxidation of SOMc to SOMa (f10. SOMc-->SOMa), reverse Michaelis-Menton Kinetics
          OXIDAT = (mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%SOMc(npt)    &
                    / (mimicsbiome%KO(2)*mimicsbiome%Km(npt,K2) + mimicspool%MICk(npt)))      &
                    + (mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%SOMc(npt) &
                    / (mimicsbiome%KO(1)*mimicsbiome%Km(npt,R2) + mimicspool%MICr(npt)))
    
          !Divide total litter inputs (mgC/cm3/day) by number of hours in the daily timestep (NHOURSf)
          dLITm = mimicsflux%ClitInput(npt,metbc)/NHOURSf * (1.0-mimicsbiome%Fi(metbc)) - LITmin(1) - LITmin(3)
          dMICr = mimicsbiome%MGE(1)*(LITmin(1)+ SOMmin(1)) + mimicsbiome%MGE(2)*(LITmin(2)) - sum(MICtrn(1:3))
          dSOMp = mimicsflux%ClitInput(npt,metbc)/NHOURSf * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4)- DEsorb 
          dLITs = mimicsflux%ClitInput(npt,struc)/NHOURSf * (1.0-mimicsbiome%Fi(struc)) - LITmin(2) - LITmin(4)
          dMICk = mimicsbiome%MGE(3)*(LITmin(3)+ SOMmin(2)) + mimicsbiome%MGE(4)*(LITmin(4)) - sum(MICtrn(4:6)) 
          dSOMc = mimicsflux%ClitInput(npt,struc)/NHOURSf * mimicsbiome%Fi(struc) + MICtrn(2) + MICtrn(5) - OXIDAT
          dSOMa = MICtrn(3) + MICtrn(6) + DEsorb + OXIDAT - SOMmin(1) - SOMmin(2)
      
          !Sum daily heterotrphic respiration flux (mgC/cm3) 
          mimicsflux%Chresp(npt) = mimicsflux%Chresp(npt)  &
                                   + (1.0 - mimicsbiome%MGE(1)) * (LITmin(1) + SOMmin(1)) &
                                   + (1.0 - mimicsbiome%MGE(2)) * (LITmin(2))             &
                                   + (1.0 - mimicsbiome%MGE(3)) * (LITmin(3) + SOMmin(2)) &
                                   + (1.0 - mimicsbiome%MGE(4)) * (LITmin(4))
    
          !Update pools (mgC/cm3)
          mimicspool%LITm(npt) = mimicspool%LITm(npt) + dLITm
          mimicspool%LITs(npt) = mimicspool%LITs(npt) + dLITs
          mimicspool%MICr(npt) = mimicspool%MICr(npt) + dMICr
          mimicspool%MICk(npt) = mimicspool%MICk(npt) + dMICk
          mimicspool%SOMa(npt) = mimicspool%SOMa(npt) + dSOMa
          mimicspool%SOMc(npt) = mimicspool%SOMc(npt) + dSOMc
          mimicspool%SOMp(npt) = mimicspool%SOMp(npt) + dSOMp

      end do

!     if (casamet%ijgcm(npt) .eq. 2945) then     ! Crop (12)
!     if (casamet%ijgcm(npt) .eq. 10354) then    ! Tundra (18)
!     if (casamet%ijgcm(npt) .eq. 10052) then    ! Grassland (10)
!     if (casamet%ijgcm(npt) .eq. 10919) then    ! Deciduous broadleaf (4)
!     if (casamet%ijgcm(npt) .eq. 11018) then    ! Evergreen needleleaf (1)
!     if (casamet%ijgcm(npt) .eq. 11569) then    ! Evergreen broadleaf (2)

      if (casamet%ijgcm(npt) .eq. iptToSave_mimics) then 
          ! Write to point file daily if mdaily==1, else write once a year on day 365. -mdh 5/14/2018
          if ((mdaily == 1) .or. (idoy==365)) then
              call WritePointMIMICS(214, sPtFileNameMIMICS, npt, mp, iYrCnt, idoy, &
                  cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
                  LITmin, MICtrn, SOMmin, DEsorb, OXIDAT, &
                  dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil)
          endif 
      endif 

  ENDIF

  end do

END SUBROUTINE mimics_soil_reverseMM

!--------------------------------------------------------------------------------

SUBROUTINE mimics_ccycle(veg,casabiome,casapool,casaflux,casamet)

! Called daily
! Update all C pool sizes that were not updated in mimics_soil_*MM
! This includes plant C pools, but not litter C or soil C pools

  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! Local variables
  INTEGER i,j,k,npt,nland

  DO npt=1,mp
      IF(casamet%iveg2(npt) == icewater) THEN
          casamet%glai(npt) = 0.0
      ELSE  
          casapool%cplant(npt,:) = casapool%cplant(npt,:) + casapool%dCplantdt(npt,:) * deltpool 
          casapool%clabile(npt)  = casapool%clabile(npt)  + casapool%dclabiledt(npt)  * deltpool 
          casamet%glai(npt)      = MAX(casabiome%glaimin(veg%iveg(npt)), &
                                       casabiome%sla(veg%iveg(npt)) * casapool%cplant(npt,leaf))
          casamet%glai(npt)      = MIN(casabiome%glaimax(veg%iveg(npt)), casamet%glai(npt))
      
          ! Update CWD pool (-mdh 5/4/2015)
          casapool%clitter(npt,:) = casapool%clitter(npt,:) + casapool%dClitterdt(npt,:) * deltpool 

          DO i=1,mplant
              IF(casapool%cplant(npt,i) < 0.0) THEN
!                 WRITE(57,*) 'Cpool: npt,ivt',npt,casamet%iveg2(npt),casapool%cplant(npt,:)
                  !call casa_poolzero(npt,1,casapool)
                  casapool%cplant(npt,i) = max(0.0, casapool%cplant(npt,i))
              ENDIF
          ENDDO
      
      ENDIF
  ENDDO !end of "npt"

END SUBROUTINE mimics_ccycle

!--------------------------------------------------------------------------------

!! Subroutine mimics_caccum accumulates daily C pool and flux values
!! Average annual values will be computed at the end of the simulation 
!!   in subroutine mimics_poolfluxout

SUBROUTINE mimics_caccum(mp,cwd2co2)

  !Subroutine arguments
  integer, INTENT(IN) :: mp	! number of grid cells
  real, dimension(mp),INTENT(IN) :: cwd2co2	! respiration from decay of CWD (gC/m2/day)

  !Local variables
  integer   :: npt
  real(r_2) :: unitConv		! mgC/cm3 * depth(cm)* (1g/10^3mg)*(10^4cm2)/m2 = gC/m2

  do npt = 1,mp

  if (casamet%iveg2(npt) /= icewater) then

      unitConv = 10.0 * mimicsbiome%depth(veg%iveg(npt))
    
!     if (npt .eq. 1003) then
!     print *
!     print *, 'caccum: unitConv = ', unitConv
!     print *, 'caccum: mimicspoolAn%ClitterAv(',npt,',METBC) = ', mimicspoolAn%ClitterAn(npt,METBC), mimicspool%LITm(npt)
!     print *, 'caccum: mimicspoolAn%ClitterAn(',npt,',STRUC) = ', mimicspoolAn%ClitterAn(npt,STRUC), mimicspool%LITs(npt)
!
!     print *, 'caccum: mimicspoolAn%CmicrobeAn(',npt,',RSEL) = ', mimicspoolAn%CmicrobeAn(npt,RSEL), mimicspool%MICr(npt)
!     print *, 'caccum: mimicspoolAn%CmicrobeAn(',npt,',KSEL) = ', mimicspoolAn%CmicrobeAn(npt,KSEL), mimicspool%MICk(npt)
!
!     print *, 'caccum: mimicspoolAn%CsoilAn(',npt,',AVAL) = ', mimicspoolAn%CsoilAn(npt,AVAL), mimicspool%SOMa(npt)
!     print *, 'caccum: mimicspoolAn%CsoilAn(',npt,',CHEM) = ', mimicspoolAn%CsoilAn(npt,CHEM), mimicspool%SOMc(npt)
!     print *, 'caccum: mimicspoolAn%CsoilAn(',npt,',PHYS) = ', mimicspoolAn%CsoilAn(npt,PHYS), mimicspool%SOMp(npt)
!     endif

      mimicspoolAn%ClitterAn(npt,METBC) = mimicspoolAn%ClitterAn(npt,METBC) + mimicspool%LITm(npt) * unitConv
      mimicspoolAn%ClitterAn(npt,STRUC) = mimicspoolAn%ClitterAn(npt,STRUC) + mimicspool%LITs(npt) * unitConv
    
      mimicspoolAn%CmicrobeAn(npt,RSEL) = mimicspoolAn%CmicrobeAn(npt,RSEL) + mimicspool%MICr(npt) * unitConv
      mimicspoolAn%CmicrobeAn(npt,KSEL) = mimicspoolAn%CmicrobeAn(npt,KSEL) + mimicspool%MICk(npt) * unitConv
    
      mimicspoolAn%CsoilAn(npt,AVAL) = mimicspoolAn%CsoilAn(npt,AVAL) + mimicspool%SOMa(npt) * unitConv
      mimicspoolAn%CsoilAn(npt,CHEM) = mimicspoolAn%CsoilAn(npt,CHEM) + mimicspool%SOMc(npt) * unitConv
      mimicspoolAn%CsoilAn(npt,PHYS) = mimicspoolAn%CsoilAn(npt,PHYS) + mimicspool%SOMp(npt) * unitConv
    
      !Litter Inputs and Heterotrophic Respiration (gC/m2)
      mimicsfluxAn%CLitInputAn(npt,METBC) = mimicsfluxAn%CLitInputAn(npt,METBC) + mimicsflux%CLitInput(npt,METBC) * unitConv 
      mimicsfluxAn%CLitInputAn(npt,STRUC) = mimicsfluxAn%CLitInputAn(npt,STRUC) + mimicsflux%CLitInput(npt,STRUC) * unitConv

      ! Include calculation for cwd2co2 in daily as well as annual CO2 flux (-mdh 4/26/2015, 1/23/2017)
      ! mimicsflux%Chresp is in units of mg/c3 whereas mimicsfluxAn%ChrespAn is gC/m2
      mimicsflux%Chresp(npt) = mimicsflux%Chresp(npt) + cwd2co2(npt)/unitConv 
      mimicsfluxAn%ChrespAn(npt) = mimicsfluxAn%ChrespAn(npt) + mimicsflux%Chresp(npt) * unitConv 

      ! Added f(T), f(W), thetaLiq, thetaFrzn to MIMICS output. -mdh 11/27/2017
      mimicspoolAn%fTAn(npt) = mimicspoolAn%fTAn(npt) + mimicspool%fT(npt) 
      mimicspoolAn%fWAn(npt) = mimicspoolAn%fWAn(npt) + mimicspool%fW(npt) 
      mimicspoolAn%thetaLiqAn(npt) = mimicspoolAn%thetaLiqAn(npt) + mimicspool%thetaLiq(npt) 
      mimicspoolAn%thetaFrznAn(npt) = mimicspoolAn%thetaFrznAn(npt) + mimicspool%thetaFrzn(npt) 

      if (cwd2co2(npt) .lt. 0.0 .or. mimicsflux%Chresp(npt) .lt. 0.0) then
          write(57,*) 'WARNING NEGATIVE RESPIRATION: ', &
                      'mimicsflux%Chresp(',npt,') = ', mimicsflux%Chresp(npt), &
                      'cwd2co2(',npt,') = ', cwd2co2(npt)
      endif
    
  endif
  end do

END SUBROUTINE mimics_caccum

!--------------------------------------------------------------------------------

END MODULE mimics_cycle_module
