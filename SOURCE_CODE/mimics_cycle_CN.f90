!--------------------------------------------------------------------------------
! FILE: mimics_cycle_CN.f90
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
!         - Implements the reverse Michaelis-Menten kinetics
!         - alternative to mimics_soil_forwardMM
!     SUBROUTINE mimics_soil_reverseMM_CN
!         - calculate daily changes in litter, microbial, and SOM pools
!         - replaces casa_coeffsoil and casa_delsoil when MIMICS-CN model is used 
!         - Implements the reverse Michaelis-Menten kinetics
!         - alternative to mimics_soil_forwardMM (currently no CN version fo forwardMM)
!     SUBROUTINE mimics_coeffplant
!         - calculate litter inputs to MIMICS model
!         - replaces casa_coeffplant when MIMICS model is used
!     SUBROUTINE mimics_delplant
!         - replaces casa_delplant when MIMICS model is used with icycle==1
!     SUBROUTINE mimics_delplant_CN
!         - replaces casa_delplant when MIMICS-CN model is used with icycle==2
!     SUBROUTINE mimics_xratesoil (-mdh 4/26/2015) 
!         - replaces casa_xratesoil and casa_coeffsoil when MIMICS or CORPSE models are used
!         - used to compute the dcomposition rate for CWD, casaflux%klitter(npt,cwd)
!     SUBROUTINE mimics_ccycle
!         - called daily
!         - update all C pools that were not updated in mimics_soil_forwardMM or mimics_soil_reverseMM 
!         - replaces casa_cnpcyle when MIMICS or CORPSE models are used with icycle==1
!     SUBROUTINE mimics_cncycle
!         - called daily
!         - update all C&N pools that were not updated in mimics_soil_forwardMM or mimics_soil_reverseMM 
!         - replaces casa_cnpcyle when MIMICS-CN is used with icycle==2
!     SUBROUTINE mimics_caccum
!         - called daily
!         - accumulate annual C fluxes and pool values for MIMICS
!     SUBROUTINE mimics_naccum
!         - called daily
!         - accumulate annual N fluxes and pool values for MIMICS-CN when icycle==2
!
! Contact: Melannie Hartman
!          melannie@ucar.edu
!
! History:
!   1/12/2015 - Created
!   8/24/2015 - Added reverse Michaelis-Menton option (subroutine mimics_soil_reverseMM)
!   4/29/2019 - Added MIMICS-CN code reverse Michaelis-Menton
!     SUBROUTINE mimics_soil_reverseMM_CN
!   1/13/2020 - Added root exudate flux to mimics_delplant (C only) and 
!               mimics_delplant_CN
!   9/27/2020 - Updated Nleaching calculation
!
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
  REAL(r_2), DIMENSION(mp,mplant) :: ratioLignintoN
  INTEGER :: npt    
  REAL(r_2), DIMENSION(mplant) :: CNplantMax

  casaflux%fromPtoL(:,:,:)      = 0.0
  casaflux%kplant(:,:)          = 0.0   

  ! When simulating N, calculate dynamic lignin:N ratio and reset mimicsbiome%ligninNratio
  ! which was initially calculated in mimics_readbiome. -mdh 1/20/2020

  ! ATTENTION: troubleshooting here becasue the lignin:N can be outrageoulsy large. 
  ! It happens when there is a lot of C but little plant N. 
  ! Restrict ratioLignintoN using maximum C:N ratio for the plant part. -mdh 2/3/2020
 
  if (icycle > 1) then
      do npt = 1,mp 
          if (casamet%iveg2(npt) /= icewater) then
              ! Dynamic ratioLignintoN calculation from casa_coeffplant. 
              ! Reset mimicsbiome%ligninNratio (initially calculated in readbiome) to ratioLignintoN. -mdh 1/20/2020
    
!             ! using max function to avoid dividing by zero, ypw 14/may/2008
!             ratioLignintoN(npt,leaf) = (casapool%Cplant(npt,leaf) &
!                                  /(max(1.0e-10,casapool%Nplant(npt,leaf)) *casabiome%ftransNPtoL(veg%iveg(npt),leaf))) &
!                                  * casabiome%fracLigninplant(veg%iveg(npt),leaf)  
!             ratioLignintoN(npt,froot)= (casapool%Cplant(npt,froot)&
!                                      /(max(1.0e-10,casapool%Nplant(npt,froot))*casabiome%ftransNPtoL(veg%iveg(npt),froot))) &
!                                  * casabiome%fracLigninplant(veg%iveg(npt),froot) 
    
              CNplantMax(leaf) = 1.0 / casabiome%ratioNCplantmin(veg%iveg(npt),leaf)
              CNplantMax(froot) = 1.0 / casabiome%ratioNCplantmin(veg%iveg(npt),froot)
              CNplantMax(wood) = 1.0 / casabiome%ratioNCplantmin(veg%iveg(npt),wood)
    
              if (casapool%Cplant(npt,leaf) / max(1.0e-10,casapool%Nplant(npt,leaf)) > CNplantMax(leaf)) then
                  ratioLignintoN(npt,leaf) = CNplantMax(leaf) / casabiome%ftransNPtoL(veg%iveg(npt),leaf) &
                                       * casabiome%fracLigninplant(veg%iveg(npt),leaf)  
              else
                  ratioLignintoN(npt,leaf) = (casapool%Cplant(npt,leaf) &
                                       /(max(1.0e-10,casapool%Nplant(npt,leaf)) *casabiome%ftransNPtoL(veg%iveg(npt),leaf))) &
                                       * casabiome%fracLigninplant(veg%iveg(npt),leaf)  
              endif
    
              if (casapool%Cplant(npt,froot) / max(1.0e-10,casapool%Nplant(npt,froot)) > CNplantMax(froot)) then
                  ratioLignintoN(npt,froot) = CNplantMax(froot) / casabiome%ftransNPtoL(veg%iveg(npt),froot) &
                                       * casabiome%fracLigninplant(veg%iveg(npt),froot)  
              else
                  ratioLignintoN(npt,froot) = (casapool%Cplant(npt,froot) &
                                       /(max(1.0e-10,casapool%Nplant(npt,froot)) *casabiome%ftransNPtoL(veg%iveg(npt),froot))) &
                                       * casabiome%fracLigninplant(veg%iveg(npt),froot)  
              endif
    
              mimicsbiome%ligninNratio(npt,leaf) = ratioLignintoN(npt,leaf)
              mimicsbiome%ligninNratio(npt,froot) = ratioLignintoN(npt,froot)
          end if

!         if (npt == 19) then
!             write(*,*)
!             write(*,*) 'ratioLignintoN(npt,leaf) =', ratioLignintoN(npt,leaf)
!             write(*,*) 'ratioLignintoN(npt,froot) =', ratioLignintoN(npt,froot)
!             write(*,*) 'froot C:N =', casapool%Cplant(npt,froot) / max(1.0e-10,casapool%Nplant(npt,froot))
!             write(*,*) 'casabiome%ftransNPtoL(veg%iveg(npt),froot) =', casabiome%ftransNPtoL(veg%iveg(npt),froot)
!             write(*,*) 'casabiome%fracLigninplant(veg%iveg(npt),froot) =', casabiome%fracLigninplant(veg%iveg(npt),froot)
!         endif
 
      end do
  endif

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
!     write(*,*)
!     write(*,*) 'ratioLignintoN leaf, froot, wood =', &
!                ratioLignintoN(npt,leaf), ratioLignintoN(npt,froot), ratioLignintoN(npt,wood)
!     write(*,*) 'mimicsbiome%ligninNratio leaf, froot, wood =', &
!                mimicsbiome%ligninNratio(npt,leaf), mimicsbiome%ligninNratio(npt,froot), mimicsbiome%ligninNratio(npt,wood)
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
!         write(*,*)
!         write(*,*) 'mimics_xratesoil:'
!         write(*,'(a42,f10.6)') 'casabiome%xkoptlitter(veg%iveg(npt)) =', casabiome%xkoptlitter(veg%iveg(npt))
!         write(*,'(a42,f10.6)') 'xktemp(npt) =', xktemp(npt)
!         write(*,'(a42,f10.6)') 'xkwater(npt) =', xkwater(npt)
!         write(*,'(a42,f10.6)') 'casabiome%litterrate(veg%iveg(npt),cwd) =', casabiome%litterrate(veg%iveg(npt),cwd)
!         write(*,'(a42,f10.6)') 'casaflux%klitter(npt,cwd) =', casaflux%klitter(npt,cwd)


      !! from casa_coeffsoil for reference
      !! casaflux%fromLtoS(:,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(:),wood)) ! CWD -> fmic
      !! casaflux%fromLtoS(:,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(:),wood)        ! CWD -> slow

      ! Fraction of cwd decomposition that goes to heterotrophic respiration
      casaflux%fromLtoCO2(npt,cwd) = 1.0 - 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(npt),wood)) &
                                         - 0.7 * casabiome%fracLigninplant(veg%iveg(npt),wood) 
      casaflux%fromLtoCO2(npt,cwd) = MAX(0.0, casaflux%fromLtoCO2(npt,cwd))
      casaflux%fromLtoCO2(npt,cwd) = MIN(1.0, casaflux%fromLtoCO2(npt,cwd))

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

  real(r_2), dimension(mp),INTENT(OUT) :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                                     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                                     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                                     cwd2co2,cwd2str
  ! Local Variables
  integer ::  npt
  real(r_2):: fracLeafMetabolic

  casaflux%FluxCtolitter = 0.0
  casaflux%FluxNtolitter = 0.0
  casaflux%FluxPtolitter = 0.0
  ! Added root exudate flux -mdh 1/13/2020
  casaflux%Cexudate = 0.0
  casaflux%Nexudate = 0.0
  casaflux%Pexudate = 0.0

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

          ! Compute root exudate C flux as a fraction of froot NPP. -mdh 1/13/2020
          casaflux%Cexudate(npt) = max(0.0, casabiome%fracRootExud(veg%iveg(npt)) * casaflux%Cnpp(npt) &
                                   * casaflux%fracCalloc(npt,froot))  
          casapool%dCplantdt(npt,froot) = casapool%dCplantdt(npt,froot) - casaflux%Cexudate(npt) 
   
          ! Calculate fraction C to labile pool as a fraction of GPP, not NPP
          casapool%dClabiledt(npt) = casaflux%Cgpp(npt) * casaflux%fracClabile(npt) - casaflux%clabloss(npt)
   
          !Compute litter inputs from casa (gC/m2/day)
          cleaf2met(npt) = casaflux%fromPtoL(npt,metb,leaf)  * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          cleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf)   * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          ! Add casaflux%Cexudate(npt) to metabolic litter. -mdh 1/13/2019
          !croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
          croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot) &
                           + casaflux%Cexudate(npt)
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

          ! mimicsflux%ClitInput(*,*) is both a flux to the MIMICS litter pools and an output variable
          ! Convert litter inputs from gC/m2/day to mgC/cm3/day.  
          ! Later this flux will be divided by the number of hourly timesteps.

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
          ! set limits on Lignin:N to keep fmet > 0.2 -ww 11/20 
          ! necessary for litter quality in boreal forests with high cwd flux
          mimicsbiome%ligninNratioAvg(npt) = min(40.0, mimicsbiome%ligninNratioAvg(npt))

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
!     if(mimicsbiome%fmet(npt) < 0.0) then
!         write(*,'(a28,i5,a2,f10.6,a8,i5,a31,f10.6)') '  WARNING: mimicsbiome%fmet(', npt, ')=', mimicsbiome%fmet(npt), &
!             '  PFT =', casamet%iveg2(npt), '  mimicsbiome%ligninNratioAvg =', mimicsbiome%ligninNratioAvg(npt)
!         write(*,*) '  Resetting fmet to zero...'
!         mimicsbiome%fmet(npt) = 0.0
!     endif

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
!         STOP
!     endif
  ENDDO


END SUBROUTINE mimics_delplant

!--------------------------------------------------------------------------------

!ATTENTION: update this to call mimics_delplant and only include additional
!N calculations!

SUBROUTINE mimics_delplant_CN(veg,casabiome,casapool,casaflux,casamet,            &
                           cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                           nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                           pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                           cwd2co2,cwd2str,nwd2str)

! Calculate the change in plant C pools when the MIMICS model is used

  TYPE (veg_parameter_type), INTENT(INOUT) :: veg        ! vegetation parameters
  TYPE (casa_biome),         INTENT(INOUT) :: casabiome
  TYPE (casa_pool),          INTENT(INOUT) :: casapool
  TYPE (casa_flux),          INTENT(INOUT) :: casaflux
  TYPE (casa_met),           INTENT(INOUT) :: casamet

  real(r_2), dimension(mp),INTENT(OUT) :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                                     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                                     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                                     cwd2co2,cwd2str,nwd2str
  ! Local Variables
  integer ::  npt
  real(r_2):: fracLeafMetabolic

  casaflux%FluxCtolitter = 0.0
  casaflux%FluxNtolitter = 0.0
  casaflux%FluxPtolitter = 0.0
  ! Added root exudate flux -mdh 1/13/2020
  casaflux%Cexudate = 0.0
  casaflux%Nexudate = 0.0
  casaflux%Pexudate = 0.0

  mimicsbiome%ligninNratioAvg(:) = 0.0
  casapool%dClitterdt(:,:) = 0.0
  casapool%dNlitterdt(:,:) = 0.0

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
  nwd2str = 0.0

  pleaf2met = 0.0
  pleaf2str = 0.0
  proot2met = 0.0
  proot2str = 0.0
  pwood2cwd = 0.0

  DO npt=1,mp
      IF (casamet%iveg2(npt) /= icewater) THEN

          casapool%dCplantdt(npt,:)= casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)   &
                                     - casaflux%kplant(npt,:) * casapool%cplant(npt,:)

          ! Compute root exudate C flux as a fraction of froot NPP. -mdh 1/13/2020
          casaflux%Cexudate(npt) = max(0.0, casabiome%fracRootExud(veg%iveg(npt)) * casaflux%Cnpp(npt) &
                                   * casaflux%fracCalloc(npt,froot))  
          casapool%dCplantdt(npt,froot) = casapool%dCplantdt(npt,froot) - casaflux%Cexudate(npt) 
   
          ! Calculate fraction C to labile pool as a fraction of GPP, not NPP
          casapool%dClabiledt(npt) = casaflux%Cgpp(npt) * casaflux%fracClabile(npt) - casaflux%clabloss(npt)
   
          !Compute litter inputs from casa (gC/m2/day)
          cleaf2met(npt) = casaflux%fromPtoL(npt,metb,leaf)  * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          cleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf)   * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
          ! Add casaflux%Cexudate(npt) to metabolic litter. -mdh 1/13/2019
          !croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
          croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot) &
                           + casaflux%Cexudate(npt)
          croot2str(npt) = casaflux%fromPtoL(npt,str,froot)  * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
          cwood2cwd(npt) = casaflux%fromPtoL(npt,cwd,wood)   * casaflux%kplant(npt,wood)  * casapool%cplant(npt,wood)

          !---------------------------------------------------------------------------------------------------- 
          ! This code is from casa_delplant.
          if (icycle > 1) then 
             if(casaflux%fracNalloc(npt,leaf)==0.0) THEN
                casapool%dNplantdt(npt,leaf)  = -casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf)
             else
                casapool%dNplantdt(npt,leaf)  = -casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf) &
                                              * casabiome%ftransNPtoL(veg%iveg(npt),leaf)
             endif
             casapool%dNplantdt(npt,wood)  = -casaflux%kplant(npt,wood) * casapool%Nplant(npt,wood) &
                                              * casabiome%ftransNPtoL(veg%iveg(npt),wood)
             casapool%dNplantdt(npt,froot)  = -casaflux%kplant(npt,froot) * casapool%Nplant(npt,froot) &
                                              * casabiome%ftransNPtoL(veg%iveg(npt),froot)

             ! Compute root exudate N flux as a fraction of froot N uptake. -mdh 1/13/2020
             if (casaflux%Cexudate(npt) > 0.0) then
               casaflux%Nexudate(npt) = max(0.0,casabiome%fracRootExud(veg%iveg(npt)) * casaflux%Nminuptake(npt) &
                                        * casaflux%fracNalloc(npt,froot))
               casapool%dNplantdt(npt,froot) = casapool%dNplantdt(npt,froot) - casaflux%Nexudate(npt)
             endif

             ! added by ypwang 5/nov/2012
      
             nleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                            * casapool%cplant(npt,leaf)       * ratioNCstrfix
             nroot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                            * casapool%cplant(npt,froot)      * ratioNCstrfix
      
             ! dNplantdt includes casaflux%Nexudate. -mdh 1/13/2020
             nleaf2met(npt) = -casapool%dNplantdt(npt,leaf)  - nleaf2str(npt)
             nroot2met(npt) = -casapool%dNplantdt(npt,froot) - nroot2str(npt)
             nwood2cwd(npt) = -casapool%dNplantdt(npt,wood)

             ! ATTENTION: is this correct? -mdh 6/23/2019
             ! casaflux%Nminuptake and casaflux%fracNalloc computed in casa_nuptake
             ! adding N uptake
             casapool%dNplantdt(npt,:) = casapool%dNplantdt(npt,:) &
                                          + casaflux%Nminuptake(npt)*casaflux%fracNalloc(npt,:) 
          endif
         !---------------------------------------------------------------------------------------------------- 

          if ((cleaf2met(npt) + cleaf2str(npt)) .gt. 0.0 ) then
              fracLeafMetabolic = cleaf2met(npt) / (cleaf2met(npt) + cleaf2str(npt))
          else
              fracLeafMetabolic = 0.0
          endif

          !For MIMICS, litter C&N inputs are converted from g/m2/day to mg/cm3/day  

          !! Will suggested that CWD does not go into structural material all at once (4/20/2015).
          !! The decay rate of CWD and fraction lost to Rh are computed in mimics_xratesoil (-mdh 4/26/2015)
          !!   casaflux%klitter(npt,cwd) - decay rate of CWD
          !!   fromLtoCO2(mp,mlitter)    - fraction of decomposed litter emitted as CO2 (fraction)

          cwd2co2(npt) =  casaflux%fromLtoCO2(npt,cwd) * casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)
          cwd2str(npt) = (1.0-casaflux%fromLtoCO2(npt,cwd)) * casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)
!         write(*,*)
!         write(*,*) 'mimics_delplant_CN:'
!         write(*,'(a14,f10.6)') 'cwd2co2(npt) =', cwd2co2(npt)
!         write(*,'(a14,f10.6)') 'cwd2str(npt) =', cwd2str(npt)
!         write(*,'(a30,f10.6)') 'casaflux%fromLtoCO2(npt,cwd) =', casaflux%fromLtoCO2(npt,cwd)
!         write(*,'(a30,f10.6)') 'casaflux%klitter(npt,cwd) =', casaflux%klitter(npt,cwd)
!         write(*,'(a30,f18.6)') 'casapool%clitter(npt,cwd) =', casapool%clitter(npt,cwd)
          
          ! ATTENTION: Is nwd2str calculation correct? -mdh 6/23/2019
          if (icycle > 1) then 
              nwd2str(npt) = (1.0-casaflux%fromLtoCO2(npt,cwd)) &
                             * casaflux%klitter(npt,cwd) * casapool%nlitter(npt,cwd)
          endif

          ! Calculate change in CWD pool.  CWD pool is updated in mimics_ccycle (-mdh 5/4/2015)
          casapool%dClitterdt(npt,cwd) = cwood2cwd(npt) - casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)
          ! Add N component of CWD. -mdh 7/1/2019
          ! ATTENTION: not really sure about this calculation. 
          casapool%dNlitterdt(npt,cwd) = nwood2cwd(npt) - casaflux%klitter(npt,cwd) * casapool%nlitter(npt,cwd)

          !ATTENTION. make sure all intermediate calculations are done correctly!
          !casapool%dNlitterdt(npt,:) =  casaflux%FluxNtolitter(npt,:)  &
          !                         - casaflux%klitter(npt,:) &
          !                         * max(0.0,casapool%Nlitter(npt,:))


          ! mimicsflux%ClitInput(*,*) is both a flux to the MIMICS litter pools and an output variable
          ! Convert litter inputs from g/m2/day to mg/cm3/day.  
          ! In mimics_soil subroutines these fluxes will be divided by the number of hourly timesteps.

          mimicsflux%ClitInput(npt,metbc) = ( cleaf2met(npt) &
                                            + croot2met(npt)) * 0.1 / mimicsbiome%depth(veg%iveg(npt))

          mimicsflux%ClitInput(npt,struc) = ( cleaf2str(npt) &
                                            + croot2str(npt) &
                                            + cwd2str(npt)) * 0.1 / mimicsbiome%depth(veg%iveg(npt))

          if (icycle > 1) then
            mimicsflux%NlitInput(npt,metbc) = ( nleaf2met(npt) &
                                            + nroot2met(npt)) * 0.1 / mimicsbiome%depth(veg%iveg(npt))

            mimicsflux%NlitInput(npt,struc) = ( nleaf2str(npt) &
                                            + nroot2str(npt) &
                                            + nwd2str(npt)) * 0.1 / mimicsbiome%depth(veg%iveg(npt))
          endif

          ! For consistency, compute the Litter Inputs that are written to casa's netcdf file. -mdh 1/23/2017
          casaflux%ClitInptMet(npt) = cleaf2met(npt) + croot2met(npt)
          casaflux%ClitInptStruc(npt) = cleaf2str(npt) + croot2str(npt)+ cwd2str(npt)

          if (icycle > 1) then
            casaflux%NlitInptMet(npt) = nleaf2met(npt) + nroot2met(npt)
            casaflux%NlitInptStruc(npt) = nleaf2str(npt) + nroot2str(npt)+ nwd2str(npt)
          endif

          ! Weighted average of lignin:N ratio for all litter inputs. Used for mimicsbiome%fmet(npt) calculation below
          mimicsbiome%ligninNratioAvg(npt) =  &
                                ( mimicsbiome%ligninNratio(npt,leaf)  * (cleaf2met(npt) + cleaf2str(npt)) &
                                + mimicsbiome%ligninNratio(npt,froot) * (croot2met(npt) + croot2str(npt)) &
                                + mimicsbiome%ligninNratio(npt,wood)  * (cwd2str(npt)) ) &
                                / max(0.001,cleaf2met(npt)+cleaf2str(npt)+croot2met(npt)+croot2str(npt)+cwd2str(npt))
          ! set limits on Lignin:N to keep fmet > 0.2 -ww 11/20 
          ! necessary for litter quality in boreal forests with high cwd flux
          mimicsbiome%ligninNratioAvg(npt) = min(40.0, mimicsbiome%ligninNratioAvg(npt))

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
!     if(mimicsbiome%fmet(npt) < 0.0) then
!         write(*,'(a28,i5,a2,f10.6,a8,i5,a31,f10.6)') '  WARNING: mimicsbiome%fmet(', npt, ')=', mimicsbiome%fmet(npt), &
!             '  PFT =', casamet%iveg2(npt), '  mimicsbiome%ligninNratioAvg =', mimicsbiome%ligninNratioAvg(npt)
!         write(*,*) '  Resetting fmet to zero...'
!         mimicsbiome%fmet(npt) = 0.0
!     endif

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
!         write(*,'(a21,f6.4,1x,f6.4)') 'fAVAL(1), fAVAL(2) = ', mimicsbiome%fAVAL(npt,1), mimicsbiome%fAVAL(npt,2)
!         write(*,'(a21,f6.4,1x,f6.4)') 'fCHEM(1), fCHEM(2) = ', mimicsbiome%fCHEM(npt,1), mimicsbiome%fCHEM(npt,2)
!         write(*,'(a21,f6.4,1x,f6.4)') 'fPHYS(1), fPHYS(2) = ', mimicsbiome%fPHYS(npt,1), mimicsbiome%fPHYS(npt,2)
!         STOP
!     endif

  ENDDO


END SUBROUTINE mimics_delplant_CN


!--------------------------------------------------------------------------------
!  SUBROUTINE mimics_soil_forward_MM
!    - calculate daily changes in litter, microbial, and SOM pools

SUBROUTINE mimics_soil_forwardMM(mp,iYrCnt,idoy,mdaily,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)

  ! Function Arguments
  integer, INTENT(IN) :: mp,iYrCnt,idoy  ! number of grid points, simulation year count, day of year
  integer, INTENT(IN) :: mdaily
  real(r_2), dimension(mp),INTENT(IN) :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd

  ! Local Variables
  integer :: npt, ihr
  integer, parameter :: NHOURS = 24
  real(r_2) :: LITmin(4), MICtrn(6), SOMmin(2), DEsorp, OXIDAT, desorption
  real(r_2) :: dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk
  real(r_2) :: NHOURSf
  real(r_2) :: Tsoil           ! average soil temperature for the day (degrees C)
  real(r_2) :: theta_liq       ! WW average liquid soil water
  real(r_2) :: theta_frzn      ! WW average frozen soil water
  real(r_2) :: air_filled_porosity !Fraction of 1.0.  Different from 1.0-theta_liq since it includes ice
  real(r_2) :: fW              ! CORPSE moisture function theta_liq^3*(1-air_filled_porosity)^2.5,
                               ! WW adjusted to give max values of 1, min = 0.01
  real(r_2) :: Cbalance
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps 
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)

! if (iptToSave_mimics > 0) then
!     open(214,file=sPtFileNameMIMICS, access='APPEND')
! endif

  NHOURSf = real(NHOURS)
  Cbalance = -9999.9

  do npt=1,mp

  IF(casamet%iveg2(npt) /= icewater) THEN

      mimicsflux%Chresp(npt) = 0.0
      mimicsflux%CSOMpInput(npt) = 0.0

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

      ! Desorption a function of soil temperature, Q10 = 1.1 w/ reference
      ! temperature of 25C. -WW 4/6/2021. Created parameter names -mdh 4/12/2021
      ! desorption = mimicsbiome%desorp(npt) * (1.1 * exp((Tsoil-25)/10))
      ! Make the Q10 modification to desorption optional. -mdh 6/30/2021
      if (mimicsbiome%desorpQ10 < 0.0) then
          desorption = mimicsbiome%desorp(npt) 
      else
          desorption = mimicsbiome%desorp(npt) * (mimicsbiome%desorpQ10 * exp((Tsoil-mimicsbiome%desorpTref)/10.0))
      endif
    
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
          !DEsorp = mimicspool%SOMp(npt) * mimicsbiome%desorp(npt) 
          DEsorp = mimicspool%SOMp(npt) * desorption 

          ! Oxidation of SOMc to SOMa (f10. SOMc-->SOMa)
          OXIDAT = (mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%SOMc(npt)    &
                    / (mimicsbiome%KO(2)*mimicsbiome%Km(npt,K2) + mimicspool%SOMc(npt)))      &
                    + (mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%SOMc(npt) &
                    / (mimicsbiome%KO(1)*mimicsbiome%Km(npt,R2) + mimicspool%SOMc(npt)))
    
          !Divide total litter inputs (mgC/cm3/day) by number of hours in the daily timestep (NHOURSf)
          dLITm = mimicsflux%ClitInput(npt,metbc)/NHOURSf * (1.0-mimicsbiome%Fi(metbc)) - LITmin(1) - LITmin(3)
          dMICr = mimicsbiome%MGE(1)*(LITmin(1)+ SOMmin(1)) + mimicsbiome%MGE(2)*(LITmin(2)) - sum(MICtrn(1:3))
          dSOMp = mimicsflux%ClitInput(npt,metbc)/NHOURSf * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4)- DEsorp 
          dLITs = mimicsflux%ClitInput(npt,struc)/NHOURSf * (1.0-mimicsbiome%Fi(struc)) - LITmin(2) - LITmin(4)
          dMICk = mimicsbiome%MGE(3)*(LITmin(3)+ SOMmin(2)) + mimicsbiome%MGE(4)*(LITmin(4)) - sum(MICtrn(4:6)) 
          dSOMc = mimicsflux%ClitInput(npt,struc)/NHOURSf * mimicsbiome%Fi(struc) + MICtrn(2) + MICtrn(5) - OXIDAT
          dSOMa = MICtrn(3) + MICtrn(6) + DEsorp + OXIDAT - SOMmin(1) - SOMmin(2)
      
          !Sum daily heterotrphic respiration flux (mgC/cm3) 
          mimicsflux%Chresp(npt) = mimicsflux%Chresp(npt)  &
                                   + (1.0 - mimicsbiome%MGE(1)) * (LITmin(1) + SOMmin(1)) &
                                   + (1.0 - mimicsbiome%MGE(2)) * (LITmin(2))             &
                                   + (1.0 - mimicsbiome%MGE(3)) * (LITmin(3) + SOMmin(2)) &
                                   + (1.0 - mimicsbiome%MGE(4)) * (LITmin(4))

          !Sum daily inputs to SOMp (mgC/cm3). -mdh 12/3/2018
          mimicsflux%CSOMpInput(npt) = mimicsflux%CSOMpInput(npt)  &
              + mimicsflux%ClitInput(npt,metbc)/NHOURSf * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4) 
    
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

!     if (casamet%ijgcm(npt) .eq. iptToSave_mimics) then 
      if (npt .eq. iptToSave_mimics) then 
          ! Write to point file daily if mdaily==1, else write once a year on day 365. -mdh 5/14/2018
          if ((mdaily == 1) .or. (idoy==365)) then
              call WritePointMIMICS(214, sPtFileNameMIMICS, npt, mp, iYrCnt, idoy, &
                  cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
                  LITmin, MICtrn, SOMmin, DEsorp, OXIDAT, &
                  dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil, Cbalance)
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
  real(r_2), dimension(mp),INTENT(IN) :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd

  ! Local Variables
  integer :: npt, ihr
  integer, parameter :: NHOURS = 24
  real(r_2) :: LITmin(4), MICtrn(6), SOMmin(2), DEsorp, OXIDAT, desorption
  real(r_2) :: dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk
  real(r_2) :: NHOURSf, NHOURSfrac
  real(r_2) :: Tsoil           ! average soil temperature for the day (degrees C)
  real(r_2) :: theta_liq       ! WW average liquid soil water
  real(r_2) :: theta_frzn      ! WW average frozen soil water
  real(r_2) :: air_filled_porosity !Fraction of 1.0.  Different from 1.0-theta_liq since it includes ice
  real(r_2) :: fW              ! CORPSE moisture function theta_liq^3*(1-air_filled_porosity)^2.5,
                               ! WW adjusted to give max values of 1, min =
                               ! 0.01
  real(r_2) :: Cbalance, CorgSum1, CorgSum2
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps 
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)

! if (iptToSave_mimics > 0) then
!     open(214,file=sPtFileNameMIMICS, access='APPEND')
! endif

  NHOURSf = real(NHOURS)
  NHOURSfrac = 1.0/NHOURSf
  Cbalance = -9999.9

  do npt=1,mp

  IF(casamet%iveg2(npt) /= icewater) THEN

      mimicsflux%Chresp(npt) = 0.0
      mimicsflux%CSOMpInput(npt) = 0.0

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
    
      CorgSum1 = mimicspool%LITm(npt) + mimicspool%LITs(npt) + &
                 mimicspool%SOMa(npt) + mimicspool%SOMc(npt) + mimicspool%SOMp(npt)

      ! Desorption a function of soil temperature, Q10 = 1.1 w/ reference
      ! temperature of 25C. -WW 4/6/2021. Created parameter names -mdh 4/12/2021
      ! desorption = mimicsbiome%desorp(npt) * (1.1 * exp((Tsoil-25)/10))
      ! Make the Q10 modification to desorption optional. -mdh 6/30/2021
      if (mimicsbiome%desorpQ10 < 0.0) then
          desorption = mimicsbiome%desorp(npt) 
      else
          desorption = mimicsbiome%desorp(npt) * (mimicsbiome%desorpQ10 * exp((Tsoil-mimicsbiome%desorpTref)/10.0))
      endif

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
          !DEsorp = mimicspool%SOMp(npt) * mimicsbiome%desorp(npt)  
          DEsorp = mimicspool%SOMp(npt) * desorption  

          ! Oxidation of SOMc to SOMa (f10. SOMc-->SOMa), reverse Michaelis-Menton Kinetics
          OXIDAT = (mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%SOMc(npt)    &
                    / (mimicsbiome%KO(2)*mimicsbiome%Km(npt,K2) + mimicspool%MICk(npt)))      &
                    + (mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%SOMc(npt) &
                    / (mimicsbiome%KO(1)*mimicsbiome%Km(npt,R2) + mimicspool%MICr(npt)))
    
          ! Divide total litter inputs (mgC/cm3/day) by number of hours in the daily timestep (NHOURSf)
          ! Optimization: Multiply by NHOURSfrac = 1/NHOURSf
          dLITm = mimicsflux%ClitInput(npt,metbc)*NHOURSfrac * (1.0-mimicsbiome%Fi(metbc)) - LITmin(1) - LITmin(3)
          dMICr = mimicsbiome%MGE(1)*(LITmin(1)+ SOMmin(1)) + mimicsbiome%MGE(2)*(LITmin(2)) - sum(MICtrn(1:3))
          dSOMp = mimicsflux%ClitInput(npt,metbc)*NHOURSfrac * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4)- DEsorp 
          dLITs = mimicsflux%ClitInput(npt,struc)*NHOURSfrac * (1.0-mimicsbiome%Fi(struc)) - LITmin(2) - LITmin(4)
          dMICk = mimicsbiome%MGE(3)*(LITmin(3)+ SOMmin(2)) + mimicsbiome%MGE(4)*(LITmin(4)) - sum(MICtrn(4:6)) 
          dSOMc = mimicsflux%ClitInput(npt,struc)*NHOURSfrac * mimicsbiome%Fi(struc) + MICtrn(2) + MICtrn(5) - OXIDAT
          dSOMa = MICtrn(3) + MICtrn(6) + DEsorp + OXIDAT - SOMmin(1) - SOMmin(2)
      
          !Sum daily heterotrophic respiration flux (mgC/cm3) 
          mimicsflux%Chresp(npt) = mimicsflux%Chresp(npt)  &
                                   + (1.0 - mimicsbiome%MGE(1)) * (LITmin(1) + SOMmin(1)) &
                                   + (1.0 - mimicsbiome%MGE(2)) * (LITmin(2))             &
                                   + (1.0 - mimicsbiome%MGE(3)) * (LITmin(3) + SOMmin(2)) &
                                   + (1.0 - mimicsbiome%MGE(4)) * (LITmin(4))

          !Sum daily inputs to SOMp (mgC/cm3). -mdh 12/3/2018
          mimicsflux%CSOMpInput(npt) = mimicsflux%CSOMpInput(npt)  &
              + mimicsflux%ClitInput(npt,metbc)*NHOURSfrac * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4)
    
          !Update pools (mgC/cm3)
          mimicspool%LITm(npt) = mimicspool%LITm(npt) + dLITm
          mimicspool%LITs(npt) = mimicspool%LITs(npt) + dLITs
          mimicspool%MICr(npt) = mimicspool%MICr(npt) + dMICr
          mimicspool%MICk(npt) = mimicspool%MICk(npt) + dMICk
          mimicspool%SOMa(npt) = mimicspool%SOMa(npt) + dSOMa
          mimicspool%SOMc(npt) = mimicspool%SOMc(npt) + dSOMc
          mimicspool%SOMp(npt) = mimicspool%SOMp(npt) + dSOMp

      end do

      CorgSum2 = mimicspool%LITm(npt) + mimicspool%LITs(npt) + &
                 mimicspool%SOMa(npt) + mimicspool%SOMc(npt) + mimicspool%SOMp(npt)
      ! CorgSum2 = CorgSum1 + gains - losses
      Cbalance = -1.0*(CorgSum1 - CorgSum2 - mimicsflux%Chresp(npt) &
                  + mimicsflux%ClitInput(npt,metbc)+mimicsflux%ClitInput(npt,struc))

!     if (casamet%ijgcm(npt) .eq. 2945) then     ! Crop (12)
!     if (casamet%ijgcm(npt) .eq. 10354) then    ! Tundra (18)
!     if (casamet%ijgcm(npt) .eq. 10052) then    ! Grassland (10)
!     if (casamet%ijgcm(npt) .eq. 10919) then    ! Deciduous broadleaf (4)
!     if (casamet%ijgcm(npt) .eq. 11018) then    ! Evergreen needleleaf (1)
!     if (casamet%ijgcm(npt) .eq. 11569) then    ! Evergreen broadleaf (2)

      if (.FALSE.) then
!     if (.TRUE.) then
          write(*,*)
          write(*,*) 'C balance > 0.0 means unexpected gain of C from start to end of day'
          write(*,'(a13,f10.6)') 'Cbalance = ', Cbalance
          write(*,'(a13,f10.6)') 'CorgSum2 = ', CorgSum2
          write(*,'(a13,f10.6)') 'CorgSum1 = ', CorgSum1
          write(*,'(a13,f10.6)') 'dCorg = ', CorgSum2 - CorgSum1
          write(*,'(a13,f10.6)') 'litCinput = ', mimicsflux%ClitInput(npt,metbc) + mimicsflux%ClitInput(npt,struc)
          write(*,'(a13,f10.6)') 'Chresp = ', mimicsflux%Chresp(npt)
      endif

!     if (casamet%ijgcm(npt) .eq. iptToSave_mimics) then 
      if (npt .eq. iptToSave_mimics) then 
          ! Write to point file daily if mdaily==1, else write once a year on day 365. -mdh 5/14/2018
          if ((mdaily == 1) .or. (idoy==365)) then
              call WritePointMIMICS(214, sPtFileNameMIMICS, npt, mp, iYrCnt, idoy, &
                  cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
                  LITmin, MICtrn, SOMmin, DEsorp, OXIDAT, &
                  dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil, Cbalance)
          endif 
      endif 

  ENDIF

  end do

END SUBROUTINE mimics_soil_reverseMM

!--------------------------------------------------------------------------------

SUBROUTINE mimics_ccycle(veg,casabiome,casapool,casaflux,casamet)

! Called daily
! Update all C pool that were not updated in mimics_soil_forwardMM or mimics_soil_reverseMM.
! This includes plant C pools, but not litter C or soil C pools.

  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! Local variables
  INTEGER i,npt

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
SUBROUTINE mimics_cncycle(veg,casabiome,casapool,casaflux,casamet)

! Called daily
! Update all C&N pools that were not updated in mimics_soil_reverseMM_CN.
! This includes plant C&N pools, and mineral N pool, but not litter C&N or soil organic C&N pools.

  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! Local variables
  INTEGER i,npt

  DO npt=1,mp
      IF (casamet%iveg2(npt) == icewater) THEN
          casamet%glai(npt) = 0.0
      ELSE  
          casapool%cplant(npt,:) = casapool%cplant(npt,:) + casapool%dCplantdt(npt,:) * deltpool 
          casapool%clabile(npt)  = casapool%clabile(npt)  + casapool%dClabiledt(npt)  * deltpool 

          IF (casapool%cplant(npt,leaf) > 0.0) THEN
            IF(icycle > 1) casapool%Nplant(npt,:) = casapool%Nplant(npt,:) &
                                 +casapool%dNplantdt(npt,:)*deltpool
          ENDIF

          casamet%glai(npt)      = MAX(casabiome%glaimin(veg%iveg(npt)), &
                                       casabiome%sla(veg%iveg(npt)) * casapool%cplant(npt,leaf))
          casamet%glai(npt)      = MIN(casabiome%glaimax(veg%iveg(npt)), casamet%glai(npt))
      
          ! Update CWD pool (-mdh 5/4/2015)
          casapool%clitter(npt,:) = casapool%clitter(npt,:) + casapool%dClitterdt(npt,:) * deltpool 
!         write(*,*)
!         write(*,*) 'mimics_cncycle:'
!         write(*,'(a32,f10.6)') 'casapool%dClitterdt(npt,metb) =', casapool%dClitterdt(npt,metb)
!         write(*,'(a32,f10.6)') 'casapool%dClitterdt(npt,str) =', casapool%dClitterdt(npt,str)
!         write(*,'(a32,f18.6)') 'casapool%dClitterdt(npt,cwd) =', casapool%dClitterdt(npt,cwd)

          IF (icycle > 1) THEN
            casapool%Nlitter(npt,:) = casapool%Nlitter(npt,:) + casapool%dNlitterdt(npt,:)* deltpool
            !casapool%Nsoil(npt,:)   = casapool%Nsoil(npt,:) + casapool%dNsoildt(npt,:)  * deltpool
            casapool%Nsoilmin(npt)  = casapool%Nsoilmin(npt) + casapool%dNsoilmindt(npt) * deltpool
          ENDIF  ! end of "icycle > 1"

          DO i=1,mplant
            IF (casapool%cplant(npt,i) < 0.0) THEN
!               WRITE(57,*) 'Cpool: npt,ivt',npt,casamet%iveg2(npt),casapool%cplant(npt,:)
                !call casa_poolzero(npt,1,casapool)
                casapool%cplant(npt,i) = max(0.0, casapool%cplant(npt,i))
            ENDIF
          ENDDO
          IF (icycle > 1) THEN
            DO i=1,mplant
              IF (casapool%nplant(npt,i) < 0.0) THEN
!                WRITE(57,*) 'Npool:', 'npt,ivt,ipool',npt,casamet%iveg2(np),casapool%nplant(npt,:)
                !call casa_poolzero(npt,2,casapool)
                casapool%nplant(npt,i) = max(0.0, casapool%nplant(npt,i))
              ENDIF
            ENDDO
          ENDIF ! end of "icycle > 1"
      
      ENDIF
  ENDDO !end of "npt"

END SUBROUTINE mimics_cncycle

!--------------------------------------------------------------------------------

!! Subroutine mimics_caccum accumulates daily C pool and flux values
!! Average annual values will be computed at the end of the simulation 
!!   in subroutine mimics_poolfluxout

SUBROUTINE mimics_caccum(mp,cwd2co2)

  !Subroutine arguments
  implicit none
  integer, INTENT(IN) :: mp                        ! number of grid cells
  real(r_2), dimension(mp),INTENT(IN) :: cwd2co2   ! respiration from decay of CWD (gC/m2/day)

  !Local variables
  integer   :: npt
  real(r_2) :: unitConv                         ! mgC/cm3 * depth(cm)* (1g/10^3 mg)*(10^4 cm2)/m2 = gC/m2

  do npt = 1,mp

    if (casamet%iveg2(npt) /= icewater) then

      unitConv = 10.0 * mimicsbiome%depth(veg%iveg(npt))
    
!     if (npt .eq. 1003) then
!     print *
!     print *, 'caccum: unitConv = ', unitConv
!     print *, 'caccum: mimicspoolAn%ClitterAn(',npt,',METBC) = ', mimicspoolAn%ClitterAn(npt,METBC), mimicspool%LITm(npt)
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
      ! mimicsflux%CSOMpInput is in units of mg/c3 whereas mimicsfluxAn%CSOMpInputAn is gC/m2. -mdh 12/3/2018
      mimicsfluxAn%CSOMpInputAn(npt) = mimicsfluxAn%CSOMpInputAn(npt) + mimicsflux%CSOMpInput(npt) * unitConv 
      mimicsfluxAn%Overflow_rAn(npt) = mimicsfluxAn%Overflow_rAn(npt) + mimicsflux%Overflow_r(npt) * unitConv 
      mimicsfluxAn%Overflow_kAn(npt) = mimicsfluxAn%Overflow_kAn(npt) + mimicsflux%Overflow_k(npt) * unitConv 

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

!! Subroutine mimics_naccum accumulates daily N pool and flux values
!! Average annual values will be computed at the end of the simulation 
!!   in subroutine mimics_poolfluxout

SUBROUTINE mimics_naccum(mp)

  !Subroutine arguments
  implicit none
  integer, INTENT(IN) :: mp                     ! number of grid cells

  !Local variables
  integer   :: npt
  real(r_2) :: unitConv                         ! mgN/cm3 * depth(cm)* (1g/10^3 mg)*(10^4 cm2)/m2 = gN/m2

  do npt = 1,mp

    if (casamet%iveg2(npt) /= icewater) then

      unitConv = 10.0 * mimicsbiome%depth(veg%iveg(npt))
    
!     if (npt .eq. 1003) then
!     print *
!     print *, 'naccum: unitConv = ', unitConv
!     print *, 'naccum: mimicspoolAn%NlitterAn(',npt,',METBC) = ', mimicspoolAn%NlitterAn(npt,METBC), mimicspool%LITmN(npt)
!     print *, 'naccum: mimicspoolAn%NlitterAn(',npt,',STRUC) = ', mimicspoolAn%NlitterAn(npt,STRUC), mimicspool%LITsN(npt)
!
!     print *, 'naccum: mimicspoolAn%NmicrobeAn(',npt,',RSEL) = ', mimicspoolAn%NmicrobeAn(npt,RSEL), mimicspool%MICrN(npt)
!     print *, 'naccum: mimicspoolAn%NmicrobeAn(',npt,',KSEL) = ', mimicspoolAn%NmicrobeAn(npt,KSEL), mimicspool%MICkN(npt)
!
!     print *, 'naccum: mimicspoolAn%NsoilAn(',npt,',AVAL) = ', mimicspoolAn%NsoilAn(npt,AVAL), mimicspool%SOMaN(npt)
!     print *, 'naccum: mimicspoolAn%NsoilAn(',npt,',CHEM) = ', mimicspoolAn%NsoilAn(npt,CHEM), mimicspool%SOMcN(npt)
!     print *, 'naccum: mimicspoolAn%NsoilAn(',npt,',PHYS) = ', mimicspoolAn%NsoilAn(npt,PHYS), mimicspool%SOMpN(npt)
!     endif

      mimicspoolAn%NlitterAn(npt,METBC) = mimicspoolAn%NlitterAn(npt,METBC) + mimicspool%LITmN(npt) * unitConv
      mimicspoolAn%NlitterAn(npt,STRUC) = mimicspoolAn%NlitterAn(npt,STRUC) + mimicspool%LITsN(npt) * unitConv
    
      mimicspoolAn%NmicrobeAn(npt,RSEL) = mimicspoolAn%NmicrobeAn(npt,RSEL) + mimicspool%MICrN(npt) * unitConv
      mimicspoolAn%NmicrobeAn(npt,KSEL) = mimicspoolAn%NmicrobeAn(npt,KSEL) + mimicspool%MICkN(npt) * unitConv
    
      mimicspoolAn%NsoilAn(npt,AVAL) = mimicspoolAn%NsoilAn(npt,AVAL) + mimicspool%SOMaN(npt) * unitConv
      mimicspoolAn%NsoilAn(npt,CHEM) = mimicspoolAn%NsoilAn(npt,CHEM) + mimicspool%SOMcN(npt) * unitConv
      mimicspoolAn%NsoilAn(npt,PHYS) = mimicspoolAn%NsoilAn(npt,PHYS) + mimicspool%SOMpN(npt) * unitConv
      mimicspoolAn%DINAn(npt)        = mimicspoolAn%DINAn(npt)        + mimicspool%DIN(npt) * unitConv
    
      !Litter Inputs and Heterotrophic Respiration (gC/m2)
      mimicsfluxAn%NLitInputAn(npt,METBC) = mimicsfluxAn%NLitInputAn(npt,METBC) + mimicsflux%NLitInput(npt,METBC) * unitConv 
      mimicsfluxAn%NLitInputAn(npt,STRUC) = mimicsfluxAn%NLitInputAn(npt,STRUC) + mimicsflux%NLitInput(npt,STRUC) * unitConv

    endif
  end do

END SUBROUTINE mimics_naccum

!--------------------------------------------------------------------------------
!  SUBROUTINE mimics_soil_reverseMM_CN
!    - calculate daily changes in litter, microbial, and SOM pools
!    - alternative to mimics_soil_forwardMM; implements the reverse Michaelis-Mention kinetics

SUBROUTINE mimics_soil_reverseMM_CN(mp,iYrCnt,idoy,mdaily,cleaf2met,cleaf2str,croot2met, &
                                    croot2str,cwd2str,cwd2co2,cwood2cwd, &
                                    nleaf2met,nleaf2str,nroot2met,nroot2str,nwd2str,nwood2cwd)

  USE define_types

  ! Function Arguments
  integer, INTENT(IN) :: mp,iYrCnt,idoy  ! number of grid points, simulation year count, day of year
  integer, INTENT(IN) :: mdaily
  real(r_2), dimension(mp),INTENT(IN) :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd
  real(r_2), dimension(mp),INTENT(IN) :: nleaf2met,nleaf2str,nroot2met,nroot2str,nwd2str,nwood2cwd

  ! Local Variables
  integer :: npt, ihr
  integer, parameter :: NHOURS = 24
  real(r_2) :: LITmin(4), MICtrn(6), SOMmin(2), DEsorp, OXIDAT
  real(r_2) :: LITminN(4), MICtrnN(6), SOMminN(2), DEsorpN, OXIDATN
  real(r_2) :: dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk
  real(r_2) :: dLITmN, dLITsN, dSOMaN, dSOMcN, dSOMpN, dMICrN, dMICkN, dDIN
  real(r_2) :: upMICrC, upMICrN, upMICkC, upMICkN
  real(r_2) :: CNup_r, CNup_k
  real(r_2) :: DINup_r, DINup_k 
  real(r_2) :: Overflow_r, Overflow_k, Nspill_r, Nspill_k 
  real(r_2) :: unitConv, unitConvCtoM
  real(r_2) :: NHOURSf, NHOURSfrac, NDAYSfrac
  real(r_2) :: MICr_recip, MICk_recip,desorption
  real(r_2) :: Tsoil           ! average soil temperature for the day (degrees C)
  real(r_2) :: theta_liq       ! WW average liquid soil water
  real(r_2) :: theta_frzn      ! WW average frozen soil water
  real(r_2) :: air_filled_porosity !Fraction of 1.0.  Different from 1.0-theta_liq since it includes ice
  real(r_2) :: fW              ! CORPSE moisture function theta_liq^3*(1-air_filled_porosity)^2.5,
                               ! WW adjusted to give max values of 1, min = 0.01
  real(r_2) :: ClitInputNPPhr  ! NPP in mg C/cm3/hour
  real(r_2) :: DINstart  ! Amount of DIN available to microbes at the beginning of the day (mg N/cm3)
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps 
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)

  ! C&N balance variables
  real(r_2) CorgSum1, CorgSum2
  real(r_2) Cbalance
  real(r_2) NorgSum1, NorgSum2
  real(r_2) Nbalance

  !ATTENTION: where are parameters mimicsbiome%fracNimport_r and mimicsbiome%fracNimport_k used? -mdh 6/21/2019

  NHOURSf = real(NHOURS)
  NHOURSfrac = 1.0/NHOURSf
  NDAYSfrac = 1.0/365.0

  do npt=1,mp

  IF(casamet%iveg2(npt) /= icewater) THEN

      ! (mg N/cm3)*(1/1000)(g/mg)*(10000)(cm2/m2)*depth(cm) = g N/m2
      unitConv = 10.0*mimicsbiome%depth(veg%iveg(npt))    ! Convert mgN/cm3 to gN/m2 (MIMICS to CASA) by multipling by this factor
      unitConvCtoM = 1.0/unitConv ! Convert gN/m2 to mgN/cm3 (CASA to MIMICS) to by multipling by this factor

      ! Convert annual NPP (gC/m2/yr) to mgC/cm3/hour
      ClitInputNPPhr = NHOURSfrac * casaflux%CnppAn(npt) * unitConvCtoM * NDAYSfrac 

      ! CASACNP fluxes computed by MIMICS-CN. MIMICS-CN fluxes (mg N/cm3) will be converted to gN/m2.
      casaflux%Nlittermin(npt) = 0.0  ! Gross N mineralization from litter decomposition (gN/m2/day)
      casaflux%Nsmin(npt) = 0.0       ! Gross N mineraliztion from SOM decomposition (gN/m2/day)
      casaflux%Nsimm(npt) = 0.0       ! Soil N immobilization (gN/m2/day)
      casaflux%Nminloss(npt) = 0.0    ! N volatilization (N2O) (gN/m2/day)
      casapool%dNsoilmindt(npt) = 0.0 ! Change in mineral N pool during 24-hour loop (gN/m2/day)

      ! Compute mineral N leaching prior to MIMICS 24-hour loop, and remove N leaching calculations from this loop. -mdh 9/27/2020
      casaflux%Nminleach(npt) = casapool%Nsoilmin(npt) * casaflux%fNminleach(npt)  ! N leached from mineral soil (gN/m2/day)
      casapool%Nsoilmin(npt) = casapool%Nsoilmin(npt) - casaflux%Nminleach(npt)

      ! mimicsbiome%fracDINavailMIC is a new parameter which allows only a fraction of soil mineral N
      ! to be available to microbes.  Emily's model did not have it since there was no N competition
      ! with plants. -mdh 6/21/2019
      mimicspool%DIN(npt) = mimicsbiome%fracDINavailMIC*casapool%Nsoilmin(npt) * unitConvCtoM  ! Convert gN/m2 to mgN/cm3
      DINstart = mimicspool%DIN(npt)

      mimicsflux%Chresp(npt) = 0.0
      mimicsflux%CSOMpInput(npt) = 0.0
      ! Accumulate daily Overflow_r and Overflow_k respiration and add to output netCDF file. -mdh 10/12/2020
      mimicsflux%Overflow_r(npt) = 0.0
      mimicsflux%Overflow_k(npt) = 0.0

      ! Vmax - temperature sensitive maximum reaction velocities (mg C (mg MIC)-1 h-1) 
      Tsoil = casamet%tsoilavg(npt) - tkzeroc
      ! Read in soil moisture data as in CORPSE
      theta_liq  = min(1.0, casamet%moistavg(npt)/soil%ssat(npt))     ! fraction of liquid water-filled pore space (0.0 - 1.0)
      theta_frzn = min(1.0, casamet%frznmoistavg(npt)/soil%ssat(npt)) ! fraction of frozen water-filled pore space (0.0 - 1.0)
      air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)

      if (mimicsbiome%fWFunction .eq. CORPSE) then
        ! CORPSE water scalar, adjusted to give maximum values of 1
        ! Optimiztion: 1.0/0.022600567942709 = 44.2466756824402
        ! fW = (theta_liq**3 * air_filled_porosity**2.5)/0.022600567942709
        fW = (theta_liq**3 * air_filled_porosity**2.5) * 44.2466756824402
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
    
      ! casaflux%CnppAn(npt) =  average annual NPP (gC/m2/yr)
      ! mimicsbiome%tauMod(npt) = SQRT(casaflux%CnppAn(npt)/100)
      mimicsbiome%tauMod(npt) = SQRT(casaflux%CnppAn(npt)/mimicsbiome%tauModDenom)
      mimicsbiome%tauMod(npt) = MAX(mimicsbiome%tauMod_MIN, mimicsbiome%tauMod(npt))
      mimicsbiome%tauMod(npt) = MIN(mimicsbiome%tauMod_MAX, mimicsbiome%tauMod(npt))

      ! Updated based on Will's (-mdh 6/1/2015).  I had been using site-level values, so only tauR changed.
      ! mimicsbiome%tauR(npt) = 5.2 * 0.0001 * exp(0.3 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      ! mimicsbiome%tauK(npt) = 2.4 * 0.0001 * exp(0.1 * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      ! 0.00052  tau_r(1)
      ! 0.3      tau_r(2)
      ! 0.00024  tau_k(1)
      ! 0.1      tau_k(2)
      mimicsbiome%CN_r(npt) = mimicsbiome%CNr * SQRT(mimicsbiome%cnModNum/mimicsbiome%fmet(npt))
      mimicsbiome%CN_k(npt) = mimicsbiome%CNk * SQRT(mimicsbiome%cnModNum/mimicsbiome%fmet(npt))

      mimicsbiome%tauR(npt) = mimicsbiome%tau_r(1) * &
                                exp(mimicsbiome%tau_r(2) * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)
      mimicsbiome%tauK(npt) = mimicsbiome%tau_k(1) * &
                                exp(mimicsbiome%tau_k(2) * mimicsbiome%fmet(npt)) * mimicsbiome%tauMod(npt)

! This code was from the oringinal version of MIMICS-CN. -mdh 2/10/2020
!     ! Modification of tauR and tauK for MIMICS-CN (-mdh 4/29/2019)
!     ! Emily's R code for reference
!     ! turnover      <<- c(5.2e-4*exp(0.3*(fMET)), 2.4e-4*exp(0.1*(fMET)))	#WORKS BETTER FOR N_RESPONSE RATIO
!     ! turnover_MOD1 <<- sqrt(ANPP_C[s]/100)  #basicaily standardize against NWT
!     ! turnover      <<- turnover * turnover_MOD1
!     ! #turnover     <<- turnover * 2 #Sulman et al. 2018
!     ! turnover      <<- turnover * 1.7
!     ! turnover      <<- turnover^2*0.55/(.45*Inputs) #  Inputs(1)=metabolic, Inputs(2)=structural (mgC/cm3/hr)
!
!     ! Multiply ClitInput by NHOURSfrac = 1/NHOURSf
!     !mimicsbiome%tauR(npt) = (1.7*mimicsbiome%tauR(npt))**2.0*0.55/(0.45*mimicsflux%ClitInput(npt,metbc)*NHOURSfrac)
!     !mimicsbiome%tauK(npt) = (1.7*mimicsbiome%tauK(npt))**2.0*0.55/(0.45*mimicsflux%ClitInput(npt,struc)*NHOURSfrac)

      ! Use code updates from Emily 10/24/2019. -mdh 2/10/2020
      ! Also use ANPP, converted to mgC/cm3/hr as the litter input instead of instantaneous litter input. -mdh 2/10/2020
      ! turnover      <<- c(5.2e-4*exp(0.3*(fMET)), 2.4e-4*exp(0.1*(fMET)))	#WORKS BETTER FOR N_RESPONSE RATIO
      ! turnover_MOD1 <<- sqrt(ANPP_C[s]/100)  #basicaily standardize against NWT
      ! turnover_MOD1[turnover_MOD1 < 0.6] <<- 0.6 # correction not used in LIDET resutls 
      ! turnover_MOD1[turnover_MOD1 > 1.3] <<- 1.3      #Sulman et al. 2018
      ! turnover      <<- turnover * turnover_MOD1
      ! turnover <<- turnover/2.2
      ! turnover <<- turnover^2*0.55/(.45*Inputs)

      ! ATTENTION: comment out density dependent modifications to tauR and tauK
      ! because they are causing some cells to produce NANs. -mdh 2/10/2020
      ! ClitInputNPPhr is annual NPP (gC/m2/yr) converted to mgC/cm3/hour
!     mimicsbiome%tauR(npt) = (mimicsbiome%tauR(npt)/2.2)**2.0*0.55/(0.45*ClitInputNPPhr)
!     mimicsbiome%tauK(npt) = (mimicsbiome%tauK(npt)/2.2)**2.0*0.55/(0.45*ClitInputNPPhr)

      ! WW also modify TAU as a function of soil moisture, so things don't
      ! collapse in frozen soils...
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

      ! Desorption a function of soil temperature, Q10 = 1.1 w/ reference
      ! temperature of 25C. -WW 4/6/2021. Created parameter names -mdh 4/12/2021
      ! desorption = mimicsbiome%desorp(npt) * (1.1 * exp((Tsoil-25)/10))
      ! Make the Q10 modification to desorption optional. -mdh 6/30/2021
      if (mimicsbiome%desorpQ10 < 0.0) then
          desorption = mimicsbiome%desorp(npt) 
      else
          desorption = mimicsbiome%desorp(npt) * (mimicsbiome%desorpQ10 * exp((Tsoil-mimicsbiome%desorpTref)/10.0))
      endif

      CorgSum1 = mimicspool%LITm(npt) + mimicspool%LITs(npt) + &
                 mimicspool%SOMa(npt) + mimicspool%SOMc(npt) + mimicspool%SOMp(npt)
      NorgSum1 = mimicspool%LITmN(npt) + mimicspool%LITsN(npt) + &
                 mimicspool%SOMaN(npt) + mimicspool%SOMcN(npt) + mimicspool%SOMpN(npt)

      do ihr = 1, NHOURS

          ! Flows to and from MICr

          ! Multiply by recipical of MICr to speed up calculations. -mdh 10/26/2020
          MICr_recip = 1.0 / (mimicspool%MICr(npt) + 1.0e-10)
  
          !MICr decomp of METABOLIC litter (f1. LITm-->MICr), reverse Michaelis-Menton Kinetics
          LITmin(1) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R1) * mimicspool%LITm(npt) &
                      / (mimicsbiome%Km(npt,R1) + mimicspool%MICr(npt))
          
          LITminN(1) = LITmin(1) * mimicspool%LITmN(npt) / (mimicspool%LITm(npt) + 1.0e-10)
  
          !MICr decomp of STRUCTURAL litter (f2. LITs-->MICr), reverse Michaelis-Menton Kinetics   
          LITmin(2) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%LITs(npt) &
                      / (mimicsbiome%Km(npt,R2) + mimicspool%MICr(npt))  

          LITminN(2) = LITmin(2) * mimicspool%LITsN(npt) / (mimicspool%LITs(npt) + 1.0e-10)
      
          ! MIMICS-CN includes density dependence for Microbial Turnover (-mdh 4/29/2019)    

          !MICr turnover to SOMp (f41. MICr-->SOMp) 
          !MICtrn(1) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fPHYS(npt,1)

          MICtrn(1) = mimicspool%MICr(npt)**(mimicsbiome%densDep) * mimicsbiome%tauR(npt) * mimicsbiome%fPHYS(npt,1)
      
          !MICtrnN(1) = MICtrn(1) * mimicspool%MICrN(npt) / (mimicspool%MICr(npt) + 1.0e-10)
          MICtrnN(1) = MICtrn(1) * mimicspool%MICrN(npt) * MICr_recip

          !MICr turnover to SOMc (f42. MICr-->SOMc)                  
          !MICtrn(2) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fCHEM(npt,1)   
          MICtrn(2) = mimicspool%MICr(npt)**(mimicsbiome%densDep) * mimicsbiome%tauR(npt) * mimicsbiome%fCHEM(npt,1)   

          !MICtrnN(2) = MICtrn(2) * mimicspool%MICrN(npt) / (mimicspool%MICr(npt) + 1.0e-10)
          MICtrnN(2) = MICtrn(2) * mimicspool%MICrN(npt) * MICr_recip
      
          !MICr turnover to SOMa (f43. MICr-->SOMa)               
          !MICtrn(3) = mimicspool%MICr(npt) * mimicsbiome%tauR(npt) * mimicsbiome%fAVAL(npt,1)    
          MICtrn(3) = mimicspool%MICr(npt)**(mimicsbiome%densDep) * mimicsbiome%tauR(npt) * mimicsbiome%fAVAL(npt,1) 

          !MICtrnN(3) = MICtrn(3) * mimicspool%MICrN(npt) / (mimicspool%MICr(npt) + 1.0e-10)   
          MICtrnN(3) = MICtrn(3) * mimicspool%MICrN(npt) * MICr_recip  
      
          !decomp of SOMa by MICr (f3. SOMa-->MICr), reverse Michaelis-Menton Kinetics
          SOMmin(1) = mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R3) * mimicspool%SOMa(npt) &
                      / (mimicsbiome%Km(npt,R3) + mimicspool%MICr(npt))   

          SOMminN(1) = SOMmin(1) * mimicspool%SOMaN(npt) / (mimicspool%SOMa(npt) + 1.0e-10)
      
          !Flows to and from MICk

          ! Multiply by recipical of MICk to speed up calculations. -mdh 10/26/2020
          MICk_recip = 1.0 / (mimicspool%MICk(npt) + 1.0e-10)
      
          !decomp of METABOLIC litter (f5. LITm-->MICk), reverse Michaelis-Menton Kinetics
          LITmin(3) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K1) * mimicspool%LITm(npt) &
                      / (mimicsbiome%Km(npt,K1) + mimicspool%MICk(npt)) 
      
          LITminN(3) = LITmin(3) * mimicspool%LITmN(npt) / (mimicspool%LITm(npt) + 1.0e-10)


          !decomp of STRUCTURAL litter (f6. LITs-->MICk), reverse Michaelis-Menton Kinetics 
          LITmin(4) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%LITs(npt) &
                      / (mimicsbiome%Km(npt,K2) + mimicspool%MICk(npt))  

          LITminN(4) = LITmin(4) * mimicspool%LITsN(npt) / (mimicspool%LITs(npt) + 1.0e-10)

          ! MIMICS-CN includes density dependence for Microbial Turnover (-mdh 4/29/2019)    
      
          !MICk turnover to SOMp (f81. MICk-->SOMp) 
          !MICtrn(4) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fPHYS(npt,2) 
          MICtrn(4) = mimicspool%MICk(npt)**(mimicsbiome%densDep) * mimicsbiome%tauK(npt) * mimicsbiome%fPHYS(npt,2)     

          !MICtrnN(4) = MICtrn(4) * mimicspool%MICkN(npt) / (mimicspool%MICk(npt) + 1.0e-10)
          MICtrnN(4) = MICtrn(4) * mimicspool%MICkN(npt) * MICk_recip
      
          !MICk turnover to SOMc (f82. MICk-->SOMc)             
          !MICtrn(5) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fCHEM(npt,2) 
          MICtrn(5) = mimicspool%MICk(npt)**(mimicsbiome%densDep) * mimicsbiome%tauK(npt) * mimicsbiome%fCHEM(npt,2) 

          !MICtrnN(5) = MICtrn(5) * mimicspool%MICkN(npt) / (mimicspool%MICk(npt) + 1.0e-10)
          MICtrnN(5) = MICtrn(5) * mimicspool%MICkN(npt) * MICk_recip
      
          !MICk turnover to SOMa (f83. MICk-->SOMa)
          !MICtrn(6) = mimicspool%MICk(npt) * mimicsbiome%tauK(npt) * mimicsbiome%fAVAL(npt,2)       
          MICtrn(6) = mimicspool%MICk(npt)**(mimicsbiome%densDep) * mimicsbiome%tauK(npt) * mimicsbiome%fAVAL(npt,2)  

          !MICtrnN(6) = MICtrn(6) * mimicspool%MICkN(npt) / (mimicspool%MICk(npt) + 1.0e-10)     
          MICtrnN(6) = MICtrn(6) * mimicspool%MICkN(npt) * MICk_recip    
      
          !decomp of SOMa by MICk (f7. SOMa-->MICk, reverse Michaelis-Menton Kinetics         
          SOMmin(2) = mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K3) * mimicspool%SOMa(npt) &
                      / (mimicsbiome%Km(npt,K3) + mimicspool%MICk(npt))   

          SOMminN(2) = SOMmin(2) * mimicspool%SOMaN(npt) / (mimicspool%SOMa(npt) + 1.0e-10)
      
          ! Desorbtion of SOMp to SOMa (function of fCLAY) (f9. SOMp-->SOMa)
          !DEsorp = mimicspool%SOMp(npt) * mimicsbiome%desorp(npt)  
          DEsorp = mimicspool%SOMp(npt) * desorption               

          DEsorpN = DEsorp * mimicspool%SOMpN(npt) / (mimicspool%SOMp(npt) + 1.0e-10)

          ! Oxidation of SOMc to SOMa (f10. SOMc-->SOMa), reverse Michaelis-Menton Kinetics
          OXIDAT = (mimicspool%MICk(npt) * mimicsbiome%Vmax(npt,K2) * mimicspool%SOMc(npt)    &
                    / (mimicsbiome%KO(2)*mimicsbiome%Km(npt,K2) + mimicspool%MICk(npt)))      &
                    + (mimicspool%MICr(npt) * mimicsbiome%Vmax(npt,R2) * mimicspool%SOMc(npt) &
                    / (mimicsbiome%KO(1)*mimicsbiome%Km(npt,R2) + mimicspool%MICr(npt)))

          OXIDATN = OXIDAT * mimicspool%SOMcN(npt) / (mimicspool%SOMc(npt) + 1.0e-10)
    

          ! Partitions available DIN between microbial pools based on relative biomass
          ! The amount of DIN uptake is not demand driven.  It is a fraction of a large pool. -mdh 7/29/2019

          ! Assume that microbes take up all mineral N that is available. The unneeded portion will be spilled later.
          DINup_r = mimicspool%DIN(npt)*mimicspool%MICr(npt)/(mimicspool%MICr(npt)+ mimicspool%MICk(npt)) 
          DINup_k = mimicspool%DIN(npt)*mimicspool%MICk(npt)/(mimicspool%MICr(npt)+ mimicspool%MICk(npt)) 

          upMICrC = mimicsbiome%MGE(1)*(LITmin(1)+SOMmin(1)) + mimicsbiome%MGE(2)*(LITmin(2))
          upMICrN = mimicsbiome%NUE(1)*(LITminN(1)+SOMminN(1)) + mimicsbiome%NUE(2)*(LITminN(2)) + DINup_r
          CNup_r = upMICrC/(upMICrN + 1e-10) !avoiding /0
          Overflow_r = upMICrC - upMICrN*min(mimicsbiome%CN_r(npt), CNup_r) 
          Nspill_r   = upMICrN - upMICrC/max(mimicsbiome%CN_r(npt), CNup_r) 
          ! Add Overflow_r and Overflow_k to output netCDF file. Units conversion occurs in subroutine mimics_caccum. -mdh 10/12/2020
          mimicsflux%Overflow_r(npt) = mimicsflux%Overflow_r(npt) + Overflow_r

          upMICkC = mimicsbiome%MGE(3)*(LITmin(3)+ SOMmin(2)) + mimicsbiome%MGE(4)*(LITmin(4))
          upMICkN = mimicsbiome%NUE(3)*(LITminN(3) + SOMminN(2)) + mimicsbiome%NUE(4)*LITminN(4) + DINup_k
          CNup_k = upMICkC/(upMICkN + 1e-10)
          Overflow_k = upMICkC - upMICkN*min(mimicsbiome%CN_k(npt), CNup_k) 
          Nspill_k = upMICkN - upMICkC/max(mimicsbiome%CN_k(npt), CNup_k)   
          ! Add Overflow_r and Overflow_k to output netCDF file. Units conversion occurs in subroutine mimics_caccum. -mdh 10/12/2020
          mimicsflux%Overflow_k(npt) = mimicsflux%Overflow_k(npt) + Overflow_k

          ! Divide total litter inputs (mgC/cm3/day) by number of hours in the daily timestep (NHOURSf)
          ! Optimization: Multiply by NHOURSfrac = 1/NHOURSf
          dLITm = mimicsflux%ClitInput(npt,metbc)*NHOURSfrac * (1.0-mimicsbiome%Fi(metbc)) - LITmin(1) - LITmin(3)
          !dMICr = mimicsbiome%MGE(1)*(LITmin(1)+ SOMmin(1)) + mimicsbiome%MGE(2)*(LITmin(2)) - sum(MICtrn(1:3))
          dMICr = upMICrC - sum(MICtrn(1:3)) - Overflow_r
          dSOMp = mimicsflux%ClitInput(npt,metbc)*NHOURSfrac * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4)- DEsorp 
          dLITs = mimicsflux%ClitInput(npt,struc)*NHOURSfrac * (1.0-mimicsbiome%Fi(struc)) - LITmin(2) - LITmin(4)
          !dMICk = mimicsbiome%MGE(3)*(LITmin(3)+ SOMmin(2)) + mimicsbiome%MGE(4)*(LITmin(4)) - sum(MICtrn(4:6))
          dMICk = upMICkC - sum(MICtrn(4:6))- Overflow_k
          dSOMc = mimicsflux%ClitInput(npt,struc)*NHOURSfrac * mimicsbiome%Fi(struc) + MICtrn(2) + MICtrn(5) - OXIDAT
          dSOMa = MICtrn(3) + MICtrn(6) + DEsorp + OXIDAT - SOMmin(1) - SOMmin(2)

          dLITmN = mimicsflux%NlitInput(npt,metbc)*NHOURSfrac * (1.0-mimicsbiome%Fi(metbc)) - LITminN(1) - LITminN(3)
          dMICrN = upMICrN - sum(MICtrnN(1:3)) - Nspill_r
          dSOMpN = mimicsflux%NlitInput(npt,metbc)*NHOURSfrac * mimicsbiome%Fi(metbc) + MICtrnN(1) + MICtrnN(4)- DEsorpN 
          dLITsN = mimicsflux%NlitInput(npt,struc)*NHOURSfrac * (1.0-mimicsbiome%Fi(struc)) - LITminN(2) - LITminN(4)
          dMICkN = upMICkN - sum(MICtrnN(4:6)) - Nspill_k
          dSOMcN = mimicsflux%NlitInput(npt,struc)*NHOURSfrac * mimicsbiome%Fi(struc) + MICtrnN(2) + MICtrnN(5) - OXIDATN
          dSOMaN = MICtrnN(3) + MICtrnN(6) + DEsorpN + OXIDATN - SOMminN(1) - SOMminN(2)
          dDIN = (1-mimicsbiome%NUE(1))*(LITminN(1)+SOMminN(1)) + (1-mimicsbiome%NUE(2))*(LITminN(2)) +  &
                 (1-mimicsbiome%NUE(3))*(LITminN(3)+SOMminN(2)) + (1-mimicsbiome%NUE(4))*(LITminN(4)) +  &
                 Nspill_r + Nspill_k - DINup_r - DINup_k
 
          !Sum daily heterotrophic respiration flux (mgC/cm3) 
          mimicsflux%Chresp(npt) = mimicsflux%Chresp(npt)  &
                                   + (1.0 - mimicsbiome%MGE(1)) * (LITmin(1) + SOMmin(1)) &
                                   + (1.0 - mimicsbiome%MGE(2)) * (LITmin(2))             &
                                   + (1.0 - mimicsbiome%MGE(3)) * (LITmin(3) + SOMmin(2)) &
                                   + (1.0 - mimicsbiome%MGE(4)) * (LITmin(4))             &
                                   + Overflow_r + Overflow_k

          !Sum daily C inputs to SOMp (mgC/cm3). -mdh 12/3/2018
          !ClitInput is initialized in mimics_delplant_CN
          mimicsflux%CSOMpInput(npt) = mimicsflux%CSOMpInput(npt)  &
              + mimicsflux%ClitInput(npt,metbc)*NHOURSfrac * mimicsbiome%Fi(metbc) + MICtrn(1) + MICtrn(4)

          ! Connect to CASA fluxes. Convert mgC/cm3 to gC/m2

          !Probably should not modify casa pools and fluxes here. Instead pass back values.
          !So these pools and fluxes can be set in CASA. 
          !Too many values to pass back. Work on this later. -mdh 6/10/2019.

          casaflux%Nlittermin(npt) = casaflux%Nlittermin(npt) &
                 + unitConv*((1-mimicsbiome%NUE(1))*LITminN(1) + (1-mimicsbiome%NUE(2))*LITminN(2) &
                 + (1-mimicsbiome%NUE(3))*LITminN(3) + (1-mimicsbiome%NUE(4))*LITminN(4))

          ! Correction 8/26/2019  
          ! casaflux%Nsmin(npt) = casaflux%Nsmin(npt) &
          !       + unitConv*((1-mimicsbiome%NUE(1))*SOMminN(1) + (1-mimicsbiome%NUE(3))*SOMminN(2) &
          !       + Nspill_r + Nspill_k) 
          ! casaflux%Nsimm(npt) = casaflux%Nsimm(npt) - unitConv*(DINup_r+DINup_k) 

          casaflux%Nsmin(npt) = casaflux%Nsmin(npt) &
                 + unitConv*((1-mimicsbiome%NUE(1))*SOMminN(1) + (1-mimicsbiome%NUE(3))*SOMminN(2))
          casaflux%Nsimm(npt) = casaflux%Nsimm(npt) - unitConv*(DINup_r + DINup_k - Nspill_r - Nspill_k)

          !Update C pools (mgC/cm3)
          mimicspool%LITm(npt) = mimicspool%LITm(npt) + dLITm
          mimicspool%LITs(npt) = mimicspool%LITs(npt) + dLITs
          mimicspool%MICr(npt) = mimicspool%MICr(npt) + dMICr
          mimicspool%MICk(npt) = mimicspool%MICk(npt) + dMICk
          mimicspool%SOMa(npt) = mimicspool%SOMa(npt) + dSOMa
          mimicspool%SOMc(npt) = mimicspool%SOMc(npt) + dSOMc
          mimicspool%SOMp(npt) = mimicspool%SOMp(npt) + dSOMp

          !Update N pools (mgN/cm3)
          mimicspool%LITmN(npt) = mimicspool%LITmN(npt) + dLITmN
          mimicspool%LITsN(npt) = mimicspool%LITsN(npt) + dLITsN
          mimicspool%MICrN(npt) = mimicspool%MICrN(npt) + dMICrN
          mimicspool%MICkN(npt) = mimicspool%MICkN(npt) + dMICkN
          mimicspool%SOMaN(npt) = mimicspool%SOMaN(npt) + dSOMaN
          mimicspool%SOMcN(npt) = mimicspool%SOMcN(npt) + dSOMcN
          mimicspool%SOMpN(npt) = mimicspool%SOMpN(npt) + dSOMpN
          mimicspool%DIN(npt) = mimicspool%DIN(npt) + dDIN

!     ! Code placed here temporarily for debugging hourly output. -mdh 2/3/2020
!     if (npt .eq. iptToSave_mimics) then 
!         ! Write to point file daily if mdaily==1, else write once a year on day 365. -mdh 5/14/2018
!         if ((mdaily == 1) .or. (idoy==365)) then
!             call WritePointMIMICS_CN(214, sPtFileNameMIMICS, npt, mp, ihr, idoy, &
!                 cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
!                 LITmin, MICtrn, SOMmin, DEsorp, OXIDAT, &
!                 dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil, &
!
!                 nleaf2met,nleaf2str,nroot2met,nroot2str,nwd2str,nwood2cwd, &
!                 LITminN, MICtrnN, SOMminN, DEsorpN, OXIDATN, &
!                 dLITmN, dLITsN, dMICrN, dMICkN, dSOMaN, dSOMcN, dSOMpN, dDIN, &
!                 DINup_r, DINup_k, upMICrC, upMICrN, upMICkC, upMICkN, & 
!                 Overflow_r, Overflow_k, Nspill_r, Nspill_k, Cbalance, Nbalance)
!         endif 
!     endif 

      end do ! End of 24-hour loop

      ! This is from casa_delsoil. Incorporate into this subroutine. -mdh 6/24/2019
      ! casaflux%Nsnet(npt)=casaflux%Nlittermin(npt)   &
      !                    + casaflux%Nsmin(npt) &
      !                    + casaflux%Nsimm(npt)
      ! casaflux%Nsnet = Net N mineralization (gN/m2/day), calculated below.
      ! casaflux%Nmindep = N deposition (gN/m2/day) from netcdf input file, set in casacnpdriver each day.
      ! casaflux%Nminfix = biome-specific N fixation rate (gN/m2/day), set in casa_readbiome.
      ! casaflux%Nminloss = portion of Nsnet lost (assumed to be N2O) (gN/m2/day), calculated below.
      ! casaflux%Nupland = N uptake by plants (gN/m2/day), calculated in casa_nuptake.
      ! casaflux%Nminleach = N leaching losses (gN/m2/day) not subtracted here since they were removed before 24-hour loop. -mdh 9/27/2020
      !
      ! casapool%dNsoilmindt(npt)= casaflux%Nsnet(npt)&
      !                          + casaflux%Nmindep(npt) 
      !                          + casaflux%Nminfix(npt)   &
      !                          - casaflux%Nminloss(npt)  &
      !                          - casaflux%Nupland(npt) 

      ! The Nsnet calculation below is from casa_delsoil.
      casaflux%Nsnet(npt) = casaflux%Nlittermin(npt) + casaflux%Nsmin(npt) + casaflux%Nsimm(npt)

      ! The Nminloss calculations below are from casa_delsoil.
      ! There is additional leaching according to CASA model. -mdh 8/5/2019
      if (casapool%Nsoilmin(npt) > 2.0 .AND. casamet%tsoilavg(npt) > 273.12) THEN
        casaflux%Nminloss(npt)   = casaflux%fNminloss(npt)  &
                                   * MAX(0.0,casaflux%Nsnet(npt))
      else
        ! ATTENTION: this code from casa_delsoil doesn't look correct.  
        ! I don't understand why the product was used in the calculation. 
        ! But changing these calculations in casa_delsoil messed up the 
        ! N cycle when CASACNP is run without MIMICS-CN. -mdh 8/12/2019
        casaflux%Nminloss(npt)   = casaflux%fNminloss(npt)  &
                                   * MAX(0.0,casaflux%Nsnet(npt)) &
                                   * MAX(0.0,casapool%Nsoilmin(npt)/2.0)
      endif

      ! dNsoilmindt is used in mimics_cncycle to update casapool%Nsoilmin. 
      ! The change in soil mineral N for the day includes microbial uptake (Nsimm), 
      ! and N mineralization (Nlittermin+Nsmin)
      casapool%dNsoilmindt(npt) = unitConv*(mimicspool%DIN(npt) - DINstart)

      ! Add components of dNsoilmindt not included in MIMICS hourly calculations above. 
      ! This is everything but Nsnet=Nlittermin+Nsmin+Nsimm and mineral N leaching
      ! casapool%dNsoilmindt(npt) is removed from casapool%Nsoilmin(npt) in subroutine mimics_cncycle.
      casapool%dNsoilmindt(npt) = casapool%dNsoilmindt(npt) &
                                  + casaflux%Nmindep(npt)   &
                                  + casaflux%Nminfix(npt)   &
                                  - casaflux%Nminloss(npt)  &
                                  - casaflux%Nupland(npt) 
      ! Temporarily record Rh in the casaflux variable for model evaluation. -mdh 8/5/2019. Remove later.
      casaflux%Crsoil = mimicsflux%Chresp(npt) * unitConv + cwd2co2

      CorgSum2 = mimicspool%LITm(npt) + mimicspool%LITs(npt) + &
                 mimicspool%SOMa(npt) + mimicspool%SOMc(npt) + mimicspool%SOMp(npt)
      NorgSum2 = mimicspool%LITmN(npt) + mimicspool%LITsN(npt) + &
                 mimicspool%SOMaN(npt) + mimicspool%SOMcN(npt) + mimicspool%SOMpN(npt)

      ! Compute C&N balance in terms of MIMICS units (N leaching losses already considered in change in DIN

      ! CorgSum2 = CorgSum1 + gains - losses
      Cbalance = -1.0*(CorgSum1 - CorgSum2 - mimicsflux%Chresp(npt) &
                  + mimicsflux%ClitInput(npt,metbc)+mimicsflux%ClitInput(npt,struc))

      ! NorgSum2 + mimicspool%DIN(npt) = NorgSum1 + DINstart + gains - losses
      Nbalance = -1.0*(NorgSum1 - NorgSum2 + DINstart - mimicspool%DIN(npt) &
                  + mimicsflux%NlitInput(npt,metbc)+mimicsflux%NlitInput(npt,struc)) 

!     if (.TRUE.) then
      if (.FALSE.) then
    
          write(*,*)
          write(*,*) 'N balance > 0.0 means unexpected gain of N from start to end of day'
          write(*,'(a13,f10.6)') 'Nbalance = ', Nbalance
          write(*,'(a13,f10.6)') 'NorgSum2 = ', NorgSum2
          write(*,'(a13,f10.6)') 'NorgSum1 = ', NorgSum1
          write(*,'(a13,f10.6)') 'dNorg = ', NorgSum2 - NorgSum1
          write(*,'(a13,f10.6)') 'DINend = ', mimicspool%DIN(npt)
          write(*,'(a13,f10.6)') 'DINstart = ', DINstart
          write(*,'(a13,f10.6)') 'dDIN = ', mimicspool%DIN(npt) - DINstart
          write(*,'(a13,f10.6)') 'litNinput = ', mimicsflux%NlitInput(npt,metbc) + mimicsflux%NlitInput(npt,struc)
          write(*,'(a13,f10.6)') 'Cbalance = ', Cbalance
          write(*,'(a13,f10.6)') 'CorgSum2 = ', CorgSum2
          write(*,'(a13,f10.6)') 'CorgSum1 = ', CorgSum1
          write(*,'(a13,f10.6)') 'dCorg = ', CorgSum2 - CorgSum1
          write(*,'(a13,f10.6)') 'litCinput = ', mimicsflux%ClitInput(npt,metbc) + mimicsflux%ClitInput(npt,struc)
          write(*,'(a13,f10.6)') 'Chresp = ', mimicsflux%Chresp(npt)
    
      endif

!     if (casamet%ijgcm(npt) .eq. iptToSave_mimics) then 
      if (npt .eq. iptToSave_mimics) then 
          ! Write to point file daily if mdaily==1, else write once a year on day 365. -mdh 5/14/2018
          if ((mdaily == 1) .or. (idoy==365)) then
              call WritePointMIMICS_CN(214, sPtFileNameMIMICS, npt, mp, iYrCnt, idoy, &
                  cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
                  LITmin, MICtrn, SOMmin, DEsorp, OXIDAT, &
                  dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil, &

                  nleaf2met,nleaf2str,nroot2met,nroot2str,nwd2str,nwood2cwd, &
                  LITminN, MICtrnN, SOMminN, DEsorpN, OXIDATN, &
                  dLITmN, dLITsN, dMICrN, dMICkN, dSOMaN, dSOMcN, dSOMpN, dDIN, &
                  DINup_r, DINup_k, upMICrC, upMICrN, upMICkC, upMICkN, & 
                  mimicsflux%Overflow_r, mimicsflux%Overflow_k, Nspill_r, Nspill_k, Cbalance, Nbalance)
          endif 
      endif 
       
  ENDIF  !IF(casamet%iveg2(npt) /= icewater)

  end do

END SUBROUTINE mimics_soil_reverseMM_CN
!--------------------------------------------------------------------------------

END MODULE mimics_cycle_module
