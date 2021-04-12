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
! Purpose: subroutines for calculating carbon, nitrogen, phosphorus cycle 
!          including plant growth
!
! Called from: biogeochem (mostly) or casa_xnp
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed by Yingping Wang (Wang et al., BG, 2011)
!          Current version uses fixed phenology.
!
!
! ==============================================================================
! casa_cnp.f90
!
! This module contains the following subroutines:
!   casa_xnp
!   casa_allocation
!   casa_rplant
!   casa_xrateplant,        
!   casa_xratesoil
!   casa_coeffplant,        
!   casa_coeffsoil
!   casa_delplant,          
!   casa_delsoil
!   avgsoil
!   casa_xkN (moved to casa_nlim.f90 on 11/11/2019 -mdh)
!   casa_nuptake,           
!   casa_puptake
!   casa_Nrequire,          
!   casa_Prequire
!   casa_cnpcycle
!   casa_poolzero
!   casa_cnpbal
!   casa_ndummy
!   phenology

MODULE casa_cnp_module
USE define_dimensions
USE define_types
USE casadimension
USE casaparm
USE casavariable
USE phenvariable
IMPLICIT NONE
CONTAINS


SUBROUTINE casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xnplimit
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xNPuptake
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  INTEGER :: np
  REAL(r_2), DIMENSION(mp)        :: xnlimit,xplimit
  REAL(r_2), DIMENSION(mp)        :: xncleaf,xpcleaf
  REAL(r_2), DIMENSION(mp)        :: xnCnpp,xpCnpp
  REAL(r_2), DIMENSION(mp,mplant) :: Nreqmax, Nreqmin, NtransPtoP
  REAL(r_2), DIMENSION(mp)        :: totNreqmax,totNreqmin
  REAL(r_2), DIMENSION(mp)        :: xNuptake,xPuptake
  REAL(r_2), DIMENSION(mp,mplant) :: Preqmax, Preqmin, PtransPtoP
  REAL(r_2), DIMENSION(mp)        :: totPreqmax,totPreqmin

  xnlimit  = 1.0
  xplimit  = 1.0
  xnplimit = 1.0
  casaflux%fracClabile(:) = 0.0

  !print *, 'xnp:icycle', icycle
  SELECT CASE(icycle)
  CASE(2)
    WHERE(casamet%iveg2/=icewater) 
      xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
      !xnlimit(:) = xncleaf(:)/(xncleaf(:)+casabiome%KminN(veg%iveg(:)))
      xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01)
      xplimit(:) = 1.0
      xnplimit(:) =min(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
    ENDWHERE
  CASE(3)
    WHERE(casamet%iveg2/=icewater) 
      xncleaf(:) = casapool%nplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
      xnlimit(:) = xncleaf(:)/(xncleaf(:)+0.01)
      xpcleaf(:) = casapool%pplant(:,leaf)/(casapool%cplant(:,leaf)+1.0e-10)
      !xplimit(:) = xpcleaf(:)/(xpcleaf(:)+casabiome%Kuplabp(veg%iveg(:)))
      xplimit(:) = xpcleaf(:)/(xpcleaf(:)+0.0006)
      !xnplimit(:) = min(1.0,casabiome%Kuptake(veg%iveg(:))*min(xnlimit(:),xplimit(:)))
      xnplimit(:) =min(xnlimit(:),xplimit(:)) * casabiome%xnpmax(veg%iveg(:))
    ENDWHERE
  END SELECT

  ! now check if soil nutrient supply can meet the plant uptake,
  ! otherwise reduce NPP
  xNuptake = 1.0
  xPuptake = 1.0 

  IF(icycle >1) THEN
    Nreqmin(:,:)    = 0.0
    Nreqmax(:,:)    = 0.0
    NtransPtoP(:,:) = 0.0
    totNreqmax = 0.0
    totNreqmin = 0.0
    xNuptake   = 1.0
    xnCnpp = max(0.0,casaflux%Cnpp)
    call casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
    DO np=1,mp
      IF(casamet%iveg2(np)/=icewater) THEN 
        totNreqmax(np) = Nreqmax(np,leaf)+Nreqmax(np,wood)+Nreqmax(np,froot)
        totNreqmin(np) = Nreqmin(np,leaf)+Nreqmin(np,wood)+Nreqmin(np,froot)
        xNuptake(np)   = MAX(0.0,MIN(1.0,casapool%Nsoilmin(np) &
                                         /(totNreqmin(np)*deltpool+1.0e-10)))
      ENDIF
    ENDDO
  ENDIF
  IF(icycle >2) THEN
    Preqmin(:,:)       = 0.0
    Preqmax(:,:)       = 0.0
    PtransPtoP(:,:)    = 0.0
    totPreqmax = 0.0
    totPreqmin = 0.0
    xPuptake   = 1.0
    xpCnpp = max(0.0,casaflux%Cnpp)
    call casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                       casabiome,casapool,casaflux,casamet)
    DO np=1,mp
      IF(casamet%iveg2(np)/=icewater) THEN 
        totPreqmax(np) = Preqmax(np,leaf)+Preqmax(np,wood)+Preqmax(np,froot)
        totPreqmin(np) = Preqmin(np,leaf)+Preqmin(np,wood)+Preqmin(np,froot)
        xPuptake(np)   = MAX(0.0,MIN(1.0,casapool%psoillab(np) &
                                         /(totPreqmin(np)*deltpool+1.0e-10)))
      ENDIF
    ENDDO
  ENDIF

  xnplimit(:)  = 1.0
  xNPuptake(:)     = min(xnuptake(:), xpuptake(:))
  
  do np =1, mp
     if(casamet%iveg2(np)/=icewater.and.casaflux%Cnpp(np) > 0.0.and.xNPuptake(np) < 1.0) then
!       write(*,*) 'NPP BEFORE REDUCTION casaflux%Cnpp(',np,')=', casaflux%Cnpp(np)
        casaflux%fracClabile(np) =min(1.0,max(0.0,(1.0- xNPuptake(np)))) * max(0.0,casaflux%Cnpp(np))/(casaflux%Cgpp(np) +1.0e-10)
        casaflux%Cnpp(np)    = casaflux%Cnpp(np) - casaflux%fracClabile(np) * casaflux%Cgpp(np)
!       write(*,*) 'NPP REDUCED casaflux%Cnpp(',np,')=', casaflux%Cnpp(np)
     endif
  enddo

!  casaflux%Cnpp(:) = xNPuptake(:) * xnplimit(:) * casaflux%Cnpp(:)
!  write(*,911) xNuptake(1), totNreqmin(1), totNreqmax(1), casapool%Nsoilmin(1)
!911  format('xnp: ',10(f8.3,2x))

END SUBROUTINE casa_xnp


SUBROUTINE casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen)
! compute fraction of net photosynthate allocated to leaf, wood and froot
!
! inputs
!   moistavg(mp)           as an argument (volumetric liquid fraction)
!   frznmoistavg(mp)       as an argument (volumetric frozen fraction)
!   tsoilavg(mp)           as an argument (K)
!   btran(mp)              as an argument (dimensionless)
! outputs:
!   fracCalloc(mp,mplant1)  
!
! modified Piere's alocation scheme 
! input: leaf stage 
!        leaf area

  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  ! local variables
  !INTEGER :: npt,ns,is,iv
  REAL(r_2), DIMENSION(mp,mplant) :: fracCallocx
  !REAL(r_2), DIMENSION(mp,mplant) :: delc
  !REAL(r_2), DIMENSION(mp)        :: ctotal
  REAL(r_2), DIMENSION(mp)        :: xLalloc,xwsalloc,xTalloc
  REAL(r_2), DIMENSION(mp)        :: xWorNalloc,xNalloc,xWalloc
  REAL(r_2), DIMENSION(mp)        :: totfracCalloc                 

  ! initlization
  casaflux%fracCalloc  = 0.0
  casaflux%fracClabile = 0.0
  fracCallocx = 0.0

  SELECT CASE (LALLOC)

  CASE(2)   
    ! calculate the allocation coefficients
    call casa_wolf(veg,casabiome,casaflux,casapool,casamet)

  CASE(1)   ! dynamic allocation
    WHERE(casamet%iveg2/=icewater) 
      xLalloc(:) = min(1.0,max(0.0,exp(-0.5*casamet%glai(:))))   ! L limiting
      ! Pseudo-nutrient limitation calculation 
      WHERE(casamet%tsoilavg > 0.0) 
        xwsalloc(:) = min( max(casamet%moistavg(:)-soil%swilt(:),0.0) &
                         /(soil%sfc(:)-soil%swilt(:)), 1.0 )
      ELSE WHERE
        xwsalloc(:) = 0.01
      END WHERE
      xTalloc(:)    = min(1.0,max(0.0,Q10alloc** &
                      ((casamet%tsoilavg(:)-TkzeroC-30.0)/10.0) )) !T limiting
      xNalloc(:)    = min(1.0,max(0.0,xwsalloc(:)*xTalloc(:)))     !N limiting
      xWalloc(:)    = min(1.0,max(0.0,casamet%btran(:)))           !W limiting
      xWorNalloc(:) = min(xWalloc(:),xNalloc(:))
      WHERE(casamet%lnonwood==0) 
        casaflux%fracCalloc(:,FROOT) = R0 * 3.0 * xLalloc(:) &
                                     / (xLalloc(:)+ 2.0*xWorNalloc(:))
        casaflux%fracCalloc(:,WOOD)  = S0 * 3.0 * xWorNalloc(:) &
                                     / (2.0*xLalloc(:)+ xWorNalloc(:))
        casaflux%fracCalloc(:,LEAF)  = 1.0 - casaflux%fracCalloc(:,FROOT) &
                                     - casaflux%fracCalloc(:,WOOD)
      ELSE WHERE
        casaflux%fracCalloc(:,FROOT) = R0 * 3.0 * xLalloc(:) &
                                     / (xLalloc(:)+2.0*xWorNalloc(:))
        casaflux%fracCalloc(:,WOOD)  = 0.0
        casaflux%fracCalloc(:,LEAF)  = 1.0 - casaflux%fracCalloc(:,FROOT)
      END WHERE
    END WHERE
  CASE (0)   ! fixed allocation
    casaflux%fracCalloc(:,:) = casabiome%fracnpptop(veg%iveg(:),:)      
  END SELECT


  ! during leaf growth phase 0 or 3, no carbon is allocated to leaf, 
  ! during maximal leaf growth phase, all C is allocated to leaf
  ! during steady growth period, C allocation is estimated in such 
  ! a way that approach the allometric relationship
  ! the relationships are:(all pools in g C/m2)
  ! for forests
  !   fineroot/totalC C=0.3192-0.0485*(totalC)^0.1755, see mokany et al. (2003)
  !   fineroot = ratiofrootleaf*cleaf
  ! for grassland
  !   root=ratiofinerootleaf*cleaf

  WHERE(casamet%iveg2/=icewater) 
    WHERE(phen%phase==0) 
      casaflux%fracCalloc(:,leaf)  = 0.0
      casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                                   /(casaflux%fracCalloc(:,froot) &
                                     +casaflux%fracCalloc(:,wood))
      casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
    END WHERE 

    WHERE(phen%phase==1)          
      casaflux%fracCalloc(:,leaf)  = 0.8
      WHERE(casamet%lnonwood==0)  !woodland or forest
        casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
        casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
      ELSEWHERE !grassland
        casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
      ENDWHERE
    END WHERE

    WHERE(phen%phase==3) 
!      casaflux%fracClabile(:)  = casaflux%fracCalloc(:,leaf)
      casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,wood) 
      casaflux%fracCalloc(:,leaf)  = 0.0
    ENDWHERE

  ! IF Prognostic LAI reached glaimax, no C is allocated to leaf
  ! Q.Zhang 17/03/2011
    WHERE(casamet%glai(:)>=casabiome%glaimax(veg%iveg(:)))
      casaflux%fracCalloc(:,leaf)  = 0.0
      casaflux%fracCalloc(:,froot) =  casaflux%fracCalloc(:,froot) &
                                   /(casaflux%fracCalloc(:,froot) &
                                     +casaflux%fracCalloc(:,wood))
      casaflux%fracCalloc(:,wood)  = 1.0 -casaflux%fracCalloc(:,froot)
    ENDWHERE

    WHERE(casamet%glai(:)<casabiome%glaimin(veg%iveg(:)))
      casaflux%fracCalloc(:,leaf)  = 0.8
      WHERE(casamet%lnonwood==0)  !woodland or forest
        casaflux%fracCalloc(:,froot) = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
        casaflux%fracCalloc(:,wood)  = 0.5*(1.0-casaflux%fracCalloc(:,leaf))
      ELSEWHERE !grassland
        casaflux%fracCalloc(:,froot) = 1.0-casaflux%fracCalloc(:,leaf)
      ENDWHERE
    ENDWHERE
    ! added in for negative NPP and one of biomass pool being zero ypw 27/jan/2014
    WHERE(casaflux%Cnpp<0.0) 
       casaflux%fracCalloc(:,leaf)  = casaflux%Crmplant(:,leaf)/sum(casaflux%Crmplant,2)   
       casaflux%fracCalloc(:,wood)  = casaflux%Crmplant(:,wood)/sum(casaflux%Crmplant,2)   
       casaflux%fracCalloc(:,froot) = casaflux%Crmplant(:,froot)/sum(casaflux%Crmplant,2)   
    ENDWHERE    

  ENDWHERE
  ! normalization the allocation fraction to ensure they sum up to 1
  totfracCalloc(:) = sum(casaflux%fracCalloc(:,:),2)
  casaflux%fracCalloc(:,leaf) = casaflux%fracCalloc(:,leaf)/totfracCalloc(:)
  casaflux%fracCalloc(:,wood) = casaflux%fracCalloc(:,wood)/totfracCalloc(:)
  casaflux%fracCalloc(:,froot) = casaflux%fracCalloc(:,froot)/totfracCalloc(:)

END SUBROUTINE casa_allocation  

SUBROUTINE casa_wolf(veg,casabiome,casaflux,casapool,casamet)
   ! carbon allocation based on
   ! Wolf,Field and Berry, 2011. Ecological Applications, p1546-1556
   ! Wolf et al. 2011. Global Biogeochemical Cycles, 25, GB3015, doi:10.1019/2010GB003917
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(IN) :: casabiome
  TYPE (casa_met),            INTENT(IN) :: casamet
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux

   real, parameter :: wolf_alpha1=6.22
   real, parameter :: wolf_beta=-1.33
   real, parameter :: wolf_c1=-wolf_alpha1/(1+wolf_beta)
   real, parameter :: wolf_c2=1.0/(1.0+wolf_beta)

   !real  totleaf,totwood,totcroot,totfroot,totnpp
   !real  fracleaf,fracwood,fraccroot,fracfroot
   !
   ! local variables
   integer   npt
   real(r_2), dimension(mp)  ::  totbmdm,ntree,nppdm
   real(r_2), dimension(mp)  ::  gleaf,gwood,gcroot,gfroot,gtot
   !
   ! input
   !  totleaf, totwood, totcroot, totfroot :    g C m-2
   !  totnpp:                                   g C m-2 d-1
   ! output
   !  fracleaf,fracwood, fraccroot, fracfroot:  fractions
   !

    do npt=1,mp
       IF(casamet%iveg2(npt)==forest .and. casaflux%Cnpp(npt)>0.0001) THEN  !forest types
          totbmdm(npt) = sum(casapool%cplant(npt,:)) *10.0 / fracCbiomass      !10.0 for convert gc/m2 to kg/ha
          totbmdm(npt) = max(30000.0, totbmdm(npt))
          ! calculate tree stocking density
           ntree(npt) = 10**(wolf_c1+wolf_c2*log10(totbmdm(npt)))   ! tree ha-1, based on eqn (4) of Wolf et al. 2011, GBC
           ntree(npt) = min(200000.0,ntree(npt))
           ! changed by ypw 23/april/2012 to avoid negative npp
           nppdm(npt)  = (abs(casaflux%Cnpp(npt)) *365.0*0.001/fracCbiomass)/(0.0001*ntree(npt))  ! in kg dm tree-1 yr-1

           gleaf(npt)  = 0.156*(nppdm(npt)**1.106)     ! Figure 2a of Wolf, Field and Berry (2011)
           gwood(npt)  = 0.232*(nppdm(npt)**1.165)     ! Figure 2b of Wolf, Field and Berry (2011)
           gcroot(npt) = 0.0348*(nppdm(npt)**1.310)    ! Figure 2d of Wolf, Field and Berry (2011)
           gfroot(npt) = 0.247*(nppdm(npt)**0.987)     ! Figure 2c of Wolf, Field and Berry (2011)
           gtot(npt)   = gleaf(npt) + gwood(npt) + gcroot(npt) + gfroot(npt)

           casaflux%fracCalloc(npt,leaf)  = gleaf(npt)/gtot(npt)
           casaflux%fracCalloc(npt,wood)  = gwood(npt)/gtot(npt)
           casaflux%fracCalloc(npt,froot) = (gcroot(npt)+gfroot(npt))/gtot(npt)

!        write(87,*) 'allocation = ',npt,casamet%iveg2(npt), totbmdm(npt),ntree(npt),nppdm(npt),casaflux%fracCalloc(npt,:)

        ELSE                ! other types
           casaflux%fracCalloc(npt,:) = casabiome%fracnpptop(veg%iveg(npt),:)
        ENDIF
    enddo

END SUBROUTINE casa_wolf


SUBROUTINE casa_rplant(veg,casabiome,casapool,casaflux,casamet)
! maintenance respiration of woody tisse and fineroots 
! see Sitch et al. (2003), GCB, reqn (23)

  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  !INTEGER :: npt

  real(r_2), dimension(mp)        :: Ygrow        ! growth efficiency Q.Zhang 22/02/2011
  real(r_2), dimension(mp,mplant) :: ratioPNplant ! Q.Zhang 22/02/2011
  real(r_2), dimension(mp)        :: delcrmwood,delcrmfroot    ! reduction in wood and root respiration when NPP <0.0   

  ratioPNplant = 0.0
  Ygrow        = 0.0

  WHERE(casapool%Nplant>0.0)
    ratioPNplant = casapool%Pplant/(casapool%Nplant+ 1.0e-10)
  ENDWHERE

  Ygrow(:) = 0.65+0.2*ratioPNplant(:,leaf)/(ratioPNplant(:,leaf)+1.0/15.0)

  casaflux%Crmplant(:,wood) = 0.0
  casaflux%Crmplant(:,froot) = 0.0
  delcrmwood   = 0.0
  delcrmfroot  = 0.0
  casaflux%Crgplant = 0.0
  casaflux%Clabloss = 0.0

  WHERE(casamet%iveg2/=icewater) 
    WHERE(casamet%tairk >250.0) 
      WHERE(casapool%cplant(:,wood)>1.0e-6)
      casaflux%Crmplant(:,wood)  = casabiome%rmplant(veg%iveg(:),wood) &
                                 * casapool%nplant(:,wood)             &
                                 * exp(308.56*(1.0/56.02-1.0           &
                                 / (casamet%tairk(:)+46.02-tkzeroc)))
      ENDWHERE
      casaflux%Clabloss(:)  =  casabiome%kclabrate(veg%iveg(:)) &
                            * max(0.0,casapool%Clabile(:))      &
                            * exp(308.56*(1.0/56.02-1.0         &
                            / (casamet%tairk(:)+46.02-tkzeroc)))
    ENDWHERE
    WHERE(casamet%tsoilavg >250.0.and.casapool%cplant(:,froot)>1.0e-6) 
      casaflux%Crmplant(:,froot) = casabiome%rmplant(veg%iveg(:),froot) &
                                 * casapool%nplant(:,froot)             &
                                 * exp(308.56*(1.0/56.02-1.0            &
                                 / (casamet%tsoilavg(:)+46.02-tkzeroc)))
    ENDWHERE
!    casaflux%Crmplant(:,leaf) = casaflux%Crmplant(:,leaf) + casaflux%Clabloss(:)

    WHERE((casaflux%Cgpp-SUM(casaflux%Crmplant,2))>0.0)
    !casaflux%crgplant(:)  = 0.25* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%Crmplant(:,:),2))
    ! Growth efficiency correlated to leaf N:P ratio. Q.Zhang @ 22/02/2011
      casaflux%crgplant(:)  = (1.0-Ygrow(:))* max(0.0,casaflux%Cgpp(:)-SUM(casaflux%Crmplant(:,:),2))

    ELSEWHERE
      casaflux%crgplant(:) = 0.0
    ENDWHERE

    ! changes made by yp wang 5 april 2013
    casaflux%Cnpp(:) = casaflux%Cgpp(:)-SUM(casaflux%Crmplant(:,:),2) - casaflux%crgplant(:) 
    !! Added plant respiration calculation here because the Crp variable is not set anywhere else (-MDH 7/21/2014)
    casaflux%Crp(:) =  SUM(casaflux%Crmplant(:,:),2) + casaflux%crgplant(:)
  ENDWHERE

!  write(57,*) 'rplant', casaflux%Cgpp(243),casaflux%Crmplant(243,:),casamet%tsoilavg(243),casamet%tairk(243)
!$$$$$$$$$$$$$$$$$$$$$$
!    WHERE(casaflux%Cnpp < 0.0)
!        delcrmwood(:)  = casaflux%Cnpp(:) * casaflux%Crmplant(:,wood)/ (1.0e-10+ casaflux%Crmplant(:,wood) + casaflux%Crmplant(:,froot))
!        delcrmfroot(:) = casaflux%Cnpp(:) * casaflux%Crmplant(:,froot)/(1.0e-10+ casaflux%Crmplant(:,wood) + casaflux%Crmplant(:,froot))
!        casaflux%Crmplant(:,wood) = casaflux%Crmplant(:,wood) + delcrmwood(:)
!        casaflux%Crmplant(:,froot) = casaflux%Crmplant(:,froot) + delcrmfroot(:)
!        casaflux%Cnpp(:) = casaflux%Cnpp(:) -delcrmwood(:)-delcrmfroot(:)
!    ENDWHERE
!$$$$$$$$$$$$$$$$$$$$$$$
!  ENDWHERE

!  print *, 'calling rplant',veg%iveg(1),casamet%tairk(1)
!,tkzeroc,casapool%nplant(1,:),casaflux%Crmplant(1,:),casaflux%crgplant(1)
!   npt =290 
!   write(77,*) npt,casaflux%Cgpp(npt),casaflux%Cnpp(npt), casaflux%Crmplant(npt,:),casaflux%Crgplant(npt),casapool%cplant(npt,:),casapool%nplant(npt,:)
!701  format('point ',i6,100(f10.3,2x))
END SUBROUTINE casa_rplant


SUBROUTINE casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                           casamet,phen)
! use xleafcold and xleafdry to account for
! cold and drought stress on death rate of leaf
! inputs:
!     ivt(mp)  :       biome type
!     phase(mp):       leaf growth stage
!     tairk(mp)    :   air temperature in K
! outputs
!     xkleafcold(mp):  cold stress induced leaf senescence rate (1/day)
!     xkleafdry(mp):   drought-induced leaf senescence rate (1/day)
!     xkleaf(mp):      set to 0.0 during maximal leaf growth phase

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleafcold
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleafdry
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xkleaf
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mp)              :: xcoldleaf
  INTEGER :: npt

  xkleafcold(:) = 0.0
  xkleafdry(:)  = 0.0
  xkleaf(:)     = 1.0        

  ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
  DO npt=1,mp
  IF(casamet%iveg2(npt)/=icewater) THEN
  !    following the formulation of Arora (2005) on the 
  !    effect of cold or drought stress on leaf litter fall
  !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
    IF(casamet%tairk(npt)>=phen%TKshed(veg%iveg(npt))) THEN
      xcoldleaf(npt) = 1.0
    ELSE
      IF(casamet%tairk(npt)<=(phen%TKshed(veg%iveg(npt))-5.0)) THEN
        xcoldleaf(npt)=0.0
      ELSE
        ! Gordon Bonan discovered missing parentheses in the following equation. -mdh 5/14/2018
        !xcoldleaf(npt) = (casamet%tairk(npt)-phen%TKshed(veg%iveg(npt))-5.0)/5.0
        xcoldleaf(npt) = (casamet%tairk(npt)-(phen%TKshed(veg%iveg(npt))-5.0))/5.0
      ENDIF
    ENDIF
    xcoldleaf(npt) = min(1.0,max(0.0,xcoldleaf(npt)))
    xkleafcold(npt) = casabiome%xkleafcoldmax(veg%iveg(npt)) &
                    * (1.0-xcoldleaf(npt)) &
                    ** casabiome%xkleafcoldexp(veg%iveg(npt))
    xkleafdry(npt)  = casabiome%xkleafdrymax(veg%iveg(npt)) &
                    * (1.0-casamet%btran(npt))&
                    ** casabiome%xkleafdryexp(veg%iveg(npt))
    IF (phen%phase(npt)==1) xkleaf(npt)= 0.0
  END IF
  END DO

!  WHERE(casamet%iveg2/=icewater) 
!  !    following the formulation of Arora (2005) on the 
!  !    effect of cold or drought stress on leaf litter fall
!  !    calculate cold stress (eqn (18), Arora 2005, GCB 11:39-59)
!    WHERE(casamet%tairk(:)>=phen%TKshed(veg%iveg(:))) 
!      xcoldleaf(:) = 1.0
!    ELSEWHERE 
!      WHERE(casamet%tairk(:)<=(phen%TKshed(veg%iveg(:))-5.0)) 
!        xcoldleaf(:)=0.0
!      ELSEWHERE
!        xcoldleaf(:) = (casamet%tairk(:)-phen%TKshed(veg%iveg(:))-5.0)/5.0
!      ENDWHERE
!    ENDWHERE
!    xcoldleaf(:) = min(1.0,max(0.0,xcoldleaf(:)))
!    xkleafcold(:) = casabiome%xkleafcoldmax(veg%iveg(:)) * (1.0-xcoldleaf(:)) &
!                 ** casabiome%xkleafcoldexp(veg%iveg(:))
!    xkleafdry(:) = casabiome%xkleafdrymax(veg%iveg(:))*(1.0-casamet%btran(:))&
!                 ** casabiome%xkleafdryexp(veg%iveg(:))
!    WHERE(phen%phase(:)==1) xkleaf(:)= 0.0  
!  ENDWHERE 

END SUBROUTINE casa_xrateplant


SUBROUTINE casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
!  to account for effects of T and W on litter decomposition: xklitter, xksoil
!  inputs:
!     ivt(mp)  :        biome type
!     tsoilavg(mp):     soil temperature in K
!     moistavg(mp):     volumetric soil liquid moisture
!     frznmoistavg(mp): volumetric soil frozen moisture
!
!  outputs
!     xk(mp):          modifier of soil litter decomposition rate (dimensionless)
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(OUT) :: xklitter,xksoil
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome

  ! local variables
  !INTEGER np         
  REAL(r_2), parameter :: wfpscoefa=0.55   ! Kelly et al. (2000) JGR, Figure 2b), optimal wfps 
  REAL(r_2), parameter :: wfpscoefb=1.70   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefc=-0.007 ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefd=3.22   ! Kelly et al. (2000) JGR, Figure 2b)
  REAL(r_2), parameter :: wfpscoefe=6.6481 ! =wfpscoefd*(wfpscoefb-wfpscoefa)/(wfpscoefa-wfpscoefc)
! Kirschbaum function parameters
  REAL(r_2), parameter :: xkalpha=-3.764   ! Kirschbaum (1995, SBB)
  REAL(r_2), parameter :: xkbeta=0.204
  REAL(r_2), parameter :: xktoptc=36.9
  REAL(r_2), DIMENSION(mp)       :: xkwater,xktemp
  REAL(r_2), DIMENSION(mp)       :: fwps,tsavg
!,tsurfavg  !!, msurfavg
  INTEGER :: npt

  xklitter(:) = 1.0
  xksoil(:)   = 1.0
  fwps(:)     =  min(1.0, casamet%moistavg(:)/soil%ssat(:))
  tsavg(:)    =  casamet%tsoilavg(:) 

  ! BP changed the WHERE construct to DO-IF for Mk3L (jun2010)
  DO npt=1,mp
  IF(casamet%iveg2(npt)/=icewater) THEN
    xktemp(npt)  = casabiome%q10soil(veg%iveg(npt))**(0.1*(tsavg(npt)-TKzeroC-35.0))
    xkwater(npt) = ((fwps(npt)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
               * ((fwps(npt)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
    IF (veg%iveg(npt) == cropland .OR. veg%iveg(npt) == croplnd2) &
               xkwater(npt)=1.0
    xklitter(npt) = casabiome%xkoptlitter(veg%iveg(npt)) * xktemp(npt) * xkwater(npt)
    xksoil(npt)   = casabiome%xkoptsoil(veg%iveg(npt))   * xktemp(npt) * xkwater(npt)

    casapool%thetaLiq(npt) = fwps(npt)
    casapool%fT(npt) = xktemp(npt)
    casapool%fW(npt) = xkwater(npt)

    !print *, 'xratesoil: ', npt, casamet%moistavg(npt),soil%ssat(npt),casamet%tsoilavg(npt),xksoil(npt)
    !print *, 'xklitter: ', npt, tsavg(npt), xklitter(npt), xktemp(npt), xkwater(npt)

  END IF
  END DO
!  WHERE(casamet%iveg2/=icewater)  
!!    ! Kirschbaum function
!!    xktemp(:) = exp(xkalpha + xkbeta*(tsavg(:)-TKzeroC) &
!!              * (1.0-0.5*(tsavg(:)-TKzeroc)/xktoptc))
!    ! add by ypwang on 3/april/2009
!!    xktemp(:) = xkoptcoeff(veg%iveg(:))*exp(xkbeta*(tsavg(:)-TKzeroC-xktoptc))
!    xktemp(:)  = q10soil**(0.1*(tsavg(:)-TKzeroC-35.0))
!    xkwater(:) = ((fwps(:)-wfpscoefb)/(wfpscoefa-wfpscoefb))**wfpscoefe    &
!               * ((fwps(:)-wfpscoefc)/(wfpscoefa-wfpscoefc))**wfpscoefd
!    WHERE(veg%iveg==12)
!      xkwater(:)=1.0
!    ENDWHERE
!    xklitter(:) = xkoptlitter(veg%iveg(:)) * xktemp(:) * xkwater(:)
!    xksoil(:)   = xkoptsoil(veg%iveg(:))   * xktemp(:) * xkwater(:)
!  ENDWHERE 


END SUBROUTINE casa_xratesoil

SUBROUTINE casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                           casaflux,casamet)
! calculate the plant litter fall rate, litter fall and sOM decomposition rate (1/day)
! and the transfer coefficients between different pools
! NOTE: this subroutine is not called when MIMICS or CORPSE are the SOM model.
!       See mimics_coeffplant.
!
! inputs:
!     xkleafcold(mp):  cold stress induced leaf senescence rate (1/day)
!     xkleafdry(mp):   drought-induced leaf senescence rate (1/day)
!     xkleaf(mp):      set to 0.0 during maximal leaf growth phase
!
! outputs:
!     kplant(mp,mplant):        senescence rate of plant pool (1/day)
!     fromPtoL(mp,mlitter,mplant): fraction of senesced plant biomass to litter pool (fraction)

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(IN)    :: xkleafcold,xkleafdry,xkleaf
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  !REAL(r_2), DIMENSION(mp)  :: xk
  REAL(r_2), DIMENSION(mp,mplant)         :: ratioLignintoN
  INTEGER npt     

  casaflux%fromPtoL(:,:,:)      = 0.0
  casaflux%kplant(:,:)          = 0.0   ! (BPjun2010)

  WHERE(casamet%iveg2/=icewater) 
  ! using max function to avoid dividing by zero, ypw 14/may/2008
    ratioLignintoN(:,leaf) = (casapool%Cplant(:,leaf) &
                             /(max(1.0e-10,casapool%Nplant(:,leaf)) *casabiome%ftransNPtoL(veg%iveg(:),leaf))) &
                             * casabiome%fracLigninplant(veg%iveg(:),leaf)  
    ratioLignintoN(:,froot)= (casapool%Cplant(:,froot)&
                             /(max(1.0e-10,casapool%Nplant(:,froot))*casabiome%ftransNPtoL(veg%iveg(:),froot))) &
                             * casabiome%fracLigninplant(veg%iveg(:),froot) 

    casaflux%fromPtoL(:,metb,leaf)    = max(0.001, (0.75*(0.85 - 0.013 *ratioLignintoN(:,leaf))))  ! WW modified to reflect MIMICS parameterization
    casaflux%fromPtoL(:,metb,froot)   = max(0.001, (0.75*(0.85 - 0.013 *ratioLignintoN(:,froot)))) ! provides lower fMET estimates, Dec 6, 2017
    casaflux%fromPtoL(:,str,leaf)    = 1.0 - casaflux%fromPtoL(:,metb,leaf)
    casaflux%fromPtoL(:,str,froot)   = 1.0 - casaflux%fromPtoL(:,metb,froot) 
    casaflux%fromPtoL(:,cwd,wood)    = 1.0

    casaflux%kplant(:,leaf)        = casabiome%plantrate(veg%iveg(:),leaf)*xkleaf(:) &
                                   + xkleafcold(:) + xkleafdry(:)
    casaflux%kplant(:,wood)        = casabiome%plantrate(veg%iveg(:),wood) 
    casaflux%kplant(:,froot)       = casabiome%plantrate(veg%iveg(:),froot) 
  ENDWHERE

  ! When glai<glaimin,leaf biomass will not decrease anymore. (Q.Zhang 10/03/2011)
  DO npt = 1,mp 
    if(casamet%glai(npt).le.casabiome%glaimin(veg%iveg(npt))) casaflux%kplant(npt,leaf) = 0.0
  ENDDO
  ! end change

END SUBROUTINE casa_coeffplant

SUBROUTINE casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
!  calculate the plant litter fall rate, litter fall and sOM decomposition rate (1/day)
!  and the transfer coefficients between different pools
!  NOTE: this subroutine is not called when MIMICS or CORPSE are the SOM model.
!        See mimics_soil_forwardMM and mimics_soil_reverseMM.
!
! inputs:
!     xk(mp):          modifier of soil litter decomposition rate (dimensionless)
!
! outputs:
!     klitter(mp,mlitter):      decomposition rate of litter pool (1/day)
!     ksoil(mp,msoil):          decomposition rate of soil pool (1/day)
!     fromLtoS(mp,mlitter,msoil):  fraction of decomposed litter to soil (fraction)
!     fromStoS(mp,msoil,msoil):    fraction of decomposed soil C to another soil pool (fraction)
!     fromLtoCO2(mp,mlitter):      fraction of decomposed litter emitted as CO2 (fraction)
!     fromStoCO2(mp,msoil):        fraction of decomposed soil C emitted as CO2 (fraction)

  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(IN) :: xklitter,xksoil
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  INTEGER j,k,kk,nland,nv       !i: for plant pool, j for litter, k for soil

  casaflux%fromLtoS(:,:,:)      = 0.0   
  casaflux%fromStoS(:,:,:)      = 0.0                                
                                          ! flow from soil to soil
  DO k = 1, msoil
     casaflux%fromStoS(:,k,k)   = -1.0
  ENDDO  ! "k"
  casaflux%fromLtoCO2(:,:) = 0.0             ! flow from L or S to CO2
  casaflux%fromStoCO2(:,:) = 0.0

  casaflux%klitter(:,:) = 0.0        !initialize klitter (Q.Zhang 03/03/2011)

  WHERE(casamet%iveg2/=icewater) 

    casaflux%klitter(:,metb)   = xklitter(:) * casabiome%litterrate(veg%iveg(:),metb) 
    casaflux%klitter(:,str)    = xklitter(:) * casabiome%litterrate(veg%iveg(:),str) &
                                 * exp(-3.0*casabiome%fracLigninplant(veg%iveg(:),leaf))
    casaflux%klitter(:,cwd)    = xklitter(:) * casabiome%litterrate(veg%iveg(:),cwd)            

    casaflux%ksoil(:,mic)      = xksoil(:) * casabiome%soilrate(veg%iveg(:),mic)   &
                               * (1.0 - 0.75 *(soil%silt(:)+soil%clay(:)))
    casaflux%ksoil(:,slow)     = xksoil(:) * casabiome%soilrate(veg%iveg(:),slow)
    casaflux%ksoil(:,pass)     = xksoil(:) * casabiome%soilrate(veg%iveg(:),pass)
    casaflux%kplab(:)          = xksoil(:) * casabiome%xkplab(casamet%isorder(:))
    casaflux%kpsorb(:)         = xksoil(:) * casabiome%xkpsorb(casamet%isorder(:))
    casaflux%kpocc(:)          = xksoil(:) * casabiome%xkpocc(casamet%isorder(:))


    WHERE(veg%iveg==cropland)      ! for cultivated land type
       casaflux%ksoil(:,mic)  = casaflux%ksoil(:,mic) * 1.25
       casaflux%ksoil(:,slow) = casaflux%ksoil(:,slow)* 1.5
       casaflux%ksoil(:,pass) = casaflux%ksoil(:,pass)* 1.5 
    ENDWHERE  ! 

    !! Allow carbon use efficiencies to be read from the casa pft parameter file. -MDH 12/23/2019 
!                                         ! flow from litter to soil 
!   casaflux%fromLtoS(:,mic,metb)   = 0.45                                  
!                                         ! metb -> mic
!   casaflux%fromLtoS(:,mic,str)   = 0.45*(1.0-casabiome%fracLigninplant(veg%iveg(:),leaf))  
!                                         ! str -> mic
!   casaflux%fromLtoS(:,slow,str)  = 0.7 * casabiome%fracLigninplant(veg%iveg(:),leaf)       
!                                         ! str -> slow
!   casaflux%fromLtoS(:,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(:),wood)) 
!                                         ! CWD -> fmic
!   casaflux%fromLtoS(:,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(:),wood)        
!                                         ! CWD -> slow
!
!! set the following two backflow to set (see Bolker 199x)
!    casaflux%fromStoS(:,mic,slow)  = 0.45 * (0.997 - 0.009 *soil%clay(:))
!    casaflux%fromStoS(:,mic,pass)  = 0.45
!
!   casaflux%fromStoS(:,slow,mic)  = (0.85 - 0.68 * (soil%clay(:)+soil%silt(:))) &
!                                    * (0.997 - 0.032*soil%clay(:))
!   casaflux%fromStoS(:,pass,mic)  = (0.85 - 0.68 * (soil%clay(:)+soil%silt(:))) &
!                                    * (0.003 + 0.032*soil%clay(:))
!   casaflux%fromStoS(:,pass,slow) = 0.45 * (0.003 + 0.009 * soil%clay(:)) 


    casaflux%fromLtoS(:,mic,metb)  = casabiome%CUEmetbmic(veg%iveg(:))                                  
                                          ! metb -> mic
    casaflux%fromLtoS(:,mic,str)   = casabiome%CUEstrmic(veg%iveg(:)) * (1.0-casabiome%fracLigninplant(veg%iveg(:),leaf))  
                                          ! str -> mic
    casaflux%fromLtoS(:,slow,str)  = casabiome%CUEstrslow(veg%iveg(:)) * casabiome%fracLigninplant(veg%iveg(:),leaf)       
                                          ! str -> slow
    casaflux%fromLtoS(:,mic,cwd)   = casabiome%CUEcwdmic(veg%iveg(:)) * (1.0-casabiome%fracLigninplant(veg%iveg(:),wood)) 
                                          ! CWD -> fmic
    casaflux%fromLtoS(:,slow,cwd)  = casabiome%CUEcwdslow(veg%iveg(:)) * casabiome%fracLigninplant(veg%iveg(:),wood)        
                                          ! CWD -> slow
 
     !! set the following two backflow to set (see Bolker 199x)
     !casaflux%fromStoS(nv,mic,slow)  = 0.0
     !casaflux%fromStoS(nv,mic,pass)  = 0.0

    casaflux%fromStoS(:,slow,mic)  = casabiome%CUEmicslow(veg%iveg(:)) * (0.85 - 0.68*(soil%clay(:)+soil%silt(:))) &
                                     * (0.997 - 0.032*soil%clay(:))
    casaflux%fromStoS(:,pass,mic)  = casabiome%CUEmicpass(veg%iveg(:)) * (0.85 - 0.68*(soil%clay(:)+soil%silt(:))) &
                                     * (0.003 + 0.032*soil%clay(:))
    casaflux%fromStoS(:,pass,slow) = casabiome%CUEpassslow(veg%iveg(:)) * (0.003 + 0.009*soil%clay(:)) 

  ENDWHERE
   
  DO nland=1,mp
    IF(casamet%iveg2(nland)/=icewater) THEN
      DO j=1,mlitter
        DO k=1,msoil
          casaflux%fromLtoCO2(nland,j) = casaflux%fromLtoCO2(nland,j)  &
                                       + casaflux%fromLtoS(nland,k,j)
        ENDDO  !"k"
        casaflux%fromLtoCO2(nland,j) = 1.0 - casaflux%fromLtoCO2(nland,j) 
      ENDDO !"j"
      DO k=1,msoil
        DO kk=1,msoil
          casaflux%fromStoCO2(nland,k) = casaflux%fromStoCO2(nland,k) &
                                       + casaflux%fromStoS(nland,kk,k)
        ENDDO  !"kk"
      ENDDO   !"k"
      casaflux%fromStoCO2(nland,:) = -casaflux%fromStoCO2(nland,:)

      ! ATTENTION: should the eqn above be casaflux%fromStoCO2(nland,:) = 1.0-casaflux%fromStoCO2(nland,:) ??
      !write(*,*)
      !write(*,*) 'casaflux%fromStoS(nland,:,1)  =', casaflux%fromStoS(nland,:,1) 
      !write(*,*) 'casaflux%fromStoS(nland,:,2)  =', casaflux%fromStoS(nland,:,2) 
      !write(*,*) 'casaflux%fromStoS(nland,:,3)  =', casaflux%fromStoS(nland,:,3) 
      !write(*,*) 'casaflux%fromStoCO2(nland,1)  =', casaflux%fromStoCO2(nland,1) 
      !write(*,*) 'casaflux%fromStoCO2(nland,2)  =', casaflux%fromStoCO2(nland,2) 
      !write(*,*) 'casaflux%fromStoCO2(nland,3)  =', casaflux%fromStoCO2(nland,3) 
      ! The StoCO2 fractions appear to be OK because casaflux%fromStoS(:,k,k) = -1.0. -mdh 4/28/2017 

    ENDIF   
  ENDDO   ! "nland"

END SUBROUTINE casa_coeffsoil

! modified by ypw following Chris Lu 5/nov/2012
SUBROUTINE casa_delplant(veg,casabiome,casapool,casaflux,casamet,            &
                         cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                         nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                         pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

!  calculate the chnage in plant C, N and P pools
!  uptake of N and P will be computed in casa_uptake
!  labile C pool will be computed casa_labile
!  NOTE: this subroutine is not called when MIMICS or CORPSE are the SOM model.
!        See mimics_delplant (C only) and mimics_delplant_CN (C&N).

  IMPLICIT NONE
  TYPE (veg_parameter_type), INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),         INTENT(INOUT) :: casabiome
  TYPE (casa_pool),          INTENT(INOUT) :: casapool
  TYPE (casa_flux),          INTENT(INOUT) :: casaflux
  TYPE (casa_met),           INTENT(INOUT) :: casamet

  ! added by ypwang following Chris Lu 5/nov/2012
  real(r_2), dimension(mp),INTENT(OUT) :: cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                                     nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                                     pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd

  INTEGER  npt,nL,nP

   casaflux%FluxCtolitter = 0.0
   casaflux%FluxNtolitter = 0.0
   casaflux%FluxPtolitter = 0.0
   ! Added root exudate flux -mdh 1/13/2020
   casaflux%Cexudate = 0.0
   casaflux%Nexudate = 0.0
   casaflux%Pexudate = 0.0

   ! added by ypwang following Chris Lu 5/nov/2012

   cleaf2met = 0.0
   cleaf2str = 0.0
   croot2met = 0.0
   croot2str = 0.0
   cwood2cwd = 0.0

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

  !MPI
  DO npt=1,mp
  IF(casamet%iveg2(npt)/=icewater) THEN
!    PRINT *, 'npt = ', npt
!    PRINT *, 'casapool%cplant(npt,:) = ', casapool%cplant(npt,:)
    casapool%dCplantdt(npt,:)  =  casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,:)     &
                               - casaflux%kplant(npt,:)  * casapool%cplant(npt,:)

    ! Compute root exudate C flux as a fraction of froot NPP. -mdh 1/13/2020
    casaflux%Cexudate(npt) = max(0.0, casabiome%fracRootExud(veg%iveg(npt)) * casaflux%Cnpp(npt) * casaflux%fracCalloc(npt,froot))  
    casapool%dCplantdt(npt,froot) = casapool%dCplantdt(npt,froot) - casaflux%Cexudate(npt)   

    ! change here made by ypw on 26august 2011
    ! calculate fraction c to labile pool as a fraction of gpp, not npp
    ! casapool%dClabiledt(npt)   = casaflux%Cnpp(npt)    * casaflux%fracClabile(npt)
    casapool%dClabiledt(npt)   =  casaflux%Cgpp(npt)  * casaflux%fracClabile(npt) - casaflux%Clabloss(npt)
    ! added by ypwang 5/nov/2012
    cleaf2met(npt) = casaflux%fromPtoL(npt,metb,leaf)  * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
    cleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf)   * casaflux%kplant(npt,leaf)  * casapool%cplant(npt,leaf)
    ! Add casaflux%Cexudate(npt) to metabolic litter. -mdh 1/13/2019
    !croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
    croot2met(npt) = casaflux%fromPtoL(npt,metb,froot) * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot) &
                     + casaflux%Cexudate(npt)
    croot2str(npt) = casaflux%fromPtoL(npt,str,froot)  * casaflux%kplant(npt,froot) * casapool%cplant(npt,froot)
    cwood2cwd(npt) = casaflux%fromPtoL(npt,cwd,wood)   * casaflux%kplant(npt,wood)  * casapool%cplant(npt,wood)

    casaflux%ClitInptMet(npt) = cleaf2met(npt) + croot2met(npt)
    casaflux%ClitInptStruc(npt) = cleaf2str(npt) + croot2str(npt) 

!    PRINT *, 'npt, mp, iveg', npt, mp, veg%iveg(npt)
    IF(icycle > 1) THEN
!    PRINT *, 'casapool%Nplant(npt,:) = ', casapool%Nplant(npt,:)
       IF(casaflux%fracNalloc(npt,leaf)==0.0) THEN
          casapool%dNplantdt(npt,leaf)  = -casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf)
       else
          casapool%dNplantdt(npt,leaf)  = -casaflux%kplant(npt,leaf) * casapool%Nplant(npt,leaf) &
                                        * casabiome%ftransNPtoL(veg%iveg(npt),leaf)
       ENDIF
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

       nleaf2met(npt) = -casapool%dNplantdt(npt,leaf)  - nleaf2str(npt)
       ! dNplantdt includes casaflux%Nexudate. -mdh 1/13/2020
       nroot2met(npt) = -casapool%dNplantdt(npt,froot) - nroot2str(npt)
       nwood2cwd(npt) = -casapool%dNplantdt(npt,wood)

       casaflux%NlitInptMet(npt) = nleaf2met(npt) + nroot2met(npt)
       casaflux%NlitInptStruc(npt) = nleaf2str(npt) + nroot2str(npt) 

    ENDIF

!    PRINT *, 'before icycle >2; npt, mp', npt, mp
    IF(icycle >2) THEN
!    PRINT *, 'casapool%Pplant(npt,:) = ', casapool%Pplant(npt,:)
       IF(casaflux%fracPalloc(npt,leaf)==0.0) THEN
          casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf)
       else 
          casapool%dPplantdt(npt,leaf)  = - casaflux%kplant(npt,leaf) * casapool%Pplant(npt,leaf) &
                                        * casabiome%ftransPPtoL(veg%iveg(npt),leaf)
       ENDIF


       casapool%dPplantdt(npt,wood)  = - casaflux%kplant(npt,wood) * casapool%Pplant(npt,wood) &
                                     * casabiome%ftransPPtoL(veg%iveg(npt),wood)
       casapool%dPplantdt(npt,froot)  = - casaflux%kplant(npt,froot) * casapool%Pplant(npt,froot) &
                                     * casabiome%ftransPPtoL(veg%iveg(npt),froot)

       ! Compute root exudate P flux as a fraction of froot P uptake. -mdh 1/13/2020
       if (casaflux%Cexudate(npt) > 0.0) then
         casaflux%Pexudate(npt) = max(0.0,casabiome%fracRootExud(veg%iveg(npt)) * casaflux%Plabuptake(npt) &
                                 * casaflux%fracPalloc(npt,froot))
         casapool%dPplantdt(npt,froot) = casapool%dPplantdt(npt,froot) - casaflux%Pexudate(npt)
       endif
       ! added by ypwang 5/nov/2012

       pleaf2str(npt) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                      * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix
       proot2str(npt) = casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                      * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix
       pleaf2met(npt) = -casapool%dPplantdt(npt,leaf)  - pleaf2str(npt)
       ! dPplantdt includes casaflux%Pexudate. -mdh 1/13/2020
       proot2met(npt) = -casapool%dPplantdt(npt,froot) - proot2str(npt)
       pwood2cwd(npt) = -casapool%dPplantdt(npt,wood)

       
    ENDIF


    DO nL=1,mlitter
       DO nP=1,mplant
          casaflux%FluxCtolitter(npt,nL) = casaflux%FluxCtolitter(npt,nL) &
                                 + casaflux%fromPtoL(npt,nL,nP) &
                                 * casaflux%kplant(npt,nP) &
                                 * casapool%cplant(npt,nP)
       ENDDO
    ENDDO
 
    ! It appears that the calculation above did not account for Cexudate from froot to metb. -mdh 1/13/2020
    casaflux%FluxCtolitter(npt,metb) = casaflux%FluxCtolitter(npt,metb) + casaflux%Cexudate(npt)

!    PRINT *, 'before 2nd icycle >1; npt, mp', npt, mp
    IF(icycle > 1) THEN
       casaflux%FluxNtolitter(npt,str) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                               * casapool%cplant(npt,leaf)       * ratioNCstrfix              &     
                               + casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                               * casapool%cplant(npt,froot)      * ratioNCstrfix  
       ! FluxNtolitter includes Nexudate in dNplantdt calculation. -mdh 1/13/2019
       casaflux%FluxNtolitter(npt,metb) = -casapool%dNplantdt(npt,leaf)-casapool%dNplantdt(npt,froot) &
                                          - casaflux%FluxNtolitter(npt,str)
       casaflux%FluxNtolitter(npt,CWD) = -casapool%dNplantdt(npt,wood)

! adding N uptake
       casapool%dNplantdt(npt,:) = casapool%dNplantdt(npt,:) &
                                 + casaflux%Nminuptake(npt)*casaflux%fracNalloc(npt,:) 
! Note: This calculation was commented out because Nsoilmin is reduced by Nminuptake (== Nupland)
! in casa_delsoil instead.
!       casapool%Nsoilmin(npt)    = casapool%Nsoilmin(npt) - casaflux%Nminuptake(npt) *deltpool
    ENDIF !end "icycle >1"

!    PRINT *, 'second icycle > 2'
    IF(icycle>2) THEN
       casaflux%FluxPtolitter(npt,str) = casaflux%fromPtoL(npt,str,leaf) * casaflux%kplant(npt,leaf)  &
                               * casapool%cplant(npt,leaf)       * ratioNCstrfix/ratioNPstrfix        &     
                               + casaflux%fromPtoL(npt,str,froot)* casaflux%kplant(npt,froot) &
                               * casapool%cplant(npt,froot)      * ratioNCstrfix/ratioNPstrfix
       ! FluxPtolitter includes Pexudate in dPplantdt calculation. -mdh 1/13/2019
       casaflux%FluxPtolitter(npt,metb) = -casapool%dPplantdt(npt,leaf)-casapool%dPplantdt(npt,froot) &
                               - casaflux%FluxPtolitter(npt,str)
       casaflux%FluxPtolitter(npt,CWD) = -casapool%dPplantdt(npt,wood)
! add P uptake
       casapool%dPplantdt(npt,:) = casapool%dPplantdt(npt,:) &
                                 + casaflux%Plabuptake(npt)*casaflux%fracPalloc(npt,:) 
!       casapool%Psoillab(npt)    = casapool%Psoillab(npt) - casaflux%Plabuptake(npt) * deltpool
    ENDIF  !of "icycle >2"
!    PRINT *, 'End of all endifs'

  ENDIF
  ENDDO
!  PRINT *, 'Done casa_delplant; npt, mp', npt, mp

END SUBROUTINE casa_delplant

SUBROUTINE casa_delsoil(veg,casapool,casaflux,casamet,casabiome)
! calculate changes in litter and soil pools
! NOTE: this subroutine not called when MIMICS or CORPSE are the SOM model.
!       See mimics_soil_forwardMM and mimics_soil_reverseMM.

  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome

  ! local variables
  REAL(r_2), DIMENSION(mp)    :: xdplabsorb, fluxptase
  REAL(r_2), DIMENSION(mp)    :: cwd2str, nwd2str
  INTEGER j,jj,k,kk,kkk,nL,nS,nSS,nland

  casaflux%FluxCtoCO2    = 0.0
  casaflux%FluxCtosoil   = 0.0
  casaflux%FluxNtosoil   = 0.0
  casaflux%FluxPtosoil   = 0.0
  casaflux%Crsoil        = 0.0  ! initialization added by BP jul2010
  casaflux%CpassInpt     = 0.0  ! -mdh 12/3/2018

  casapool%dClitterdt = 0.0
  casapool%dCsoildt   = 0.0

  casapool%dNlitterdt = 0.0
  casapool%dNsoildt   = 0.0
  casapool%dNsoilmindt= 0.0
  casaflux%Nsmin = 0.0
  casaflux%Nsimm = 0.0
  casaflux%Nsnet = 0.0
  casaflux%Nminloss = 0.0
  casaflux%Nminleach = 0.0
  casaflux%Nlittermin=0.0

  casapool%dPlitterdt   = 0.0
  casapool%dPsoildt     = 0.0
  casapool%dPsoillabdt  = 0.0
  casapool%dPsoilsorbdt = 0.0
  casapool%dPsoiloccdt  = 0.0
  casaflux%Psmin = 0.0
  casaflux%Psimm = 0.0
  casaflux%Psnet = 0.0
  casaflux%Pleach = 0.0
  casaflux%Ploss  = 0.0
  casaflux%Plittermin=0.0
  fluxptase = 0.0
  cwd2str = 0.0      ! Added for output only (-mdh 1/23/2017)

DO nland=1,mp
IF(casamet%iveg2(nland)/=icewater) THEN
   DO nL=1,mlitter
      casaflux%FluxCtoCO2(nland) = casaflux%FluxCtoCO2(nland)  &
                        + casaflux%fromLtoCO2(nland,nL)  &
                        * casaflux%klitter(nland,nL) &
                        * casapool%clitter(nland,nL)
   ENDDO
      
   DO nS=1,msoil
      DO nL=1,mlitter
         casaflux%FluxCtosoil(nland,nS) = casaflux%FluxCtosoil(nland,nS) &
                               + casaflux%fromLtoS(nland,nS,nL) &
                               * casaflux%klitter(nland,nL) &
                               * casapool%clitter(nland,nL)
      ENDDO
      DO nSS=1,msoil
         IF(nSS/=nS) THEN
            casaflux%FluxCtosoil(nland,nS) = casaflux%FluxCtosoil(nland,nS) &
                                  + casaflux%fromStoS(nland,nS,nSS) &
                                  * casaflux%ksoil(nland,nSS) &
                                  * casapool%csoil(nland,nSS)
            ! nS is destination pool. Save inputs to the passive pool. -mdh 12/3/2018
            if (nS == PASS) then
               casaflux%CpassInpt(nland) = casaflux%FluxCtosoil(nland,nS)
            endif

         ENDIF
      ENDDO
      casaflux%FluxCtoCO2(nland) = casaflux%FluxCtoCO2(nland)  &
                        + casaflux%fromStoCO2(nland,nS) &
                        * casaflux%ksoil(nland,nS) &
                        * casapool%csoil(nland,nS)
      !write(*,*) 'soilresp(',nS,')=', casaflux%fromStoCO2(nland,nS) &
      !                 * casaflux%ksoil(nland,nS) &
      !                 * casapool%csoil(nland,nS)

   ENDDO

   !-----------------------------------------------------------------------------------------------
   ! The MIMICS and CORPSE models add decayed CWD to structural litter instead of to soil directly.
   ! For consistency, record that flux in the output variable casaflux%ClitInptStruc(npt). -mdh 1/23/2017

   DO nS=1,msoil
      cwd2str(nland) = cwd2str(nland) + casaflux%fromLtoS(nland,nS,cwd) &
                               * casaflux%klitter(nland,cwd) &
                               * casapool%clitter(nland,cwd)

      ! Added nwd2str. -mdh 7/15/2019
      IF(icycle>1) THEN
        nwd2str(nland) = nwd2str(nland) + casaflux%fromLtoS(nland,nS,cwd) &
                               * casaflux%klitter(nland,cwd) &
                               * casapool%Nlitter(nland,cwd)
      ENDIF
   ENDDO
   casaflux%ClitInptStruc(nland) = casaflux%ClitInptStruc(nland) + cwd2str(nland)
   ! Added nwd2str. -mdh 7/15/2019
   casaflux%NlitInptStruc(nland) = casaflux%NlitInptStruc(nland) + nwd2str(nland)
   !-----------------------------------------------------------------------------------------------

   IF(icycle>1) THEN
      DO j=1,mlitter
         casaflux%Nlittermin(nland) = casaflux%Nlittermin(nland) &
                                    + casaflux%klitter(nland,j) * casapool%Nlitter(nland,j)
      ENDDO
      DO k=1,msoil
         casaflux%Nsmin(nland)   = casaflux%Nsmin(nland)   &
                                    + casaflux%ksoil(nland,k)   * casapool%Nsoil(nland,k)
      ENDDO    !gross mineralisation

      DO kk=1,msoil
         DO jj=1,mlitter    ! immobilisation from litter to soil
            casaflux%Nsimm(nland) = casaflux%Nsimm(nland) &
                                     - casaflux%fromLtoS(nland,kk,jj) &
                                     * casaflux%klitter(nland,jj)     &
                                     * casapool%Clitter(nland,jj)     &
                                     * casapool%ratioNCsoilnew(nland,kk)
         ENDDO
         DO kkk=1,msoil      ! immobilisation from soil to soil
            IF(kkk.ne.kk) THEN
               casaflux%Nsimm(nland) = casaflux%Nsimm(nland) &
                                        - casaflux%fromStoS(nland,kk,kkk)  &
                                        * casaflux%ksoil(nland,kkk) &
                                        * casapool%Csoil(nland,kkk) &
                                        * casapool%ratioNCsoilnew(nland,kk)
            ENDIF
         ENDDO
      ENDDO  ! immobilization
 
      casaflux%Nsnet(nland)=casaflux%Nlittermin(nland) &
                                 +casaflux%Nsmin(nland)   &
                                 +casaflux%Nsimm(nland)
                                      ! net mineralization 
      IF(casapool%Nsoilmin(nland)>2.0.AND.casamet%tsoilavg(nland)>273.12) THEN
        casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
                                   * MAX(0.0,casaflux%Nsnet(nland))
        casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                                   * MAX(0.0,casapool%Nsoilmin(nland))
      ELSE
        casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
                                   * MAX(0.0,casaflux%Nsnet(nland)) &
                                   * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
        casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
                                   * MAX(0.0,casapool%Nsoilmin(nland)) &
                                   * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
      ENDIF

!     IF(casapool%Nsoilmin(nland)>2.0.AND.casamet%tsoilavg(nland)>273.12) THEN
!       casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
!                                  * MAX(0.0,casaflux%Nsnet(nland))
!       casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
!                                  * MAX(0.0,casapool%Nsoilmin(nland))
!     ELSE
!!      ! ATTENTION: these calculations doesn't look right.  Why the products? 
!!      ! I am updating the code. The Nminloss calculation also resides in subroutine
!!      ! mimics_soil_reverseMM_CN (mimics_cycle_CN.f90) -mdh 7/29/2019
!!      casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
!!                                 * MAX(0.0,casaflux%Nsnet(nland)) &
!!                                 * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
!!      casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
!!                                 * MAX(0.0,casapool%Nsoilmin(nland)) &
!!                                 * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
! 
!       casaflux%Nminloss(nland)   = casaflux%fNminloss(nland)  &
!                                  * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
!       casaflux%Nminleach(nland)  = casaflux%fNminleach(nland) &
!                                  * MAX(0.0,casapool%Nsoilmin(nland)/2.0)
!     ENDIF
      DO k=1,msoil
         DO j=1,mlitter
            casaflux%FluxNtosoil(nland,k) =  casaflux%FluxNtosoil(nland,k)  &
                                 + casaflux%fromLtoS(nland,k,j) &
                                 * casaflux%klitter(nland,j)    &
                                 * casapool%Clitter(nland,j)    &
                                 * casapool%ratioNCsoilnew(nland,k)
         ENDDO  ! end of "j"
         DO kk=1,msoil
            IF(kk.ne.k) THEN
               casaflux%FluxNtosoil(nland,k) = casaflux%FluxNtosoil(nland,k)  &
                                    + casaflux%fromStoS(nland,k,kk) &
                                    * casaflux%ksoil(nland,kk)      &
                                    * casapool%Csoil(nland,kk)      &
                                    * casapool%ratioNCsoilnew(nland,k)
            ENDIF
         ENDDO ! end of "kk"
      ENDDO    ! end of "k"

   ENDIF !end of icycle >1

   IF(icycle >2) THEN
      DO j=1,mlitter
         casaflux%Plittermin(nland) = casaflux%Plittermin(nland) &
                                    + casaflux%klitter(nland,j) * casapool%Plitter(nland,j)
      ENDDO
      DO k=1,msoil
         casaflux%Psmin(nland)   = casaflux%Psmin(nland)   &
                                    + casaflux%ksoil(nland,k)   * casapool%Psoil(nland,k)
      ENDDO    !gross mineralisation

      DO kk=1,msoil
         DO jj=1,mlitter    ! immobilisation from litter to soil
!            casaflux%Psimm(nland) = casaflux%Psimm(nland) &
!                                     - casaflux%fromLtoS(nland,kk,jj) &
!                                     * casaflux%klitter(nland,jj)     &
!                                     * casapool%Clitter(nland,jj)     &
!                                     * casapool%ratioPCsoil(nland,kk)
            casaflux%Psimm(nland) = casaflux%Psimm(nland) &
                                     - casaflux%fromLtoS(nland,kk,jj) &
                                     * casaflux%klitter(nland,jj)     &
                                     * casapool%Nlitter(nland,jj)     &
                                     /casapool%ratioNPsoil(nland,kk)
!                                     * casapool%ratioPCsoil(nland,kk)/casapool%ratioNCsoil(nland,kk)
         ENDDO
         DO kkk=1,msoil      ! immobilisation from soil to soil
            IF(kkk.ne.kk) THEN
!               casaflux%Psimm(nland) = casaflux%Psimm(nland) &
!                                        - casaflux%fromStoS(nland,kk,kkk)  &
!                                        * casaflux%ksoil(nland,kkk) &
!                                        * casapool%Csoil(nland,kkk) &
!                                        * casapool%ratioPCsoil(nland,kk)
               casaflux%Psimm(nland) = casaflux%Psimm(nland) &
                                        - casaflux%fromStoS(nland,kk,kkk)  &
                                        * casaflux%ksoil(nland,kkk) &
                                        * casapool%Nsoil(nland,kkk) &
                                        /casapool%ratioNPsoil(nland,kk)
!                                        * casapool%ratioPCsoil(nland,kk)/casapool%ratioNCsoil(nland,kk)
            ENDIF
         ENDDO
      ENDDO  ! immobilization
 
      casaflux%Psnet(nland)=casaflux%Plittermin(nland) &
                                 +casaflux%Psmin(nland)   &
                                 +casaflux%Psimm(nland)
                                      ! net mineralization 

!      casaflux%Pleach(nland)  =  (1.0e-4) &
!                                 * max(0.0,casapool%Psoillab(nland))

      casaflux%Pleach(nland)  =  casaflux%fPleach(nland) &
                                 * max(0.0,casapool%Psoillab(nland))

      DO k=1,msoil
         DO j=1,mlitter
!            casaflux%FluxPtosoil(nland,k) =  casaflux%FluxPtosoil(nland,k)  &
!                                 + casaflux%fromLtoS(nland,k,j) &
!                                 * casaflux%klitter(nland,j)    &
!                                 * casapool%Clitter(nland,j)    &
!                                 * casapool%ratioPCsoil(nland,k)
            casaflux%FluxPtosoil(nland,k) =  casaflux%FluxPtosoil(nland,k)  &
                                 + casaflux%fromLtoS(nland,k,j) &
                                 * casaflux%klitter(nland,j)    &
                                 * casapool%Nlitter(nland,j)    &
                                 /casapool%ratioNPsoil(nland,k)
!                                 * casapool%ratioPCsoil(nland,k)/casapool%ratioNCsoil(nland,k)
         ENDDO  ! end of "j"
         DO kk=1,msoil
            IF(kk.ne.k) THEN
!               casaflux%FluxPtosoil(nland,k) = casaflux%FluxPtosoil(nland,k)  &
!                                    + casaflux%fromStoS(nland,k,kk) &
!                                    * casaflux%ksoil(nland,kk)      &
!                                    * casapool%Csoil(nland,kk)      &
!                                    * casapool%ratioPCsoil(nland,k)
               casaflux%FluxPtosoil(nland,k) = casaflux%FluxPtosoil(nland,k)  &
                                    + casaflux%fromStoS(nland,k,kk) &
                                    * casaflux%ksoil(nland,kk)      &
                                    * casapool%Nsoil(nland,kk)      &
                                    /casapool%ratioNPsoil(nland,k)
!                                    * casapool%ratioPCsoil(nland,k)/casapool%ratioNCsoil(nland,k)
            ENDIF
         ENDDO ! end of "kk"
      ENDDO    ! end of "k"
! need to account for flow from sorbed to occluded pool
   ENDIF
ENDIF  ! end of /=icewater
ENDDO  ! end of nland



!   write(*,991) casaflux%psnet(1883),casaflux%Plittermin(1883), casaflux%Psmin(1883),casaflux%Psimm(1883), &
!                casaflux%klitter(1883,:), casaflux%ksoil(1883,:),casaflux%fromLtoS(1883,:,:), casaflux%fromStoS(1883,:,:), &
!                casapool%clitter(1883,:),casapool%csoil(1883,:),casapool%nlitter(1883,:), casapool%nsoil(1883,:), &
!                casapool%plitter(1883,:),casapool%psoil(1883,:),  &
!                casapool%ratioPCsoil(1883,:),casapool%ratioNCsoil(1883,:)/casapool%ratioPCsoil(1883,:)
!991  format(' delsoil point 1883 ', 100(f15.8,1x))

DO nland=1,mp
IF(casamet%iveg2(nland)/=icewater) THEN
   casapool%dClitterdt(nland,:) =  casaflux%FluxCtolitter(nland,:) - casaflux%klitter(nland,:) * casapool%clitter(nland,:)
   casapool%dCsoildt(nland,:)   =  casaflux%FluxCtosoil(nland,:)   - casaflux%ksoil(nland,:)   * casapool%csoil(nland,:)
   casaflux%Crsoil(nland)       =  casaflux%FluxCtoCO2(nland)
   IF(icycle > 1) THEN
      casapool%dNlitterdt(nland,:) =  casaflux%FluxNtolitter(nland,:)  &
                                   - casaflux%klitter(nland,:) &
                                   * max(0.0,casapool%Nlitter(nland,:))

      casapool%dNsoildt(nland,:) = casaflux%FluxNtosoil(nland,:) &
                                 - casaflux%ksoil(nland,:) * casapool%Nsoil(nland,:)
      casapool%dNsoilmindt(nland)= casaflux%Nsnet(nland)&
                                 + casaflux%Nmindep(nland) + casaflux%Nminfix(nland)   &
                                 - casaflux%Nminloss(nland)    &
                                 - casaflux%Nminleach(nland)   &
                                 - casaflux%Nupland(nland)  
                            
   ENDIF

   IF(icycle >2) THEN

!      fluxptase(nland) =  casabiome%prodptase(veg%iveg(nland))  &
!                       *  max(0.0,(casapool%Psoil(nland,2)+casapool%Psoil(nland,3))) &
!                       *  max(0.0,(casabiome%costNpup(veg%iveg(nland))-15.0))/(max(0.0,(casabiome%costNpup(veg%iveg(nland))-15.0)) + 150.0)

      fluxptase(nland) =  casabiome%prodptase(veg%iveg(nland))  &
                       *  max(0.0,(casapool%Psoil(nland,2)*casaflux%ksoil(nland,2) &
                          +casapool%Psoil(nland,3)*casaflux%ksoil(nland,3))) &
                       *  max(0.0,(casabiome%costNpup(veg%iveg(nland))-15.0)) &
                          /(max(0.0,(casabiome%costNpup(veg%iveg(nland))-15.0)) + 150.0)
      !fluxptase(nland)  = 0.0
      xdplabsorb(nland) = 1.0+ casaflux%Psorbmax(nland)*casaflux%kmlabp(nland) &
                        /((casaflux%kmlabp(nland)+casapool%Psoillab(nland))**2)
      casapool%dPlitterdt(nland,:) = casaflux%fluxPtolitter(nland,:)  &
                                   - casaflux%klitter(nland,:)                 &
                                   * max(0.0,casapool%Plitter(nland,:))

      casapool%dPsoildt(nland,1) = casaflux%FluxPtosoil(nland,1)                        &
                                 - casaflux%ksoil(nland,1) * casapool%Psoil(nland,1)

      casapool%dPsoildt(nland,2) = casaflux%FluxPtosoil(nland,2)                        &
                                 - casaflux%ksoil(nland,2) * casapool%Psoil(nland,2)    &
                                 - fluxptase(nland) * casaflux%ksoil(nland,2)*casapool%Psoil(nland,2) &
                                  /(casaflux%ksoil(nland,2)*casapool%Psoil(nland,2)+casaflux%ksoil(nland,3)*casapool%Psoil(nland,3))

      casapool%dPsoildt(nland,3) = casaflux%FluxPtosoil(nland,3)                        &
                                 - casaflux%ksoil(nland,3) * casapool%Psoil(nland,3)    &
                                 - fluxptase(nland) * casaflux%ksoil(nland,3)*casapool%Psoil(nland,3) &
                                  /(casaflux%ksoil(nland,2)*casapool%Psoil(nland,2)+casaflux%ksoil(nland,3)*casapool%Psoil(nland,3))

      casapool%dPsoillabdt(nland)= casaflux%Psnet(nland) + fluxptase(nland)         &
                                 + casaflux%Pdep(nland) + casaflux%Pwea(nland)      &
                                 - casaflux%Pleach(nland)-casaflux%pupland(nland)   &
                                 - casaflux%kpsorb(nland)*casapool%Psoilsorb(nland) &
                                 + casaflux%kpocc(nland) * casapool%Psoilocc(nland)
      ! here the dPsoillabdt =(dPsoillabdt+dPsoilsorbdt)
      ! dPsoilsorbdt  = xdplabsorb
      casapool%dPsoillabdt(nland)  = casapool%dPsoillabdt(nland)/xdplabsorb(nland)
      casapool%dPsoilsorbdt(nland) = 0.0

      casapool%dPsoiloccdt(nland)  = casaflux%kpsorb(nland)* casapool%Psoilsorb(nland) &
                                   - casaflux%kpocc(nland) * casapool%Psoilocc(nland)
      ! P loss to non-available P pools
!      casaflux%Ploss(nland)        = casaflux%kpocc(nland) * casapool%Psoilocc(nland)

!      casaflux%Ploss(nland)       = casaflux%fPleach(nland) &
!                                 * max(0.0,casapool%Psoillab(nland))
      casaflux%Ploss(nland)       = 0.0                                   
   ENDIF
ENDIF
ENDDO

!  nland=4104
!  write(76,701) nland, casapool%Nsoilmin(nland),casaflux%Nlittermin(nland),casaflux%nsmin(nland),casaflux%nsimm(nland),casaflux%Nmindep(nland), &
!                casaflux%Nminfix(nland),casaflux%Nminloss(nland),casaflux%Nminleach(nland),casaflux%Nupland(nland),casaflux%Cnpp(nland), &
!                casapool%ratioNCsoilmin(nland,:),casapool%ratioNCsoilmax(nland,:),casapool%ratioNCsoilnew(nland,:)
!701 format(i6,100(f10.5,2x))
END SUBROUTINE casa_delsoil

SUBROUTINE avgsoil(veg,soil,casamet)
! Get avg soil moisture, avg soil temperature 
! need to estimate the land cell mean soil temperature and moisture weighted by the area fraction
! of each tile within the land cell

  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  INTEGER                     :: ns,nland

  casamet%tsoilavg      = 0.0
  casamet%moistavg      = 0.0
  casamet%frznmoistavg  = 0.0
  casamet%btran         = 0.0
!  print *, 'avgsoil: froot', veg%froot

  DO ns = 1, ms
  DO nland=1,mp
    casamet%tsoilavg(nland)  = casamet%tsoilavg(nland)+veg%froot(nland,ns)  &
                               * casamet%tsoil(nland,ns)

    ! Just noticed that this calculation prevents soil moisture > field capacity. -mdh 3/13/2017
    casamet%moistavg(nland)  = casamet%moistavg(nland)+ veg%froot(nland,ns) &
                               * min(soil%sfc(nland),casamet%moist(nland,ns)) 

    ! frznmoistavg calculation added. -mdh 3/13/2017
    casamet%frznmoistavg(nland) = casamet%frznmoistavg(nland)+ veg%froot(nland,ns) &
                                  * casamet%frznmoist(nland,ns) 

    casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
            * (min(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland)) &
            /(soil%sfc(nland)-soil%swilt(nland))

    !print *, 'veg%froot: ', nland, ns, veg%froot(nland,ns)
    !print *, 'tsoilavg: ', nland, casamet%tsoilavg(nland)
  ENDDO
  ENDDO

END SUBROUTINE avgsoil

SUBROUTINE casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
! (1) compute (1)N uptake by plants; 
! (2) allocation of uptaken N to plants 
!
  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  REAL(r_2), DIMENSION(mp),     INTENT(IN)    :: xkNlimiting
  !INTEGER,                      INTENT(IN)    :: mp

  ! local variables
  INTEGER                              :: np
  REAL(r_2), DIMENSION(mp,mplant)      :: Nreqmax, Nreqmin, NtransPtoP, xnuptake
  REAL(r_2), DIMENSION(mp)             :: totNreqmax,totNreqmin
  REAL(r_2), DIMENSION(mp)             :: xnCnpp

  Nreqmin(:,:)       = 0.0
  Nreqmax(:,:)       = 0.0
  NtransPtoP(:,:)    = 0.0
  totNreqmax = 0.0
  totNreqmin = 0.0

  casaflux%Nminuptake(:)     = 0.0
  casaflux%fracNalloc(:,:)   = 0.0
  xnCnpp = max(0.0,casaflux%Cnpp)
  call casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
  
  DO np=1,mp
  IF(casamet%iveg2(np)/=icewater) THEN 
    totNreqmax(np) = Nreqmax(np,leaf)+Nreqmax(np,wood)+Nreqmax(np,froot)
    totNreqmin(np) = Nreqmin(np,leaf)+Nreqmin(np,wood)+Nreqmin(np,froot)

    xnuptake(np,leaf) = Nreqmin(np,leaf) + xkNlimiting(np)* (Nreqmax(np,leaf)-Nreqmin(np,leaf))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))
    xnuptake(np,wood) = Nreqmin(np,wood) + xkNlimiting(np)* (Nreqmax(np,wood)-Nreqmin(np,wood))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))
    xnuptake(np,froot) = Nreqmin(np,froot) + xkNlimiting(np)* (Nreqmax(np,froot)-Nreqmin(np,froot))     &
                      *  casapool%Nsoilmin(np)/(casapool%Nsoilmin(np)+casabiome%kminN(veg%iveg(np)))

    casaflux%Nminuptake(np) = xnuptake(np,leaf) + xnuptake(np,wood) + xnuptake(np,froot)+1.0e-10
    casaflux%fracNalloc(np,leaf)  = xnuptake(np,leaf)/casaflux%Nminuptake(np)
    casaflux%fracNalloc(np,wood)  = xnuptake(np,wood)/casaflux%Nminuptake(np)
    casaflux%fracNalloc(np,froot) = xnuptake(np,froot)/casaflux%Nminuptake(np)
  ENDIF
  ENDDO

!  np=1
!  write(*,911) casapool%nsoilmin(np),casaflux%Nminuptake(np),xnuptake(np,leaf), xnuptake(np,wood), xnuptake(np,froot), &
!               casaflux%fracNalloc(np,leaf),casaflux%fracNalloc(np,wood),casaflux%fracNalloc(np,froot)
!911 format('N uptake:',100(f8.3,2x))

  casaflux%Nupland = casaflux%Nminuptake

END SUBROUTINE casa_nuptake

SUBROUTINE casa_Nrequire(xnCnpp,Nreqmin,Nreqmax,NtransPtoP,veg, &
                         casabiome,casapool,casaflux,casamet)
!
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp),        INTENT(IN)    :: xnCnpp
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: Nreqmax, Nreqmin
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: NtransPtoP
  TYPE (veg_parameter_type),             INTENT(INOUT) :: veg
  TYPE (casa_biome),                     INTENT(INOUT) :: casabiome
  TYPE (casa_pool),                      INTENT(INOUT) :: casapool
  TYPE (casa_flux),                      INTENT(INOUT) :: casaflux
  TYPE (casa_met),                       INTENT(INOUT) :: casamet

  ! local variable
  INTEGER :: np
  REAL(r_2), DIMENSION(mp,mplant)     :: ncplantmax
  REAL(r_2) :: fNC

  Nreqmin(:,:)    = 0.0
  Nreqmax(:,:)    = 0.0
  NtransPtoP(:,:) = 0.0
  
  ! It appears that the new parameters casabiome%xkNlimitmin and casabiome%xkNlimitmax should
  ! replace 0.5 and 2.0 below. The multiplier is 0.0 <= max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0) <= 1.0,
  ! and therefore casabiome%ratioNCplantmin <= ncplantmax <= casabiome%ratioNCplantmax,
  ! but changing values 0.5 and 2.0 alters the range of this multiplier.
  ! This 0.0-1.0 multiplier is a near linear function so I replaced it with a linear
  ! function that uses casabiome%xkNlimitmin and casabiome%xkNlimitmax.  -mdh 12/31/2019.

  DO np=1,mp
  IF(casamet%iveg2(np)/=icewater) THEN
    !if(casapool%Nsoilmin(np)<2.0) then
    if (casapool%Nsoilmin(np)< casabiome%xkNlimitmax(veg%iveg(np))) then
     !  ncplantmax(np,leaf) =casabiome%ratioNCplantmin(veg%iveg(np),leaf)  &
     !                      +(casabiome%ratioNCplantmax(veg%iveg(np),leaf)-casabiome%ratioNCplantmin(veg%iveg(np),leaf)) &
     !                        * min(1.0,max(0.0,casapool%Nsoilmin(np)*0.5)) 
     !  ncplantmax(np,wood) =casabiome%ratioNCplantmin(veg%iveg(np),wood)  &
     !                      +(casabiome%ratioNCplantmax(veg%iveg(np),wood)-casabiome%ratioNCplantmin(veg%iveg(np),wood)) &
     !                        * min(1.0,max(0.0,casapool%Nsoilmin(np)*0.5)) 
     !  ncplantmax(np,froot) =casabiome%ratioNCplantmin(veg%iveg(np),froot)  &
     !                      +(casabiome%ratioNCplantmax(veg%iveg(np),froot)-casabiome%ratioNCplantmin(veg%iveg(np),froot)) &
     !                        * min(1.0,max(0.0,casapool%Nsoilmin(np)*0.5)) 

     ! Code replaced. -mdh 12/31/2019
     ! ncplantmax(np,leaf) =casabiome%ratioNCplantmin(veg%iveg(np),leaf)  &
     !                     +(casabiome%ratioNCplantmax(veg%iveg(np),leaf)-casabiome%ratioNCplantmin(veg%iveg(np),leaf)) &
     !                       * min(1.0,max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0)) 
     ! ncplantmax(np,wood) =casabiome%ratioNCplantmin(veg%iveg(np),wood)  &
     !                     +(casabiome%ratioNCplantmax(veg%iveg(np),wood)-casabiome%ratioNCplantmin(veg%iveg(np),wood)) &
     !                       * min(1.0,max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0)) 
     ! ncplantmax(np,froot) =casabiome%ratioNCplantmin(veg%iveg(np),froot)  &
     !                     +(casabiome%ratioNCplantmax(veg%iveg(np),froot)-casabiome%ratioNCplantmin(veg%iveg(np),froot)) &
     !                       * min(1.0,max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0)) 

       ! N limitation, so use plant N:C ratio that is between the min N:C and max N:C to compute N demand. 
       ! Linear ramp from 0.0 to 1.0 for casabiome%xkNlimitmin < casapool%Nsoilmin(:) < casabiome%xkNlimitmax
       ! y = (y2 - y1) / (x2 - x1) * (x - x2) + y2
       ! max(0.0,fNC) = 0 when casapool%Nsoilmin <= casabiome%xkNlimitmin, low N 
       !                  ncplantmax(np,:) = casabiome%ratioNCplantmin(veg%iveg(np),:)
       ! max(0.0,fNC) = 1 when casapool%Nsoilmin >= casabiome%xkNlimitmax, ample N 
       !                  ncplantmax(np,:) = casabiome%ratioNCplantmin(veg%iveg(np),:)) + NC ratio difference

       fNC = (1.0 - 0.0) / (casabiome%xkNlimitmax(veg%iveg(np)) - casabiome%xkNlimitmin(veg%iveg(np))) &
                       * (casapool%Nsoilmin(np) - casabiome%xkNlimitmax(veg%iveg(np))) + 1.0

       ncplantmax(np,leaf) =casabiome%ratioNCplantmin(veg%iveg(np),leaf)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),leaf)-casabiome%ratioNCplantmin(veg%iveg(np),leaf)) &
                             * min(1.0,max(0.0,fNC)) 
       ncplantmax(np,wood) =casabiome%ratioNCplantmin(veg%iveg(np),wood)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),wood)-casabiome%ratioNCplantmin(veg%iveg(np),wood)) &
                             * min(1.0,max(0.0,fNC)) 
       ncplantmax(np,froot) =casabiome%ratioNCplantmin(veg%iveg(np),froot)  &
                           +(casabiome%ratioNCplantmax(veg%iveg(np),froot)-casabiome%ratioNCplantmin(veg%iveg(np),froot)) &
                             * min(1.0,max(0.0,fNC)) 
    else
      ! No N limitation, so use maximum plant N:C ratio to compute N demand. 
      ncplantmax(np,leaf)  = casabiome%ratioNCplantmax(veg%iveg(np),leaf) 
      ncplantmax(np,wood)  = casabiome%ratioNCplantmax(veg%iveg(np),wood) 
      ncplantmax(np,froot) = casabiome%ratioNCplantmax(veg%iveg(np),froot) 
    endif

    Nreqmax(np,leaf)  = xnCnpp(np)* casaflux%fracCalloc(np,leaf) *ncplantmax(np,leaf)
    Nreqmax(np,wood)  = xnCnpp(np)* casaflux%fracCalloc(np,wood) *ncplantmax(np,wood)
    Nreqmax(np,froot) = xnCnpp(np)* casaflux%fracCalloc(np,froot)*ncplantmax(np,froot)

    Nreqmin(np,leaf) =  xnCnpp(np)* casaflux%fracCalloc(np,leaf) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),leaf)
    Nreqmin(np,wood) =  xnCnpp(np)* casaflux%fracCalloc(np,wood) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),wood)
    Nreqmin(np,froot) =  xnCnpp(np)* casaflux%fracCalloc(np,froot) &
                       * casabiome%ratioNCplantmin(veg%iveg(np),froot)

    NtransPtoP(np,leaf) = casaflux%kplant(np,leaf)*casapool%Nplant(np,leaf) &
                       * (1.0-casabiome%ftransNPtoL(veg%iveg(np),leaf))
    NtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Nplant(np,wood) &
                       * (1.0-casabiome%ftransNPtoL(veg%iveg(np),wood))
    NtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Nplant(np,froot) &
                       * (1.0-casabiome%ftransNPtoL(veg%iveg(np),froot))

    Nreqmax(np,leaf)  = max(0.0,Nreqmax(np,leaf) - NtransPtoP(np,leaf))
    Nreqmax(np,wood)  = max(0.0,Nreqmax(np,wood) - NtransPtoP(np,wood))
    Nreqmax(np,froot) = max(0.0,Nreqmax(np,froot) - NtransPtoP(np,froot))
    Nreqmin(np,leaf)  = max(0.0,Nreqmin(np,leaf) - NtransPtoP(np,leaf))
    Nreqmin(np,wood)  = max(0.0,Nreqmin(np,wood) - NtransPtoP(np,wood))
    Nreqmin(np,froot) = max(0.0,Nreqmin(np,froot) - NtransPtoP(np,froot))

    if(casapool%nplant(np,leaf)/(casapool%cplant(np,leaf)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),leaf)) then 
       Nreqmax(np,leaf) = 0.0 
       Nreqmin(np,leaf) =0.0
    endif
    if(casapool%nplant(np,wood)/(casapool%cplant(np,wood)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),wood)) then
       Nreqmax(np,wood) = 0.0
       Nreqmin(np,wood) =0.0
    endif
    if(casapool%nplant(np,froot)/(casapool%cplant(np,froot)+1.0e-10)>casabiome%ratioNCplantmax(veg%iveg(np),froot)) then
       Nreqmax(np,froot) = 0.0
       Nreqmin(np,froot) =0.0
    endif

  ENDIF
  ENDDO

!  np=4104  
!  write(77,701) np,1.0/casabiome%rationcplantmin(veg%iveg(np),:),1.0/casabiome%rationcplantmax(veg%iveg(np),:),casapool%Nsoilmin(np), &
!                   1.0/ncplantmax(np,:),Nreqmin(np,:),Nreqmax(np,:),casapool%cplant(np,:)/casapool%nplant(np,:),  &
!                       casapool%cplant(np,:),casapool%nplant(np,:), xnCnpp(np),casaflux%fracCalloc(np,:), min(1.0,max(0.0,2.0**(0.5*casapool%Nsoilmin(np))-1.0))
!701 format(i6,100(f12.5,2x))
END SUBROUTINE casa_Nrequire

SUBROUTINE casa_puptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
! (1) compute  P uptake by plants; 
! (2) allocation of uptaken P to plants 
!
  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  REAL(r_2), DIMENSION(mp),     INTENT(IN)    :: xkNlimiting


  ! local variables
  !INTEGER                        :: np
  REAL(r_2), DIMENSION(mp,mplant) :: Preqmax,Preqmin,PtransPtoP,xPuptake
  REAL(r_2), DIMENSION(mp)        :: totPreqmax,totPreqmin
  REAL(r_2), DIMENSION(mp)        :: xpCnpp
   
  Preqmin(:,:)             = 0.0
  Preqmax(:,:)             = 0.0
  PtransPtoP(:,:)          = 0.0
  casaflux%Plabuptake(:)   = 0.0
  casaflux%fracPalloc(:,:) = 0.0
  totPreqmax               = 0.0
  totPreqmin               = 0.0

  xpCnpp = max(0.0,casaflux%Cnpp)
  call casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                     casabiome,casapool,casaflux,casamet)
  WHERE(casamet%iveg2/=icewater) 
    totPreqmax(:) = Preqmax(:,leaf)+Preqmax(:,wood)+Preqmax(:,froot)
    totPreqmin(:) = Preqmin(:,leaf)+Preqmin(:,wood)+Preqmin(:,froot)

    xpuptake(:,leaf) = Preqmin(:,leaf) + xkNlimiting(:)* (Preqmax(:,leaf)-Preqmin(:,leaf))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))
    xpuptake(:,wood) = Preqmin(:,wood) + xkNlimiting(:)* (Preqmax(:,wood)-Preqmin(:,wood))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))
    xpuptake(:,froot) = Preqmin(:,froot) + xkNlimiting(:)* (Preqmax(:,froot)-Preqmin(:,froot))     &
                      *  casapool%Psoillab(:)/(casapool%Psoillab(:)+casabiome%KuplabP(veg%iveg(:)))

    casaflux%Plabuptake(:) = xpuptake(:,leaf) + xpuptake(:,wood) + xpuptake(:,froot)+1.0e-10
    casaflux%fracPalloc(:,leaf)  = xpuptake(:,leaf)/casaflux%Plabuptake(:)
    casaflux%fracPalloc(:,wood)  = xpuptake(:,wood)/casaflux%Plabuptake(:)
    casaflux%fracPalloc(:,froot) = xpuptake(:,froot)/casaflux%Plabuptake(:)

  ENDWHERE

  casaflux%Pupland = casaflux%Plabuptake 

!  np=1921
!  write(77,911) casapool%Psoillab(np),casaflux%Plabuptake(np), &
!               xpuptake(np,leaf), xpuptake(np,wood), xpuptake(np,froot), &
!               casaflux%fracPalloc(np,leaf),casaflux%fracPalloc(np,wood), &
!               casaflux%fracPalloc(np,froot), casabiome%ratioNPplantmin(veg%iveg(np),wood),&
!               xkNlimiting(np),Preqmax(np,wood),Preqmin(np,wood), &
!               (casapool%Psoillab(np),casabiome%KuplabP(veg%iveg(np)))
!911 format('P uptake:',100(f8.3,2x))

  ! only used in spinning up the model
!  DO  np=1,mp
!    casaflux%Plabuptake(np) = TotPreqmax(np) 
!    casaflux%Pupland(np)    = TotPreqmax(np)
!    casaflux%Pwea(np)       = TotPreqmax(np)
!  ENDDO

END SUBROUTINE casa_puptake

SUBROUTINE casa_Prequire(xpCnpp,Preqmin,Preqmax,PtransPtoP,veg, &
                         casabiome,casapool,casaflux,casamet)
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp),        INTENT(IN)    :: xpCnpp
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: Preqmax, Preqmin
  REAL(r_2), DIMENSION(mp,mplant), INTENT(INOUT) :: PtransPtoP
  TYPE (veg_parameter_type),             INTENT(INOUT) :: veg
  TYPE (casa_biome),                     INTENT(INOUT) :: casabiome
  TYPE (casa_pool),                      INTENT(INOUT) :: casapool
  TYPE (casa_flux),                      INTENT(INOUT) :: casaflux
  TYPE (casa_met),                       INTENT(INOUT) :: casamet

  ! local variables
  !INTEGER :: nland,ip
  INTEGER :: np,mp

  Preqmin(:,:)       = 0.0
  Preqmax(:,:)       = 0.0
  PtransPtoP(:,:)    = 0.0
  do np=1,mp
  IF(casamet%iveg2(np)/=icewater) then
    Preqmax(np,leaf) = xpCnpp(np)* casaflux%fracCalloc(np,leaf) &
                    * (casapool%Nplant(np,leaf)/(casapool%Cplant(np,leaf)+1.0e-10))/casabiome%ratioNPplantmin(veg%iveg(np),leaf)
    Preqmax(np,wood) = xpCnpp(np)* casaflux%fracCalloc(np,wood) &
                    * (casapool%Nplant(np,wood)/(casapool%Cplant(np,wood)+1.0e-10))/casabiome%ratioNPplantmin(veg%iveg(np),wood)
    Preqmax(np,froot) = xpCnpp(np)* casaflux%fracCalloc(np,froot) &
                    * (casapool%Nplant(np,froot)/(casapool%Cplant(np,froot)+1.0e-10))/casabiome%ratioNPplantmin(veg%iveg(np),froot)

    Preqmin(np,leaf) = xpCnpp(np) * casaflux%fracCalloc(np,leaf) &
                    * (casapool%Nplant(np,leaf)/(casapool%Cplant(np,leaf)+1.0e-10))/casabiome%ratioNPplantmax(veg%iveg(np),leaf)
    Preqmin(np,wood) = xpCnpp(np) * casaflux%fracCalloc(np,wood) &
                    * (casapool%Nplant(np,wood)/(casapool%Cplant(np,wood)+1.0e-10))/casabiome%ratioNPplantmax(veg%iveg(np),wood)
    Preqmin(np,froot) = xpCnpp(np) * casaflux%fracCalloc(np,froot) &
                    * (casapool%Nplant(np,froot)/(casapool%Cplant(np,froot)+1.0e-10))/casabiome%ratioNPplantmax(veg%iveg(np),froot)

    PtransPtoP(np,leaf) = casaflux%kplant(np,leaf)*casapool%Pplant(np,leaf) &
                       * (1.0-casabiome%ftransPPtoL(veg%iveg(np),leaf))
    PtransPtoP(np,wood) = casaflux%kplant(np,wood)*casapool%Pplant(np,wood) &
                       * (1.0-casabiome%ftransPPtoL(veg%iveg(np),wood))
    PtransPtoP(np,froot) = casaflux%kplant(np,froot)*casapool%Pplant(np,froot) &
                       * (1.0-casabiome%ftransPPtoL(veg%iveg(np),froot))

    Preqmax(np,leaf)    = max(0.0,Preqmax(np,leaf) - PtransPtoP(np,leaf))
    Preqmax(np,wood)    = max(0.0,Preqmax(np,wood) - PtransPtoP(np,wood))
    Preqmax(np,froot)    = max(0.0,Preqmax(np,froot) - PtransPtoP(np,froot))

    Preqmin(np,leaf)    = max(0.0,Preqmin(np,leaf) - PtransPtoP(np,leaf))
    Preqmin(np,wood)    = max(0.0,Preqmin(np,wood) - PtransPtoP(np,wood))
    Preqmin(np,froot)    = max(0.0,Preqmin(np,froot) - PtransPtoP(np,froot))


    if(casapool%pplant(np,leaf)/(casapool%nplant(np,leaf)+1.0e-10)> 1.0/casabiome%ratioNPplantmin(veg%iveg(np),leaf)) then
       Preqmax(np,leaf) = 0.0
       Preqmin(np,leaf) =0.0
    endif
    if(casapool%pplant(np,wood)/(casapool%nplant(np,wood)+1.0e-10)> 1.0/casabiome%ratioNPplantmin(veg%iveg(np),wood)) then
       Preqmax(np,wood) = 0.0
       Preqmin(np,wood) =0.0
    endif
    if(casapool%pplant(np,froot)/(casapool%nplant(np,froot)+1.0e-10)> 1.0/casabiome%ratioNPplantmin(veg%iveg(np),froot)) then
       Preqmax(np,froot) = 0.0
       Preqmin(np,froot) =0.0
    endif

  endif
  ENDDO    

! yow
!    ip=1921
!    write(77,791)  ip,veg%iveg(ip),xpCnpp(ip) , casaflux%fracCalloc(ip,wood),  &
!                   casapool%Nplant(ip,wood), casapool%Cplant(ip,wood), casabiome%ratioNPplantmax(veg%iveg(ip),wood)
!
!791 format('P require =', 2(i6,2x),100(f10.4,2x))
END SUBROUTINE casa_Prequire


SUBROUTINE casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)
! update all pool sizes
!
  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  ! local variables
  !REAL(r_2), DIMENSION(mp)   :: plabsorb,deltap
  INTEGER i,j,k,np

  DO np=1,mp
  IF(casamet%iveg2(np) == icewater) THEN
    casamet%glai(np)   = 0.0
  ELSE  
    casapool%cplant(np,:)  = casapool%cplant(np,:)  &
                           + casapool%dCplantdt(np,:)  * deltpool 
    casapool%clabile(np)   = casapool%clabile(np)   &
                           + casapool%dclabiledt(np)   * deltpool  
    IF(casapool%cplant(np,leaf) > 0.0) THEN
      IF(icycle >1) casapool%Nplant(np,:) = casapool%Nplant(np,:) &
                                 +casapool%dNplantdt(np,:)*deltpool
      IF(icycle >2) casapool%Pplant(np,:) = casapool%Pplant(np,:) &
                                 +casapool%dPplantdt(np,:)*deltpool
    ENDIF
!    casamet%glai(np)   = MIN(0.0, casabiome%sla(veg%iveg(np))  &
!                                  * casapool%cplant(np,leaf))
    casamet%glai(np)   = MAX(casabiome%glaimin(veg%iveg(np)), &
                               casabiome%sla(veg%iveg(np)) * casapool%cplant(np,leaf))
    casamet%glai(np)   = MIN(casabiome%glaimax(veg%iveg(np)), casamet%glai(np))
    casapool%clitter(np,:) = casapool%clitter(np,:) &
                           + casapool%dClitterdt(np,:) * deltpool 
    casapool%csoil(np,:)   = casapool%csoil(np,:)   &
                           + casapool%dCsoildt(np,:)   * deltpool

    IF(icycle >1) THEN
      casapool%Nlitter(np,:) = casapool%Nlitter(np,:) &
                             + casapool%dNlitterdt(np,:)* deltpool
      casapool%Nsoil(np,:)   = casapool%Nsoil(np,:)   &
                             + casapool%dNsoildt(np,:)  * deltpool
      casapool%Nsoilmin(np)  = casapool%Nsoilmin(np)  &
                             + casapool%dNsoilmindt(np) * deltpool
    ENDIF

    IF(icycle >2) THEN
      casapool%Plitter(np,:) = casapool%Plitter(np,:) &
                             + casapool%dPlitterdt(np,:)* deltpool
      casapool%Psoil(np,:)   = casapool%Psoil(np,:)   &
                             + casapool%dPsoildt(np,:)  * deltpool
      casapool%Psoillab(np)  = casapool%Psoillab(np)  &
                             + casapool%dPsoillabdt(np) * deltpool
      casapool%Psoilsorb(np) = casaflux%Psorbmax(np)*casapool%Psoillab(np) &
                             /(casaflux%kmlabp(np)+casapool%Psoillab(np))
!      casapool%Psoilsorb(np) = casapool%Psoilsorb(np)  &
!                             + casapool%dPsoilsorbdt(np) * deltpool
      casapool%Psoilocc(np)   = casapool%Psoilocc(np)  &
                              + casapool%dPsoiloccdt(np)  * deltpool
    ENDIF

    DO i=1,mplant
      IF(casapool%cplant(np,i) < 0.0)  THEN
        WRITE(57,*)  'Cpool: np,ivt',np,casamet%iveg2(np),casapool%cplant(np,:)
        call casa_poolzero(np,1,casapool)
        casapool%cplant(np,i) = max(0.0, casapool%cplant(np,i))
      ENDIF
    ENDDO
    IF(icycle >1) THEN
      DO i=1,mplant
        IF(casapool%nplant(np,i) < 0.0) THEN
          WRITE(57,*) 'Npool:', 'np,ivt,ipool',np,casamet%iveg2(np),casapool%nplant(np,:)
          call casa_poolzero(np,2,casapool)
          casapool%nplant(np,i) = max(0.0, casapool%nplant(np,i))
        ENDIF
      ENDDO
    ENDIF ! end of "icycle >1"

    DO j=1,mlitter
      IF(casapool%clitter(np,j) < 0.0)  THEN
        WRITE(57,*)  'Clitter: np,ivt2',np,casamet%iveg2(np),casapool%clitter(np,:)
        call casa_poolzero(np,3,casapool) 
        casapool%clitter(np,j) = max(0.0, casapool%clitter(np,j))
      ENDIF
    ENDDO

    DO k=1,msoil
      IF(casapool%csoil(np,k) < 0.0)    THEN
        WRITE(57,*)  'Csoil: np,ivt2',np,casamet%iveg2(np),casapool%csoil(np,:)
        call casa_poolzero(np,5,casapool) 
        casapool%csoil(np,k) = max(0.0, casapool%csoil(np,k))
      ENDIF
    ENDDO

!  check if any pool size, and terminate model run if any pool size is negative!!
    IF(icycle >1) THEN
      DO j=1,mlitter
        IF(casapool%nlitter(np,j) < 0.0)  THEN
          WRITE(57,*)  'Nlitter: np,ivt2',np,casamet%iveg2(np),casapool%Nlitter(np,:)
          call casa_poolzero(np,4,casapool) 
          casapool%nlitter(np,j) = max(0.0, casapool%nlitter(np,j))
        ENDIF
      ENDDO
      DO k=1,msoil
        IF(casapool%nsoil(np,k) < 0.0) THEN
          WRITE(57,*)  'Nsoil: np,ivt2',np,casamet%iveg2(np),casapool%nsoil(np,:)
          call casa_poolzero(np,6,casapool) 
          casapool%nsoil(np,k) = max(0.0, casapool%nsoil(np,k))
        ENDIF
      ENDDO
    ENDIF  !end of "icycle >1"
  ENDIF
  ENDDO !end of "np"

END SUBROUTINE casa_cnpcycle

SUBROUTINE casa_poolzero(n,ipool,casapool)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, ipool
  TYPE (casa_pool), INTENT(INOUT) :: casapool

  WRITE(57,*) ' WARNING: negative pools are reset to ZERO!!'
  SELECT CASE(ipool)
  CASE(1)
     WRITE(57,*) 'plant carbon pool size negative!!'
     WRITE(57,*) 'plant C pools: ', n,casapool%cplant(n,:)
  CASE(2)
     WRITE(57,*) 'plant nitrogen pool size negative!!'
     WRITE(57,*) 'plant C pools: ',n,casapool%cplant(n,:)
     WRITE(57,*) 'plant N pools: ',n,casapool%nplant(n,:)
  CASE(3)
     WRITE(57,*) 'litter carbon pool size negative!!'
     WRITE(57,*) 'litter C pools: ',n,casapool%clitter(n,:)
  CASE(4)
     WRITE(57,*) 'litter nitrogen pool size negative!!'
     WRITE(57,*) 'carbon pool: ',n,casapool%clitter(n,:)
     WRITE(57,*) 'nitrogen pools: ',n,casapool%nlitter(n,:)
  CASE(5)
     WRITE(57,*) 'soil carbon pool size negative!!'
     WRITE(57,*) 'soil C pools: ',n,casapool%csoil(n,:)
  CASE(6)
     WRITE(57,*) 'soil nitrogen pool size negative!!'
     WRITE(57,*) 'soil C pools: ', n,casapool%csoil(n,:)
     WRITE(57,*) 'soil N pools: ', n,casapool%nsoil(n,:)
  END SELECT 

END SUBROUTINE casa_poolzero

SUBROUTINE casa_cnpbal(casapool,casaflux,casabal)

  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_balance),          INTENT(INOUT) :: casabal

  ! local variables
  !INTEGER :: npt
  REAL(r_2), DIMENSION(mp) :: cbalplant,  nbalplant,  pbalplant
  REAL(r_2), DIMENSION(mp) :: cbalsoil,   nbalsoil,   pbalsoil
  !REAL(r_2), DIMENSION(mp) :: cbalplantx, nbalplantx, pbalplantx
        

  cbalplant(:) = 0.0
  cbalsoil(:)  = 0.0
  nbalplant(:) = 0.0
  nbalsoil(:)  = 0.0
  pbalplant(:) = 0.0
  pbalsoil(:)  = 0.0

  casabal%cbalance(:)  = 0.0
  casabal%nbalance(:)  = 0.0
  casabal%pbalance(:)  = 0.0

!C balance
   Cbalplant(:)  = sum(casabal%cplantlast,2) -sum(casapool%cplant,2)            &
                 + casabal%Clabilelast(:)-casapool%clabile(:)                   &        
                 +(casaflux%Cnpp(:) - SUM((casaflux%kplant*casabal%cplantlast),2))*deltpool &               
                 + casapool%dClabiledt(:)* deltpool
   Cbalsoil(:)   = sum(casabal%clitterlast,2) - sum(casapool%clitter,2)         &
                 + sum(casabal%csoillast,2)   - sum(casapool%csoil,2)           & 
                 +(SUM((casaflux%kplant*casabal%cplantlast),2)-casaflux%Crsoil(:))*deltpool

   casabal%cbalance(:) = Cbalplant(:) + Cbalsoil(:)

   ! test
!   casabal%cbalance(:) = sum(casapool%cplant,2)+casapool%clabile+sum(casapool%clitter,2)+sum(casapool%csoil,2) &
!                       - sum(casabal%cplantlast,2)-casabal%clabilelast-sum(casabal%clitterlast,2)-sum(casabal%csoillast,2) &
!                       - (casaflux%Cnpp(:)+casapool%dClabiledt(:)-casaflux%Crsoil(:))*deltpool
!
!   npt=243   
!
!   write(77,91) casabal%cbalance(npt),Cbalplant(npt),Cbalsoil(npt), &
!               casapool%cplant(npt,:),casabal%cplantlast(npt,:),casapool%dCplantdt(npt,:), &
!               casaflux%kplant(npt,:)*casabal%cplantlast(npt,:),                            &
!               casaflux%fraccalloc(npt,:), &
!               casaflux%Cnpp(npt)*casaflux%fraccalloc(npt,1)-casaflux%kplant(npt,1)*casabal%cplantlast(npt,1), &
!               casaflux%Cnpp(npt)*casaflux%fraccalloc(npt,2)-casaflux%kplant(npt,2)*casabal%cplantlast(npt,2), &
!               casaflux%Cnpp(npt)*casaflux%fraccalloc(npt,3)-casaflux%kplant(npt,3)*casabal%cplantlast(npt,3), &
!               casapool%clitter(npt,:),casabal%clitterlast(npt,:),casapool%dClitterdt(npt,:), &
!               sum(casaflux%kplant(npt,:)*casabal%cplantlast(npt,:)),sum(casaflux%FluxCtolitter(npt,:)), &
!               sum(casaflux%FluxCtolitter(npt,:)-casaflux%klitter(npt,:)*casabal%clitterlast(npt,:)),   &             
!               casapool%csoil(npt,:),casabal%csoillast(npt,:),casapool%dCsoildt(npt,:),                  &
!               sum(casaflux%FluxCtosoil(npt,:))-casaflux%Crsoil(npt),          &
!               sum(casapool%dCsoildt(npt,:)),casapool%csoil(npt,2)-casabal%csoillast(npt,2)-casapool%dCsoildt(npt,2), &
!               casaflux%ksoil(npt,2)*casabal%csoillast(npt,2), &
!               casapool%clabile(npt), casabal%clabilelast(npt),casapool%dClabiledt(npt), &
!               casaflux%Cnpp(npt), casaflux%Crsoil(npt)
!91 format('balance= ',100(f12.5,2x))

   casabal%cplantlast  = casapool%cplant
   casabal%clabilelast = casapool%clabile
   casabal%clitterlast = casapool%clitter
   casabal%csoillast   = casapool%csoil
   casabal%sumcbal     = casabal%sumcbal + casabal%cbalance
 
   IF(icycle >1) THEN
      Nbalplant(:) = sum(casabal%nplantlast,2) -sum(casapool%nplant,2)                  & 
                    +casaflux%Nminuptake(:) *deltpool
      Nbalsoil(:)  = -sum(casapool%nlitter,2)-sum(casapool%nsoil,2)                     &
                     -casapool%nsoilmin(:)+ casabal%nsoilminlast(:)                     &
                     + sum(casabal%nlitterlast,2)    + sum(casabal%nsoillast,2)         &
                     +(casaflux%Nmindep(:) + casaflux%Nminfix(:)- casaflux%Nminloss(:)                        &
                       -casaflux%Nminleach(:)-casaflux%Nupland(:)) * deltpool

      casabal%nbalance(:) = Nbalplant(:) + Nbalsoil(:)

      casabal%nplantlast  = casapool%nplant
      casabal%nlitterlast = casapool%nlitter
      casabal%nsoillast   = casapool%nsoil
      casabal%nsoilminlast= casapool%nsoilmin
      casabal%sumnbal     = casabal%sumnbal + casabal%nbalance

   ENDIF

   IF(icycle >2) THEN
      Pbalplant(:) = sum(casabal%Pplantlast,2) -sum(casapool%Pplant,2)                      &
                   + casaflux%Plabuptake(:) *deltpool
      Pbalsoil(:)  = -sum(casapool%Plitter,2)        - sum(casapool%Psoil,2)                      &
                     + sum(casabal%Plitterlast,2)    + sum(casabal%Psoillast,2)                   &
                   -casapool%psoillab(:)-casapool%psoilsorb(:)-casapool%psoilocc(:)               &
                   + casabal%psoillablast(:) + casabal%psoilsorblast(:) + casabal%psoilocclast(:) &
                   +(casaflux%Pdep(:) + casaflux%Pwea(:)                                          &            
                     -casaflux%Pleach(:)-casaflux%Pupland(:)                                      &
                     -casaflux%Ploss(:)) * deltpool

      casabal%pbalance(:) = pbalplant(:) + pbalsoil(:) 
       
      casabal%pplantlast   = casapool%pplant
      casabal%plitterlast  = casapool%plitter
      casabal%psoillast    = casapool%psoil
      casabal%psoillablast = casapool%psoillab
      casabal%psoilsorblast= casapool%psoilsorb
      casabal%psoilocclast = casapool%psoilocc
      casabal%sumpbal  = casabal%sumpbal + casabal%pbalance
   ENDIF

END SUBROUTINE casa_cnpbal

SUBROUTINE casa_ndummy(casapool)
  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool

  casapool%Nplant(:,:) = casapool%Cplant(:,:) * casapool%ratioNCplant(:,:)

END SUBROUTINE casa_ndummy

SUBROUTINE casa_pdummy(casapool)
  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool

  casapool%Pplant(:,:) = casapool%Nplant(:,:) / casapool%ratioNPplant(:,:)

END SUBROUTINE casa_pdummy

SUBROUTINE phenology(iday,veg,phen)
  IMPLICIT NONE
  INTEGER,              INTENT(IN)    :: iday
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (phen_variable), INTENT(INOUT) :: phen

  ! local variables (temprary)
  INTEGER :: np
  INTEGER, DIMENSION(mp)  :: days,days1to2, days2to3, days3to4, days4to1 

!  PRINT *, 'Within SUBROUTINE phenology, mp = ', mp 
  DO np=1,mp
    days1to2(np) = phen%doyphase(np,2) - phen%doyphase(np,1)
    days2to3(np) = phen%doyphase(np,3) - phen%doyphase(np,2)
    days3to4(np) = phen%doyphase(np,4) - phen%doyphase(np,3)
    days4to1(np) = phen%doyphase(np,1) - phen%doyphase(np,4)
    IF(days1to2(np) < 0) days1to2(np) = days1to2(np) +365
    IF(days2to3(np) < 0) days2to3(np) = days2to3(np) +365
    IF(days3to4(np) < 0) days3to4(np) = days3to4(np) +365
    IF(days4to1(np) < 0) days4to1(np) = days4to1(np) +365 
  ENDDO
  ! compute leaf phenology
  DO np=1,mp
    SELECT CASE(phen%phase(np))
      CASE(0)
        days(np) = iday - phen%doyphase(np,4)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days4to1(np)) phen%phase(np) =1
      CASE(1)
        days(np) = iday - phen%doyphase(np,1)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days1to2(np)) phen%phase(np) =2
      CASE(2)
        days(np) = iday - phen%doyphase(np,2)
        IF(days(np) <0) days(np) = days(np) + 365
        IF(days(np) > days2to3(np)) phen%phase(np) =3
      CASE(3)
        days(np) = iday - phen%doyphase(np,3)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days3to4(np)) phen%phase(np) =0
    END SELECT
  ENDDO

  WHERE(veg%iveg==1 .or. veg%iveg ==2 )
       phen%phase = 2
  ENDWHERE

END SUBROUTINE phenology

END MODULE casa_cnp_module
