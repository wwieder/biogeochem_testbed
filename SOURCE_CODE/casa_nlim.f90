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
! Purpose: subroutines for calculating N limitation on litter decay and
!          plant growth
!
! Called from: biogeochem (mostly) or casa_xnp
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed by Yingping Wang (Wang et al., BG, 2011)
!
! ==============================================================================
! casa_nlim.f90
!
! This module contains the following subroutines:
!   casa_xkN (moved from casa_cnp.f90 on 11/11/2019 -mdh)
!   casa_xkN2 (alternative to casa_xkN, -mdh 11/11/2019)

MODULE casa_nlim_module
USE define_dimensions
USE define_types
USE casadimension
USE casaparm
USE casavariable
IMPLICIT NONE
CONTAINS

SUBROUTINE casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
! computing the reduction in litter and SOM decomposition 
! when decomposition rate is N-limiting
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(INOUT) :: xkNlimiting
  TYPE (casa_pool),         INTENT(INOUT) :: casapool
  TYPE (casa_flux),         INTENT(INOUT) :: casaflux
  TYPE (casa_met),          INTENT(INOUT) :: casamet
  TYPE (casa_biome),        INTENT(INOUT) :: casabiome
!    
  TYPE (veg_parameter_type),   INTENT(IN) :: veg  ! vegetation parameters

  ! local variables
  INTEGER j,k,kk,nland
  REAL(r_2), DIMENSION(mp)         :: xFluxNlittermin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilmin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilimm
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilminnet
! A maximum Clitter set to avoid positive feedback for litter accumulation
! when N mode is activated. (Q.Zhang 23/05/2011)
!  real(r_2), dimension(17)         :: xClitter
!  data xClitter/100.0,100.0,100.0,100.0,50.0,150.0,150.0,100.0,&
!                150.0,150.0,100.0, 20.0,20.0, 20.0, 20.0, 20.0,20.0/

  xkNlimiting(:)  = 1.0
!  set N mineral N fluxes to zero
  xFluxNlittermin(:)  = 0.0
  xFluxNsoilmin(:)    = 0.0
  xFluxNsoilimm(:)    = 0.0  !negative for microbial uptake and positive for release of mineral N
  xFluxNsoilminnet(:) = 0.0
!   PRINT *, 'within casa_xkN'

!  calculate gross mineralisation
  DO nland=1,mp
  IF (casamet%iveg2(nland)/=icewater) THEN

    ! calculate C:N ratio of newly formed SOM as function of soil mineral N pool
    IF (casapool%Nsoilmin(nland) < 2.0) THEN
      casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmin(nland,:)  &
                                       + (casapool%ratioNCsoilmax(nland,:) &
                                         -casapool%ratioNCsoilmin(nland,:)) &
                                       * max(0.0,casapool%Nsoilmin(nland)) / 2.0
    ELSE
      casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmax(nland,:)
    ENDIF

    DO j=1,mlitter
      xFluxNlittermin(nland) = xFluxNlittermin(nland) &
                         + casaflux%klitter(nland,j) * casapool%Nlitter(nland,j)
    ENDDO
    DO k=1,msoil
      xFluxNsoilmin(nland)   = xFluxNsoilmin(nland)   &
                         + casaflux%ksoil(nland,k)   * casapool%Nsoil(nland,k)
    ENDDO

    ! calculate N immobilisation from L to S and S to S
    DO kk=1,msoil
      DO j=1,mlitter    ! immobilisation from litter to soil
        xFluxNsoilimm(nland) = xFluxNsoilimm(nland) &
             - casaflux%fromLtoS(nland,kk,j) * casaflux%klitter(nland,j) &
             * casapool%Clitter(nland,j) * casapool%ratioNCsoilnew(nland,kk)
      ENDDO
      DO k=1,msoil      ! immobilisation from soil to soil
        IF(k.ne.kk) THEN
          xFluxNsoilimm(nland) = xFluxNsoilimm(nland) &
               - casaflux%fromStoS(nland,kk,k) * casaflux%ksoil(nland,k) &
               * casapool%Csoil(nland,k) * casapool%ratioNCsoilnew(nland,kk)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  ENDDO

  ! now check if there is sufficient mineral N 
  xFluxNsoilminnet(:) = xFluxNlittermin(:) + xFluxNsoilmin(:) + xFluxNsoilimm(:)
!   PRINT *, 'casamet%iveg2 = ', casamet%iveg2
!   PRINT *, 'deltpool = ',deltpool
!   PRINT *, 'xFluxNsoilminnet = ', xFluxNsoilminnet
! WHERE(casamet%iveg2(:)/=icewater) 
!    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0) 
!      xkNlimiting(:) =1.0
!    ELSEWHERE
!      xkNlimiting(:) =max(0.0, - (casapool%Nsoilmin(:)-2.0)/(deltpool*xFluxNsoilminnet(:))) 
!      xkNlimiting(:) =MIN(1.0,xkNlimiting(:))
!    ENDWHERE
! ENDWHERE 

! Q.Zhang 23/05/2011 test code according to YPW
  WHERE(casamet%iveg2(:)/=icewater)
!! Not getting any N limitation.  Drop the OR condition for testing. -mdh 10/7/2019
!!    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0 &
!!          .OR. xFluxNsoilminnet(:) > 0.0)
    WHERE((xFluxNsoilminnet(:)*deltpool + (casapool%Nsoilmin(:)-2.0)) > 0.0) 
      xkNlimiting(:) =1.0
    ELSEWHERE
      xkNlimiting(:) =MAX(0.0, - (casapool%Nsoilmin(:)-0.5) &
                                /(deltpool*xFluxNsoilminnet(:)))
      xkNlimiting(:) =MIN(1.0,xkNlimiting(:))
    ENDWHERE
! Q.Zhang 23/05/2011 test
! If pool size larger than xClitter, turnover rate will not constrained by Nsoilmin.
!    where(casapool%clitter(:,1) > xClitter(veg%iveg(:)))
!     xkNlimiting(:) = 1.0
!    end where
! end (Q.Zhang 23/05/2011)
    where(sum(casapool%clitter,2) > casabiome%maxfinelitter(veg%iveg(:)) + casabiome%maxcwd(veg%iveg(:)))
     xkNlimiting(:) = 1.0
    end where
  ENDWHERE

  ! Added casaflux%xkNlimiting for output. -mdh 10/7/2019
! write(*,*)
! write(*,*) 'xFluxNsoilminnet(:)*deltpool =', xFluxNsoilminnet*deltpool 
! write(*,*) 'casapool%Nsoilmin(:)-2.0 =', casapool%Nsoilmin-2.0
! write(*,*) 'xFluxNsoilminnet(:) =', xFluxNsoilminnet
  casaflux%xkNlimiting(:) = xkNlimiting(:)

END SUBROUTINE casa_xkN

!======================================================================================================
! Subroutine casa_xkN2 is similar to subroutine casa_xkN, but it does a simple linear ramp 
! of xkNlimiting between 0.0 and 1.0 when xkNlimitmin < casapool%Nsoilmin < xkNlimitmax. 
! (formerly when 0.5 < casapool%Nsoilmin < 2.0). 
! Melannie Hartman 11/11/2019

! Added isomModel argument to the casa_xkN2 function call. -mdh 3/23/2020.
!!SUBROUTINE casa_xkN2(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
SUBROUTINE casa_xkN2(xkNlimiting,casapool,casaflux,casamet,casabiome,veg,isomModel)
! computing the reduction in litter and SOM decomposition 
! when decomposition rate is N-limiting
  IMPLICIT NONE
  REAL(r_2), DIMENSION(mp), INTENT(INOUT) :: xkNlimiting
  TYPE (casa_pool),         INTENT(INOUT) :: casapool
  TYPE (casa_flux),         INTENT(INOUT) :: casaflux
  TYPE (casa_met),          INTENT(INOUT) :: casamet
  TYPE (casa_biome),        INTENT(INOUT) :: casabiome
  INTEGER,                  INTENT(IN)    :: isomModel
!    
  TYPE (veg_parameter_type),   INTENT(IN) :: veg  ! vegetation parameters

  ! local variables
  INTEGER j,k,kk,nland
  REAL(r_2), DIMENSION(mp)         :: xFluxNlittermin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilmin
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilimm
  REAL(r_2), DIMENSION(mp)         :: xFluxNsoilminnet
! A maximum Clitter set to avoid positive feedback for litter accumulation
! when N mode is activated. (Q.Zhang 23/05/2011)

  xkNlimiting(:) = 1.0
!   PRINT *, 'within casa_xkN2'

  DO nland=1,mp
    IF (casamet%iveg2(nland)/=icewater) THEN
      ! calculate C:N ratio of newly formed SOM as function of soil mineral N pool
      ! This code is only relevant to CASA SOM model. -mdh 3/23/2020.
      if (isomModel == CASACNP) then

        !!IF (casapool%Nsoilmin(nland) < 2.0) THEN
        !!  ratioNCnew = ratioNCmin + (ratioNCmax-ratioNCmin)*max(0.0,Nsoilmin)/2.0
        !!  casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmin(nland,:)  &
        !!                                     + (casapool%ratioNCsoilmax(nland,:) &
        !!                                     - casapool%ratioNCsoilmin(nland,:)) &
        !!                                     * max(0.0,casapool%Nsoilmin(nland)) / 2.0

        ! Substitute 2.0 with xkNlimitmax. -mdh 3/23/2020.
        IF (casapool%Nsoilmin(nland) < casabiome%xkNlimitmax(veg%iveg(nland))) THEN
          !ratioNCnew = ratioNCmin + (ratioNCmax-ratioNCmin)*max(0.0,Nsoilmin)/xkNlimitmax
          casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmin(nland,:)  &
                                             + (casapool%ratioNCsoilmax(nland,:) &
                                             - casapool%ratioNCsoilmin(nland,:)) &
                                             * max(0.0,casapool%Nsoilmin(nland)) / casabiome%xkNlimitmax(veg%iveg(nland))
        ELSE
          casapool%ratioNCsoilnew(nland,:) = casapool%ratioNCsoilmax(nland,:)
        ENDIF
      endif
    ENDIF
  ENDDO

  ! Now check if there is sufficient mineral N 
  ! Use parameters from CASA PFT parameter file. -mdh 12/23/2019
  ! casabiome%xkNlimitmin(veg%iveg(:)) replaces 0.5
  ! casabiome%xkNlimitmax(veg%iveg(:)) replaces 2.0
  WHERE(casamet%iveg2(:)/=icewater)
    !WHERE((casapool%Nsoilmin(:)-2.0) > 0.0) 
    WHERE((casapool%Nsoilmin(:)-casabiome%xkNlimitmax(veg%iveg(:))) > 0.0) 
      xkNlimiting(:) = 1.0
    ELSEWHERE
      ! Linear ramp from 0.0 to 1.0 for 0.5 < casapool%Nsoilmin(:) < 2.0
      ! y = (y2 - y1) / (x2 - x1) * (x - x2) + y2

      !xkNlimiting(:) = (1.0 - 0.0) / (2.0 - 0.5) * (casapool%Nsoilmin(:) - 2.0) + 1.0
      xkNlimiting(:) = (1.0 - 0.0) / (casabiome%xkNlimitmax(veg%iveg(:)) - casabiome%xkNlimitmin(veg%iveg(:))) &
                       * (casapool%Nsoilmin(:) - casabiome%xkNlimitmax(veg%iveg(:))) + 1.0
      xkNlimiting(:) = MAX(0.0, xkNlimiting(:)) 
      xkNlimiting(:) = MIN(1.0, xkNlimiting(:))
    ENDWHERE
  ENDWHERE

  ! The following code is only relevant to the CASA SOM model. -mdh 3/23/2020.
  if (isomModel == CASACNP) then
    ! If Clitter pool (metb+struc+cwd) size larger > maximum, turnover rate (klitter) will not constrained by Nsoilmin.
    WHERE(sum(casapool%clitter,2) > casabiome%maxfinelitter(veg%iveg(:)) + casabiome%maxcwd(veg%iveg(:)))
      xkNlimiting(:) = 1.0
    ENDWHERE
  endif

  ! Added casaflux%xkNlimiting for output. -mdh 10/7/2019
  casaflux%xkNlimiting(:) = xkNlimiting(:)

END SUBROUTINE casa_xkN2
!======================================================================================================

END MODULE casa_nlim_module
