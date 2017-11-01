!----------------------------------------------------------------------------------------------------
! FILE: corpse_cycle.f90
!
! Purpose: 
!   Subroutines for C cycling in the CORPSE model
!
!     SUBROUTINE corpse_delplant
!         - replaces casa_delplant when CORPSE model is are used
!     SUBROUTINE corpse_soil 
!         - calculate daily changes in litter, microbial, and SOM pools
!         - replaces casa_coeffsoil and casa_delsoil when corpse model is used
!
! Contact: Melannie Hartman
!          melannie@ucar.edu
!
! History:
!   01/04/2016 - Created
!   02/22/2016 - Added subroutine corpse_delplant
!   03/21/2016 - Run daily instead of hourly timestep
!----------------------------------------------------------------------------------------------------

MODULE corpse_cycle_module

USE define_types
USE casavariable
USE casaparm
USE corpsedimension
USE corpsevariable
IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------------

SUBROUTINE corpse_delplant(mp,veg,casabiome,casapool,casaflux,casamet,            &
                           cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                           nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                           pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                           cwd2co2,cwd2str)


! Calculate the change in plant C pools when the CORPSE model is used

  IMPLICIT NONE
  integer, INTENT(IN)                      :: mp         ! number of grid cells
  TYPE (veg_parameter_type), INTENT(IN)    :: veg        ! vegetation parameters
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

  casaflux%FluxCtolitter = 0.0
  casaflux%FluxNtolitter = 0.0
  casaflux%FluxPtolitter = 0.0
  casapool%dClitterdt(:,:) = 0.0;

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
          cwd2co2(npt) =  casaflux%fromLtoCO2(npt,cwd) * casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)
          cwd2str(npt) = (1.0-casaflux%fromLtoCO2(npt,cwd)) * casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)
          ! Calculate change in CWD pool.  CWD pool is updated in mimics_ccyle (-mdh 5/4/2015)
          casapool%dClitterdt(npt,cwd) = cwood2cwd(npt) - casaflux%klitter(npt,cwd) * casapool%clitter(npt,cwd)

          casaflux%ClitInptMet(npt) = cleaf2met(npt) + croot2met(npt)
          ! Add cwd2str instead of cwood2cwd to structural litter input. -mdh 1/23/2017
          !casaflux%ClitInptStruc(npt) = cleaf2str(npt) + croot2str(npt)+ cwood2cwd(npt)
          casaflux%ClitInptStruc(npt) = cleaf2str(npt) + croot2str(npt)+ cwd2str(npt)

      ENDIF
  ENDDO
 
END SUBROUTINE corpse_delplant


!----------------------------------------------------------------------------------------------------
! Calculate daily changes in litter, microbial, and SOM pools
SUBROUTINE corpse_soil(mp,idoy,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2)

  IMPLICIT NONE
  ! Function Arguments
  integer, INTENT(IN) :: mp,idoy  ! number of grid points, day of year
  !Daily Litter Inputs (gC/m2)
  real, dimension(mp),INTENT(IN) :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2

  ! Local Variables
  integer :: npt, ihr, jj, doy
  integer :: badLitter
  real    :: zeroThreshold
  real :: air_filled_porosity, theta_frzn
  real,dimension(mp) :: theta

  ! Variables to hold the daily litter inputs from the CASACNP model
  real,dimension(nspecies):: daily_leaflitter_input
  real,dimension(nspecies):: daily_rootlitter_input
  real,dimension(nspecies):: daily_exudate_input

  zeroThreshold = -10.0e-5
  theta(:) = 0.0

  if (iptToSave_corpse > 0) then
      open(215,file=sPtFileNameCORPSE, access='APPEND')
      close(215)
  endif

  do npt=1,mp
  
      badLitter = 0
      ! Only the daily_exudate_input(LABILE) value will be reset to non-zero value. -mdh 2/8/2017
      daily_exudate_input(:) = 0

      ! Ben said that theta is actually fraction of water-filled pore space, not volumetric swc. -mdh 2/27/2017
      !!theta = casamet%moistavg(npt) ! mean volumetric soil water content (0.0 - 1.0)
      theta(npt) = casamet%moistavg(npt)/soil%ssat(npt) ! fraction of liquid water-filled pore space (0.0 - 1.0)

      IF(casamet%iveg2(npt) /= icewater) THEN

          ! New output variable to track daily CO2 respiration losses when cwd decomposes to structural litter. -mdh 2/6/2017
          if (idoy == 1) then
              !pt(npt)%litterlayer_outputs%cwd2co2(1..366) 
              pt(npt)%litterlayer_outputs%cwd2co2(:) = 0.0
          endif

          ! ONCE A DAY
  
          T = casamet%tsoilavg(npt)     ! degrees K

          !! Ben Sulman suggested that exudation is a function of NPP rather than a fixed input (-mdh 11/14/2016)
          !! daily_exudate_input = annual_exudate_input*dt
          !! Divide by 1000 to convert gC/m2/yr to kgC/m2/yr. -mdh 11/18/2016
          daily_exudate_input(LABILE) = exudate_npp_frac(LABILE)*casaflux%CnppAn(npt)*dt/1000.0
   
          ! Convert litter inputs from gC/m2 to kgC/m2
          daily_leaflitter_input(LABILE) = cleaf2met(npt)/1000.0
          daily_leaflitter_input(RECALCTRNT) = (cleaf2str(npt) + cwd2str(npt))/1000.0
          daily_leaflitter_input(DEADMICRB) = 0.0
  
          daily_rootlitter_input(LABILE) = croot2met(npt)/1000.0
          daily_rootlitter_input(RECALCTRNT) = croot2str(npt)/1000.0
          daily_rootlitter_input(DEADMICRB) = 0.0

          !Reduce labile root litter inputs by the amount of daily_exudate_input.  
          !Reduce daily_exudate_input if the flux exceeds this labile litter input. - mdh 2/8/2017
          if ((daily_rootlitter_input(LABILE) - daily_exudate_input(LABILE)) > 0.0) then
              daily_rootlitter_input(LABILE) = daily_rootlitter_input(LABILE) - daily_exudate_input(LABILE)
              !write(*,*) 'npt:', npt, ' daily_exudate_input(LABILE) =', daily_exudate_input(LABILE)
          else
              !write(*,*) 'Reducing daily_exudate_input by: ', daily_exudate_input(LABILE) - daily_rootlitter_input(LABILE)
              !write(*,*) '  npt:', npt, ' daily_exudate_input =', daily_exudate_input(LABILE)
              !write(*,*) '  npt:', npt, ' daily_rootlitter_input(LABILE) =', daily_rootlitter_input(LABILE)
              daily_exudate_input(LABILE) = daily_rootlitter_input(LABILE)
              daily_rootlitter_input(LABILE) = 0.0
          endif

          if (daily_leaflitter_input(LABILE) < zeroThreshold) then
              write(*,*) 'npt:', npt, ' daily_leaflitter_input(LABILE) =', daily_leaflitter_input(LABILE)
              daily_leaflitter_input(LABILE) = 0.0
              badLitter = badLitter + 1
          endif
          if (daily_leaflitter_input(RECALCTRNT) < zeroThreshold) then
              write(*,*) 'npt:', npt, ' daily_leaflitter_input(RECALCTRNT) =', daily_leaflitter_input(RECALCTRNT)
              daily_leaflitter_input(RECALCTRNT) = 0.0
              badLitter = badLitter + 1
          endif
          if (daily_rootlitter_input(LABILE) < zeroThreshold) then
              write(*,*) 'npt:', npt, ' daily_rootlitter_input(LABILE) =', daily_rootlitter_input(LABILE)
              daily_rootlitter_input(LABILE) = 0.0
              badLitter = badLitter + 1
          endif
          if (daily_rootlitter_input(RECALCTRNT) < zeroThreshold) then
              write(*,*) 'npt:', npt, ' daily_rootlitter_input(RECALCTRNT) =', daily_rootlitter_input(RECALCTRNT)
              daily_rootlitter_input(RECALCTRNT) = 0.0
              badLitter = badLitter + 1
          endif       
  
          if (badLitter >= 1) then
              !!write(*,*) 'npt:', npt, ' Bad Litter Input on time/day ', timestep/(365*24)+1.0, idoy
              write(*,*) 'npt:', npt, ' Bad Litter Input on time/day ', timestep/(365)+1.0, idoy
              !!STOP
          endif
  
          call add_litter2(pt(npt)%litterlayer, daily_leaflitter_input)
          do jj=1,num_lyr
              !Add root litter to soil
              call add_litter2(pt(npt)%soil(jj), &
                              daily_rootlitter_input, &
                              rhizosphere_frac)
          enddo

          ! HOURLY
          !!do ihr = 1, NHOURS - commented out to run daily (-mdh 3/21/2016)
  
          !Loop through soil layers, calling update_pool to run decomposition etc
  
              do jj=1,num_lyr
      
                  !Add exudate carbon to rhizosphere (more continuous in time).
  
                  call add_carbon_to_rhizosphere(pt(npt)%soil(jj),daily_exudate_input)
          
                  ! Update using casamet%frznmoistavg(npt) variable. -mdh 3/13/2107
                  !air_filled_porosity = 1.0-theta(npt)
                  theta_frzn = casamet%frznmoistavg(npt)/soil%ssat(npt) ! fraction of frozen water-filled pore space (0.0 - 1.0)
                  air_filled_porosity = max(0.0, 1.0-theta(npt)-theta_frzn)

                  call update_pool(pool=pt(npt)%soil(jj), &
                               T=T, &
                               theta=theta(npt),&
                               air_filled_porosity=air_filled_porosity, &
                               liquid_water=theta(npt), &
                               frozen_water=theta_frzn, &
                               dt=dt,&
                               layerThickness=pt(npt)%dz(jj), &
                               fast_C_loss_rate=pt(npt)%fast_C_loss_rate, &
                               slow_C_loss_rate=pt(npt)%slow_C_loss_rate, &
                               deadmic_C_loss_rate=pt(npt)%deadmic_C_loss_rate, &
                               CO2prod=pt(npt)%CO2prod, &
                               deadmic_produced=pt(npt)%deadmic_produced, &
                               protected_produced=pt(npt)%protected_produced, &
                               protected_turnover_rate=pt(npt)%protected_turnover_rate, &
                               C_dissolved=pt(npt)%C_dissolved, &
                               deposited_C=pt(npt)%deposited_C, &
                               npt=npt, &
                               doy=idoy, &
                               hr=ihr)
              enddo
      
              ! Update using casamet%frznmoistavg(npt) variable. -mdh 3/13/2107
              !air_filled_porosity = 1.0-theta(npt)
              theta_frzn = casamet%frznmoistavg(npt)/soil%ssat(npt) ! fraction of frozen water-filled pore space (0.0 - 1.0)
              air_filled_porosity = max(0.0, 1.0-theta(npt)-theta_frzn)

              !Do the decomposition etc for the litter layer
              call update_pool(pool=pt(npt)%litterlayer,&
                           T=T, &
                           theta=theta(npt), &
                           air_filled_porosity=air_filled_porosity, &
                           liquid_water=theta(npt), &
                           frozen_water=theta_frzn, &
                           dt=dt,&
                           layerThickness=pt(npt)%dz(1),&
                           fast_C_loss_rate=pt(npt)%fast_C_loss_rate, &
                           slow_C_loss_rate=pt(npt)%slow_C_loss_rate, &
                           deadmic_C_loss_rate=pt(npt)%deadmic_C_loss_rate, &
                           CO2prod=pt(npt)%CO2prod, &
                           deadmic_produced=pt(npt)%deadmic_produced, &
                           protected_produced=pt(npt)%protected_produced, &
                           protected_turnover_rate=pt(npt)%protected_turnover_rate, &
                           C_dissolved=pt(npt)%C_dissolved, &
                           deposited_C=pt(npt)%deposited_C, &
                           npt=npt, &
                           doy=idoy, &
                           hr=ihr)
      
          !!enddo   ! End hourly loop
  

          ! Save the CO2 flux from the transformation of CWD to structural litter. This flux 
          ! is not included to the decomposition of litter since it comes from the CASA model. 
          ! The cwd2co2 values will be added to the output heterotrophic respiration fluxes 
          ! in WritePoolFluxNcFile_corpse. Convert cwd2co2 to kg/m2. -mdh 2/6/2017
          pt(npt)%litterlayer_outputs%cwd2co2(idoy) = cwd2co2(npt)/1000.0

      ENDIF

  enddo   ! End loop through npt grid cells


  !!timestep = timestep + NHOURS  - commented out for daily timestep (-mdh 3/21/2016)
  timestep = timestep + 1
  doy = MOD(timestep, 365)   ! Day of Year
  if (doy == 0) then
      doy = 365
  end if

  do npt=1,mp

      !Record the output if time step matches recording step
      if (mod(timestep,recordtime) .eq. 0) then
          do jj=1,num_lyr
              !call save_output_line(pt(npt)%soil(jj), pt(npt)%soil_outputs(jj), casamet%moistavg(npt), casamet%tsoilavg(npt))
              call save_output_line(pt(npt)%soil(jj), pt(npt)%soil_outputs(jj), theta(npt), casamet%tsoilavg(npt))
          enddo
          !call save_output_line(pt(npt)%litterlayer, pt(npt)%litterlayer_outputs, casamet%moistavg(npt), casamet%tsoilavg(npt))
          call save_output_line(pt(npt)%litterlayer, pt(npt)%litterlayer_outputs, theta(npt), casamet%tsoilavg(npt))

          if (casamet%ijgcm(npt) .eq. iptToSave_corpse) then
              !ATTENTION: TO DO
              !write(*,*) 'litterlayer_totalC(', npt, ')= ', &
              !            pt(npt)%litterlayer_outputs%totalC(pt(npt)%litterlayer_outputs%linesWritten)
              !write(*,*) 'soil_totalC(', npt, ')= ', &
              !            pt(npt)%soil_outputs(1)%totalC(pt(npt)%soil_outputs(1)%linesWritten)
          endif
      endif

      ! Accumulate C pools and fluxes every day
      do jj=1,num_lyr
          call corpse_caccum(pt(npt)%soil(jj), pt(npt)%soil_outputs(jj), theta(npt), casamet%tsoilavg(npt), &
                             doy,casamet%ijgcm(npt))
      enddo
      call corpse_caccum(pt(npt)%litterlayer, pt(npt)%litterlayer_outputs, theta(npt), casamet%tsoilavg(npt), &
                         doy,casamet%ijgcm(npt))

  enddo

end SUBROUTINE corpse_soil

!--------------------------------------------------------------------------------

END MODULE corpse_cycle_module
