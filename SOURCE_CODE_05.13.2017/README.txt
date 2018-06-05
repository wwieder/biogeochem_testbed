module load tool/netcdf/4.3.2/gcc
rm *.o *.mod
make -f Makefile.txt casaclm_mimics_corpse_test

Code Updates 11/07/2016
Melannie Hartman

----------------------------------------------------------------------------------------------------
I was getting this message when running ncrcat: 
  "ERROR: units attribute "simulation time in years" is not listed in UDUnits2 SI system database"
  I set attr_units = '1..ntimes' for the time variable (as I originally had it) instead of 'simulation time in years'.
  This update was made to casa_input.f90 and mimics_inout.f90, but corpse_inout.f90 did not need to be updated.

----------------------------------------------------------------------------------------------------

Will noticed that the protected soil pool was growing unrealistically large for grid cellid = 01110 (middle of Africa).
This is a very dry grid cell.  After consulting over email with Ben Sulman, we decided to modify the respiration calculation 
in subroutine update_cohort (corpse_soil_carbon.f90):

    ! Put a lower limit of 0.001 on theta^3 (w.wieder, 11/7/2016).
!   Resp=Vmax(T)*vmax_multiplier*theta**3*(Cavail)*enz/(sum(Cavail)*kC+enz) &
!        *max((air_filled_porosity)**gas_diffusion_exp,min_anaerobic_resp_factor)
    Resp=Vmax(T)*vmax_multiplier*(theta**3+0.001)*(Cavail)*enz/(sum(Cavail)*kC+enz) &
         *max((air_filled_porosity)**gas_diffusion_exp,min_anaerobic_resp_factor)

----------------------------------------------------------------------------------------------------
Will found a bug in mimics_coeffplant.  The fmet term was not based on the fmet_p(*) parameters as it is in
mimics_delplant (file mimics_cycle.f90).

fmet = mimicsbiome%fmet_p(1)*(mimicsbiome%fmet_p(2)-mimicsbiome%fmet_p(3) * mimicsbiome%ligninNratio(:,leaf))

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

----------------------------------------------------------------------------------------------------
Code Updates 11/14/2016
Melannie Hartman

CORPSE now reads exudate_npp_frac from corpse_params.nml instead of annual_exudation_input (kgC/m2/yr).
CORPSE uses exudate_npp_frac to compute root_exudation as a fraction of annual NPP from the previous 
year insted of using a fixed amount which was 0.01, 0.0, and 0.0 for labile, recalcitrant, and dead 
microbes, respectively.  

I'm in the process of outputing structural and metabolic litterfall (inputs) to the casa netcdf files.
These inputs apply to all models, so it is best to put it here rather than in MIMICS and CORPSE netcdf files 
seperately. MIMICS does have litter inputs in its netcdf files so I need to check for consistency, then I can 
potentially delete them from MIMICS netcdf files.

----------------------------------------------------------------------------------------------------
Code Updates 11/21/2016
Melannie Hartman

Will showed me graphs of CORPSE litter layer (above-ground) C at the beginning of the transient.  
The LitterLayer_C2 pool (the majority of litter layer C) declines rapidly at the beginning of the 
transient run, but soil pools are more stable.  Will wondered if there was initialization problem 
at restart.   The total Litter Layer C pool is globally the same as at the end of the spinup and the 
beginning of the transient. There was no change in litter inputs between the too simulations either.  
However, respiration from the litter layer is high at the beginning of the transient run.  

There was an initialization problem with the litter layer at restart.  I was portioning the litter layer into 
bulk and rhizosphere after reading the restart file, and for the litter layer(above-ground) there is no rhizosphere!
Rhizoshere litter would decompose faster.


**********
Code before (subroutine corpse_init in corpse_inout.f90)

          do jj = 1, nspecies
              pt(ii)%litterlayer%litterCohorts(RHIZ)%litterC(jj)= litter_litterC(jj)*rhizosphere_frac
              pt(ii)%litterlayer%litterCohorts(RHIZ)%protectedC(jj) = litter_protectedC(jj)*rhizosphere_frac
              pt(ii)%soil(lyr)%litterCohorts(RHIZ)%litterC(jj) = soil_litterC(jj)*rhizosphere_frac
              pt(ii)%soil(lyr)%litterCohorts(RHIZ)%protectedC(jj) = soil_protectedC(jj)*rhizosphere_frac
              pt(ii)%litterlayer%litterCohorts(BULK)%litterC(jj)= litter_litterC(jj)*(1.0-rhizosphere_frac)
              pt(ii)%litterlayer%litterCohorts(BULK)%protectedC(jj) = litter_protectedC(jj)*(1.0-rhizosphere_frac)
              pt(ii)%soil(lyr)%litterCohorts(BULK)%litterC(jj) = soil_litterC(jj)*(1.0-rhizosphere_frac)
              pt(ii)%soil(lyr)%litterCohorts(BULK)%protectedC(jj) = soil_protectedC(jj)*(1.0-rhizosphere_frac)
          end do 

          pt(ii)%litterlayer%litterCohorts(RHIZ)%livingMicrobeC =  litter_livingMicrobeC*rhizosphere_frac
          pt(ii)%soil(lyr)%litterCohorts(RHIZ)%livingMicrobeC =  soil_livingMicrobeC*rhizosphere_frac
          pt(ii)%litterlayer%litterCohorts(BULK)%livingMicrobeC =  litter_livingMicrobeC*(1.0-rhizosphere_frac)
          pt(ii)%soil(lyr)%litterCohorts(BULK)%livingMicrobeC =  soil_livingMicrobeC*(1.0-rhizosphere_frac)


**********
Code after (subroutine corpse_init in corpse_inout.f90)

          ! Rhizosphere fraction = 0.0 for the litterlayer. -mdh 11/21/2106
          do jj = 1, nspecies
              pt(ii)%litterlayer%litterCohorts(RHIZ)%litterC(jj)= 0.0
              pt(ii)%litterlayer%litterCohorts(RHIZ)%protectedC(jj) = 0.0
              pt(ii)%litterlayer%litterCohorts(BULK)%litterC(jj)= litter_litterC(jj)
              pt(ii)%litterlayer%litterCohorts(BULK)%protectedC(jj) = litter_protectedC(jj)

              pt(ii)%soil(lyr)%litterCohorts(RHIZ)%litterC(jj) = soil_litterC(jj)*rhizosphere_frac
              pt(ii)%soil(lyr)%litterCohorts(RHIZ)%protectedC(jj) = soil_protectedC(jj)*rhizosphere_frac
              pt(ii)%soil(lyr)%litterCohorts(BULK)%litterC(jj) = soil_litterC(jj)*(1.0-rhizosphere_frac)
              pt(ii)%soil(lyr)%litterCohorts(BULK)%protectedC(jj) = soil_protectedC(jj)*(1.0-rhizosphere_frac)
          end do 

          pt(ii)%litterlayer%litterCohorts(RHIZ)%livingMicrobeC =  0.0
          pt(ii)%litterlayer%litterCohorts(BULK)%livingMicrobeC =  litter_livingMicrobeC

          pt(ii)%soil(lyr)%litterCohorts(RHIZ)%livingMicrobeC =  soil_livingMicrobeC*rhizosphere_frac
          pt(ii)%soil(lyr)%litterCohorts(BULK)%livingMicrobeC =  soil_livingMicrobeC*(1.0-rhizosphere_frac)

**********
The above code update reduced the litterlayer_C2 losses at the beginning of the transient, but there is 
still another startup issue.

----------------------------------------------------------------------------------------------------

Code Updates 11/28/2016
Melannie Hartman

Ignore CORPSE updates from 11/21/2016. Read/Write end-of-simulation CORPSE pools directly by cohort

          read(107,*) ipt, & 
                      ijgcm, &
                      tcnt, &
                      ftime, &
                      latz, &
                      lonz, &
                      ivtz, &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%litterC(LABILE), &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%litterC(RECALCTRNT), &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%litterC(DEADMICRB), &
                      pt(ii)%litterlayer%litterCohorts(BULK)%litterC(LABILE), &
                      pt(ii)%litterlayer%litterCohorts(BULK)%litterC(RECALCTRNT), &
                      pt(ii)%litterlayer%litterCohorts(BULK)%litterC(DEADMICRB), &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%protectedC(LABILE), &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%protectedC(RECALCTRNT), &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%protectedC(DEADMICRB), &
                      pt(ii)%litterlayer%litterCohorts(BULK)%protectedC(LABILE), &
                      pt(ii)%litterlayer%litterCohorts(BULK)%protectedC(RECALCTRNT), &
                      pt(ii)%litterlayer%litterCohorts(BULK)%protectedC(DEADMICRB), &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%livingMicrobeC, &
                      pt(ii)%litterlayer%litterCohorts(BULK)%livingMicrobeC, &
                      pt(ii)%litterlayer%litterCohorts(RHIZ)%CO2, &
                      pt(ii)%litterlayer%litterCohorts(BULK)%CO2, &
                      pt(ii)%litterlayer%protection_rate, &
                      pt(ii)%litterlayer%Qmax, &
                      pt(ii)%litterlayer%dissolved_carbon(LABILE), &
                      pt(ii)%litterlayer%dissolved_carbon(RECALCTRNT), &
                      pt(ii)%litterlayer%dissolved_carbon(RECALCTRNT), &

                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%litterC(LABILE), &
                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%litterC(RECALCTRNT), &
                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%litterC(DEADMICRB), &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%litterC(LABILE), &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%litterC(RECALCTRNT), &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%litterC(DEADMICRB), &
                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%protectedC(LABILE), &
                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%protectedC(RECALCTRNT), &
                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%protectedC(DEADMICRB), &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%protectedC(LABILE), &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%protectedC(RECALCTRNT), &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%protectedC(DEADMICRB), &
                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%livingMicrobeC, &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%livingMicrobeC, &
                      pt(ii)%soil(lyr)%litterCohorts(RHIZ)%CO2, &
                      pt(ii)%soil(lyr)%litterCohorts(BULK)%CO2, &
                      pt(ii)%soil(lyr)%protection_rate, &
                      pt(ii)%soil(lyr)%Qmax, &
                      pt(ii)%soil(lyr)%dissolved_carbon(LABILE), &
                      pt(ii)%soil(lyr)%dissolved_carbon(RECALCTRNT), &
                      pt(ii)%soil(lyr)%dissolved_carbon(RECALCTRNT)



----------------------------------------------------------------------------------------------------
Code updates 1/23/2017
Melannie Hartman

I am reconciling the differences between CASA, MIMIMS, and CORPSE litter inputs and CO2 fluxes.
There are complications because of the way the models treat coarse woody debri (CWD).
MIMICS and CORPSE do not have CWD litter pools like CASA does. 

In CASA: the decompostion of CWD results in direct inputs of C to the soil with CO2 loss.
MIMICS & CORPSE: we explicitly compute the composition of CWD.  The associated CO2 loss is
cwd2co2.  This flux must be added to the daily and annual Rh fluxes of the models.  The remaining 
carbon is transferred to sturctural litter in the variable cwd2str.  


----------------------------------------------------------------------------------------------------
Email from Ben Sulman 2/2/2017:

After looking at the code in corpse_inout.f90, I'm not sure the things I pointed out are the source 
of the problem because the CO2prod output is not being used to calculate the RH that the model saves 
the way I thought it was. But I might have found something else:

In corpse_cycle.f90 line 250, cwd2co2 gets added to pt(npt)%litterlayer_outputs%CO2 every time step.
However, the same pt(npt)%litterlayer_outputs is used in the subroutine save_output_line (line 270 in 
corpse_cycle.f90). Looking at that subroutine in corpse_inout.f90, I think it will overwrite 
litterlayer_outputs%CO2 with sum(pool%litterCohorts(:)%CO2) [line 482], which means you lose the 
cwd2co2 flux that was added.

In general, I think it would be safer to use the CO2prod output from update_pool to keep track of 
cumulative CO2 rather than the CO2 field of the soil cohort. The issue is that the cohort's 
cumulative CO2 is always growing, so if you run the model long enough you could eventually reach 
a point where the total cumulative CO2 value is so much bigger than the input in a particular day 
that the machine precision makes your ability to estimate CO2 flux from changes in the total value 
less accurate. For example, if your machine precision is 9 significant digits and you've been running 
for 1000 years (3.65e5 days) then you can only resolve daily CO2 flux to 4 significant digits if you're 
relying on differences in cumulative CO2. The same thing applies to the originalLitterC field, which 
will keep growing whenever carbon is added, and never decreases. This is a problem with the way I 
coded it in general, and something I probably have to fix if I want to run very long simulations and 
have the internal carbon balance check still work right. Maybe the cumulative CO2 and originalLitterC 
could be reset once a year or something. Although I guess whether this is worth worrying about would 
depend on what numerical precision you're using to compile the fortran code.

My response 2/6/2017:

Thanks!  It seems you fond the issue!  I considered several ways to correct it, including your suggestion 
for using the CO2prod variable below.  I ultimately decided to create a separate accumulator for the 
cwd2co2 flux and add it to the output when the litterlayer Rh is written to NetCDF.  I get nervous 
about changing the ways the CORPSE does C accounting for fear of introducing unintended side effects.  

Thanks also about pointing out the issue with relying on differences in the cumulative CO2 to compute 
daily CO2.  That is clearly something we need to be concerned with if outputting daily fluxes during 
a long simulation.  For now I think we are OK with the way we have been using the testbed.  When I 
run a spinup (> 1000 years), I specify annual output rather than daily output, so I am computing 
that difference once every 365 days, not every day.  I only save daily output for transient runs 
that are < 200 years in duration, and I use double precision for CORPSE.  Those CO2 accumulators 
are reset to zero whenever the model is restarted.  

----------------------------------------------------------------------------------------------------
Update to corpse_cycle.f90
Melannie Hartman
2/8/2017

CORPSE has an additional C inputs that the other models don't have: root exudation.
To keep the C balance consistent between the models we will reduce the amount
of labile root inputs to CORPSE.

  do npt=1,mp
  
      badLitter = 0
      ! Only the daily_exudate_input(LABILE) value will be reset to non-zero value. -mdh 2/8/2017
      daily_exudate_input(:) = 0
      IF(casamet%iveg2(npt) /= icewater) THEN

          ! New output variable to track daily CO2 respiration losses when cwd decomposes to structural litter. -mdh 2/6/2017
          if (idoy == 1) then
              !pt(npt)%litterlayer_outputs%cwd2co2(1..366) 
              pt(npt)%litterlayer_outputs%cwd2co2(:) = 0.0
          endif

          ! ONCE A DAY
  
          theta = casamet%moistavg(npt) ! mean volumetric soil water content (0.0 - 1.0)
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
              write(*,*) 'npt:', npt, ' daily_exudate_input(LABILE) =', daily_exudate_input(LABILE)
          else
              write(*,*) 'Reducing daily_exudate_input by: ', daily_exudate_input(LABILE) - daily_rootlitter_input(LABILE)
              write(*,*) '  npt:', npt, ' daily_exudate_input =', daily_exudate_input(LABILE)
              write(*,*) '  npt:', npt, ' daily_rootlitter_input(LABILE) =', daily_rootlitter_input(LABILE)
              daily_exudate_input(LABILE) = daily_rootlitter_input(LABILE)
              daily_rootlitter_input(LABILE) = 0.0
          endif

----------------------------------------------------------------------------------------------------
Melannie Hartman
2/27/2017

Important CORPSE theta update in corpse_cycle.csv

          ! Ben said that theta is actually fraction of water-filled pore space, not volumetric swc. -mdh 2/27/2017
          !!theta = casamet%moistavg(npt) ! mean volumetric soil water content (0.0 - 1.0)
          theta = casamet%moistavg(npt)/soil%ssat(npt) ! fraction of water-filled porespace (0.0 - 1.0)
          T = casamet%tsoilavg(npt)     ! degrees K



Write a restart file every 100 years instead of at the end of the simulation.
I did not change to frequency that the poolfluxout subroutines are called.
For CASA and MIMICS, these functions do other caculations.

Added writeToRestartCSVfile as an argument to corpse_poolfluxout (corpse_inout.f90)
SUBROUTINE corpse_poolfluxout(filename_corpseepool,mp,writeToRestartCSVfile)


casacnpdriver in casa_inout.f90:

          !! Compute annual means in casa_poolout and casa_fluxout but only write 
          !! to output CSV files every 100 years.
          if (MOD((nloop-1)*myear+iyear, 100) == 0) then
              writeToRestartCSVfile = .true.
          else
              writeToRestartCSVfile = .false.
          endif
          call casa_poolout(filename_cnpepool,iYrCnt,myear,writeToRestartCSVfile)
          call casa_fluxout(filename_cnpflux,myear,clitterinput,csoilinput,writeToRestartCSVfile)
          :
          if (isomModel == MIMICS) then
               call mimics_poolfluxout(filename_mimicsepool,mp,iYrCnt,myear,writeToRestartCSVfile)
               :
          else if (isomModel == CORPSE) then
              ! Output current year's CORPSE results for non-transient run(-mdh 5/16/2016)
              call corpse_poolfluxout(filename_corpseepool,mp,writeToRestartCSVfile)


----------------------------------------------------------------------------------------------------
Melannie Hartman
3/9/2017

Added lots of "implicit none" to corpse subroutines
Dimension theta(mp) and write that to CORPSE output instead of casamet%moistavg(npt) in subroutine corpse_soil:

    theta(npt) = casamet%moistavg(npt)/soil%ssat(npt) ! fraction of water-filled pore space (0.0 - 1.0)

    !call save_output_line(pt(npt)%litterlayer, pt(npt)%litterlayer_outputs, casamet%moistavg(npt), casamet%tsoilavg(npt))
    call save_output_line(pt(npt)%litterlayer, pt(npt)%litterlayer_outputs, theta(npt), casamet%tsoilavg(npt))

    call corpse_caccum(pt(npt)%litterlayer, pt(npt)%litterlayer_outputs, theta(npt), casamet%tsoilavg(npt), &
                       doy,casamet%ijgcm(npt))


----------------------------------------------------------------------------------------------------
Melannie Hartman
3/15/2017

The variable xfrznsoil was added to new met.nc files.  If this variable exists, read it in and use it
to compute air_filled_porosity in CORPSE.

NOTE:
This calculation in SUBROUTINE avgsoil(veg,soil,casamet) prevents vswc from being > field capacity.
    casamet%moistavg(nland)  = casamet%moistavg(nland)+ veg%froot(nland,ns) &
                               * min(soil%sfc(nland),casamet%moist(nland,ns)) 

----------------------------------------------------------------------------------------------------
Melannie Hartman
3/20/2017

If xfrzmoist is not found in the met.nc file, set the frozen moisture content to 0.0 and assume that
the soil moisture is 100% liquid.

Subroutine read_soil_carbon_namelist in corpse_soil_carbon .f90 was never actually called, so the first 
set of parameters in the namelist file (soil_carbon_nml) retained their internal hard-coded default values.  
Now, read_soil_carbon_namelist is called from casa_driver_clm.f90 right after the CORPSE_casa_nml parameters 
are read. This was never an issue before because we never changed any of those top parameters.

----------
Example: corpse_params_new.nml

&soil_carbon_nml
        vmaxref=4500e0,25e0,600e0     !Vmax at reference temperature (yr-1)
        Ea=37e3,54e3,50e3           !Activation energy (kJ/mol)
        kC=3*.01                    !Michaelis-Menton C parameter (dimensionless microbe fraction of total C)
        minMicrobeC=1e-3            !Minimum microbial biomass (fraction of total C)
        Tmic=0.25                    !Microbial turnover rate (yr-1)
        eup=.6,.05,.6                !Microbial uptake efficiency (dimensionless fraction)
        protection_rate=1.5         !Rate that carbon becomes protected (yr-1 or yr-1 kg-microbial-biomass-1 depending on microbe_driven_protection)
        microbe_driven_protection=.FALSE.   !Whether to use microbial biomass in protection rate
        protection_species=1.0,0.001,1.0     !Relative protection rate of each carbon species (between 0 and 1)
        tProtected=45.0             !Turnover time for protected carbon transition back to unprotected pool (years)
        protected_carbon_decomp_factor=0.0 !vmaxref for protected carbon is multiplied by this (0.0 if protected C is inaccessible to microbial decomposition)
        soilMaxCohorts=2            !Maximum number of cohorts in each soil carbon pool
        gas_diffusion_exp=2.5       !Exponent for gas diffusion power law dependence on theta
                                                           !See Meslin et al 2010, SSAJ
        et=0.6                     !Fraction of microbial turnover not converted to CO2
        leaching_solubility=0.0     !Rate carbon dissolves in soil water at saturated moisture (yr-1)
        DOC_deposition_rate=1.0e10  !Rate carbon is deposited from DOC (yr-1) -- currently set very high so there is no persistent DOC
        flavor_relative_solubility=1.0,0.1,1.0  !Relative solubility of each C species, between 0 and 1
        protected_relative_solubility=0.0       !Relative solubility of protected carbon, between 0 and 1
        litterDensity=22.0             !C density of litter layer (kg/m3)
                                        !22.0 roughly from Gaudinsky et al 2000


    /

&CORPSE_casa_nml
    initial_C=0.0,2.0,0.0
    exudate_npp_frac=0.02,0.0,0.0
    rhizosphere_frac=0.3

    /
----------

----------------------------------------------------------------------------------------------------
3/23/2017
Update CORPSE's Qmax function so that it is a function of clay

Email from Ben 3/22/2017:

Here is the clay function I used for CORPSE. It’s based on Melanie Mayes’ 2012 paper where she measured 
sorption isotherms for a bunch of soils. I’m using the Qmax parameter calculated from that study. 
I think in those isotherm measurements it was actually a maximum sorption amount, but I’m using it 
as a multiplier for protection_rate in each grid cell. If clay content is zero, then it sets Qmax 
to zero (there is no protection occurring at all in that case).

In this calculation, clay is clay % (out of 100) and porosity is soil porosity out of 1.0. I attached 
a plot of how the function should look across a range of clay content.

  !Qmax in mgC/kg soil from Mayes et al 2012, converted to g/m3 using solid density of 2650 kg/m3
  if(clay .le. 0) then
      Qmax= 0
  else
      Qmax = max(0.0,10**(.4833*log10(clay)+2.3282)*(1.0-porosity)*2650*1e-6)
  endif

So in LM3 I calculate this function for each grid cell using a map of clay content, and then 
multiply protection_rate by Qmax when calculating the actual protection rate in the CORPSE code. 
In the end I’m not sure the units make that much sense but it does give values in a sensible 
range for using as a multiplier.

----------------------------------------------------------------------------------------------------
4/28/2017
Melannie Hartman
I was documenting CASACNP code and questioned the fromStoCO2 calculations below, but after outputting 
the values they look fine. It was just not obvious to me why the calculations worked.

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

----------------------------------------------------------------------------------------------------
5/15/2017

I added a new option to fcasacnp_clm_testbed.lst: initcasa=3 (repeated transient) to allow a sequence 
of transient weather files to be read over and over mloop times. This can be used as an alterative
way of doing spinups.  The option initcasa=0 or initcasa=1 read a single met.nc file over and
over, and the size of this file was limited to 5-7 years.

----------------------------------------------------------------------------------------------------
11/6/2017

Started with code in /project/tss/wwieder/CASACLM/SOURCE_CODE_11.02.2017_mimTHETA2

1) Streamline netCDF tools so only necessary fields are written to met files (omit NPP, LAI)

2) Echo CASA, MIMICS, and CORPSE parameter files as they are read
   
3) Echo CASA, MIMICS, and CORPSE restart files being read upon initialization.  

----------------------------------------------------------------------------------------------------
06/02/2018

Since April 26, 2018 I've updated the point-level input/output for all three models.
If point-level output is turned on in the .lst file, then the output interval will depend
on mdaily. If mdaily=1, output is each day.  If mdaily=0, output is the state of the 
model on day 365 only.

Also updated an equation in casa_cnp.f90, SUBROUTINE casa_xrateplant. 
On 5/14/2018 Gordon discovered the a set of missing parentheses.

        !xcoldleaf(npt) = (casamet%tairk(npt)-phen%TKshed(veg%iveg(npt))-5.0)/5.0
        xcoldleaf(npt) = (casamet%tairk(npt)-(phen%TKshed(veg%iveg(npt))-5.0))/5.0

----------------------------------------------------------------------------------------------------


