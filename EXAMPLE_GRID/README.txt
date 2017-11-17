

The version of casaclm_corpse_mimics in this directory reads N deposition from met.nc files now instead of 
from the gridinfo_igbpz.csv file.  

Melannie Hartman
9/24/2016

====================================================================================================
Run spinups and 20th century transient runs for all three models
Melannie Hartman
October 17-24, 2016

----------------------------------------------------------------------------------------------------
CASACNP (prespin, spin, transient)

1) Run the prespin to initialize vegetation pools for all models.  This prespin is actaully only needed for
   MIMICS and CORPSE, but I'm doing it for CASACNP also for consistency.

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_prespin_1901_x100_STEP_1of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &
 
2) Run the CASACNP accelerated spinup, which accelerated decomposition of passive pool 10x
   * Note, the restart file has to be manually updated, with passive SOC multiplied 10x for step 3

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_casacnp_ADspin_1901_1920_x499_STEP_2of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &

   * Concatenate the annual spinup files to be used with NCL scripts later to view time series data.

     cd GRID/OUTPUT_GRID_CRU_NCEP/CASACNP_ADSPIN_test3
     ncrcat casaclm_pool_flux_9???.nc casaclm_pool_flux_9000_9980.nc &
     ls -lrt *_pool_flux_????_????.nc
     rm *_pool_flux_????.nc

   * Run NCL scripts to look at spinup results

     cd /GRID/NCL_GRID/CASACNP_SPIN
     ncl globalCTimeseries.ncl

   * Check that SOM pools are spunup, 
     1) dSOC < 0.01 Pg C globally,
     2) > 98% of grid cells changed less than 1 g C m-2
     3) > 98% of grid cells changed less than 0.1%  

   * these scripts also provide a visual look at results, but have not been updated in some time
     ncl casaclm_byVegCat.ncl
     ncl casaclm_map.ncl
     ncl globalGPPTimeseries.ncl

3) Run the CASACNP post accelerated spinup
   returns to default parameterization
   * Note, the restart file has to be manually updated, with passive SOC (column u) multiplied 10x
   * I did this in excel, but save the new files a a wondow, formatted .csv

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_casacnp_postADspin_1901_1920_x200_STEP_3of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &

   * Concatenate the annual spinup files to be used with NCL scripts later to view time series data.

     cd GRID/OUTPUT_GRID_CRU_NCEP/CASACNP_SPIN_test3
     ncrcat casaclm_pool_flux_3???.nc casaclm_pool_flux_3000_3999.nc &
     ls -lrt *_pool_flux_????_????.nc
     rm *_pool_flux_????.nc

   * Run NCL scripts to look at spinup results

     cd GRID/NCL_GRID/CASACNP_SPIN
     ncl globalCTimeseries.ncl

4) Run the CASACNP transient simulation 1901-2010

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_casacnp_trans_1901_2010_STEP_4of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &

   * Concatenate the netcdf files for time series analysis.
     cacluclate annual averages of stocks & fluxes for times series analyses

     cd GRID/NCL_GRID/
     ncl daily_to_annual_3.ncl

   * instructions below are depreciated, but retained  
     Note that daily files are too large to concatenate for all 110 years.
     This example shows how to look at a subset of years.

     cd GRID/OUTPUT_GRID/CASACNP_TRANS
     rm *yyyy*.nc
     ncrcat casaclm_pool_flux_2???_daily.nc casaclm_pool_flux_2000_2010_daily.nc

   * Run NCL scripts to look at transient results

     cd GRID/NCL_GRID/CASACNP_TRANS
     ncl globalCTimeseries_daily.ncl
     ncl globalGPPTimeseries_daily.ncl


----------------------------------------------------------------------------------------------------
MIMICS (prespin, spin, transient)

ncrcat mimics_pool_flux_????.nc mimics_pool_flux_0001_2100.nc

1) Run the prespin to initialize vegetation pools for all models.  This same prespin is 
   used for all three SOM models.

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_prespin_1901_x100_STEP_1of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &

2) Run the MIMICS spinup

   cd GRID/EXAMPLE_GRID
   cp Rfcsacnp_clm_mimics_spin_1901_1920_x499_STEP_2of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &
   
   * Concatenate the netcdf files for time series analysis.

     cd GRID/OUTPUT_GRID/MIMICS_SPIN_test3
     ncrcat casaclm_pool_flux_9???.nc casaclm_pool_flux_9000_9980.nc &
     ncrcat mimics_pool_flux_9???.nc mimics_pool_flux_9000_9980.nc &
     ls -lrt *_pool_flux_????_????.nc
     rm *_pool_flux_????.nc

   * Run NCL scripts to look at spinup results

     cd /GRID/NCL_GRID/CASACNP_SPIN
     ncl globalCTimeseries.ncl

   * Check that SOM pools are spunup,
     1) dSOC < 0.01 Pg C globally,
     2) > 98% of grid cells changed less than 1 g C m-2
     3) > 98% of grid cells changed less than 0.1%

   * Users may need to run multiple spinup simulations, pointing to old restart files for all models

3) Run the MIMICS transient simulation 1901-2010

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_mimics_trans_1901_2010_STEP_4of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &

   * Concatenate the netcdf files for time series analysis as with CASA, above.

----------------------------------------------------------------------------------------------------
CORPSE (prespin, spin, transient)

1) Run the prespin to initialize vegetation pools for all models.  This same prespin is 
   used for all three SOM models.

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_prespin_1901_x100_STEP_1of4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &

  
2) Run the CORPSE spinup

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_corpse_spin_1901_1920_x499_STEP_2aof4.lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &

   * Concatenate the netcdf files for time series analysis.

     cd GRID/OUTPUT_GRID/CORPSE_SPIN_test3a
     ncrcat casaclm_pool_flux_9???.nc casaclm_pool_flux_9000_9980.nc &
     ncrcat corpse_pool_flux_9???.nc corpse_pool_flux_9000_9980.nc

     ls -lrt c*_pool_flux_????_????.nc
     rm c*_pool_flux_????.nc

   * Run NCL scripts to look at spinup results

     cd GRID/NCL_GRID/CORPSE_SPIN
     ncl globalCTimeseries.ncl

   * Check that SOM pools are spunup,
     1) dSOC < 0.01 Pg C globally,
     2) > 98% of grid cells changed less than 1 g C m-2
     3) > 98% of grid cells changed less than 0.1%
     
   * Users may need to run multiple spinup simulations, pointing to old restart files for all models
   
   3) Run the MIMICS transient simulation 1901-2010

   cd GRID/EXAMPLE_GRID
   cp fcasacnp_clm_corpse_trans_1901_2010_STEP_4of4lst fcasacnp_clm_testbed.lst
   ./casaclm_mimics_corpse &
   
   * Concatenate the netcdf files for time series analysis as with CASA, above.
