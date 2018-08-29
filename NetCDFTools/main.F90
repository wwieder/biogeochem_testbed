program main

   ! Read in climate data from CLM h0 (monthly) and h1 (daily or hourly) NetCDF files 
   !
   ! Melannie Hartman
   !   Last update:
   !     October 7, 2013
   !     September 12, 2016
   !
   ! h0 files - monthly values (one file per month)
   !   each grid cell has 400+ variables
   !     each variable has one monthly value 
   !   time in these files is days since Jan 01, 1901 on the last day of the month
   !   time = 40150 in the 12th month of 2010
   !
   !   Example files:
   !     /project/tss/slevis/forMelannie/monthly/clm45sp_2deg4506_hist.clm2.h0.yyyy-mm.nc 
   !
   !     /project/tss/slevis/forMelannie/monthly/clm45sp_2deg4506_hist.clm2.h0.2000-01.nc 
   !     /project/tss/slevis/forMelannie/monthly/clm45sp_2deg4506_hist.clm2.h0.2000-02.nc 
   !         :
   !     /project/tss/slevis/forMelannie/monthly/clm45sp_2deg4506_hist.clm2.h0.2000-12.nc 
   !     /project/tss/slevis/forMelannie/monthly/clm45sp_2deg4506_hist.clm2.h0.2001-01.nc 
   !         :
   !     /project/tss/slevis/forMelannie/monthly/clm45sp_2deg4506_hist.clm2.h0.2001-12.nc 
   !
   !
   ! h1 hourly files - hourly values for one year (365 days = 8760 hours)
   !   ignore the first time slice of file 1901.
   !   the last hour in 1901 is found in the first timeslice of 1902.  Similarly for other years.
   !   time in these files is days since Jan 01, 1901 on the 2nd last hour of the last day of the year.
   !   1/24 = 0.041667  (40150/365 = 110 years)
   !   110 years * 8760 hours/year = 963,600
   !   The last 24 hours stored in the 2010 file:
   !     40149.00, 40149.04, 40149.08, 40149.12, 40149.17, 40149.21, 40149.25,
   !     40149.29, 40149.33, 40149.38, 40149.42, 40149.46, 40149.5, 40149.54,
   !     40149.58, 40149.62, 40149.67, 40149.71, 40149.75, 40149.79, 40149.83,
   !     40149.88, 40149.92, 40149.96
   !
   !   Example files:
   !     /project/tss/slevis/forMelannie/hourly/clm45sp_2deg4506_hist.clm2.h1.yyyy-mm-dd-00000.nc
   !
   !     /project/tss/slevis/forMelannie/hourly/clm45sp_2deg4506_hist.clm2.h1.2000-01-01-00000.nc
   !     /project/tss/slevis/forMelannie/hourly/clm45sp_2deg4506_hist.clm2.h1.2001-01-01-00000.nc
   !         :
   !     /project/tss/slevis/forMelannie/hourly/clm45sp_2deg4506_hist.clm2.h1.2010-01-01-00000.nc
   !
   !
   ! h1 daily files - daily values for one year (365 days)
   !   Example files:
   !     /project/bgc01/bonan/casa-cnp/clm/clm45sp_2deg4506_hist.clm2.h1.1901-01-01-00000.nc
   !         :
   !     /project/bgc01/bonan/casa-cnp/clm/clm45sp_2deg4506_hist.clm2.h1.2010-01-01-00000.nc
   !     
   !     OR
   !
   !     /project/bgc01/bonan/casa-cnp/clm/clm45sp_2deg4506_rcp85.clm2.h1.2006-01-01-00000.nc
   !         :
   !     /project/bgc01/bonan/casa-cnp/clm/clm45sp_2deg4506_rcp85.clm2.h1.2100-01-01-00000.nc
   !
  
   ! USES
   use clm_netcdf_tools
   implicit none

   ! LOCAL VARIABLES
   character*100    :: filesin           ! files.ini

!  ! Dimension lengths in NetCDF files
!  integer :: nlon
!  integer :: nlat
!  integer :: ntimes
!  ! Meteorological variables stored in NetCDF files (full time series)
!  ! All these variables are the mean over the timestep 
!  real(4), allocatable :: lat1d(:)         ! lat(lat) latitude (degrees_north)
!  real(4), allocatable :: lon1d(:)         ! lon(lon) longitude (degrees_east)
!  real(4), allocatable :: time(:)          ! days since 1901-01-01 00:00:00
!  real(4), allocatable :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
!  real(4), allocatable :: fpsn(:,:,:)      ! fpsn(lon, lat, time) photosyntheis (umol/m2/sec)
!  real(4), allocatable :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
!  real(4), allocatable :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
!  real(4), allocatable :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
!  real(4), allocatable :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation (W/m^2)
!  real(4), allocatable :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
!  real(4), allocatable :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
!  real(4), allocatable :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
!  real(4), allocatable :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
!  real(4), allocatable :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
!  real(4), allocatable :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
!  real(4), allocatable :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
!  real(4), allocatable :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
!
!  Some utilities that are available:
!  call printNCfileInfo('/project/tss/slevis/forMelannie/monthly/clm45sp_2deg4506_hist.clm2.h0.2000-01.nc')

   filesin = 'files.ini'

   !Read and process monthly, daily, or hourly climate files one year at a time
   !This function also runs the aggregated canopy model each year when hourly climate is read.

   call readClimateSequence(filesin)

   write(*,*) 'Simulation Complete.'

end program main

