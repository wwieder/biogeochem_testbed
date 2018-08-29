module clm_netcdf_tools

!-----------------------------------------------------------------------
!BOP
!
!  !MODULE: clm_netcdf_tools
!
!  !DESCRIPTION:
!
!  !USES:
!   use clm_common
!   use canopy_model
    use casacnp_model
!
!  !PUBLIC TYPES:
   implicit none

   !Parameters are defined in clm_common

! !PUBLIC MEMBER FUNCTIONS:
   public :: allocateMonthlyVariables
!  public :: allocateHourlyVariables
   public :: allocateDailyHourlyVariables !! Updated 8/11/2014
   public :: readMonthlyHistoryFile
   public :: readDailyHistoryFile         !! Added 8/11/2014
!  public :: readHourlyHistoryFile
   public :: readClimateSequence
   public :: readClimateSequenceMerge
   public :: readInitFile
   public :: readNdepositionFile
   public :: readSurfaceData
   public :: soilWaterPotential
   public :: writeClimateToCSVFile

! !PUBLIC DATA MEMBERS:
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !PRIVATE DATA MEMBERS:
!
! !I/O text files
!  unit=10      Input (required)  : Initialization file (files.ini)
!  unit=14      Input (required)  : transient CO2 concentrations (year  CO2ppm)
!  unit=11      Output (optional) : hourly climate text file for LATINDX, LONINDX
!  unit=12      Output (optional) : daily climate text file for LATINDX, LONINDX
!  unit=16      Output (optional) : daily gpp text file for LATINDX, LONINDX

! !REVISION HISTORY:
!  Created by Melannie Hartman
!  melannie@ucar.edu
!  July 2013-March 2014
!  Remove NPP and LAI from met.nc files -mdh 11/6/2017
!
!  NOTES:
!    
!  NetCDF float are real(4), double is real(8). If types are not correct SIGSEGV occurs:
!    Program received signal 11 (SIGSEGV): Segmentation fault.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readSurfaceData
!
! !INTERFACE:
   subroutine readSurfaceData(ncsrfcfilename, lsmlon, lsmlat, lsmpft, nlevsoi, &
                              pctpft, pctsand, pctclay, orgdens, watsat, watfc)
! WW added nlevgrnd, watsat & watfc

! !DESCRIPTION:
!  Read PCT_PFT, PCT_SAND, PCT_CLAY, and ORGANIC from surface data set.
!  WW also WATFC, WATSAT, copied from h0 files & renamed
!  Return those arrays and their dimensions. 
!  Note: Code to read LATIXY and LONGXY is available but commented out to save execution time.
!
! !USES:
!  none
!
! !ARGUMENTS:
   character(len=*), intent(in) :: ncsrfcfilename       ! netcdf filename for surface data set
   integer, intent(out) :: lsmlon                       ! length of longitude dimensions 
   integer, intent(out) :: lsmlat                       ! length of latitude dimension 
   integer, intent(out) :: lsmpft                       ! length of pft dimension 
   integer, intent(out) :: nlevsoi                      ! number of soil levels

   !PCT_PFT, PCT_CLAY and PCT_CLAY are defined as DOUBLE in surface data netcdf file. Use real(8).
   real(8), allocatable, intent(out) :: pctpft(:,:,:)   ! percent plant functional type of gridcell
   real(8), allocatable, intent(out) :: pctsand(:,:,:)  ! percent sand (lon,lat,nlevsoi) at soil levels 
   real(8), allocatable, intent(out) :: pctclay(:,:,:)  ! percent clay (lon,lat,nlevsoi) at soil levels (lon,lat,nlevsoi)
   real(8), allocatable, intent(out) :: orgdens(:,:,:)  ! organic matter density (lon,lat,nlevsoi) at soil levels (kg/m3) 
   real(8), allocatable, intent(out) :: watsat(:,:,:)   ! saturated water content(lon,lat,levgrnd) at soil levels 
   real(8), allocatable, intent(out) :: watfc(:,:,:)    ! water field capacity(lon,lat,levgrnd) at soil levels 
!! WW added, WATSAT and WATFC both floats

! !CALLED FROM:
!  subroutine readHourlyHistoryFile
!  subroutine readDailyHistoryFile
!
! !REVISION HISTORY:
!  Created by Melannie Hartman
!
!EOP

! !LOCAL VARIABLES:
      integer :: ncid                           ! netcdf file ID
      integer :: status                         ! function return status
      integer :: lat_dimid                      ! netcdf dimension id
      integer :: lon_dimid                      ! netcdf dimension id
      integer :: pft_dimid                      ! netcdf dimension id
      integer :: nlevsoi_dimid                  ! netcdf dimension id
!     integer :: latixy_varid                   ! netcdf variable id
!     integer :: longxy_varid                   ! netcdf variable id
!     integer :: start2(2), count2(2)           ! start and count arrays for reading 2-D data from netcdf files
!     !LONGXY and LATIXY are defined as DOUBLE in netcdf file. Use real(8).
!     real(8), allocatable :: longxy(:,:)       ! longitude in 2-dimensions (lsmlon,lsmlat)
!     real(8), allocatable :: latixy(:,:)       ! latitude in 2-dimensions (lsmlon,lsmlat)
      integer :: pctpft_varid                   ! netcdf variable id
      integer :: pctsand_varid                  ! netcdf variable id
      integer :: pctclay_varid                  ! netcdf variable id
      integer :: orgdens_varid                  ! netcdf variable id
! WW added
      integer :: watsat_varid                   ! netcdf variable id
      integer :: watfc_varid                    ! netcdf variable id
      integer :: start3(3), count3(3),count4(3)  ! start and count arrays for reading 3-D data from netcdf files

!-----------------------------------------------------------------------
      if (verbose .ge. 0) print *, "Reading surface data from file ", trim(ncsrfcfilename), "..."
     
      ! Get dimension ids
      ! There may be more dimensions in the file - if so use the examples below to read them in

      status = nf_open(ncsrfcfilename, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, "")
         
      status = nf_inq_dimid(ncid, "lsmlat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lsmlat")

      status = nf_inq_dimid(ncid, "lsmlon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lsmlon")

      status = nf_inq_dimid(ncid, "lsmpft", pft_dimid)
      if (status /= nf_noerr) call handle_err(status, "lsmpft")

      status = nf_inq_dimid(ncid, "nlevsoi", nlevsoi_dimid)
      if (status /= nf_noerr) call handle_err(status, "nlevsoi")

      ! Get dimension sizes

      status = nf_inq_dimlen(ncid, lat_dimid, lsmlat)
      if (status /= nf_noerr) call handle_err(status, "lsmlat")

      status = nf_inq_dimlen(ncid, lon_dimid, lsmlon)
      if (status /= nf_noerr) call handle_err(status, "lsmlon")

      status = nf_inq_dimlen(ncid, pft_dimid, lsmpft)
      if (status /= nf_noerr) call handle_err(status, "lsmpft")

      status = nf_inq_dimlen(ncid, nlevsoi_dimid, nlevsoi)
      if (status /= nf_noerr) call handle_err(status, "nlevsoi")

      if (verbose .ge. 1) then
         print *, "  readSurfaceData: lsmlat: ", lsmlat
         print *, "  readSurfaceData: lsmlon: ", lsmlon
         print *, "  readSurfaceData: lsmpft: ", lsmpft
         print *, "  readSurfaceData: nlevsoi: ", nlevsoi
      endif


      ! Get variable ids

!     status = nf_inq_varid(ncid, "LATIXY", latixy_varid)
!     if (status /= nf_noerr) call handle_err(status, "LATIXY")
!
!     status = nf_inq_varid(ncid, "LONGXY", longxy_varid)
!     if (status /= nf_noerr) call handle_err(status, "LONGXY")

      status = nf_inq_varid(ncid, "PCT_PFT", pctpft_varid)
      if (status /= nf_noerr) call handle_err(status, "PCT_PFT")

      status = nf_inq_varid(ncid, "PCT_SAND", pctsand_varid)
      if (status /= nf_noerr) call handle_err(status, "PCT_SAND")

      status = nf_inq_varid(ncid, "PCT_CLAY", pctclay_varid)
      if (status /= nf_noerr) call handle_err(status, "PCT_CLAY")

      status = nf_inq_varid(ncid, "ORGANIC", orgdens_varid)
      if (status /= nf_noerr) call handle_err(status, "ORGANIC")

! WW added
      status = nf_inq_varid(ncid, "WATSAT", watsat_varid)
      if (status /= nf_noerr) call handle_err(status, "WATSAT")

      status = nf_inq_varid(ncid, "WATFC", watfc_varid)
      if (status /= nf_noerr) call handle_err(status, "WATFC" )
! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: lsmlat, lsmlon (2-D) 
!     allocate(latixy(1:lsmlon, 1:lsmlat))
!     allocate(longxy(1:lsmlon, 1:lsmlat))
!
!     latixy(:,:) = 0.0
!     longxy(:,:) = 0.0
!
!     start2 = (/ 1, 1 /)
!     count2 = (/ lsmlat, lsmlon /)
!
!     status = nf_get_var(ncid, latixy_varid, latixy, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "LATIXY")
!
!     status = nf_get_var(ncid, longxy_varid, longxy, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "LONGXY")

      
! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: nlevsoi or lsmpft, lsmlat, lsmlon (3-D) 
      allocate(pctpft(1:lsmlon, 1:lsmlat, 1:lsmpft))
      allocate(pctsand(1:lsmlon, 1:lsmlat, 1:nlevsoi))
      allocate(pctclay(1:lsmlon, 1:lsmlat, 1:nlevsoi))
      allocate(orgdens(1:lsmlon, 1:lsmlat, 1:nlevsoi))
! WW added
      allocate(watsat(1:lsmlon, 1:lsmlat, 1:nlevsoi))
      allocate(watfc(1:lsmlon, 1:lsmlat, 1:nlevsoi))

      pctpft(:,:,:) = 0.0
      pctsand(:,:,:) = 0.0
      pctclay(:,:,:) = 0.0
      orgdens(:,:,:) = 0.0
! WW added
      watsat(:,:,:) = 0.0
      watfc(:,:,:) = 0.0

      start3 = (/ 1, 1, 1 /)
      count3 = (/ lsmpft, lsmlat, lsmlon /)

      status = nf_get_var(ncid, pctpft_varid, pctpft, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "get PCT_PFT")


      start3 = (/ 1, 1, 1 /)
      count3 = (/ nlevsoi, lsmlat, lsmlon /)

      status = nf_get_var(ncid, pctsand_varid, pctsand, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "get PCT_SAND")

      status = nf_get_var(ncid, pctclay_varid, pctclay, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "get PCT_CLAY")

      status = nf_get_var(ncid, orgdens_varid, orgdens, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "get ORGANIC")

! WW added 
     status = nf_get_var(ncid, watsat_varid, watsat, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "get WATSAT")

     status = nf_get_var(ncid, watfc_varid, watfc, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "get WATFC")

!     Echo values for testing...
!     do i = 1, lsmlat
!        do j = 1, lsmlon
!           print *, j, ",", i, ": ", latixy(j,i), longxy(j,i), pctsand(j,i,1), pctclay(j,i,1), orgdens(j,i,1)
!        end do
!     end do

      status = nf_close(ncid)

      if (verbose .ge. 0) print *, "Done reading surface data from ", trim(ncsrfcfilename)

   end subroutine readSurfaceData

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateMonthlyVariables
!
! !INTERFACE:
   subroutine allocateMonthlyVariables(ncfilename, nFiles, nlon, nlat, ntime, lon1d, lat1d, time, &
                                       flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                                       rain, snow, tbot, thbot, wind)
!
! !DESCRIPTION:
!  Allocate memory for monthly 3-D meterological variables for multiple time periods (nFiles * ntime)
!
! !USES:
!  none
!
!  !ARGUMENTS:
    character(len=*), intent(in)    :: ncfilename        ! netcdf filename
    integer, intent(in)             :: nFiles            ! number of NetCDF files to read
    integer, intent(out) :: nlon                         ! number of longitudes 
    integer, intent(out) :: nlat                         ! number of latitudes
    integer, intent(out) :: ntime                        ! number of monthly time periods per file 
    real(4), allocatable :: lon1d(:)                     ! lon(lon) (degrees_east)
    real(4), allocatable :: lat1d(:)                     ! lat(lat) (degrees_north)
    real(4), allocatable :: time(:)                       ! days since 1901-01-01 00:00:00
    ! All these variables are the mean over the timestep except where indicated as "total for month"
    real(4), allocatable, intent(out) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
    real(4), allocatable, intent(out) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
    real(4), allocatable, intent(out) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
    real(4), allocatable, intent(out) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
    real(4), allocatable, intent(out) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
    real(4), allocatable, intent(out) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
    real(4), allocatable, intent(out) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
    real(4), allocatable, intent(out) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
    real(4), allocatable, intent(out) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
    real(4), allocatable, intent(out) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
    real(4), allocatable, intent(out) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
    real(4), allocatable, intent(out) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
    real(4), allocatable, intent(out) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)

! !CALLED FROM:
!  subroutine ??
!
! !REVISION HISTORY:
!  Created by Melannie Hartman
!
! !LOCAL VARIABLES:
!EOP

      !integer :: lat, lon                      ! length of latitude and longitude dimensions (global grid sizes)
      integer :: ncid                           ! netcdf file ID
      integer :: status                         ! function return status
      integer :: lat_dimid                      ! netcdf dimension id
      integer :: lon_dimid                      ! netcdf dimension id
      integer :: time_dimid                     ! netcdf dimension id
      integer :: lat_varid                      ! netcdf variable id
      integer :: lon_varid                      ! netcdf variable id
      integer :: time_varid                     ! netcdf variable id
      !integer :: start3(3), count3(3)           ! start and count arrays for reading data from netcdf files
      !real(4), allocatable :: var3d(:,:,:)      ! variable with 3-dimensions (lon,lat,time)
!-----------------------------------------------------------------------
      ! Determine dimensions and if grid file is 2d or 1d

      status = nf_open(ncfilename, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, ncfilename)
         
      status = nf_inq_dimid(ncid, "time", time_dimid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_dimid(ncid, "lat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimid(ncid, "lon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_dimlen(ncid, time_dimid, ntime)
      if (status /= nf_noerr) call handle_err(status, "ntime")

      status = nf_inq_dimlen(ncid, lat_dimid, nlat)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimlen(ncid, lon_dimid, nlon)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_varid(ncid, "lat", lat_varid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_varid(ncid, "lon", lon_varid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_varid(ncid, "time", time_varid)
      if (status /= nf_noerr) call handle_err(status, "time")

      if (verbose .ge. 1) then
         print *, "  allocateMonthlyVariables: length(lat): ", nlat
         print *, "  allocateMonthlyVariables: length(lon): ", nlon
         print *, "  allocateMonthlyVariables: length(time): ", ntime
      endif

      allocate(lon1d(1:nlon))
      allocate(lat1d(1:nlat))
      allocate(time(1:nFiles*ntime))
! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: time, lat, lon
      allocate(flds(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsds(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsnd(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsni(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsvd(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsvi(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(pbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(qbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(rain(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(snow(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(tbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(thbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(wind(1:nlon, 1:nlat, 1:nFiles*ntime))

      lon1d(:) = 0.0
      lat1d(:) = 0.0
      time(:) = 0.0

      flds(:,:,:) = 0.0
      fsds(:,:,:) = 0.0
      fsdsnd(:,:,:) = 0.0
      fsdsni(:,:,:) = 0.0
      fsdsvd(:,:,:) = 0.0
      fsdsvi(:,:,:) = 0.0
      pbot(:,:,:) = 0.0
      qbot(:,:,:) = 0.0
      rain(:,:,:) = 0.0
      snow(:,:,:) = 0.0
      tbot(:,:,:) = 0.0
      thbot(:,:,:) = 0.0
      wind(:,:,:) = 0.0

      status = nf_close(ncid)

   end subroutine allocateMonthlyVariables


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readMonthlyHistoryFile
!
! !INTERFACE:
   subroutine readMonthlyHistoryFile(ncfilename, nlon, nlat, ntime, lon1d, lat1d, time, &
                                     flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                                     rain, snow, tbot, thbot, wind)
!
! !DESCRIPTION:
!  Read longitude, latitude, and meteorological variables from monthly history file.
!
! !USES:
!  none
!
!  !ARGUMENTS:
   character(len=*), intent(in)    :: ncfilename  ! netcdf filename
   integer, intent(out) :: nlon                          ! number of longitudes 
   integer, intent(out) :: nlat                          ! number of latitudes
   integer, intent(out) :: ntime                         ! number of monthly time periods per file 
   real(4), allocatable, intent(out) :: lon1d(:)         ! lon(lon) (degrees_east)
   real(4), allocatable, intent(out) :: lat1d(:)         ! lat(lat) (degrees_north)
   real(4), allocatable, intent(out) :: time(:)          ! days since 1901-01-01 00:00:00
   ! All these variables are the mean over the timestep except where indicated as "total for month"
   real(4), allocatable, intent(out) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), allocatable, intent(out) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), allocatable, intent(out) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), allocatable, intent(out) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), allocatable, intent(out) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), allocatable, intent(out) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), allocatable, intent(out) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)

!
! !CALLED FROM:
! subroutine ??
!
! !REVISION HISTORY:
! Created by Melannie Hartman
!
!
! !LOCAL VARIABLES:
!EOP

      !integer :: lat, lon                      ! length of latitude and longitude dimensions (global grid sizes)
      integer :: ncid                           ! netcdf file ID
      integer :: status                         ! function return status
      integer :: lat_dimid                      ! netcdf dimension id
      integer :: lon_dimid                      ! netcdf dimension id
      integer :: time_dimid                     ! netcdf dimension id
      integer :: lat_varid                      ! netcdf variable id
      integer :: lon_varid                      ! netcdf variable id
      integer :: time_varid                     ! netcdf variable id
      integer :: varid                          ! netcdf variable id
      integer :: start3(3), count3(3)           ! start and count arrays for reading data from netcdf files
      !real(4), allocatable :: var3d(:,:,:)     ! variable with 3-dimensions (lon,lat,time)
!-----------------------------------------------------------------------

      if (verbose .ge. 0) print *, "Reading monthly climate data from ", trim(ncfilename), "..."

      status = nf_open(ncfilename, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, ncfilename)
         
      status = nf_inq_dimid(ncid, "time", time_dimid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_dimid(ncid, "lat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimid(ncid, "lon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_dimlen(ncid, time_dimid, ntime)
      if (status /= nf_noerr) call handle_err(status, "ntime")

      status = nf_inq_dimlen(ncid, lat_dimid, nlat)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimlen(ncid, lon_dimid, nlon)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_varid(ncid, "lat", lat_varid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_varid(ncid, "lon", lon_varid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_varid(ncid, "time", time_varid)
      if (status /= nf_noerr) call handle_err(status, "time")

      if (verbose .ge. 1) then
         print *, "  readMonthlyHistoryFile: length(lat): ", nlat
         print *, "  readMonthlyHistoryFile: length(lon): ", nlon
         print *, "  readMonthlyHistoryFile: length(time): ", ntime
      endif

      allocate(lon1d(1:nlon))
      allocate(lat1d(1:nlat))
      allocate(time(1:ntime))
! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: time, lat, lon
      allocate(flds(1:nlon, 1:nlat, 1:ntime))
      allocate(fsds(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsnd(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsni(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsvd(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsvi(1:nlon, 1:nlat, 1:ntime))
      allocate(pbot(1:nlon, 1:nlat, 1:ntime))
      allocate(qbot(1:nlon, 1:nlat, 1:ntime))
      allocate(rain(1:nlon, 1:nlat, 1:ntime))
      allocate(snow(1:nlon, 1:nlat, 1:ntime))
      allocate(tbot(1:nlon, 1:nlat, 1:ntime))
      allocate(thbot(1:nlon, 1:nlat, 1:ntime))
      allocate(wind(1:nlon, 1:nlat, 1:ntime))
 
      lon1d(:) = 0.0
      lat1d(:) = 0.0
      time(:) = 0.0
      flds(:,:,:) = 0.0
      fsds(:,:,:) = 0.0
      fsdsnd(:,:,:) = 0.0
      fsdsni(:,:,:) = 0.0
      fsdsvd(:,:,:) = 0.0
      fsdsvi(:,:,:) = 0.0
      pbot(:,:,:) = 0.0
      qbot(:,:,:) = 0.0
      rain(:,:,:) = 0.0
      snow(:,:,:) = 0.0
      tbot(:,:,:) = 0.0
      thbot(:,:,:) = 0.0
      wind(:,:,:) = 0.0

      status = nf_get_var(ncid, lat_varid, lat1d)
      if (status /= nf_noerr) call handle_err(status, "LAT1D")

      status = nf_get_var(ncid, lon_varid, lon1d)
      if (status /= nf_noerr) call handle_err(status, "LON1D")

      status = nf_get_var(ncid, time_varid, time)
      if (status /= nf_noerr) call handle_err(status, "time")


      start3 = (/ 1, 1 ,1 /)
      count3 = (/ ntime, nlat, nlon /)

      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
 
      if (verbose .ge. 2) print *, "FLDS"
      status = nf_inq_varid(ncid, "FLDS", varid)
      if (status /= nf_noerr) call handle_err(status, "FLDS")

      status = nf_get_var(ncid, varid, flds, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "flds")


      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDS"
      status = nf_inq_varid(ncid, "FSDS", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDS")

      status = nf_get_var(ncid, varid, fsds, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsds")


      ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDSND"
      status = nf_inq_varid(ncid, "FSDSND", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSND")

      status = nf_get_var(ncid, varid, fsdsnd, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsnd")


      ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDSNI"
      status = nf_inq_varid(ncid, "FSDSNI", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSNI")

      status = nf_get_var(ncid, varid, fsdsni, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsni")


      ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)

      if (verbose .ge. 2) print *, "FSDSVD"
      status = nf_inq_varid(ncid, "FSDSVD", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSVD")

      status = nf_get_var(ncid, varid, fsdsvd, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsvd")


      ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDSVI"
      status = nf_inq_varid(ncid, "FSDSVI", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSVI")

      status = nf_get_var(ncid, varid, fsdsvi, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsvi")


      ! pbot(lon, lat, time) atmospheric pressure (Pa)

      if (verbose .ge. 2) print *, "PBOT"
      status = nf_inq_varid(ncid, "PBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "PBOT")

      status = nf_get_var(ncid, varid, pbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "pbot")


      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)

      if (verbose .ge. 2) print *, "QBOT"
      status = nf_inq_varid(ncid, "QBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "QBOT")

      status = nf_get_var(ncid, varid, qbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "qbot")


      ! rain(lon, lat, time) atmospheric rain (mm/sec)

      if (verbose .ge. 2) print *, "RAIN"
      status = nf_inq_varid(ncid, "RAIN", varid)
      if (status /= nf_noerr) call handle_err(status, "RAIN")

      status = nf_get_var(ncid, varid, rain, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "rain")


      ! snow(lon, lat, time) atmospheric snow (mm/sec)

      if (verbose .ge. 2) print *, "SNOW"
      status = nf_inq_varid(ncid, "SNOW", varid)
      if (status /= nf_noerr) call handle_err(status, "SNOW")

      status = nf_get_var(ncid, varid, snow, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "snow")


      ! tbot(lon, lat, time) atmospheric air temperature (K)

      if (verbose .ge. 2) print *, "TBOT"
      status = nf_inq_varid(ncid, "TBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "TBOT")

      status = nf_get_var(ncid, varid, tbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "tbot")


      ! thbot(lon, lat, time) atmospheric air potential temperature (K)

      if (verbose .ge. 2) print *, "THBOT"
      status = nf_inq_varid(ncid, "THBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "THBOT")

      status = nf_get_var(ncid, varid, thbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "thbot")


      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)

      if (verbose .ge. 2) print *, "WIND"
      status = nf_inq_varid(ncid, "WIND", varid)
      if (status /= nf_noerr) call handle_err(status, "WIND")

      status = nf_get_var(ncid, varid, wind, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "wind")

      status = nf_close(ncid)

      if (verbose .ge. 0) print *, "Done reading monthly climate data from ", trim(ncfilename)

   end subroutine readMonthlyHistoryFile

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readNdepositionFile
!
! !INTERFACE:
   subroutine readNdepositionFile(ncfilename, nlon, nlat, ntime, lon1d, lat1d, time, calyear, ndep)
!
! !DESCRIPTION:
!  Read annual N deposition from netcdf file.  
!  The file I was testing with contained annual values for years 1849-2006.
!
! !USES:
!  none
!
!  !ARGUMENTS:
   character(len=*), intent(in)    :: ncfilename         ! netcdf filename
   integer, intent(out) :: nlon                          ! number of longitudes 
   integer, intent(out) :: nlat                          ! number of latitudes
   integer, intent(out) :: ntime                         ! number of years of data
   ! lon1d, lat1d, and time must be DOUBLE, but ndep is real(4)
   real(8), allocatable, intent(out) :: lon1d(:)         ! lon(lon) (degrees_east)
   real(8), allocatable, intent(out) :: lat1d(:)         ! lat(lat) (degrees_north)
   real(8), allocatable, intent(out) :: time(:)          ! days since 0000-01-01 00:00
   integer, allocatable, intent(out) :: calyear(:)       ! calendar year for each time
   real(4), allocatable, intent(out) :: ndep(:,:,:)      ! ndep(lon, lat, time) Sum of NOy and NHx deposition (gN/m2/yr)

!
! !CALLED FROM:
! subroutine ??
!
! !REVISION HISTORY:
! Created by Melannie Hartman
!
!
! !LOCAL VARIABLES:
!EOP

      integer :: ncid                           ! netcdf file ID
      integer :: status                         ! function return status
      integer :: lat_dimid                      ! netcdf dimension id
      integer :: lon_dimid                      ! netcdf dimension id
      integer :: time_dimid                     ! netcdf dimension id
      integer :: lat_varid                      ! netcdf variable id
      integer :: lon_varid                      ! netcdf variable id
      integer :: time_varid                     ! netcdf variable id
      integer :: year_varid                     ! netcdf variable id
      integer :: ndep_varid                     ! netcdf variable id
      integer :: start3(3), count3(3)           ! start and count arrays for reading data from netcdf files
!-----------------------------------------------------------------------

      if (verbose .ge. 0) print *, "Reading N deposition file ", trim(ncfilename), "..."

      status = nf_open(ncfilename, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, ncfilename)
         
      status = nf_inq_dimid(ncid, "time", time_dimid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_dimid(ncid, "lat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimid(ncid, "lon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_dimlen(ncid, time_dimid, ntime)
      if (status /= nf_noerr) call handle_err(status, "ntime")

      status = nf_inq_dimlen(ncid, lat_dimid, nlat)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimlen(ncid, lon_dimid, nlon)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_varid(ncid, "lat", lat_varid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_varid(ncid, "lon", lon_varid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_varid(ncid, "time", time_varid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_varid(ncid, "YEAR", year_varid)
      if (status /= nf_noerr) call handle_err(status, "YEAR")

      if (verbose .ge. 1) then
         print *, "  fndep length(lat): ", nlat
         print *, "  fndep length(lon): ", nlon
         print *, "  fndep length(time): ", ntime
      endif

      allocate(lon1d(1:nlon))
      allocate(lat1d(1:nlat))
      allocate(time(1:ntime))
      allocate(calyear(1:ntime))
! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: time, lat, lon
      allocate(ndep(1:nlon, 1:nlat, 1:ntime))

      lon1d(:) = 0.0
      lat1d(:) = 0.0
      time(:) = 0.0
      calyear(:) = 0
      ndep(:,:,:) = 0.0

      status = nf_get_var(ncid, lat_varid, lat1d)
      if (status /= nf_noerr) call handle_err(status, "nf_get_var(LAT1D)")

      status = nf_get_var(ncid, lon_varid, lon1d)
      if (status /= nf_noerr) call handle_err(status, "nf_get_var(LON1D)")

      status = nf_get_var(ncid, time_varid, time)
      if (status /= nf_noerr) call handle_err(status, "nf_get_var(time)")

      status = nf_inq_varid(ncid, "NDEP_year", ndep_varid)
      if (status /= nf_noerr) call handle_err(status, "NDEP_year")


      !Retrieve variable values 

      !! Note, calyear could be read directly from variable "YEAR"
      !! calyear(:) = time(:) / DAYS_PER_YEAR
      status = nf_get_var(ncid, year_varid, calyear)
      if (status /= nf_noerr) call handle_err(status, "nf_get_var(YEAR)")

      start3 = (/ 1, 1 ,1 /)
      count3 = (/ ntime, nlat, nlon /)

      ! ndep(lon, lat, time) atmospheric N deposition

      status = nf_get_var(ncid, ndep_varid, ndep, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "nf_get_var(ndep)")

      status = nf_close(ncid)

      if (verbose .ge. 0) print *, "Done reading N deposition from file ", trim(ncfilename)

   end subroutine readNdepositionFile

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateDailyHourlyVariables
!
! !INTERFACE:
   subroutine allocateDailyHourlyVariables(ncfilename, nFiles, nlon, nlat, ntime, lon1d, lat1d, time, &
                                           flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                                           rain, snow, tbot, thbot, wind)
!
! !DESCRIPTION:
!  Allocate memory for daily 3-D meterological variables for multiple time periods (nFiles * ntime)
!
! !USES:
!  none
!
!     !ARGUMENTS:
      character(len=*), intent(in)    :: ncfilename      ! netcdf filename
      integer, intent(in)             :: nFiles          ! number of NetCDF files that will be read
   integer, intent(out) :: nlon                          ! number of longitudes 
   integer, intent(out) :: nlat                          ! number of latitudes
   integer, intent(out) :: ntime                         ! number of hourly time periods per file 
   real(4), allocatable :: lon1d(:)                      ! lon(lon) (degrees_east)
   real(4), allocatable :: lat1d(:)                      ! lat(lat) (degrees_north)
   real(4), allocatable :: time(:)                       ! days since 1901-01-01 00:00:00
   ! All these variables are the mean over the timestep 
   real(4), allocatable, intent(out) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), allocatable, intent(out) :: fpsn(:,:,:)      ! fpsn(lon, lat, time) photosynthesis (umol/m2/sec)
   real(4), allocatable, intent(out) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), allocatable, intent(out) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), allocatable, intent(out) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), allocatable, intent(out) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), allocatable, intent(out) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), allocatable, intent(out) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), allocatable, intent(out) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)

!
! !CALLED FROM:
! subroutine ??
!
! !REVISION HISTORY:
! Created by Melannie Hartman
!
!
! !LOCAL VARIABLES:
!EOP

      integer :: ncid				! netcdf file ID
      integer :: status				! function return status
      integer :: lat_dimid			! netcdf dimension id
      integer :: lon_dimid			! netcdf dimension id
      integer :: time_dimid			! netcdf dimension id
      integer :: lat_varid			! netcdf variable id
      integer :: lon_varid			! netcdf variable id
      integer :: time_varid			! netcdf variable id
!-----------------------------------------------------------------------

      status = nf_open(ncfilename, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, ncfilename)
         
      status = nf_inq_dimid(ncid, "time", time_dimid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_dimid(ncid, "lat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimid(ncid, "lon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_dimlen(ncid, time_dimid, ntime)
      if (status /= nf_noerr) call handle_err(status, "ntime")

      status = nf_inq_dimlen(ncid, lat_dimid, nlat)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimlen(ncid, lon_dimid, nlon)
      if (status /= nf_noerr) call handle_err(status, "lon")


      status = nf_inq_varid(ncid, "lat", lat_varid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_varid(ncid, "lon", lon_varid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_varid(ncid, "time", time_varid)
      if (status /= nf_noerr) call handle_err(status, "time")

      if (verbose .ge. 1) then
         print *, "  allocateDailyHourlyVariables: length(lat): ", nlat
         print *, "  allocateDailyHourlyVariables: length(lon): ", nlon
         print *, "  allocateDailyHourlyVariables: length(time): ", ntime
      endif

      allocate(lon1d(1:nlon))
      allocate(lat1d(1:nlat))
      allocate(time(1:nFiles*ntime))
! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: time, lat, lon
      allocate(flds(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fpsn(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsds(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsnd(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsni(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsvd(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(fsdsvi(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(pbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(qbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(rain(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(snow(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(tbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(thbot(1:nlon, 1:nlat, 1:nFiles*ntime))
      allocate(wind(1:nlon, 1:nlat, 1:nFiles*ntime))

      lon1d(:) = 0.0
      lat1d(:) = 0.0
      time(:) = 0.0
      flds(:,:,:) = 0.0
      fpsn(:,:,:) = 0.0
      fsds(:,:,:) = 0.0
      fsdsnd(:,:,:) = 0.0
      fsdsni(:,:,:) = 0.0
      fsdsvd(:,:,:) = 0.0
      fsdsvi(:,:,:) = 0.0
      pbot(:,:,:) = 0.0
      qbot(:,:,:) = 0.0
      rain(:,:,:) = 0.0
      snow(:,:,:) = 0.0
      tbot(:,:,:) = 0.0
      thbot(:,:,:) = 0.0
      wind(:,:,:) = 0.0

      status = nf_close(ncid)

   end subroutine allocateDailyHourlyVariables

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readDailyHistoryFile
!
! !INTERFACE:
!  subroutine readDailyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, lon1d, lat1d, &
!                                  area, topo, landfrac, landmask, pftmask, time, levgrnd, &
!                                  flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
!                                  rain, snow, tbot, thbot, wind, tlai, h2osoi, soilliq, tsoi)
   subroutine readDailyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, lon1d, lat1d, &
                                   area, topo, landfrac, landmask, pftmask, time, levgrnd, &
                                   flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                                   rain, snow, tbot, thbot, wind, h2osoi, soilliq, tsoi)
!
! !DESCRIPTION:
!  Read latitude, longitude and meteorlogical variables from daily netcdf history files. 
!  It is assumed that each netcdf file contains one year of data.
!
! !USES:
!  none
!
! !ARGUMENTS:
   character(len=*), intent(in)    :: ncfilename  ! netcdf filename
   integer, intent(out) :: nlon				 ! number of longitudes 
   integer, intent(out) :: nlat				 ! number of latitudes
   integer, intent(out) :: ntime			 ! number of hourly time periods per file 
   integer, intent(out) :: nlevgrnd			 ! number of soil levels
   real(4), allocatable, intent(out) :: lon1d(:)	 ! lon(lon) (degrees_east)
   real(4), allocatable, intent(out) :: lat1d(:)  	 ! lat(lat) (degrees_north)
   real(4), allocatable, intent(out) :: area(:,:)	 ! grid cell areas (km^2)
   real(4), allocatable, intent(out) :: topo(:,:)	 ! grid cell topography (m)
   real(4), allocatable, intent(out) :: landfrac(:,:)    ! land fraction
   integer, allocatable, intent(out) :: landmask(:,:)    ! land/ocean mask (0=ocean and 1=land)
   integer, allocatable, intent(out) :: pftmask(:,:)     ! pft real/fake mask (0=fake and 1=real)
   real(4), allocatable, intent(out) :: time(:)          ! days since 1901-01-01 00:00:00
   real(4), allocatable, intent(out) :: levgrnd(:)       ! coordinate soil levels (m)
   ! All these variables are the mean over the timestep
   real(4), allocatable, intent(out) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), allocatable, intent(out) :: fpsn(:,:,:)      ! flds(lon, lat, time) photosynthesis (umol/m2/sec)
   real(4), allocatable, intent(out) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), allocatable, intent(out) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), allocatable, intent(out) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), allocatable, intent(out) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), allocatable, intent(out) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), allocatable, intent(out) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), allocatable, intent(out) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
   !real(4), allocatable, intent(out) :: tlai(:,:,:)      ! tlai(lon, lat, time) total projected LAI
   real(4), allocatable, intent(out) :: h2osoi(:,:,:,:)  ! h2osoi(lon, lat, levgrnd, time) volumetric soil water (mm3/mm3)(vegetated landunits only) 
   real(4), allocatable, intent(out) :: soilliq(:,:,:,:) ! soilliq(lon, lat, levgrnd, time) soil liquid water (kg/m2)(vegetated landunits only) 
   real(4), allocatable, intent(out) :: tsoi(:,:,:,:)    ! tsoi(lon, lat, levgrnd, time) soil temperature (K) (vegetated landunits only) 

! !CALLED FROM:
! subroutine ??
!
! !REVISION HISTORY:
! Created by Melannie Hartman
!
!
! !LOCAL VARIABLES:
!EOP

      integer :: ncid				! netcdf file ID
      integer :: status				! function return status
      integer :: lat_dimid			! netcdf dimension id
      integer :: lon_dimid			! netcdf dimension id
      integer :: time_dimid			! netcdf dimension id
      integer :: levgrnd_dimid			! netcdf dimension id
      integer :: lat_varid			! netcdf variable id
      integer :: lon_varid			! netcdf variable id
      integer :: area_varid			! netcdf variable id
      integer :: topo_varid			! netcdf variable id
      integer :: levgrnd_varid			! netcdf variable id
      integer :: landfrac_varid			! netcdf variable id
      integer :: landmask_varid			! netcdf variable id
      integer :: pftmask_varid			! netcdf variable id
      integer :: time_varid			! netcdf variable id
      integer :: varid	                	! netcdf variable id
      integer :: start2(2), count2(2)           ! start and count arrays for reading 2-D data from netcdf files
      integer :: start3(3), count3(3)           ! start and count arrays for reading 3-D data from netcdf files
      integer :: start4(4), count4(4)           ! start and count arrays for reading 4-D data from netcdf files
!-----------------------------------------------------------------------

      if (verbose .ge. 0) print *, "Reading daily climate data file ", trim(ncfilename), "..."

      status = nf_open(ncfilename, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, ncfilename)
         

      status = nf_inq_dimid(ncid, "time", time_dimid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_dimid(ncid, "lat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimid(ncid, "lon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_dimid(ncid, "levgrnd", levgrnd_dimid)
      if (status /= nf_noerr) call handle_err(status, "levgrnd")


      status = nf_inq_dimlen(ncid, time_dimid, ntime)
      if (status /= nf_noerr) call handle_err(status, "ntime")

      status = nf_inq_dimlen(ncid, lat_dimid, nlat)
      if (status /= nf_noerr) call handle_err(status, "nlat")

      status = nf_inq_dimlen(ncid, lon_dimid, nlon)
      if (status /= nf_noerr) call handle_err(status, "nlon")

      status = nf_inq_dimlen(ncid, levgrnd_dimid, nlevgrnd)
      if (status /= nf_noerr) call handle_err(status, "nlevgrnd")


      status = nf_inq_varid(ncid, "lat", lat_varid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_varid(ncid, "lon", lon_varid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_varid(ncid, "area", area_varid)
      if (status /= nf_noerr) call handle_err(status, "area")

! topo is no longer included in CLM files. -mdh 12/12/2016
!     status = nf_inq_varid(ncid, "topo", topo_varid)
!     if (status /= nf_noerr) call handle_err(status, "topo")

      status = nf_inq_varid(ncid, "landfrac", landfrac_varid)
      if (status /= nf_noerr) call handle_err(status, "landfrac")

      status = nf_inq_varid(ncid, "landmask", landmask_varid)
      if (status /= nf_noerr) call handle_err(status, "landmask")

      status = nf_inq_varid(ncid, "pftmask", pftmask_varid)
      if (status /= nf_noerr) call handle_err(status, "pftmask")

      status = nf_inq_varid(ncid, "time", time_varid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_varid(ncid, "levgrnd", levgrnd_varid)
      if (status /= nf_noerr) call handle_err(status, "levgrnd")

      if (verbose .ge. 1) then
         print *, "  readDailyHistoryFile: length(lat): ", nlat
         print *, "  readDailyHistoryFile: length(lon): ", nlon
         print *, "  readDailyHistoryFile: length(time): ", ntime
         print *, "  readDailyHistoryFile: length(levgrnd): ", nlevgrnd
      endif

! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: lat, lon (2-D); time, lat, lon (3-D); time, levgrnd, lat, lon (4-D)
      allocate(lon1d(1:nlon))
      allocate(lat1d(1:nlat))
      allocate(time(1:ntime))
      allocate(levgrnd(1:nlevgrnd))
      allocate(area(1:nlon, 1:nlat))
      allocate(topo(1:nlon, 1:nlat))
      allocate(landfrac(1:nlon, 1:nlat))
      allocate(landmask(1:nlon, 1:nlat))
      allocate(pftmask(1:nlon, 1:nlat))
      allocate(flds(1:nlon, 1:nlat, 1:ntime))
      allocate(fpsn(1:nlon, 1:nlat, 1:ntime))
      allocate(fsds(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsnd(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsni(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsvd(1:nlon, 1:nlat, 1:ntime))
      allocate(fsdsvi(1:nlon, 1:nlat, 1:ntime))
      allocate(pbot(1:nlon, 1:nlat, 1:ntime))
      allocate(qbot(1:nlon, 1:nlat, 1:ntime))
      allocate(rain(1:nlon, 1:nlat, 1:ntime))
      allocate(snow(1:nlon, 1:nlat, 1:ntime))
      allocate(tbot(1:nlon, 1:nlat, 1:ntime))
      allocate(thbot(1:nlon, 1:nlat, 1:ntime))
      allocate(wind(1:nlon, 1:nlat, 1:ntime))
      !allocate(tlai(1:nlon, 1:nlat, 1:ntime))
      allocate(h2osoi(1:nlon, 1:nlat, 1:nlevgrnd, 1:ntime))
      allocate(soilliq(1:nlon, 1:nlat, 1:nlevgrnd, 1:ntime))
      allocate(tsoi(1:nlon, 1:nlat, 1:nlevgrnd, 1:ntime))
 
      lon1d(:) = 0.0
      lat1d(:) = 0.0
      time(:) = 0.0
      levgrnd(:) = 0.0

      area(:,:) = 0.0
      topo(:,:) = 0.0
      landfrac(:,:) = 0.0
      landmask(:,:) = 0
      pftmask(:,:) = 0
 
      flds(:,:,:) = 0.0
      fpsn(:,:,:) = 0.0
      fsds(:,:,:) = 0.0
      fsdsnd(:,:,:) = 0.0
      fsdsni(:,:,:) = 0.0
      fsdsvd(:,:,:) = 0.0
      fsdsvi(:,:,:) = 0.0
      pbot(:,:,:) = 0.0
      qbot(:,:,:) = 0.0
      rain(:,:,:) = 0.0
      snow(:,:,:) = 0.0
      tbot(:,:,:) = 0.0
      thbot(:,:,:) = 0.0
      wind(:,:,:) = 0.0
      !tlai(:,:,:) = 0.0

      h2osoi(:,:,:,:) = 0.0
      soilliq(:,:,:,:) = 0.0
      tsoi(:,:,:,:) = 0.0
 

      status = nf_get_var(ncid, lat_varid, lat1d)
      if (status /= nf_noerr) call handle_err(status, "lat1d")

      status = nf_get_var(ncid, lon_varid, lon1d)
      if (status /= nf_noerr) call handle_err(status, "lon1d")

      status = nf_get_var(ncid, time_varid, time)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_get_var(ncid, levgrnd_varid, levgrnd)
      if (status /= nf_noerr) call handle_err(status, "levgrnd")

!     do i = 1, nlevgrnd
!         print *, i, levgrnd(i)
!     end do

      start2 = (/ 1 ,1 /)
      count2 = (/ nlat, nlon /)

      ! area(nlon, nlat)

      status = nf_inq_varid(ncid, "area", varid)
      if (status /= nf_noerr) call handle_err(status, "area")

      status = nf_get_var(ncid, varid, area, start2, count2)
      if (status /= nf_noerr) call handle_err(status, "area")


!     ! topo(nlon, nlat)
!
!     status = nf_inq_varid(ncid, "topo", varid)
!     if (status /= nf_noerr) call handle_err(status, "topo")
!
!     status = nf_get_var(ncid, varid, topo, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "topo")


      ! landfrac(nlon, nlat)

      status = nf_inq_varid(ncid, "landfrac", varid)
      if (status /= nf_noerr) call handle_err(status, "landfrac")

      status = nf_get_var(ncid, varid, landfrac, start2, count2)
      if (status /= nf_noerr) call handle_err(status, "landfrac")


      ! landmask(nlon, nlat)

      status = nf_inq_varid(ncid, "landmask", varid)
      if (status /= nf_noerr) call handle_err(status, "landmask")

      status = nf_get_var(ncid, varid, landmask, start2, count2)
      if (status /= nf_noerr) call handle_err(status, "landmask")


      ! pftmask(nlon, nlat)

      status = nf_inq_varid(ncid, "pftmask", varid)
      if (status /= nf_noerr) call handle_err(status, "pftmask")

      status = nf_get_var(ncid, varid, pftmask, start2, count2)
      if (status /= nf_noerr) call handle_err(status, "pftmask")


      start3 = (/ 1, 1 ,1 /)
      count3 = (/ ntime, nlat, nlon /)

      ! flds(nlon, nlat, ntime) atmospheric longwave radiation (W/m^2)
 
      if (verbose .ge. 2) print *, "FLDS"
      status = nf_inq_varid(ncid, "FLDS", varid)
      if (status /= nf_noerr) call handle_err(status, "FLDS")

      status = nf_get_var(ncid, varid, flds, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "flds")


      ! fpsn(nlon, nlat, ntime) photosynthesis (umol/m2/sec)
      ! Added 8/11/2014
 
      if (verbose .ge. 2) print *, "FPSN"
      status = nf_inq_varid(ncid, "FPSN", varid)
      if (status /= nf_noerr) call handle_err(status, "FPSN")

      status = nf_get_var(ncid, varid, fpsn, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fpsn")


      ! fsds(nlon, nlat, ntime) atmospheric incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDS"
      status = nf_inq_varid(ncid, "FSDS", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDS")

      status = nf_get_var(ncid, varid, fsds, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsds")


      ! fsdsnd(nlon, nlat, ntime) direct nir incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDSND"
      status = nf_inq_varid(ncid, "FSDSND", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSND")

      status = nf_get_var(ncid, varid, fsdsnd, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsnd")


      ! fsdsni(nlon, nlat, ntime) diffuse nir incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDSNI"
      status = nf_inq_varid(ncid, "FSDSNI", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSNI")

      status = nf_get_var(ncid, varid, fsdsni, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsni")


      ! fsdsvd(nlon, nlat, ntime)direct vis incident solar radiation" (W/m^2)

      if (verbose .ge. 2) print *, "FSDSVD"
      status = nf_inq_varid(ncid, "FSDSVD", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSVD")

      status = nf_get_var(ncid, varid, fsdsvd, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsvd")


      ! fsdsvi(nlon, nlat, ntime) diffuse vis incident solar radiation (W/m^2)

      if (verbose .ge. 2) print *, "FSDSVI"
      status = nf_inq_varid(ncid, "FSDSVI", varid)
      if (status /= nf_noerr) call handle_err(status, "FSDSVI")

      status = nf_get_var(ncid, varid, fsdsvi, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "fsdsvi")


      ! pbot(nlon, nlat, ntime) atmospheric pressure (Pa)

      if (verbose .ge. 2) print *, "PBOT"
      status = nf_inq_varid(ncid, "PBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "PBOT")

      status = nf_get_var(ncid, varid, pbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "pbot")


      ! qbot(nlon, nlat, ntime) atmospheric specific humidity (kg/kg)

      if (verbose .ge. 2) print *, "QBOT"
      status = nf_inq_varid(ncid, "QBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "QBOT")

      status = nf_get_var(ncid, varid, qbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "qbot")


      ! rain(nlon, nlat, ntime) atmospheric rain (mm/sec)

      if (verbose .ge. 2) print *, "RAIN"
      status = nf_inq_varid(ncid, "RAIN", varid)
      if (status /= nf_noerr) call handle_err(status, "RAIN")

      status = nf_get_var(ncid, varid, rain, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "rain")


      ! snow(nlon, nlat, ntime) atmospheric snow (mm/sec)

      if (verbose .ge. 2) print *, "SNOW"
      status = nf_inq_varid(ncid, "SNOW", varid)
      if (status /= nf_noerr) call handle_err(status, "SNOW")

      status = nf_get_var(ncid, varid, snow, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "snow")


      ! tbot(nlon, nlat, ntime) atmospheric air temperature (K)

      if (verbose .ge. 2) print *, "TBOT"
      status = nf_inq_varid(ncid, "TBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "TBOT")

      status = nf_get_var(ncid, varid, tbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "tbot")


      ! thbot(nlon, nlat, ntime) atmospheric air potential temperature (K)

      if (verbose .ge. 2) print *, "THBOT"
      status = nf_inq_varid(ncid, "THBOT", varid)
      if (status /= nf_noerr) call handle_err(status, "THBOT")

      status = nf_get_var(ncid, varid, thbot, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "thbot")


      ! wind(nlon, nlat, ntime) atmospheric wind velocity magnitude (m/s)

      if (verbose .ge. 2) print *, "WIND"
      status = nf_inq_varid(ncid, "WIND", varid)
      if (status /= nf_noerr) call handle_err(status, "WIND")

      status = nf_get_var(ncid, varid, wind, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "wind")


      ! tlai(nlon, nlat, ntime) total projected LAI

      !if (verbose .ge. 2) print *, "TLAI"
      !status = nf_inq_varid(ncid, "TLAI", varid)
      !if (status /= nf_noerr) call handle_err(status, "TLAI")
      !
      !status = nf_get_var(ncid, varid, tlai, start3, count3)
      !if (status /= nf_noerr) call handle_err(status, "tlai")


      start4 = (/ 1, 1 ,1, 1 /)
      count4 = (/ ntime, nlevgrnd, nlat, nlon /)

      ! h2osoi(lon, lat, nlevgrnd, time) volumetric soil water content (mm3/mm3)

      if (verbose .ge. 2) print *, "H2OSOI"
      status = nf_inq_varid(ncid, "H2OSOI", varid)
      if (status /= nf_noerr) call handle_err(status, "H2OSOI")

      status = nf_get_var(ncid, varid, h2osoi, start4, count4)
      if (status /= nf_noerr) call handle_err(status, "h2osoi")


      ! soilliq(lon, lat, nlevgrnd, time) soil liquid water (kg/m2)

      if (verbose .ge. 2) print *, "SOILLIQ"
      status = nf_inq_varid(ncid, "SOILLIQ", varid)
      if (status /= nf_noerr) call handle_err(status, "SOILLIQ")

      status = nf_get_var(ncid, varid, soilliq, start4, count4)
      if (status /= nf_noerr) call handle_err(status, "soilliq")


      ! tsoi(lon, lat, nlevgrnd, time) soil temperature (K)

      if (verbose .ge. 2) print *, "TSOI"
      status = nf_inq_varid(ncid, "TSOI", varid)
      if (status /= nf_noerr) call handle_err(status, "TSOI")

      status = nf_get_var(ncid, varid, tsoi, start4, count4)
      if (status /= nf_noerr) call handle_err(status, "tsoi")


      status = nf_close(ncid)

      if (verbose .ge. 0) print *, "Done reading daily climate data from ", trim(ncfilename)

   end subroutine readDailyHistoryFile

!-----------------------------------------------------------------------
! This subroutine is commented out because it is not currently being used
! and I didn't want to confuse this code with other code I am modifying.
! Melannie Hartman 2/20/2017
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readHourlyHistoryFile
!
! !INTERFACE:
!  subroutine readHourlyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, lon1d, lat1d, &
!                                   area, topo, landfrac, landmask, pftmask, time, levgrnd, &
!                                   flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
!                                   rain, snow, tbot, thbot, wind, tlai, h2osoi, tsoi)
!
! !DESCRIPTION:
!  Read latitude, longitude and meteorlogical variables from hourly netcdf history files. 
!  It is assumed that each netcdf file contains one year of data.
!
! !USES:
!  none
!
! !ARGUMENTS:
!  character(len=*), intent(in)    :: ncfilename  ! netcdf filename
!  integer, intent(out) :: nlon				 ! number of longitudes 
!  integer, intent(out) :: nlat				 ! number of latitudes
!  integer, intent(out) :: ntime			 ! number of hourly time periods per file 
!  integer, intent(out) :: nlevgrnd			 ! number of soil levels
!  real(4), allocatable, intent(out) :: lon1d(:)	 ! lon(lon) (degrees_east)
!  real(4), allocatable, intent(out) :: lat1d(:)  	 ! lat(lat) (degrees_north)
!  real(4), allocatable, intent(out) :: area(:,:)	 ! grid cell areas (km^2)
!  real(4), allocatable, intent(out) :: topo(:,:)	 ! grid cell topography (m)
!  real(4), allocatable, intent(out) :: landfrac(:,:)    ! land fraction
!  integer, allocatable, intent(out) :: landmask(:,:)    ! land/ocean mask (0=ocean and 1=land)
!  integer, allocatable, intent(out) :: pftmask(:,:)     ! pft real/fake mask (0=fake and 1=real)
!  real(4), allocatable, intent(out) :: time(:)          ! days since 1901-01-01 00:00:00
!  real(4), allocatable, intent(out) :: levgrnd(:)       ! coordinate soil levels (m)
!  ! All these variables are the mean over the timestep
!  real(4), allocatable, intent(out) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
!  real(4), allocatable, intent(out) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
!  real(4), allocatable, intent(out) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
!  real(4), allocatable, intent(out) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
!  real(4), allocatable, intent(out) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
!  real(4), allocatable, intent(out) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
!  real(4), allocatable, intent(out) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
!  real(4), allocatable, intent(out) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
!  real(4), allocatable, intent(out) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
!  real(4), allocatable, intent(out) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
!  real(4), allocatable, intent(out) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
!  real(4), allocatable, intent(out) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
!  real(4), allocatable, intent(out) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
!  real(4), allocatable, intent(out) :: tlai(:,:,:)      ! tlai(lon, lat, time) total projected LAI
!  real(4), allocatable, intent(out) :: h2osoi(:,:,:,:)  ! h2osoi(lon, lat, levgrnd, time) volumetric soil water (mm3/mm3)(vegetated landunits only) 
!  real(4), allocatable, intent(out) :: tsoi(:,:,:,:)    ! tsoi(lon, lat, levgrnd, time) soil temperature (K) (vegetated landunits only) 
!
! !CALLED FROM:
! subroutine ??
!
! !REVISION HISTORY:
! Created by Melannie Hartman
!
!
! !LOCAL VARIABLES:
!EOP
!
!     integer :: ncid				! netcdf file ID
!     integer :: status				! function return status
!     integer :: lat_dimid			! netcdf dimension id
!     integer :: lon_dimid			! netcdf dimension id
!     integer :: time_dimid			! netcdf dimension id
!     integer :: levgrnd_dimid			! netcdf dimension id
!     integer :: lat_varid			! netcdf variable id
!     integer :: lon_varid			! netcdf variable id
!     integer :: area_varid			! netcdf variable id
!     integer :: topo_varid			! netcdf variable id
!     integer :: levgrnd_varid			! netcdf variable id
!     integer :: landfrac_varid			! netcdf variable id
!     integer :: landmask_varid			! netcdf variable id
!     integer :: pftmask_varid			! netcdf variable id
!     integer :: time_varid			! netcdf variable id
!     integer :: varid	                	! netcdf variable id
!     integer :: start2(2), count2(2)           ! start and count arrays for reading 2-D data from netcdf files
!     integer :: start3(3), count3(3)           ! start and count arrays for reading 3-D data from netcdf files
!     integer :: start4(4), count4(4)           ! start and count arrays for reading 4-D data from netcdf files
!-----------------------------------------------------------------------
!
!     if (verbose .ge. 0) print *, "Reading hourly climate data file ", trim(ncfilename), "..."
!
!     status = nf_open(ncfilename, nf_nowrite, ncid)
!     if (status /= nf_noerr) call handle_err(status, ncfilename)
!        
!
!     status = nf_inq_dimid(ncid, "time", time_dimid)
!     if (status /= nf_noerr) call handle_err(status, "time")
!
!     status = nf_inq_dimid(ncid, "lat", lat_dimid)
!     if (status /= nf_noerr) call handle_err(status, "lat")
!
!     status = nf_inq_dimid(ncid, "lon", lon_dimid)
!     if (status /= nf_noerr) call handle_err(status, "lon")
!
!     status = nf_inq_dimid(ncid, "levgrnd", levgrnd_dimid)
!     if (status /= nf_noerr) call handle_err(status, "levgrnd")
!
!
!     status = nf_inq_dimlen(ncid, time_dimid, ntime)
!     if (status /= nf_noerr) call handle_err(status, "ntime")
!
!     status = nf_inq_dimlen(ncid, lat_dimid, nlat)
!     if (status /= nf_noerr) call handle_err(status, "nlat")
!
!     status = nf_inq_dimlen(ncid, lon_dimid, nlon)
!     if (status /= nf_noerr) call handle_err(status, "nlon")
!
!     status = nf_inq_dimlen(ncid, levgrnd_dimid, nlevgrnd)
!     if (status /= nf_noerr) call handle_err(status, "nlevgrnd")
!
!
!     status = nf_inq_varid(ncid, "lat", lat_varid)
!     if (status /= nf_noerr) call handle_err(status, "lat")
!
!     status = nf_inq_varid(ncid, "lon", lon_varid)
!     if (status /= nf_noerr) call handle_err(status, "lon")
!
!     status = nf_inq_varid(ncid, "area", area_varid)
!     if (status /= nf_noerr) call handle_err(status, "area")
!
!     status = nf_inq_varid(ncid, "topo", topo_varid)
!     if (status /= nf_noerr) call handle_err(status, "topo")
!
!     status = nf_inq_varid(ncid, "landfrac", landfrac_varid)
!     if (status /= nf_noerr) call handle_err(status, "landfrac")
!
!     status = nf_inq_varid(ncid, "landmask", landmask_varid)
!     if (status /= nf_noerr) call handle_err(status, "landmask")
!
!     status = nf_inq_varid(ncid, "pftmask", pftmask_varid)
!     if (status /= nf_noerr) call handle_err(status, "pftmask")
!
!     status = nf_inq_varid(ncid, "time", time_varid)
!     if (status /= nf_noerr) call handle_err(status, "time")
!
!     status = nf_inq_varid(ncid, "levgrnd", levgrnd_varid)
!     if (status /= nf_noerr) call handle_err(status, "levgrnd")
!
!     if (verbose .ge. 1) then
!        print *, "  readHourlyHistoryFile: length(lat): ", nlat
!        print *, "  readHourlyHistoryFile: length(lon): ", nlon
!        print *, "  readHourlyHistoryFile: length(time): ", ntime
!        print *, "  readHourlyHistoryFile: length(levgrnd): ", nlevgrnd
!     endif
!
! Dimensions in FORTRAN are in Column Major Order: the first array index varies the most rapidly.
! Note, in NetCDF file the dimensions appear in the opposite order: lat, lon (2-D); time, lat, lon (3-D); time, levgrnd, lat, lon (4-D)
!     allocate(lon1d(1:nlon))
!     allocate(lat1d(1:nlat))
!     allocate(time(1:ntime))
!     allocate(levgrnd(1:nlevgrnd))
!     allocate(area(1:nlon, 1:nlat))
!     allocate(topo(1:nlon, 1:nlat))
!     allocate(landfrac(1:nlon, 1:nlat))
!     allocate(landmask(1:nlon, 1:nlat))
!     allocate(pftmask(1:nlon, 1:nlat))
!     allocate(flds(1:nlon, 1:nlat, 1:ntime))
!     allocate(fsds(1:nlon, 1:nlat, 1:ntime))
!     allocate(fsdsnd(1:nlon, 1:nlat, 1:ntime))
!     allocate(fsdsni(1:nlon, 1:nlat, 1:ntime))
!     allocate(fsdsvd(1:nlon, 1:nlat, 1:ntime))
!     allocate(fsdsvi(1:nlon, 1:nlat, 1:ntime))
!     allocate(pbot(1:nlon, 1:nlat, 1:ntime))
!     allocate(qbot(1:nlon, 1:nlat, 1:ntime))
!     allocate(rain(1:nlon, 1:nlat, 1:ntime))
!     allocate(snow(1:nlon, 1:nlat, 1:ntime))
!     allocate(tbot(1:nlon, 1:nlat, 1:ntime))
!     allocate(thbot(1:nlon, 1:nlat, 1:ntime))
!     allocate(wind(1:nlon, 1:nlat, 1:ntime))
!     allocate(tlai(1:nlon, 1:nlat, 1:ntime))
!     allocate(h2osoi(1:nlon, 1:nlat, 1:nlevgrnd, 1:ntime))
!     allocate(tsoi(1:nlon, 1:nlat, 1:nlevgrnd, 1:ntime))
!
!     lon1d(:) = 0.0
!     lat1d(:) = 0.0
!     time(:) = 0.0
!     levgrnd(:) = 0.0
!
!     area(:,:) = 0.0
!     topo(:,:) = 0.0
!     landfrac(:,:) = 0.0
!     landmask(:,:) = 0
!     pftmask(:,:) = 0
!
!     flds(:,:,:) = 0.0
!     fsds(:,:,:) = 0.0
!     fsdsnd(:,:,:) = 0.0
!     fsdsni(:,:,:) = 0.0
!     fsdsvd(:,:,:) = 0.0
!     fsdsvi(:,:,:) = 0.0
!     pbot(:,:,:) = 0.0
!     qbot(:,:,:) = 0.0
!     rain(:,:,:) = 0.0
!     snow(:,:,:) = 0.0
!     tbot(:,:,:) = 0.0
!     thbot(:,:,:) = 0.0
!     wind(:,:,:) = 0.0
!     tlai(:,:,:) = 0.0
!
!     h2osoi(:,:,:,:) = 0.0
!     tsoi(:,:,:,:) = 0.0
!
!
!     status = nf_get_var(ncid, lat_varid, lat1d)
!     if (status /= nf_noerr) call handle_err(status, "lat1d")
!
!     status = nf_get_var(ncid, lon_varid, lon1d)
!     if (status /= nf_noerr) call handle_err(status, "lon1d")
!
!     status = nf_get_var(ncid, time_varid, time)
!     if (status /= nf_noerr) call handle_err(status, "time")
!
!     status = nf_get_var(ncid, levgrnd_varid, levgrnd)
!     if (status /= nf_noerr) call handle_err(status, "levgrnd")
!
!     do i = 1, nlevgrnd
!         print *, i, levgrnd(i)
!     end do
!
!     start2 = (/ 1 ,1 /)
!     count2 = (/ nlat, nlon /)
!
!     ! area(nlon, nlat)
!
!     status = nf_inq_varid(ncid, "area", varid)
!     if (status /= nf_noerr) call handle_err(status, "area")
!
!     status = nf_get_var(ncid, varid, area, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "area")
!
!
!     ! topo(nlon, nlat)
!
!     status = nf_inq_varid(ncid, "topo", varid)
!     if (status /= nf_noerr) call handle_err(status, "topo")
!
!     status = nf_get_var(ncid, varid, topo, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "topo")
!
!
!     ! landfrac(nlon, nlat)
!
!     status = nf_inq_varid(ncid, "landfrac", varid)
!     if (status /= nf_noerr) call handle_err(status, "landfrac")
!
!     status = nf_get_var(ncid, varid, landfrac, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "landfrac")
!
!
!     ! landmask(nlon, nlat)
!
!     status = nf_inq_varid(ncid, "landmask", varid)
!     if (status /= nf_noerr) call handle_err(status, "landmask")
!
!     status = nf_get_var(ncid, varid, landmask, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "landmask")
!
!
!     ! pftmask(nlon, nlat)
!
!     status = nf_inq_varid(ncid, "pftmask", varid)
!     if (status /= nf_noerr) call handle_err(status, "pftmask")
!
!     status = nf_get_var(ncid, varid, pftmask, start2, count2)
!     if (status /= nf_noerr) call handle_err(status, "pftmask")
!
!
!     start3 = (/ 1, 1 ,1 /)
!     count3 = (/ ntime, nlat, nlon /)
!
!     ! flds(nlon, nlat, ntime) atmospheric longwave radiation (W/m^2)
!
!     if (verbose .ge. 2) print *, "FLDS"
!     status = nf_inq_varid(ncid, "FLDS", varid)
!     if (status /= nf_noerr) call handle_err(status, "FLDS")
!
!     status = nf_get_var(ncid, varid, flds, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "flds")
!
!
!     ! fsds(nlon, nlat, ntime) atmospheric incident solar radiation (W/m^2)
!
!     if (verbose .ge. 2) print *, "FSDS"
!     status = nf_inq_varid(ncid, "FSDS", varid)
!     if (status /= nf_noerr) call handle_err(status, "FSDS")
!
!     status = nf_get_var(ncid, varid, fsds, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "fsds")
!
!
!     ! fsdsnd(nlon, nlat, ntime) direct nir incident solar radiation (W/m^2)
!
!     if (verbose .ge. 2) print *, "FSDSND"
!     status = nf_inq_varid(ncid, "FSDSND", varid)
!     if (status /= nf_noerr) call handle_err(status, "FSDSND")
!
!     status = nf_get_var(ncid, varid, fsdsnd, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "fsdsnd")
!
!
!     ! fsdsni(nlon, nlat, ntime) diffuse nir incident solar radiation (W/m^2)
!
!     if (verbose .ge. 2) print *, "FSDSNI"
!     status = nf_inq_varid(ncid, "FSDSNI", varid)
!     if (status /= nf_noerr) call handle_err(status, "FSDSNI")
!
!     status = nf_get_var(ncid, varid, fsdsni, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "fsdsni")
!
!
!     ! fsdsvd(nlon, nlat, ntime)direct vis incident solar radiation" (W/m^2)
!
!     if (verbose .ge. 2) print *, "FSDSVD"
!     status = nf_inq_varid(ncid, "FSDSVD", varid)
!     if (status /= nf_noerr) call handle_err(status, "FSDSVD")
!
!     status = nf_get_var(ncid, varid, fsdsvd, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "fsdsvd")
!
!
!     ! fsdsvi(nlon, nlat, ntime) diffuse vis incident solar radiation (W/m^2)
!
!     if (verbose .ge. 2) print *, "FSDSVI"
!     status = nf_inq_varid(ncid, "FSDSVI", varid)
!     if (status /= nf_noerr) call handle_err(status, "FSDSVI")
!
!     status = nf_get_var(ncid, varid, fsdsvi, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "fsdsvi")
!
!
!     ! pbot(nlon, nlat, ntime) atmospheric pressure (Pa)
!
!     if (verbose .ge. 2) print *, "PBOT"
!     status = nf_inq_varid(ncid, "PBOT", varid)
!     if (status /= nf_noerr) call handle_err(status, "PBOT")
!
!     status = nf_get_var(ncid, varid, pbot, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "pbot")
!
!
!     ! qbot(nlon, nlat, ntime) atmospheric specific humidity (kg/kg)
!
!     if (verbose .ge. 2) print *, "QBOT"
!     status = nf_inq_varid(ncid, "QBOT", varid)
!     if (status /= nf_noerr) call handle_err(status, "QBOT")
!
!     status = nf_get_var(ncid, varid, qbot, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "qbot")
!
!
!     ! rain(nlon, nlat, ntime) atmospheric rain (mm/sec)
!
!     if (verbose .ge. 2) print *, "RAIN"
!     status = nf_inq_varid(ncid, "RAIN", varid)
!     if (status /= nf_noerr) call handle_err(status, "RAIN")
!
!     status = nf_get_var(ncid, varid, rain, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "rain")
!
!
!     ! snow(nlon, nlat, ntime) atmospheric snow (mm/sec)
!
!     if (verbose .ge. 2) print *, "SNOW"
!     status = nf_inq_varid(ncid, "SNOW", varid)
!     if (status /= nf_noerr) call handle_err(status, "SNOW")
!
!     status = nf_get_var(ncid, varid, snow, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "snow")
!
!
!     ! tbot(nlon, nlat, ntime) atmospheric air temperature (K)
!
!     if (verbose .ge. 2) print *, "TBOT"
!     status = nf_inq_varid(ncid, "TBOT", varid)
!     if (status /= nf_noerr) call handle_err(status, "TBOT")
!
!     status = nf_get_var(ncid, varid, tbot, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "tbot")
!
!
!     ! thbot(nlon, nlat, ntime) atmospheric air potential temperature (K)
!
!     if (verbose .ge. 2) print *, "THBOT"
!     status = nf_inq_varid(ncid, "THBOT", varid)
!     if (status /= nf_noerr) call handle_err(status, "THBOT")
!
!     status = nf_get_var(ncid, varid, thbot, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "thbot")
!
!
!     ! wind(nlon, nlat, ntime) atmospheric wind velocity magnitude (m/s)
!
!     if (verbose .ge. 2) print *, "WIND"
!     status = nf_inq_varid(ncid, "WIND", varid)
!     if (status /= nf_noerr) call handle_err(status, "WIND")
!
!     status = nf_get_var(ncid, varid, wind, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "wind")
!
!
!     ! tlai(nlon, nlat, ntime) total projected LAI
!
!     if (verbose .ge. 2) print *, "TLAI"
!     status = nf_inq_varid(ncid, "TLAI", varid)
!     if (status /= nf_noerr) call handle_err(status, "TLAI")
!
!     status = nf_get_var(ncid, varid, tlai, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "tlai")
!
!
!     start4 = (/ 1, 1 ,1, 1 /)
!     count4 = (/ ntime, nlevgrnd, nlat, nlon /)
!
!     ! h2osoi(lon, lat, nlevgrnd, time) volumetric soil water content (mm3/mm3)
!
!     if (verbose .ge. 2) print *, "H2OSOI"
!     status = nf_inq_varid(ncid, "H2OSOI", varid)
!     if (status /= nf_noerr) call handle_err(status, "H2OSOI")
!
!     status = nf_get_var(ncid, varid, h2osoi, start4, count4)
!     if (status /= nf_noerr) call handle_err(status, "h2osoi")
!
!
!     ! tsoi(lon, lat, nlevgrnd, time) soil temperature (K)
!
!     if (verbose .ge. 2) print *, "TSOI"
!     status = nf_inq_varid(ncid, "TSOI", varid)
!     if (status /= nf_noerr) call handle_err(status, "TSOI")
!
!     status = nf_get_var(ncid, varid, tsoi, start4, count4)
!     if (status /= nf_noerr) call handle_err(status, "tsoi")
!
!
!     status = nf_close(ncid)
!
!     if (verbose .ge. 0) print *, "Done reading hourly climate data from ", trim(ncfilename)
!
!  end subroutine readHourlyHistoryFile
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readClimateSequenceMerge
!
! !INTERFACE:
   subroutine readClimateSequenceMerge(filesin, nlon, nlat, ntimes, lon1d, lat1d, time, &
                                       flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, &
                                       pbot, qbot, rain, snow, tbot, thbot, wind)
!
! !DESCRIPTION:
!  Read a monthly, daily, or hourly sequence of climate files and store meterological  
!  variables in 3-D arrays defined in function argument list.  
!  Climate data for all years is concatenated into long multi-year arrays.
!  WARNING: For hourly climate files, more than a few years can cause arrays to grow so large
!  that machine memory is exceeded.   
!
! !USES:
!  none
!
!  !ARGUMENTS:
   character(len=*), intent(in)    :: filesin            ! name of file with climate information
   integer, intent(out) :: nlon				 ! number of longitudes 
   integer, intent(out) :: nlat				 ! number of latitudes
   integer, intent(out) :: ntimes			 ! total number of monthly/hourly time periods returned
   real(4), allocatable, intent(out) :: lon1d(:)	 ! lon(lon) (degrees_east)
   real(4), allocatable, intent(out) :: lat1d(:)  	 ! lat(lat) (degrees_north)
   real(4), allocatable, intent(out) :: time(:)          ! days since 1901-01-01 00:00:00
   ! All these variables are the mean over the timestep 
   real(4), allocatable, intent(out) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), allocatable, intent(out) :: fpsn(:,:,:)      ! fpsn(lon, lat, time) photosynthesis (umol/m2/sec)
   real(4), allocatable, intent(out) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), allocatable, intent(out) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), allocatable, intent(out) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), allocatable, intent(out) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), allocatable, intent(out) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), allocatable, intent(out) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), allocatable, intent(out) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), allocatable, intent(out) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)

   ! LOCAL VARIBLES
   integer       :: imo1, imo2
   integer       :: ntime		    	! number of time periods per file
   ! Variables used to create file names
   integer       :: iyr, imo
   integer       :: i1, i2
   character*200 :: ncfilename  	    	! netcdf climate filename
   character*200 :: ncsrfcfilename  	    	! netcdf surface data file name
   character*200 :: ncndepfilename  	    	! netcdf N deposition file name
   character*200 :: filetransCO2  	 	! name of transient CO2 file (.txt)
   character*2   :: filetype		    	! h0 | h1
   integer       :: h			        ! start position of last occurrence of "h0" or "h1" in ncfilename (count starts at 1)
   integer       :: nFiles			! number of NetCDF files to be read
   integer       :: iyear1             		! first year of climate data to read from ncfilename
   integer       :: imonth1            		! first month of climate data to read from ncfilename
   integer       :: iday1             		! first day of climate data to read from ncfilename
   integer       :: iyear2             		! last year of climate data to read from ncfilename
   integer       :: imonth2            		! last month of climate data to read from ncfilename
   integer       :: iday2              		! last day of climate data to read from ncfilename
   character*2   :: xmonth
   character*4   :: xyear
   real(4)       :: time1			! time associated with iyear1, imonth1, iday1 (days since 1901-01-01 00:00:00)
   real(4)       :: time2			! time associated with iyear2, imonth2, iday2 (days since 1901-01-01 00:00:00)
   ! Meteorological variables stored in NetCDF files (partial time series)
   ! All these variables are the mean over the timestep 
   real(4), allocatable :: lon1d_loc(:)         ! lon(lon) longitude (degrees_east)
   real(4), allocatable :: lat1d_loc(:)         ! lat(lat) latitude (degrees_north)
   real(4), allocatable :: area_loc(:,:)	! grid cell areas (km^2)
   real(4), allocatable :: topo_loc(:,:)	! grid cell topography (m)
   real(4), allocatable :: landfrac_loc(:,:)    ! land fraction
   integer, allocatable :: landmask_loc(:,:)    ! land/ocean mask (0=ocean and 1=land)
   integer, allocatable :: pftmask_loc(:,:)     ! pft real/fake mask (0=fake and 1=real)
   real(4), allocatable :: time_loc(:)          ! days since 1901-01-01 00:00:00
   real(4), allocatable :: flds_loc(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), allocatable :: fpsn_loc(:,:,:)      ! fpsn(lon, lat, time) photosynthesis (umol/m2/sec)
   real(4), allocatable :: fsds_loc(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), allocatable :: fsdsnd_loc(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), allocatable :: fsdsni_loc(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), allocatable :: fsdsvd_loc(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), allocatable :: fsdsvi_loc(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), allocatable :: pbot_loc(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), allocatable :: qbot_loc(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), allocatable :: rain_loc(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), allocatable :: snow_loc(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), allocatable :: tbot_loc(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), allocatable :: thbot_loc(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), allocatable :: wind_loc(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
   !The values below are not currently being used but are needed as a function argument to readDailyHistoryFile()
   integer              :: nlevgrnd		! number of soil levels
   real(4), allocatable :: levgrnd_loc(:)       ! coordinate soil levels (m)
   !real(4), allocatable :: tlai_loc(:,:,:)      ! tlai(lon, lat, time) total projected LAI
   real(4), allocatable :: h2osoi_loc(:,:,:,:)  ! h2osoi(lon, lat, levgrnd, time) volumetric soil water (mm3/mm3)(vegetated landunits only) 
   real(4), allocatable :: soilliq_loc(:,:,:,:) ! soilliq(lon, lat, levgrnd, time) soil liquid water (kg/m2)(vegetated landunits only) 
   real(4), allocatable :: tsoi_loc(:,:,:,:)    ! tsoi(lon, lat, levgrnd, time) soil temperature (K) (vegetated landunits only) 

   call readInitFile(filesin, filetype, ncfilename, ncsrfcfilename, ncndepfilename, &
                     filetransco2, h, &
                     iyear1, imonth1, iday1, iyear2, imonth2, iday2, time1, time2)

   ! Read data from all the file names in the monthly time series (there is one file per month)
   if (filetype .eq. 'h0') then

       nFiles = (iyear2 - iyear1 + 1) * 12 - (imonth1 - 1) - (12 - imonth2)
       if (verbose .ge. 2) print *, "nFiles =", nFiles

       call allocateMonthlyVariables(ncfilename, nFiles, nlon, nlat, ntime, lon1d, lat1d, time, &
          flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, rain, snow, tbot, thbot, wind)

       ntimes = nFiles * ntime

       i1 = 1
       i2 = ntime

       do iyr = iyear1, iyear2
          imo1 = imonth1
          if (iyr .gt. iyear1) then 
             imo1 = 1
          endif
          if ((iyear1 .lt. iyear2) .and. (iyr .ne. iyear2)) then 
             imo2 = 12
          else
             imo2 = imonth2
          endif

          do imo = imo1, imo2
             write(xyear,'(i4)') iyr	! copy iyr to string xyear
             write(xmonth,'(i2)') imo	! copy imo to string xmonth
             if (imo .lt. 10) xmonth(1:1) = '0'
             ncfilename(h+3:h+6) = xyear
             ncfilename(h+8:h+9) = xmonth
             if (verbose .ge. 2) print *, iyr, imo, trim(ncfilename)

             call readMonthlyHistoryFile(ncfilename, nlon, nlat, ntime, lon1d_loc, lat1d_loc, time_loc, &
                flds_loc, fsds_loc, fsdsnd_loc, fsdsni_loc, fsdsvd_loc, fsdsvi_loc, pbot_loc, qbot_loc, &
                rain_loc, snow_loc, tbot_loc, thbot_loc, wind_loc)

             ! Concatenate the monthly weather arrays

             lon1d(1:nlon) = lon1d_loc(1:nlon)
             lat1d(1:nlat) = lat1d_loc(1:nlat)
             time(i1:i2) = time_loc(1:ntime)
             flds(1:nlon, 1:nlat, i1:i2) = flds_loc(1:nlon, 1:nlat, 1:ntime)
             fsds(1:nlon, 1:nlat, i1:i2) = fsds_loc(1:nlon, 1:nlat, 1:ntime)
             fsdsnd(1:nlon, 1:nlat, i1:i2) = fsdsnd_loc(1:nlon, 1:nlat, 1:ntime)
             fsdsni(1:nlon, 1:nlat, i1:i2) = fsdsni_loc(1:nlon, 1:nlat, 1:ntime)
             fsdsvd(1:nlon, 1:nlat, i1:i2) = fsdsvd_loc(1:nlon, 1:nlat, 1:ntime)
             fsdsvi(1:nlon, 1:nlat, i1:i2) = fsdsvi_loc(1:nlon, 1:nlat, 1:ntime)
             pbot(1:nlon, 1:nlat, i1:i2) = pbot_loc(1:nlon, 1:nlat, 1:ntime)
             qbot(1:nlon, 1:nlat, i1:i2) = qbot_loc(1:nlon, 1:nlat, 1:ntime)
             rain(1:nlon, 1:nlat, i1:i2) = rain_loc(1:nlon, 1:nlat, 1:ntime)
             snow(1:nlon, 1:nlat, i1:i2) = snow_loc(1:nlon, 1:nlat, 1:ntime)
             tbot(1:nlon, 1:nlat, i1:i2) = tbot_loc(1:nlon, 1:nlat, 1:ntime)
             thbot(1:nlon, 1:nlat, i1:i2) = thbot_loc(1:nlon, 1:nlat, 1:ntime)
             wind(1:nlon, 1:nlat, i1:i2) = wind_loc(1:nlon, 1:nlat, 1:ntime)
             i1 = i1 + ntime
             i2 = i2 + ntime

             deallocate(lon1d_loc, lat1d_loc, time_loc, &
                        flds_loc, fsds_loc, fsdsnd_loc, fsdsni_loc, fsdsvd_loc, fsdsvi_loc, &
                        pbot_loc, qbot_loc, rain_loc, snow_loc, tbot_loc, thbot_loc, wind_loc)


          end do
       end do
   endif

   ! Read data from all the file names in the hourly time series (there is one file per year)
   if (filetype .eq. 'h1') then

       nFiles = (iyear2 - iyear1 + 1) 

       call allocateDailyHourlyVariables(ncfilename, nFiles, nlon, nlat, ntime, lon1d, lat1d, time, &
          flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, rain, snow, tbot, thbot, wind)

       ntimes = nFiles * ntime

       i1 = 1
       i2 = ntime

       do iyr = iyear1, iyear2

          write(xyear,'(i4)') iyr
          ncfilename(h+3:h+6) = xyear
          if (verbose .ge. 2) print *, iyr, trim(ncfilename)

          if (h1Daily) then
              !call readDailyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, lon1d_loc, lat1d_loc, &
              !   area_loc, topo_loc, landfrac_loc, landmask_loc, pftmask_loc, time_loc, levgrnd_loc, &
              !   flds_loc, fpsn_loc, fsds_loc, fsdsnd_loc, fsdsni_loc, fsdsvd_loc, fsdsvi_loc, pbot_loc, qbot_loc, &
              !   rain_loc, snow_loc, tbot_loc, thbot_loc, wind_loc, tlai_loc, h2osoi_loc, soilliq_loc, tsoi_loc)
              call readDailyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, lon1d_loc, lat1d_loc, &
                 area_loc, topo_loc, landfrac_loc, landmask_loc, pftmask_loc, time_loc, levgrnd_loc, &
                 flds_loc, fpsn_loc, fsds_loc, fsdsnd_loc, fsdsni_loc, fsdsvd_loc, fsdsvi_loc, pbot_loc, qbot_loc, &
                 rain_loc, snow_loc, tbot_loc, thbot_loc, wind_loc, h2osoi_loc, soilliq_loc, tsoi_loc)
          else
              write(*,*) 'netcdfTools is not currently configured to read hourly CLM history files.'
              STOP

!             call readHourlyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, lon1d_loc, lat1d_loc, &
!                area_loc, topo_loc, landfrac_loc, landmask_loc, pftmask_loc, time_loc, levgrnd_loc, &
!                flds_loc, fsds_loc, fsdsnd_loc, fsdsni_loc, fsdsvd_loc, fsdsvi_loc, pbot_loc, qbot_loc, &
!                rain_loc, snow_loc, tbot_loc, thbot_loc, wind_loc, tlai_loc, h2osoi_loc, tsoi_loc)
          endif

          ! Concatenate the hourly weather arrays
          ! There is a limit to the size of this array: on tramhill I found it to be 10 years of hourly values

          lon1d(1:nlon) = lon1d_loc(1:nlon)
          lat1d(1:nlat) = lat1d_loc(1:nlat)
          time(i1:i2) = time_loc(1:ntime)
          flds(1:nlon, 1:nlat, i1:i2) = flds_loc(1:nlon, 1:nlat, 1:ntime)
          fsds(1:nlon, 1:nlat, i1:i2) = fsds_loc(1:nlon, 1:nlat, 1:ntime)
          fsdsnd(1:nlon, 1:nlat, i1:i2) = fsdsnd_loc(1:nlon, 1:nlat, 1:ntime)
          fsdsni(1:nlon, 1:nlat, i1:i2) = fsdsni_loc(1:nlon, 1:nlat, 1:ntime)
          fsdsvd(1:nlon, 1:nlat, i1:i2) = fsdsvd_loc(1:nlon, 1:nlat, 1:ntime)
          fsdsvi(1:nlon, 1:nlat, i1:i2) = fsdsvi_loc(1:nlon, 1:nlat, 1:ntime)
          pbot(1:nlon, 1:nlat, i1:i2) = pbot_loc(1:nlon, 1:nlat, 1:ntime)
          qbot(1:nlon, 1:nlat, i1:i2) = qbot_loc(1:nlon, 1:nlat, 1:ntime)
          rain(1:nlon, 1:nlat, i1:i2) = rain_loc(1:nlon, 1:nlat, 1:ntime)
          snow(1:nlon, 1:nlat, i1:i2) = snow_loc(1:nlon, 1:nlat, 1:ntime)
          tbot(1:nlon, 1:nlat, i1:i2) = tbot_loc(1:nlon, 1:nlat, 1:ntime)
          thbot(1:nlon, 1:nlat, i1:i2) = thbot_loc(1:nlon, 1:nlat, 1:ntime)
          wind(1:nlon, 1:nlat, i1:i2) = wind_loc(1:nlon, 1:nlat, 1:ntime)
          i1 = i1 + ntime
          i2 = i2 + ntime

          deallocate(lon1d_loc, lat1d_loc, area_loc, topo_loc, &
                     landfrac_loc, landmask_loc, pftmask_loc, time_loc, &
                     flds_loc, fsds_loc, fsdsnd_loc, fsdsni_loc, fsdsvd_loc, fsdsvi_loc, &
                     pbot_loc, qbot_loc, rain_loc, snow_loc, tbot_loc, thbot_loc, wind_loc)


       end do
   endif

   end subroutine readClimateSequenceMerge

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readClimateSequence
!
! !INTERFACE:
   subroutine readClimateSequence(filesin)
!
! !DESCRIPTION:
!  Read a sequence of monthly, daily, or hourly climate files and store meterological 
!  variables in arrays to drive and create input files for other models.
!  Call readInitFile: Read initialization file, files.ini
!  Call readSurfaceData: Read netcdf surface data set
!  Call readNdepositionFile: Read netcdf N deposition file
!  For hourly meteorology:
!  o Call readHourlyHistoryFile: Read hourly meteorological variables
!  o Call summarizeDailyFromHourlyWeather: Hourly meteorology is aggregated to daily meteorology
!  o Call RunCanopyModel: Runs the aggregated canopy model each year to compute daily GPP.
!  o Call writeCASACNPfiles: Create input files for the CASACNP model  
!  For daily meteorology:
!  o Call readDailyHistoryFile: Read hourly meteorological variables
!  o Call summarizeDailyVariables: Daily meteorology is aggregated as needed
!  o Call writeCASACNPfiles: Create input files for the CASACNP model  
!  For monthly meteorology:
!  o Call readMonthlyHistoryFile: Read monthly meteorological variables
!  o Currently monthly meteorology is not used for other purposes.
!  After each year, the arrays are deallocated before reading in the next year's climate.
!
! !USES:
!  use canopy_model
!
!  !ARGUMENTS:
   character(len=*), intent(in)    :: filesin            ! name of file with climate information


   ! LOCAL VARIBLES
   integer       :: imo1, imo2
   integer       :: nlat, nlon, ntime, ndays
   integer       :: nlevgrnd                ! number of soil levels in climate file
   ! Variables used to create file names
   integer       :: iyr, imo
   character*200 :: ncfilename              ! netcdf climate filename
   character*200 :: ncsrfcfilename          ! netcdf surface data filename
   character*200 :: ncndepfilename          ! netcdf annual N deposition filename
   character*200 :: filetransCO2            ! name of transient CO2 file (.txt)
   character*200 :: fileClimOutputCSV       ! Name of climate output text file (.csv)
   character*100 :: fileGPPOutputCSV        ! Name of gpp output text file (.csv)
   character*2   :: filetype                ! 'h0' | 'h1'
   character*2   :: xmonth
   character*4   :: xyear
   integer       :: h                       ! start position of last occurrence of "h0" or "h1" in ncfilename (count starts at 1)
   integer       :: iyear1                  ! first year of climate data to read from ncfilename
   integer       :: imonth1                 ! first month of climate data to read from ncfilename
   integer       :: iday1                   ! first day of climate data to read from ncfilename
   integer       :: iyear2                  ! last year of climate data to read from ncfilename
   integer       :: imonth2                 ! last month of climate data to read from ncfilename
   integer       :: iday2                   ! last day of climate data to read from ncfilename
   integer       :: nFiles                  ! number of NetCDF files to be read
   real(4)       :: time1                   ! time associated with iyear1, imonth1, iday1 (days since 1901-01-01 00:00:00)
   real(4)       :: time2                   ! time associated with iyear2, imonth2, iday2 (days since 1901-01-01 00:00:00)
   ! Meteorological variables stored in NetCDF files (partial time series)
   real(4), allocatable :: lon1d(:)         ! lon(lon) longitude (degrees_east)
   real(4), allocatable :: lat1d(:)         ! lat(lat) latitude (degrees_north)
   real(4), allocatable :: area(:,:)        ! grid cell areas (km^2)
   real(4), allocatable :: topo(:,:)        ! grid cell topography (m)
   real(4), allocatable :: landfrac(:,:)    ! land fraction
   integer, allocatable :: landmask(:,:)    ! land/ocean mask (0=ocean and 1=land)
   integer, allocatable :: pftmask(:,:)     ! pft real/fake mask (0=fake and 1=real)
   real(4), allocatable :: time(:)          ! days since 1901-01-01 00:00:00
   real(4), allocatable :: levgrnd(:)       ! coordinate soil levels (m)
   ! All these variables are the mean over the timestep, and time is in hours or days.
   real(4), allocatable :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), allocatable :: fpsn(:,:,:)      ! fpsn(lon, lat, time) photosynthesis (umol/m2/sec)
   real(4), allocatable :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), allocatable :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), allocatable :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), allocatable :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), allocatable :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), allocatable :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), allocatable :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), allocatable :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), allocatable :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), allocatable :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), allocatable :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), allocatable :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
   !real(4), allocatable :: tlai(:,:,:)      ! tlai(lon, lat, time) total projected LAI
   real(4), allocatable :: h2osoi(:,:,:,:)  ! h2osoi(lon, lat, levgrnd, time) volumetric soil water (mm3/mm3)(vegetated landunits only) 
   real(4), allocatable :: soilliq(:,:,:,:) ! soilliq(lon, lat, levgrnd, time) soil liquid water (kg/m2) (vegetated landunits only) 
   real(4), allocatable :: tsoi(:,:,:,:)    ! tsoi(lon, lat, levgrnd, time) soil temperature BY LAYER (K) (vegetated landunits only) 

   ! Surface data set
   integer :: lsmlon                        ! length of longitude dimensions 
   integer :: lsmlat                        ! length of latitude dimension 
   integer :: lsmpft                        ! length of PFT dimension 
   integer :: nlevsoi                       ! number of soil levels in surface data set
   ! Note: pctpft, pctsand, pctclay, omfrac from surface data set are DOUBLE. Use real(8).
   real(8), allocatable :: pctpft(:,:,:)    ! Percent of each PFT in grid cell (lsmlon,lsmlat,lsmpft)
   real(8), allocatable :: pctsand(:,:,:)   ! Percent Sand (lsmlon,lsmlat,nlevsoi)
   real(8), allocatable :: pctclay(:,:,:)   ! Percent Clay (lsmlon,lsmlat,nlevsoi)
   real(8), allocatable :: orgdens(:,:,:)   ! organic matter density at soil levels (kg/m3) (lsmlon,lsmlat,nlevsoi)
!! WW added, floats
   real(8), allocatable :: watsat(:,:,:)   ! (lsmlon,lsmlat,levgrnd)
   real(8), allocatable :: watfc(:,:,:)   ! (lsmlon,lsmlat,levgrnd)

   ! Daily weather and soil water variables
   real(4), allocatable :: tminday(:,:,:)   ! tminday(lon, lat, day) minimum daily air temperature (C)
   real(4), allocatable :: tmaxday(:,:,:)   ! tmaxday(lon, lat, day) maximum daily air temperature (C)
   real(4), allocatable :: rainday(:,:,:)   ! rainday(lon, lat, day) daily atmospheric rain (mm/day) 
   real(4), allocatable :: snowday(:,:,:)   ! snowday(lon, lat, day) daily atmospheric snow (mm/day) 
   real(4), allocatable :: fsdsday(:,:,:)   ! fsdsday(lon, lat, day) 24-hour average atmospheric incident solar radiation (W/m^2)
   real(4), allocatable :: windday(:,:,:)   ! windday(lon, lat, day) 24-hour average atmospheric wind velocity magnitude (m/s)
   !real(4), allocatable :: tlaiday(:,:,:)   ! tlaiday(lon, lat, day) daily total projected LAI
   real(4), allocatable :: gppday(:,:,:)    ! gppday(lon, lat, day) Gross Primary Production (gC/m2/day)
   real(4), allocatable :: soipsiday(:,:,:) ! soipsiday(lon,lat,day) daily depth-weighted soil water potential(vegetated landunits only) (MPa)
   real(4), allocatable :: h2osoiday(:,:,:,:)  ! h2osoiday(lon, lat, nlevsoi, day) average volumetric liquid soil water BY LAYER (mm3/mm3) (vegetated landunits only) 
   real(4), allocatable :: h2ofrznsoiday(:,:,:,:) ! h2ofrznsoiday(lon, lat, nlevsoi, day) average volumetric frzn soil water BY LAYER (mm3/mm3) (vegetated landunits only) 
!  real(4), allocatable :: soilliqday(:,:,:,:) ! soilliqday(lon, lat, nlevsoi, day) soil liquid water BY LAYER (kg/m2) (vegetated landunits only) 
   real(4), allocatable :: tsoiday(:,:,:,:)    ! tsoiday(lon, lat, nlevsoi, day) soil temperature BY LAYER (K) (vegetated landunits only) 

   ! Transient CO2
   integer :: y, io
   integer :: ystart                        ! index where co2year(ystart) = iyear, the first year of climate file
   integer :: co2year(MAXCO2YRS)            ! calendar years read from transient CO2 file
   real(4) :: co2ppm(MAXCO2YRS)             ! co2 concentrations (ppm) read from transient CO2 file

   ! Annual N deposition. Time is in years.
   integer :: ndlat, ndlon                  ! ndlat and ndlon should equal nlat, nlon from other NetCDF files
   integer :: ndtime                        ! Number of years in N deposition file
   ! ndep_lon1d, ndep_lat1d, and ndep_time must be DOUBLE, but ndep is real(4)
   real(8), allocatable :: ndep_lon1d(:)    ! lon(lon) longitude (degrees_east)
   real(8), allocatable :: ndep_lat1d(:)    ! lat(lat) latitude (degrees_north)
   real(8), allocatable :: ndep_time(:)     ! days since 0000-01-01 00:00
   integer, allocatable :: ndep_calyear(:)  ! calendar year for each time
   real(4), allocatable :: ndep(:,:,:)      ! ndep(lon, lat, time) annual N deposition (gN/m2/yr)

   !Output .csv file for testing purposes (see subroutine writeClimateToCSVFile)
   if (doWriteHourlyClimateTestFile .eq. 1) then
      fileClimOutputCSV = 'climtest.csv'
      open(unit=11, file=trim(fileClimOutputCSV))
      !write(unit=11,fmt='(23(a7))') 'lon,', 'lat,', 'time,', 'tindex,', 'year,', 'day,', 'hour,', & 
      !   'flds,', 'fpsn', 'fsds,', 'fsdsnd,', 'fsdsni,', 'fsdsvd,', 'fsdsvi,', &
      !   'pbot,', 'qbot,', 'rain,', 'snow,', 'tbot,', 'thbot,', 'wind,', 'tlai,', & 
      !   'time2'
      write(unit=11,fmt='(22(a7))') 'lon,', 'lat,', 'time,', 'tindex,', 'year,', 'day,', 'hour,', & 
         'flds,', 'fpsn', 'fsds,', 'fsdsnd,', 'fsdsni,', 'fsdsvd,', 'fsdsvi,', &
         'pbot,', 'qbot,', 'rain,', 'snow,', 'tbot,', 'thbot,', 'wind,', 'time2'
   endif

   !Output .csv file for testing purposes (see subroutine RunCanopyModel)
   if (doWriteDailyGPPTestFile .eq. 1) then
      fileGPPOutputCSV = 'gpp.csv'
      open(unit=16, file=fileGPPOutputCSV)
      !write(unit=16, fmt='(12(a7))') 'lat,', 'lon,', 'year,', 'day,', 'daylen,', &
      !                               'tmin,', 'tmax,', 'rad,', 'lwpmin,', 'swp,', 'tlai,', 'gpp,'
      write(unit=16, fmt='(11(a7))') 'lat,', 'lon,', 'year,', 'day,', 'daylen,', &
                                     'tmin,', 'tmax,', 'rad,', 'lwpmin,', 'swp,', 'gpp,'
   endif

   call readInitFile(filesin, filetype, ncfilename, ncsrfcfilename, ncndepfilename, &
                     filetransco2, h, iyear1, imonth1, iday1, iyear2, imonth2, iday2, &
                     time1, time2)
! WW added watsat,watfc  
   call readSurfaceData(ncsrfcfilename, lsmlon, lsmlat, lsmpft, nlevsoi, &
                        pctpft, pctsand, pctclay, orgdens, watsat, watfc)

   call readNdepositionFile(ncndepfilename, ndlon, ndlat, ndtime, ndep_lon1d, ndep_lat1d, &
                            ndep_time, ndep_calyear, ndep)

   ! Check for consistent grid sizes
   if (ndlat .ne. lsmlat) then
      write(*,*) 'The size of the latitude dimension in ', ncndepfilename, ':', ndlat
      write(*,*) 'does not equal the size of the latitude dimension in ', ncsrfcfilename, ':', lsmlat
      STOP
   endif
   if (ndlon .ne. lsmlon) then
      write(*,*) 'The size of the longitude dimension in ', ncndepfilename, ':', ndlon
      write(*,*) 'does not equal the size of the longitude dimension in ', ncsrfcfilename, ':', lsmlon
      STOP
   endif

   !Read transient CO2 file
   open(unit=14, file=trim(filetransco2))
   rewind(unit=14)
   y = 0
   ystart = 0
   co2year(:) = 0
   co2ppm(:) = 0
   do while (y .lt. MAXCO2YRS)
       y = y + 1
       read(unit=14,fmt='(i4,2x,f8.4)',IOSTAT=io) co2year(y), co2ppm(y)
       if (io < 0) exit
       if (co2year(y) .eq. iyear1) ystart = y      ! find position of first climate year
!      print *, y, co2year(y), co2ppm(y)
   end do
!  print*, 'ystart = ', ystart
   close(unit=14)
   if (ystart .eq. 0) then
       print*, iyear1, "was not found in file ", filetransco2
       stop
   endif

   ! Read data from all the file names in the monthly time series (there is one file per month)
   if (filetype .eq. 'h0') then

       nFiles = (iyear2 - iyear1 + 1) * 12 - (imonth1 - 1) - (12 - imonth2)
       if (verbose .ge. 2) print *, "nFiles =", nFiles

       do iyr = iyear1, iyear2

          imo1 = imonth1
          if (iyr .gt. iyear1) then 
             imo1 = 1
          endif
          if ((iyear1 .lt. iyear2) .and. (iyr .ne. iyear2)) then 
             imo2 = 12
          else
             imo2 = imonth2
          endif

          do imo = imo1, imo2
             write(xyear,'(i4)') iyr    ! copy iyr to string xyear
             write(xmonth,'(i2)') imo   ! copy imo to string xmonth
             if (imo .lt. 10) xmonth(1:1) = '0'
             ncfilename(h+3:h+6) = xyear
             ncfilename(h+8:h+9) = xmonth
             if (verbose .ge. 2) print *, iyr, imo, trim(ncfilename)

             call readMonthlyHistoryFile(ncfilename, nlon, nlat, ntime, lon1d, lat1d, time, &
                flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                rain, snow, tbot, thbot, wind)

             ! Check for consistent grid sizes
             if (nlat .ne. lsmlat) then
                write(*,*) 'The size of the latitude dimension in ', ncfilename, ':', nlat
                write(*,*) 'does not equal the size of the latitude dimension in ', ncsrfcfilename, ':', lsmlat
                STOP
             endif
             if (nlon .ne. lsmlon) then
                write(*,*) 'The size of the longitude dimension in ', ncfilename, ':', nlon
                write(*,*) 'does not equal the size of the longitude dimension in ', ncsrfcfilename, ':', lsmlon
                STOP
             endif

! Write monthly weather to .csv file for testing purposes
!            call writeClimateToCSVFile(fileClimOutput, nlon, nlat, ntime, lon1d, lat1d, time, &
!               flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
!               rain, snow, tbot, thbot, wind)

             deallocate(lon1d, lat1d, time, &
                flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, &
                pbot, qbot, rain, snow, tbot, thbot, wind)

          end do
       end do
   endif

   ! Read data from all the file names in the hourly time series (there is one file per year)
   if (filetype .eq. 'h1') then
       nFiles = (iyear2 - iyear1 + 1) 

       do iyr = iyear1, iyear2

          write(xyear,'(i4)') iyr
          ncfilename(h+3:h+6) = xyear
          if (verbose .ge. 2) print *, iyr, trim(ncfilename)

          if (h1Daily) then
              !topo is no longer read. Set to zero in this call. -mdh 12/12/2016
              !call readDailyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, &
              !    lon1d, lat1d, area, topo, landfrac, landmask, pftmask, time, levgrnd, &
              !    flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
              !    rain, snow, tbot, thbot, wind, tlai, h2osoi, soilliq, tsoi)
              call readDailyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, &
                  lon1d, lat1d, area, topo, landfrac, landmask, pftmask, time, levgrnd, &
                  flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                  rain, snow, tbot, thbot, wind, h2osoi, soilliq, tsoi)
          else
              write(*,*) 'netcdfTools is not currently configured to read hourly CLM history files.'
              STOP

!             call readHourlyHistoryFile(ncfilename, nlon, nlat, ntime, nlevgrnd, &
!                 lon1d, lat1d, area, topo, landfrac, landmask, pftmask, time, levgrnd, &
!                 flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
!                 rain, snow, tbot, thbot, wind, tlai, h2osoi, tsoi)
          endif


             ! Check for consistent grid sizes
             if (nlat .ne. lsmlat) then
                write(*,*) 'The size of the latitude dimension in ', ncfilename, ':', nlat
                write(*,*) 'does not equal the size of the latitude dimension in ', ncsrfcfilename, ':', lsmlat
                STOP
             endif
             if (nlon .ne. lsmlon) then
                write(*,*) 'The size of the longitude dimension in ', ncfilename, ':', nlon
                write(*,*) 'does not equal the size of the longitude dimension in ', ncsrfcfilename, ':', lsmlon
                STOP
             endif

          if (h1Daily) then
              ! Summarize the daily variables in the clm daily history file (mostly soil variables).
              ! Daily clm history files contain fpsn (used to compute daily GPP).  
              ! No need to call aggregated canopy model to compute GPP.

              !call summarizeDailyVariables(nlon, nlat, ntime, nlevgrnd, nlevsoi, &
              !    lon1d, lat1d, time, levgrnd, &
              !    flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
              !    rain, snow, tbot, thbot, wind, tlai, h2osoi, soilliq, tsoi, &
              !    pctsand, pctclay, orgdens, &
              !    tminday, tmaxday, rainday, snowday, fsdsday, windday, &
              !    tlaiday, gppday, soipsiday, h2osoiday, h2ofrznsoiday, tsoiday)
              call summarizeDailyVariables(nlon, nlat, ntime, nlevgrnd, nlevsoi, &
                  lon1d, lat1d, time, levgrnd, &
                  flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                  rain, snow, tbot, thbot, wind, h2osoi, soilliq, tsoi, &
                  pctsand, pctclay, orgdens, &
                  tminday, tmaxday, rainday, snowday, fsdsday, windday, &
                  gppday, soipsiday, h2osoiday, h2ofrznsoiday, tsoiday)

                  ndays = ntime
          else

              write(*,*) 'netcdfTools is not currently configured to read hourly CLM history files.'
              STOP

!             call summarizeDailyFromHourlyWeather(nlon, nlat, ntime, nlevgrnd, nlevsoi, &
!                 lon1d, lat1d, time, levgrnd, &
!                 flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
!                 rain, snow, tbot, thbot, wind, tlai, h2osoi, tsoi, &
!                 pctsand, pctclay, orgdens, &
!                 tminday, tmaxday, rainday, snowday, fsdsday, windday, &
!                 tlaiday, soipsiday, h2osoiday, tsoiday)
!
!             ! Run the aggregated canopy model
!   
!             if (verbose .ge. 0) print *, "Running Canopy Model for year ", iyr
!             if (co2year(ystart) .ne. iyr) then
!                 print*, "There is no CO2 concentration for year ", iyr
!                 stop
!             endif
!             ndays = 365
!             call RunCanopyModel(iyr, nlon, nlat, ndays, lon1d, lat1d, &
!                 area, topo, landfrac, landmask, pftmask, &
!                 tminday, tmaxday, fsdsday, tlaiday, soipsiday, co2ppm(ystart), gppday)
!             ystart = ystart+1
!   
           endif

          if (doWriteHourlyClimateTestFile .eq. 1) then
              ! Write hourly weather to .csv file for testing purposes
              !call writeClimateToCSVFile(fileClimOutputCSV, nlon, nlat, ntime, lon1d, lat1d, time, &
              !    flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
              !    rain, snow, tbot, thbot, wind, tlai)
              call writeClimateToCSVFile(fileClimOutputCSV, nlon, nlat, ntime, lon1d, lat1d, time, &
                  flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, qbot, &
                  rain, snow, tbot, thbot, wind)
          endif

          ! Create input files for the CASACNP model
     
          !call writeCASACNPfiles(iyr, iyear1, iyear2, nlon, nlat, lsmpft, &
          !    nlevgrnd, nlevsoi, ndays, ndtime, ndep_calyear, &
          !    lon1d, lat1d, levgrnd, area, landfrac, landmask, pftmask, &
          !    ndep, tminday, tmaxday, tlaiday, gppday, &
          !    h2osoiday, h2ofrznsoiday, tsoiday, pctpft, pctsand, pctclay, &
          !    watsat, watfc)
          call writeCASACNPfiles(iyr, iyear1, iyear2, nlon, nlat, lsmpft, &
              nlevgrnd, nlevsoi, ndays, ndtime, ndep_calyear, &
              lon1d, lat1d, levgrnd, area, landfrac, landmask, pftmask, &
              ndep, tminday, tmaxday, gppday, &
              h2osoiday, h2ofrznsoiday, tsoiday, pctpft, pctsand, pctclay, &
              watsat, watfc)
!!WW added

          !Deallocate all variables that are read from the hourly climate file each year 
          !and the daily variables derived from these hourly variables.

          !deallocate(lon1d, lat1d, area, topo, landfrac, landmask, pftmask, time, levgrnd, &
          !   flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, &
          !   pbot, qbot, rain, snow, tbot, thbot, wind, tlai, h2osoi, soilliq, tsoi, &
          !   tminday, tmaxday, rainday, snowday, fsdsday, windday, tlaiday, gppday, &
          !   soipsiday, h2osoiday, h2ofrznsoiday, tsoiday)
          deallocate(lon1d, lat1d, area, topo, landfrac, landmask, pftmask, time, levgrnd, &
             flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, &
             pbot, qbot, rain, snow, tbot, thbot, wind, h2osoi, soilliq, tsoi, &
             tminday, tmaxday, rainday, snowday, fsdsday, windday, gppday, &
             soipsiday, h2osoiday, h2ofrznsoiday, tsoiday)

       end do
   endif
  
   if (doWriteHourlyClimateTestFile .eq. 1) close(unit=11)
   if (doWriteDailyGPPTestFile .eq. 1) close(unit=16)

   deallocate(pctsand, pctclay, orgdens, watsat, watfc)
! WW added

   end subroutine readClimateSequence

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: writeClimateToCSVFile
!
! !INTERFACE:
   !subroutine writeClimateToCSVFile(fileClimOutput, nlon, nlat, ntimes, lon1d, lat1d, time, &
   !                              flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, &
   !                              qbot, rain, snow, tbot, thbot, wind, tlai)
   subroutine writeClimateToCSVFile(fileClimOutput, nlon, nlat, ntimes, lon1d, lat1d, time, &
                                 flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, &
                                 qbot, rain, snow, tbot, thbot, wind)
!
! !DESCRIPTION:
!  Write a series of meterological variables in 3-D arrays to a comma delimted file.
!
! !USES:
!  none
!
!  !ARGUMENTS:
   character(len=*), intent(in)    :: fileClimOutput     ! output file name (.csv)
   integer, intent(in) :: nlon                           ! number of longitudes 
   integer, intent(in) :: nlat                           ! number of latitudes
   integer, intent(in) :: ntimes                         ! total number of monthly/hourly time periods 
   real(4), intent(in) :: lon1d(:)         ! lon(lon) (degrees_east)
   real(4), intent(in) :: lat1d(:)         ! lat(lat) (degrees_north)
   real(4), intent(in) :: time(:)          ! days since 1901-01-01 00:00:00
   ! All these variables are the mean over the timestep 
   real(4), intent(in) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), intent(in) :: fpsn(:,:,:)      ! fpsn(lon, lat, time) photosynthesis (umol/m2/sec)
   real(4), intent(in) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), intent(in) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), intent(in) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), intent(in) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), intent(in) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), intent(in) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), intent(in) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), intent(in) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), intent(in) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), intent(in) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), intent(in) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), intent(in) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
   !real(4), intent(in) :: tlai(:,:,:)      ! tlai(lon, lat, time) total projected LAI

   ! LOCAL VARIBLES
   integer       :: i, j, t
   character*1   :: c           ! delimeter (a comma)
   integer :: calyear           ! calendar year
   integer :: day               ! day of year (1..365)
   integer :: month             ! month (1..12)
   integer :: dayofmo           ! day of month (1..31)
   integer :: hour              ! hour of the day (1..24), hour 1 is 1:00am, hour 24 is midnight of the next day
   real(4) :: time2

   c = ','

!  This code goes into the calling subroutine
!  if (doWriteHourlyClimateTestFile .eq. 1) then
!     open(unit=11, file=trim(fileClimOutput))
!     write(unit=11,fmt='(23(a7))') 'lon,', 'lat,', 'time,', 'tindex,', 'year,', 'day,', 'hour,', & 
!        'flds,', 'fpsn', 'fsds,', 'fsdsnd,', 'fsdsni,', 'fsdsvd,', 'fsdsvi,', &
!        'pbot,', 'qbot,', 'rain,', 'snow,', 'tbot,', 'thbot,', 'wind,', 'tlai,', & 
!        'time2'
!  endif

    ! Write weather variables to .csv file for testing purposes
    do j = LONINDX,LONINDX
       do i = LATINDX,LATINDX
          do t = 1, ntimes
                call getYearMonthDayHour(time(t), calyear, day, month, dayofmo, hour)
                call getTime(calyear, day, hour, time2)        ! Test: retrieving time2 to compare to time(t). They should be equal.
  
                write(unit=11, fmt='(3(f0.2,a1), 4(i6,a1), 9(f0.6,a1), 2(e14.8,a1), 5(f0.6,a1))') &
                    lon1d(j),c, lat1d(i),c, time(t),c, &
                    t,c, calyear,c, day,c, hour,c,  &
                    flds(j,i,t),c, fpsn(j,i,t),c, fsds(j,i,t),c, fsdsnd(j,i,t),c, fsdsni(j,i,t),c, fsdsvd(j,i,t),c,  &
                    fsdsvi(j,i,t),c, pbot(j,i,t),c, qbot(j,i,t),c, rain(j,i,t),c, snow(j,i,t),c,  &
                    tbot(j,i,t),c, thbot(j,i,t),c, wind(j,i,t),c, time2
                    !tbot(j,i,t),c, thbot(j,i,t),c, wind(j,i,t),c, tlai(j,i,t),c, time2

             end do
        end do
    end do

    end subroutine writeClimateToCSVFile


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readInitFile
!
! !INTERFACE:
   subroutine readInitFile(filesin, filetype, ncfilename, ncsrfcfilename, ncndepfilename, &
                           filetransCO2, hindex, &
                           iyear1, imonth1, iday1, iyear2, imonth2, iday2, time1, time2)
!
! !DESCRIPTION:
!  Read the initialization file (filesin).  Retrieve the following from filesin:
!    climate file file type (line 1)
!    starting and ending dates (lines 2 and 3)
!    path and file name of one of the netcdf climate files (line 4)*  
!    exact path and name of surface data file (line 5)
!    exact path and name of N deposition file (line 6)
!    exact path and name of transient CO2 file name (line 7)
!
!  *The year and month in the climate file name will be substituted according to start and end dates.
!
!  Example file format:
!
!  h1        ! file type (h0=monthly, h1=hourly)
!  20000101  ! start (yyyymmdd)
!  20020101  ! end (yyyymmdd)
!  /project/bgc01/melannie/clm45sp_2deg4506_hist.clm2.h1.1902-01-01-00000.nc
!  /project/bgc01/melannie/surfdata_1.9x2.5_simyr1850_c130421.nc
!  /project/bgc01/melannie/fndep_clm_hist_simyr1849-2006_1.9x2.5_c100428.nc
!  ./CO2_1768-2010.txt
!
! !USES:
!  none
!
!  !ARGUMENTS:
   character(len=*), intent(in)    :: filesin            ! name of file with climate information
   character(len=*), intent(out)   :: filetype           ! type of climate file to read: 'h0' (monthly) or 'h1' (daily or hourly)
   character(len=*), intent(out)   :: ncfilename         ! name of first climate netcdf climate file to read
   character(len=*), intent(out)   :: ncsrfcfilename     ! name of surface data netcdf climate file to read
   character(len=*), intent(out)   :: ncndepfilename     ! name of annual N deposition netcdf file
   character(len=*), intent(out)   :: filetransCO2       ! name of transient CO2 file (.txt)
   integer, intent(out)            :: hindex             ! start position of last occurrence of "h0" or "h1" in ncfilename (count starts at 1)
   integer, intent(out)            :: iyear1             ! first year of climate data to read from ncfilename
   integer, intent(out)            :: imonth1            ! first month of climate data to read from ncfilename
   integer, intent(out)            :: iday1              ! first day of climate data to read from ncfilename
   integer, intent(out)            :: iyear2             ! last year of climate data to read from ncfilename
   integer, intent(out)            :: imonth2            ! last month of climate data to read from ncfilename
   integer, intent(out)            :: iday2              ! last day of climate data to read from ncfilename
   real(4), intent(out)            :: time1              ! time associated with iyear1, imonth1, iday1 (days since 1901-01-01 00:00:00)
   real(4), intent(out)            :: time2              ! time associated with iyear2, imonth2, iday2 (days since 1901-01-01 00:00:00)

!  !LOCAL VARIABLES

   character*8   :: start_date, end_date        ! yyyymmdd
   character*4   :: start_year, end_year        ! yyyy
   character*2   :: start_month, end_month      ! mm
   character*2   :: start_day, end_day          ! dd
   logical       :: back                        ! used as argument to index function
   integer       :: ihour1, ihour2              ! 1..24
   integer       :: idoy1, idoy2                ! 1 - 365
   integer       :: doy(0:12)                   ! last day of each month (day of year)

   if (verbose .ge. 0) print *, "Reading initialization file ", trim(filesin), "..."
   open(unit=10, file=trim(filesin))
   rewind(unit=10)
   read(unit=10, fmt='(a2)') filetype           ! 'h0' or 'h1'
   read(unit=10, fmt='(a8)') start_date         ! yyyymmdd
   read(unit=10, fmt='(a8)') end_date           ! yyyymmdd
   read(unit=10, fmt='(a200)') ncfilename
   read(unit=10, fmt='(a200)') ncsrfcfilename
   read(unit=10, fmt='(a200)') ncndepfilename
   read(unit=10, fmt='(a200)') filetransCO2
   ncfilename = trim(ncfilename)
   ncsrfcfilename = trim(ncsrfcfilename)
   ncndepfilename = trim(ncndepfilename)
   filetransCO2 = trim(filetransCO2)
   close(unit=10)

   back = .true.
   hindex = index(ncfilename, filetype, back)   ! start position of last occurrence of "h0" or "h1" in ncfilename (count starts at 1)
   if (hindex .eq. 0) then
       write(*,*) 'Substring ', filetype, ' was not found in netCDF file name ' , ncfilename
       STOP
   endif
   if (verbose .ge. 2) then
      print *, 'index = ', hindex
      print *, 'start_date = ', start_date
      print *, 'end_date = ', end_date
   endif
  
   doy(0) = 0
   doy(1) = doy(0) + 31
   doy(2) = doy(1) + 28
   doy(3) = doy(2) + 31
   doy(4) = doy(3) + 30
   doy(5) = doy(4) + 31
   doy(6) = doy(5) + 30
   doy(7) = doy(6) + 31
   doy(8) = doy(7) + 31
   doy(9) = doy(8) + 30
   doy(10) = doy(9) + 31
   doy(11) = doy(10) + 30
   doy(12) = doy(11) + 31

   ! INSERT RANGE CHECKING ON imonth1, imonth2

   start_year = start_date(1:4)
   start_month = start_date(5:6)
   start_day = start_date(7:8)
   end_year = end_date(1:4)
   end_month = end_date(5:6)
   end_day = end_date(7:8)
   if (verbose .ge. 1) then
      print *, '  start_year = ', start_year
      print *, '  start_month = ', start_month
      print *, '  start_day = ', start_day
      print *, '  end_year = ', end_year
      print *, '  end_month = ', end_month
      print *, '  end_day = ', end_day
   endif

   ! Copy string variables (start_*, end_*) into integer variables
   read(start_year, '(i4)') iyear1
   read(start_month, '(i4)') imonth1
   read(start_day,  '(i4)') iday1
   read(end_year, '(i4)') iyear2
   read(end_month, '(i4)') imonth2
   read(end_day,  '(i4)') iday2

   ihour1 = 1
   ihour2 = 24
   idoy1 = doy(imonth1-1) + iday1 
   idoy2 = doy(imonth2-1) + iday2 
   if (filetype .eq. 'h0') then
      ihour1 = 24
      ihour2 = 24
      idoy1 = doy(imonth1) 
      idoy2 = doy(imonth2) 
   endif
   
   call getTime(iyear1, idoy1, ihour1, time1)
   call getTime(iyear2, idoy2, ihour2, time2)

   if (verbose .ge. 2) then
      print *, 'time1 = ', time1
      print *, 'time2 = ', time2
   endif

   if (verbose .ge. 0) print *, "Done reading initialization file ", trim(filesin)
  
   end subroutine readInitFile

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: summarizeDailyVariables
!
! !INTERFACE:
   !subroutine summarizeDailyVariables(nlon, nlat, ntimes, nlevgrnd, nlevsoi, &
   !                              lon1d, lat1d, time, levgrnd, &
   !                              flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, &
   !                              qbot, rain, snow, tbot, thbot, wind, &
   !                              tlai, h2osoi, soilliq, tsoi, pctsand, pctclay, orgdens, &
   !                              tminday, tmaxday, rainday, snowday, fsdsday, windday, &
   !                              tlaiday, gppday, soipsiday, h2osoiday, h2ofrznsoiday, tsoiday)
   subroutine summarizeDailyVariables(nlon, nlat, ntimes, nlevgrnd, nlevsoi, &
                                 lon1d, lat1d, time, levgrnd, &
                                 flds, fpsn, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, &
                                 qbot, rain, snow, tbot, thbot, wind, &
                                 h2osoi, soilliq, tsoi, pctsand, pctclay, orgdens, &
                                 tminday, tmaxday, rainday, snowday, fsdsday, windday, &
                                 gppday, soipsiday, h2osoiday, h2ofrznsoiday, tsoiday)
!
! !DESCRIPTION:
!  Average and aggregate daily weather variables to compute daily weather variables.
!  Compute daily depth-weighted soil water potential from soilliq (not h2osoi), pctsand, pctclay, and orgdens.
!
! !USES:
!  none
!
!  !ARGUMENTS:
   integer, intent(in) :: nlon		   ! number of longitudes 
   integer, intent(in) :: nlat		   ! number of latitudes
   integer, intent(in) :: ntimes	   ! total number of hourly time periods (8760 hours - 1 year)
   integer, intent(in) :: nlevgrnd	   ! number of soil levels from history file
   integer, intent(in) :: nlevsoi	   ! number of soil levels from surface data set
   real(4), intent(in) :: lon1d(:)	   ! longitude(lon) (degrees_east)
   real(4), intent(in) :: lat1d(:)  	   ! latitude(lat) (degrees_north)
   real(4), intent(in) :: time(:)          ! days since 1901-01-01 00:00:00
   real(4), intent(in) :: levgrnd(:)       ! coordinate soil levels (m)
   ! All the hourly variables are the mean over the timestep 
   real(4), intent(in) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
   real(4), intent(in) :: fpsn(:,:,:)      ! fpsn(lon, lat, time) photosynthesis (umol/m2/sec)
   real(4), intent(in) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
   real(4), intent(in) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
   real(4), intent(in) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
   real(4), intent(in) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
   real(4), intent(in) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
   real(4), intent(in) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
   real(4), intent(in) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
   real(4), intent(in) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
   real(4), intent(in) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
   real(4), intent(in) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
   real(4), intent(in) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
   real(4), intent(in) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
   !real(4), intent(in) :: tlai(:,:,:)      ! tlai(lon, lat, time) total projected LAI
   real(4), intent(in) :: h2osoi(:,:,:,:)  ! h2osoi(lon, lat, levgrnd, time) daily volumetric soil water (mm3/mm3)(vegetated landunits only) 
   real(4), intent(in) :: soilliq(:,:,:,:) ! soilliq(lon, lat, levgrnd, time) daily soil liquid water (kg/m2)(vegetated landunits only) 
   real(4), intent(in) :: tsoi(:,:,:,:)    ! tsoi(lon, lat, levgrnd, time) daily soil temperature (K) (vegetated landunits only) 
   ! Note: pctpft, pctsand, pctclay, omfrac from surface data set are DOUBLE. Use real(8).
   real(8), intent(in) :: pctsand(:,:,:)   ! pctsand(lon, lat, nlevsoi) percent sand
   real(8), intent(in) :: pctclay(:,:,:)   ! pctclay(lon, lat, nlevsoi) percent clay
   real(8), intent(in) :: orgdens(:,:,:)   ! orgdens(lon, lat, nlevsoi) organic matter density (kg/m3)

   !Daily weather and soil variables
   real(4), allocatable, intent(out) :: tminday(:,:,:)      ! tminday(lon, lat, day) minimum daily air temperature (C)
   real(4), allocatable, intent(out) :: tmaxday(:,:,:)      ! tmaxday(lon, lat, day) maximum daily air temperature (C)
   real(4), allocatable, intent(out) :: rainday(:,:,:)      ! rainday(lon, lat, day) daily atmospheric rain (mm/day)
   real(4), allocatable, intent(out) :: snowday(:,:,:)      ! snowday(lon, lat, day) daily atmospheric snow (mm/day)
   real(4), allocatable, intent(out) :: fsdsday(:,:,:)      ! fsdsday(lon, lat, day) 24-hour average atmospheric incident solar radiation (W/m^2)
   real(4), allocatable, intent(out) :: windday(:,:,:)      ! windday(lon, lat, day) 24-hour average atmospheric wind velocity magnitude (m/s)
   !real(4), allocatable, intent(out) :: tlaiday(:,:,:)      ! tlaiday(lon, lat, day) daily total projected LAI
   real(4), allocatable, intent(out) :: gppday(:,:,:)       ! gppday(lon, lat, day) daily GPP (gC/m2/day)
   real(4), allocatable, intent(out) :: soipsiday(:,:,:)    ! soipsiday(lon, lat, day) depth-weighted soil water potential (MPa) (vegetated landunits only) 
   real(4), allocatable, intent(out) :: h2osoiday(:,:,:,:)  ! h2osoiday(lon, lat, nlevsoi, day) average liquid volumetric soil water BY LAYER (mm3/mm3) (vegetated landunits only) 
   real(4), allocatable, intent(out) :: h2ofrznsoiday(:,:,:,:) ! h2ofrznsoiday(lon, lat, nlevsoi, day) average frozen volumetric soil water BY LAYER (mm3/mm3) (vegetated landunits only) 
!  real(4), allocatable, intent(out) :: soilliqday(:,:,:,:) ! soilliqday(lon, lat, nlevsoi, day) average soil liquid water BY LAYER (kg/m2) (vegetated landunits only) 
   real(4), allocatable, intent(out) :: tsoiday(:,:,:,:)    ! tsoiday(lon, lat, nlevsoi, day) soil temperature BY LAYER (K) (vegetated landunits only) 

   ! LOCAL VARIABLES
   integer :: i, j, k, day
   integer ::  calyear, dayofyr, month, dayofmo, hour  ! needed only to compute calyear for test output file
   character*200 :: fileDailyClimOutputCSV      ! Name of output file (.csv)
   character*1   :: c                           ! delimeter (a comma)
   integer :: tbotMissing                       ! set to 1 if missing hourly tbot value, reset to 0 each day
   integer :: rainMissing                       ! set to 1 if missing hourly rain value, reset to 0 each day
   integer :: snowMissing                       ! set to 1 if missing hourly snow value, reset to 0 each day
   integer :: fsdsMissing                       ! set to 1 if missing hourly fsds value, reset to 0 each day
   integer :: windMissing                       ! set to 1 if missing hourly wind value, reset to 0 each day
   !integer :: tlaiMissing                       ! set to 1 if missing hourly tlai value, reset to 0 each day
   integer :: gppMissing                        ! set to 1 if missing hourly fpsn value, reset to 0 each day
   integer, allocatable :: h2osoiMissing(:)     ! BY LAYER, set to 1 if missing daily h2osoi value, reset to 0 each day
   integer, allocatable :: soilliqMissing(:)    ! BY LAYER, set to 1 if missing daily soilliq value, reset to 0 each day
   integer, allocatable :: tsoiMissing(:)       ! BY LAYER, set to 1 if missing daily tsoi value, reset to 0 each day
   character*4   :: xyear                       ! string version of calyear
   integer :: lyr50cm                           ! soil layer with lower boundary ~50cm (0.50 m)
   real(4) :: levlower, levupper, levlower50cm  ! lower and upper boundaries of soil layer with center at levgrnd(i) (m)
   real(4) :: soipsi                            ! Intermediate calculation for soil water potential (MPa)
   real(4) :: sand, clay, org                   ! Temporary variables for soipsi calculation
   real(4), allocatable :: levgrndthickness(:)  ! levgrndthickness(nlevgrnd) thickness of soil layers (m)
  

   if (verbose .ge. 0) print *, "Summarizing daily variables..."

   !----------------------------------------------------------------------------------------------------
   ! Write weather variables to .csv file for testing purposes
   !----------------------------------------------------------------------------------------------------
   calyear = 0
   if (ntimes .gt. 1) call getYearMonthDayHour(time(ntimes-1), calyear, dayofyr, month, dayofmo, hour)
   write(xyear,'(i4)') calyear
   if (doWriteDailyClimateTestFile .eq. 1) then
      fileDailyClimOutputCSV = 'dailyClim_' // xyear // '.csv'
      open(unit=12, file=trim(fileDailyClimOutputCSV))
      c = ','
      !write(unit=12,fmt='(13(a12))') 'lon,', 'lat,', 'year,', 'day,', & 
      !   'tmin(C),', 'tmax(C),', 'rain(mm),', 'snow(mm),', 'fsds(W/m2),', &
      !   'wind(m/s),', 'tlai,', 'gpp(gC/m2)', 'soipsi(MPa)'
      write(unit=12,fmt='(13(a12))') 'lon,', 'lat,', 'year,', 'day,', & 
         'tmin(C),', 'tmax(C),', 'rain(mm),', 'snow(mm),', 'fsds(W/m2),', &
         'wind(m/s),', 'gpp(gC/m2)', 'soipsi(MPa)'
   endif
   !----------------------------------------------------------------------------------------------------

   allocate(tminday(1:nlon, 1:nlat, 1:365))
   allocate(tmaxday(1:nlon, 1:nlat, 1:365))
   allocate(rainday(1:nlon, 1:nlat, 1:365))
   allocate(snowday(1:nlon, 1:nlat, 1:365))
   allocate(fsdsday(1:nlon, 1:nlat, 1:365))
   allocate(windday(1:nlon, 1:nlat, 1:365))
   !allocate(tlaiday(1:nlon, 1:nlat, 1:365))
   allocate(gppday(1:nlon, 1:nlat, 1:365))

   tminday(:,:,:) = 1000
   tmaxday(:,:,:) = -1000
   rainday(:,:,:) = 0
   snowday(:,:,:) = 0
   fsdsday(:,:,:) = 0
   windday(:,:,:) = 0
   !tlaiday(:,:,:) = 0
   gppday(:,:,:) = 0

   !------------------------------------------------------------------------------------------
   ! Compute  soil water potential by layer and  depth-weighted average of soil water potential 
   ! from volumetric soil water content, pctsand, pctclay, and orgdens.
   ! ATTENTION - calculate for top 50 cm soil levels only, use variable from surface data set
   ! Make sure lyr50cm=6 < nlevgrnd=15 and lyr50cm=6 < nlevsoi=10.
   !------------------------------------------------------------------------------------------

   lyr50cm = 6

   !tsoiday, h2osoiday, and h2ofrznsoiday need to be dimensioned for 1:nlevsoi (1:10) for casacnp
   allocate(h2osoiday(1:nlon, 1:nlat, 1:nlevsoi, 1:365))        !This calculation is done by layer
   allocate(h2ofrznsoiday(1:nlon, 1:nlat, 1:nlevsoi, 1:365))       !This calculation is done by layer
!  allocate(soilliqday(1:nlon, 1:nlat, 1:nlevsoi, 1:365))       !This calculation is done by layer
   allocate(tsoiday(1:nlon, 1:nlat, 1:nlevsoi, 1:365))          !This calculation is done by layer
   allocate(soipsiday(1:nlon, 1:nlat, 1:365))                   !This calculation is depth-weighted average 
   allocate(levgrndthickness(1:nlevgrnd))
   allocate(h2osoiMissing(1:nlevsoi))
   allocate(soilliqMissing(1:nlevsoi))
   allocate(tsoiMissing(1:nlevsoi))

   levgrndthickness(:) = 0.0
   h2osoiday(:,:,:,:) = 0
   h2ofrznsoiday(:,:,:,:) = 0
!  soilliqday(:,:,:,:) = 0
   tsoiday(:,:,:,:) = 0
   soipsiday(:,:,:) = 0

   levupper = 0
   levlower50cm = 0.5  ! Depth of soil (m) for the soil layer closest to 50cm. This will be reset below. -mdh 2/21/2017
   !!do k = 1, lyr50cm ! Computer levgrndthickness(k) below 50cm -mdh 2/21/2017
   do k = 1, nlevsoi
      levlower = levgrnd(k) + (levgrnd(k+1) - levgrnd(k))/2
      levgrndthickness(k) = levlower - levupper
      levupper = levlower
      if (verbose .ge. 2) print *, k, levgrndthickness(k)
      if (k == lyr50cm) levlower50cm = levlower
   end do

   !------------------------------------------------------------------------------------------
   ! Compute daily values
   !------------------------------------------------------------------------------------------

   do i = 1, nlat
      do j = 1, nlon
         do day = 1, 365
   
            tbotMissing = 0
            rainMissing = 0
            snowMissing = 0
            fsdsMissing = 0
            windMissing = 0
            !tlaiMissing = 0
            h2osoiMissing(:) = 0
            soilliqMissing(:) = 0
            tsoiMissing(:) = 0
            gppMissing = 0
   
            if (tbot(j,i,day) .lt. MISSING_VALUE) then
               if (tbot(j,i,day) .lt. tminday(j,i,day) + KELVIN) tminday(j,i,day) = tbot(j,i,day) - KELVIN
               if (tbot(j,i,day) .gt. tmaxday(j,i,day) + KELVIN) tmaxday(j,i,day) = tbot(j,i,day) - KELVIN
            else
               tbotMissing = 1
            endif
   
            if (rain(j,i,day) .lt. MISSING_VALUE) then
               rainday(j,i,day) = rain(j,i,day) * SEC_PER_DAY 
            else
               rainMissing = 1
            endif
   
            if (snow(j,i,day) .lt. MISSING_VALUE) then
               snowday(j,i,day) = snow(j,i,day) * SEC_PER_DAY
            else
               snowMissing = 1
            endif
   
            if (fsds(j,i,day) .lt. MISSING_VALUE) then
               fsdsday(j,i,day) = fsds(j,i,day)
            else
               fsdsMissing = 1
            endif
   
            if (wind(j,i,day) .lt. MISSING_VALUE) then
               windday(j,i,day) = wind(j,i,day) 
            else
               windMissing = 1
            endif
   
            !if (tlai(j,i,day) .lt. MISSING_VALUE) then
            !   tlaiday(j,i,day) = tlai(j,i,day)
            !else
            !   tlaiMissing = 1
            !endif
   
            if (fpsn(j,i,day) .lt. MISSING_VALUE) then
               ! Convert umol/m2/sec to gC/m2/day
               gppday(j,i,day) = fpsn(j,i,day) * 12.0 * MOL_PER_UMOL * SEC_PER_DAY
            else
               gppMissing = 1
            endif
   
            do k = 1, nlevsoi 
               !if (h2osoi(j,i,k,day) .lt. MISSING_VALUE) then
               if (soilliq(j,i,k,day) .lt. MISSING_VALUE) then
                  !!h2osoiday(j,i,k,day) = h2osoi(j,i,k,day)
                  ! Compute volumentric soil water content (mm3/mm3) from soil liquid water (kg/m2) -mdh 2/20/2017
                  h2osoiday(j,i,k,day) = soilliq(j,i,k,day) * 0.001 / levgrndthickness(k)
                  ! Compute the frozen fraction also. -mdh 3/13/2017
                  h2ofrznsoiday(j,i,k,day) = h2osoi(j,i,k,day) - h2osoiday(j,i,k,day)
                  if (h2ofrznsoiday(j,i,k,day) < 0.0) then
                     !write(*,*) 'WARNING: h2ofrznsoiday(',j,i,k,day,') < 0:', h2ofrznsoiday(j,i,k,day)
                     !write(*,*) '  Resetting h2ofrznsoiday to 0.0'
                     h2ofrznsoiday(j,i,k,day) = 0.0
                  endif
                  !write(*,*)
                  !write(*,*) 'soilliq(',j,i,k,day, ') =', soilliq(j,i,k,day)
                  !write(*,*) 'h2osoiday(',j,i,k,day, ') =', h2osoiday(j,i,k,day)
                  !write(*,*) 'h2ofrznsoiday(',j,i,k,day, ') =', h2ofrznsoiday(j,i,k,day)
                  !write(*,*) 'levgrndthickness(', k, ') =', levgrndthickness(k)
               else
                  h2osoiMissing(k) = 1
                  soilliqMissing(k) = 1
               end if
            end do
   
            do k = 1, nlevsoi 
               if (tsoi(j,i,k,day) .lt. MISSING_VALUE) then
                  tsoiday(j,i,k,day) = tsoi(j,i,k,day) 
               else
                  tsoiMissing(k) = 1
               end if
            end do
   
            if (tbotMissing .eq. 1) tminday(j,i,day) = MISSING_VALUE
            if (tbotMissing .eq. 1) tmaxday(j,i,day) = MISSING_VALUE
            if (rainMissing .eq. 1) rainday(j,i,day) = MISSING_VALUE
            if (snowMissing .eq. 1) snowday(j,i,day) = MISSING_VALUE
            if (fsdsMissing .eq. 1) fsdsday(j,i,day) = MISSING_VALUE
            if (windMissing .eq. 1) windday(j,i,day) = MISSING_VALUE
            !if (tlaiMissing .eq. 1) tlaiday(j,i,day) = MISSING_VALUE
            if (gppMissing .eq. 1) gppday(j,i,day) = MISSING_VALUE
            do k = 1, nlevsoi
               if (h2osoiMissing(k) .eq. 1) h2osoiday(j,i,k,day) = MISSING_VALUE
               if (soilliqMissing(k) .eq. 1) h2osoiday(j,i,k,day) = MISSING_VALUE
               if (h2osoiMissing(k) .eq. 1 .or. soilliqMissing(k) .eq. 1) h2ofrznsoiday(j,i,k,day) = MISSING_VALUE
               if (tsoiMissing(k) .eq. 1) tsoiday(j,i,k,day) = MISSING_VALUE
            end do
         enddo
      enddo
   enddo

   !------------------------------------------------------------------------------------------
   ! Compute soil water potential by layer and depth-weighted average of soil water  
   ! potential from volumetric soil water content, pctsand, pctclay, and orgdens.
   ! Depth-weighted averages include the top 50 cm only.
   !------------------------------------------------------------------------------------------

   do day = 1, 365
      do k = 1, lyr50cm
         do i = 1, nlat
            do j = 1, nlon
               if (h2osoiday(j,i,k,day) .lt. MISSING_VALUE) then
                  sand = REAL(pctsand(j,i,k))
                  clay = REAL(pctclay(j,i,k))
                  org = REAL(orgdens(j,i,k))
                  ! Compute volumentric soil water content (mm3/mm3) from soil liquid water (kg/m2) -mdh 2/20/2017
                  ! h2osoiday was already computed from soilliq above, so don't need to recalculate? -mdh 2/20/2017
                  ! h2osoiday(j,i,k,day) = soilliq(j,i,k,day) * 0.001 / levgrndthickness(k)
                  ! h2ofrznsoiday(j,i,k,day) = h2osoi(j,i,k,day) - h2osoiday(j,i,k,day)
                  !write(*,*) 'soilliq(',j,i,k,day, ') =', soilliq(j,i,k,day)
                  !write(*,*) 'h2osoiday(',j,i,k,day, ') =', h2osoiday(j,i,k,day)
                  !write(*,*) 'h2ofrznsoiday(',j,i,k,day, ') =', h2ofrznsoiday(j,i,k,day)
                  !write(*,*) 'levgrndthickness(', k, ') =', levgrndthickness(k)
                  soipsi = soilWaterPotential(h2osoiday(j,i,k,day), sand, clay, org, levgrnd(k))
                  soipsiday(j,i,day) = soipsiday(j,i,day) + soipsi * levgrndthickness(k)
               else
                  soipsiday(j,i,day) = MISSING_VALUE
               end if
            end do
         end do
      end do
   end do

   do day = 1, 365
      do i = 1, nlat
         do j = 1, nlon
            if (soipsiday(j,i,day) .lt. MISSING_VALUE) then
               soipsiday(j,i,day) = soipsiday(j,i,day) / levlower50cm
            end if
         end do
      end do
   end do

   !------------------------------------------------------------------------------------------
   ! Write to daily weather to .csv file for testing purposes.
   if (doWriteDailyClimateTestFile .eq. 1) then
      do j = LONINDX,LONINDX
         do i = LATINDX,LATINDX
            do day = 1,365
               write(unit=12, fmt='(2(f0.2,a1), 2(i6,a1), 9(f0.6,a1))') &
                   lon1d(j),c, lat1d(i),c, calyear,c, day,c, &
                   tminday(j,i,day),c, tmaxday(j,i,day),c, rainday(j,i,day),c, snowday(j,i,day),c, &
                   fsdsday(j,i,day),c, windday(j,i,day),c, gppday(j,i,day),c, soipsiday(j,i,day),c
                   !fsdsday(j,i,day),c, windday(j,i,day),c, tlaiday(j,i,day),c, gppday(j,i,day),c, soipsiday(j,i,day),c
            enddo
         enddo
      enddo
      close(unit=12)
   endif

   if (verbose .ge. 0) print *, "Done summarizing daily variables"

   end subroutine summarizeDailyVariables


!-----------------------------------------------------------------------
! This subroutine is commented out because it is not currently being used
! and I didn't want to confuse this code with other code I am modifying.
! Melannie Hartman 2/20/2017
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: summarizeDailyFromHourlyWeather
!
! !INTERFACE:
!  subroutine summarizeDailyFromHourlyWeather(nlon, nlat, ntimes, nlevgrnd, nlevsoi, &
!                                lon1d, lat1d, time, levgrnd, &
!                                flds, fsds, fsdsnd, fsdsni, fsdsvd, fsdsvi, pbot, &
!                                qbot, rain, snow, tbot, thbot, wind, &
!                                tlai, h2osoi, tsoi, pctsand, pctclay, orgdens, &
!                                tminday, tmaxday, rainday, snowday, fsdsday, windday, &
!                                tlaiday, soipsiday, h2osoiday, tsoiday)
!
! !DESCRIPTION:
!  Average and aggregate hourly weather variables to compute daily weather variables.
!  Compute daily depth-weighted soil water potential from h2osoi, pctsand, pctclay, and orgdens.
!
! !USES:
!  none
!
!  !ARGUMENTS:
!  integer, intent(in) :: nlon             ! number of longitudes 
!  integer, intent(in) :: nlat             ! number of latitudes
!  integer, intent(in) :: ntimes           ! total number of hourly time periods (8760 hours - 1 year)
!  integer, intent(in) :: nlevgrnd         ! number of soil levels from history file
!  integer, intent(in) :: nlevsoi          ! number of soil levels from surface data set
!  real(4), intent(in) :: lon1d(:)         ! longitude(lon) (degrees_east)
!  real(4), intent(in) :: lat1d(:)         ! latitude(lat) (degrees_north)
!  real(4), intent(in) :: time(:)          ! days since 1901-01-01 00:00:00
!  real(4), intent(in) :: levgrnd(:)       ! coordinate soil levels (m)
!  ! All the hourly variables are the mean over the timestep 
!  real(4), intent(in) :: flds(:,:,:)      ! flds(lon, lat, time) atmospheric longwave radiation (W/m^2)
!  real(4), intent(in) :: fsds(:,:,:)      ! fsds(lon, lat, time) atmospheric incident solar radiation (W/m^2)
!  real(4), intent(in) :: fsdsnd(:,:,:)    ! fsdsnd(lon, lat, time) direct nir incident solar radiation (W/m^2)
!  real(4), intent(in) :: fsdsni(:,:,:)    ! fsdsni(lon, lat, time) diffuse nir incident solar radiation (W/m^2)
!  real(4), intent(in) :: fsdsvd(:,:,:)    ! fsdsvd(lon, lat, time) direct vis incident solar radiation" (W/m^2)
!  real(4), intent(in) :: fsdsvi(:,:,:)    ! fsdsvi(lon, lat, time) diffuse vis incident solar radiation (W/m^2)
!  real(4), intent(in) :: pbot(:,:,:)      ! pbot(lon, lat, time) atmospheric pressure (Pa)
!  real(4), intent(in) :: qbot(:,:,:)      ! qbot(lon, lat, time) atmospheric specific humidity (kg/kg)
!  real(4), intent(in) :: rain(:,:,:)      ! rain(lon, lat, time) atmospheric rain (mm/sec) 
!  real(4), intent(in) :: snow(:,:,:)      ! snow(lon, lat, time) atmospheric snow (mm/sec) 
!  real(4), intent(in) :: tbot(:,:,:)      ! tbot(lon, lat, time) atmospheric air temperature (K)
!  real(4), intent(in) :: thbot(:,:,:)     ! thbot(lon, lat, time) atmospheric air potential temperature (K)
!  real(4), intent(in) :: wind(:,:,:)      ! wind(lon, lat, time) atmospheric wind velocity magnitude (m/s)
!  real(4), intent(in) :: tlai(:,:,:)      ! tlai(lon, lat, time) total projected LAI
!  real(4), intent(in) :: h2osoi(:,:,:,:)  ! h2osoi(lon, lat, levgrnd, time) hourly volumetric soil water (mm3/mm3)(vegetated landunits only) 
!  real(4), intent(in) :: tsoi(:,:,:,:)    ! tsoi(lon, lat, levgrnd, time) hourly soil temperature (K) (vegetated landunits only) 
!  ! Note: pctpft, pctsand, pctclay, omfrac from surface data set are DOUBLE. Use real(8).
!  real(8), intent(in) :: pctsand(:,:,:)   ! pctsand(lon, lat, nlevsoi) percent sand
!  real(8), intent(in) :: pctclay(:,:,:)   ! pctclay(lon, lat, nlevsoi) percent clay
!  real(8), intent(in) :: orgdens(:,:,:)   ! orgdens(lon, lat, nlevsoi) organic matter density (kg/m3)
!
!  !Daily weather and soil variables
!  real(4), allocatable, intent(out) :: tminday(:,:,:)      ! tminday(lon, lat, day) minimum daily air temperature (C)
!  real(4), allocatable, intent(out) :: tmaxday(:,:,:)      ! tmaxday(lon, lat, day) maximum daily air temperature (C)
!  real(4), allocatable, intent(out) :: rainday(:,:,:)      ! rainday(lon, lat, day) daily atmospheric rain (mm/day)
!  real(4), allocatable, intent(out) :: snowday(:,:,:)      ! snowday(lon, lat, day) daily atmospheric snow (mm/day)
!  real(4), allocatable, intent(out) :: fsdsday(:,:,:)      ! fsdsday(lon, lat, day) 24-hour average atmospheric incident solar radiation (W/m^2)
!  real(4), allocatable, intent(out) :: windday(:,:,:)      ! windday(lon, lat, day) 24-hour average atmospheric wind velocity magnitude (m/s)
!  real(4), allocatable, intent(out) :: tlaiday(:,:,:)      ! tlaiday(lon, lat, day) daily total projected LAI
!  real(4), allocatable, intent(out) :: soipsiday(:,:,:)    ! soipsiday(lon, lat, day) depth-weighted soil water potential (MPa) (vegetated landunits only) 
!  real(4), allocatable, intent(out) :: h2osoiday(:,:,:,:)  ! h2osoiday(lon, lat, nlevsoi, day) average volumetric soil water BY LAYER (mm3/mm3) (vegetated landunits only) 
!  real(4), allocatable, intent(out) :: tsoiday(:,:,:,:)    ! tsoiday(lon, lat, nlevsoi, day) soil temperature BY LAYER (K) (vegetated landunits only) 
!
!  ! LOCAL VARIABLES
!  integer :: i, j, k, tcnt, day
!  integer ::  calyear, dayofyr, month, dayofmo, hour   ! needed only to compute calyear for test output file
!  character*200 :: fileDailyClimOutputCSV      ! Name of output file (.csv)
!  character*1   :: c                           ! delimeter (a comma)
!  integer :: tbotMissing                       ! set to 1 if missing hourly tbot value, reset to 0 each day
!  integer :: rainMissing                       ! set to 1 if missing hourly rain value, reset to 0 each day
!  integer :: snowMissing                       ! set to 1 if missing hourly snow value, reset to 0 each day
!  integer :: fsdsMissing                       ! set to 1 if missing hourly fsds value, reset to 0 each day
!  integer :: windMissing                       ! set to 1 if missing hourly wind value, reset to 0 each day
!  integer :: tlaiMissing                       ! set to 1 if missing hourly tlai value, reset to 0 each day
!  integer, allocatable :: h2osoiMissing(:)     ! BY LAYER, set to 1 if missing hourly h2osoi value, reset to 0 each day
!  integer, allocatable :: tsoiMissing(:)       ! BY LAYER, set to 1 if missing hourly tsoi value, reset to 0 each day
!  character*4   :: xyear                       ! string version of calyear
!  integer :: lyr50cm                           ! soil layer with lower boundary ~50cm (0.50 m)
!  real(4) :: levlower, levupper                ! lower and upper boundaries of soil layer with center at levgrnd(i) (m)
!  real(4) :: soipsi                            ! Intermediate calculation for soil water potential (MPa)
!  real(4) :: sand, clay, org                   ! Temporary variables for soipsi calculation
!  real(4), allocatable :: levgrndthickness(:)  ! levgrndthickness(nlevgrnd) thickness of soil layers (m)
! 
!
!  if (verbose .ge. 0) print *, "Summarizing daily weather from hourly data..."
!
!  !----------------------------------------------------------------------------------------------------
!  ! Write weather variables to .csv file for testing purposes
!  !----------------------------------------------------------------------------------------------------
!  calyear = 0
!  if (ntimes .gt. 1) call getYearMonthDayHour(time(ntimes-1), calyear, dayofyr, month, dayofmo, hour)
!  write(xyear,'(i4)') calyear
!  if (doWriteDailyClimateTestFile .eq. 1) then
!     fileDailyClimOutputCSV = 'dailyClim_' // xyear // '.csv'
!     open(unit=12, file=trim(fileDailyClimOutputCSV))
!     c = ','
!     write(unit=12,fmt='(12(a12))') 'lon,', 'lat,', 'year,', 'day,', & 
!        'tmin(C),', 'tmax(C),', 'rain(mm),', 'snow(mm),', 'fsds(W/m2),', &
!        'wind(m/s),', 'tlai,', 'soipsi(MPa)'
!  endif
!  !----------------------------------------------------------------------------------------------------
!
!  allocate(tminday(1:nlon, 1:nlat, 1:365))
!  allocate(tmaxday(1:nlon, 1:nlat, 1:365))
!  allocate(rainday(1:nlon, 1:nlat, 1:365))
!  allocate(snowday(1:nlon, 1:nlat, 1:365))
!  allocate(fsdsday(1:nlon, 1:nlat, 1:365))
!  allocate(windday(1:nlon, 1:nlat, 1:365))
!  allocate(tlaiday(1:nlon, 1:nlat, 1:365))
!
!  tminday(:,:,:) = 1000
!  tmaxday(:,:,:) = -1000
!  rainday(:,:,:) = 0
!  snowday(:,:,:) = 0
!  fsdsday(:,:,:) = 0
!  windday(:,:,:) = 0
!  tlaiday(:,:,:) = 0
!
!  !------------------------------------------------------------------------------------------
!  ! Compute  soil water potential by layer and  depth-weighted average of soil water potential 
!  ! from volumetric soil water content, pctsand, pctclay, and orgdens.
!  ! ATTENTION - calculate for top 50 cm soil levels only, use variable from surface data set
!  ! Make sure lyr50cm=6 < nlevgrnd=15 and lyr50cm=6 < nlevsoi=10.
!  !------------------------------------------------------------------------------------------
!
!  lyr50cm = 6
!
!  !tsoiday and h2osoiday need to be dimensioned for 1:nlevsoi (1:10) for casacnp
!  allocate(h2osoiday(1:nlon, 1:nlat, 1:nlevsoi, 1:365))        !This calculation is done by layer
!  allocate(tsoiday(1:nlon, 1:nlat, 1:nlevsoi, 1:365))          !This calculation is done by layer
!  allocate(soipsiday(1:nlon, 1:nlat, 1:365))                   !This calculation is depth-weighted average 
!  allocate(levgrndthickness(1:nlevgrnd))
!  allocate(h2osoiMissing(1:nlevsoi))
!  allocate(tsoiMissing(1:nlevsoi))
!
!  levgrndthickness(:) = 0.0
!  h2osoiday(:,:,:,:) = 0
!  tsoiday(:,:,:,:) = 0
!  soipsiday(:,:,:) = 0
!
!  levupper = 0
!  do k = 1, lyr50cm
!     levlower = levgrnd(k) + (levgrnd(k+1) - levgrnd(k))/2
!     levgrndthickness(k) = levlower - levupper
!     levupper = levlower
!     if (verbose .ge. 2) print *, k, levgrndthickness(k)
!  end do
!
!  !------------------------------------------------------------------------------------------
!  ! Compute daily values
!  !------------------------------------------------------------------------------------------
!
!  do i = 1, nlat
!     do j = 1, nlon
!
!        tcnt = 1
!        day = 1
!        tbotMissing = 0
!        rainMissing = 0
!        snowMissing = 0
!        fsdsMissing = 0
!        windMissing = 0
!        tlaiMissing = 0
!        h2osoiMissing(:) = 0
!        tsoiMissing(:) = 0
!
!        !If any hourly data for a variable is invalid, the daily value for that variable is also invalid
!        do while (tcnt .le. ntimes .AND. day .le. 365) 
!
!           if (tbot(j,i,tcnt) .lt. MISSING_VALUE) then
!              if (tbot(j,i,tcnt) .lt. tminday(j,i,day) + KELVIN) tminday(j,i,day) = tbot(j,i,tcnt) - KELVIN
!              if (tbot(j,i,tcnt) .gt. tmaxday(j,i,day) + KELVIN) tmaxday(j,i,day) = tbot(j,i,tcnt) - KELVIN
!           else
!              tbotMissing = 1
!           endif
!
!           if (rain(j,i,tcnt) .lt. MISSING_VALUE) then
!              rainday(j,i,day) = rainday(j,i,day) + rain(j,i,tcnt) * SEC_PER_HOUR 
!           else
!              rainMissing = 1
!           endif
!
!           if (snow(j,i,tcnt) .lt. MISSING_VALUE) then
!              snowday(j,i,day) = snowday(j,i,day) + snow(j,i,tcnt) * SEC_PER_HOUR
!           else
!              snowMissing = 1
!           endif
!
!           if (fsds(j,i,tcnt) .lt. MISSING_VALUE) then
!              fsdsday(j,i,day) = fsdsday(j,i,day) + fsds(j,i,tcnt) * 0.041666667       ! 0.041666667 = (1.0/HOURS_PER_DAY)
!           else
!              fsdsMissing = 1
!           endif
!
!           if (wind(j,i,tcnt) .lt. MISSING_VALUE) then
!              windday(j,i,day) = windday(j,i,day) + wind(j,i,tcnt) * 0.041666667       ! 0.041666667 = (1.0/HOURS_PER_DAY)
!           else
!              windMissing = 1
!           endif
!
!           if (tlai(j,i,tcnt) .lt. MISSING_VALUE) then
!              tlaiday(j,i,day) = tlaiday(j,i,day) + tlai(j,i,tcnt) * 0.041666667       ! 0.041666667 = (1.0/HOURS_PER_DAY)
!           else
!              tlaiMissing = 1
!           endif
!
!           do k = 1, nlevsoi 
!              if (h2osoi(j,i,k,tcnt) .lt. MISSING_VALUE) then
!                 h2osoiday(j,i,k,day) = h2osoiday(j,i,k,day) + h2osoi(j,i,k,tcnt) * 0.041666667    ! 0.041666667 = (1.0/HOURS_PER_DAY)
!              else
!                 h2osoiMissing(k) = 1
!              end if
!           end do
!
!           do k = 1, nlevsoi 
!              if (tsoi(j,i,k,tcnt) .lt. MISSING_VALUE) then
!                 tsoiday(j,i,k,day) = tsoiday(j,i,k,day) + tsoi(j,i,k,tcnt) * 0.041666667      ! 0.041666667 = (1.0/HOURS_PER_DAY)
!              else
!                 tsoiMissing(k) = 1
!              end if
!           end do
!
!           if (MODULO(tcnt,24) .eq. 0) then
!              !The last hour for the day has been read
!              if (tbotMissing .eq. 1) tminday(j,i,day) = MISSING_VALUE
!              if (tbotMissing .eq. 1) tmaxday(j,i,day) = MISSING_VALUE
!              if (rainMissing .eq. 1) rainday(j,i,day) = MISSING_VALUE
!              if (snowMissing .eq. 1) snowday(j,i,day) = MISSING_VALUE
!              if (fsdsMissing .eq. 1) fsdsday(j,i,day) = MISSING_VALUE
!              if (windMissing .eq. 1) windday(j,i,day) = MISSING_VALUE
!              if (tlaiMissing .eq. 1) tlaiday(j,i,day) = MISSING_VALUE
!              do k = 1, nlevsoi
!                 if (h2osoiMissing(k) .eq. 1) h2osoiday(j,i,k,day) = MISSING_VALUE
!              end do
!              do k = 1, nlevsoi
!                 if (tsoiMissing(k) .eq. 1) tsoiday(j,i,k,day) = MISSING_VALUE
!              end do
!
!              ! Reinitialize for the next day
!              day = day + 1
!              if (day .le. 365) then
!                 tbotMissing = 0
!                 rainMissing = 0
!                 snowMissing = 0
!                 fsdsMissing = 0
!                 windMissing = 0
!                 tlaiMissing = 0
!                 h2osoiMissing(:) = 0
!                 tsoiMissing(:) = 0
!              endif
!           endif
!           tcnt = tcnt + 1
!
!         enddo
!     enddo
!  enddo
!
!
!  !------------------------------------------------------------------------------------------
!  ! Compute soil water potential by layer and depth-weighted average of soil water  
!  ! potential from volumetric soil water content, pctsand, pctclay, and orgdens.
!  ! Depth-weighted averages include the top 50 cm only.
!  !------------------------------------------------------------------------------------------
!
!  do day = 1, 365
!     do k = 1, lyr50cm
!        do i = 1, nlat
!           do j = 1, nlon
!              if (h2osoiday(j,i,k,day) .lt. MISSING_VALUE) then
!                 sand = REAL(pctsand(j,i,k))
!                 clay = REAL(pctclay(j,i,k))
!                 org = REAL(orgdens(j,i,k))
!                 soipsi = soilWaterPotential(h2osoiday(j,i,k,day), sand, clay, org, levgrnd(k))
!                 soipsiday(j,i,day) = soipsiday(j,i,day) + soipsi * levgrndthickness(k)
!              else
!                 soipsiday(j,i,day) = MISSING_VALUE
!              end if
!           end do
!        end do
!     end do
!  end do
!
!  do day = 1, 365
!     do i = 1, nlat
!        do j = 1, nlon
!           if (soipsiday(j,i,day) .lt. MISSING_VALUE) then
!              soipsiday(j,i,day) = soipsiday(j,i,day) / levlower
!           end if
!        end do
!     end do
!  end do
!
!  !------------------------------------------------------------------------------------------
!  ! Write to daily weather to .csv file for testing purposes.
!  if (doWriteDailyClimateTestFile .eq. 1) then
!     do j = LONINDX,LONINDX
!        do i = LATINDX,LATINDX
!           do day = 1,365
!              write(unit=12, fmt='(2(f0.2,a1), 2(i6,a1), 8(f0.6,a1))') &
!                  lon1d(j),c, lat1d(i),c, calyear,c, day,c, &
!                  tminday(j,i,day),c, tmaxday(j,i,day),c, rainday(j,i,day),c, snowday(j,i,day),c, &
!                  fsdsday(j,i,day),c, windday(j,i,day),c, tlaiday(j,i,day),c, soipsiday(j,i,day),c
!           enddo
!        enddo
!     enddo
!     close(unit=12)
!  endif
!
!  if (verbose .ge. 0) print *, "Done summarizing daily weather from hourly data"
!
!  end subroutine summarizeDailyFromHourlyWeather

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilWaterPotential
!
! !INTERFACE:
   real function soilWaterPotential(h2osoi, sand, clay, org, zsoi)
!
! !DESCRIPTION:
!  Computes soil water potential (MPa).
!  Reference: See CLM4 Technical Notes, eq. 7.85, p. 141 (some modifications have been made)
!
! !USES:
!  none
!
!  !ARGUMENTS:
   real(4), intent(in) :: h2osoi    ! volumetric soil water content (mm3/mm3)
   real(4), intent(in) :: sand      ! percent sand (0-100), from surface dataset
   real(4), intent(in) :: clay      ! percent clay (0-100), from surface dataset
   real(4), intent(in) :: org       ! organic matter density (kg/m3), from surface dataset
   real(4), intent(in) :: zsoi      ! soil level (m), from surface dataset

!  !LOCAL VARIABLES:

   ! Mineral soil:
   real(4) :: Watsat           ! volumetric soil water at saturation (porosity)
   real(4) :: bsw              ! Clapp and Hornberger "b"
   real(4) :: sucsat           ! Saturated soil water suction (mm)

   ! Organic soil:
   real(4) :: om_watsat        ! porosity of organic soil
   real(4) :: om_b             ! Clapp Hornberger paramater for organic soil
   real(4) :: om_sucsat        ! saturated suction for organic matter (mm)
   real(4) :: zsapric          ! depth (m) that organic matter takes on
   real(4) :: om_frac          ! organic matter fraction

   ! Soil water potential:
   real(4) :: head             ! Head of pressure (MPa/m)
   real(4) :: s_node           ! water content relative to saturation
   real(4) :: smp_mm           ! soil water potential (mm)
   real(4) :: smp_mpa          ! soil water potential (MPa)

   Watsat = 0.489 - 0.00126 * sand
   bsw    = 2.91 + 0.159 * clay
   sucsat = 10.0 * ( 10.0**(1.88 - 0.0131 * sand) )

   !Characteristics of sapric peat
   zsapric    = 0.5 
   om_watsat  = max(0.93 - 0.1 * (zsoi/zsapric), 0.83)
   om_b       = min(2.7  + 9.3 * (zsoi/zsapric), 12.0)
   om_sucsat  = min(10.3 - 0.2 * (zsoi/zsapric), 10.1)

   !This is the combined mineral and organic values:
  
   om_frac = org / ORGANIC_MAX
   Watsat  = (1.0 - om_frac)*Watsat + om_watsat*om_frac
   bsw     = (1.0 - om_frac)*bsw    + om_b*om_frac
   sucsat  = (1.0 - om_frac)*sucsat + om_sucsat*om_frac

   s_node = h2osoi / Watsat
   smp_mm = -sucsat * s_node**(-bsw)

   !Convert soil water potential from mm to MPa

   head = DENH2O * GRAV * 1.e-06        ! Head of pressure  (MPa/m)
   smp_mpa = smp_mm * 1.e-03 * head     ! mm -> m -> MPa

   soilWaterPotential = smp_mpa

   end function soilWaterPotential

end module clm_netcdf_tools
