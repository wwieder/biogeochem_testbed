module casacnp_model

!-----------------------------------------------------------------------
!BOP
!
!  !MODULE: casacnp_model
!
!  !DESCRIPTION:
!   Create gridded input files needed to run the casacnp_clm model.
!   The casacnp_clm model is CASACNP with some of the I/O changed so that
!   CLM derived meteorology and soil properties can be used to drive the model.
!   The casacnp_clm model is a separate executable and is not run here.
!
! !USES:
   use clm_common

! !OUTPUT text files
!  unit=201     gridinfo_igbpz.csv (required by casacnp_clm model)
!  unit=202     gridinfo_soil.csv (required by casacnp_clm model)
!  unit=203     gridinfo_pft.csv (optional, for testing purposes)
!  unit=204      met.txt (optional, for testing purposes)
!
! !OUTPUT NetCDF files
!  met.nc (required by casacnp_clm model)

! !CONSTANTS
   implicit none
!  include 'netcdf.inc'
   integer, parameter :: XSOIL = 6              ! number of soil layers used by casacnp model
!  real(4), parameter :: KELVIN = 273.15        ! Equivalent of 0 degrees Celsius
   integer, parameter :: seed = 86456
   integer, parameter :: MISSING_INT_SM = -99
   real(4), parameter :: MISSING_FLT_SM = -99.9

! !PUBLIC MEMBER FUNCTIONS:
   public :: readCLMgridIGBPGrid
   public :: writeCASACNPfiles
   public :: createCASACNPMetNcFile
   public :: writeCASACNPMetNcFile
   public :: createCASACNPMetNcFilePt
   public :: writeCASACNPMetNcFilePt

! !REVISION HISTORY:
!  Created by Melannie Hartman
!  melannie@ucar.edu
!  December 2013 - February 2014
!  Updated to create single point files - Aug. 15, 2016
!  Updated to write daily N deposition to met.nc files - Aug. 22, 2016
!  Remove NPP and LAI from met.nc files -mdh 11/6/2017
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: createCASACNPMetNcFile
!
! !INTERFACE:
   subroutine createCASACNPMetNcFile(ncfilename,iyear1, iyear2, ndays, nlon, nlat, &
                                     nsoilyrs, lon1d, lat1d, landfrac, cellMissing, cellid)

! !DESCRIPTION:
! ! Define the dimensions and variables in the met.nc file that contains the multiple
! ! years of daily meteorological inputs and soil temperature and moisture conditions 
! ! for the CASACNP model. Write missing values to all arrays.  These arrays will be 
! ! filled later by calling subroutine writeCASACNPMetNcFile.
!
! !USES:
!
! !REVISION HISTORY:
!
! !CALLED BY: subroutine writeCASACNPfiles
!
! !ARGUMENTS:
   character(len=*), intent(in) :: ncfilename           ! name of output ncfile
   integer, intent(in) :: iyear1, iyear2                ! first and last calendar years
   integer, intent(in) :: ndays                         ! number of days in the year
   integer, intent(in) :: nlon                          ! number of longitudes 
   integer, intent(in) :: nlat                          ! number of latitudes
   integer, intent(in) :: nsoilyrs                      ! number of soil layers
   real(4), intent(in) :: lon1d(:)                      ! lat(lat) (degrees_east)
   real(4), intent(in) :: lat1d(:)                      ! lat(lat) (degrees_north)
   real(4), intent(in) :: landfrac(:,:)                 ! landfrac(nlon,nlat) land fraction (0-1)
   integer, intent(in) :: cellMissing(:,:)              ! cellMissing(nlon,nlat) 0=no missing data, 1=missing data
   integer, intent(in) :: cellid(:,:)                   ! cellid(nlon,nlat) grid cell IDs (1..nlat*nlon)


! !LOCAL VARIABLES:
   integer :: iyr
   integer :: status
   integer :: ncid                              ! NetCDF file ID
   integer :: dimid_lon, dimid_lat              ! NetCDF dimension IDs for latitude and longitude
   integer :: dimid_myear, dimid_time           ! NetCDF dimension IDs for time variables
   integer :: dimid_xsoil                       ! NetCDF dimention ID for soil layers
   integer :: dims(4)                           ! Array of NetCDF dimension IDs for defining variables
   integer :: nyears                            ! number of years of data
   integer :: ntimes                            ! total number of timesteps (nyears * ndays)
   integer :: varid_lon, varid_lat              ! NetCDF variable ID for latitude and longitude
   integer :: varid_landfrac                    ! NetCDF variable ID for land fraction
   integer :: varid_year                        ! NetCDF variable ID for year (calendar years)
   integer :: varid_xtairk                      ! NetCDF variable ID for air temperature
   integer :: varid_ndep                        ! NetCDF variable ID for N deposition
   integer :: varid_xtsoil                      ! NetCDF variable ID for soil temperature 
   integer :: varid_xmoist, varid_xfrznmoist    ! NetCDF variable ID for soil moisture
   !integer :: varid_xcnpp,                     ! NetCDF variable ID for NPP 
   !integer :: varid_xlai                       ! NetCDF variable ID for LAI
   integer :: varid_xcgpp                       ! NetCDF variable ID for GPP
   integer :: varid_mask                        ! NetCDF variable ID for cellMissing(nlon,nlat)
   integer :: varid_cellid                      ! NetCDF variable ID for cellid(nlon,nlat)
   integer, allocatable :: year(:)              ! year(nyears) calendar years
   !real(4), allocatable :: xcnpp(:,:,:)        ! xcnpp(nlon,nlat,ntimes) daily NPP (gC/m2/day)
   !real(4), allocatable :: xlai(:,:,:)         ! xlai(nlon,nlat,ntimes) daily LAI (m2/m2)
   real(4), allocatable :: xcgpp(:,:,:)         ! xcgpp(nlon,nlat,ntimes) daily GPP (gC/m2/day)
   real(4), allocatable :: xtairk(:,:,:)        ! xtairk(nlon,nlat,ntimes) daily average air temperature (K) 
   real(4), allocatable :: ndepDay(:,:,:)       ! ndepDay(nlonpt,nlatpt,ntimes) daily N deposition (gN/m2/day)
   real(4), allocatable :: xtsoil(:,:,:,:)      ! xtsoil(nlon,nlat,XSOIL,ntimes) daily soil temperature (K)
   real(4), allocatable :: xmoist(:,:,:,:)      ! xmoist(nlon,nlat,XSOIL,ntimes) daily volumetric soil liquid moisture (m3/m3)
   real(4), allocatable :: xfrznmoist(:,:,:,:)  ! xfrznmoist(nlon,nlat,XSOIL,ntimes) daily volumetric soil frozen moisture (m3/m3)
   character*100 :: attr_name                   ! String for assigning global and variable attributes
   character*100 :: attr_units                  ! String for assigning global and variable attributes
   character*10 :: date_string                  ! String for assigning date to global attributes
   character*8 :: time_string                   ! String for assigning time to global attributes

   if (verbose .ge. 0) print *, 'Creating NetCDF to file ', trim(ncfilename)

   nyears = iyear2 - iyear1 + 1
   ntimes = ndays * nyears

   if (verbose .ge. 0) then
      print *, "  createCASACNPMetNcFile: nyears: ", nyears
      print *, "  createCASACNPMetNcFile: ndays: ", ndays
      print *, "  createCASACNPMetNcFile: nlon: ", nlon
      print *, "  createCASACNPMetNcFile: nlat: ", nlat
      print *, "  createCASACNPMetNcFile: nsoilyrs: ", nsoilyrs
      print *, "  createCASACNPMetNcFile: ntimes: ", ntimes
   endif

   allocate(year(nyears))
   !allocate(xcnpp(nlon,nlat,ntimes))
   !allocate(xlai(nlon,nlat,ntimes))
   allocate(xcgpp(nlon,nlat,ntimes))
   allocate(xtairk(nlon,nlat,ntimes))
   allocate(ndepDay(nlon,nlat,ntimes))
   allocate(xtsoil(nlon,nlat,nsoilyrs,ntimes))
   allocate(xmoist(nlon,nlat,nsoilyrs,ntimes))
   allocate(xfrznmoist(nlon,nlat,nsoilyrs,ntimes))

   !Create an array of calendar years to write to the NetCDF file
   do iyr = 1, nyears
       year(iyr) = iyear1 + iyr - 1
   enddo
   
   !xcnpp(:,:,:) = MISSING_VALUE 
   !xlai(:,:,:) = MISSING_VALUE 
   xcgpp(:,:,:) = MISSING_VALUE 
   xtairk(:,:,:) = MISSING_VALUE 
   ndepDay(:,:,:) = MISSING_VALUE 
   xtsoil(:,:,:,:) = MISSING_VALUE 
   xmoist(:,:,:,:) = MISSING_VALUE 
   xfrznmoist(:,:,:,:) = MISSING_VALUE 

   dims(1) = 0
   dims(2) = 0
   dims(3) = 0
   dims(4) = 0

   ! Create the netcdf file 

   status = nf_create(ncfilename, NF_CLOBBER, ncid)
   if (status /= nf_noerr) call handle_err(status, "ncfilename")


   ! Define file dimensions

   status = nf_def_dim(ncid, 'lon', nlon, dimid_lon)
   if (status /= nf_noerr) call handle_err(status, "lon")

   status = nf_def_dim(ncid, 'lat', nlat, dimid_lat)
   if (status /= nf_noerr) call handle_err(status, "lat")

   status = nf_def_dim(ncid, 'time', ntimes, dimid_time)
   if (status /= nf_noerr) call handle_err(status, "time")

   status = nf_def_dim(ncid, 'nsoilyrs', nsoilyrs, dimid_xsoil)
   if (status /= nf_noerr) call handle_err(status, "nsoilyrs")

   status = nf_def_dim(ncid, 'myear', nyears, dimid_myear)
   if (status /= nf_noerr) call handle_err(status, "myear")


   ! Define variables

   dims(1) = dimid_lon
   status = nf_def_var(ncid, 'lon', NF_REAL, 1, dims, varid_lon)
   if (status /= nf_noerr) call handle_err(status, "def_var(lon)")

   dims(1) = dimid_lat
   status = nf_def_var(ncid, 'lat', NF_REAL, 1, dims, varid_lat)
   if (status /= nf_noerr) call handle_err(status, "def_var(lat)")

   dims(1) = dimid_myear
   status = nf_def_var(ncid, 'year', NF_INT, 1, dims, varid_year)
   if (status /= nf_noerr) call handle_err(status, "def_var(year)")


   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions are in the opposite 
   ! order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat

   status = nf_def_var(ncid, 'landfrac', NF_REAL, 2, dims, varid_landfrac)
   if (status /= nf_noerr) call handle_err(status, "landfrac")

   status = nf_def_var(ncid, 'cellMissing', NF_INT, 2, dims, varid_mask)
   if (status /= nf_noerr) call handle_err(status, "cellMissing")

   status = nf_def_var(ncid, 'cellid', NF_INT, 2, dims, varid_cellid)
   if (status /= nf_noerr) call handle_err(status, "cellid")


   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions are in the opposite 
   ! order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat
   dims(3) = dimid_time

   status = nf_def_var(ncid, 'xtairk', NF_REAL, 3, dims, varid_xtairk)
   if (status /= nf_noerr) call handle_err(status, "xtairk")

   status = nf_def_var(ncid, 'ndep', NF_REAL, 3, dims, varid_ndep)
   if (status /= nf_noerr) call handle_err(status, "ndep")

   !status = nf_def_var(ncid, 'xlai', NF_REAL, 3, dims, varid_xlai)
   !if (status /= nf_noerr) call handle_err(status, "xlai")

   !status = nf_def_var(ncid, 'xcnpp', NF_REAL, 3, dims, varid_xcnpp)
   !if (status /= nf_noerr) call handle_err(status, "xcnpp")

   status = nf_def_var(ncid, 'xcgpp', NF_REAL, 3, dims, varid_xcgpp)
   if (status /= nf_noerr) call handle_err(status, "xcgpp")
 

   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions are in the opposite 
   ! order that they appear in the NetCDF file with ncdump. 

   dims(1) = dimid_lon
   dims(2) = dimid_lat
   dims(3) = dimid_xsoil
   dims(4) = dimid_time
 
   status = nf_def_var(ncid, 'xtsoil', NF_REAL, 4, dims, varid_xtsoil)
   if (status /= nf_noerr) call handle_err(status, "def_var(xtsoil)")
 
   status = nf_def_var(ncid, 'xmoist', NF_REAL, 4, dims, varid_xmoist)
   if (status /= nf_noerr) call handle_err(status, "def_var(xmoist)")
 
   status = nf_def_var(ncid, 'xfrznmoist', NF_REAL, 4, dims, varid_xfrznmoist)
   if (status /= nf_noerr) call handle_err(status, "def_var(xfrznmoist)")
 
   ! Global attributes
   attr_name = 'CLM History file information for CASACNP model'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'title', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "title")
 
   attr_name = 'NOTE: None of the variables are weighted by land fraction!'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'comment', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "comment")
 
   call get_time_and_date(date_string, time_string)
   attr_name = 'created on ' // date_string // ' ' // time_string
   status = nf_put_att_text(ncid, NF_GLOBAL, 'history', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "history")
 
   attr_name = 'CLM Model'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'source', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "source")
 
 
   ! Attributes of the variables

   ! Attributes of year variable
   attr_name = 'calendar years'
   status = nf_put_att_text(ncid, varid_year, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   ! Attributes of lon variable
   attr_name = 'coordinate longitude'
   attr_units = 'degrees_east'
   call PutVariableAttributeReal(ncid, varid_lon, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of lat variable
   attr_name = 'coordinate latitude'
   attr_units = 'degrees_north'
   call PutVariableAttributeReal(ncid, varid_lat, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of landfrac variable
   attr_name = 'land fraction from pft dataset'
   attr_units = 'unitless'
   call PutVariableAttributeReal(ncid, varid_landfrac, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cellMissing variable
   attr_name = 'Missing Data Mask'
   status = nf_put_att_text(ncid, varid_mask, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   attr_units = '0=no missing data, 1=missing data'
   status = nf_put_att_text(ncid, varid_mask, 'units', len(trim(attr_units)), trim(attr_units))

   ! Attributes of cellid variable
   attr_name = 'Grid Cell ID'
   status = nf_put_att_text(ncid, varid_cellid, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   attr_units = '1..nlat*nlon'
   status = nf_put_att_text(ncid, varid_cellid, 'units', len(trim(attr_units)), trim(attr_units))

   ! Attributes of xcnpp variable
   !attr_name = 'net primary production'
   !attr_units = 'gC m-2 day-1'
   !call PutVariableAttributeReal(ncid, varid_xcnpp, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of xcgpp variable
   attr_name = 'gross primary production'
   attr_units = 'gC m-2 day-1'
   call PutVariableAttributeReal(ncid, varid_xcgpp, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xlai variable
   !attr_name = 'Leaf Area Index'
   !attr_units = 'm2/m2'
   !call PutVariableAttributeReal(ncid, varid_xlai, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of xtairk variable
   attr_name = 'average daily air temperature'
   attr_units = 'K'
   call PutVariableAttributeReal(ncid, varid_xtairk, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of ndepDay variable
   attr_name = 'daily N deposition derived from annual N deposition'
   attr_units = 'gN/m2/day'
   call PutVariableAttributeReal(ncid, varid_ndep, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xtsoil variable
   attr_name = 'average daily soil temperature by layer'
   attr_units = 'K'
   call PutVariableAttributeReal(ncid, varid_xtsoil, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xmoist variable
   attr_name = 'volumetric soil liquid water content by layer'
   attr_units = 'm3/m3'
   call PutVariableAttributeReal(ncid, varid_xmoist, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xfrznmoist variable
   attr_name = 'volumetric soil frozen water content by layer'
   attr_units = 'm3/m3'
   call PutVariableAttributeReal(ncid, varid_xfrznmoist, attr_name, attr_units, MISSING_VALUE)

   !--------------------------------------------------------------------------------------
   ! End the definition phase so that variables can be written to the file

   status = nf_enddef(ncid)
   if (status /= nf_noerr) call handle_err(status, "enddef")
   !--------------------------------------------------------------------------------------

   ! Write variable values to ncfilename

   status =  nf_put_var(ncid, varid_year, year)
   if (status /= nf_noerr) call handle_err(status, "put_var(year)")

   status =  nf_put_var(ncid, varid_lon, lon1d)
   if (status /= nf_noerr) call handle_err(status, "put_var(lon)")

   status =  nf_put_var(ncid, varid_lat, lat1d)
   if (status /= nf_noerr) call handle_err(status, "put_var(lat)")

   status =  nf_put_var(ncid, varid_landfrac, landfrac)
   if (status /= nf_noerr) call handle_err(status, "put_var(landfrac)")

   status =  nf_put_var(ncid, varid_mask, cellMissing)
   if (status /= nf_noerr) call handle_err(status, "put_var(cellMissing)")

   status =  nf_put_var(ncid, varid_cellid, cellid)
   if (status /= nf_noerr) call handle_err(status, "put_var(cellid)")

   status =  nf_put_var(ncid, varid_xtairk, xtairk)
   if (status /= nf_noerr) call handle_err(status, "put_var(xtairk)")

   status =  nf_put_var(ncid, varid_ndep, ndepDay)
   if (status /= nf_noerr) call handle_err(status, "put_var(ndepDay)")

   !status =  nf_put_var(ncid, varid_xlai, xlai)
   !if (status /= nf_noerr) call handle_err(status, "put_var(xlai)")

   !status =  nf_put_var(ncid, varid_xcnpp, xcnpp)
   !if (status /= nf_noerr) call handle_err(status, "put_var(xcnpp)")

   status =  nf_put_var(ncid, varid_xcgpp, xcgpp)
   if (status /= nf_noerr) call handle_err(status, "put_var(xcgpp)")

   status =  nf_put_var(ncid, varid_xtsoil, xtsoil)
   if (status /= nf_noerr) call handle_err(status, "put_var(xtsoil)")

   status =  nf_put_var(ncid, varid_xmoist, xmoist)
   if (status /= nf_noerr) call handle_err(status, "put_var(xmoist)")

   status =  nf_put_var(ncid, varid_xfrznmoist, xfrznmoist)
   if (status /= nf_noerr) call handle_err(status, "put_var(xfrznmoist)")


   status = nf_close(ncid)

   !deallocate(year,xcnpp,xcgpp,xlai,xtairk,ndepDay,xtsoil,xmoist,xfrznmoist)
   deallocate(year,xcgpp,xtairk,ndepDay,xtsoil,xmoist,xfrznmoist)

   if (verbose .ge. 0) print *, 'Done creating file ', trim(ncfilename)

end subroutine createCASACNPMetNcFile


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: writeCASACNPMetNcFile
!
! !INTERFACE:
!  subroutine writeCASACNPMetNcFile(ncFileMet, calyear, iyear1, iyear2, nlon, nlat, &
!      ndays, nsoilyrs, xcnpp, xcgpp, xlai, xtairk, xtsoil, xmoist, xfrznmoist, ndepDay)
   subroutine writeCASACNPMetNcFile(ncFileMet, calyear, iyear1, iyear2, nlon, nlat, &
       ndays, nsoilyrs, xcgpp, xtairk, xtsoil, xmoist, xfrznmoist, ndepDay)

! !DESCRIPTION:
!  Write a year's worth of daily values to met.nc file.
!
! !USES:
!
! !CALLED BY: subroutine writeCASACNPfiles
!
   !ARGUMENTS
    character(len=*), intent(in) :: ncFileMet   ! name of output ncfile
    integer, intent(in) :: calyear              ! current calendar year
    integer, intent(in) :: iyear1, iyear2       ! first and last calendar years in the file
    integer, intent(in) :: nlon                 ! number of longitudes 
    integer, intent(in) :: nlat                 ! number of latitudes
    integer, intent(in) :: ndays                ! number of days in the year
    integer, intent(in) :: nsoilyrs             ! number of soil layers
    !real(4), intent(in) :: xlai(:,:,:)         ! xlai(nlon,nlat,ndays) total projected LAI
    !real(4), intent(in) :: xcnpp(:,:,:)        ! xcnpp(nlon,nlat,ndays) Net Primary Production (gC/m2/day)
    real(4), intent(in) :: xcgpp(:,:,:)         ! xcgpp(nlon,nlat,ndays) Gross Primary Production (gC/m2/day)
    real(4), intent(in) :: xtairk(:,:,:)        ! xtairk(nlon,nlat,ndays) daily average air temperature (K) 
    real(4), intent(in) :: ndepDay(:,:,:)       ! ndepDay*nlon,nlat,ndays) daily N deposition (gN/m2/day)
    real(4), intent(in) :: xtsoil(:,:,:,:)      ! xtsoil(nlon,nlat,XSOIL,ndays), CASACNP soil temperature by layer (K)
    real(4), intent(in) :: xmoist(:,:,:,:)      ! xmoist(nlon,nlat,XSOIL,ndays), CASACNP volumetric soil liquid water content by layer 
    real(4), intent(in) :: xfrznmoist(:,:,:,:)  ! xfrznmoist(nlon,nlat,XSOIL,ndays), CASACNP volumetric soil frozen water content by layer 
   
   !LOCAL VARIABLES
    integer :: status
    integer :: ncid
    integer :: yrCnt                            ! Year of output relative to the starting year (0..nyears-1)
    !integer :: xcnpp_varid, xcgpp_varid, xlai_varid
    integer :: xcgpp_varid 
    integer :: ndep_varid, xtairk_varid, xtsoil_varid, xmoist_varid, xfrznmoist_varid
    integer :: start3(3), count3(3)             ! start and count arrays for reading 3-D data from netcdf files
    integer :: start4(4), count4(4)             ! start and count arrays for reading 4-D data from netcdf files


    if (verbose .ge. 0) print *, 'writeCASACNPMetNcFile: Updating file ', trim(ncFileMet), ' for year ', calyear
    
    yrCnt = calyear - iyear1

    ! Open the netcdf file for writing

    status = nf_open(ncFileMet, NF_WRITE, ncid)
    if (status /= nf_noerr) call handle_err(status, "ncFileMet")

    !Retrieve variable IDs
 
    !status = nf_inq_varid(ncid, "xcnpp", xcnpp_varid)
    !if (status /= nf_noerr) call handle_err(status, "inq_varid(xcnpp)")

    status = nf_inq_varid(ncid, "xcgpp", xcgpp_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xcgpp)")
 
    !status = nf_inq_varid(ncid, "xlai", xlai_varid)
    !if (status /= nf_noerr) call handle_err(status, "inq_varid(xlai)")
  
    status = nf_inq_varid(ncid, "xtairk", xtairk_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xtairk)")

    status = nf_inq_varid(ncid, "ndep", ndep_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(ndep)")

    status = nf_inq_varid(ncid, "xtsoil", xtsoil_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xtsoil)")

    status = nf_inq_varid(ncid, "xmoist", xmoist_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xmoist)")

    status = nf_inq_varid(ncid, "xfrznmoist", xfrznmoist_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xfrznmoist)")


    ! Write the 3-D variables.  Dimensions are in the opposite order they will
    ! appear in the NetCDF file because of Column Major Order FORTRAN arrays.

    start3(1) = 1
    start3(2) = 1
    start3(3) = yrCnt * ndays + 1

    !count3 = (/ ndays, nlat, nlon /)   ! netCDF file dimension order
    count3 = (/ nlon, nlat, ndays /)    ! FORTRAN dimension order

    !status = nf_put_vara(ncid, xcnpp_varid, start3, count3, xcnpp)
    !if (status /= nf_noerr) call handle_err(status, "put_vara(xcnpp)")

    status = nf_put_vara(ncid, xcgpp_varid, start3, count3, xcgpp)
    if (status /= nf_noerr) call handle_err(status, "put_vara(xcgpp)")

    !status = nf_put_vara(ncid, xlai_varid, start3, count3, xlai)
    !if (status /= nf_noerr) call handle_err(status, "put_vara(xlai)")

    status = nf_put_vara(ncid, xtairk_varid, start3, count3, xtairk)
    if (status /= nf_noerr) call handle_err(status, "put_vara(xtairk)")

    status = nf_put_vara(ncid, ndep_varid, start3, count3, ndepDay)
    if (status /= nf_noerr) call handle_err(status, "put_vara(ndepDay)")


    ! Write the 4-D variables. Dimensions are in the opposite order they will
    ! appear in the NetCDF file because of Column Major Order FORTRAN arrays.

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = yrCnt * ndays + 1

    !count4 = (/ ndays, nsoilyrs, nlat, nlon /) ! netCDF file dimension order
    count4 = (/ nlon, nlat, nsoilyrs, ndays /)  ! FORTRAN dimension order

    status = nf_put_vara(ncid, xtsoil_varid, start4, count4, xtsoil)
    if (status /= nf_noerr) call handle_err(status, "put_vara(xtsoil)")

    status = nf_put_vara(ncid, xmoist_varid, start4, count4, xmoist)
    if (status /= nf_noerr) call handle_err(status, "put_vara(xmoist)")

    status = nf_put_vara(ncid, xfrznmoist_varid, start4, count4, xfrznmoist)
    if (status /= nf_noerr) call handle_err(status, "put_vara(xfrznmoist)")


    status = nf_close(ncid)

    if (verbose .ge. 0) print *, 'Done updating file ', trim(ncFileMet), ' for year ', calyear

    end subroutine writeCASACNPMetNcFile




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: createCASACNPMetNcFilePt
!
! !INTERFACE:
   subroutine createCASACNPMetNcFilePt(ncfilename,iyear1, iyear2, ndays, ilon, ilat, nlon, nlat, &
                                     nsoilyrs, lon1d, lat1d, landfrac, cellMissing, cellid)

! !DESCRIPTION:
! ! Define the dimensions and variables in the met.nc file that contains the multiple
! ! years of daily meteorological inputs and soil temperature and moisture conditions 
! ! for the CASACNP model. Write missing values to all arrays.  These arrays will be 
! ! filled later by calling subroutine writeCASACNPMetNcFilePt.
!
! !USES:
!
! !REVISION HISTORY:
!
! !CALLED BY: subroutine writeCASACNPfiles
!
! !ARGUMENTS:
   character(len=*), intent(in) :: ncfilename           ! name of output ncfile
   integer, intent(in) :: iyear1, iyear2                ! first and last calendar years
   integer, intent(in) :: ndays                         ! number of days in the year
   integer, intent(in) :: ilon                          ! longitude index 
   integer, intent(in) :: ilat                          ! latitude index
   integer, intent(in) :: nlon                          ! number of longitudes in the entire grid 
   integer, intent(in) :: nlat                          ! number of latitudes in the entire grid
   integer, intent(in) :: nsoilyrs                      ! number of soil layers
   real(4), intent(in) :: lon1d(:)                      ! lat(lat) (degrees_east)
   real(4), intent(in) :: lat1d(:)                      ! lat(lat) (degrees_north)
   real(4), intent(in) :: landfrac(:,:)                 ! landfrac(nlon,nlat) land fraction (0-1)
   integer, intent(in) :: cellMissing(:,:)              ! cellMissing(nlon,nlat) 0=no missing data, 1=missing data
   integer, intent(in) :: cellid(:,:)                   ! cellid(nlon,nlat) grid cell IDs (1..nlat*nlon)


! !LOCAL VARIABLES:
   integer :: iyr
   integer :: nlonpt, nlatpt
   integer :: status
   integer :: ncid                              ! NetCDF file ID
   integer :: dimid_lon, dimid_lat              ! NetCDF dimension IDs for latitude and longitude
   integer :: dimid_myear, dimid_time           ! NetCDF dimension IDs for time variables
   integer :: dimid_xsoil                       ! NetCDF dimention ID for soil layers
   integer :: dims(4)                           ! Array of NetCDF dimension IDs for defining variables
   integer :: nyears                            ! number of years of data
   integer :: ntimes                            ! total number of timesteps (nyears * ndays)
   integer :: varid_lon, varid_lat              ! NetCDF variable ID for latitude and longitude
   integer :: varid_landfrac                    ! NetCDF variable ID for land fraction
   integer :: varid_year                        ! NetCDF variable ID for year (calendar years)
   integer :: varid_xtairk                      ! NetCDF variable ID for air temperature
   integer :: varid_ndep                        ! NetCDF variable ID for N deposition
   integer :: varid_xtsoil                      ! NetCDF variable ID for soil temperature 
   integer :: varid_xmoist, varid_xfrznmoist    ! NetCDF variable ID for soil moisture
   !integer :: varid_xcnpp                      ! NetCDF variable ID for NPP 
   !integer :: varid_xlai                       ! NetCDF variable ID for LAI
   integer :: varid_xcgpp                       ! NetCDF variable ID for GPP
   integer :: varid_mask                        ! NetCDF variable ID for cellMissing(nlonpt,nlatpt)
   integer :: varid_cellid                      ! NetCDF variable ID for cellid(nlonpt,nlatpt)
   integer, allocatable :: year(:)              ! year(nyears) calendar years
   !real(4), allocatable :: xcnpp(:,:,:)        ! xcnpp(nlonpt,nlatpt,ntimes) daily NPP (gC/m2/day)
   !real(4), allocatable :: xlai(:,:,:)         ! xlai(nlonpt,nlatpt,ntimes) daily LAI (m2/m2)
   real(4), allocatable :: xcgpp(:,:,:)         ! xcgpp(nlonpt,nlatpt,ntimes) daily GPP (gC/m2/day)
   real(4), allocatable :: xtairk(:,:,:)        ! xtairk(nlonpt,nlatpt,ntimes) daily average air temperature (K) 
   real(4), allocatable :: ndepDay(:,:,:)       ! ndepDay(nlonpt,nlatpt,ntimes) daily N deposition (gN/m2/day)
   real(4), allocatable :: xtsoil(:,:,:,:)      ! xtsoil(nlonpt,nlatpt,XSOIL,ntimes) daily soil temperature (K)
   real(4), allocatable :: xmoist(:,:,:,:)      ! xmoist(nlonpt,nlatpt,XSOIL,ntimes) daily volumetric liquid soil moisture (m3/m3)
   real(4), allocatable :: xfrznmoist(:,:,:,:)  ! xfrznmoist(nlonpt,nlatpt,XSOIL,ntimes) daily volumetric frozen soil moisture (m3/m3)
   character*100 :: attr_name                   ! String for assigning global and variable attributes
   character*100 :: attr_units                  ! String for assigning global and variable attributes
   character*10 :: date_string                  ! String for assigning date to global attributes
   character*8 :: time_string                   ! String for assigning time to global attributes

   if (verbose .ge. 0) print *, 'Creating NetCDF to file ', trim(ncfilename)

   nyears = iyear2 - iyear1 + 1
   ntimes = ndays * nyears
   nlonpt = 1
   nlatpt = 1

   if (verbose .ge. 0) then
      print *, "  createCASACNPMetNcFilePt: nyears: ", nyears
      print *, "  createCASACNPMetNcFilePt: ndays: ", ndays
      print *, "  createCASACNPMetNcFilePt: ilon: ", ilon
      print *, "  createCASACNPMetNcFilePt: ilat: ", ilat
      print *, "  createCASACNPMetNcFilePt: nlon: ", nlon
      print *, "  createCASACNPMetNcFilePt: nlat: ", nlat
      print *, "  createCASACNPMetNcFilePt: nsoilyrs: ", nsoilyrs
      print *, "  createCASACNPMetNcFilePt: ntimes: ", ntimes
   endif

   allocate(year(nyears))
   !allocate(xcnpp(nlonpt,nlatpt,ntimes))
   !allocate(xlai(nlonpt,nlatpt,ntimes))
   allocate(xcgpp(nlonpt,nlatpt,ntimes))
   allocate(xtairk(nlonpt,nlatpt,ntimes))
   allocate(ndepDay(nlonpt,nlatpt,ntimes))
   allocate(xtsoil(nlonpt,nlatpt,nsoilyrs,ntimes))
   allocate(xmoist(nlonpt,nlatpt,nsoilyrs,ntimes))
   allocate(xfrznmoist(nlonpt,nlatpt,nsoilyrs,ntimes))

   !Create an array of calendar years to write to the NetCDF file
   do iyr = 1, nyears
       year(iyr) = iyear1 + iyr - 1
   enddo
   
   !xcnpp(:,:,:) = MISSING_VALUE 
   !xlai(:,:,:) = MISSING_VALUE 
   xcgpp(:,:,:) = MISSING_VALUE 
   xtairk(:,:,:) = MISSING_VALUE 
   ndepDay(:,:,:) = MISSING_VALUE 
   xtsoil(:,:,:,:) = MISSING_VALUE 
   xmoist(:,:,:,:) = MISSING_VALUE 
   xfrznmoist(:,:,:,:) = MISSING_VALUE 

   dims(1) = 0
   dims(2) = 0
   dims(3) = 0
   dims(4) = 0

   ! Create the netcdf file 

   status = nf_create(ncfilename, NF_CLOBBER, ncid)
   if (status /= nf_noerr) call handle_err(status, "ncfilename")


   ! Define file dimensions

   status = nf_def_dim(ncid, 'lon', nlonpt, dimid_lon)
   if (status /= nf_noerr) call handle_err(status, "lon")

   status = nf_def_dim(ncid, 'lat', nlatpt, dimid_lat)
   if (status /= nf_noerr) call handle_err(status, "lat")

   status = nf_def_dim(ncid, 'time', ntimes, dimid_time)
   if (status /= nf_noerr) call handle_err(status, "time")

   status = nf_def_dim(ncid, 'nsoilyrs', nsoilyrs, dimid_xsoil)
   if (status /= nf_noerr) call handle_err(status, "nsoilyrs")

   status = nf_def_dim(ncid, 'myear', nyears, dimid_myear)
   if (status /= nf_noerr) call handle_err(status, "myear")


   ! Define variables

   dims(1) = dimid_lon
   status = nf_def_var(ncid, 'lon', NF_REAL, 1, dims, varid_lon)
   if (status /= nf_noerr) call handle_err(status, "def_var(lon)")

   dims(1) = dimid_lat
   status = nf_def_var(ncid, 'lat', NF_REAL, 1, dims, varid_lat)
   if (status /= nf_noerr) call handle_err(status, "def_var(lat)")

   dims(1) = dimid_myear
   status = nf_def_var(ncid, 'year', NF_INT, 1, dims, varid_year)
   if (status /= nf_noerr) call handle_err(status, "def_var(year)")


   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions are in the opposite 
   ! order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat

   status = nf_def_var(ncid, 'landfrac', NF_REAL, 2, dims, varid_landfrac)
   if (status /= nf_noerr) call handle_err(status, "landfrac")

   status = nf_def_var(ncid, 'cellMissing', NF_INT, 2, dims, varid_mask)
   if (status /= nf_noerr) call handle_err(status, "cellMissing")

   status = nf_def_var(ncid, 'cellid', NF_INT, 2, dims, varid_cellid)
   if (status /= nf_noerr) call handle_err(status, "cellid")


   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions are in the opposite 
   ! order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat
   dims(3) = dimid_time

   status = nf_def_var(ncid, 'xtairk', NF_REAL, 3, dims, varid_xtairk)
   if (status /= nf_noerr) call handle_err(status, "xtairk")

   status = nf_def_var(ncid, 'ndep', NF_REAL, 3, dims, varid_ndep)
   if (status /= nf_noerr) call handle_err(status, "ndep")

   !status = nf_def_var(ncid, 'xlai', NF_REAL, 3, dims, varid_xlai)
   !if (status /= nf_noerr) call handle_err(status, "xlai")

   !status = nf_def_var(ncid, 'xcnpp', NF_REAL, 3, dims, varid_xcnpp)
   !if (status /= nf_noerr) call handle_err(status, "xcnpp")

   status = nf_def_var(ncid, 'xcgpp', NF_REAL, 3, dims, varid_xcgpp)
   if (status /= nf_noerr) call handle_err(status, "xcgpp")
 

   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions are in the opposite 
   ! order that they appear in the NetCDF file with ncdump. 

   dims(1) = dimid_lon
   dims(2) = dimid_lat
   dims(3) = dimid_xsoil
   dims(4) = dimid_time
 
   status = nf_def_var(ncid, 'xtsoil', NF_REAL, 4, dims, varid_xtsoil)
   if (status /= nf_noerr) call handle_err(status, "def_var(xtsoil)")
 
   status = nf_def_var(ncid, 'xmoist', NF_REAL, 4, dims, varid_xmoist)
   if (status /= nf_noerr) call handle_err(status, "def_var(xmoist)")
 
   status = nf_def_var(ncid, 'xfrznmoist', NF_REAL, 4, dims, varid_xfrznmoist)
   if (status /= nf_noerr) call handle_err(status, "def_var(xfrznmoist)")
 
 
   ! Global attributes
   attr_name = 'CLM History file information for CASACNP model'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'title', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "title")
 
   attr_name = 'NOTE: None of the variables are weighted by land fraction!'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'comment', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "comment")
 
   call get_time_and_date(date_string, time_string)
   attr_name = 'created on ' // date_string // ' ' // time_string
   status = nf_put_att_text(ncid, NF_GLOBAL, 'history', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "history")
 
   attr_name = 'CLM Model'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'source', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "source")
 
 
   ! Attributes of the variables

   ! Attributes of year variable
   attr_name = 'calendar years'
   status = nf_put_att_text(ncid, varid_year, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   ! Attributes of lon variable
   attr_name = 'coordinate longitude'
   attr_units = 'degrees_east'
   call PutVariableAttributeReal(ncid, varid_lon, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of lat variable
   attr_name = 'coordinate latitude'
   attr_units = 'degrees_north'
   call PutVariableAttributeReal(ncid, varid_lat, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of landfrac variable
   attr_name = 'land fraction from pft dataset'
   attr_units = 'unitless'
   call PutVariableAttributeReal(ncid, varid_landfrac, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cellMissing variable
   attr_name = 'Missing Data Mask'
   status = nf_put_att_text(ncid, varid_mask, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   attr_units = '0=no missing data, 1=missing data'
   status = nf_put_att_text(ncid, varid_mask, 'units', len(trim(attr_units)), trim(attr_units))

   ! Attributes of cellid variable
   attr_name = 'Grid Cell ID'
   status = nf_put_att_text(ncid, varid_cellid, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   attr_units = '1..nlat*nlon'
   status = nf_put_att_text(ncid, varid_cellid, 'units', len(trim(attr_units)), trim(attr_units))

   ! Attributes of xcnpp variable
   !attr_name = 'net primary production'
   !attr_units = 'gC m-2 day-1'
   !call PutVariableAttributeReal(ncid, varid_xcnpp, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of xcgpp variable
   attr_name = 'gross primary production'
   attr_units = 'gC m-2 day-1'
   call PutVariableAttributeReal(ncid, varid_xcgpp, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xlai variable
   !attr_name = 'Leaf Area Index'
   !attr_units = 'm2/m2'
   !call PutVariableAttributeReal(ncid, varid_xlai, attr_name, attr_units, MISSING_VALUE)
 
   ! Attributes of xtairk variable
   attr_name = 'average daily air temperature'
   attr_units = 'K'
   call PutVariableAttributeReal(ncid, varid_xtairk, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of ndep variable
   attr_name = 'daily N deposition derived from annual N deposition'
   attr_units = 'gN/m2/day'
   call PutVariableAttributeReal(ncid, varid_ndep, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xtsoil variable
   attr_name = 'average daily soil temperature by layer'
   attr_units = 'K'
   call PutVariableAttributeReal(ncid, varid_xtsoil, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xmoist variable
   attr_name = 'volumetric soil liquid water content by layer'
   attr_units = 'm3/m3'
   call PutVariableAttributeReal(ncid, varid_xmoist, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of xfrznmoist variable
   attr_name = 'volumetric soil frozen water content by layer'
   attr_units = 'm3/m3'
   call PutVariableAttributeReal(ncid, varid_xfrznmoist, attr_name, attr_units, MISSING_VALUE)

   !--------------------------------------------------------------------------------------
   ! End the definition phase so that variables can be written to the file

   status = nf_enddef(ncid)
   if (status /= nf_noerr) call handle_err(status, "enddef")
   !--------------------------------------------------------------------------------------

   ! Write variable values to ncfilename

   status =  nf_put_var(ncid, varid_year, year)
   if (status /= nf_noerr) call handle_err(status, "put_var(year)")

   status =  nf_put_var(ncid, varid_lon, lon1d(ilon))
   if (status /= nf_noerr) call handle_err(status, "put_var(lon1d(ilon))")

   status =  nf_put_var(ncid, varid_lat, lat1d(ilat))
   if (status /= nf_noerr) call handle_err(status, "put_var(lat1d(ilat))")

   status =  nf_put_var(ncid, varid_landfrac, landfrac(ilon,ilat))
   if (status /= nf_noerr) call handle_err(status, "put_var(landfrac(ilon,ilat))")

   status =  nf_put_var(ncid, varid_mask, cellMissing(ilon,ilat))
   if (status /= nf_noerr) call handle_err(status, "put_var(cellMissing)")

   status =  nf_put_var(ncid, varid_cellid, cellid(ilon,ilat))
   if (status /= nf_noerr) call handle_err(status, "put_var(cellid)")

   status =  nf_put_var(ncid, varid_xtairk, xtairk)
   if (status /= nf_noerr) call handle_err(status, "put_var(xtairk)")

   status =  nf_put_var(ncid, varid_ndep, ndepDay)
   if (status /= nf_noerr) call handle_err(status, "put_var(ndepDay)")

   !status =  nf_put_var(ncid, varid_xlai, xlai)
   !if (status /= nf_noerr) call handle_err(status, "put_var(xlai)")

   !status =  nf_put_var(ncid, varid_xcnpp, xcnpp)
   !if (status /= nf_noerr) call handle_err(status, "put_var(xcnpp)")

   status =  nf_put_var(ncid, varid_xcgpp, xcgpp)
   if (status /= nf_noerr) call handle_err(status, "put_var(xcgpp)")

   status =  nf_put_var(ncid, varid_xtsoil, xtsoil)
   if (status /= nf_noerr) call handle_err(status, "put_var(xtsoil)")

   status =  nf_put_var(ncid, varid_xmoist, xmoist)
   if (status /= nf_noerr) call handle_err(status, "put_var(xmoist)")

   status =  nf_put_var(ncid, varid_xfrznmoist, xfrznmoist)
   if (status /= nf_noerr) call handle_err(status, "put_var(xfrznmoist)")


   status = nf_close(ncid)

   !deallocate(year,xcnpp,xcgpp,xlai,xtairk,ndepDay,xtsoil,xmoist,xfrznmoist)
   deallocate(year,xcgpp,xtairk,ndepDay,xtsoil,xmoist,xfrznmoist)

   if (verbose .ge. 0) print *, 'Done creating file ', trim(ncfilename)

end subroutine createCASACNPMetNcFilePt


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: writeCASACNPMetNcFilePt
!
! !INTERFACE:
!  subroutine writeCASACNPMetNcFilePt(ncFileMet, calyear, iyear1, iyear2, ilon, ilat, nlon, nlat, &
!      ndays, nsoilyrs, xcnpp, xcgpp, xlai, xtairk, xtsoil, xmoist, xfrznmoist, ndepDay)
   subroutine writeCASACNPMetNcFilePt(ncFileMet, calyear, iyear1, iyear2, ilon, ilat, nlon, nlat, &
       ndays, nsoilyrs, xcgpp, xtairk, xtsoil, xmoist, xfrznmoist, ndepDay)

! !DESCRIPTION:
!  Write a year's worth of daily values to met.nc file.
!
! !USES:
!
! !CALLED BY: subroutine writeCASACNPfiles
!
   !ARGUMENTS
    character(len=*), intent(in) :: ncFileMet   ! name of output ncfile
    integer, intent(in) :: calyear              ! current calendar year
    integer, intent(in) :: iyear1, iyear2       ! first and last calendar years in the file
    integer, intent(in) :: ilon                 ! longitude index 
    integer, intent(in) :: ilat                 ! latitude index
    integer, intent(in) :: nlon                 ! number of longitudes 
    integer, intent(in) :: nlat                 ! number of latitudes
    integer, intent(in) :: ndays                ! number of days in the year
    integer, intent(in) :: nsoilyrs             ! number of soil layers
    !real(4), intent(in) :: xlai(:,:,:)          ! xlai(nlon,nlat,ndays) total projected LAI
    !real(4), intent(in) :: xcnpp(:,:,:)        ! xcnpp(nlon,nlat,ndays) Net Primary Production (gC/m2/day)
    real(4), intent(in) :: xcgpp(:,:,:)         ! xcgpp(nlon,nlat,ndays) Gross Primary Production (gC/m2/day)
    real(4), intent(in) :: xtairk(:,:,:)        ! xtairk(nlon,nlat,ndays) daily average air temperature (K) 
    real(4), intent(in) :: ndepDay(:,:,:)       ! ndepDay(nlon,nlat,ndays) daily average N deposition (gN/m2/day)
    real(4), intent(in) :: xtsoil(:,:,:,:)      ! xtsoil(nlon,nlat,XSOIL,ndays), CASACNP soil temperature by layer (K)
    real(4), intent(in) :: xmoist(:,:,:,:)      ! xmoist(nlon,nlat,XSOIL,ndays), CASACNP volumetric soil liquid water content by layer 
    real(4), intent(in) :: xfrznmoist(:,:,:,:)  ! xfrznmoist(nlon,nlat,XSOIL,ndays), CASACNP volumetric soil frozen water content by layer 
   
   !LOCAL VARIABLES
    integer :: nlonpt, nlatpt
    integer :: status
    integer :: ncid
    integer :: yrCnt                            ! Year of output relative to the starting year (0..nyears-1)
    !integer :: xcnpp_varid, xcgpp_varid, xlai_varid
    integer ::  xcgpp_varid
    integer :: ndep_varid, xtairk_varid, xtsoil_varid, xmoist_varid, xfrznmoist_varid
    integer :: start3(3), count3(3)             ! start and count arrays for reading 3-D data from netcdf files
    integer :: start4(4), count4(4)             ! start and count arrays for reading 4-D data from netcdf files


    if (verbose .ge. 0) print *, 'writeCASACNPMetNcFilePt: Updating file ', trim(ncFileMet), ' for year ', calyear
    
    yrCnt = calyear - iyear1
    nlonpt = 1
    nlatpt = 1

    ! Open the netcdf file for writing

    status = nf_open(ncFileMet, NF_WRITE, ncid)
    if (status /= nf_noerr) call handle_err(status, "ncFileMet")

    !Retrieve variable IDs
 
    !status = nf_inq_varid(ncid, "xcnpp", xcnpp_varid)
    !if (status /= nf_noerr) call handle_err(status, "inq_varid(xcnpp)")

    status = nf_inq_varid(ncid, "xcgpp", xcgpp_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xcgpp)")
 
    !status = nf_inq_varid(ncid, "xlai", xlai_varid)
    !if (status /= nf_noerr) call handle_err(status, "inq_varid(xlai)")
  
    status = nf_inq_varid(ncid, "xtairk", xtairk_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xtairk)")

    status = nf_inq_varid(ncid, "ndep", ndep_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(ndep)")

    status = nf_inq_varid(ncid, "xtsoil", xtsoil_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xtsoil)")

    status = nf_inq_varid(ncid, "xmoist", xmoist_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xmoist)")

    status = nf_inq_varid(ncid, "xfrznmoist", xfrznmoist_varid)
    if (status /= nf_noerr) call handle_err(status, "inq_varid(xfrznmoist)")


    ! Write the 3-D variables.  Dimensions are in the opposite order they will
    ! appear in the NetCDF file because of Column Major Order FORTRAN arrays.

    start3(1) = 1
    start3(2) = 1 
    start3(3) = yrCnt * ndays + 1

    !count3 = (/ ndays, nlatpt, nlonpt /)       ! netCDF file dimension order
    count3 = (/ nlonpt, nlatpt, ndays /)        ! FORTRAN dimension order

    !status = nf_put_vara(ncid, xcnpp_varid, start3, count3, xcnpp(ilon,ilat,:))
    !if (status /= nf_noerr) call handle_err(status, "put_vara(xcnpp(ilon,ilat,:))")

    status = nf_put_vara(ncid, xcgpp_varid, start3, count3, xcgpp(ilon,ilat,:))
    if (status /= nf_noerr) call handle_err(status, "put_vara(xcgpp(ilon,ilat,:))")

    !status = nf_put_vara(ncid, xlai_varid, start3, count3, xlai(ilon,ilat,:))
    !if (status /= nf_noerr) call handle_err(status, "put_vara(xlai(ilon,ilat,:))")

    status = nf_put_vara(ncid, xtairk_varid, start3, count3, xtairk(ilon,ilat,:))
    if (status /= nf_noerr) call handle_err(status, "put_vara(xtairk(ilon,ilat,:))")

    status = nf_put_vara(ncid, ndep_varid, start3, count3, ndepDay(ilon,ilat,:))
    if (status /= nf_noerr) call handle_err(status, "put_vara(ndepDay(ilon,ilat,:))")


    ! Write the 4-D variables. Dimensions are in the opposite order they will
    ! appear in the NetCDF file because of Column Major Order FORTRAN arrays.

    start4(1) = 1
    start4(2) = 1
    start4(3) = 1
    start4(4) = yrCnt * ndays + 1

    !count4 = (/ ndays, nsoilyrs, nlatpt, nlonpt /)     ! netCDF file dimension order
    count4 = (/ nlonpt, nlatpt, nsoilyrs, ndays /)      ! FORTRAN dimension order

    status = nf_put_vara(ncid, xtsoil_varid, start4, count4, xtsoil(ilon,ilat,:,:))
    if (status /= nf_noerr) call handle_err(status, "put_vara(xtsoil(ilon,ilat,:,:))")

    status = nf_put_vara(ncid, xmoist_varid, start4, count4, xmoist(ilon,ilat,:,:))
    if (status /= nf_noerr) call handle_err(status, "put_vara(xmoist(ilon,ilat,:,:))")

    status = nf_put_vara(ncid, xfrznmoist_varid, start4, count4, xfrznmoist(ilon,ilat,:,:))
    if (status /= nf_noerr) call handle_err(status, "put_vara(xfrznmoist(ilon,ilat,:,:))")


    status = nf_close(ncid)

    if (verbose .ge. 0) print *, 'Done updating file ', trim(ncFileMet), ' for year ', calyear

    end subroutine writeCASACNPMetNcFilePt


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readCLMgridIGBPGrid
!
! !INTERFACE:
   subroutine readCLMgridIGBPGrid(fileCLMgridIGBP, nlon, nlat, lon1d, lat1d, vegIGBP, wrtcell, cellid, wrtSingleCells)

! !DESCRIPTION:
!  Read the IGBP veg classes based on CLM grid from text file fileCLMgridIGBP 
!  (clmGrid_IGBP.csv).  This file should have 13825 lines (the first line is
!  a column header).  The points should be sorted first by LONGXY is ascending
!  order then by LATIXY is ascending order.  
!
! !USES:
!
! !CALLED BY: subroutine writeCASACNPfiles
!
!  !ARGUMENTS:
   character(len=*), intent(in) :: fileCLMgridIGBP      ! clmGrid_IGBP.csv (IGBP veg types on CLM grid). Input file.
   integer, intent(in) :: nlon                          ! number of longitudes 
   integer, intent(in) :: nlat                          ! number of latitudes
   real(4), intent(in) :: lon1d(:)                      ! longitude(lon) (degrees_east)
   real(4), intent(in) :: lat1d(:)                      ! latitude(lat) (degrees_north)
   integer, intent(out) :: vegIGBP(:,:)                 ! vegIGBP(nlon,nlat) IGBP veg classification
   integer, intent(out) :: wrtcell(:,:)                 ! wrtcell(nlon,nlat) save individual cell results
   integer, intent(out) :: cellid(:,:)                  ! cellid(nlon,nlat) grid cell ID
   integer, intent(out) :: wrtSingleCells               ! =1 if write output to single cell files instead of to gridded files

  !LOCAL VARIABLES
   integer :: IOstatus
   integer :: i,j, cellCnt
   integer :: lsmlon, lsmlat
   real(4) :: landfrac, longitude, latitude

   if (verbose .ge. 0) print *, 'Reading vegetation file ', trim(fileCLMgridIGBP), '...'
   wrtSingleCells = 0

   open(unit=205,file=trim(fileCLMgridIGBP))
   read(205,*)          !Read past column header
   do i = 1, nlon
      do j = 1, nlat

         read(205,*,IOSTAT=IOstatus) lsmlon, lsmlat, landfrac, longitude, latitude, vegIGBP(i,j), cellCnt, wrtcell(i,j)
         write(*,206) lsmlon, lsmlat, landfrac, longitude, latitude, vegIGBP(i,j), cellid(i,j), wrtcell(i,j)
         206 format(i4, ',', i4, ',', f4.2, ',', f6.2, ',', f6.2, ',', i4, ',', i6, ',', i2, ',', a20)

         if (cellid(i,j) /= cellCnt) then
            write(*,*) 'Error in readCLMgridIGBPGrid'
            write(*,*) 'cellids in ', trim(fileCLMgridIGBP), ' are not specified correctly'
            write(*,*) 'cellid = ', cellid(i,j), ' cellCnt = ', cellCnt
            STOP
         endif

         if (wrtcell(i,j) == 1) then
               wrtSingleCells = 1
         endif

         if (IOstatus < 0) then
            print *, "Unexpected EOF while reading file ", TRIM(fileCLMgridIGBP)
            STOP
         elseif (IOstatus > 0) then
            print *, "An error occurred while reading file ", TRIM(fileCLMgridIGBP)
            STOP
         endif
         if (abs(longitude - lon1d(i)) .gt. 0.05) then
             print *, TRIM(fileCLMgridIGBP), ": Expected longitude ", lon1d(i), " got longitude ", longitude
             STOP
         endif
         if (abs(latitude - lat1d(j)) .gt. 0.05) then
             print *, TRIM(fileCLMgridIGBP), ": Expected latitude ", lat1d(j), " got latitude ", latitude
             STOP
         endif
      enddo
   enddo

   close(205)
   if (verbose .ge. 0) print *, 'Done reading vegetation file ', trim(fileCLMgridIGBP)


   end subroutine readCLMgridIGBPGrid


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: writeCASACNPfiles
!
! !INTERFACE:
!  subroutine writeCASACNPfiles(calyear, iyear1, iyear2, nlon, nlat, npft, &
!                        nlevgrnd, nlevsoi, ndays, ndtime, ndep_calyear, &
!                        lon1d, lat1d, levgrnd, area, landfrac, landmask, pftmask, &
!                        ndepAn, tminday, tmaxday, tlaiday, gpp, &
!                        h2osoiday, h2ofrznsoiday, tsoiday, pctpft, pctsand, pctclay, watsat, watfc)
   subroutine writeCASACNPfiles(calyear, iyear1, iyear2, nlon, nlat, npft, &
                         nlevgrnd, nlevsoi, ndays, ndtime, ndep_calyear, &
                         lon1d, lat1d, levgrnd, area, landfrac, landmask, pftmask, &
                         ndepAn, tminday, tmaxday, gpp, &
                         h2osoiday, h2ofrznsoiday, tsoiday, pctpft, pctsand, pctclay, watsat, watfc)
!! WW added

! !DESCRIPTION:
!  This subroutine is called once a year to create input files for the CASACNP model.  
!  For all years write a year's worth of daily values (GPP, air temperature, 
!  soil temperature, and soil moisture) to the met.nc file.
!  Create grid files for CASACNP model the first year only (calyear = iyear1).
!
! !CALLS:
!    call createCASACNPMetNcFile - create and define grid met.nc file
!    call createCASACNPMetNcFilePt - create and define point met.nc file
!    call readCLMgridIGBPGrid - read the IGBP veg classes based on CLM grida(fileCLMgridIGBP = "clmGrid_IGBP.csv")
!    call hydraulicProperties - determine soil hydraulic properties from soil texture
!    call writeCASACNPMetNcFile - write to grid met.nc file
!    call writeCASACNPMetNcFilePt - write to point met.nc file
!
! !USES:
!
!  !ARGUMENTS:
   integer, intent(in) :: calyear               ! current calendar year (YYYY)
   integer, intent(in) :: iyear1, iyear2        ! first and last calendar years (YYYY)
   integer, intent(in) :: nlon                  ! number of longitudes 
   integer, intent(in) :: nlat                  ! number of latitudes
   integer, intent(in) :: npft                  ! number of PFTs
   integer, intent(in) :: nlevgrnd              ! number of soil levels from history file (~15)
   integer, intent(in) :: nlevsoi               ! number of soil levels from surface data set (~10)
   integer, intent(in) :: ndays                 ! total number of daily time periods (usually 365)
   integer, intent(in) :: ndtime                ! number of years of N deposition
   integer, intent(in) :: ndep_calyear(:)       ! calendar years in N deposition file
   real(4), intent(in) :: lon1d(:)              ! longitude(lon) (degrees_east)
   real(4), intent(in) :: lat1d(:)              ! latitude(lat) (degrees_north)
   real(4), intent(in) :: levgrnd(:)            ! levgrnd(nlevgrnd), coordinate soil levels (m)
   real(4), intent(in) :: area(:,:)             ! grid cell areas (km^2)
   real(4), intent(in) :: landfrac(:,:)         ! landfrac(nlon,nlat) land fraction
   integer, intent(in) :: landmask(:,:)         ! landmask(nlon,nlat) land/ocean mask (0=ocean and 1=land)
   integer, intent(in) :: pftmask(:,:)          ! pftmask(nlon,nlat) pft real/fake mask (0=fake and 1=real)

   !Variables needed to create gridinfo_igbpz.csv
   !Use ndepAn to compute average N deposition for each grid cell (year = 1,ndtime)
   real(4), intent(in) :: ndepAn(:,:,:)         ! ndepAn(lon, lat, year) annual N deposition (gN/m2/yr)

   !Variables needed to create met.txt
   real(4), intent(in) :: tminday(:,:,:)        ! tminday(nlon,nlat,ndays) minimum daily air temperature (C)
   real(4), intent(in) :: tmaxday(:,:,:)        ! tmaxday(nlon,nlat,ndays) maximum daily air temperature (C)
   !real(4), intent(in) :: tlaiday(:,:,:)       ! tlai(nlon,nlat,ndays) total projected LAI
   real(4), intent(in) :: gpp(:,:,:)            ! Gross Primary Production (gC/m2/day)
   real(4), intent(in) :: h2osoiday(:,:,:,:)    ! h2osoiday(nlon,nlat,nlevsoi,ndays) daily average liquid volumetric soil water BY LAYER (mm3/mm3)
   real(4), intent(in) :: h2ofrznsoiday(:,:,:,:)! h2ofrznsoiday(nlon,nlat,nlevsoi,ndays) daily average frozen volumetric soil water BY LAYER (mm3/mm3)
                                                !   (vegetated landunits only) 
   real(4), intent(in) :: tsoiday(:,:,:,:)      ! tsoiday(nlon,nlat,nlevsoi,ndays) daily soil temperature r BY LAYER (K)
                                                !   (vegetated landunits only) 
    
   !Note: pctpft, pctsand, pctclay are from surface data set and are DOUBLE. Use real(8).
   real(8), intent(in) :: pctpft(:,:,:)         ! pctpft(nlon,nlat,npft) percent of each PFT (1-17) in grid cell
   real(8), intent(in) :: pctsand(:,:,:)        ! pctsand(nlon,nlat,nlevsoi) percent sand
   real(8), intent(in) :: pctclay(:,:,:)        ! pctclay(nlon,nlat,nlevsoi) percent clay
!WW added watsat & watfield to surface data set 
!orig. included in h0 files
    real(8), intent(in) :: watsat(:,:,:)         ! watsat(nlon,nlat,nlevsoi)water saturated (m3/m3)
    real(8), intent(in) :: watfc(:,:,:)          ! watfc(nlon,nlat,nlevsoi)water field capacity (m3/m3)

  !LOCAL VARIABLES
   character(len=5 )    :: sCellName            ! name appended to files for individual cell results
   integer, allocatable :: cellMissing(:,:)     ! cellMissing(nlon,nlat), 0=simulate, 1=missing data, don't simulate
   integer, allocatable :: cellid(:,:)          ! cell IDs read from fileCLMgridIGBP, written to fileGridInfo and fileGridSoil
   integer, allocatable :: vegIGBP(:,:)         ! vegIGBP(nlon,nlat) IGBP veg classification on CLM grid
   !real(4), allocatable :: xcnpp(:,:,:)        ! xcnpp(nlon,nlat,ndays) Net Primary Production (gC/m2/yr)
   real(4), allocatable :: xtairk(:,:,:)        ! xtairk(nlon,nlat,ndays) daily average air temperature (K) 
   real(4), allocatable :: ndepDay(:,:,:)       ! ndepDay(nlon,nlat,ndays) daily average N deposition (gN/m2/day)
   real(4), allocatable :: xtsoil(:,:,:,:)      ! xtsoil(nlon,nlat,XSOIL,ndays), CASACNP soil temperature by layer (K)
   real(4), allocatable :: xmoist(:,:,:,:)      ! xmoist(nlon,nlat,XSOIL,ndays), CASACNP volumetric soil liquid water content by layer 
   real(4), allocatable :: xfrznmoist(:,:,:,:)  ! xfrznmoist(nlon,nlat,XSOIL,ndays), CASACNP volumetric soil frozen water content by layer 
   real(4), allocatable :: levgrndthickness(:)  ! levgrndthickness(nlevsoi) thickness of soil layers (m)
   real levlower, levupper                      ! used to compute CLM layer thickness from levgrnd
   integer :: i,j,k,ndyr                        ! loop indices
   integer :: day                               ! day of year (1 - 365)
   !integer :: nyears                           ! current calendar year, number of years
   integer :: ndyr1, ndyr2                      ! time indices in ndepAn(lon, lat, ndyr) corresponding to calendar years iyear1 and iyear2
   real(4) :: ndepAvg                           ! average N deposition from ndyr1 to ndyr1 (gN/m2/yr)

   integer :: cellCnt                           ! grid cell count
   integer :: ist, iso, pft, vtypex             ! variables in fileGridInfo
   integer :: doy1, doy2, doy3, doy4, phase     ! place holder variables in fileGridInfo
   real(4) :: pwea, pdust, nfix                 ! variables in fileGridInfo
   real(4) :: sand, clay, silt                  ! soil texture variables in fileGridSoil (%sand,%clay,%silt)
   real(4) :: wwilt, wfield, wsat               ! volumetric percentages in fileGridSoil (wilting pt, field capacity, saturation)
   real(4) :: landarea

   character(len=100) :: ncfileMet              ! met.nc. Output file.
!  character(len=100) :: fileMet                ! met.txt. Output file (one cell only).
   character(len=100) :: fileGridInfo           ! gridinfo_igpbz.csv. Output file.
   character(len=100) :: fileGridSoil           ! gridinfo_soil.csv (Not an original casacnp file). Output file.
!  character(len=100) :: fileGridPFT            ! gridinfo_pft.csv (Not an original casacnp file). Output file.
   character(len=100) :: fileCLMgridIGBP        ! clmGrid_IGBP.csv (IGBP veg types on CLM grid). Input file.
   character*1::c                               ! text file delimeter
   character*4 :: xyear, xyear2
   c = ','

   !------------------------------------------------------------------------------------------
   ! ATTENTION  
   ! Make sure nlevsoi=10 < nlevgrnd=15 
   ! nlevsoi comes from the surface data, nlevgrnd comes from the history files.
   !------------------------------------------------------------------------------------------
 
   !nyears = iyear2 - iyear1 + 1

   allocate(cellMissing(nlon,nlat))
   allocate(cellid(nlon,nlat))
   allocate(vegIGBP(nlon,nlat))
   !allocate(xcnpp(nlon,nlat,ndays))
   allocate(xtairk(nlon,nlat,ndays))
   allocate(ndepDay(nlon,nlat,ndays))
   allocate(xtsoil(nlon,nlat,XSOIL,ndays))
   allocate(xmoist(nlon,nlat,XSOIL,ndays))
   allocate(xfrznmoist(nlon,nlat,XSOIL,ndays))
   allocate(levgrndthickness(1:nlevsoi+1))

   if (iyear1 .eq. calyear) then
      ! First year of the annual sequence
      allocate(wrtcell(nlon,nlat))
   endif

   cellMissing(:,:) = 0
   cellid(:,:) = 0
   xtairk(:,:,:) = MISSING_VALUE 
   ndepDay(:,:,:) = MISSING_VALUE 
   xtsoil(:,:,:,:) = MISSING_VALUE 
   xmoist(:,:,:,:) = MISSING_VALUE 
   xfrznmoist(:,:,:,:) = MISSING_VALUE 
   !xcnpp(:,:,:) = MISSING_VALUE

   if (iyear1 .le. ndep_calyear(1)) then
      ndyr1 = 1
   else if (iyear1 .ge. ndep_calyear(ndtime)) then
      ndyr1 = ndtime
   else
      ndyr1 = iyear1 - ndep_calyear(1) + 1
   endif
   if (iyear2 .le. ndep_calyear(1)) then
      ndyr2 = 1
   else if (iyear2 .ge. ndep_calyear(ndtime)) then
      ndyr2 = ndtime
   else
      ndyr2 = iyear2 - ndep_calyear(1) + 1
   endif

   ! ndyr is the position in ndepAn array that contains N dep for the current year
   if (calyear .le. ndep_calyear(1)) then
      ndyr = 1
   else if (calyear .ge. ndep_calyear(ndtime)) then
      ndyr = ndtime
   else
      ndyr = calyear - ndep_calyear(1) + 1
   endif

   if (verbose .ge. 0) then
      write(*,*) '  writeCASACNPfiles: calyear =', calyear, ' ndyr =', ndyr
      write(*,*) '  writeCASACNPfiles: ndep_calyear(ndyr1) =', ndep_calyear(ndyr1)
      write(*,*) '  writeCASACNPfiles: ndep_calyear(ndyr2) =', ndep_calyear(ndyr2)
      write(*,*) '  writeCASACNPfiles: ndep_calyear(ndyr) =', ndep_calyear(ndyr)
   endif
  
   cellCnt = 0

   do i = 1, nlon
      do j = 1, nlat
         cellCnt = cellCnt + 1
         cellid(i,j) = cellCnt

         do day = 1, ndays
            if (gpp(i,j,day) .lt. MISSING_VALUE) then
               !xcnpp(i,j,day) = 0.5 * gpp(i,j,day)
            else
               cellMissing(i,j) = 1
            endif
            if ((tminday(i,j,day) .lt. MISSING_VALUE) .and. (tmaxday(i,j,day) .lt. MISSING_VALUE)) then
               !Air temperature (convert C to K)
               xtairk(i,j,day) = 0.5 * (tminday(i,j,day) + tmaxday(i,j,day)) + KELVIN
            else
               cellMissing(i,j) = 1
            endif
         enddo

         if (ndepAn(i,j,ndyr).lt. MISSING_VALUE) then
            do day = 1, ndays
               ndepDay(i,j,day) = ndepAn(i,j,ndyr) / 365
            enddo
         else
            cellMissing(i,j) = 1             
         endif

      enddo
   enddo
    
   levupper = 0
   do k = 1, nlevsoi
      levlower = levgrnd(k) + (levgrnd(k+1) - levgrnd(k))/2
      levgrndthickness(k) = levlower - levupper
      levupper = levlower
   end do

   !Compute the soil depth of the 10th soil layer.
   !Assuming nlevsoi=10 and levgrnd has 15=levgrnd elements.
   levlower = levgrnd(nlevsoi) + (levgrnd(nlevsoi+1) - levgrnd(nlevsoi))/2

   !Find cells that have missing soils data
   do day = 1, ndays
      do k = 1, nlevsoi
         do j = 1, nlat
            do i = 1, nlon
               if (tsoiday(i,j,k,day) .ge. MISSING_VALUE) then
                  cellMissing(i,j) = 1
               endif
               if (h2osoiday(i,j,k,day) .ge. MISSING_VALUE .or. &
                   h2ofrznsoiday(i,j,k,day) .ge. MISSING_VALUE) then
                  cellMissing(i,j) = 1
               endif
            enddo
         enddo
      enddo
   enddo


   !ATTENTION: select correct year, or average ndep below. Area units?


   ! Write grid files only once, during the first year`
   if (calyear .eq. iyear1) then
    
!!     !Define dimensions and variables in the met.nc file. 
!!     call createCASACNPMetNcFile(ncFileMet, iyear1, iyear2, ndays, nlon, nlat, XSOIL, &
!!                                 lon1d, lat1d, landfrac, cellMissing, cellid)

       !Read the IGBP veg classes based on CLM grid
       fileCLMgridIGBP = "clmGrid_IGBP.csv"
       call readCLMgridIGBPGrid(fileCLMgridIGBP, nlon, nlat, lon1d, lat1d, vegIGBP, wrtcell, cellid, wrtSingleCells)

       if (wrtSingleCells == 1) then
          !Write comma-delimted gridinfo_igbpz.csv file for INDIVIDUAL cells
          do i = 1, nlon
             do j = 1, nlat

                if (wrtcell(i,j) == 1 .AND. cellMissing(i,j) == 0) then

                   call CreateMetFileNamePt(ncFileMet, cellid(i,j), iyear1, iyear2)

                   !Define dimensions and variables in the single-cell met.nc file. 
                   call createCASACNPMetNcFilePt(ncFileMet, iyear1, iyear2, ndays, i, j, nlat, nlon, XSOIL, &
                                   lon1d, lat1d, landfrac, cellMissing, cellid)

                   call GetCellName(sCellName, cellid(i,j))

                   !Write comma-delimted gridinfo_igbpz.csv file for the grid cell
                   fileGridInfo = "gridinfo_igbpz_pt" // trim(sCellName) // ".csv"
                   open(unit=201,file=trim(fileGridInfo))
                   print *, "Writing ", trim(fileGridInfo), " file for year ", calyear
                   write(201,fmt='(a110)') "ijcam,lat,lon,ivt_igbp,ist,iso,landarea,Ndep,Nfix,Pwea, &
                                            Pdust,doyP1,doyP2,doyP3,doyP4,Phase(1),ivcasa,ilat,ilon"
                   close(unit=201)
                   
                   !Write comma-delimted gridinfo_soil.csv file for the grid cell
                   fileGridSoil = "gridinfo_soil_pt" // trim(sCellName) // ".csv"
                   open(unit=202,file=trim(fileGridSoil))
                   print *, "Writing ", trim(fileGridSoil), " for year ", calyear
                   write(202,fmt='(a46)') "ijcam,lat,lon,sand,clay,silt,wwilt,wfield,wsat"
                   close(unit=202)
                endif
                if (wrtcell(i,j) == 1 .AND. cellMissing(i,j) == 1) then
                   write(*,111) ' No files written for cell ID = ', cellid(i,j), &
                                'latitude = ', lat1d(j), &
                                'longitude = ', lon1d(i), &
                                ': one or more CLM variables are missing.'
111 format(a32,i6,2x,a11,f7.2,2x,a12,f7.2,a40)
                endif
             end do
          end do
       else

          !Insert first and last calendar years in the met.nc file name.
          write(xyear,'(i4)') iyear1
          write(xyear2,'(i4)') iyear2
          ncFileMet = 'met_' // xyear // '_' // xyear2 // '.nc'

          !Define dimensions and variables in the met.nc file for ALL cells. 
          call createCASACNPMetNcFile(ncFileMet, iyear1, iyear2, ndays, nlon, nlat, XSOIL, &
                                   lon1d, lat1d, landfrac, cellMissing, cellid)

          !Write comma-delimted gridinfo_igbpz.csv file for ALL cells
          fileGridInfo = "gridinfo_igbpz.csv"
          open(unit=201,file=trim(fileGridInfo))
          print *, "Writing gridinfo_igbpz.csv file for year ", calyear
          write(201,fmt='(a110)') "ijcam,lat,lon,ivt_igbp,ist,iso,landarea,Ndep,Nfix,Pwea, &
                                   Pdust,doyP1,doyP2,doyP3,doyP4,Phase(1),ivcasa,ilat,ilon"
          
          fileGridSoil = "gridinfo_soil.csv"
          open(unit=202,file=trim(fileGridSoil))
          print *, "Writing gridinfo_soil.csv file for year ", calyear
          write(202,fmt='(a46)') "ijcam,lat,lon,sand,clay,silt,wwilt,wfield,wsat"
       endif 
    
!      fileGridPFT = "gridinfo_pft.csv"
!      open(unit=203,file=trim(fileGridPFT))
!      print *, "Writing gridinfo_pft.csv file for year ", calyear
!      write(203,fmt='(a83)') "ijcam,lat,lon,landfrac,landmask,pftmask,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,"
    
       ! CLM PFTs (as of Dec. 2013)
       ! 0      bare ground
       ! 1      Needleleaf evergreen tree – temperate NET Temperate
       ! 2      Needleleaf evergreen tree - boreal NET Boreal
       ! 3      Needleleaf deciduous tree – boreal NDT Boreal
       ! 4      Broadleaf evergreen tree – tropical BET Tropical
       ! 5      Broadleaf evergreen tree – temperate BET Temperate
       ! 6      Broadleaf deciduous tree – tropical BDT Tropical
       ! 7      Broadleaf deciduous tree – temperate BDT Temperate
       ! 8      Broadleaf deciduous tree – boreal BDT Boreal
       ! 9      Broadleaf evergreen shrub - temperate BES Temperate
       ! 10     Broadleaf deciduous shrub – temperate BDS Temperate
       ! 11     Broadleaf deciduous shrub – boreal BDS Boreal
       ! 12     C3 arctic grass -
       ! 13     C3 grass -
       ! 14     C4 grass -
       ! 15     Crop1 -
       ! 16     Crop2
    
       call srand(seed)
       vtypex = MISSING_INT_SM  ! Not used by casacnp when IGBP vegtypes are used
       doy1 = MISSING_INT_SM    ! Variable in file that is not used
       doy2 = MISSING_INT_SM    ! Variable in file that is not used
       doy3 = MISSING_INT_SM    ! Variable in file that is not used
       doy4 = MISSING_INT_SM    ! Variable in file that is not used
       phase = MISSING_INT_SM   ! Variable in file that is not used
       pwea = 0.0               ! P weathering
       pdust = 0.0              ! P dust input
       nfix = 0.1               ! N fixation

       print *
       print *, "WARNING: if casaclm is to be run for C, N, and P (icycle=3), you will "
       print *, "need to correctly define soil order (iso) in file ", TRIM(fileGridInfo)
       print *, "for P calculations. The current assignment of iso is not correct."
       print *

       do i = 1, nlon
          do j = 1, nlat

             !ist = INT(7*rand()) + 1   ! Soil type classification (1-7). Not needed if using gridinfo_soil.csv.
             ist = MISSING_INT_SM

             !iso = INT(12*rand()) + 1  ! Temporary assignment for soil order (1-12). Determine this if running P (icycle=3)
             iso = MISSING_INT_SM       ! An undefined value causes an error. So using a random number for now. 

             pft = vegIGBP(i,j)

             sand = MISSING_FLT_SM 
             silt = MISSING_FLT_SM 
             clay = MISSING_FLT_SM 
             wwilt = MISSING_FLT_SM 
             wfield = MISSING_FLT_SM 
             wsat = MISSING_FLT_SM 

             if (cellMissing(i,j) .eq. 0) then
                sand = 0.0
                clay = 0.0
! WW added
                wsat = 0.0
                wfield = 0.0
                !ATTENTION: average top 50 cm only (1..lyr50cm) instead of 1..10?
                do k = 1, nlevsoi
                   sand   = sand   + real(pctsand(i,j,k))*levgrndthickness(k)
                   clay   = clay   + real(pctclay(i,j,k))*levgrndthickness(k)
! WW added
                   wsat   = wsat   + real(watsat(i,j,k))*levgrndthickness(k)
                   wfield = wfield + real(watfc(i,j,k))*levgrndthickness(k)
                enddo
                sand = sand / levlower
                clay = clay / levlower
                silt = 100.0 - sand - clay
! WW added
                wsat   = wsat   / levlower
                wfield = wfield / levlower

! WW modified   call hydraulicProperties(sand,clay,silt,wwilt,wfield,wsat)
                call hydraulicProperties(sand,clay,silt,wwilt,wfield)

                if ((landfrac(i,j) .lt. MISSING_VALUE) .AND. (area(i,j) .lt. MISSING_VALUE)) then
                   !Convert km^2 to m^2
                   landarea = landfrac(i,j)*area(i,j)*1000000.0
                else
                   landarea = MISSING_FLT_SM
                endif
       
                ! Note: ndepAvg, the value written to the gridinfo file, is the average for all years,  
                !  while ndepDay, written to met.nc file, is daily deposition for the current year only. -mdh 8/22/2016
                ndepAvg = 0.0
                do ndyr = ndyr1, ndyr2
                   ndepAvg = ndepAvg + ndepAn(i,j,ndyr)
                enddo
                ndepAvg = ndepAvg / (ndyr2 - ndyr1 + 1)

               if (wrtSingleCells == 1) then
                 !Create these "gridinfo" files for each cell in a subset of cells 
                 if (wrtcell(i,j) == 1) then
                       call GetCellName(sCellName, cellid(i,j))
                       !Make sure file is named the same as it is above
                       fileGridInfo = "gridinfo_igbpz_pt" // trim(sCellName) // ".csv"
                       open(unit=201,file=trim(fileGridInfo),access='append')
                       write(unit=201, fmt='(i6,a1,2(f7.2,a1),3(i3,a1),f20.4,a1,4(f8.5,a1),9(i4,a1))') &
                          cellid(i,j),c, lat1d(j),c, lon1d(i),c, pft,c, ist,c, iso,c, landarea,c, &
                          ndepAvg,c, nfix,c, pwea,c, pdust,c, &
                          doy1,c, doy2,c, doy3,c, doy4,c, phase,c, vtypex,c, 1,c, 1,c
                          !!doy1,c, doy2,c, doy3,c, doy4,c, phase,c, vtypex,c, j,c, i,c
                          !!ilat and ilon in gridinfo_igbpz*.csv should = 1 for point files
                        close(unit=201)
       
                       !Make sure file is named the same as it is above
                       fileGridSoil = "gridinfo_soil_pt" // trim(sCellName) // ".csv"
                       open(unit=202,file=trim(fileGridSoil),access='append')
                       write(unit=202, fmt='(i6,a1,2(f7.2,a1),6(f8.5,a1))') &
                          cellid(i,j),c, lat1d(j),c, lon1d(i),c, sand/100,c, clay/100,c, silt/100,c, &
!! WWW modified           wwilt/100,c, wfield/100,c, wsat/100,c
                          wwilt/100,c, wfield,c, wsat,c
                       close(unit=202)
                  endif
               else
                   !All cells are written to a single file
                   write(unit=201, fmt='(i6,a1,2(f7.2,a1),3(i3,a1),f20.4,a1,4(f8.5,a1),9(i4,a1))') &
                      cellid(i,j),c, lat1d(j),c, lon1d(i),c, pft,c, ist,c, iso,c, landarea,c, &
                      ndepAvg,c, nfix,c, pwea,c, pdust,c, &
                      doy1,c, doy2,c, doy3,c, doy4,c, phase,c, vtypex,c, j,c, i,c
      
                   write(unit=202, fmt='(i6,a1,2(f7.2,a1),6(f8.5,a1))') &
                      cellid(i,j),c, lat1d(j),c, lon1d(i),c, sand/100,c, clay/100,c, silt/100,c, &
!! WWW modified       wwilt/100,c, wfield/100,c, wsat/100,c
                      wwilt/100,c, wfield,c, wsat,c

               endif


   
!               write(unit=203, fmt='(i6,a1,3(f7.2,a1),2(i6,a1))',advance='no') &
!                  cellid(i,j),c, lat1d(j),c, lon1d(i),c, landfrac(i,j),c, landmask(i,j),c, pftmask(i,j),c 
!                  do p = 1, npft
!                     write(unit=203, fmt='(f7.2,a1)',advance='no') &
!                        pctpft(i,j,p),c
!                  enddo
!               write(unit=203,fmt='(a1)') " "

             endif
          enddo
       enddo
       close(unit=201)
       close(unit=202)
!      close(unit=203)
       print *, "Done writing gridinfo_igbpz.csv file(s)..."
       print *, "Done writing gridinfo_soil.csv file(s)..."
!      print *, "Done writing gridinfo_pft.csv file(s)..."
   endif  !endif (calyear .eq. iyear1) then
    
   !Calculate CASACNP soil temperature and moisture to variables from CLM variables
   !ATTENTION: Document CLM to CASA soil layer mapping here.
   do j = 1, nlat
      do i = 1, nlon
         if (cellMissing(i,j) .eq. 0) then
            do day = 1, ndays
               !soil temperature
               xtsoil(i,j,1,day) = tsoiday(i,j,1,day)
               xtsoil(i,j,2,day) = (tsoiday(i,j,2,day)*levgrndthickness(2) + tsoiday(i,j,3,day)*levgrndthickness(3)) &
                                   / (levgrndthickness(2) + levgrndthickness(3))
               xtsoil(i,j,3,day) = (tsoiday(i,j,4,day)*levgrndthickness(4) + tsoiday(i,j,5,day)*levgrndthickness(5)) &
                                   / (levgrndthickness(4) + levgrndthickness(5))
               xtsoil(i,j,4,day) = (tsoiday(i,j,6,day)*levgrndthickness(6) + tsoiday(i,j,7,day)*levgrndthickness(7)) &
                                   / (levgrndthickness(6) + levgrndthickness(7))
               xtsoil(i,j,5,day) = (tsoiday(i,j,8,day)*levgrndthickness(8) + tsoiday(i,j,9,day)*levgrndthickness(9)) &
                                   / (levgrndthickness(8) + levgrndthickness(9))
               xtsoil(i,j,6,day) = tsoiday(i,j,10,day)
   
               ! volumetric soil liquid water content
               xmoist(i,j,1,day) = h2osoiday(i,j,1,day)
               xmoist(i,j,2,day) = (h2osoiday(i,j,2,day)*levgrndthickness(2) + h2osoiday(i,j,3,day)*levgrndthickness(3)) &
                                   / (levgrndthickness(2) + levgrndthickness(3)) 
               xmoist(i,j,3,day) = (h2osoiday(i,j,4,day)*levgrndthickness(4) + h2osoiday(i,j,5,day)*levgrndthickness(5)) &
                                   / (levgrndthickness(4) + levgrndthickness(5)) 
               xmoist(i,j,4,day) = (h2osoiday(i,j,6,day)*levgrndthickness(6) + h2osoiday(i,j,7,day)*levgrndthickness(7)) &
                                   / (levgrndthickness(6) + levgrndthickness(7)) 
               xmoist(i,j,5,day) = (h2osoiday(i,j,8,day)*levgrndthickness(8) + h2osoiday(i,j,9,day)*levgrndthickness(9)) &
                                   / (levgrndthickness(8) + levgrndthickness(9)) 
               xmoist(i,j,6,day) = h2osoiday(i,j,10,day) 

               ! volumetric soil frozen water content
               xfrznmoist(i,j,1,day) = h2ofrznsoiday(i,j,1,day)
               xfrznmoist(i,j,2,day) = (h2ofrznsoiday(i,j,2,day)*levgrndthickness(2) &
                                        + h2ofrznsoiday(i,j,3,day)*levgrndthickness(3)) &
                                        / (levgrndthickness(2) + levgrndthickness(3)) 
               xfrznmoist(i,j,3,day) = (h2ofrznsoiday(i,j,4,day)*levgrndthickness(4) &
                                        + h2ofrznsoiday(i,j,5,day)*levgrndthickness(5)) &
                                        / (levgrndthickness(4) + levgrndthickness(5)) 
               xfrznmoist(i,j,4,day) = (h2ofrznsoiday(i,j,6,day)*levgrndthickness(6) &
                                        + h2ofrznsoiday(i,j,7,day)*levgrndthickness(7)) &
                                        / (levgrndthickness(6) + levgrndthickness(7)) 
               xfrznmoist(i,j,5,day) = (h2ofrznsoiday(i,j,8,day)*levgrndthickness(8) &
                                        + h2ofrznsoiday(i,j,9,day)*levgrndthickness(9)) &
                                        / (levgrndthickness(8) + levgrndthickness(9)) 
               xfrznmoist(i,j,6,day) = h2ofrznsoiday(i,j,10,day)

               !write(*,*)
               !write(*,*) 'writeCASACNPfiles: xmoist(', i,j,1,day, ') =', xmoist(i,j,1,day)
               !write(*,*) 'writeCASACNPfiles: xmoist(', i,j,2,day, ') =', xmoist(i,j,2,day)
               !write(*,*) 'writeCASACNPfiles: xmoist(', i,j,3,day, ') =', xmoist(i,j,3,day)
               !write(*,*) 'writeCASACNPfiles: xmoist(', i,j,4,day, ') =', xmoist(i,j,4,day)
               !write(*,*) 'writeCASACNPfiles: xmoist(', i,j,5,day, ') =', xmoist(i,j,5,day)
               !write(*,*) 'writeCASACNPfiles: xmoist(', i,j,6,day, ') =', xmoist(i,j,6,day)
            enddo
         endif
      enddo
   enddo


!  !For testing...
!  !Write tab-delimted met file for one cell and one year
!  print *, "Writing met.txt file for year ", calyear
!  iyear=1
!  write(xyear,'(i4)') calyear
!  fileMet = 'met_' // xyear // '.txt'
!  i = LONINDX
!  j = LATINDX
!  open(unit=204,file=trim(fileMet))
!  write(unit=204, fmt='(i1)') iyear
!  do day = 1, ndays
!     write(unit=204, fmt='(T9,i1,T17,i3,T25,f4.2,T33,f5.3,T41,f8.4,T57,6(F8.4,2x),6(F6.4,2x))') &
!        iyear,day,tlaiday(i,j,day),xcnpp(i,j,day),xtairk(i,j,day),(xtsoil(i,j,k,day),k=1,6),(xmoist(i,j,k,day),k=1,6)
!  enddo
!  close(unit=204)
!  print *, "Done writing met.txt file..."

   if (wrtSingleCells /= 1) then

      !Recreate the name of the met.nc file
      write(xyear,'(i4)') iyear1
      write(xyear2,'(i4)') iyear2
      ncFileMet = 'met_' // xyear // '_' // xyear2 // '.nc'

      !call writeCASACNPMetNcFile(ncFileMet, calyear, iyear1, iyear2, nlon, nlat, ndays, &
      !                           XSOIL, xcnpp, gpp, tlaiday, xtairk, xtsoil, xmoist, xfrznmoist, ndepDay)
      call writeCASACNPMetNcFile(ncFileMet, calyear, iyear1, iyear2, nlon, nlat, ndays, &
                                 XSOIL, gpp, xtairk, xtsoil, xmoist, xfrznmoist, ndepDay)
   else
      !write(*,*) 'wrtSingleCells = ', wrtSingleCells, 'calyear = ', calyear
      do i = 1, nlon
         do j = 1, nlat
            if (wrtcell(i,j) == 1 .AND. cellMissing(i,j) == 0) then
               call CreateMetFileNamePt(ncFileMet, cellid(i,j), iyear1, iyear2)
               !call writeCASACNPMetNcFilePt(ncFileMet, calyear, iyear1, iyear2, i, j, nlat, nlon, ndays, &
               !                             XSOIL, xcnpp, gpp, tlaiday, xtairk, &
               !                             xtsoil, xmoist, xfrznmoist, ndepDay)
               call writeCASACNPMetNcFilePt(ncFileMet, calyear, iyear1, iyear2, i, j, nlat, nlon, ndays, &
                                            XSOIL, gpp, xtairk, xtsoil, xmoist, xfrznmoist, ndepDay)
            endif
         end do 
      end do
   endif

   !deallocate(cellMissing, cellid, vegIGBP, xcnpp, xtairk, ndepDay, xtsoil, &
   !           xmoist, xfrznmoist, levgrndthickness)
   deallocate(cellMissing, cellid, vegIGBP, xtairk, ndepDay, xtsoil, &
              xmoist, xfrznmoist, levgrndthickness)

   end subroutine writeCASACNPfiles

!----------------------------------------------------------------------------------------------------

   SUBROUTINE CreateMetFileNamePt(ncFileMet, cellid, iyear1, iyear2)
   
       character(len=100), intent(INOUT) :: ncFileMet
       integer, INTENT(IN)               :: cellid
       integer, INTENT(IN)               :: iyear1, iyear2
   
       !Local Variables
       character(len=4) syear1, syear2
       character(len=5) sCellName 
     
       call GetCellName(sCellName, cellid)
     
       !Insert first and last calendar years in the met.nc file name.
       write(syear1,'(i4)') iyear1
       write(syear2,'(i4)') iyear2
       ncFileMet = 'met_pt' // trim(sCellName) // '_' // syear1 // '_' // syear2 // '.nc'
       !write(*,*) 'Created file name ', trim(ncFileMet)
   
   END SUBROUTINE CreateMetFileNamePt

!----------------------------------------------------------------------------------------------------

   SUBROUTINE GetCellName(sCellName, cellid)
   
       character(len=5), INTENT(INOUT)   :: sCellName 
       integer, INTENT(IN)               :: cellid
   
       write(sCellName,'(i5)') cellid    ! copy the value of cellid to string sCellName
       if (cellid .lt. 10000) then       ! make sure cellid is 5 digits, add leading zeros if necessary
       sCellName(1:1) = '0'
       endif
       if (cellid .lt. 1000) then
       sCellName(2:2) = '0'
       endif
       if (cellid .lt. 100) then
       sCellName(3:3) = '0'
       endif
       if (cellid .lt. 10) then
       sCellName(4:4) = '0'
       endif
     
   END SUBROUTINE GetCellName

!----------------------------------------------------------------------------------------------------

end module casacnp_model
