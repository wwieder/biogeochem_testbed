module canopy_model

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
!  !MODULE: canopy_model
!
!  !DESCRIPTION:
!   This module contains subroutines to run the Aggregated Canopy Model
!   for a CLM grid and to save daily GPP output in NetCDF format.
!
! !USES:
   use clm_common

! !I/O text files
!  unit=16	Output (optional) : daily gpp text file for LATINDX, LONINDX

   implicit none
!  include 'netcdf.inc'

! !PUBLIC MEMBER FUNCTIONS:
   public :: Photosyn
   public :: Daylength
   public :: RunCanopyModel
   public :: NcWriteDailyGPPOutput

!  SEE clm_common.F90 for the following parameter values!
!  real(4), parameter :: SEC_PER_HOUR = 3600
!  real(4), parameter :: HOURS_PER_DAY = 24
!  real(4), parameter :: DAYS_PER_YEAR = 365
!  real(4), parameter :: MISSING_VALUE = 1.e+36
!  real(4), parameter :: PI =  3.14159265
!  integer, parameter :: LATINDX = 70		! Testing location
!  integer, parameter :: LONINDX = 20		! Testing location
!  integer, parameter :: baseyear = 1901	! Times in the climate files are days since 1901-01-01 00:00:00
!
!  Aggregated Canopy Model reference (Note: there have been code updates since this publication):  
! 
!  Williams, M., E. B. Rastetter, D. Fernandes, M. L. Goulden, G. Shaver,
!  and L. Johnson. 1997. Predicting gross primary productivity in
!  terrestrial ecosystems. Ecological Applications 7:882–894.

!  type paramType
!     sequence
!     real(4) :: a1
!     real(4) :: a2
!     real(4) :: a3
!     real(4) :: a4
!     real(4) :: a5
!     real(4) :: a6
!     real(4) :: a7
!     real(4) :: a8
!     real(4) :: a9
!     real(4) :: a10
!  end type paramType
!
!  type climType
!     sequence
!     real(4) :: tmax		! Daily maximum temperature (C)
!     real(4) :: tmin		! Daily minimum temperature (C)
!     real(4) :: co2		! atmospheric CO2 (umol/mol)
!     real(4) :: rad		! Total daily irradiance (MJ m-2 d-1)
!  end type climType
!
!  type stateType
!     sequence
!     real(4) :: DayLength	! Daylength (hours)
!     real(4) :: psid		! Difference between min. leaf water potential and soil water potential (MPa)
!     real(4) :: rtot		! Total plant-soil hydraulic resistance (MPa m-2 s mmol-1) (set to 0.1)  
!  end type stateType
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Aggregrated Canopy Model
!
! !INTERFACE:
   real(4) function Photosyn(day, latitude, tmin, tmax, rad, co2, nit, lai, rtot, lwpmin, swp, daylen)

! !DESCRIPTION:
!  Return the gross primary production, GPP (gC/m2/day), for the day.
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: day   	! Day of year [1, 365]
   real(4), intent(in) :: latitude     	! latitude (decimal degrees) [-90, 90]
   real(4), intent(in) :: tmin		! Daily minimum temperature (C)
   real(4), intent(in) :: tmax		! Daily maximum temperature (C)
   real(4), intent(in) :: rad		! Total daily irradiance (MJ m-2 d-1)
   real(4), intent(in) :: co2		! atmospheric CO2 (umol/mol)
   real(4), intent(in) :: nit		! foliar N concentration (gN/m2 leaf area)
   real(4), intent(in) :: lai		! Leaf Area Index (m2/m2)
   real(4), intent(in) :: rtot		! Total plant-soil hydraulic resistance (MPa m-2 s mmol-1) (set to 0.1)  
   real(4), intent(in) :: lwpmin	! Minimum leaf water potential (MPa) (for example, -2.0)
   real(4), intent(in) :: swp          	! Soil water potential (MPa)
   real(4), intent(out) :: daylen       ! Daylength (hours)


! !REVISION HISTORY:
!
! !LOCAL VARIABLES:

   real(4) :: ci		! concentration at the site of carboxylation (µmol/mole)
   real(4) :: cps		! light limitation ??
   real(4) :: e0		! canopy-level quantum yield (gC/MJ/m2/day)
   real(4) :: gs		! daily canopy conductance (units?)
   real(4) :: pp		! nit / gs * a1 * exp(a8 * tmax)
   real(4) :: qq 		! a3 - a4
   real(4) :: tav		! daily average temperature (C)
   real(4) :: trange		! daily temperature range (C)
   real(4) :: psid		! Difference between min. leaf water potential and soil water potential (MPa)
   real(4) :: gpp          	! gross primary production (gC/m2/day)

!  type (climType), pointer :: clim
!  type (stateType), pointer :: state
!  type (paramType), pointer :: param

   character*1 :: c

   !ACM parameters
   real(4) :: a1		! NUE parameter (gC GPP/gN in leaves)
   real(4) :: a2		! Day Length constant 
   real(4) :: a3		! Canopy CO2 compensation point
   real(4) :: a4		! Canopy CO2 half-saturation point
   real(4) :: a5		! Midsummer coefficient
   real(4) :: a6		! Temperature range constant of canopy conductance
   real(4) :: a7		! Maximum canopy quantum yield
   real(4) :: a8		! Temperature coefficient of NUE parameter 
   real(4) :: a9		! LAI-canopy quantum yield coefficient
   real(4) :: a10		! Temperature coefficient of canopy conductance(?)
  
!EOP
!-----------------------------------------------------------------------

   !ACM parameters (a1..a10)
!  a1 = 6.0			! NUE parameter (gC GPP/gN in leaves)
!  a2 = 0.01569			! Day Length constant
!  a3 =  4.223 			! Canopy CO2 compensation point (umol/mole)
!  a4 = 208.86			! Canopy CO2 half-saturation point (umol/mole)
!  a5 =  0.04531		! Midsummer coefficient
!  a6 = 0.378			! Temperature range constant of canopy conductance
!  a7 = 7.193			! Maximum canopy quantum yield
!  a8 = 0.0111			! Temperature coefficient of NUE parameter 
!  a9 = 2.10 			! LAI-canopy quantum yield coefficient
!  a10 = 0.7897			! Temperature coefficient of canopy conductance(?)

! Updated parameters 4/7/2014 to increase GPP
   a1 = 6.0			! NUE parameter (gC GPP/gN in leaves)
   a2 = 0.01800			! Day Length constant
   a3 =  4.223 			! Canopy CO2 compensation point (umol/mole)
   a4 = 208.86			! Canopy CO2 half-saturation point (umol/mole)
   a5 =  0.09000		! Midsummer coefficient
   a6 = 0.31500			! Temperature range constant of canopy conductance
   a7 = 7.193			! Maximum canopy quantum yield
   a8 = 0.01800			! Temperature coefficient of NUE parameter 
   a9 = 0.98900			! LAI-canopy quantum yield coefficient
   a10 = 0.8730			! Temperature coefficient of canopy conductance(?)

   gpp = 0

!  c = ','
!  write(*, fmt='(1(f0.2,a1), 1(i6,a1), 6(f0.6,a1))') latitude,c, day,c, &
!        tmin,c, tmax,c, rad,c, lwpmin,c, swp,c, gpp,c

   ! tav: average daily temperature (C)
   tav = (tmax + tmin) / 2.0
   daylen = Daylength(day, latitude)
   psid = lwpmin - swp

   if(tav .gt. 0.) then       
      trange = tmax - tmin	! daily temperature range (C)
      gs = abs(psid)**a10 / ((a6* rtot + 0.5*trange))
      pp = nit / gs * a1 * exp(a8 * tmax)
      qq = a3 - a4
      ci = 0.5 * (co2 + qq - pp + ((co2 + qq - pp)**2 - 4 * &
          (co2 * qq - pp * a3)) **0.5)
      e0 = a7 * lai**2/(lai**2 + a9)
      cps = e0* rad * gs * (co2 - ci)/(e0 * rad + &
            gs * (co2-ci))
      
      !GPP for the day (gC/m2/day)
      gpp = cps * (a2 * daylen + a5)
   else
      gpp = 0.
   endif

   Photosyn = gpp

   end function Photosyn

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Daylength function
!
! !INTERFACE:
   real(4) function Daylength(jday, latitude)

! !DESCRIPTION
!  Returns daylength (hours)
   
! !ARGUMENTS
   integer :: jday		! day of year (1 - 365)
   real(4) :: latitude		! latitude in decimal degrees

  !LOCAL VARIABLES
   real(4) :: rlat		! latitude in radians
   real(4) :: adelt
   real(4) :: ahou
   real(4) :: temp1, temp2

   ! Convert latitude from decimal degrees to radians
   if (abs(abs(latitude) - 90.0) .lt. 0.01) then
      !Adjust for latitude +90, -90 because tangent is undefined there.
      rlat = (latitude + 0.01) * PI / 180.0
   else
      rlat = latitude * PI / 180.0
   endif
    
   temp1 = 2.0 * PI * (jday - 77.0) / 365.0
   adelt = 0.4014 * sin(temp1)
   temp1 = 1.0 - (-tan(rlat) * (adelt))**2.0
   if (temp1 .lt. 0.0) temp1 = 0.0
   temp1 = sqrt(temp1)
   temp2 = -tan(rlat) * tan(adelt)
   ahou = atan2(temp1,temp2)
   Daylength = (ahou / PI) * 24.0

   end function Daylength

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!  RunCanopyModel()
!
! !INTERFACE:
   subroutine RunCanopyModel(year, nlon, nlat, ndays, lon1d, lat1d, &
                             area, topo, landfrac, landmask, pftmask, &
                             tminday, tmaxday, fsdsday, tlaiday, soipsiday, co2ppm, gpp)

! !DESCRIPTIONi
!  Run the aggregated canopy model for a year and write daily GPP to a netCDF file.
   
!  !ARGUMENTS:
   integer, intent(in) :: year		   ! calendar year (YYYY)
   integer, intent(in) :: nlon		   ! number of longitudes 
   integer, intent(in) :: nlat		   ! number of latitudes
   integer, intent(in) :: ndays		   ! total number of daily time periods 
   real(4), intent(in) :: lon1d(:)	   ! longitude(lon) (degrees_east)
   real(4), intent(in) :: lat1d(:)  	   ! latitude(lat) (degrees_north)
   real(4), intent(in) :: area(:,:)	   ! grid cell areas (km^2)
   real(4), intent(in) :: topo(:,:)	   ! grid cell topography (m)
   real(4), intent(in) :: landfrac(:,:)    ! land fraction
   integer, intent(in) :: landmask(:,:)    ! land/ocean mask (0=ocean and 1=land)
   integer, intent(in) :: pftmask(:,:)     ! pft real/fake mask (0=fake and 1=real)
   real(4), intent(in) :: tminday(:,:,:)   ! tminday(lon, lat, day) minimum daily air temperature (C)
   real(4), intent(in) :: tmaxday(:,:,:)   ! tmaxday(lon, lat, day) maximum daily air temperature (C)
   real(4), intent(in) :: fsdsday(:,:,:)   ! fsds(lon,lat,day) 24-hour average atmospheric incident solar radiation (W/m^2)
   real(4), intent(in) :: tlaiday(:,:,:)   ! tlai(lon,lat,day) total projected LAI
   real(4), intent(in) :: soipsiday(:,:,:) ! soipsiday(lon,lat,day) daily depth-weighted soil water potential(vegetated landunits only) (MPa)
   real(4), intent(in) :: co2ppm	   ! atmospheric CO2 (umol/mol)
   real(4), intent(out), allocatable :: gpp(:,:,:)      ! Gross Primary Production (gC/m2/day)

  !LOCAL VARIABLES
   integer :: i,j
   integer :: day			   ! day of year (1 - 365)
   character*1::c			   ! file delimeter
   
  !Aggregated canopy model
   real(4), allocatable :: time(:)         ! time (days since 1901-01-01 00:00:00)
   integer, allocatable :: mcdate(:)       ! current date (YYYYMMDD)
   real(4) :: rad			   ! Total daily irradiance (MJ m-2 d-1)
   real(4) :: nit			   ! foliar N concentration (gN/m2 leaf area)
   real(4) :: lai			   ! Leaf Area Index (m2/m2)
   real(4) :: rtot			   ! Total plant-soil hydraulic resistance (MPa m-2 s mmol-1) (set to 0.1)  
   real(4) :: lwpmin			   ! Minimum leaf water potential (for example, -2.5)
   real(4) :: swp          		   ! Soil water potential (MPa) (for example, 0.0)
   real(4) :: daylen          		   ! day length (hours)
   real(4) :: tmin			   ! Temporary variable for minimum temperature of the grid cell
   real(4) :: tmax			   ! Temporary variable for maximum temperature of the grid cell
   real(4) :: latitude			   ! Temporary variable for latitude of the grid cell
   real(4) :: longitude			   ! Temporary variable for longitude of the grid cell
   integer :: nYears			   ! used to compute "time" array
   character*4   :: xyear		   ! Calendar year converted to a string
   integer :: h
   character*100 :: ncfilename		   ! Name of output netcdf file

   allocate(gpp(1:nlon, 1:nlat, 1:ndays))
   allocate(time(1:ndays))
   allocate(mcdate(1:ndays))
   c = ','

!  !Put this code in calling subroutine to create a .csv output file for the canopy model 
!  if (doWriteDailyGPPTestFile .eq. 1) then
!     fileGPPOutput = 'gpp.csv'
!     open(unit=16, file=fileGPPOutput)
!     write(unit=16, fmt='(12(a7))') 'lat,', 'lon,', 'year,', 'day,', 'daylen,', &
!                                    'tmin,', 'tmax,', 'rad,', 'lwpmin,', 'swp,', 'tlai,', 'gpp,'
!  endif

   ! Run the aggregated canopy model
   do j = 1, nlon
      do i = 1, nlat
         do day = 1, ndays

            !nit = 1.92		! Use a N leaf conc. value from a deciduous forest (Gordon Bonan 10/28/2013)
            nit = 3.0		! Use a N leaf conc. increased to boost GPP (4/7/2014)
 				! 1.92 is the published value for Harvard Forest (Williams et al. 1997)
            rtot = 0.1		! Set to 0.1 (Quinn Thomas [rqthomas@vt.edu] 09/09/2013)
            lwpmin = -2.0	! Use -2.0 for minimum leaf water potential (Gordon Bonan 10/14/2013)
            swp = soipsiday(j,i,day)
            tmin = tminday(j,i,day)
            tmax = tmaxday(j,i,day)
            lai = tlaiday(j,i,day)
            longitude = lon1d(j)
            latitude = lat1d(i)
            gpp(j,i,day) = 0
      
            if (tmin .lt. MISSING_VALUE .AND. tmax .lt. MISSING_VALUE .AND. &
                fsdsday(j,i,day) .lt. MISSING_VALUE .AND. swp .lt. MISSING_VALUE) then
      
               rad = fsdsday(j,i,day) * HOURS_PER_DAY * SEC_PER_HOUR * 0.000001		! Convert W/m2 to MJ/m2/day
      
               gpp(j,i,day) = Photosyn(day, latitude, tmin, tmax, rad, co2ppm, nit, lai, rtot, lwpmin, swp, daylen)
                
               !Write output to a .csv file for testing purposes.
               if (doWriteDailyGPPTestFile .eq. 1) then
                  if (j .eq. LONINDX .AND. i .eq. LATINDX) then
                     write(unit=16, fmt='(2(f0.2,a1), 2(i6,a1), 8(f0.6,a1))') latitude,c, longitude,c, year,c, day,c, &
                        daylen,c, tmin,c, tmax,c, rad,c, lwpmin,c, swp,c, lai,c, gpp(j,i,day),c
                  endif
               endif
      
            else
               gpp(j,i,day) = MISSING_VALUE
            endif

            nYears = year - baseyear
            time(day) = nYears * DAYS_PER_YEAR + day 
            mcdate(day) = getmcdate(time(day))

         end do
      end do
   end do
 
   write(xyear,'(i4)') year	! copy year to string xyear
   ncfilename = 'GPP_Daily_Output'
   h = len(trim(ncfilename))
   ncfilename(h+1:h+1) = '.'
   ncfilename(h+2:h+5) = xyear
   ncfilename(h+6:h+20) = '-01-01-00000.nc'

   call NcWriteDailyGPPOutput(ncfilename, nlon, nlat, ndays, lon1d, lat1d, &
      area, topo, landfrac, landmask, pftmask, time, mcdate, gpp)

   end subroutine RunCanopyModel

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NcWriteDailyGPPOutput
!
! !INTERFACE:
   subroutine NcWriteDailyGPPOutput(ncfilename, nlon, nlat, ntimes, lon1d, lat1d, &
      area, topo, landfrac, landmask, pftmask, time, mcdate, gpp)

! !DESCRIPTION:
!  Create a NetCDF output file of daily GPP with attributes similar to 'h0' CSM files

   !ARGUMENTS
   character(len=*), intent(in) :: ncfilename        	! name of output ncfile
   integer, intent(in) :: nlon				! number of longitudes 
   integer, intent(in) :: nlat				! number of latitudes
   integer, intent(in) :: ntimes			! total number of monthly/hourly time periods returned
   real(4), intent(in) :: lon1d(:)	 		! lat(lat) (degrees_east)
   real(4), intent(in) :: lat1d(:)  	 		! lat(lat) (degrees_north)
   real(4), intent(in) :: area(:,:)	   		! grid cell areas (km^2)
   real(4), intent(in) :: topo(:,:)	   		! grid cell topography (m)
   real(4), intent(in) :: landfrac(:,:)    		! land fraction
   integer, intent(in) :: landmask(:,:)    		! land/ocean mask (0=ocean and 1=land)
   integer, intent(in) :: pftmask(:,:)     		! pft real/fake mask (0=fake and 1=real)
   real(4), intent(in) :: time(:)          		! days since 1901-01-01 00:00:00
   integer, intent(in) :: mcdate(:)          		! current date (YYYYMMDD)
   real(4), intent(in) :: gpp(:,:,:)          		! gpp(lon, lat, time) gross primary production (gC/m2/day)

   !LOCAL VARIABLES
   integer :: status
   integer :: ncid
   integer :: dimid_lon, dimid_lat, dimid_time	! NetCDF dimension IDs
   integer :: dims(3)
   ! NetCDF variable IDs
   integer :: varid_lon, varid_lat, varid_time, varid_mcdate, varid_gpp
   integer :: varid_area, varid_topo, varid_landfrac, varid_landmask, varid_pftmask
   character*100 :: attr_str			! String for assigning global and variable attributes
   character*10 :: date_string			! String for assigning date to global attributes
   character*8 :: time_string			! String for assigning time to global attributes

   if (verbose .ge. 0) print *, 'Writing GPP to file ', trim(ncfilename)

   dims(1) = 0
   dims(2) = 0
   dims(3) = 0

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

!  status = nf_def_dim(ncid, 'time', NF_UNLIMITED, dimid_time)
!  if (status /= nf_noerr) call handle_err(status, "time")

   ! Define variables

   dims(1) = dimid_lon
   status = nf_def_var(ncid, 'lon', NF_REAL, 1, dims, varid_lon)
   if (status /= nf_noerr) call handle_err(status, "lon")

   dims(1) = dimid_lat
   status = nf_def_var(ncid, 'lat', NF_REAL, 1, dims, varid_lat)
   if (status /= nf_noerr) call handle_err(status, "lat")

   dims(1) = dimid_time
   status = nf_def_var(ncid, 'time', NF_REAL, 1, dims, varid_time)
   if (status /= nf_noerr) call handle_err(status, "time")

   dims(1) = dimid_time
   status = nf_def_var(ncid, 'mcdate', NF_INT, 1, dims, varid_mcdate)
   if (status /= nf_noerr) call handle_err(status, "mcdate")

   ! Because dimensions in FORTRAN are in Column Major Order (the first array index varies the most rapidly)
   ! dimensions are in the opposite order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat

   status = nf_def_var(ncid, 'area', NF_REAL, 2, dims, varid_area)
   if (status /= nf_noerr) call handle_err(status, "area")

   status = nf_def_var(ncid, 'topo', NF_REAL, 2, dims, varid_topo)
   if (status /= nf_noerr) call handle_err(status, "topo")

   status = nf_def_var(ncid, 'landfrac', NF_REAL, 2, dims, varid_landfrac)
   if (status /= nf_noerr) call handle_err(status, "landfrac")

   status = nf_def_var(ncid, 'landmask', NF_INT, 2, dims, varid_landmask)
   if (status /= nf_noerr) call handle_err(status, "landmask")

   status = nf_def_var(ncid, 'pftmask', NF_INT, 2, dims, varid_pftmask)
   if (status /= nf_noerr) call handle_err(status, "pftmask")


! ATTENTION: these variables are not included.  Should they be?
!int mcsec(time) ;
!	mcsec:long_name = "current seconds of current date" ;
!	mcsec:units = "s" ;
!int mdcur(time) ;
!	mdcur:long_name = "current day (from base day)" ;
!int mscur(time) ;
!	mscur:long_name = "current seconds of current day" ;

   ! Because dimensions in FORTRAN are in Column Major Order (the first array index varies the most rapidly)
   ! dimensions are in the opposite order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat
   dims(3) = dimid_time
   status = nf_def_var(ncid, 'gpp', NF_REAL, 3, dims, varid_gpp)
   if (status /= nf_noerr) call handle_err(status, "gpp")

   ! Global attributes
   attr_str = 'CLM History file information'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'title', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "title")

   attr_str = 'NOTE: None of the variables are weighted by land fraction!'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'comment', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "comment")

   call get_time_and_date(date_string, time_string)
   attr_str = 'created on ' // date_string // ' ' // time_string
   status = nf_put_att_text(ncid, NF_GLOBAL, 'history', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "history")

   attr_str = 'Aggregated Canopy Model'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'source', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "source")


   ! Attributes of the variables

   ! Attributes of lon
   attr_str = 'coordinate longitude'
   status = nf_put_att_text(ncid, varid_lon, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   attr_str = 'degrees_east'
   status = nf_put_att_text(ncid, varid_lon, 'units', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "units")

   status = nf_put_att_real(ncid, varid_lon, '_FillValue', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_real(ncid, varid_lon, 'missing_value', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "missing_value")


   ! Attributes of lat
   attr_str = 'coordinate latitude'
   status = nf_put_att_text(ncid, varid_lat, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   attr_str = 'degrees_north'
   status = nf_put_att_text(ncid, varid_lat, 'units', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "units")

   status = nf_put_att_real(ncid, varid_lat, '_FillValue', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_real(ncid, varid_lat, 'missing_value', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "missing_value")


   ! Attributes of time
   attr_str = 'time'
   status = nf_put_att_text(ncid, varid_time, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   attr_str = 'days since 1901-01-01 00:00:00'
   status = nf_put_att_text(ncid, varid_time, 'units', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "units")

   attr_str = 'noleap'
   status = nf_put_att_text(ncid, varid_time, 'calendar', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "calendar")

   attr_str = 'time_bounds'
   status = nf_put_att_text(ncid, varid_time, 'bounds', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "bounds")


   ! Attributes of mcdate
   attr_str = 'current date (YYYYMMDD)'
   status = nf_put_att_text(ncid, varid_mcdate, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")


   !Attributes of area
   attr_str = 'grid cell areas'
   status = nf_put_att_text(ncid, varid_area, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   attr_str = 'km^2'
   status = nf_put_att_text(ncid, varid_area, 'units', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "units")

   status = nf_put_att_real(ncid, varid_area, '_FillValue', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_real(ncid, varid_area, 'missing_value', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "missing_value")


   !Attributes of topo
   attr_str = 'grid cell topography'
   status = nf_put_att_text(ncid, varid_topo, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   attr_str = 'm'
   status = nf_put_att_text(ncid, varid_topo, 'units', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "units")

   status = nf_put_att_real(ncid, varid_topo, '_FillValue', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_real(ncid, varid_topo, 'missing_value', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "missing_value")


   !Attributes of landfrac
   attr_str = 'land fraction'
   status = nf_put_att_text(ncid, varid_landfrac, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   status = nf_put_att_real(ncid, varid_landfrac, '_FillValue', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_real(ncid, varid_landfrac, 'missing_value', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "missing_value")


   !Attributes of landmask
   attr_str = 'land/ocean mask (0.=ocean and 1.=land)'
   status = nf_put_att_text(ncid, varid_landmask, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")
 
   status = nf_put_att_int(ncid, varid_landmask, '_FillValue', NF_INT, 1, MISSING_INT)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_int(ncid, varid_landmask, 'missing_value', NF_INT, 1, MISSING_INT)
   if (status /= nf_noerr) call handle_err(status, "missing_value")
 
 
   !Attributes of pftmask
   attr_str = 'pft real/fake mask (0.=fake and 1.=real)'
   status = nf_put_att_text(ncid, varid_pftmask, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")
 
   status = nf_put_att_int(ncid, varid_pftmask, '_FillValue', NF_INT, 1, MISSING_INT)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_int(ncid, varid_pftmask, 'missing_value', NF_INT, 1, MISSING_INT)
   if (status /= nf_noerr) call handle_err(status, "missing_value")


   ! Attributes of gpp
   attr_str = 'gross primary production'
   status = nf_put_att_text(ncid, varid_gpp, 'long_name', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "long_name")

   attr_str = 'gC m-2 day-1'
   status = nf_put_att_text(ncid, varid_gpp, 'units', len(trim(attr_str)), trim(attr_str))
   if (status /= nf_noerr) call handle_err(status, "units")

   status = nf_put_att_real(ncid, varid_gpp, '_FillValue', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
   status = nf_put_att_real(ncid, varid_gpp, 'missing_value', NF_REAL, 1, MISSING_VALUE)
   if (status /= nf_noerr) call handle_err(status, "missing_value")


   ! End the definition phase so that variables can be written to the file

   status = nf_enddef(ncid)
   if (status /= nf_noerr) call handle_err(status, "enddef")


   ! Write variable values to ncfilename

   status =  nf_put_var(ncid, varid_lon, lon1d)
   if (status /= nf_noerr) call handle_err(status, "lon")

   status =  nf_put_var(ncid, varid_lat, lat1d)
   if (status /= nf_noerr) call handle_err(status, "lat")

   status =  nf_put_var(ncid, varid_time, time)
   if (status /= nf_noerr) call handle_err(status, "time")

   status =  nf_put_var(ncid, varid_mcdate, mcdate)
   if (status /= nf_noerr) call handle_err(status, "mcdate")

   status =  nf_put_var(ncid, varid_area, area)
   if (status /= nf_noerr) call handle_err(status, "area")

   status =  nf_put_var(ncid, varid_topo, topo)
   if (status /= nf_noerr) call handle_err(status, "topo")

   status =  nf_put_var(ncid, varid_landfrac, landfrac)
   if (status /= nf_noerr) call handle_err(status, "landfrac")

   status =  nf_put_var(ncid, varid_landmask, landmask)
   if (status /= nf_noerr) call handle_err(status, "landmask")

   status =  nf_put_var(ncid, varid_pftmask, pftmask)
   if (status /= nf_noerr) call handle_err(status, "pftmask")

   status =  nf_put_var(ncid, varid_gpp, gpp)
   if (status /= nf_noerr) call handle_err(status, "gpp")


   status = nf_close(ncid)

   if (verbose .ge. 0) print *, 'Done writing GPP to file ', trim(ncfilename)

   end subroutine NcWriteDailyGPPOutput


end module canopy_model
