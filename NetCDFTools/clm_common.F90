module clm_common

!-----------------------------------------------------------------------
!BOP
!
!  !MODULE: clm_common
!
!  !DESCRIPTION:
!  !  Contains constants and subroutines common to multiple modules
!
!  !USES:
   use shr_kind_mod
!
!  !PUBLIC TYPES:
   implicit none
   include 'netcdf.inc'
   public :: handle_err             ! Subroutine to print error messages related to netCDF reads/writes
   public :: getYearMonthDayHour    ! Return the calendar year, day of year, month, day of month, and hour of day associated with a time from CLM h0 or h1 file.
   public :: getTime                ! Return the time (days since 1901-01-01 00:00:00) associated with the calendar year, day of year, and hour of day.
   public :: getmcdate              ! Return the mcdate (YYYYMMDD) associated with a time (days since 1901-01-01 00:00:00) read from a CLM h0 or h1 file.
   public :: hydraulicProperties    ! Compute volumetric soil water content for wilting point, field capacity, and saturation

   real(4), parameter :: SEC_PER_HOUR = 3600
   real(4), parameter :: SEC_PER_DAY = 86400
   real(4), parameter :: HOURS_PER_DAY = 24
   real(4), parameter :: DAYS_PER_YEAR = 365
   real(4), parameter :: MISSING_VALUE = 1.e+36
   integer, parameter :: MISSING_INT = -9999
   integer, parameter :: MAXCO2YRS = 400	! maximum number of years in transient CO2 file
   real(4), parameter :: KELVIN = 273.15	! Equivalent of 0 degrees Celsius
   real(4), parameter :: PI =  3.14159265
   integer, parameter :: baseyear = 1901	! Times in the climate files are days since 1901-01-01 00:00:00
   real(4), parameter :: ORGANIC_MAX = 130	! Maximum soil organic matter density, for peat (kg/m3)
   real(4), parameter :: DENH2O = 1000.0        ! Density of liquid water (kg/m3)
   real(4), parameter :: GRAV = 9.80665         ! Gravitational acceleration (m/s2)
   real(4), parameter :: MOL_PER_UMOL = 1.0e-6  ! Moles per micromole

!! integer, parameter :: LATINDX = 70		! Testing location (40.74N, near Fort Collins, CO)
!! integer, parameter :: LONINDX = 103		! Testing location (255, near Fort Collins, CO)
!! integer, parameter :: LATINDX = 49		! Testing location (near equator, Amazon, South America)
!! integer, parameter :: LONINDX = 117		! Testing location (70W/290, Amazon, South America)
   integer, parameter :: LATINDX = 71		! Testing location (42.43N, Harvard Forest)
   integer, parameter :: LONINDX = 116		! Testing location (-72.19W/287.81, Harvard Forest)

   integer, parameter :: doWriteHourlyClimateTestFile = 0   	! Write hourly climate text file for LATINDX, LONINDX (unit=11)
   integer, parameter :: doWriteDailyClimateTestFile = 0	! Write daily climate text file for LATINDX, LONINDX (unit=12)
   integer, parameter :: doWriteDailyGPPTestFile = 1		! Write daily GPP text file for LATINDX, LONINDX (unit=16)

   ! verbose controls the amount of text to stdout: -1 = none, 0 = progress, 1 = debug1, 2 = debug2
   integer :: verbose = 1  
   logical :: h1Daily = .TRUE.	! True if h1 files are daily, they are hourly otherwise
   integer :: wrtSingleCells    ! =1 if write output to single cell files instead of to gridded files
   integer, allocatable :: wrtcell(:,:)     	! wrtcell(nlon,nlat) save individual cell results

   save

!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine handle_err(status, errmsg)
      integer, intent(in) :: status
      character(len=*), intent(in)    :: errmsg  ! append error message
     
      if (status /= nf_noerr) then
         print *, trim(nf_strerror(status)), ": ", errmsg
         stop "Stopped"
      endif
   end subroutine handle_err

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getYearMonthDayHour
!
! !INTERFACE:
   subroutine getYearMonthDayHour(time, calyear, dayofyr, month, dayofmo, hour)

! !DESCRIPTION:
!  Return the calendar year, day of year (1-365), month, day of month (1..31), and hour of day (0-23)  
!  associated with a time (days since 1901-01-01 00:00:00) from an h0 or h1 file.
!  Assumes no leap years.

   !ARGUMENTS
   real(4), intent(in)  :: time 		! days since 1901-01-01 00:00:00
   integer, intent(out) :: calyear		! 4-digit calendar year
   integer, intent(out) :: dayofyr	        ! day of year (1..365)
   integer, intent(out) :: month	        ! month (1..12)
   integer, intent(out) :: dayofmo	        ! day of mo (1..31)
   integer, intent(out) :: hour	        	! hour of the day (1..24), hour 1 is 1:00am, hour 24 is midnight of the next day

   !LOCAL VARIABLES
   character*2 :: filetype 	! 'h0' or 'h1'
   integer :: nYears    	! number of years since base year
   real(4) :: timeleft 		! intermediate calculation to compute day of year
   real(4) :: timeadj		! adjusted time to handle h1 file exceptions (see below)
   integer :: lastday(1:12)     ! day of year corresponding to the last day of the month 
   
   filetype = 'h1'		! filetype can be either 'h1' or 'h0', the time calculation works for both.
 
   calyear = 0
   dayofyr = 0
   month = 0
   dayofmo = 0
   hour = 0

   ! For h1 (hourly) files
   !   Ignore the first time slice of file 1901.
   !   The last hour in 1901 is found in the first timeslice of the 1902 file.  
   !   The last hour in 2000 is found in the first timeslice of the 2001 file.  Similarly for other years.
   if ((filetype .eq. 'h0') .or. (filetype .eq. 'h1')) then
      timeadj = time - (1.0/HOURS_PER_DAY)		! Subtract 1 hour from time (as fraction of a day)
      nYears = INT(timeadj / DAYS_PER_YEAR)
      calyear = baseyear + nYears 
      timeleft = timeadj - nYears * DAYS_PER_YEAR 
      dayofyr = INT(timeleft) + 1
      hour = NINT(HOURS_PER_DAY * (timeleft - (dayofyr-1))) + 1
   endif

   lastday(1) = 31
   lastday(2) = 59
   lastday(3) = 90
   lastday(4) = 120
   lastday(5) = 151
   lastday(6) = 181
   lastday(7) = 212
   lastday(8) = 243
   lastday(9) = 273
   lastday(10) = 304
   lastday(11) = 334
   lastday(12) = 365

   if (dayofyr .le. lastday(1)) then
      month = 1
      dayofmo = dayofyr
   else if (dayofyr .le. lastday(2)) then
      month = 2
      dayofmo = dayofyr - lastday(1)
   else if (dayofyr .le. lastday(3)) then
      month = 3
      dayofmo = dayofyr - lastday(2)
   else if (dayofyr .le. lastday(4)) then
      month = 4
      dayofmo = dayofyr - lastday(3)
   else if (dayofyr .le. lastday(5)) then
      month = 5
      dayofmo = dayofyr - lastday(4)
   else if (dayofyr .le. lastday(6)) then
      month = 6
      dayofmo = dayofyr - lastday(5)
   else if (dayofyr .le. lastday(7)) then
      month = 7
      dayofmo = dayofyr - lastday(6)
   else if (dayofyr .le. lastday(8)) then
      month = 8
      dayofmo = dayofyr - lastday(7)
   else if (dayofyr .le. lastday(9)) then
      month = 9
      dayofmo = dayofyr - lastday(8)
   else if (dayofyr .le. lastday(10)) then
      month = 10
      dayofmo = dayofyr - lastday(9)
   else if (dayofyr .le. lastday(11)) then
      month = 11
      dayofmo = dayofyr - lastday(10)
   else
      month = 12
      dayofmo = dayofyr - lastday(11)
   endif
   
   end subroutine getYearMonthDayHour

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getTime
!
! !INTERFACE:
   subroutine getTime(calyear, day, hour, time)
!
! !DESCRIPTION:
!  Return the time (days since 1901-01-01 00:00:00) associated with the calendar year, 
!  day of year (1-365), and hour of day (1-24) 
!
!  For h1 (hourly) files
!  Ignore the first time slice of file 1901.
!  The last hour in 1901 is found in the first timeslice of the 1902 file.  
!  The last hour in 2000 is found in the first timeslice of the 2001 file.  Similarly for other years.

   !ARGUMENTS
   integer, intent(in) :: calyear	! 4-digit calendar year
   integer, intent(in) :: day	        ! day of year (1..365)
   integer, intent(in) :: hour	        ! hour of the day (1..24), hour 1 is 1:00am, hour 24 is midnight of the next day
   real(4), intent(out) :: time 	! days since 1901-01-01 00:00:00

   !LOCAL VARIABLES
   integer :: nYears    	! number of years since base year
   real(4) :: timeleft 		! intermediate calculation to compute day of year
   real(4) :: timeadj		! adjusted time to handle h1 file exceptions (see below)
   
 
   nYears = calyear - baseyear
   time = nYears * DAYS_PER_YEAR + (day - 1) + hour * (1.0/HOURS_PER_DAY)

   end subroutine getTime

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getmcdate
!
! !INTERFACE:
   integer function getmcdate(time)

! !DESCRIPTION:
!  Return the mcdate (YYYYMMDD) associated 
!  with a time (days since 1901-01-01 00:00:00) read from an h0 or h1 file.

   !ARGUMENTS
   real(4), intent(in) :: time 			! days since 1901-01-01 00:00:00


   !LOCAL VARIABLES
   integer :: calyear		! 4-digit calendar year
   integer :: dayofyr	        ! day of year (1..365)
   integer :: month	        ! month (1..12)
   integer :: dayofmo	        ! day of year (1..31)
   integer :: hour	        ! hour of the day (1..24), hour 1 is 1:00am, hour 24 is midnight of the next day

   call getYearMonthDayHour(time, calyear, dayofyr, month, dayofmo, hour)

   getmcdate = calyear * 10000 + month * 100 + dayofmo

   end function getmcdate

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_time_and_date
!
! !INTERFACE:
   subroutine get_time_and_date(date_string, time_string)

! !DESCRIPTION:

      !ARGUMENTS
      character*10, intent(out) :: date_string
      character*8, intent(out) :: time_string

      !LOCAL VARIABLES
!     character(8)  :: date
!     character(10) :: time
!     character(5)  :: zone
      character*4 :: syear
      character*2 :: smonth, sday
      character*2 :: shour, smin, ssec
   
      integer,dimension(8) :: values
!     call date_and_time(date,time,zone,values)
!     call date_and_time(DATE=date,ZONE=zone)
!     call date_and_time(TIME=time)
!     print '(a,2x,a,2x,a)', date, time, zone
!     print '(8i5))', values

      call date_and_time(VALUES=values)

      write(syear,  '(i4)') values(1)
      write(smonth, '(i2)') values(2)
      write(sday,   '(i2)') values(3)
      write(shour,  '(i2)') values(5)
      write(smin,   '(i2)') values(6)
      write(ssec,   '(i2)') values(7)

      if (values(2) .lt. 10) smonth(1:1) = '0'
      if (values(3) .lt. 10) sday(1:1) = '0'
      if (values(5) .lt. 10) shour(1:1) = '0'
      if (values(6) .lt. 10) smin(1:1) = '0'
      if (values(7) .lt. 10) ssec(1:1) = '0'
      
      date_string(1:2) = smonth
      date_string(3:3) = "/"
      date_string(4:5) = sday
      date_string(6:6) = "/"
      date_string(7:10) = syear

      time_string(1:2) = shour
      time_string(3:3) = ":"
      time_string(4:5) = smin
      time_string(6:6) = ":"
      time_string(7:8) = ssec

!     date_string = trim(smonth) // '/' // trim(sday) '/' // trim(syear)
!     time_string = trim(shour) // ':' // trim(smin) ':' // trim(ssec)

!     print *, date_string
!     print *, time_string

   end subroutine get_time_and_date

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: hydraulicProperties
!
! !INTERFACE:
!WWsubroutine hydraulicProperties(sand,clay,silt,wwilt,wfield,wsat)
   subroutine hydraulicProperties(sand,clay,silt,wwilt,wfield)

! !DESCRIPTION:
!  Given soil texture properties (%sand, %clay, %silt) compute
!  volumetric percentages for wilting point, field capacity, and saturation.
!  The algorithms are based on Jabloun and Sahli 2006:
!     "Development and comparative analysis of pedotransfer functions for 
!      predicting soil water characteristic content for Tunisian soil.
!      Proceedings of the 7th Edition of TJASSST 2006."
!  I selected this method because I could calculate the 3 soil properties
!  from sand, clay, and silt alone. 
!  Consider adding Saxton equations also.

   !ARGUMENTS
   real(4), intent(in) :: sand 		! percent sand (%)
   real(4), intent(in) :: clay 		! percent clay (%)
   real(4), intent(in) :: silt 		! percent silt (%)
   real(4), intent(in) :: wfield        ! WW added so wwilt !< wfield
   real(4), intent(out) :: wwilt 	! wilting point (%)
! WW
!   real(4), intent(out) :: wfield 	! field capacity (%)
!   real(4), intent(out) :: wsat 	! saturation (%)
!   real(4)  :: wfield       ! field capacity (%)
!   real(4)  :: wsat         ! saturation (%)
 
!   wsat = -99.9
!   wfield = -99.9
   wwilt = -99.9

   ! Check for missing values.  No error.
   if ((sand .gt. 100) .or. (silt .gt. 100) .or. (clay .gt. 100)) then
      return
   endif

   if (abs(sand + clay + silt - 100.0) .gt. 0.01) then
      print *, "WARNING in hydraulicProperties: sand + silt + clay <> 100 ", sand+clay+silt
      return
   endif

! WW commented
!   wsat    = 0.6658*silt + 0.1567*sand - 0.0079*silt*silt - 12.31121/sand  &
!             - 6.4756*alog(sand) - 0.0038*clay*silt + 0.0038*clay*sand &
!             - 0.0042*silt*sand + 52.7526 
!   wfield  = 118.932*clay + 119.0866*silt + 119.1104*sand + 162.31731/clay  &
!             - 46.21921/silt-5.12991/sand  + 18.1733*alog(clay) + 0.0013*clay*silt &
!             + 0.0022*silt*sand - 11939.3493
   wwilt   = - 1.5722*silt - 0.5423*sand - 0.0072*clay*clay + 0.0072*silt*silt  &
             - 0.0059*sand*sand + 160.14591/clay  +  6.60011/sand  &
             + 0.0022*clay*silt - 0.0039*clay*sand + 92.3851 

   ! Found a bug when silt is high (sand, clay, silt) = (29.398, 1.009, 69.593). 
   !   wwilt = 155.744%
   !   wfiel = 134.11%
   !   wsat =  34.374%	
   ! If bounds on calculations aren't correct, use average values. -MDH 8/25/2014

!  if (wwilt .gt. 30.0) then
!      wwilt = 15.0
!  end if
! WW added above to be consistent below

   if ((wwilt .ge. wfield*100) .or. (wfield*100 - wwilt .le. 5.0)) then
      ! Keep minimum difference between field capacity and wilting point. -mdh 4/17/2017
      !wwilt = (wfield*100) - 10.0
      wwilt = (wfield*100) - 5.0
   endif
! WW added to keep wwilt<wfield

   if (wwilt .le. 0.0) then
      !This condition (wfield <= 5%) is unlikely. -mdh 4/17/2017
      write(*,*) 'WARNING in hydraulicProperties: field capacity % below excepted range:', wfield*100
      write(*,*) '  Resetting wilting point % to ', wwilt
      wwilt = 0.5*(wfield*100)
   end if

!  if ((wwilt .gt. wfield) .or. (wfield .gt. wsat)) then
!     wsat = 43.0
!     wfield = 32.0
!     wwilt = 15.0
!  endif
!  if ((wwilt .gt. 30.0) .or. (wfield .gt. 45.0) .or. (wsat .gt. 50.0)) then
!     wsat = 43.0
!     wfield = 32.0
!     wwilt = 15.0
!  endif

   end subroutine hydraulicProperties

!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PutVariableAttributeReal
!
! !INTERFACE:
   subroutine PutVariableAttributeReal(ncid, varid, attr_name, attr_units, missing_value)

! !DESCRIPTION:
  !  Write variable attributes for a real variable to NetCDF file

!   ARGUMENTS
      integer, intent(in) :: ncid			! netcdf file id
      integer, intent(in) :: varid			! netcdf variable id
      character*100, intent(in) :: attr_name		! String for assigning variable "long_name" attribute
      character*100, intent(in) :: attr_units		! String for assigning variable "units" attribute
      real(4), intent(in) :: missing_value
      
!   LOCAL VARIABLES
      integer status	! NetCDF error status

      status = nf_put_att_text(ncid, varid, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
 
      status = nf_put_att_text(ncid, varid, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
 
      status = nf_put_att_real(ncid, varid, '_FillValue', NF_REAL, 1, missing_value)
      if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
      status = nf_put_att_real(ncid, varid, 'missing_value', NF_REAL, 1, missing_value)
      if (status /= nf_noerr) call handle_err(status, "missing_value")

   end subroutine PutVariableAttributeReal

end module clm_common
