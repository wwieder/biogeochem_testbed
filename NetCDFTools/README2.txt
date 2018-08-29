module load compiler/gnu/8.1.0
module load tool/netcdf/4.6.1/gcc


Update 4/17/2017:

Don't allow the difference between field capacity and wilting point to be < 5%

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

----------------------------------------------------------------------------------------------------
Update 2/20/2017

Use SOILLIQ(time, levgrnd, lat, lon) (kg/m2) instead of H2OSOI(time, levgrnd, lat, lon) (mm3/mm3) 
to compute volumetric soil water content for the met.nc files.

Here is the conversion from kg/m2 to volumetric swc (mm3/mm3):

kg    1000 cm3       m2             1             m      cm3     mm3
--- * --------  * ---------  * -------------- * ----- =  ---  =  ---
m2       kg       100^2 cm2    thickness (m)    100 cm   cm3     mm3


Also read in SOILPSI(time, levgrnd, lat, lon).



NetCDF variable metadata:

/project/tss/wwieder/CASACLM/clm_forcing/CRU_hist/clm4_5_12_r191_CLM45spHIST_CRU.clm2.h1.1850-01-01-00000.nc

netcdf clm4_5_12_r191_CLM45spHIST_CRU.clm2.h1.1850-01-01-00000 {
dimensions:
	lon = 144 ;
	lat = 96 ;
	levgrnd = 15 ;
	time = UNLIMITED ; // (365 currently)

variables:
	float levgrnd(levgrnd) ;
		levgrnd:long_name = "coordinate soil levels" ;
		levgrnd:units = "m" ;

	float H2OSOI(time, levgrnd, lat, lon) ;
		H2OSOI:long_name = "volumetric soil water (vegetated landunits only)" ;
		H2OSOI:units = "mm3/mm3" ;
		H2OSOI:cell_methods = "time: mean" ;
		H2OSOI:_FillValue = 1.e+36f ;
		H2OSOI:missing_value = 1.e+36f ;

	float SOILLIQ(time, levgrnd, lat, lon) ;
		SOILLIQ:long_name = "soil liquid water (vegetated landunits only)" ;
		SOILLIQ:units = "kg/m2" ;
		SOILLIQ:cell_methods = "time: mean" ;
		SOILLIQ:_FillValue = 1.e+36f ;
		SOILLIQ:missing_value = 1.e+36f ;

	float SOILPSI(time, levgrnd, lat, lon) ;
		SOILPSI:long_name = "soil water potential in each soil layer" ;
		SOILPSI:units = "MPa" ;
		SOILPSI:cell_methods = "time: mean" ;
		SOILPSI:_FillValue = 1.e+36f ;
		SOILPSI:missing_value = 1.e+36f ;


----------------------------------------------------------------------------------------------------
Updates from 2/20/2017 to 3/19/2017

Read CLM variable ??? and calculate xfrznmoist (volumentric frozen soil water content) in addition to 
xmoist (volumetric liquid soil water content).

Code from clm_netcdf_tools.f90:
            do k = 1, nlevsoi 
               !if (h2osoi(j,i,k,day) .lt. MISSING_VALUE) then
               if (soilliq(j,i,k,day) .lt. MISSING_VALUE) then
                  !!h2osoiday(j,i,k,day) = h2osoi(j,i,k,day)
                  ! Compute volumentric soil water content (mm3/mm3) from soil liquid water (kg/m2) -mdh 2/20/2017
                  h2osoiday(j,i,k,day) = soilliq(j,i,k,day) * 0.001 / levgrndthickness(k)
                  ! Compute the frozen fraction also. -mdh 3/13/2017
                  h2ofrznsoiday(j,i,k,day) = h2osoi(j,i,k,day) - h2osoiday(j,i,k,day)

Updates from WWeider 3/17-3/19:
The subroutine hydraulicProperties is used to compute volumetric soil water content for wilting point, 
field capacity, and saturation. WWeider addded WATSAT and WATFC to the CLM surface data file, and now 
these variables are output to grid_info_soil.csv.  Wilting point is still calculated internally in 
subroutine hydraulicProperties.
