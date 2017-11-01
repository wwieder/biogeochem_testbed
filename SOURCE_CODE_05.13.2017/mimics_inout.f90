!--------------------------------------------------------------------------------
! FILE: mimics_inout.f90
!
! Purpose: 
!   Input/Output subroutines for the MIMICS model
!
!     SUBROUTINE mimics_readbiome - read pftlookup_mimics.csv (mimics parameters)
!     SUBROUTINE mimics_init - initialize mimics pools. If initcasa >= 1, read .csv file 
!     SUBROUTINE mimics_poolfluxout - write mimics pools to restart .csv output file
!     SUBROUTINE WritePoolFluxNcFile_mimics_annual - write annual mimics pools and fluxes to netCDF file
!     SUBROUTINE WritePoolFluxNcFile_mimics_daily - write daily mimics pools and fluxes to netCDF file
!
! Contact: Melannie Hartman
!          melannie@ucar.edu
!
! History:
!   1/5/2015 - Created
!--------------------------------------------------------------------------------

SUBROUTINE mimics_readbiome(fname_mimicsbiome, mp, mvtype)
  use casadimension
  use define_types
  use casavariable
  use mimicsdimension
  use mimicsparam
  use mimicsvariable
  implicit none

  ! Function arguments
  character(len=100), INTENT(IN) :: fname_mimicsbiome 
  integer, INTENT(IN) :: mp     ! number of grid points
  integer, INTENT(IN) :: mvtype ! number of vegetation (biome) types

  ! Local variables
  integer   :: npt, nP, nv, nv1
  real(r_2) :: Pscalar, kmod_val(6)

  open(101,file=fname_mimicsbiome)

  !Skip past vegetation descriptions
  read(101,*) 
  do nv=1,mvtype
      read(101,*) 
  end do

  ! Fixed parameters
  ! TO DO - check input file variable name to assure corrert placement of parameter values
  read(101,*) 
 
  read(101,*) mimicsbiome%Vslope(R1) 
  read(101,*) mimicsbiome%Vslope(R2) 
  read(101,*) mimicsbiome%Vslope(R3)
  read(101,*) mimicsbiome%Vslope(K1)
  read(101,*) mimicsbiome%Vslope(K2)
  read(101,*) mimicsbiome%Vslope(K3) 

  read(101,*) mimicsbiome%Vint(R1) 
  read(101,*) mimicsbiome%Vint(R2) 
  read(101,*) mimicsbiome%Vint(R3)
  read(101,*) mimicsbiome%Vint(K1)
  read(101,*) mimicsbiome%Vint(K2)
  read(101,*) mimicsbiome%Vint(K3) 

  read(101,*) mimicsbiome%av(R1) 
  read(101,*) mimicsbiome%av(R2) 
  read(101,*) mimicsbiome%av(R3)
  read(101,*) mimicsbiome%av(K1)
  read(101,*) mimicsbiome%av(K2)
  read(101,*) mimicsbiome%av(K3) 

  read(101,*) mimicsbiome%Kslope(R1) 
  read(101,*) mimicsbiome%Kslope(R2) 
  read(101,*) mimicsbiome%Kslope(R3)  
  read(101,*) mimicsbiome%Kslope(K1)  
  read(101,*) mimicsbiome%Kslope(K2) 
  read(101,*) mimicsbiome%Kslope(K3) 

  read(101,*) mimicsbiome%Kint(R1) 
  read(101,*) mimicsbiome%Kint(R2) 
  read(101,*) mimicsbiome%Kint(R3)  
  read(101,*) mimicsbiome%Kint(K1)  
  read(101,*) mimicsbiome%Kint(K2) 
  read(101,*) mimicsbiome%Kint(K3)

  read(101,*) mimicsbiome%ak(R1) 
  read(101,*) mimicsbiome%ak(R2) 
  read(101,*) mimicsbiome%ak(R3)
  read(101,*) mimicsbiome%ak(K1)
  read(101,*) mimicsbiome%ak(K2)
  read(101,*) mimicsbiome%ak(K3) 

  read(101,*) mimicsbiome%Vmod(R1)
  read(101,*) mimicsbiome%Vmod(R2)
  read(101,*) mimicsbiome%Vmod(R3)
  read(101,*) mimicsbiome%Vmod(K1)
  read(101,*) mimicsbiome%Vmod(K2)
  read(101,*) mimicsbiome%Vmod(K3)

  read(101,*) kmod_val(R1)
  read(101,*) kmod_val(R2)
  read(101,*) kmod_val(R3)
  read(101,*) kmod_val(K1)
  read(101,*) kmod_val(K2)
  read(101,*) kmod_val(K3)

  read(101,*) mimicsbiome%KO(1)
  read(101,*) mimicsbiome%KO(2)

  read(101,*) mimicsbiome%MGE(1)
  read(101,*) mimicsbiome%MGE(2)
  read(101,*) mimicsbiome%MGE(3)
  read(101,*) mimicsbiome%MGE(4)

  read(101,*) mimicsbiome%tau_r(1)
  read(101,*) mimicsbiome%tau_r(2)

  read(101,*) mimicsbiome%tau_k(1)
  read(101,*) mimicsbiome%tau_k(2)

  read(101,*) mimicsbiome%tauModDenom
  read(101,*) mimicsbiome%tauMod_MIN
  read(101,*) mimicsbiome%tauMod_MAX

  read(101,*) mimicsbiome%fPHYS_r(1)
  read(101,*) mimicsbiome%fPHYS_r(2)

  read(101,*) mimicsbiome%fPHYS_K(1)
  read(101,*) mimicsbiome%fPHYS_K(2)

  read(101,*) mimicsbiome%fCHEM_r(1)
  read(101,*) mimicsbiome%fCHEM_r(2)
  read(101,*) mimicsbiome%fCHEM_r(3)

  read(101,*) mimicsbiome%fCHEM_K(1)
  read(101,*) mimicsbiome%fCHEM_K(2)
  read(101,*) mimicsbiome%fCHEM_K(3)

  read(101,*) mimicsbiome%fSOM_p(1)
  read(101,*) mimicsbiome%fSOM_p(2)

  read(101,*) mimicsbiome%phys_scalar(1)
  read(101,*) mimicsbiome%phys_scalar(2)

  read(101,*) mimicsbiome%Fi(metbc)
  read(101,*) mimicsbiome%Fi(struc)

  read(101,*) mimicsbiome%fmet_p(1)
  read(101,*) mimicsbiome%fmet_p(2)
  read(101,*) mimicsbiome%fmet_p(3)

! write(*,*) 'mimicsbiome%MGE(1) =', mimicsbiome%MGE(1)
! write(*,*) 'mimicsbiome%MGE(2) =', mimicsbiome%MGE(2)
! write(*,*) 'mimicsbiome%MGE(3) =', mimicsbiome%MGE(3)
! write(*,*) 'mimicsbiome%MGE(4) =', mimicsbiome%MGE(4)

  write(*,*) 'mimicsbiome%fmet_p(1) =', mimicsbiome%fmet_p(1)
  write(*,*) 'mimicsbiome%fmet_p(2) =', mimicsbiome%fmet_p(2)
  write(*,*) 'mimicsbiome%fmet_p(3) =', mimicsbiome%fmet_p(3)


  !Biome-specific parameters
  read(101,*) 
  read(101,*) 
  do nv=1,mvtype
     read(101,*) nv1,mimicsbiome%depth(nv)
     if (mimicsbiome%depth(nv) <= 0.0) then
        write(*,*) 'Error in ', trim(filename_mimicsbiome),' (depth <= 0): depth(',nv,')=',mimicsbiome%depth(nv)
        STOP
     endif
  end do

  close(101)

  ! Calculate cell-specific parameters
  do npt=1,mp

      ! Use site-level value for Pscalar (-mdh 4/20/2015)
      ! Pscalar = 1.0 / (2.0 * exp(-2.0*SQRT(soil%clay(npt))))
      ! Use global value for Pscalar (-mdh 4/6/2015)
      ! Pscalar = 1.0 / (0.8 * exp(-3.0*SQRT(soil%clay(npt))))
      ! Updated based on Will's (-mdh 6/1/2015)
      ! Pscalar = 1.0 / (2.0 * exp(-3.0*SQRT(soil%clay(npt))))

      Pscalar = mimicsbiome%phys_scalar(1) * exp(mimicsbiome%phys_scalar(2)*SQRT(soil%clay(npt)))
      mimicsbiome%Kmod(npt,R1) = kmod_val(R1)           ! modifies Km[r1] for fluxes from LITm to MICr 
      mimicsbiome%Kmod(npt,R2) = kmod_val(R2)           ! modifies Km[r2] for fluxes from LITs to MICr 
      mimicsbiome%Kmod(npt,R3) = kmod_val(R3) * Pscalar ! modifies Km[r3] for fluxes from SOMa to MICr 
      mimicsbiome%Kmod(npt,K1) = kmod_val(K1)           ! modifies Km[k1] for fluxes from LITm to MICk 
      mimicsbiome%Kmod(npt,K2) = kmod_val(K2)           ! modifies Km[k2] for fluxes from LITs MICk 
      mimicsbiome%Kmod(npt,K3) = kmod_val(K3) * Pscalar ! modifies Km[k3] for fluxes from SOMa to MICk 

      do nP=1,mplant
          ! fmet = fraction of plant residue transferred to metabolic litter (0.0 - 1.0)
          ! It is computed from the biome-specific lignin:N ratio for plant pool nP 
          !   (g lignin / g N) =  (g lignin / g C) / (g N / g C) 
          mimicsbiome%ligninNratio(npt,nP) =  casabiome%fracLigninPlant(veg%iveg(npt), nP) &
                                              / max(0.001,casapool%ratioNCplant(npt,nP))
          ! write(*,*) 'lignin:N = ', mimicsbiome%ligninNratio(npt,nP)
      end do


      ! mimicsbiome%desorb(npt) - desorbsion rate from SOMp to SOMa (hr-1)
      ! mimicsbiome%desorb(npt) = 1.5 * 0.00001 * exp(-1.5 * soil%clay(npt))
      mimicsbiome%desorb(npt) = mimicsbiome%fSOM_p(1) * exp(mimicsbiome%fSOM_p(2) * soil%clay(npt))

  end do

END SUBROUTINE mimics_readbiome

!--------------------------------------------------------------------------------

SUBROUTINE mimics_init(filename_mimicsipool,mp,ms,mst)
! Initialize mimics litter, microbe, and SOM pools
! Reset casa litter and SOM pools to zero because they aren't needed by the MIMICS model
! This subroutine should be called after subroutine init_casa (not in place of)

  use casadimension
  use casaparm
  use casavariable
  use mimicsdimension
  use mimicsparam
  use mimicsvariable
  implicit none

  !Subroutine arguments
  character(len=100), INTENT(IN) :: filename_mimicsipool
  integer,            INTENT(IN) :: mp, ms, mst
  
  !Local Variables
  integer   :: np,npt,npz,nl,ns,nland,nlandz
  real(r_2) :: nyearz,ivtz,latz,lonz,areacellz

  print *, 'initial MIMICS pool file: ',filename_mimicsipool

  mimicspool%LITm(:) = 0.0
  mimicspool%LITs(:) = 0.0 
  mimicspool%MICr(:) = 0.0 
  mimicspool%MICk(:) = 0.0 
  mimicspool%SOMa(:) = 0.0 
  mimicspool%SOMc(:) = 0.0 
  mimicspool%SOMp(:) = 0.0 

  WHERE(casamet%iveg2 /= icewater)
      mimicspool%LITm(:) = 1.0
      mimicspool%LITs(:) = 1.0
      mimicspool%MICr(:) = 0.015
      mimicspool%MICk(:) = 0.025
      mimicspool%SOMa(:) = 1.0
      mimicspool%SOMc(:) = 1.0
      mimicspool%SOMp(:) = 1.0
  END WHERE

  !If not a spinup run (initcasa .ne. 0) read initial pool values from file
  if (initcasa >= 1) then
      open(105,file=filename_mimicsipool)
      read(105,*)  ! Skip past header line
      do npt =1, mp
          read(105,*) nyearz,npz,ivtz,latz,lonz,areacellz,         &
                      mimicspool%LITm(npt), mimicspool%LITs(npt),  &
                      mimicspool%MICr(npt), mimicspool%MICk(npt),  &
                      mimicspool%SOMa(npt), mimicspool%SOMc(npt), mimicspool%SOMp(npt)

          !ATTENTION: Check npz, ivtz, latz, lonz, areacellz against values read by casa_init
          !TO DO
      end do
      close(105)
  endif

  ! check for negative pool sizes
  mimicspool%LITm = max(0.0, mimicspool%LITm)
  mimicspool%LITs = max(0.0, mimicspool%LITs)
  mimicspool%MICr = max(0.0, mimicspool%MICr)
  mimicspool%MICk = max(0.0, mimicspool%MICk)
  mimicspool%SOMa = max(0.0, mimicspool%SOMa)
  mimicspool%SOMc = max(0.0, mimicspool%SOMc)
  mimicspool%SOMp = max(0.0, mimicspool%SOMp)

! ATTENTION: need corresponding assignments for MIMICS?
! casabal%clitterlast = casapool%clitter
! casabal%csoillast   = casapool%csoil
! casabal%sumcbal     = 0.0

! ATTENTION: check if these assignments will create divide by zero errors
! I commented these assignments out to see if it makes a difference in
! the output. (-MDH 2/16/2015)
! casapool%Csoil(:,:)   = 0.0
! casapool%Clitter(:,:) = 0.0
! casapool%Nsoil(:,:)   = 0.0
! casapool%Psoil(:,:)   = 0.0
! casapool%Nsoilmin(:)  = 0.0

END SUBROUTINE mimics_init

!--------------------------------------------------------------------------------

SUBROUTINE mimics_poolfluxout(filename_mimicsepool,mp,iYrCnt,myear,writeToRestartCSVfile)

  use define_types
! use casadimension
  use casaparm
  use casavariable
  use mimicsdimension
  use mimicsparam
  use mimicsvariable
  implicit none

  !Subroutine arguments
  character(len=100), INTENT(IN) :: filename_mimicsepool
  integer, INTENT(IN)            :: mp, iYrCnt, myear
  logical, INTENT(IN)            :: writeToRestartCSVfile

  !Local Variables
  integer   :: npt,nout,nso
  real(r_2) :: xyear, xyear2

!--------------------------------------------------------------------------
! Calculate average daily pool over the years and average annual fluxes
! These pools and fluxes were accumulated each day of the simulation.
! See subroutine mimics_caccum in mimics_cycle.f90.

! xyear=1.0/(real(myear)*365)
! xyear2 = 1.0/(real(myear))

  ! Output is every year now, not every myear years. -mdh 10/17/2016
  xyear=1.0/365.0
  xyear2 = 1.0

  !Divide by number of simulation days to get average daily pool value
  mimicspoolAn%ClitterAn(:,:)  = mimicspoolAn%ClitterAn(:,:)  * xyear
  mimicspoolAn%CmicrobeAn(:,:) = mimicspoolAn%CmicrobeAn(:,:) * xyear
  mimicspoolAn%CsoilAn(:,:)    = mimicspoolAn%CsoilAn(:,:)    * xyear

  !Divide by number of simulation years to get average annual flux
  mimicsfluxAn%ChrespAn(:)      = mimicsfluxAn%ChrespAn(:) * xyear2    
  mimicsfluxAn%CLitInputAn(:,:) = mimicsfluxAn%CLitInputAn(:,:) * xyear2   

!--------------------------------------------------------------------------

  if (writeToRestartCSVfile) then

      nout=106
      open(nout,file=filename_mimicsepool)

      ! mimicsbal%sumcbal=min(9999.0,max(-9999.0,mimicsbal%sumcbal))

      !Write header line
      write(nout,91) "iYrCnt,npt,veg%iveg,",                &
                     "casamet%lat,casamet%lon,",          &
                     "casamet%areacell,",                 &
                     "mimicspool%LITm,mimicspool%LITs,",  &
                     "mimicspool%MICr,mimicspool%MICk,",  &
                     "mimicspool%SOMa,mimicspool%SOMc,",  &
                     "mimicspool%SOMp"
    91 format(a20,a24,a17,a32,a32,a32,a15)

      do npt =1, mp

! Attention...If initial pool values for MIMICS transient run are not set to 0.0 for icewater, 
! daily global mean output drops at the end of the first year when these pools are set to 0.0 below. 
! Rely on initialization to set icewater cells correctly since MIMICS does not simulate them. -mdh 7/11/2016
!         if (casamet%iveg2(npt) == icewater) then
!             mimicspool%LITm(npt) = 0.0
!             mimicspool%LITs(npt) = 0.0
!             mimicspool%MICr(npt) = 0.0
!             mimicspool%MICk(npt) = 0.0
!             mimicspool%SOMa(npt) = 0.0
!             mimicspool%SOMc(npt) = 0.0
!             mimicspool%SOMp(npt) = 0.0
!!            mimicsbal%sumcbal(npt) = 0.0 
!         endif

          write(nout,92) iYrCnt,npt,veg%iveg(npt),            &
                         casamet%lat(npt),casamet%lon(npt), &
                         casamet%areacell(npt)*(1.0e-9),    &
                         mimicspool%LITm(npt), mimicspool%LITs(npt),  &
                         mimicspool%MICr(npt), mimicspool%MICk(npt),  &
                         mimicspool%SOMa(npt), mimicspool%SOMc(npt), mimicspool%SOMp(npt)
      end do

      CLOSE(nout)

  endif

92  format(3(i6,',',2x),10(f18.10,',',2x))

!--------------------------------------------------------------------------
  nout=107
  open(nout,file="mimics_diagnostic.csv")

  !Write header line
  write(nout,93) "iYrCnt,npt,veg%iveg,",                &
                 "casamet%lat,casamet%lon,",          &
                 "casamet%areacell,",                 &
                 "lignin:N,fmet,"
93 format(a20,a24,a17,a13)


  do npt =1, mp

      write(nout,94) iYrCnt,npt,veg%iveg(npt),            &
                     casamet%lat(npt),casamet%lon(npt), &
                     casamet%areacell(npt)*(1.0e-9),    &
                     mimicsbiome%ligninNratioAvg(npt), mimicsbiome%fmet(npt)
  end do

  CLOSE(nout)

94  format(3(i6,',',2x),5(f18.10,',',2x))

END SUBROUTINE mimics_poolfluxout

!----------------------------------------------------------------------------------------------------
    SUBROUTINE WritePoolFluxNcFile_mimics_annual(filename_ncOut, mp, year)

!   DESCRIPTION
!   Define and write average annual output variables from MIMICS model to netcdf file filename_ncOut.
!   Called by casacnpdriver once a year and again at the end of the simulation.
!   The file will contain values for a single year.
!
!   Melannie Hartman. January 26, 2015 
 
    USE define_types
    USE casaparm
    USE casavariable
    USE clmgridvariable
    USE mimicsdimension
    USE mimicsparam
    USE mimicsvariable
    implicit none
    include 'netcdf.inc'
    real(4), parameter :: MISSING_VALUE = 1.e+36
    integer, parameter :: MISSING_INT = -9999

!   ARGUMENTS
      character(len=*), intent(in) :: filename_ncOut	! NetCDF output file name 
      integer, intent(in) :: mp			        ! number of points with valid data
      integer, intent(in) :: year			! output year
!   LOCAL VARIABLES
      integer :: i
      integer :: ncid				! netcdf file ID
      integer :: status				! function return status
!     integer :: dimid_mp			! netcdf dimension id
      integer :: dimid_lat			! netcdf dimension id
      integer :: dimid_lon			! netcdf dimension id
      integer :: dimid_time			! netcdf dimension id
      integer :: nlon, nlat, ntimes	 	! Dimension sizes for NetCDf file
      integer :: dims(3)			! Array of NetCDF dimension IDs for defining variables
      integer :: start1(1), count1(1)           ! start and count arrays for writing 1-D data from netcdf files
      integer :: start2(2), count2(2)           ! start and count arrays for writing 2-D data from netcdf files
      integer :: start3(3), count3(3)           ! start and count arrays for writing 3-D data from netcdf files
      integer :: varid_lon, varid_lat		! NetCDF variable ID for latitude and longitude
      integer :: varid_time			! NetCDF variable ID for time
!     integer :: varid_year			! NetCDF variable ID for year
      integer :: varid_mask			! NetCDF variable ID for cellMissing(nlon,nlat)
      integer :: varid_cellid			! NetCDF variable ID for cellid(nlon,nlat)
      integer :: varid_igbp			! NetCDF variable ID for IGBP_PFT(nlon,nlat)
      integer :: varid_landarea			! NetCDF variable ID for land area
      integer :: varid_cLITm			! NetCDF variable ID for metablic litter C
      integer :: varid_cLITs			! NetCDF variable ID for structural litter C
      integer :: varid_cMICr			! NetCDF variable ID for microbe r-selected pool C
      integer :: varid_cMICk			! NetCDF variable ID for microbe k-selected pool C
      integer :: varid_cSOMa			! NetCDF variable ID for available soil C
      integer :: varid_cSOMc			! NetCDF variable ID for chemically protected soil C
      integer :: varid_cSOMp			! NetCDF variable ID for physically protected soil C
      integer :: varid_cHresp			! NetCDF variable ID for soil heterotrophic respiration C
      integer :: varid_cLitIn_m			! NetCDF variable ID for metabolic litter Input C
      integer :: varid_cLitIn_s			! NetCDF variable ID for structural litter Input C
      character*100 :: attr_name		! String for assigning global and variable attributes
      character*100 :: attr_units		! String for assigning global and variable attributes
      character*10 :: date_string		! String for assigning date to global attributes
      character*8 :: time_string		! String for assigning time to global attributes
      integer :: verbose=0
      integer :: npt, ilon, ilat, itime
      integer, allocatable :: IGBP_PFT(:,:)	! IGBP_PFT(nlon,nlat) IGBP PFT classification (1-18)
      real(4), allocatable :: landarea(:,:)	! landarea(nlon,nlat) km^2
      real(4), allocatable :: var1(:,:,:)	! gridded output variable
      real(4), allocatable :: var2(:,:,:)	! gridded output variable
      real(4), allocatable :: var3(:,:,:)	! gridded output variable
      real(4), allocatable :: var4(:,:,:)	! gridded output variable
      real(4), allocatable :: var5(:,:,:)	! gridded output variable
      real(4), allocatable :: var6(:,:,:)	! gridded output variable
      real(4), allocatable :: var7(:,:,:)	! gridded output variable
      real(4), allocatable :: var8(:,:,:)	! gridded output variable
      real(4), allocatable :: var9(:,:,:)	! gridded output variable
      real(4), allocatable :: var10(:,:,:)	! gridded output variable
      real(4) :: time                          ! time in years (1..ntimes)

      if (verbose .ge. 0) print *, "Writing output to file ", trim(filename_ncOut), "..."

! Output these variables in NetCDF file:
!     veg%iveg(npt)				- IGBP PFTs
!     casamet%lat(npt)				- latitudes	
!     casamet%lon(npt)				- longitudes
!     casamet%areacell(npt)*(1.0e-6) 		- landarea (km^2)
!   Average annual values for fluxes and pools      
!     mimicsfluxAn%CLitInputAn(npt,METBC)	- metabolic litter inputs (gC/m2/yr)
!     mimicsfluxAn%CLitInputAn(npt,STRUC)	- structural litter inputs (gC/m2/yr)
!     mimicsfluxAn%ChrespAn(npt)		- heterotrphic respiration flux (gC/m2/yr)
!     mimicspoolAn%ClitterAn(npt,METBC)		- metabolic litter pool (gC/m2)
!     mimicspoolAn%ClitterAn(npt,STRUC)		- structural litter pool (gC/m2)
!     mimicspoolAn%CmicrobeAn(npt,RSEL) 	- r-selected microbe pool (gC/m2)
!     mimicspoolAn%CmicrobeAn(npt,KSEL) 	- k-selected microbe pool (gC/m2)
!     mimicspoolAn%CsoilAn(npt,AVAL)		- available SOM pool (gC/m2)
!     mimicspoolAn%CsoilAn(npt,CHEM)		- chemically protected SOM pool (gC/m2)
!     mimicspoolAn%CsoilAn(npt,PHYS)		- physically protected SOM pool (gC/m2)

   dims(1) = 0
   dims(2) = 0
   ntimes = 1
   nlat = clmgrid%nlat
   nlon = clmgrid%nlon

   ! Create the netcdf file 

   status = nf_create(filename_ncOut, NF_CLOBBER, ncid)
   if (status /= nf_noerr) call handle_err(status, trim(filename_ncOut))


   ! Define file dimensions

   status = nf_def_dim(ncid, 'lon', nlon, dimid_lon)
   if (status /= nf_noerr) call handle_err(status, "lon")

   status = nf_def_dim(ncid, 'lat', nlat, dimid_lat)
   if (status /= nf_noerr) call handle_err(status, "lat")

!  status = nf_def_dim(ncid, 'time', ntimes, dimid_time)
   status = nf_def_dim(ncid, 'time', nf_unlimited, dimid_time)
   if (status /= nf_noerr) call handle_err(status, "time")

!  status = nf_def_dim(ncid, 'mp', mp, dimid_mp)
!  if (status /= nf_noerr) call handle_err(status, "mp")


   ! Define variables

   dims(1) = dimid_lon
   status = nf_def_var(ncid, 'lon', NF_REAL, 1, dims, varid_lon)
   if (status /= nf_noerr) call handle_err(status, "def_var(lon)")

   dims(1) = dimid_lat
   status = nf_def_var(ncid, 'lat', NF_REAL, 1, dims, varid_lat)
   if (status /= nf_noerr) call handle_err(status, "def_var(lat)")

   dims(1) = dimid_time
   status = nf_def_var(ncid, 'time', NF_FLOAT, 1, dims, varid_time)
   if (status /= nf_noerr) call handle_err(status, "def_var(time)")

!  dims(1) = dimid_time
!  status = nf_def_var(ncid, 'year', NF_INT, 1, dims, varid_year)
!  if (status /= nf_noerr) call handle_err(status, "def_var(year)")

   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions of FORTRAN arrays 
   ! are in the opposite order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat

   status = nf_def_var(ncid, 'IGBP_PFT', NF_INT, 2, dims, varid_igbp)
   if (status /= nf_noerr) call handle_err(status, "IGBP_PFT")

   status = nf_def_var(ncid, 'landarea', NF_REAL, 2, dims, varid_landarea)
   if (status /= nf_noerr) call handle_err(status, "landarea")

   status = nf_def_var(ncid, 'cellMissing', NF_INT, 2, dims, varid_mask)
   if (status /= nf_noerr) call handle_err(status, "cellMissing")

   status = nf_def_var(ncid, 'cellid', NF_INT, 2, dims, varid_cellid)
   if (status /= nf_noerr) call handle_err(status, "cellid")


   ! Because dimensions in FORTRAN are in Column Major Order (the first 
   ! array index varies the most rapidly) dimensions of FORTRAN arrays 
   ! are in the opposite order that they appear in the NetCDF file with ncdump. 
   dims(1) = dimid_lon
   dims(2) = dimid_lat
   dims(3) = dimid_time

   status = nf_def_var(ncid, 'cLITm', NF_REAL, 3, dims, varid_cLITm)
   if (status /= nf_noerr) call handle_err(status, "cLITm")

   status = nf_def_var(ncid, 'cLITs', NF_REAL, 3, dims, varid_cLITs)
   if (status /= nf_noerr) call handle_err(status, "cLITs")

   status = nf_def_var(ncid, 'cMICr', NF_REAL, 3, dims, varid_cMICr)
   if (status /= nf_noerr) call handle_err(status, "cMICr")

   status = nf_def_var(ncid, 'cMICk', NF_REAL, 3, dims, varid_cMICk)
   if (status /= nf_noerr) call handle_err(status, "cMICk")

   status = nf_def_var(ncid, 'cSOMa', NF_REAL, 3, dims, varid_cSOMa)
   if (status /= nf_noerr) call handle_err(status, "cSOMa")

   status = nf_def_var(ncid, 'cSOMc', NF_REAL, 3, dims, varid_cSOMc)
   if (status /= nf_noerr) call handle_err(status, "cSOMc")

   status = nf_def_var(ncid, 'cSOMp', NF_REAL, 3, dims, varid_cSOMp)
   if (status /= nf_noerr) call handle_err(status, "cSOMp")

   status = nf_def_var(ncid, 'cHresp', NF_REAL, 3, dims, varid_cHresp)
   if (status /= nf_noerr) call handle_err(status, "cHresp")

  status = nf_def_var(ncid, 'cLitInput_metbc', NF_REAL, 3, dims, varid_cLitIn_m)
   if (status /= nf_noerr) call handle_err(status, "cLitInput_metbc")

  status = nf_def_var(ncid, 'cLitInput_struc', NF_REAL, 3, dims, varid_cLitIn_s)
   if (status /= nf_noerr) call handle_err(status, "cLitInput_struc")


   ! Global attributes
   attr_name = 'MIMICS model output'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'title', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "title")
 
   attr_name = 'NOTE: None of the variables are weighted by land fraction!'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'comment', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "comment")
 
   call get_time_and_date(date_string, time_string)
   attr_name = 'created on ' // date_string // ' ' // time_string
   status = nf_put_att_text(ncid, NF_GLOBAL, 'history', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "history")
 
   attr_name = 'MIMICS Model'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'source', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "source")

   attr_name = trim(filename_mimicsbiome)
   status = nf_put_att_text(ncid, NF_GLOBAL, 'parameters', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "parameters")

   ! ------------------------------ Attributes of the variables ------------------------------

   ! Attributes of time variable
   attr_name = 'coordinate time'
   status = nf_put_att_text(ncid, varid_time, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   !!attr_units = 'simulation time in years'
   attr_units = '1..ntimes'
   status = nf_put_att_text(ncid, varid_time, 'units', len(trim(attr_units)), trim(attr_units))
   if (status /= nf_noerr) call handle_err(status, "units")

!  ! Attributes of year variable
!  attr_name = 'calendar year'
!  status = nf_put_att_text(ncid, varid_year, 'long_name', len(trim(attr_name)), trim(attr_name))
!  if (status /= nf_noerr) call handle_err(status, "long_name")
!  attr_units = 'year'
!  status = nf_put_att_text(ncid, varid_year, 'units', len(trim(attr_units)), trim(attr_units))
!  if (status /= nf_noerr) call handle_err(status, "units")

   ! Attributes of lon variable
   attr_name = 'coordinate longitude'
   attr_units = 'degrees_east'
   !!call PutVariableAttributeReal(ncid, varid_lon, attr_name, attr_units, MISSING_VALUE)
   status = nf_put_att_text(ncid, varid_lon, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   status = nf_put_att_text(ncid, varid_lon, 'units', len(trim(attr_units)), trim(attr_units))
   if (status /= nf_noerr) call handle_err(status, "units")
 
   ! Attributes of lat variable
   attr_name = 'coordinate latitude'
   attr_units = 'degrees_north'
   !!call PutVariableAttributeReal(ncid, varid_lat, attr_name, attr_units, MISSING_VALUE)
   status = nf_put_att_text(ncid, varid_lat, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   status = nf_put_att_text(ncid, varid_lat, 'units', len(trim(attr_units)), trim(attr_units))
   if (status /= nf_noerr) call handle_err(status, "units")

   ! Attributes of IGBP_PFT variable
   attr_name = 'IGBP PFT classification'
   status = nf_put_att_text(ncid, varid_igbp, 'long_name', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "long_name")
   attr_units = '1-18'
   status = nf_put_att_text(ncid, varid_igbp, 'units', len(trim(attr_units)), trim(attr_units))
   if (status /= nf_noerr) call handle_err(status, "units")
   status = nf_put_att_int(ncid, varid_igbp, '_FillValue', NF_INT, 1, MISSING_INT)
   if (status /= nf_noerr) call handle_err(status, "_FillValue")
   status = nf_put_att_int(ncid, varid_igbp, 'missing_value', NF_INT, 1, MISSING_INT)
   if (status /= nf_noerr) call handle_err(status, "missing_value")

   ! Attributes of landarea variable
   attr_name = 'land area, icewater set to 0.0'
   attr_units = 'km^2'
   call PutVariableAttributeReal(ncid, varid_landarea, attr_name, attr_units, MISSING_VALUE)

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

   ! Attributes of cHresp variable
   attr_name = 'soil heterotrophic respiration'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cHresp, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLitInput_metbc variable
   attr_name = 'metabolic litter inputs'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cLitIn_m, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLitInput_struc variable
   attr_name = 'structural litter inputs'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cLitIn_s, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLITm variable
   attr_name = 'metabolic litter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cLITm, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLITs variable
   attr_name = 'structural litter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cLITs, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cMICr variable
   attr_name = 'r-selected microbial soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cMICr, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cMICk variable
   attr_name = 'k-selected microbial soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cMICk, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cSOMa variable
   attr_name = 'active soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cSOMa, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cSOMc variable
   attr_name = 'chemically protected soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cSOMc, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cSOMp variable
   attr_name = 'physically protected soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cSOMp, attr_name, attr_units, MISSING_VALUE)



   ! --------------- End the definition phase so that variables can be written to the file ---------------
   status = nf_enddef(ncid)
   if (status /= nf_noerr) call handle_err(status, "enddef")


   ! ------------------------------  Write variable values to filename_ncOut ------------------------------

   !! Write time using nf_put_vara_real instead of nf_put_var so
   !! that time gets written correctly. -mdh 12/28/2016

   time = real(year)
   start1 = (/ 1 /)
   count1 = (/ ntimes /)    ! ntimes = 1 in this subroutine

!! status =  nf_put_var(ncid, varid_time, time)
   status =  nf_put_vara_real(ncid, varid_time, start1, count1, time)
   if (status /= nf_noerr) call handle_err(status, "put_vara_real(time)")

!  status =  nf_put_var(ncid, varid_year, year)
!  if (status /= nf_noerr) call handle_err(status, "put_var(year)")
 
   status =  nf_put_var(ncid, varid_lon, clmgrid%lon1d)
   if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%lon1d)")
 
   status =  nf_put_var(ncid, varid_lat, clmgrid%lat1d)
   if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%lat1d)")
 
   status =  nf_put_var(ncid, varid_mask, clmgrid%cellMissing)
   if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%cellMissing)")
 
   status =  nf_put_var(ncid, varid_cellid, clmgrid%cellid)
   if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%cellid)")

!  Define start3 and count3 for record variables (those with unlimited time dimension)
!  Write up to 6 output variables at a time


   start3 = (/ 1, 1, 1 /)
   count3 = (/nlon, nlat, ntimes/)

   allocate(IGBP_PFT(1:clmgrid%nlon,1:clmgrid%nlat))
   allocate(landarea(1:clmgrid%nlon,1:clmgrid%nlat))
   allocate(var1(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var2(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var3(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var4(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var5(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var6(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var7(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var8(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var9(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var10(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))

!  Average annual values for fluxes and pools
!  mimicsfluxAn%CLitInputAn(npt,METBC)	- metabolic litter inputs (gC/m2/yr)
!  mimicsfluxAn%CLitInputAn(npt,STRUC)	- structural litter inputs (gC/m2/yr)
!  mimicsfluxAn%ChrespAn(npt)	 	- heterotrophic respiration flux (gC/m2/yr)
!  mimicspoolAn%ClitterAn(npt,METBC)	- metabolic litter pool (gC/m2)
!  mimicspoolAn%ClitterAn(npt,STRUC)	- structural litter pool (gC/m2)
!  mimicspoolAn%CmicrobeAn(npt,RSEL) 	- r-selected microbe pool (gC/m2)
!  mimicspoolAn%CmicrobeAn(npt,KSEL) 	- k-selected microbe pool (gC/m2)
!  mimicspoolAn%CsoilAn(npt,AVAL)	- available SOM pool (gC/m2)
!  mimicspoolAn%CsoilAn(npt,CHEM)	- chemically protected SOM pool (gC/m2)
!  mimicspoolAn%CsoilAn(npt,PHYS)	- physically protected SOM pool (gC/m2)

   IGBP_PFT(:,:) = MISSING_INT
   landarea(:,:) = MISSING_VALUE
   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE
   var8(:,:,:) = MISSING_VALUE
   var9(:,:,:) = MISSING_VALUE
   var10(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_mimics_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
!     else
!        print *, 'mimicspoolAn%ClitterAn(',npt,',METBC) = ', mimicspoolAn%ClitterAn(npt,METBC)
!        print *, 'mimicspoolAn%ClitterAn(',npt,',STRUC) = ', mimicspoolAn%ClitterAn(npt,STRUC)
!
!        print *, 'mimicspoolAn%CmicrobeAn(',npt,',RSEL) = ', mimicspoolAn%CmicrobeAn(npt,RSEL)
!        print *, 'mimicspoolAn%CmicrobeAn(',npt,',KSEL) = ', mimicspoolAn%CmicrobeAn(npt,KSEL)
!
!        print *, 'mimicspoolAn%CsoilAn(',npt,',AVAL) = ', mimicspoolAn%CsoilAn(npt,AVAL)
!        print *, 'mimicspoolAn%CsoilAn(',npt,',CHEM) = ', mimicspoolAn%CsoilAn(npt,CHEM)
!        print *, 'mimicspoolAn%CsoilAn(',npt,',PHYS) = ', mimicspoolAn%CsoilAn(npt,PHYS)
      endif

      IGBP_PFT(ilon,ilat) = veg%iveg(npt)
      !Set land area for icewater cells to 0.0 so global mean timeseries will be correct. -mdh 7/11/2016
      if (casamet%iveg2(npt) /= icewater) then
          landarea(ilon,ilat) = casamet%areacell(npt)*(1.0e-6)   !Convert m^2 to km^2
      else
          landarea(ilon,ilat) = 0.0
      endif
      var9(ilon,ilat,itime)  = mimicsfluxAn%CLitInputAn(npt,METBC)
      var10(ilon,ilat,itime) = mimicsfluxAn%CLitInputAn(npt,STRUC)
      var1(ilon,ilat,itime)  = mimicsfluxAn%ChrespAn(npt)
      var2(ilon,ilat,itime)  = mimicspoolAn%ClitterAn(npt,METBC)
      var3(ilon,ilat,itime)  = mimicspoolAn%ClitterAn(npt,STRUC)
      var4(ilon,ilat,itime)  = mimicspoolAn%CmicrobeAn(npt,RSEL)
      var5(ilon,ilat,itime)  = mimicspoolAn%CmicrobeAn(npt,KSEL)
      var6(ilon,ilat,itime)  = mimicspoolAn%CsoilAn(npt,AVAL)
      var7(ilon,ilat,itime)  = mimicspoolAn%CsoilAn(npt,CHEM)
      var8(ilon,ilat,itime)  = mimicspoolAn%CsoilAn(npt,PHYS)
   end do

   ! IGBP PFT
   status =  nf_put_var(ncid, varid_igbp, IGBP_PFT)
   if (status /= nf_noerr) call handle_err(status, "put_var(IGBP_PFT)")

   ! Land area 
   status =  nf_put_var(ncid, varid_landarea, landarea)
   if (status /= nf_noerr) call handle_err(status, "put_var(landarea)")

!! ATTENTION: I could not get the variables with a "time" dimension to write 
!! to the netcdf file when the time dimension was unlimited, UNLESS I substituted 
!! one "nf_put_vara_real" for a "nf_put_var".  I DON'T UNDERSTAND!
!! Otherwise nf_put_var seemed to ignore start3 and count3.
!! Melannie 6/3/2014
!! status =  nf_put_var(ncid, varid_cLitIn_m, var9, start3, count3)

   status =  nf_put_vara_real(ncid, varid_cLitIn_m, start3, count3, var9)
   if (status /= nf_noerr) call handle_err(status, "put_vara_real(cLitIn_m)")

   status =  nf_put_var(ncid, varid_cLitIn_s, var10, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cLitIn_s)")

   status =  nf_put_var(ncid, varid_cHresp, var1, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cHresp)")

   status =  nf_put_var(ncid, varid_cLITm, var2, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cLITm)")

   status =  nf_put_var(ncid, varid_cLITs, var3, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cLITs)")

   status =  nf_put_var(ncid, varid_cMICr, var4, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cMICr)")

   status =  nf_put_var(ncid, varid_cMICk, var5, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cMICk)")

   status =  nf_put_var(ncid, varid_cSOMa, var6, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cSOMa)")

   status =  nf_put_var(ncid, varid_cSOMc, var7, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cSOMc)")

   status =  nf_put_var(ncid, varid_cSOMp, var8, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cSOMp)")


   deallocate(var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,IGBP_PFT,landarea)

   status = nf_close(ncid)

   if (verbose .ge. 0) print *, "Done writing output to ", trim(filename_ncOut), "..."

END SUBROUTINE WritePoolFluxNcFile_mimics_annual


!----------------------------------------------------------------------------------------------------
    SUBROUTINE WritePoolFluxNcFile_mimics_daily(filename_ncOut, mp, year, iday)

!   DESCRIPTION
!   Write daily output variables from MIMICS model to netcdf file filename_ncOut
!   If iday=1, the file is create and define the file first.
!   Called once a day from casacnpdriver when mdaily=1
!   The file will contain output from a single year, 365 days.
!
!   Melannie Hartman. August 24, 2015 
 
    USE define_types
    USE casaparm
    USE casavariable
    USE clmgridvariable
    USE mimicsdimension
    USE mimicsparam
    USE mimicsvariable
    implicit none
    include 'netcdf.inc'
    real(4), parameter :: MISSING_VALUE = 1.e+36
    integer, parameter :: MISSING_INT = -9999

!   ARGUMENTS
      character(len=*), intent(in) :: filename_ncOut	! NetCDF output file name 
      integer, intent(in) :: mp			        ! number of points with valid data
      integer, intent(in) :: year			! output year
      integer, intent(in) :: iday			! output day (1..365)
!   LOCAL VARIABLES
      integer :: i
      integer :: ncid				! netcdf file ID
      integer :: status				! function return status
!     integer :: dimid_mp			! netcdf dimension id
      integer :: dimid_lat			! netcdf dimension id
      integer :: dimid_lon			! netcdf dimension id
      integer :: dimid_time			! netcdf dimension id
      integer :: nlon, nlat, ntimes	 	! Dimension sizes for NetCDf file
      integer :: nwrtimes                       ! Number of times that will be written when this subroutine is called
      integer :: dims(3)			! Array of NetCDF dimension IDs for defining variables
      integer :: start1(1), count1(1)           ! start and count arrays for writing 1-D data from netcdf files
      integer :: start2(2), count2(2)           ! start and count arrays for writing 2-D data from netcdf files
      integer :: start3(3), count3(3)           ! start and count arrays for writing 3-D data from netcdf files
      integer :: varid_lon, varid_lat		! NetCDF variable ID for latitude and longitude
      integer :: varid_time			! NetCDF variable ID for time
      integer :: varid_day			! NetCDF variable ID for day
      integer :: varid_mask			! NetCDF variable ID for cellMissing(nlon,nlat)
      integer :: varid_cellid			! NetCDF variable ID for cellid(nlon,nlat)
      integer :: varid_igbp			! NetCDF variable ID for IGBP_PFT(nlon,nlat)
      integer :: varid_landarea			! NetCDF variable ID for land area
      integer :: varid_cLITm			! NetCDF variable ID for metablic litter C
      integer :: varid_cLITs			! NetCDF variable ID for structural litter C
      integer :: varid_cMICr			! NetCDF variable ID for microbe r-selected pool C
      integer :: varid_cMICk			! NetCDF variable ID for microbe k-selected pool C
      integer :: varid_cSOMa			! NetCDF variable ID for available soil C
      integer :: varid_cSOMc			! NetCDF variable ID for chemically protected soil C
      integer :: varid_cSOMp			! NetCDF variable ID for physically protected soil C
      integer :: varid_cHresp			! NetCDF variable ID for soil heterotrophic respiration C
      integer :: varid_cLitIn_m			! NetCDF variable ID for metabolic litter Input C
      integer :: varid_cLitIn_s			! NetCDF variable ID for structural litter Input C
      character*100 :: attr_name		! String for assigning global and variable attributes
      character*100 :: attr_units		! String for assigning global and variable attributes
      character*10 :: date_string		! String for assigning date to global attributes
      character*8 :: time_string		! String for assigning time to global attributes
      integer :: verbose=0
      integer :: npt, ilon, ilat, itime
      integer, allocatable :: IGBP_PFT(:,:)	! IGBP_PFT(nlon,nlat) IGBP PFT classification (1-18)
      integer, allocatable :: days(:)          	! day array (1..ntimes)
      real(4), allocatable :: time(:)          	! time array (time in years)
      real(4), allocatable :: landarea(:,:)	! landarea(nlon,nlat) km^2
      real(4), allocatable :: var1(:,:,:)	! gridded output variable
      real(4), allocatable :: var2(:,:,:)	! gridded output variable
      real(4), allocatable :: var3(:,:,:)	! gridded output variable
      real(4), allocatable :: var4(:,:,:)	! gridded output variable
      real(4), allocatable :: var5(:,:,:)	! gridded output variable
      real(4), allocatable :: var6(:,:,:)	! gridded output variable
      real(4), allocatable :: var7(:,:,:)	! gridded output variable
      real(4), allocatable :: var8(:,:,:)	! gridded output variable
      real(4), allocatable :: var9(:,:,:)	! gridded output variable
      real(4), allocatable :: var10(:,:,:)	! gridded output variable
      real(r_2) :: unitConv		        ! mgC/cm3 * depth(cm)* (1g/10^3mg)*(10^4cm2)/m2 = gC/m2

      if (verbose .gt. 0) print *, iday, "Writing output to file ", trim(filename_ncOut), "..."

! Output these variables in NetCDF file:
!     veg%iveg(npt)				- IGBP PFTs
!     casamet%lat(npt)				- latitudes	
!     casamet%lon(npt)				- longitudes
!     casamet%areacell(npt)*(1.0e-6) 		- landarea (km^2)
!   Daily values for fluxes and pools      
!     mimicsflux%ClitInput(npt,metbc)	- metabolic litter inputs (mgC/cm3/day)
!     mimicsflux%ClitInput(npt,struc)	- structural litter inputs (mgC/cm3/day)
!     mimicsflux%Chresp(npt)		- heterotrphic respiration flux (mgC/cm3/day)
!     mimicspool%LITm(npt)   	        - metabolic litter pool (mgC/cm3)
!     mimicspool%LITs(npt)	        - structural litter pool (mgC/cm3)
!     mimicspool%MICr(npt) 	        - r-selected microbe pool (mgC/cm3)
!     mimicspool%MICk(npt) 	        - k-selected microbe pool (mgC/cm3)
!     mimicspool%SOMa(npt)		- available SOM pool (mgC/cm3)
!     mimicspool%SOMc(npt)		- chemically protected SOM pool (mgC/cm3)
!     mimicspool%SOMp(npt)		- physically protected SOM pool (mgC/cm3)

   dims(1) = 0
   dims(2) = 0
   ntimes = 365
   nlat = clmgrid%nlat
   nlon = clmgrid%nlon

   if (iday == 1) then

      ! Create the netcdf file 
   
      status = nf_create(filename_ncOut, NF_CLOBBER, ncid)
      if (status /= nf_noerr) call handle_err(status, trim(filename_ncOut))
   
   
      ! Define file dimensions
   
      status = nf_def_dim(ncid, 'lon', nlon, dimid_lon)
      if (status /= nf_noerr) call handle_err(status, "lon")
   
      status = nf_def_dim(ncid, 'lat', nlat, dimid_lat)
      if (status /= nf_noerr) call handle_err(status, "lat")
   
   !  status = nf_def_dim(ncid, 'time', ntimes, dimid_time)
      status = nf_def_dim(ncid, 'time', nf_unlimited, dimid_time)
      if (status /= nf_noerr) call handle_err(status, "time")
   
   !  status = nf_def_dim(ncid, 'mp', mp, dimid_mp)
   !  if (status /= nf_noerr) call handle_err(status, "mp")
   
   
      ! Define variables
   
      dims(1) = dimid_lon
      status = nf_def_var(ncid, 'lon', NF_REAL, 1, dims, varid_lon)
      if (status /= nf_noerr) call handle_err(status, "def_var(lon)")
   
      dims(1) = dimid_lat
      status = nf_def_var(ncid, 'lat', NF_REAL, 1, dims, varid_lat)
      if (status /= nf_noerr) call handle_err(status, "def_var(lat)")
   
      dims(1) = dimid_time
      status = nf_def_var(ncid, 'time', NF_REAL, 1, dims, varid_time)
      if (status /= nf_noerr) call handle_err(status, "def_var(time)")
   
      dims(1) = dimid_time
      status = nf_def_var(ncid, 'day', NF_INT, 1, dims, varid_day)
      if (status /= nf_noerr) call handle_err(status, "def_var(day)")
   
      ! Because dimensions in FORTRAN are in Column Major Order (the first 
      ! array index varies the most rapidly) dimensions of FORTRAN arrays 
      ! are in the opposite order that they appear in the NetCDF file with ncdump. 
      dims(1) = dimid_lon
      dims(2) = dimid_lat
   
      status = nf_def_var(ncid, 'IGBP_PFT', NF_INT, 2, dims, varid_igbp)
      if (status /= nf_noerr) call handle_err(status, "IGBP_PFT")
   
      status = nf_def_var(ncid, 'landarea', NF_REAL, 2, dims, varid_landarea)
      if (status /= nf_noerr) call handle_err(status, "landarea")
   
      status = nf_def_var(ncid, 'cellMissing', NF_INT, 2, dims, varid_mask)
      if (status /= nf_noerr) call handle_err(status, "cellMissing")
   
      status = nf_def_var(ncid, 'cellid', NF_INT, 2, dims, varid_cellid)
      if (status /= nf_noerr) call handle_err(status, "cellid")
   
   
      ! Because dimensions in FORTRAN are in Column Major Order (the first 
      ! array index varies the most rapidly) dimensions of FORTRAN arrays 
      ! are in the opposite order that they appear in the NetCDF file with ncdump. 
      dims(1) = dimid_lon
      dims(2) = dimid_lat
      dims(3) = dimid_time
   
      status = nf_def_var(ncid, 'cLITm', NF_REAL, 3, dims, varid_cLITm)
      if (status /= nf_noerr) call handle_err(status, "cLITm")
   
      status = nf_def_var(ncid, 'cLITs', NF_REAL, 3, dims, varid_cLITs)
      if (status /= nf_noerr) call handle_err(status, "cLITs")
   
      status = nf_def_var(ncid, 'cMICr', NF_REAL, 3, dims, varid_cMICr)
      if (status /= nf_noerr) call handle_err(status, "cMICr")
   
      status = nf_def_var(ncid, 'cMICk', NF_REAL, 3, dims, varid_cMICk)
      if (status /= nf_noerr) call handle_err(status, "cMICk")
   
      status = nf_def_var(ncid, 'cSOMa', NF_REAL, 3, dims, varid_cSOMa)
      if (status /= nf_noerr) call handle_err(status, "cSOMa")
   
      status = nf_def_var(ncid, 'cSOMc', NF_REAL, 3, dims, varid_cSOMc)
      if (status /= nf_noerr) call handle_err(status, "cSOMc")
   
      status = nf_def_var(ncid, 'cSOMp', NF_REAL, 3, dims, varid_cSOMp)
      if (status /= nf_noerr) call handle_err(status, "cSOMp")
   
      status = nf_def_var(ncid, 'cHresp', NF_REAL, 3, dims, varid_cHresp)
      if (status /= nf_noerr) call handle_err(status, "cHresp")
   
      status = nf_def_var(ncid, 'cLitInput_metbc', NF_REAL, 3, dims, varid_cLitIn_m)
      if (status /= nf_noerr) call handle_err(status, "cLitInput_metbc")
   
      status = nf_def_var(ncid, 'cLitInput_struc', NF_REAL, 3, dims, varid_cLitIn_s)
      if (status /= nf_noerr) call handle_err(status, "cLitInput_struc")
   
   
      ! Global attributes
      attr_name = 'MIMICS model output'
      status = nf_put_att_text(ncid, NF_GLOBAL, 'title', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "title")
    
      attr_name = 'NOTE: None of the variables are weighted by land fraction!'
      status = nf_put_att_text(ncid, NF_GLOBAL, 'comment', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "comment")
    
      call get_time_and_date(date_string, time_string)
      attr_name = 'created on ' // date_string // ' ' // time_string
      status = nf_put_att_text(ncid, NF_GLOBAL, 'history', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "history")
    
      attr_name = 'MIMICS Model'
      status = nf_put_att_text(ncid, NF_GLOBAL, 'source', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "source")
 
      attr_name = trim(filename_mimicsbiome)
      status = nf_put_att_text(ncid, NF_GLOBAL, 'parameters', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "parameters")
 
      ! ------------------------------ Attributes of the variables ------------------------------
   
      ! Attributes of time variable
      attr_name = 'coordinate time'
      status = nf_put_att_text(ncid, varid_time, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
      !!attr_units = 'simulation time in years'
      attr_units = '1..ntimes'
      status = nf_put_att_text(ncid, varid_time, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
   
      ! Attributes of day variable
      attr_name = 'day of year'
      status = nf_put_att_text(ncid, varid_day, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
      attr_units = '1..365'
      status = nf_put_att_text(ncid, varid_day, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
   
      ! Attributes of lon variable
      attr_name = 'coordinate longitude'
      attr_units = 'degrees_east'
      !!call PutVariableAttributeReal(ncid, varid_lon, attr_name, attr_units, MISSING_VALUE)
      status = nf_put_att_text(ncid, varid_lon, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
      status = nf_put_att_text(ncid, varid_lon, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
    
      ! Attributes of lat variable
      attr_name = 'coordinate latitude'
      attr_units = 'degrees_north'
      !!call PutVariableAttributeReal(ncid, varid_lat, attr_name, attr_units, MISSING_VALUE)
      status = nf_put_att_text(ncid, varid_lat, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
      status = nf_put_att_text(ncid, varid_lat, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
   
      ! Attributes of IGBP_PFT variable
      attr_name = 'IGBP PFT classification'
      status = nf_put_att_text(ncid, varid_igbp, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
      attr_units = '1-18'
      status = nf_put_att_text(ncid, varid_igbp, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
      status = nf_put_att_int(ncid, varid_igbp, '_FillValue', NF_INT, 1, MISSING_INT)
      if (status /= nf_noerr) call handle_err(status, "_FillValue")
      status = nf_put_att_int(ncid, varid_igbp, 'missing_value', NF_INT, 1, MISSING_INT)
      if (status /= nf_noerr) call handle_err(status, "missing_value")
   
      ! Attributes of landarea variable
      attr_name = 'land area, icewater set to 0.0'
      attr_units = 'km^2'
      call PutVariableAttributeReal(ncid, varid_landarea, attr_name, attr_units, MISSING_VALUE)
   
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
   
      ! Attributes of cHresp variable
      attr_name = 'soil heterotrophic respiration'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cHresp, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cLitInput_metbc variable
      attr_name = 'metabolic litter inputs'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cLitIn_m, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cLitInput_struc variable
      attr_name = 'structural litter inputs'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cLitIn_s, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cLITm variable
      attr_name = 'metabolic litter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cLITm, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cLITs variable
      attr_name = 'structural litter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cLITs, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cMICr variable
      attr_name = 'r-selected microbial soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cMICr, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cMICk variable
      attr_name = 'k-selected microbial soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cMICk, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cSOMa variable
      attr_name = 'active soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cSOMa, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cSOMc variable
      attr_name = 'chemically protected soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cSOMc, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cSOMp variable
      attr_name = 'physically protected soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cSOMp, attr_name, attr_units, MISSING_VALUE)
   
   
   
      ! --------------- End the definition phase so that variables can be written to the file ---------------
      status = nf_enddef(ncid)
      if (status /= nf_noerr) call handle_err(status, "enddef")
   

      ! ------------------------------  Write 2-D variable values to filename_ncOut ------------------------------
   
      allocate(time(ntimes))
      allocate(days(ntimes))
      do i = 1, ntimes
         days(i) = i
         time(i) = real(year) + real(i)/365.0
      enddo
   
      !! Write time using nf_put_vara_real instead of nf_put_var so
      !! that time gets written correctly. -mdh 12/28/2016

      start1 = (/ 1 /)
      count1 = (/ ntimes /)    ! ntimes = number of days in the year

      !!status =  nf_put_var(ncid, varid_time, time)
      status =  nf_put_vara_real(ncid, varid_time, start1, count1, time)
      if (status /= nf_noerr) call handle_err(status, "put_vara_real(time)")
   
      !!status =  nf_put_var(ncid, varid_day, days)
      status =  nf_put_vara_int(ncid, varid_day, start1, count1, days)
      if (status /= nf_noerr) call handle_err(status, "put_vara_int(days)")
    
      status =  nf_put_var(ncid, varid_lon, clmgrid%lon1d)
      if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%lon1d)")
    
      status =  nf_put_var(ncid, varid_lat, clmgrid%lat1d)
      if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%lat1d)")
    
      status =  nf_put_var(ncid, varid_mask, clmgrid%cellMissing)
      if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%cellMissing)")
    
      status =  nf_put_var(ncid, varid_cellid, clmgrid%cellid)
      if (status /= nf_noerr) call handle_err(status, "put_var(clmgrid%cellid)")

      allocate(IGBP_PFT(1:clmgrid%nlon,1:clmgrid%nlat))
      allocate(landarea(1:clmgrid%nlon,1:clmgrid%nlat))

      IGBP_PFT(:,:) = MISSING_INT
      landarea(:,:) = MISSING_VALUE

      do npt = 1, mp
         ilon = casamet%ilon(npt)
         ilat = casamet%ilat(npt)
         if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
            print *, 'WritePoolFluxNcFile_mimics_daily: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
            print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
            STOP
         endif

         IGBP_PFT(ilon,ilat) = veg%iveg(npt)
         !Set land area for icewater cells to 0.0 so global mean timeseries will be correct. -mdh 7/11/2016
         if (casamet%iveg2(npt) /= icewater) then
             landarea(ilon,ilat) = casamet%areacell(npt)*(1.0e-6)   !Convert m^2 to km^2
         else
             landarea(ilon,ilat) = 0.0
         endif
      end do

      ! IGBP PFT
      status =  nf_put_var(ncid, varid_igbp, IGBP_PFT)
      if (status /= nf_noerr) call handle_err(status, "put_var(IGBP_PFT)")
   
      ! Land area 
      status =  nf_put_var(ncid, varid_landarea, landarea)
      if (status /= nf_noerr) call handle_err(status, "put_var(landarea)")

      deallocate(IGBP_PFT,landarea,time,days)

   else

      status = nf_open(filename_ncOut, NF_WRITE, ncid)
      if (status /= nf_noerr) call handle_err(status, "open(filename_ncOut)")

      ! get variable ids when iday > 1

      status = nf_inq_varid(ncid, 'cLITm',varid_cLITm)
      if (status /= nf_noerr) call handle_err(status, "cLITm")
   
      status = nf_inq_varid(ncid, 'cLITs', varid_cLITs)
      if (status /= nf_noerr) call handle_err(status, "cLITs")
   
      status = nf_inq_varid(ncid, 'cMICr', varid_cMICr)
      if (status /= nf_noerr) call handle_err(status, "cMICr")
   
      status = nf_inq_varid(ncid, 'cMICk', varid_cMICk)
      if (status /= nf_noerr) call handle_err(status, "cMICk")
   
      status = nf_inq_varid(ncid, 'cSOMa', varid_cSOMa)
      if (status /= nf_noerr) call handle_err(status, "cSOMa")
   
      status = nf_inq_varid(ncid, 'cSOMc', varid_cSOMc)
      if (status /= nf_noerr) call handle_err(status, "cSOMc")
   
      status = nf_inq_varid(ncid, 'cSOMp', varid_cSOMp)
      if (status /= nf_noerr) call handle_err(status, "cSOMp")
   
      status = nf_inq_varid(ncid, 'cHresp', varid_cHresp)
      if (status /= nf_noerr) call handle_err(status, "cHresp")
   
      status = nf_inq_varid(ncid, 'cLitInput_metbc', varid_cLitIn_m)
      if (status /= nf_noerr) call handle_err(status, "cLitInput_metbc")
   
      status = nf_inq_varid(ncid, 'cLitInput_struc', varid_cLitIn_s)
      if (status /= nf_noerr) call handle_err(status, "cLitInput_struc")

   endif      ! if (iday == 1)

!  Define start3 and count3 for record variables (those with unlimited time dimension)
!  Write up to 10 output variables at a time

   nwrtimes = 1

   start3 = (/ 1, 1, iday /)
   count3 = (/nlon, nlat, nwrtimes/)

   allocate(var1(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var2(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var3(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var4(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var5(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var6(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var7(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var8(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var9(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var10(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))

!  Daily values for fluxes and pools
!  mimicsflux%ClitInput(npt,metbc)	- metabolic litter inputs (mgC/cm3/day)
!  mimicsflux%ClitInput(npt,struc)	- structural litter inputs (mgC/cm3/day)
!  mimicsflux%Chresp(npt)		- heterotrphic respiration flux (mgC/cm3/day)
!  mimicspool%LITm(npt)   	        - metabolic litter pool (mgC/cm3)
!  mimicspool%LITs(npt)	        	- structural litter pool (mgC/cm3)
!  mimicspool%MICr(npt) 	        - r-selected microbe pool (mgC/cm3)
!  mimicspool%MICk(npt) 	        - k-selected microbe pool (mgC/cm3)
!  mimicspool%SOMa(npt)			- available SOM pool (mgC/cm3)
!  mimicspool%SOMc(npt)			- chemically protected SOM pool (mgC/cm3)
!  mimicspool%SOMp(npt)			- physically protected SOM pool (mgC/cm3)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE
   var8(:,:,:) = MISSING_VALUE
   var9(:,:,:) = MISSING_VALUE
   var10(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      unitConv = 10.0 * mimicsbiome%depth(veg%iveg(npt))     ! mgC/cm3 * depth(cm)* (1g/10^3mg)*(10^4cm2)/m2 = gC/m2
      var9(ilon,ilat,itime)  = mimicsflux%CLitInput(npt,METBC) * unitConv
      var10(ilon,ilat,itime) = mimicsflux%CLitInput(npt,STRUC) * unitConv
      var1(ilon,ilat,itime)  = mimicsflux%Chresp(npt) * unitConv
      var2(ilon,ilat,itime)  = mimicspool%LITm(npt)   * unitConv
      var3(ilon,ilat,itime)  = mimicspool%LITs(npt)   * unitConv
      var4(ilon,ilat,itime)  = mimicspool%MICr(npt)   * unitConv
      var5(ilon,ilat,itime)  = mimicspool%MICk(npt)   * unitConv
      var6(ilon,ilat,itime)  = mimicspool%SOMa(npt)   * unitConv
      var7(ilon,ilat,itime)  = mimicspool%SOMc(npt)   * unitConv
      var8(ilon,ilat,itime)  = mimicspool%SOMp(npt)   * unitConv
   end do


!! ATTENTION: I could not get the variables with a "time" dimension to write 
!! to the netcdf file when the time dimension was unlimited, UNLESS I substituted 
!! one "nf_put_vara_real" for a "nf_put_var".  I DON'T UNDERSTAND!
!! Otherwise nf_put_var seemed to ignore start3 and count3.
!! Melannie 6/3/2014
!! Had to use "nf_put_vara_real" for all variables to avoid segmentation violation
!! Melannie 8/31/2015

!! status =  nf_put_var(ncid, varid_cLitIn_m, var9, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cLitIn_m, start3, count3, var9)
   if (status /= nf_noerr) call handle_err(status, "put_vara_real(cLitIn_m)")

!! status =  nf_put_var(ncid, varid_cLitIn_s, var10, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cLitIn_s, start3, count3, var10)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cLitIn_s)")

!! status =  nf_put_var(ncid, varid_cHresp, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cHresp, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cHresp)")

!! status =  nf_put_var(ncid, varid_cLITm, var2, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cLITm, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cLITm)")

!! status =  nf_put_var(ncid, varid_cLITs, var3, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cLITs, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cLITs)")

!! status =  nf_put_var(ncid, varid_cMICr, var4, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cMICr, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cMICr)")

!! status =  nf_put_var(ncid, varid_cMICk, var5, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cMICk, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cMICk)")

!! status =  nf_put_var(ncid, varid_cSOMa, var6, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cSOMa, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cSOMa)")

!! status =  nf_put_var(ncid, varid_cSOMc, var7, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cSOMc, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cSOMc)")

!! status =  nf_put_var(ncid, varid_cSOMp, var8, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cSOMp, start3, count3, var8)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cSOMp)")


   deallocate(var1,var2,var3,var4,var5,var6,var7,var8,var9,var10)

   status = nf_close(ncid)

   if (verbose .gt. 0) print *, iday, "Done writing output to ", trim(filename_ncOut), "..."

END SUBROUTINE WritePoolFluxNcFile_mimics_daily

!-------------------------------------------------------------------------------------



