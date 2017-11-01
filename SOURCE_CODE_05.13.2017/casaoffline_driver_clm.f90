!--------------------------------------------------------------------------------
! FILE: casa_driver_clm.f90
!
! This program has been modified in order to drive the CASACNP model 
! with surface and meteorological data derived from CLM history files.
! Melannie D. Hartman, March 2014
!
! Added the MIMICS model as optional SOM model (isomModel=2)
! Melannie D. Hartman, January-November 2015
!
! Added the CORPSE model as another optional SOM model (isomModel=3)
! Melannie D. Hartman, December 2015
!
!--------------------------------------------------------------------------------

PROGRAM offline_casacnp
!  calling sequence
!  1: allocate variables used in casacnpls
!  2: read in tile information
!  3: read in biome information
!  4: initilization
!  5: casacnpdriver->biogeochem->all C, N and P cycle routines
!  6: nloop=mloop, casa_poolout to write pool size out
!
  USE define_dimensions
  USE define_types
  USE casadimension
  USE casaparm
  USE casavariable
  USE mimicsvariable
  USE phenvariable
  USE corpsedimension
  USE corpsevariable

  IMPLICIT NONE
  integer mst, mvt, mloop, nloop, mdaily
  integer tyear1, tyear2, ityear, idx, idx2
  integer mreps, irep, iYrCnt, nctime
  character(len=100) :: filename_cnppoint,filename_phen, &
                        filename_cnpipool,filename_cnpmet, &
                        filename_cnpepool,filename_cnpflux, &
                        filename_soilprop, filename_ncOut
  !!character(len=100) :: filename_cnpbiome !! file for biome-specific CASACCNP parameters (see casavariable)
  character(len=100) :: filename_lst

  ! MIMICS UPDATE (-MDH 1/5/2015)
  !!character(len=100) :: filename_mimicsbiome  !! file for biome-specific MIMICS parameters (see mimicsvariable)
  character(len=100) :: filename_mimicsipool    !! file for initial MIMICS pool amounts
  character(len=100) :: filename_mimicsepool    !! file for end-of-run MIMICS pool amounts
  character(len=100) :: filename_ncOut_mimics   !! NetCDF file for annual average pools and fluxes

  ! CORPSE UPDATE (-MDH 12/28/2015)
  character(len=100) :: filename_corpseipool    !! file for initial CORPSE pool amounts
  character(len=100) :: filename_corpseepool    !! file for end-of-run CORPSE pool amounts
  !!character(len=100) :: filename_corpsenamelist !! file CORPSE namelist parameter file (see corpsevariable)
  character(len=100) :: filename_ncOut_corpse   !! NetCDF file for annual average pools and fluxes
  integer:: maxSteps                            !! Maximum number of CORPSE timesteps during the simulation

  character(len=100) :: filename_co2delta       !! File specifying deltatair, deltatsoil, NPPscalar, year (formerly co2delta.txt)

  character(len=100) :: dirPtFile               !! Directory where point file is to be written. If empty assum "./".

  character(len=4) syear
  logical back
  real(r_2)    co2air		! Atmospheric concentration (ppm) (Not used currently -mdh 2/22/2016) 
  real(r_2)    deltsoil		! Increment in soil temperature from met file (deg. C)
  real(r_2)    deltair		! Increment in air temperature from met file (deg. C)
  integer deltYr		! Year count when deltsoil or deltair begins (-mdh 2/22/2016).
  real(r_2)    nppMult		! Multiplier on NPP to increase or decrease litter inputs

  mst = 9  ! casaclm does not use soil types, but the dimension is assigned 
           ! to minimze code changes to the rest of the model.
           
  open(57,file='errors.txt')

  filename_lst = 'fcasacnp_clm_testbed.lst'
  open(10,file= filename_lst)

  read(10,*) mp      		!! number of points to simulate. Must = # of non-zero cellMissing values met.nc file.
  print *, 'mp=', mp            !!   Consistency is checked when ReadMetNcFile is called.
  read(10,*) mvt		!! number of vegetation types  !=13 for CASA; =18 for IGBP+tundra
  print *, 'mvt=', mvt
  read(10,*) mloop		!! number of times to cycle through met file 
  print *, 'mloop=', mloop
  read(10,*) mdaily 		!! =1 output daily fluxes, =0 output yearly fluxes
  print *, 'mdaily=', mdaily

  read(10,*) initcasa 		!! =0 spin; =1 restart file; =2 transient; 3=repeated transient
  print *, 'initcasa=', initcasa
  if (initcasa == 2) then
     read(10, *) tyear1, tyear2
     write(*,*) 'transient year start =', tyear1
     write(*,*) 'transient year end =', tyear2
     if (mloop .ne. 1) then
         write(*,*) 'Resetting mloop to 1 for transient simulation'
     endif
  else if (initcasa == 3) then
     ! Allow transient years to be repeated.  This is an alternative way to do a spinup run. -mdh 5/13/2017
     read(10, *) tyear1, tyear2
     write(*,*) 'This sequence of years will be repeated ', mloop, ' times:'
     write(*,*) '  year start =', tyear1
     write(*,*) '  year end =', tyear2
  else
     !Skip this line
     read(10, *)
  endif

  read(10,*) isomModel 		!! SOM cycling model: 1=CASACNP, 2=MIMICS, 3=CORPSE
  print *, 'isomModel=', isomModel
  if (isomModel < CASACNP .or. isomModel > CORPSE) then
      write(*,*) 'Unexpected value for isomModel =', isomModel
      write(*,*) 'Expecting 1, 2, or 3 in file ', trim(filename_lst)
      STOP
  endif

  read(10,*) icycle 		!! =1 for C only, =2 for C+N; =3 for C+N+P
  print *, 'icycle=', icycle

  if ((isomModel >= MIMICS) .AND. (icycle /= 1)) then
     print *, 'For MIMICS (isomModel=2) and CORPSE (isomModel=3), icycle must = 1'
     print *, 'Update options in file ', trim(filename_lst)
     STOP
  endif

  read(10,101) filename_cnppoint	!! gridinfo_igbpz.csv (unit=101)
  call rmcomments(filename_cnppoint)
  read(10,101) filename_cnpbiome	!! pftlookup_igbp.csv (unit=101)
  call rmcomments(filename_cnpbiome)
  read(10,101) filename_phen		!! modis_phenology.txt (unit=101)
  call rmcomments(filename_phen)
  read(10,101) filename_soilprop	!! gridinfo_soil.csv (unit=101)
  call rmcomments(filename_soilprop)
  read(10,101) filename_cnpmet		!! NetCDF Met file (formerly met.txt (unit=111))
  call rmcomments(filename_cnpmet)
  read(10,101) filename_cnpipool	!! cnppool_igbp_start_q.txt (unit=99)
  call rmcomments(filename_cnpipool)
  read(10,101) filename_cnpepool	!! cnppool_igbp_end_q.txt (unit=103)
  call rmcomments(filename_cnpepool)
  read(10,101) filename_cnpflux		!! cnpflux_igbp_end_q.txt (unit=104)
  call rmcomments(filename_cnpflux)
  read(10,101) filename_ncOut		!! CASACNP NetCDF Output file
  call rmcomments(filename_ncOut)

  !! MIMICS Files
  read(10,101) filename_mimicsbiome     !! pftlookup_mimics.csv (unit=101)
  call rmcomments(filename_mimicsbiome)
  read(10,101) filename_mimicsipool	!! mimics_cpool_start.txt (unit=105)
  call rmcomments(filename_mimicsipool)
  read(10,101) filename_mimicsepool	!! mimics_cpool_end.txt (unit=106)
  call rmcomments(filename_mimicsepool)
  read(10,101) filename_ncOut_mimics	!! MIMICS NetCDF Output file
  call rmcomments(filename_ncOut_mimics)

  !! CORPSE Files
  read(10,101) filename_corpseipool     !! corpse_cpool_start.txt (unit=105)
  call rmcomments(filename_corpseipool)
  read(10,101) filename_corpseepool     !! corpse_cpool_end.txt (unit=106)
  call rmcomments(filename_corpseepool)
  read(10,101) filename_corpsenamelist	!! CORPSE namelist parameters (namelistunit = 123)
  call rmcomments(filename_corpsenamelist)
  read(10,101) filename_ncOut_corpse	!! CORPSE NetCDF Output file
  call rmcomments(filename_ncOut_corpse)

  !! Read the delta file that controls changes in air/soil temperature and NPP during the simulation.
  read(10,101) filename_co2delta
  call rmcomments(filename_co2delta)
  open(1,file=filename_co2delta)
  read(1,*)
  read(1,*) co2air,deltsoil,deltair,deltYr,nppMult
  print *, 'co2air deltsoil deltair deltYr nppMult icycle =',co2air,deltsoil,deltair,deltYr,nppMult,icycle

  read(10,*) casafile%iptToSave         !! If > 0, save daily output for this point in .csv file
  if (casafile%iptToSave > 0) then
      read(10,101) dirPtFile            !! Directory where point file will be written to
      call rmcomments(dirPtFile)
      write(*,*) 'Directory for Point Files = ', trim(dirPtFile)
  else
      dirPtFile = ""
  endif

  call alloc_casavariable(mp,mvt,mst,ms)
  call alloc_mimicsvariable(mp,mvt,mplant)  !! MIMICS UPDATE (-MDH 1/5/2015)
  call alloc_phenvariable(mp,mvt)
  call alloc_casavegsoil(mp,ms)
  call alloc_corpsevariable(mp)             !! CORPSE UPDATE (-MDH 12/28/2015)

  print *,'calling casa_readpoint'
  call casa_readpoint(filename_cnppoint,mvt)

  !Added filename_soilprop to casa_readbiome function list (-MDH 2/17/2014)
  print *, 'calling casa_readbiome'
  call casa_readbiome(filename_cnpbiome,filename_soilprop,mvt,mst)

  if (isomModel == MIMICS) then

      print *, 'calling mimics_readbiome'
      call mimics_readbiome(filename_mimicsbiome,mp,mvt)

  else if (isomModel == CORPSE) then

      !Read parameters from namelist
      print *, 'Reading CORPSE parameters from namelist...', filename_corpsenamelist
      OPEN(unit=namelistunit,file=filename_corpsenamelist)
      READ(unit=namelistunit,NML=CORPSE_casa_nml)
      CLOSE(unit=namelistunit)
      ! This subroutine was never actually called before, so namelist parameters had default values -mdh 3/20/2017
      call read_soil_carbon_namelist(filename_corpsenamelist)
      if (mdaily == 1) then
         !!recordtime = 24	! output interval is every 24 hours
         recordtime = 1	        ! output interval is once a day (daily timestep - mdh 3/21/2016)
      else
         !!recordtime = 8760	! output interval is at the end of every year
         recordtime = 365	! output interval is at the end of every year (daily timestep - mdh 3/21/2016)
      endif

  endif

  print *, 'calling casa_readphen'
  call casa_readphen(filename_phen,mvt)

  print *, 'calling casa_init'
  call casa_init(filename_cnpipool,mp,ms,mst)

  if (casafile%iptToSave > 0) then
      !casa_init must be called to initialize casamet before calling this subroutine
      call WritePointFileHeaders(dirPtFile,mp)
  endif 

  if (isomModel == MIMICS) then

      print *, 'calling mimics_init'
      call mimics_init(filename_mimicsipool,mp,ms,mst)

  else if (isomModel == CORPSE) then

      print *, 'calling corpse_init'
      !! ATTENTION: The allocation of output variables may need to be moved after met.nc file is read to get the exact # simulation years.
      if (initcasa < 2) then
          maxSteps = mloop * 365 * 7                          ! 7 is the maximum number of years in a spinup file (daily time step - mdh 3/21/2016)
      else
          maxSteps = mloop * 365 * abs(tyear2 - tyear1 + 1)   ! Assumes each transient year met.nc file has 365 days  (daily time step - mdh 3/21/2016)
      endif
      print *, 'maxSteps =', maxSteps
      call corpse_init(maxSteps, filename_corpseipool, mp, initcasa, rhizosphere_frac )

  endif

  if (initcasa < 2) then

     !Non-transient run.  One met.nc file is used for the entire simulation.

     print *, 'calling casacnpdriver'
     call casacnpdriver(filename_cnpmet,filename_cnpepool,filename_cnpflux, filename_ncOut, &
                        filename_mimicsepool, filename_ncOut_mimics, &
                        filename_corpseepool, filename_ncOut_corpse, &
                        mloop, mdaily, co2air, deltsoil, deltair, deltYr, nppMult, 0)
  else

     !Transient run. Read a different met.nc file each year, output to different NetCDF files each year.

     back = .true.
     idx = index(filename_cnpmet, '.nc', back)	! start position of last occurrence of '.nc' (count starts at 1)
     if (idx == 0) then
        write(*,*) 'Error in transient met.nc file name: ', trim(filename_cnpmet) 
        write (*,*)'The name does not end with .nc extension.  See file ', trim(filename_lst)
        STOP
     endif

     if (tyear1 > 9999 .or. tyear1 < 1 .or. tyear2 > 9999 .or. tyear2 < 1) then
         write(*,*) 'Error specifying transient years.  The range is 1-9999.' 
         write (*,*)'See file ', trim(filename_lst)
         STOP
     endif

     ! Allow transient years to be repeated.  This is an alternative way to do a spinup run. -mdh 5/13/2017
     if (initcasa == 3) then
         mreps = mloop
     else
        mreps = 1
     endif
     mloop = 1

     do irep = 1, mreps
         do ityear = tyear1, tyear2

            iYrCnt = (irep-1)*(tyear2-tyear1+1) + (ityear-tyear1+1)

            ! Determine the name of the met.nc file to read by replacing the last 9 characters before ".nc" with 
            ! yyyy_yyyy, where yyyy is the value of ityear.

            write(syear,'(i4)') ityear	! copy the value of ityear to string syear
            if (ityear .lt. 1000) then      ! make sure syear is 4 digits, add leading zeros if necessary
                syear(1:1) = '0'
            endif
            if (ityear .lt. 100) then
                syear(2:2) = '0'
            endif
            if (ityear .lt. 10) then
                syear(3:3) = '0'
            endif
            filename_cnpmet(idx-9:idx-6) = syear
            filename_cnpmet(idx-4:idx-1) = syear
            
            !For transient runs, different output .nc files will be generated each year.
            !Replace the last 4 characters before .nc, usually yyyy, with the value of ityear the NetCDF file name.
            !Do this for filename_ncOut, filename_ncOut_mimics, filename_ncOut_corpse
    
            if (mreps == 1) then
                !This is a normal transient run, so insert the calendar year in the output file names
                call InsertYearInNcFileName(filename_ncOut, ityear)
                call InsertYearInNcFileName(filename_ncOut_mimics, ityear)
                call InsertYearInNcFileName(filename_ncOut_corpse, ityear)
                nctime = ityear
            else
                !This is a repeated transient run, so insert the year count in the output file names
                call InsertYearInNcFileName(filename_ncOut, iYrCnt)
                call InsertYearInNcFileName(filename_ncOut_mimics, iYrCnt)
                call InsertYearInNcFileName(filename_ncOut_corpse, iYrCnt)
                nctime = iYrCnt
            endif

    !!      Since N dep is read from met.nc, no longer need annual filename_cnppoint file. -mdh 8/22/2016
    !!      !For transient runs, read N deposition from filename_cnppoint each year!
    !!      !Replace the last 4 characters before .csv, usually yyyy, with the value of ityear in the output .csv file
    !!      idx2 = index(filename_cnppoint, '.csv', back)    ! idx2 = start position of last occurrence of '.csv' (count starts at 1)        
    !!      if (idx2 == 0) then
    !!          write(*,*) 'Error in transient file name: ', trim(filename_cnppoint) 
    !!          write (*,*)'The name does not end with .csv extension.  See file ', trim(filename_lst)
    !!          STOP
    !!      endif
    !!      filename_cnppoint(idx2-4:idx2-1) = syear

            print *, 'calling casa_readpoint for year ', syear
            call casa_readpoint(filename_cnppoint,mvt)
    
            print *, 'calling casacnpdriver for transient run'
            print *, '  filename_cnpmet =', filename_cnpmet
            print *, '  filename_ncOut =', filename_ncOut
            if (isomModel .eq. MIMICS) then
                print *, '  filename_ncOut_mimics = ', filename_ncOut_mimics
            endif
            if (isomModel .eq. CORPSE) then
                print *, '  filename_ncOut_corpse = ', filename_ncOut_corpse
            endif
            print *, '  filename_cnppoint = ', filename_cnppoint

            call casacnpdriver(filename_cnpmet,filename_cnpepool,filename_cnpflux,filename_ncOut, &
                               filename_mimicsepool, filename_ncOut_mimics, &
                               filename_corpseepool, filename_ncOut_corpse, &
                               mloop, mdaily, co2air, deltsoil, deltair, deltYr, nppMult, nctime)
         enddo
     enddo
  endif

  close(10)

 
  print *, 'simulation completed'
101 format(a100)

END PROGRAM offline_casacnp
