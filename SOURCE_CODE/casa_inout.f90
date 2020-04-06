!--------------------------------------------------------------------------------
! FILE: casa_inout.f90
!
! Subroutines in this file:
!   casa_readbiome() - reads pftlookup.csv file and gridinfo_soil.csv
!   casa_readphen()  - reads tabulated modis-derived leaf phenology data file
!   casa_readpoint() - reads gridinfo_ibgpz.csv 
!   casa_init(filename_cnpipool,mp,ms,mst) - read filename_cnpipool to initialize CASACNP pools
!   casa_poolout()   - writes end of simulation pool values for restart file
!   casa_fluxout()   - writes simulation average annual flux values 
!   casa_cnpflux()   - calculates annual fluxes by summing daily fluxes
!   casa_cnppool()   - calculates average annual C, N, P litter, soil, plant pools
!   biogeochem()     - called once a day. Calls all plant, litter, soil subroutines
!                      for all points in the grid
!   casacnpdriver()  - main driver for CASACNP, calls other subroutines
!
! The following functions/subroutines have been added in order to run the
! CASACNP model with surface and meteorological data derived from CLM
! history files.  I/O routines for MIMICS and CORPSE output are in other source files.
!
!   handle_err        - Handle NetCDF read/write errors
!   get_time_and_date - Compute date and time string for global NetCDF attributes
!   GetMetNcFileDim   - Retrieve the variable dimensions from the met.nc file.
!   ReadMetNcFile     - Read in daily values from the met.nc file
!   PutVariableAttributeReal           - Write variable attributes to NetCDF file
!   PutVariableAttributeDouble         - Write variable attributes to NetCDF file
!   WritePoolFluxNcFile_casacnp_annual - Write annual mean CASACNP pool and flux output 
!                                        to NetCDF file
!   WritePoolFluxNcFile_casacnp_daily  - Write daily CASACNP pool and flux output to
!                                        NetCDF file (IN PROGRESS 6/6/2016)
!   InsertYearInNcFileName    - Insert the simulation year into the .nc filename
!   CreateDailyNcFileName     - Create a daily .nc filename from the annual filesname 
!                               and insert simulation year into the filename
!   rmcomments - replace any non-space characers to the right of a ' ' with a ' '.
!                Used to strip comments from a line in fcasacnp_clm.lst, elaving only the file name
!   The following two subroutines write to a daily point file (.csv format) for testing:
!     * WritePointFileHeaders - Create headers for .csv files that will contain daily ouput
!                               for the test point. Includes headers for all SOM models (isomModel = 1,2,3).
!     * WritePointCASA        - Write to CASACNP daily point file (isomModel = 1)

!
! The assignment of soil properties in casa_readbiome was changed. The values
! soil%sand(*),soil%clay(*),soil%silt(*),soil%swilt(*),soil%sfc(*),soil%ssat(*)
! for each grid cell are read directly from a file rather than assigned
! from soil type (1-7).
!
! Melannie D. Hartman, March 2014 - December 2019
! Removed NPP and LAI from met.nc files. -mdh 11/6/2017
! Echo CASA parameter values as they are read in (subroutine casa_readbiome). -mdh 11/6/2017
! Echo restart file contents (subroutine casa_init). -mdh 11/6/2017
! Add N fluxes to netCDF output. -mdh 11/25/2019
! Added new set of parameters to the botton of the CASA PFT parameter file. -mdh 12/23/2019
!--------------------------------------------------------------------------------

SUBROUTINE casa_readbiome(fname_cnpbiome,filename_soilprop,mvt,mst)
  use define_dimensions
  use define_types
  use casadimension
  use casaparm
  use casavariable
  use phenvariable
  implicit none
  character(len=100)  fname_cnpbiome
  character(len=100)  filename_soilprop
  character(len=300)  buffer
  integer mvt,mst
  integer  i,iv1,nv,ns,nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12,npt,iv,is,iso
  real(r_2),    dimension(mvt)          :: leafage,frootage,woodage
  real(r_2),    dimension(mvt)          :: totroot
  real(r_2),    dimension(mvt)          :: cwdage,metage,strage
  real(r_2),    dimension(mvt)          :: micage,slowage,passage,clabileage,slax
  real(r_2),    dimension(mvt,mplant)   :: ratioCNplant
  !! Added (-MDH 6/9/2014)
  !!real(r_2),  dimension(mvt,msoil)    :: ratioCNsoil
  real(r_2),    dimension(mvt,msoil)    :: ratioCNsoil,ratioCNsoilmin,ratioCNsoilmax
  real(r_2),    dimension(ms)           :: depthsoila,depthsoilb
  real(r_2),    dimension(mvt)          :: xfNminloss, xfNminleach, xnfixrate
  real(r_2),    dimension(mvt)          :: cleaf,cwood,cfroot, &
                                      cmet,cstr,ccwd, &
                                      cmic,cslow,cpass
  real(r_2),    dimension(mvt)          :: nleaf,nwood,nfroot, &
                                      nmet,nstr,ncwd, &
                                      nmic,nslow,npass,xnsoilmin
  real(r_2),    dimension(mvt)          :: xpleaf, xpwood, xpfroot, &
                                      xpmet, xpstr, xpcwd, &
                                      xpmic,xpslow,xppass,xplab,xpsorb,xpocc

  real(r_2),    dimension(mso)          :: xkmlabp,xpsorbmax,xfPleach
  real(r_2),    dimension(mso,msoil)    :: ratioNPsoil
  real(r_2)     xratioNPleafmin,xratioNPleafmax,                 &
                xratioNPwoodmin,xratioNPwoodmax,                 &
                xratioNPfrootmin,xratioNPfrootmax 

  real(r_2),    dimension(mvt)          :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase, &
                                           xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope
  real(r_2),    dimension(mso)          :: xxkplab,xxkpsorb,xxkpocc

  real(r_2),    dimension(mvt)          :: xfherbivore,xxkleafcoldmax, xxkleafdrymax
  real(r_2),    dimension(mvt)          :: xkuplabp
  real(r_2),    dimension(mvt,ms)       :: fracroot 
  real(r_2) :: lat, lon
  integer   :: ipt, IOstatus


!
! ==== the following to be commented out when coupled with CABLE 
  real(r_2),    dimension(9)          :: silt,clay,sand,wwilt,wfield,wsat,dsoil
  real(r_2),    dimension(6)           :: dzsoil
    
   data silt/  .08,    .33, .17, .2, .06, .25, .15, .70, .33/
   data clay/  .09,    .3, .67, .2, .42, .48, .27, .17, .30/
   data sand/  .83,    .37, .16, .6, .52, .27, .58, .13, .37/
   data wwilt/ .072,   .216, .286, .135, .219, .283, .175, .395, .216/
   data wfield/.143,   .301, .367, .218, .31 , .37 , .255, .45,  .301/
   data wsat/ .398,    .479, .482, .443, .426, .482, .420, .451, .47/
   data dsoil/7*1.6e6, 1.3e6,  0.91e6/                                 ! soil density (g/m3)
   data dzsoil/.022, .058, .154, .409, 1.085, 2.872/                   ! layer thickness (m)
!  assign some data to vegetation and soil
   soil%dzsoil(:)= dzsoil(:)

!------------------------------------------------------------------------------------------------------------------------------------------
!  Replaced the assignments below with values read from gridinfo_soil.csv file (-MDH 2/17/2014)
!  soil%silt(:)  = silt(soil%isoilm(:))
!  soil%clay(:)  = clay(soil%isoilm(:))
!  soil%sand(:)  = sand(soil%isoilm(:))
!  soil%sfc(:)   = wfield(soil%isoilm(:))
!  soil%swilt(:) = wwilt(soil%isoilm(:))
!  soil%ssat(:)  = wsat(soil%isoilm(:))
   
   write(*,*)
   write(*,*) "Reading soil properties from file ", filename_soilprop, "..."
   open(121,file=filename_soilprop)

   read(121,'(a)') buffer    ! Read past column header
!  write(*,*) trim(buffer)

   do i=1,mp
      read(121,*,IOSTAT=IOstatus) ipt,lat,lon,soil%sand(i),soil%clay(i),soil%silt(i),soil%swilt(i),soil%sfc(i),soil%ssat(i)
      if (IOstatus .lt. 0) then
         write(*,*)
         write(*,*) 'Unexpected EOF while reading file ', TRIM(filename_soilprop)
         write(*,*) 'Expected ', mp, ' gridcells.  Found ', i-1, 'gridcells.'
         STOP
      elseif (IOstatus .gt. 0) then
         print *, "An error occurred while reading file ", TRIM(filename_soilprop)
         STOP
      endif
      write(*,122) ipt,lat,lon,soil%sand(i),soil%clay(i),soil%silt(i),soil%swilt(i),soil%sfc(i),soil%ssat(i)
      122 format(i6,',',2(f10.4,','),6(f8.4,','))
   enddo
   write(*,*) "Done reading soil properties from file ", filename_soilprop, "..."

!------------------------------------------------------------------------------------------------------------------------------------------

! calculate veg%froot(mp,ms) here

!=========================================
      write(*,*)
      write(*,*) "Reading CASA biome-specific parameters from file ", fname_cnpbiome, "..."
      open(101,file=fname_cnpbiome)
      do i=1,3
         read(101,'(a)') buffer 
         write(*,*) trim(buffer)
      enddo
  
      do nv=1,mvt
         read(101,*) nv0,casabiome%ivt2(nv)
         write(*,*) nv0,casabiome%ivt2(nv)
         ! print *, nv,nv0,casabiome%ivt2(nv)
      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
! Added clabileage(nv),slax(nv) to end of list (-MDH 6/29/2014)
      do nv=1,mvt
         read(101,*) nv1,casabiome%kroot(nv),casabiome%rootdepth(nv),      &
                     casabiome%kuptake(nv),casabiome%krootlen(nv),         &
                     casabiome%kminN(nv), casabiome%kuplabP(nv),           &
                     xfherbivore(nv),leafage(nv),woodage(nv),frootage(nv), &
                     metage(nv),strage(nv),cwdage(nv),  &
                     micage(nv),slowage(nv),passage(nv), & 
                     clabileage(nv),slax(nv)
         write(*,'(i2,1x,18(f12.4,1x))') nv1,casabiome%kroot(nv),casabiome%rootdepth(nv),      &
                     casabiome%kuptake(nv),casabiome%krootlen(nv),         &
                     casabiome%kminN(nv), casabiome%kuplabP(nv),           &
                     xfherbivore(nv),leafage(nv),woodage(nv),frootage(nv), &
                     metage(nv),strage(nv),cwdage(nv),  &
                     micage(nv),slowage(nv),passage(nv), & 
                     clabileage(nv),slax(nv)
         !print *, 'nv1',nv,nv1
      enddo  

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      do nv=1,mvt
         read(101,*) nv2,casabiome%fracnpptoP(nv,leaf),casabiome%fracnpptoP(nv,wood), &
                         casabiome%fracnpptoP(nv,froot),casabiome%rmplant(nv,leaf),   &
                         casabiome%rmplant(nv,wood),casabiome%rmplant(nv,froot)
         write(*,'(i2,1x,6(f12.4,1x))') nv2,casabiome%fracnpptoP(nv,leaf),casabiome%fracnpptoP(nv,wood), &
                         casabiome%fracnpptoP(nv,froot),casabiome%rmplant(nv,leaf),   &
                         casabiome%rmplant(nv,wood),casabiome%rmplant(nv,froot)
         !print *, 'nv2', nv2
      enddo 

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
!     Added ratioCNsoilmin and ratioCNsoilmax for mic, slow, and pass (-MDH 6/29/2014)
!     Added glaimax, glaimin (-MDH 6/29/2014)
      do nv=1,mvt
         read(101,*) nv3, ratioCNplant(nv,leaf),ratioCNplant(nv,wood),   &
             ratioCNplant(nv,froot),                                         &
             casabiome%ftransNPtoL(nv,leaf), casabiome%ftransNPtoL(nv,wood), &
             casabiome%ftransNPtoL(nv,froot),                                & 
             casabiome%fracligninplant(nv,leaf),                             &
             casabiome%fracligninplant(nv,wood),                             &
             casabiome%fracligninplant(nv,froot),                            &
             ratioCNsoil(nv,mic),ratioCNsoil(nv,slow),ratioCNsoil(nv,pass),  &
             ratioCNsoilmin(nv,mic),ratioCNsoilmin(nv,slow),ratioCNsoilmin(nv,pass),  &
             ratioCNsoilmax(nv,mic),ratioCNsoilmax(nv,slow),ratioCNsoilmax(nv,pass),  &
   !         xfherbivore(nv),casabiome%ratiofrootleaf(nv),                  &
             casabiome%glaimax(nv),casabiome%glaimin(nv)

         write(*,'(i2,1x,20(f12.4,1x))') nv3, ratioCNplant(nv,leaf),ratioCNplant(nv,wood),   &
             ratioCNplant(nv,froot),                                         &
             casabiome%ftransNPtoL(nv,leaf), casabiome%ftransNPtoL(nv,wood), &
             casabiome%ftransNPtoL(nv,froot),                                & 
             casabiome%fracligninplant(nv,leaf),                             &
             casabiome%fracligninplant(nv,wood),                             &
             casabiome%fracligninplant(nv,froot),                            &
             ratioCNsoil(nv,mic),ratioCNsoil(nv,slow),ratioCNsoil(nv,pass),  &
             ratioCNsoilmin(nv,mic),ratioCNsoilmin(nv,slow),ratioCNsoilmin(nv,pass),  &
             ratioCNsoilmax(nv,mic),ratioCNsoilmax(nv,slow),ratioCNsoilmax(nv,pass),  &
   !         xfherbivore(nv),casabiome%ratiofrootleaf(nv),                  &
             casabiome%glaimax(nv),casabiome%glaimin(nv)

         !print *, 'nv3',nv3

      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      do nv=1,mvt
         read(101,*) nv4,                                              &
             cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),cstr(nv),ccwd(nv), &
             cmic(nv), cslow(nv),cpass(nv)
         write(*,'(i2,1x,9(f12.4,1x))') nv4,                                              &
             cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),cstr(nv),ccwd(nv), &
             cmic(nv), cslow(nv),cpass(nv)
         !print *, 'nv4',nv4

      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      do nv=1,mvt
         read(101,*) nv5,phen%TKshed(nv),xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv),   &
             xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
         write(*,'(i2,1x,5(f12.4,1x))') nv5,phen%TKshed(nv),xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv),   &
             xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
         !print *, 'nv5',nv5
      enddo

!      read(101,*)
!      read(101,*)
!      do nv=1,mvt
!         read(101,*) nv5, &
!         xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv),   &
!         xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
!
!      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      do nv=1,mvt
         read(101,*) nv6, &
             casabiome%ratioNCplantmin(nv,leaf),casabiome%ratioNCplantmax(nv,leaf), &
             casabiome%ratioNCplantmin(nv,wood),casabiome%ratioNCplantmax(nv,wood), &
             casabiome%ratioNCplantmin(nv,froot),casabiome%ratioNCplantmax(nv,froot), &
             xfNminloss(nv), xfNminleach(nv),xnfixrate(nv)
         write(*,'(i2,1x,9(f12.4,1x))') nv6, &
             casabiome%ratioNCplantmin(nv,leaf),casabiome%ratioNCplantmax(nv,leaf), &
             casabiome%ratioNCplantmin(nv,wood),casabiome%ratioNCplantmax(nv,wood), &
             casabiome%ratioNCplantmin(nv,froot),casabiome%ratioNCplantmax(nv,froot), &
             xfNminloss(nv), xfNminleach(nv),xnfixrate(nv)
      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      do nv=1,mvt
         read(101,*) nv7,nleaf(nv),nwood(nv),nfroot(nv), &
                     nmet(nv),nstr(nv), ncwd(nv), &
                     nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
         write(*,'(i2,1x,10(f12.4,1x))') nv7,nleaf(nv),nwood(nv),nfroot(nv), &
                     nmet(nv),nstr(nv), ncwd(nv), &
                     nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      do nv=1,mvt
         read(101,*) nv8,xratioNPleafmin,xratioNPleafmax,      &
             xratioNPwoodmin,xratioNPwoodmax,                      &
             xratioNPfrootmin,xratioNPfrootmax,                    &
             casabiome%ftransPPtoL(nv,leaf), casabiome%ftransPPtoL(nv,wood), &
             casabiome%ftransPPtoL(nv,froot)
         write(*,'(i2,1x,9(f12.4,1x))') nv8,xratioNPleafmin,xratioNPleafmax,      &
             xratioNPwoodmin,xratioNPwoodmax,                      &
             xratioNPfrootmin,xratioNPfrootmax,                    &
             casabiome%ftransPPtoL(nv,leaf), casabiome%ftransPPtoL(nv,wood), &
             casabiome%ftransPPtoL(nv,froot)

!!       !! ratioPcpplantmin/ratioPcpplantmin are not defined for casabiome in this version (-MDH 6/9/2014)
!!       casabiome%ratioPcplantmin(nv,leaf)  = 1.0/(xratioNPleafmin*ratioCNplant(nv,leaf))
!!       casabiome%ratioPcplantmax(nv,leaf)  = 1.0/(xratioNPleafmax*ratioCNplant(nv,leaf))
!!       casabiome%ratioPcplantmin(nv,wood)  = 1.0/(xratioNPwoodmin*ratioCNplant(nv,wood))
!!       casabiome%ratioPcplantmax(nv,wood)  = 1.0/(xratioNPwoodmax*ratioCNplant(nv,wood))
!!       casabiome%ratioPcplantmin(nv,froot) = 1.0/(xratioNPfrootmin*ratioCNplant(nv,froot))
!!       casabiome%ratioPcplantmax(nv,froot) = 1.0/(xratioNPfrootmax*ratioCNplant(nv,froot))
         !print *, 'nv8',nv8

         !! Added (-MDH 6/9/2014)
         casabiome%ratioNPplantmin(nv,leaf)  = xratioNPleafmin
         casabiome%ratioNPplantmax(nv,leaf)  = xratioNPleafmax
         casabiome%ratioNPplantmin(nv,wood)  = xratioNPwoodmin
         casabiome%ratioNPplantmax(nv,wood)  = xratioNPwoodmax
         casabiome%ratioNPplantmin(nv,froot) = xratioNPfrootmin
         casabiome%ratioNPplantmax(nv,froot) = xratioNPfrootmax

      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
! Added xxkplab(iso),xxkpsorb(iso),xxkpocc(iso) to end of list (-MDH 6/29/2014)
      do iso=1,mso
         read(101,*) nv9,xkmlabp(iso),xpsorbmax(iso),xfPleach(iso), &
                     ratioNPsoil(iso,mic),ratioNPsoil(iso,slow),ratioNPsoil(iso,pass), &
                     xxkplab(iso),xxkpsorb(iso),xxkpocc(iso)
         write(*,'(i2,1x,9(f12.4,1x))')  nv9,xkmlabp(iso),xpsorbmax(iso),xfPleach(iso), &
                     ratioNPsoil(iso,mic),ratioNPsoil(iso,slow),ratioNPsoil(iso,pass), &
                     xxkplab(iso),xxkpsorb(iso),xxkpocc(iso)

         !print *, 'nv9',nv9
      enddo

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      do nv=1,mvt
         read(101,*) nv10, xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
                          xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
         write(*,'(i2,1x,12(f12.4,1x))') nv10, xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
                          xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
         !print *, 'nv10',nv10
      enddo

! Added this new section (-MDH 6/29/2014)
 !@@@@@@@@@@@@@@@@@@@@@@@@@

      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      read(101,'(a)') buffer
      write(*,*) trim(buffer)
      DO nv=1,mvt
        read(101,*) nv11, &
             xxnpmax(nv),xq10soil(nv),xxkoptlitter(nv),xxkoptsoil(nv),xprodptase(nv), &
             xcostnpup(nv),xmaxfinelitter(nv),xmaxcwd(nv),xnintercept(nv),xnslope(nv)                   
        write(*,'(i2,1x,10(f12.4,1x))') nv11, &
             xxnpmax(nv),xq10soil(nv),xxkoptlitter(nv),xxkoptsoil(nv),xprodptase(nv), &
             xcostnpup(nv),xmaxfinelitter(nv),xmaxcwd(nv),xnintercept(nv),xnslope(nv)                   
        !print *, 'nv11',nv11

      ENDDO
!@@@@@@@@@@@@@@@@@@@@@

!!-------------------------------------------------------------------------------------------
!! Added this section (-MDH 12/23/2019)

      ! Set default values
      casabiome%xkNlimitmin(:) = 0.50
      casabiome%xkNlimitmax(:) = 2.0
      casabiome%fracRootExud(:) = 0.0 
!     casabiome%CNRootExud(:) = 15.0 
      casabiome%CUEmetbmic(:) = 0.45 
      casabiome%CUEstrmic(:) = 0.45 
      casabiome%CUEstrslow(:) = 0.70 
      casabiome%CUEcwdmic(:) = 0.40 
      casabiome%CUEcwdslow(:) = 0.70 
      casabiome%CUEmicslow(:) = 1.0 
      casabiome%CUEmicpass(:) = 1.0 
      casabiome%CUEpassslow(:) = 0.45 

      ! Handle unexpected end of file
      read(101,'(a)', end=999) buffer
      write(*,*) trim(buffer)
      read(101,'(a)',end=999) buffer
      write(*,*) trim(buffer)
      DO nv=1,mvt
        read(101,*,end=999) nv12, &
             casabiome%xkNlimitmin(nv),casabiome%xkNlimitmax(nv),casabiome%fracRootExud(nv), &
             casabiome%CUEmetbmic(nv),casabiome%CUEstrmic(nv), &
             casabiome%CUEstrslow(nv),casabiome%CUEcwdmic(nv),casabiome%CUEcwdslow(nv), &
             casabiome%CUEmicslow(nv),casabiome%CUEmicpass(nv),casabiome%CUEpassslow(nv)                  
        write(*,'(i2,1x,12(f7.4,1x))') nv12, &
             casabiome%xkNlimitmin(nv),casabiome%xkNlimitmax(nv),casabiome%fracRootExud(nv), &
             casabiome%CUEmetbmic(nv),casabiome%CUEstrmic(nv), &
             casabiome%CUEstrslow(nv),casabiome%CUEcwdmic(nv),casabiome%CUEcwdslow(nv), &
             casabiome%CUEmicslow(nv),casabiome%CUEmicpass(nv),casabiome%CUEpassslow(nv)                    
        !print *, 'nv12',nv12

      ENDDO
!!-------------------------------------------------------------------------------------------

999   continue
      close(101)
      write(*,*) "Done reading CASA biome-specific parameters from file ", fname_cnpbiome, "..."

      fracroot   = 0.0
      depthsoila = 0.0
      depthsoilb = 0.0
      do ns=1,ms
         depthsoilb(ns) = depthsoilb(ns) + soil%dzsoil(ns)
         if(ns==1) then
            depthsoila(ns) = 0.0
         else
            depthsoila(ns) = depthsoilb(ns-1)
         endif       
      enddo


! Commented out this section (-MDH 6/29/2014)
!!     do nv=1,mvt
!!        casabiome%sla(nv)                = 0.025 * (leafage(nv)**(-0.5))            ! see eqn A1 of Arora and Boer, GCB, 2005
!!!         casabiome%fherbivore(nv)         = deltcasa*xfherbivore(nv)
!!        casabiome%fraclabile(nv,leaf)    = deltcasa*0.6    !1/day   !what is that for????
!!        casabiome%fraclabile(nv,froot)   = deltcasa*0.4    !1/day
!!        casabiome%fraclabile(nv,wood)    = deltcasa*0.0
!!        casabiome%plantrate(nv,leaf)     = deltcasa/(leafage(nv)*(1.0-xfherbivore(nv)))
!!        casabiome%plantrate(nv,froot)    = deltcasa/frootage(nv)
!!        casabiome%plantrate(nv,wood)     = deltcasa/woodage(nv)
!!        casabiome%litterrate(nv,metb)    = deltcasa/metage(nv)
!!        casabiome%litterrate(nv,str)     = deltcasa/strage(nv)
!!        casabiome%litterrate(nv,cwd)     = deltcasa/cwdage(nv)
!!        casabiome%soilrate(nv,mic)       = deltcasa/micage(nv)
!!        casabiome%soilrate(nv,slow)      = deltcasa/slowage(nv)
!!        casabiome%soilrate(nv,pass)      = deltcasa/passage(nv)
!!        casabiome%xkleafcoldmax(nv)      = deltcasa * xxkleafcoldmax(nv)
!!        casabiome%xkleafdrymax(nv)       = deltcasa * xxkleafdrymax(nv)
!!!        casabiome%kuplabp(nv)            = xkuplabp(nv)
!!        casabiome%rmplant(nv,:)          = casabiome%rmplant(nv,:)*deltcasa 
!!
!!        totroot(nv) = (1.0-exp(-casabiome%kroot(nv)*casabiome%rootdepth(nv)))
!!        do ns=1,ms
!!           fracroot(nv,ns) = (exp(-casabiome%kroot(nv)*min(casabiome%rootdepth(nv),depthsoila(ns)))  & 
!!                             -exp(-casabiome%kroot(nv)*min(casabiome%rootdepth(nv),depthsoilb(ns)))) &
!!                             /totroot(nv)
!!            write(*,992)  nv,ns,casabiome%kroot(nv),casabiome%rootdepth(nv), &
!!                          depthsoila(ns),depthsoilb(ns),totroot(nv),fracroot(nv,ns)
!!        enddo             
!!
!!     enddo
  
!!-------------------------------------------------------------------------------------------
!! Added this section (-MDH 6/29/2014)

  DO nv=1,mvt
!    casabiome%sla(nv)             = 0.025 * (leafage(nv)**(-0.5)) ! see eqn A1 of Arora and Boer, GCB, 2005
    casabiome%sla(nv)             = slax(nv) 
!    casabiome%sla(nv)             = 2.0E-4 * exp(6.15)/((12*leafage(nv))**0.46) ! see eqn 6 of Sitch, GCB, 2003
!    casabiome%fherbivore(nv)     = deltcasa*xfherbivore(nv)
    casabiome%fraclabile(nv,leaf) = deltcasa*0.6    !1/day
    casabiome%fraclabile(nv,froot)= deltcasa*0.4    !1/day
    casabiome%fraclabile(nv,wood) = deltcasa*0.0
    casabiome%plantrate(nv,leaf)  = deltcasa/(leafage(nv)*(1.0-xfherbivore(nv)))
    casabiome%plantrate(nv,froot) = deltcasa/frootage(nv)
    casabiome%plantrate(nv,wood)  = deltcasa/woodage(nv)
    casabiome%litterrate(nv,metb) = deltcasa/metage(nv)
    casabiome%litterrate(nv,str)  = deltcasa/strage(nv)
    casabiome%litterrate(nv,cwd)  = deltcasa/cwdage(nv)
    casabiome%soilrate(nv,mic)    = deltcasa/micage(nv)
    casabiome%soilrate(nv,slow)   = deltcasa/slowage(nv)
    casabiome%soilrate(nv,pass)   = deltcasa/passage(nv)
    casabiome%xkleafcoldmax(nv)   = deltcasa * xxkleafcoldmax(nv)
    casabiome%xkleafdrymax(nv)    = deltcasa * xxkleafdrymax(nv)
!    casabiome%kuplabp(nv)         = xkuplabp(nv)
    casabiome%rmplant(nv,:)       = casabiome%rmplant(nv,:)*deltcasa 
    casabiome%kclabrate(nv)       = deltcasa/clabileage(nv)

!! Added back this section (-MDH 6/30/2014)
          totroot(nv) = (1.0-exp(-casabiome%kroot(nv)*casabiome%rootdepth(nv)))
          do ns=1,ms
             fracroot(nv,ns) = (exp(-casabiome%kroot(nv)*min(casabiome%rootdepth(nv),depthsoila(ns)))  & 
                               -exp(-casabiome%kroot(nv)*min(casabiome%rootdepth(nv),depthsoilb(ns)))) &
                               /totroot(nv)
             !write(*,*)  nv,ns,casabiome%kroot(nv),casabiome%rootdepth(nv), &
             !            depthsoila(ns),depthsoilb(ns),totroot(nv),fracroot(nv,ns)

          enddo   

!@@@@@@@@@@@@@@@@@
    casabiome%xnpmax(nv)          = xxnpmax(nv)
    casabiome%q10soil(nv)         = xq10soil(nv)
    casabiome%xkoptlitter(nv)     = xxkoptlitter(nv)
    casabiome%xkoptsoil(nv)       = xxkoptsoil(nv)
    casabiome%prodptase(nv)       = xprodptase(nv)/365.0   ! convert from yearly to daily
    casabiome%costnpup(nv)        = xcostnpup(nv)
    casabiome%maxfinelitter(nv)   = xmaxfinelitter(nv)
    casabiome%maxcwd(nv)          = xmaxcwd(nv)
    casabiome%nintercept(nv)      = xnintercept(nv)
    casabiome%nslope(nv)          = xnslope(nv)
!@@@@@@@@@@@@@@
  ENDDO

!@@@@@@@@@@@@@@
  DO ns=1,mso
    casabiome%xkplab(ns)          =  xxkplab(ns)
    casabiome%xkpsorb(ns)         =  xxkpsorb(ns)
    casabiome%xkpocc(ns)          =  xxkpocc(ns)
  ENDDO
 
!@@@@@@@@@@@@@@

!  PRINT *, 'casabiome%ivt2 = ', casabiome%ivt2

!!-------------------------------------------------------------------------------------------

      do npt=1,mp
         iv1=veg%iveg(npt)
         iso=casamet%isorder(npt)
         veg%froot(npt,:) =fracroot(iv1,:)
!         print *, 'npt,iv1 ', npt, iv1
         casamet%iveg2(npt) =casabiome%ivt2(iv1)
         casamet%lnonwood(npt) = 1
         casapool%cplant(npt,wood)  = 0.0
         casapool%clitter(npt,cwd)  = 0.0
         casapool%nplant(npt,wood)  = 0.0
         casapool%nlitter(npt,cwd)  = 0.0
         casapool%pplant(npt,wood)  = 0.0
         casapool%plitter(npt,cwd)  = 0.0
         if(casamet%iveg2(npt)==forest.or.casamet%iveg2(npt)==shrub) then 
            casamet%lnonwood(npt) = 0
            casapool%cplant(npt,wood)  = cwood(iv1) 
            casapool%clitter(npt,cwd)  = ccwd(iv1)
            casapool%nplant(npt,wood)  = nwood(iv1) 
            casapool%nlitter(npt,cwd)  = ncwd(iv1)
            casapool%pplant(npt,wood)  = xpwood(iv1)
            casapool%plitter(npt,cwd)  = xpcwd(iv1)
         endif
         casapool%cplant(npt,leaf)     = cleaf(iv1)
         casapool%cplant(npt,froot)    = cfroot(iv1)
         casapool%clabile(npt)         = 0.0
         casapool%clitter(npt,metb)     = cmet(iv1)
         casapool%clitter(npt,str)     = cstr(iv1)
         casapool%csoil(npt,mic)       = cmic(iv1)
         casapool%csoil(npt,slow)      = cslow(iv1)
         casapool%csoil(npt,pass)      = cpass(iv1)
         if(icycle==1) then
            casapool%ratioNCplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
         endif

 !! Added this section (-MDH 6/29/2014)
    ! initializing glai in case not reading pool file (eg. during spin)
    casamet%glai(npt) = MAX(casabiome%glaimin(iv1), &
                        casabiome%sla(iv1) * casapool%cplant(npt,leaf))

         casaflux%fNminloss(npt)   = xfNminloss(iv1) 
        ! comment out by ypw 12/07/2009
         casaflux%fNminleach(npt)  = 10.0*xfNminleach(iv1) * deltcasa
        ! casaflux%fNminleach(npt)  = xfNminleach(iv1) 
         casapool%nplant(npt,leaf) = nleaf(iv1)
         casapool%nplant(npt,froot)= nfroot(iv1)
         casapool%nlitter(npt,metb) = nmet(iv1)
!         casapool%nlitter(npt,str) = nstr(iv1)
         casapool%nlitter(npt,str) = cstr(iv1)*ratioNCstrfix
         casapool%nsoil(npt,mic)   = nmic(iv1)
         casapool%nsoil(npt,slow)  = nslow(iv1)
         casapool%nsoil(npt,pass)  = npass(iv1) 
         casapool%nsoilmin(npt)    = xnsoilmin(iv1) 
         casapool%pplant(npt,leaf) = xpleaf(iv1)
         casapool%pplant(npt,froot)= xpfroot(iv1) 
         casapool%plitter(npt,metb) = xpmet(iv1)
!         casapool%plitter(npt,str) = xpstr(iv1)

!!       ! Replace the plitter calculation (-MDH 6/9/2014)
!!       casapool%plitter(npt,str) = cstr(iv1)*ratioPCstrfix
         casapool%plitter(npt,str) = casapool%nlitter(npt,str)/ratioNPstrfix

         casapool%psoil(npt,mic)   = xpmic(iv1)
         casapool%psoil(npt,slow)  = xpslow(iv1)
         casapool%psoil(npt,pass)  = xppass(iv1)
         casapool%psoillab(npt)    = xplab(iv1)
         casapool%psoilsorb(npt)   = xpsorb(iv1)
         casapool%psoilocc(npt)    = xpocc(iv1)
         casaflux%kmlabp(npt)      = xkmlabp(iso)
         casaflux%psorbmax(npt)    = xpsorbmax(iso)
         casaflux%fpleach(npt)     = xfPleach(iso)
         casaflux%Nminfix(npt)     = xnfixrate(iv1)/365.0  
!        write(*,*) 'casa_readbiome: casaflux%Nminfix(',npt,') =', casaflux%Nminfix(npt)

!!       !! Commented out this section (-MDH 6/9/2014)
!!       casapool%rationcplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
!!       casapool%ratiopcplant(npt,:)  = casabiome%ratioPcplantmax(iv1,:)
!!       casapool%rationclitter(npt,:) = casapool%nlitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
!!       casapool%ratiopclitter(npt,:) = casapool%plitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
!!       casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
!!       casapool%ratioPCsoil(npt,:)   = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iso,:))
!!
!!
!!       !casapool%rationcplant(npt,:) = casapool%nplant(npt,:)/(casapool%cplant(npt,:)+1.0e-10)
!!       !casapool%ratiopcplant(npt,:) = casapool%pplant(npt,:)/(casapool%cplant(npt,:)+1.0e-10)
!!       !casapool%ratioNCsoil(npt,:)    = 1.0/ratioCNsoil(iv1,:)
!!       !casapool%ratioPCsoil(npt,:)    = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iv1,:))

         !! ADDED (-MDH 6/9/2014)
         casapool%ratioNCplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
         casapool%ratioNPplant(npt,:)  = casabiome%ratioNPplantmin(iv1,:)
         casapool%ratioNClitter(npt,:) = casapool%nlitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
         casapool%ratioNPlitter(npt,:) = casapool%nlitter(npt,:)/(casapool%plitter(npt,:)+1.0e-10)
         casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
         casapool%ratioNPsoil(npt,:)   = ratioNPsoil(iso,:)
         casapool%ratioNCsoilmin(npt,:)   = 1.0/ratioCNsoilmax(iv1,:)
         casapool%ratioNCsoilmax(npt,:)   = 1.0/ratioCNsoilmin(iv1,:)
         casapool%ratioNCsoilnew(npt,:)   = casapool%ratioNCsoilmax(npt,:)
      enddo

!!    !! Commented out this section (-MDH 6/9/2014)
!!    if(icycle==1) then
!!       casapool%nplant(:,:)       = casapool%cplant(:,:) * casapool%rationcplant(:,:)
!!    else
!!       casapool%Nsoil(:,:)        = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
!!       casapool%Psoil(:,:)        = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
!!       casapool%psoilsorb(:)     =  casaflux%psorbmax(:) * casapool%psoillab(:) &
!!                                   /(casaflux%kmlabp(:)+casapool%psoillab(:))
!!    endif

      !! Added this section (-MDH 6/9/2014)
      if(icycle<2) then
         casapool%Nplant(:,:)  = casapool%Cplant(:,:) * casapool%ratioNCplant(:,:)
         casapool%Nsoil(:,:)   = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
      endif
      if(icycle<3) then
         casapool%Psoil(:,:)   = casapool%Nsoil(:,:)/ casapool%ratioNPsoil(:,:)
         casapool%Psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
                               /(casaflux%kmlabp(:)+casapool%psoillab(:))
      endif

      
!     do npt=1,mp
!        if(veg%iveg(npt)==12) then
!         print *, npt, veg%iveg(npt), casapool%Psoil(npt,:),casapool%psoilsorb(npt), &
!                       casaflux%psorbmax(npt),casapool%psoillab(npt),casaflux%kmlabp(npt)
!        endif
!     enddo

      close(101)
991 format(3(i5,2x),100(f8.2,2x))

END SUBROUTINE casa_readbiome

!--------------------------------------------------------------------------------
SUBROUTINE casa_readphen(filename_phen,mvt)
  ! read in the tabulated modis-derived leaf phenology data for latitude bands of 79.75 to -55.25
  USE define_dimensions
  USE define_types
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  CHARACTER(LEN=100), INTENT(IN)         :: filename_phen
  INTEGER(i_d), INTENT(IN)               :: mvt


  ! local variables
  integer np,nx,ilat
  integer, dimension(271,mvt)     :: greenup,fall,  phendoy1
!!integer, dimension(10)          :: greenupx,fallx,phendoy1x
!!integer, dimension(10)          :: ivtx
  integer, dimension(11)          :: greenupx,fallx,phendoy1x
  integer, dimension(11)          :: ivtx
  real(r_2), dimension(271)       :: xlat
  

  ! initialize for evergreen PFTs
  greenup(:,:) = -50
  fall(:,:)    = 367
  phendoy1(:,:)= 2


! do ilat=271,1,-1
!    read(101,*) xlat(ilat),(greenupx(nx),nx=1,10), &
!                           (fallx(nx),nx=1,10),    &
!                           (phendoy1x(nx),nx=1,10)
!    do nx=1,10
!       greenup(ilat,ivtx(nx)) = greenupx(nx)
!       fall(ilat,ivtx(nx))    = fallx(nx)
!       phendoy1(ilat,ivtx(nx))= phendoy1x(nx)
!    enddo
! enddo
!
! do np=1,mp
!
!    ilat=(casamet%lat(np)+55.25)/0.5+1
!    ilat= min(271,max(1,ilat))
!
!    phen%phase(np) = phendoy1(ilat,veg%iveg(np))
!    phen%doyphase(np,1) = greenup(ilat,veg%iveg(np))  ! DOY for greenup
!    phen%doyphase(np,2) = phen%doyphase(np,1) +14     ! DOY fo steady LAI
!    phen%doyphase(np,3) = fall(ilat,veg%iveg(np))     ! DOY for leaf senescence
!    phen%doyphase(np,4) = phen%doyphase(np,3) +14     ! DOY for minimal LAI season
!    if(phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
!    if(phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365
 
! Modifications for CLM (-melannie hartman 3/3/2014)
! The latitudes in modis_phenology start at 79.75 and count down by 0.5
! until reaching -55.25 (271 latitudes)

  open(101,file=filename_phen)
  read(101,*)

!!Found a bug on 1/26/2015 (-MDH).  There are 11 veg types in file  
!!modis_phenology_wtundra.txt since tundra (18) columns were added July 7, 2015.
!!
!!read(101,*) (ivtx(nx),nx=1,10)   ! fixed at 10, as only 10 of 17 IGBP PFT have seasonal leaf phenology
!!do ilat=1,271
!!
!!   read(101,*) xlat(ilat),(greenupx(nx),nx=1,10), &
!!                          (fallx(nx),nx=1,10),    &
!!                          (phendoy1x(nx),nx=1,10)
!!   do nx=1,10
!!      greenup(ilat,ivtx(nx)) = greenupx(nx)
!!      fall(ilat,ivtx(nx))    = fallx(nx)
!!      phendoy1(ilat,ivtx(nx))= phendoy1x(nx)
!!   enddo

  read(101,*) (ivtx(nx),nx=1,11)   ! fixed at 11, as only 11 of 18 IGBP PFT have seasonal leaf phenology
  do ilat=1,271
     read(101,*) xlat(ilat),(greenupx(nx),nx=1,11), &
                            (fallx(nx),nx=1,11),    &
                            (phendoy1x(nx),nx=1,11)
     do nx=1,11
        greenup(ilat,ivtx(nx)) = greenupx(nx)
        fall(ilat,ivtx(nx))    = fallx(nx)
        phendoy1(ilat,ivtx(nx))= phendoy1x(nx)
     enddo
  enddo

  do np=1,mp
     ilat=INT((79.75 - casamet%lat(np) + 0.25)/0.5) + 1
     !Any value of ilat < 1 or > 271 is out of range.
     ilat = max(ilat, 1)
     ilat = min(ilat, 271)

     phen%phase(np) = phendoy1(ilat,veg%iveg(np))
     phen%doyphase(np,1) = greenup(ilat,veg%iveg(np))  ! DOY for greenup
     phen%doyphase(np,2) = phen%doyphase(np,1) +14     ! DOY fo steady LAI
     phen%doyphase(np,3) = fall(ilat,veg%iveg(np))     ! DOY for leaf senescence
     phen%doyphase(np,4) = phen%doyphase(np,3) +14     ! DOY for minimal LAI season
     if(phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
     if(phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365

  enddo
END SUBROUTINE casa_readphen

!--------------------------------------------------------------------------------
SUBROUTINE casa_readpoint(filename_cnppoint,mvt)
   use define_dimensions
   use define_types
   use casaparm
   use casadimension
   use casavariable
   implicit none

   character(len=100) filename_cnppoint
   integer np,nland,mvt
   real(r_2)  annNdep,annNfix,annPwea,annPdust
   integer, dimension(mp) :: vtypex,stypex
   integer ivtigbp,doyp1,doyp2,doyp3,doyp4,phase
   integer IOstatus
     
   !! read gridinfo_ibgpz.csv

   open(101,file=filename_cnppoint)

   read(101,*) 

!  Added ilat(*) and ilon(*) to casamet so that casacnp output can be
!  remapped to a 2D grid (-Melannie Hartman 3/15/2014).

   do nland=1,mp
      !print *, 'nland = ', nland
      read(101,*,IOSTAT=IOstatus) casamet%ijgcm(nland),casamet%lat(nland),casamet%lon(nland), &
                  ivtigbp,stypex(nland),casamet%isorder(nland),casamet%areacell(nland), &
                  annNdep,annNfix,annPwea,annPdust,doyp1,doyp2,doyp3,doyp4,phase,vtypex(nland), &
                  casamet%ilat(nland), casamet%ilon(nland)
      if (IOstatus .lt. 0) then
         write(*,*)
         write(*,*) 'Unexpected EOF while reading file ', TRIM(filename_cnppoint)
         write(*,*) 'Expected ', mp, ' gridcells.  Found ', nland-1, 'gridcells.'
         STOP
      elseif (IOstatus .gt. 0) then
         print *, "An error occurred while reading file ", TRIM(filename_cnppoint)
         STOP
      endif
      if (casamet%isorder(nland) < 0) then
        write(*,*) 'Warning, iso has a negative value.'
        write(*,*) 'Check file ', trim(filename_cnppoint)
        STOP
      endif

      !! doyp1,doyp2,doyp3,doyp4,phase not used???
      !! initializing phenology from modis_phenology.txt instead

      casaflux%Nmindep(nland) = annNdep/365.0      ! This value will be over written by ndep in met.nc (-mdh 8/22/2016)
      casaflux%Pdep(nland)    = annPdust/365.0     ! gP/m2/day
      casaflux%Pwea(nland)    = annPwea/365.0      ! gP/m2/day
!!    Using 18 IGBP+tundra vegetation types (-mdh 7/7/2014)
!!    if(mvt==17) then
      if(mvt==18) then
         vtypex(nland)  = ivtigbp  ! for running IGBP veg type only
      endif
      if(vtypex(nland)==0) vtypex(nland)=iceland     
!! No longer reassigning tundra vegetation type (18) - MDH 7/7/2014
!!    if(vtypex(nland)==tundra) then
!!       vtypex(nland)=barren     ! tundra(18) is not an IGBP type defined in casacnp (-MDH 2/17/2014)
!!       print *, "Warning, reclassifing tundra (18) as barren (16) for point ", nland
!!     endif
!      print *, 'nland,vtypex,stypex ', nland,vtypex(nland),stypex(nland)
   enddo

   do np=1,mp
      veg%iveg(np)         = vtypex(np)
      soil%isoilm(np)      = stypex(np)
!      casaflux%Nminfix(np) = annNfix/365.0
! Commented out the statement "casaflux%Nminfix(np) = 0.0".  Nminfix is read from CASA PFT file (casa_readbiome). 
! This subroutine gets called multiple times in repeated transient runs (spinups), thus was
! eliminating Nfixation, contributing to plant N limitation and severely reduced NPP in global runs. -mdh 4/6/2020
!      casaflux%Nminfix(np) = 0.0
!      write(*,*) 'WARNING readpoint: casaflux%Nminfix(',np,') =', casaflux%Nminfix(np)
      if(veg%iveg(np)==12.or.veg%iveg(np)==14) casaflux%Pdep(np)=casaflux%Pdep(np)+0.7/365.0    ! for P fertilizer =13 Mt P globally in 1994
   enddo 
!   print * ,'veg type', veg%iveg   
   close(101)
   return

END SUBROUTINE casa_readpoint


!--------------------------------------------------------------------------------
SUBROUTINE casa_init(filename_cnpipool,mp,ms,mst)
!  initilize some values in phenology parameters and leaf growth phase 
use casadimension
use casaparm
use casavariable
use phenvariable
  implicit none
  character(len=100)                 filename_cnpipool
  integer  np,npt,npz,nl,ns,nland,nlandz,mp,ms,mst
  real(r_2) nyearz,ivtz,istz,latz,lonz,areacellz,glaiz,slaz,isoz

  print *, 'initcasa = ', initcasa
  !phen%phase = 2
IF (initcasa>=1) then
      write(*,*) 
      write(*,*) "Reading initial CASACNP pool file: ", filename_cnpipool, "..."
      open(99,file=filename_cnpipool)
      read(99,*)   ! Skip past file header. -mdh 1/30/2018
      do npt =1, mp
!! Commented out this section (-MDH 6/14/2014)
!!       Select Case(icycle)
!!       Case(1)
!!       !! Carbon only
!!          read(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,     &
!!          glaiz,slaz, casapool%clabile(npt),                        &
!!          casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:)
!!
!!       Case(2)
!!       !! C and N 
!!          read(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!          glaiz,slaz,casapool%clabile(npt),                     &
!!          casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:),      &
!!          casapool%nplant(npt,:),casapool%nlitter(npt,:), casapool%nsoil(npt,:),     &
!!          casapool%nsoilmin(npt)
!!       Case(3)
!!       !! C, N and P 
!!         read(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,                   &
!!          glaiz,slaz,casapool%clabile(npt),                                          &
!!          casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:),      &
!!          casapool%nplant(npt,:),casapool%nlitter(npt,:), casapool%nsoil(npt,:),     &
!!          casapool%nsoilmin(npt),                                                    &
!!          casapool%pplant(npt,:),casapool%plitter(npt,:), casapool%psoil(npt,:),     &
!!          casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt)
!!
!!       END Select 

!! Added this section (-MDH 6/14/2014)
    SELECT CASE(icycle)
      CASE(1)
        READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),  &
                   casapool%clitter(npt,:),casapool%csoil(npt,:)
!       write(*,123) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!                  casamet%glai(npt),slaz,phen%phase(npt), &
!                  casapool%clabile(npt),casapool%cplant(npt,:),  &
!                  casapool%clitter(npt,:),casapool%csoil(npt,:)
!       123 format(f4.0,',',i4,',',3(f4.0,','),5(f18.6,','),i6,',',4(f18.6,','))
 
      CASE(2)
        READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt)
!       write(*,124) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!                  casamet%glai(npt),slaz,phen%phase(npt), &
!                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!                  casapool%nsoil(npt,:),casapool%nsoilmin(npt)
!       124 format(f4.0,',',i4,',',3(f4.0,','),5(f18.6,','),i6,',',8(f18.6,','))

      CASE(3)
        READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
                   casapool%pplant(npt,:),casapool%plitter(npt,:),      &
                   casapool%psoil(npt,:),casapool%psoillab(npt),        &
                   casapool%psoilsorb(npt),casapool%psoilocc(npt)
!       write(*,125) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!                  casamet%glai(npt),slaz,phen%phase(npt), &
!                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!                  casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
!                  casapool%pplant(npt,:),casapool%plitter(npt,:),      &
!                  casapool%psoil(npt,:),casapool%psoillab(npt),        &
!                  casapool%psoilsorb(npt),casapool%psoilocc(npt)
!       125 format(f4.0,',',i4,',',3(f4.0,','),5(f18.6,','),i6,',',14(f18.6,','))

      END SELECT 
      enddo
    close(99)
    write(*,*) "Done reading initial CASACNP pool file: ", filename_cnpipool, "..."
endif
!  reset labile C pool
    casapool%clabile = 0.0    
!  check pool sizes
    casapool%cplant  = max(0.0,casapool%cplant)
    casapool%clitter = max(0.0,casapool%clitter)
    casapool%csoil   = max(0.0,casapool%csoil)
    casabal%cplantlast  = casapool%cplant
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil
    casabal%clabilelast = casapool%clabile
    casabal%sumcbal     = 0.0

!!  !! Commented out this section (-MDH 6/9/2014)
!!  if(icycle==1) then
!!      casapool%nplant(:,:)        = casapool%cplant(:,:) * casapool%rationcplant(:,:)
!!      casapool%Nsoil(:,:)        = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
!!      casapool%Psoil(:,:)        = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
!!  endif

    !! Added this section (-MDH 6/9/2014)
    IF (icycle==1) THEN
      casapool%Nplant(:,:) = casapool%cplant(:,:) * casapool%ratioNCplant(:,:)
      casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
      casapool%Psoil(:,:)  = casapool%Nsoil(:,:)/casapool%ratioNPsoil(:,:)
      casapool%Nsoilmin(:) = 2.5
    ENDIF 
    
    if(icycle >1) then
        casapool%nplant  = max(1.e-6,casapool%nplant)
        casapool%nlitter = max(1.e-6,casapool%nlitter)
        casapool%nsoil   = max(1.e-6,casapool%nsoil)
        casapool%nsoilmin= max(1.e-6,casapool%nsoilmin)
        casabal%nplantlast  = casapool%nplant
        casabal%nlitterlast = casapool%nlitter
        casabal%nsoillast   = casapool%nsoil       
        casabal%nsoilminlast= casapool%nsoilmin
        casabal%sumnbal     = 0.0
    endif

    if(icycle >2) then
        casapool%pplant = max(1.0e-7,casapool%pplant)
        casapool%plitter= max(1.0e-7,casapool%plitter)
        casapool%psoil  = max(1.0e-7,casapool%psoil)
        casapool%Psoillab  = max(2.0,casapool%psoillab)
        casapool%psoilsorb = max(10.0,casapool%psoilsorb)
        casapool%psoilocc  = max(50.0,casapool%psoilocc)
        casabal%pplantlast    = casapool%pplant
        casabal%plitterlast   = casapool%plitter
        casabal%psoillast     = casapool%psoil       
        casabal%psoillablast  = casapool%psoillab
        casabal%psoilsorblast = casapool%psoilsorb
        casabal%psoilocclast  = casapool%psoilocc
        casabal%sumpbal       = 0.0

    endif

end SUBROUTINE casa_init

!--------------------------------------------------------------------------------
! Write a header to the casacnp end-of-simulation pool (restart) file.
! This subroutine assumes that the calling routine has already opened the output
! file with unit number nout. The filename is for reference only.
!

SUBROUTINE write_cnpepool_header(icycle, nout, filename_cnpepool)
    implicit none
    integer, intent(in) :: icycle, nout
    character(len=100), intent(in) :: filename_cnpepool

    !                              10        20        30        40        50        60        70        80
    !                     123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    WRITE(nout,'(a72$)') 'iYrCnt,npt,veg%iveg,soil%isoilm,casamet%isorder,casamet%lat,casamet%lon,'
    WRITE(nout,'(a82$)') 'casamet%areacell,casamet%glai,casabiome%sla(veg%iveg),phen%phase,casapool%clabile,'
    WRITE(nout,'(a67$)') 'casapool%cplant(LEAF),casapool%cplant(WOOD),casapool%cplant(FROOT),'
    WRITE(nout,'(a67$)') 'casapool%clitter(METB),casapool%clitter(STR),casapool%clitter(CWD),'
    WRITE(nout,'(a62$)') 'casapool%csoil(MIC),casapool%csoil(SLOW),casapool%csoil(PASS),'
    WRITE(nout,'(a67$)') 'casapool%nplant(LEAF),casapool%nplant(WOOD),casapool%nplant(FROOT),'
    WRITE(nout,'(a67$)') 'casapool%nlitter(METB),casapool%nlitter(STR),casapool%nlitter(CWD),'
    WRITE(nout,'(a80$)') 'casapool%nsoil(MIC),casapool%nsoil(SLOW),casapool%nsoil(PASS),casapool%nsoilmin,'
    WRITE(nout,'(a67$)') 'casapool%pplant(LEAF),casapool%pplant(WOOD),casapool%pplant(FROOT),'
    WRITE(nout,'(a67$)') 'casapool%plitter(METB),casapool%plitter(STR),casapool%plitter(CWD),'
    WRITE(nout,'(a62$)') 'casapool%psoil(MIC),casapool%psoil(SLOW),casapool%psoil(PASS),'  
    WRITE(nout,'(a55$)') 'casapool%psoillab,casapool%psoilsorb,casapool%psoilocc,'
    WRITE(nout,'(a47)')  'casabal%sumcbal,casabal%sumnbal,casabal%sumpbal'

end SUBROUTINE write_cnpepool_header

!--------------------------------------------------------------------------------
! This subroutine assumes that the calling routine has already opened the output 
! file with unit number nout. The filename is for reference only.
!

SUBROUTINE write_cnpflux_header(icycle, nout, filename_cnpflux)
    implicit none
    integer, intent(in) :: icycle, nout
    character(len=100), intent(in) :: filename_cnpflux


! LEAF    = 1
! WOOD    = 2
! FROOT   = 3
! 
! METB    = 1
! STR     = 2
! CWD     = 3
! 
! MIC     = 1
! SLOW    = 2
! PASS    = 3
! 
! PLAB    = 1
! PSORB   = 2
! POCC    = 3

    Select Case(icycle)
    Case(1)
    !                                  10        20        30        40        50        60        70        80
    !                         123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        WRITE(nout,'(a47$)') 'myear,npt,veg%iveg,soil%isoilm,casamet%isorder,'
        WRITE(nout,'(a41$)') 'casamet%lat,casamet%lon,casamet%areacell,'
        WRITE(nout,'(a70$)') 'casabal%Fcnppyear,casabal%FCrsyear,casabal%FCneeyear,casabal%FCrpyear,'
        WRITE(nout,'(a58$)') 'clitterinput(LEAF),clitterinput(WOOD),clitterinput(FROOT),'
        WRITE(nout,'(a50)')  'csoilinput(METB),csoilinput(METB),csoilinput(METB)'
 
   Case(2)
    !                                  10        20        30        40        50        60        70        80
    !                         123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        WRITE(nout,'(a47$)') 'myear,npt,veg%iveg,soil%isoilm,casamet%isorder,'
        WRITE(nout,'(a41$)') 'casamet%lat,casamet%lon,casamet%areacell,'
        WRITE(nout,'(a70$)') 'casabal%Fcnppyear,casabal%FCrsyear,casabal%FCneeyear,casabal%FCrpyear,'
        WRITE(nout,'(a58$)') 'clitterinput(LEAF),clitterinput(WOOD),clitterinput(FROOT),'
        WRITE(nout,'(a51$)') 'csoilinput(METB),csoilinput(METB),csoilinput(METB),'
        WRITE(nout,'(a55$)') 'casabal%FNdepyear,casabal%FNfixyear,casabal%FNsnetyear,'
        WRITE(nout,'(a55)')  'casabal%FNupyear,casabal%FNleachyear,casabal%FNlossyear'
 
    !ATTENTION: Update for P
    ! Case(3)
        ! write(nout,??) myear,npt,veg%iveg(npt),soil%isoilm(npt),casamet%isorder(npt), &
        ! casamet%lat(npt),casamet%lon(npt),casamet%areacell(npt)*(1.0e-9), &
        ! casabal%FCnppyear(npt),casabal%FCrsyear(npt),casabal%FCneeyear(npt),casabal%FCrpyear(npt),&
        ! clitterinput(npt,:),csoilinput(npt,:), &
        ! casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt), &
        ! casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt), &
        ! casabal%FPweayear(npt),casabal%FPdustyear(npt), casabal%FPsnetyear(npt), &
        ! casabal%FPupyear(npt), casabal%FPleachyear(npt),casabal%FPlossyear(npt)
 
    END Select 

end SUBROUTINE write_cnpflux_header


!--------------------------------------------------------------------------------
SUBROUTINE casa_poolout(filename_cnpepool,iYrCnt,myear,writeToRestartCSVfile)
use define_dimensions
use define_types
use casadimension
use casaparm
use casavariable
use phenvariable
implicit none

  real(r_2), dimension(mso)      :: Psorder,pweasoil,fracPlab,fracPsorb,fracPocc,fracPorg,xpsoil50
  real(r_2), dimension(mp)       :: totpsoil


!Soiltype     soilnumber soil P(g P/m2)
!Alfisol    1    61.3
!Andisol    2    103.9
!Aridisol    3    92.8
!Entisol    4    136.9
!Gellisol    5    98.2
!Histosol    6    107.6
!Inceptisol    7    84.1
!Mollisol    8    110.1
!Oxisol        9    35.4    
!Spodosol    10    41.0    
!Ultisol    11    51.5    
!Vertisol    12    190.6
!Soil order specific parameters...
   data psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
   data pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/  !P weathering rate (gP/m2/yr)
   data fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
   data fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
   data fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
   data fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
   data xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

   integer  npt,nout,iYrCnt,nso,myear
   character(len=100) filename_cnpepool
   real(r_2) xyear
   logical writeToRestartCSVfile

!--------------------------------------------------------------------------
!  Added calculations for average annual pool values
   !xyear=1.0/(real(myear)*365)

   ! Output is every year now, not every myear years. -mdh 10/17/2016
   xyear = 1.0/365.0

   casapoolAn%CsoilAn=casapoolAn%CsoilAn * xyear
   casapoolAn%CplantAn=casapoolAn%CplantAn * xyear
   casapoolAn%ClitterAn=casapoolAn%ClitterAn * xyear

!  Added calculations for average temperatures. -MDH 08/17/2015
   casapoolAn%tairAn=casapoolAn%tairAn * xyear
   casapoolAn%tsoilAn=casapoolAn%tsoilAn * xyear

!  Added calculations for f(T) and f(W). -MDH 11/20/2017
   casapoolAn%fTAn=casapoolAn%fTAn * xyear
   casapoolAn%fWAn=casapoolAn%fWAn * xyear
   casapoolAn%thetaLiqAn=casapoolAn%thetaLiqAn * xyear

   if (icycle > 1) then
      casapoolAn%NsoilAn=casapoolAn%NsoilAn * xyear
      casapoolAn%NplantAn=casapoolAn%NplantAn * xyear
      casapoolAn%NlitterAn=casapoolAn%NlitterAn * xyear
      casapoolAn%NsoilminAn=casapoolAn%NsoilminAn * xyear
   endif

   if (icycle > 2) then
      casapoolAn%PsoilAn=casapoolAn%PsoilAn * xyear
      casapoolAn%PplantAn=casapoolAn%PplantAn * xyear
      casapoolAn%PlitterAn=casapoolAn%PlitterAn * xyear
   endif
!--------------------------------------------------------------------------

  if (writeToRestartCSVfile) then
     nout=103
     open(nout,file=filename_cnpepool)

    ! Write a header to the casacnp end-of-simulation pool (restart) file. -mdh 1/30/2018
    call write_cnpepool_header(icycle, nout, filename_cnpepool)

!      write(*,91) nyear,cplantsum,clittersum,csoilsum 
      casabal%sumcbal=min(9999.0,max(-9999.0,casabal%sumcbal))
      casabal%sumnbal=min(9999.0,max(-9999.0,casabal%sumnbal))
      casabal%sumpbal=min(9999.0,max(-9999.0,casabal%sumpbal))

!!    !! Commented out this section (-MDH 6/9/2014)
!!    do npt =1, mp
!!       nso = casamet%isorder(npt)
!!       totpsoil(npt) = psorder(nso) *xpsoil50(nso)
!!
!!       if(icycle<2) then
!!         casapool%nplant(npt,:)   = casapool%rationcplant(npt,:)  * casapool%cplant(npt,:)
!!         casapool%nlitter(npt,:)  = casapool%rationclitter(npt,:) * casapool%clitter(npt,:)
!!         casapool%nsoil(npt,:)    = casapool%ratioNCsoil(npt,:)   * casapool%Csoil(npt,:)
!!         casapool%nsoilmin(npt)   = 2.0
!!         casabal%sumnbal(npt)     = 0.0 
!!        endif
!!
!!        if(icycle<3) then
!!         casabal%sumpbal(npt)     = 0.0
!!         casapool%pplant(npt,:)   = casapool%ratiopcplant(npt,:)  * casapool%cplant(npt,:)
!!         casapool%plitter(npt,:)  = casapool%ratiopclitter(npt,:) * casapool%clitter(npt,:)
!!         casapool%psoil(npt,:)    = casapool%ratioPCsoil(npt,:)   * casapool%Csoil(npt,:)
!!         casapool%psoillab(npt)   = totpsoil(npt) *fracpLab(nso)
!!         casapool%psoilsorb(npt)  = casaflux%psorbmax(npt) * casapool%psoillab(npt) &
!!                                   /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
!!         casapool%psoilocc(npt)   = totpsoil(npt) *fracPocc(nso)
!!       endif
!!

  !!--------------------------------------------------------------------------------------
  !! Added this section (-MDH 6/9/2014)
  DO npt =1, mp
    nso = casamet%isorder(npt)
    totpsoil(npt) = psorder(nso) *xpsoil50(nso)
  if(casamet%iveg2(npt)>0 ) then
    IF (icycle<2) THEN
      casapool%Nplant(npt,:) = casapool%ratioNCplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%Nlitter(npt,:)= casapool%ratioNClitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%Nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0 
      if(casamet%iveg2(npt)==grass) then
         casapool%nplant(npt,wood) = 0.0
         casapool%nlitter(npt,cwd) = 0.0
      endif
    ENDIF 

    IF (icycle<3) THEN
      casabal%sumpbal(npt)    = 0.0
      casapool%pplant(npt,:)  = casapool%Nplant(npt,:)/casapool%ratioNPplant(npt,:)
      casapool%plitter(npt,:) = casapool%Nlitter(npt,:)/casapool%ratioNPlitter(npt,:)
      casapool%psoil(npt,:)   = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
      casapool%psoillab(npt)  = totpsoil(npt) *fracpLab(nso)
      casapool%psoilsorb(npt) = casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%psoilocc(npt)  = totpsoil(npt) *fracPocc(nso)
      if(casamet%iveg2(npt)==grass) then
         casapool%pplant(npt,wood) = 0.0
         casapool%plitter(npt,cwd) = 0.0
      endif
   ENDIF 
  else
     casapool%cplant(npt,:)=0.0; casapool%clitter(npt,:)=0.0; casapool%csoil(npt,:) = 0.0; casapool%clabile(npt) = 0.0
     casapool%nplant(npt,:)=0.0; casapool%nlitter(npt,:)=0.0; casapool%nsoil(npt,:) = 0.0; casapool%nsoilmin(npt) = 0.0
     casapool%pplant(npt,:)=0.0; casapool%plitter(npt,:)=0.0; casapool%psoil(npt,:) = 0.0
     casapool%psoillab(npt) = 0.0; casapool%psoilsorb(npt) = 0.0; casapool%psoilocc(npt) = 0.0
     casabal%sumcbal(npt) =0.0; casabal%sumnbal(npt) =0.0; casabal%sumpbal(npt) = 0.0
  endif
  !!--------------------------------------------------------------------------------------

!! ATTENTION: Check for consistency here! (-MDH 6/9/2014)
!!       write(nout,92) iYrCnt,npt,veg%iveg(npt),soil%isoilm(npt),casamet%isorder(npt), &
!!          casamet%lat(npt),casamet%lon(npt),casamet%areacell(npt)*(1.0e-9),          &
!!          casamet%glai(npt),casabiome%sla(veg%iveg(npt)), casapool%clabile(npt),     &
!!          casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:),      &
!!          casapool%nplant(npt,:),casapool%nlitter(npt,:), casapool%nsoil(npt,:),     &
!!          casapool%nsoilmin(npt),                                                    &
!!          casapool%pplant(npt,:),casapool%plitter(npt,:), casapool%psoil(npt,:),     &
!!          casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt),     &
!!          casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)

 !! phen%phase(npt) was added to the output list (-MDH 6/14/2014)
    WRITE(nout,92) iYrCnt,npt,veg%iveg(npt),soil%isoilm(npt),     &
          casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
          casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
          casabiome%sla(veg%iveg(npt)), phen%phase(npt), casapool%clabile(npt), &
          casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
          casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
          casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
          casapool%plitter(npt,:), casapool%psoil(npt,:),         &
          casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
          casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
      ENDDO

    CLOSE(nout)

  endif

91    format(i6,100(g9.3,2x))
!92    format(5(i6,',',2x),100(f15.6,',',2x))
!92  format(5(i6,3x),5(f18.6,3x),i6,3x,100(f18.6,3x))
92  format(5(i6,',',2x),5(f18.6,',',2x),i6,',',100(f18.6,',',2x))
end SUBROUTINE casa_poolout


!--------------------------------------------------------------------------------
SUBROUTINE casa_fluxout(filename_cnpflux,myear,clitterinput,csoilinput,writeToRestartCSVfile)
use define_dimensions
use define_types
use casadimension
use casaparm
use casavariable
use phenvariable
implicit none
  real(r_2), dimension(mp,3)   :: clitterinput,csoilinput
  integer  npt,nout,myear
  real(r_2) xyear
  character(len=100) filename_cnpflux
  logical :: writeToRestartCSVfile

   !xyear=1.0/real(myear)

   ! Output is every year now, not every myear years. -mdh 10/17/2016
   xyear = 1.0

   casabal%FCnppyear=casabal%FCnppyear * xyear
   casabal%FCgppyear=casabal%FCgppyear * xyear
   casabal%FCrsyear=casabal%FCrsyear * xyear
   casabal%FCrpyear=casabal%FCrpyear * xyear
   casabal%FCneeyear=casabal%FCneeyear * xyear
   casabal%FNdepyear=casabal%FNdepyear * xyear
   casabal%FNfixyear=casabal%FNfixyear * xyear
   casabal%FNsnetyear=casabal%FNsnetyear * xyear
   casabal%FNupyear=casabal%FNupyear * xyear
   casabal%FNleachyear=casabal%FNleachyear * xyear
   casabal%FNlossyear=casabal%FNlossyear * xyear
   casabal%FPweayear=casabal%FPweayear * xyear
   casabal%FPdustyear=casabal%FPdustyear * xyear
   casabal%FPsnetyear=casabal%FPsnetyear * xyear
   casabal%FPupyear=casabal%FPupyear * xyear
   casabal%FPleachyear=casabal%FPleachyear * xyear
   casabal%FPlossyear=casabal%FPlossyear * xyear
   clitterinput = clitterinput * xyear
   csoilinput   = csoilinput   * xyear


  if (writeToRestartCSVfile) then
      nout=104
      open(nout,file=filename_cnpflux)
      ! Write a header to the casacnp end-of-simulation flux file. -mdh 1/30/2018
      call write_cnpflux_header(icycle, nout, filename_cnpflux)

      do npt =1,mp
         Select Case(icycle)
         Case(1)
            write(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),casamet%isorder(npt), &
            casamet%lat(npt),casamet%lon(npt),casamet%areacell(npt)*(1.0e-9), &
            casabal%Fcnppyear(npt),casabal%FCrsyear(npt),casabal%FCneeyear(npt),casabal%FCrpyear(npt),&
            clitterinput(npt,:),csoilinput(npt,:)

         Case(2)
            write(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),casamet%isorder(npt), &
            casamet%lat(npt),casamet%lon(npt),casamet%areacell(npt)*(1.0e-9), &
            casabal%FCnppyear(npt),casabal%FCrsyear(npt),casabal%FCneeyear(npt),casabal%FCrpyear(npt),&
            clitterinput(npt,:),csoilinput(npt,:), &
            casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt), &
            casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)

         Case(3)
            write(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),casamet%isorder(npt), &
            casamet%lat(npt),casamet%lon(npt),casamet%areacell(npt)*(1.0e-9), &
            casabal%FCnppyear(npt),casabal%FCrsyear(npt),casabal%FCneeyear(npt),casabal%FCrpyear(npt),&
            clitterinput(npt,:),csoilinput(npt,:), &
            casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt), &
            casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt), &
            casabal%FPweayear(npt),casabal%FPdustyear(npt), casabal%FPsnetyear(npt), &
            casabal%FPupyear(npt), casabal%FPleachyear(npt),casabal%FPlossyear(npt)

         END Select 
      enddo

   close(nout)

  endif

92    format(5(i6,',',2x),100(f15.6,',',2x))
end SUBROUTINE casa_fluxout


!--------------------------------------------------------------------------------
SUBROUTINE casa_cnpflux()
use define_dimensions
use define_types
use casadimension
use casaparm
use casavariable
implicit none
!!real(r_2), dimension(mp,3)   :: clitterinput,csoilinput
!integer n

  casabal%FCgppyear(:)   = casabal%FCgppyear(:)   + casaflux%Cgpp(:)       * deltpool
  casabal%FCrpyear    = casabal%FCrpyear    + casaflux%Crp        * deltpool
  casabal%FCnppyear   = casabal%FCnppyear   + casaflux%Cnpp       * deltpool
  casabal%FCrsyear    = casabal%FCrsyear    + casaflux%Crsoil     * deltpool
  casabal%FCneeyear   = casabal%FCneeyear   + (casaflux%Cnpp-casaflux%Crsoil)* deltpool
 
! Commented out this section (-MDH 6/14/2014)
! do n=1,3
!    clitterinput(:,n)= clitterinput(:,n) + casaflux%kplant(:,n) * casapool%cplant(:,n) * deltpool
!    csoilinput(:,n)  = csoilinput(:,n)   + casaflux%fluxCtosoil(:,n) * deltpool
!    !csoilinput(:,n)  = csoilinput(:,n)   +  casaflux%fluxCtolitter(:,n) * deltpool
! enddo

! Added this section (-MDH 11/14/2016)
  casaflux%ClitInptMetAn = casaflux%ClitInptMetAn + casaflux%ClitInptMet
  casaflux%ClitInptStrucAn = casaflux%ClitInptStrucAn + casaflux%ClitInptStruc
  casaflux%CpassInptAn = casaflux%CpassInptAn + casaflux%CpassInpt

  ! Added this along with MIMICS-CN (-MDH 6/22/2019)
  casaflux%NlitInptMetAn = casaflux%NlitInptMetAn + casaflux%NlitInptMet
  casaflux%NlitInptStrucAn = casaflux%NlitInptStrucAn + casaflux%NlitInptStruc

  ! Added annual N fluxes (-MDH 11/25/2019)
  casaflux%NmindepAn = casaflux%NmindepAn + casaflux%Nmindep
  casaflux%NminfixAn = casaflux%NminfixAn + casaflux%Nminfix
  casaflux%NminleachAn = casaflux%NminleachAn + casaflux%Nminleach
  casaflux%NminlossAn = casaflux%NminlossAn + casaflux%Nminloss
  casaflux%NminuptakeAn = casaflux%NminuptakeAn + casaflux%Nminuptake
  casaflux%NlitterminAn = casaflux%NlitterminAn + casaflux%Nlittermin
  casaflux%NsminAn = casaflux%NsminAn + casaflux%Nsmin
  casaflux%NsimmAn = casaflux%NsimmAn + casaflux%Nsimm
  casaflux%NsnetAn = casaflux%NsnetAn + casaflux%Nsnet

  if(icycle >1) then
     casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep    * deltpool
     casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix    * deltpool
     casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet      * deltpool
     casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake * deltpool
     casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach  * deltpool
     casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss   * deltpool

  endif

  if(icycle >2) then
     casabal%FPweayear   = casabal%FPweayear   + casaflux%Pwea       * deltpool
     casabal%FPdustyear  = casabal%FPdustyear  + casaflux%Pdep       * deltpool
     casabal%FPsnetyear  = casabal%FPsnetyear  + casaflux%Psnet      * deltpool
     casabal%FPupyear    = casabal%FPupyear    + casaflux%Plabuptake * deltpool
     casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach     * deltpool  
     casabal%FPlossyear  = casabal%FPlossyear  + casaflux%Ploss      * deltpool 

  endif
END SUBROUTINE casa_cnpflux

!--------------------------------------------------------------------------------
!! casa_cnppool added to compute average annual pool values (-MDH 9/29/2014)
SUBROUTINE casa_cnppool()
use define_dimensions
use define_types
use casadimension
use casaparm
use casavariable
implicit none
integer n

casapoolAn%tairAn(:) = casapoolAn%tairAn(:) + casamet%tairk(:) - tkzeroc
casapoolAn%tsoilAn(:) = casapoolAn%tsoilAn(:) + casamet%tsoilavg(:) - tkzeroc
casapoolAn%fTAn(:) = casapoolAn%fTAn(:) + casapool%fT(:) 
casapoolAn%fWAn(:) = casapoolAn%fWAn(:) + casapool%fW(:)
casapoolAn%thetaLiqAn(:) = casapoolAn%thetaLiqAn(:) + casapool%thetaLiq(:)
casapoolAn%NsoilminAn(:) = casapoolAn%NsoilminAn(:) + casapool%Nsoilmin(:)

do n=1,3
  casapoolAn%CsoilAn(:,n) = casapoolAn%CsoilAn(:,n) + casapool%Csoil(:,n)
  casapoolAn%CplantAn(:,n) = casapoolAn%CplantAn(:,n) + casapool%Cplant(:,n)
  casapoolAn%ClitterAn(:,n) = casapoolAn%ClitterAn(:,n) + casapool%Clitter(:,n)

  if (icycle > 1) then
      casapoolAn%NsoilAn(:,n) = casapoolAn%NsoilAn(:,n) + casapool%Nsoil(:,n)
      casapoolAn%NplantAn(:,n) = casapoolAn%NplantAn(:,n) + casapool%Nplant(:,n)
      casapoolAn%NlitterAn(:,n) = casapoolAn%NlitterAn(:,n) + casapool%Nlitter(:,n)
  endif
   
  if (icycle > 2) then
      casapoolAn%PsoilAn(:,n) = casapoolAn%PsoilAn(:,n) + casapool%Psoil(:,n)
      casapoolAn%PplantAn(:,n) = casapoolAn%PplantAn(:,n) + casapool%Pplant(:,n)
      casapoolAn%PlitterAn(:,n) = casapoolAn%PlitterAn(:,n) + casapool%Plitter(:,n)
  endif

end do

END SUBROUTINE casa_cnppool

!-------------------------------------------------------------------------------- 
! Run all vegetation, litter, and soil subroutines for each point in the
! grid for one day. 
! Added arguments cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd
!   so the call to WritePointCASA call could be moved out of this subroutine and into 
!   casacnddriver. -mdh 5/14/2018.

SUBROUTINE biogeochem(iYrCnt,idoy,mdaily,nppScalar,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)

  USE define_dimensions
  USE define_types
  USE casadimension
  USE casa_cnp_module
  USE casa_nlim_module
  USE mimics_cycle_module
  USE corpse_cycle_module
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: iYrCnt, idoy, mdaily
  REAL(r_2), INTENT(IN)  :: nppScalar

  real(r_2), dimension(mp),INTENT(OUT)   :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd
  real(r_2), dimension(mp)               :: nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                                            pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,nwd2str

  ! local variables
  REAL(r_2),    DIMENSION(mp) :: xnplimit,xNPuptake
  REAL(r_2),    DIMENSION(mp) :: xklitter,xksoil,xkNlimiting
  REAL(r_2),    DIMENSION(mp) :: xkleafcold,xkleafdry,xkleaf
  INTEGER  j 

  xKNlimiting = 1.0
  call phenology(idoy,veg,phen)
  call avgsoil(veg,soil,casamet)
  call casa_rplant(veg,casabiome,casapool,casaflux,casamet)
  casaflux%Cnpp(:) = casaflux%Cnpp(:) * nppScalar      ! Change litter inputs if nppMult != 1 (nppMult set in co2delta.txt)
  call casa_allocation(veg,soil,casabiome,casaflux,casapool,casamet,phen)
  call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                       casamet,phen)

  casapool%fT(:) = 0.0 
  casapool%fW(:) = 0.0 
  casapool%thetaLiq(:) = 0.0 

  if (isomModel == CASACNP) then

      ! Run CASACNP SOM model

      call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                           casaflux,casamet)

      call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

      call casa_xratesoil(xklitter,xksoil,veg,soil,casamet,casabiome)
      call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)

      IF (icycle>1) THEN
          !call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
          ! Subroutine casa_xkN2 is similar to subroutine casa_xkN, but
          ! it does a simple linear ramp of xkNlimiting between 0.0 and 1.0
          ! when xkNlimitmin < casapool%Nsoilmin < xkNlimitmax. 
          ! (formerly when 0.5 < casapool%Nsoilmin < 2.0). 

          ! Added isomModel argument to the casa_xkN2 function call. -mdh 3/23/2020.
          !!call casa_xkN2(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
          call casa_xkN2(xkNlimiting,casapool,casaflux,casamet,casabiome,veg,isomModel)

          DO j=1,mlitter
              casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
          ENDDO
          call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
          IF (icycle >2) call casa_puptake(veg,xkNlimiting,casabiome, &
                                         casapool,casaflux,casamet)
      ENDIF 
    
      ! changed by ypwang following Chris Lu on 5/nov/2012
      ! Calculate the deltas for plant CNP, but don't update pools
      call casa_delplant(veg,casabiome,casapool,casaflux,casamet,            &
                         cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                         nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                         pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

      ! Calculate the deltas for soil CNP, but don't update pools
      call casa_delsoil(veg,casapool,casaflux,casamet,casabiome)
    
      ! Update plant and soil CNP pools
      call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)

  else if (isomModel == MIMICS) then
      ! Run MIMICS code

      mimicspool%fT(:) = 0.0 
      mimicspool%fW(:) = 0.0
      mimicspool%thetaLiq(:) = 0.0
      mimicspool%thetaFrzn(:) = 0.0

      call mimics_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                             casaflux,casamet)

      if (icycle > 1) then
          call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)
      endif 

      ! Compute klitter(mp,cwd) (-mdh 4/26/2015)
      call mimics_xratesoil(veg,soil,casamet,casabiome)

!! mimics_xratesoil replaces casa_xratesoil and casa_coeffsoil when MIMICS or CORPSE are used.
!! This call to casa_coeffsoil was not removed as it should have been, 
!! and was resetting casaflux%klitter(npt,cwd). -mdh 9/30/2019
!! Go over this code and determine what is needed for MIMICS-CN. -mdh 6/23/2019
!!    call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)
  
      IF (icycle>1) THEN
          !call casa_xkN(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
          ! Subroutine casa_xkN2 is similar to subroutine casa_xkN, but
          ! it does a simple linear ramp of xkNlimiting between 0.0 and 1.0
          ! when xkNlimitmin < casapool%Nsoilmin < xkNlimitmax. 
          ! (formerly when 0.5 < casapool%Nsoilmin < 2.0). 

          ! Added isomModel argument to the casa_xkN2 function call. -mdh 3/23/2020.
          !!call casa_xkN2(xkNlimiting,casapool,casaflux,casamet,casabiome,veg)
          call casa_xkN2(xkNlimiting,casapool,casaflux,casamet,casabiome,veg,isomModel)

          DO j=1,mlitter
              casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
!             write(*,*) 'biogeochem:'
!             write(*,'(a24,3(f6.4,2x))') 'xkNlimiting(1) = ', xkNlimiting(1)
!             write(*,'(a17,3(f10.6,2x))') 'casaflux%klitter(1,:) = ', casaflux%klitter(1,:)
          ENDDO
          call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
      ENDIF

!     write(*,*) 'iYrCnt = ', iYrCnt, '  DOY = ', idoy
      if (icycle == 1) then
        call mimics_delplant(veg,casabiome,casapool,casaflux,casamet,          &
                           cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                           nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                           pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                           cwd2co2,cwd2str)
      else
        call mimics_delplant_CN(veg,casabiome,casapool,casaflux,casamet,       &
                           cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                           nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                           pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                           cwd2co2,cwd2str,nwd2str)
      endif

      if (icycle == 1) then
          call mimics_soil_reverseMM(mp,iYrCnt,idoy,mdaily,cleaf2met,cleaf2str,&
                                     croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)
      else
          call mimics_soil_reverseMM_CN(mp,iYrCnt,idoy,mdaily,cleaf2met,cleaf2str, &
                                        croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
                                        nleaf2met,nleaf2str,nroot2met,nroot2str,nwd2str,nwood2cwd)
      endif

      if (icycle == 1) then
          call mimics_ccycle(veg,casabiome,casapool,casaflux,casamet)
      else
          call mimics_cncycle(veg,casabiome,casapool,casaflux,casamet)
      endif

      ! Accumulate C output variables
      call mimics_caccum(mp,cwd2co2)
      if (icycle > 1) then
          ! Accumulate N output variables
          call mimics_naccum(mp)
      endif

  else if (isomModel == CORPSE) then

      ! Run CORPSE soil code after vegetation subroutine, some are the same for
      ! mimics and corpse.

      call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                           casaflux,casamet)

      ! Compute klitter(mp,cwd) 
      call mimics_xratesoil(veg,soil,casamet,casabiome)

      call corpse_delplant(mp,veg,casabiome,casapool,casaflux,casamet,            &
                           cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,  &
                           nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,  &
                           pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd,  &
                           cwd2co2,cwd2str)

      call corpse_soil(mp,idoy,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2)

      ! Call mimics_ccycle for corpse as well
      call mimics_ccycle(veg,casabiome,casapool,casaflux,casamet)

  else 
      write(*,*) 'Unexpected value for isomModel =', isomModel
      write(*,*) 'Expecting 1, 2, or 3 in .lst file'
      STOP
  endif    

  ! Update Nplant and Pplant from Cplant and C:N and N:P ratios
  IF (icycle<3) then
     call casa_pdummy(casapool)
     IF (icycle<2) call casa_ndummy(casapool)
  ENDIF

  call casa_cnpbal(casapool,casaflux,casabal)

  call casa_cnpflux()

  ! Compute average annual soil, plant, and litter pools (-MDH 9/29/2014)
  call casa_cnppool()


END SUBROUTINE biogeochem
    
!----------------------------------------------------------------------------------------------------
SUBROUTINE GetMetNcFileDim(filename_cnpmet, ms, myear)

!  !DESCRIPTION
!   Get the dimensions of variables in the netCDF file to use for array allocation

      USE casadimension
      implicit none
      include 'netcdf.inc'

!  !ARGUMENTS
      character(len=100), intent(in) :: filename_cnpmet    ! NetCDF file name meteorological file
      integer, intent(in) :: ms                ! number of soil layers
      integer, intent(out) :: myear            ! number of years in the met file


! !LOCAL VARIABLES:
      integer :: ncid                ! netcdf file ID
      integer :: status                ! function return status
      integer :: myear_dimid            ! netcdf dimension id
      integer :: lat_dimid            ! netcdf dimension id
      integer :: lon_dimid            ! netcdf dimension id
      integer :: time_dimid            ! netcdf dimension id
      integer :: nsoil_dimid            ! netcdf dimension id
      integer :: nlon, nlat, ndays, nsoilyrs    ! Dimension sizes read from NetCDf file
      integer :: verbose=2

      if (verbose .ge. 0) print *, "Reading dimensions from met file ", trim(filename_cnpmet), "..."
     
      ! Get dimension ids

      status = nf_open(filename_cnpmet, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, "")
         
      status = nf_inq_dimid(ncid, "myear", myear_dimid)
      if (status /= nf_noerr) call handle_err(status, "myear")

      status = nf_inq_dimid(ncid, "lat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimid(ncid, "lon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_dimid(ncid, "time", time_dimid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_dimid(ncid, "nsoilyrs", nsoil_dimid)
      if (status /= nf_noerr) call handle_err(status, "nsoilyrs")


      ! Get dimension sizes

      status = nf_inq_dimlen(ncid, myear_dimid, myear)
      if (status /= nf_noerr) call handle_err(status, "myear")

      status = nf_inq_dimlen(ncid, lat_dimid, nlat)
      if (status /= nf_noerr) call handle_err(status, "nlat")

      status = nf_inq_dimlen(ncid, lon_dimid, nlon)
      if (status /= nf_noerr) call handle_err(status, "nlon")

      status = nf_inq_dimlen(ncid, time_dimid, ndays)
      if (status /= nf_noerr) call handle_err(status, "ndays")

      status = nf_inq_dimlen(ncid, nsoil_dimid, nsoilyrs)
      if (status /= nf_noerr) call handle_err(status, "nsoilyrs")

      if (verbose .ge. 1) then
         print *, "  GetMetNcFileDim: myear: ", myear
         print *, "  GetMetNcFileDim: nlat: ", nlat
         print *, "  GetMetNcFileDim: nlon: ", nlon
         print *, "  GetMetNcFileDim: time: ", ndays
         print *, "  GetMetNcFileDim: nsoilyrs: ", nsoilyrs
      endif

      !Make sure dimensions in netCDF file are consistent with casacnp

      if (ndays .ne. myear*mdyear) then
         print *, "Error in GetMetNcFileDim: ndays = ", ndays, "; myear*365 = ", myear*mdyear
         STOP
      endif

      if (nsoilyrs .ne. ms) then
         print *, "Error in GetMetNcFileDim: nsoilyrs = ", nsoilyrs, "; ms = ", ms
         STOP
      endif

      status = nf_close(ncid)


      if (verbose .ge. 0) print *, "Done reading dimensions from met file ", trim(filename_cnpmet)

END SUBROUTINE GetMetNcFileDim

!----------------------------------------------------------------------------------------------------
!SUBROUTINE ReadMetNcFile(filename_cnpmet, mp, ms, myear, xlai, xcnpp, xcgpp, xtairk, &
!                         xndepDay, xtsoil, xmoist, xfrznmoist)
SUBROUTINE ReadMetNcFile(filename_cnpmet, mp, ms, myear, xcgpp, xtairk, &
                         xndepDay, xtsoil, xmoist, xfrznmoist)
!  !DESCRIPTION
!   Read variables xcgpp, xtairk, xtsoil, xmoist, and xfrznmoist (if it exists) 
!     from the met.nc file.  
!   Dimensions of arrays must be known by first calling GetMetNcFileDim.
!     Note: I had lots of trouble with segmentation faults when I either tried to
!     allocate these arrays inside this function, or passed in arrays of unidicated 
!     sizes. Defining the dimensions in the argument list resolved the issue.
!     Melannie Hartman. Feb. 17, 2014

      USE casadimension
      USE clmgridvariable
      implicit none
      include 'netcdf.inc'

!  !ARGUMENTS
      character(len=100), intent(in) :: filename_cnpmet    ! NetCDF file name meteorological file
      integer, intent(in) :: mp                  ! number of points with valid data
      integer, intent(in) :: ms                  ! number of soil layers
      integer, intent(inout) :: myear            ! number of years in the met file
!     real(r_2), intent(inout), dimension(mp,mdyear*myear)    :: xlai       ! xlai(mp,ndays) daily LAI (m2/m2)
!     real(r_2), intent(inout), dimension(mp,mdyear*myear)    :: xcnpp      ! xcnpp(mp,ndays) daily NPP (gC/m2/day)
      real(r_2), intent(inout), dimension(mp,mdyear*myear)    :: xcgpp      ! xcgpp(mp,ndays) daily GPP (gC/m2/day)
      real(r_2), intent(inout), dimension(mp,mdyear*myear)    :: xndepDay   ! xndepDay(mp,ndays) daily N deposition (gN/m2/day)
      real(r_2), intent(inout), dimension(mp,mdyear*myear)    :: xtairk     ! xtairk(mp,ndays) daily average air temperature (K)
      real(r_2), intent(inout), dimension(mp,ms,mdyear*myear) :: xtsoil     ! xtsoil(mp,ms,ndays) daily soil temperature (K)
      real(r_2), intent(inout), dimension(mp,ms,mdyear*myear) :: xmoist     ! xmoist(mp,ms,ndays) daily volumetric liquid soil moisture (m3/m3)
      real(r_2), intent(inout), dimension(mp,ms,mdyear*myear) :: xfrznmoist ! xfrznmoist(mp,ms,ndays) daily volumetric frozen soil moisture (m3/m3)


! !LOCAL VARIABLES:
      integer :: i, j, cellCnt, npt, ims      ! Loop indices
      integer :: iday, ndoy1, ndoy2, nyear    ! Loop indices
      integer :: ncid                ! netcdf file ID
      integer :: status              ! function return status
      integer :: myear_dimid         ! netcdf dimension id
      integer :: lat_dimid           ! netcdf dimension id
      integer :: lon_dimid           ! netcdf dimension id
      integer :: time_dimid          ! netcdf dimension id
      integer :: nsoil_dimid         ! netcdf dimension id
      integer :: nlon, nlat, ndays, nsoilyrs     ! Dimension sizes read from NetCDf file
      integer :: varid                           ! netcdf variable id
      integer :: start1(1), count1(1)            ! start and count arrays for reading 1-D data from netcdf files
      integer :: start2(2), count2(2)            ! start and count arrays for reading 2-D data from netcdf files
      integer :: start3(3), count3(3)            ! start and count arrays for reading 3-D data from netcdf files
      integer :: start4(4), count4(4)            ! start and count arrays for reading 4-D data from netcdf files
      !Arrays read from NetCDF file (including cells with missing data)
      real(4), allocatable :: lai(:,:,:)         ! lai(nlon,nlat,ndays) total projected LAI
      real(4), allocatable :: Cnpp(:,:,:)        ! Cnpp(nlon,nlat,ndays) Net Primary Production (gC/m2/day)
      real(4), allocatable :: Cgpp(:,:,:)        ! Cgpp(nlon,nlat,ndays) Gross Primary Production (gC/m2/day)
      real(4), allocatable :: tairk(:,:,:)       ! tairk(nlon,nlat,ndays) daily average air temperature (K) 
      real(4), allocatable :: ndep(:,:,:)        ! ndep(nlon,nlat,ndays) daily N deposition (gN/m2/day)
      real(4), allocatable :: tsoil(:,:,:,:)     ! tsoil(nlon,nlat,ms,ndays), CASACNP soil temperature by layer (K)
      real(4), allocatable :: moist(:,:,:,:)     ! moist(nlon,nlat,ms,ndays), CASACNP volumetric soil liquid water content by layer 
      real(4), allocatable :: frznmoist(:,:,:,:) ! frznmoist(nlon,nlat,ms,ndays), CASACNP volumetric soil frozen water content by layer 
      integer :: verbose=2

      !Save the value of myear for error checking later
      nyear = myear

      if (verbose .ge. 0) print *, "Reading met file ", trim(filename_cnpmet), "..."
     
      ! Get dimension ids

      status = nf_open(filename_cnpmet, nf_nowrite, ncid)
      if (status /= nf_noerr) call handle_err(status, "")
         
      status = nf_inq_dimid(ncid, "myear", myear_dimid)
      if (status /= nf_noerr) call handle_err(status, "myear")

      status = nf_inq_dimid(ncid, "lat", lat_dimid)
      if (status /= nf_noerr) call handle_err(status, "lat")

      status = nf_inq_dimid(ncid, "lon", lon_dimid)
      if (status /= nf_noerr) call handle_err(status, "lon")

      status = nf_inq_dimid(ncid, "time", time_dimid)
      if (status /= nf_noerr) call handle_err(status, "time")

      status = nf_inq_dimid(ncid, "nsoilyrs", nsoil_dimid)
      if (status /= nf_noerr) call handle_err(status, "nsoilyrs")


      ! Get dimension sizes

      status = nf_inq_dimlen(ncid, myear_dimid, myear)
      if (status /= nf_noerr) call handle_err(status, "myear")

      status = nf_inq_dimlen(ncid, lat_dimid, nlat)
      if (status /= nf_noerr) call handle_err(status, "nlat")

      status = nf_inq_dimlen(ncid, lon_dimid, nlon)
      if (status /= nf_noerr) call handle_err(status, "nlon")

      status = nf_inq_dimlen(ncid, time_dimid, ndays)
      if (status /= nf_noerr) call handle_err(status, "ndays")

      status = nf_inq_dimlen(ncid, nsoil_dimid, nsoilyrs)
      if (status /= nf_noerr) call handle_err(status, "nsoilyrs")

      if (verbose .ge. 1) then
         print *, "  ReadMetNcFile: myear: ", myear
         print *, "  ReadMetNcFile: nlat: ", nlat
         print *, "  ReadMetNcFile: nlon: ", nlon
         print *, "  ReadMetNcFile: time: ", ndays
         print *, "  ReadMetNcFile: nsoilyrs: ", nsoilyrs
      endif

      !Make sure dimensions in netCDF file are consistent with casacnp

      if (nyear .ne. myear) then
         print *, "Error in ReadMetNcFile: nyear = ", nyear, "; myear = ", myear
         STOP
      endif

      if (ndays .ne. myear*mdyear) then
         print *, "Error in ReadMetNcFile: ndays = ", ndays, "; myear*365 = ", myear*mdyear
         STOP
      endif

      if (nsoilyrs .ne. ms) then
         print *, "Error in ReadMetNcFile: nsoilyrs = ", nsoilyrs, "; ms = ", ms
         STOP
      endif
 
      !Allocate CLM grid variables: 
      !  clmgrid%lon1d(nlon) coordinate longitude (degrees east)
      !  clmgrid%lat1d(nlat) coordinate latitude (degrees north)
      !  clmgrid%cellMissing(nlon,nlat) 0=no missing data, 1=missing data
      !  clmgrid%cellid(nlon,nlat) grid cell ids (1..nlat*nlon)
      call alloc_CLMgridVariable(nlon, nlat)

      !Allocate local variables read from NetCDF file (full grid including cells with missing data)
      allocate(lai(1:nlon,1:nlat,1:ndays))
      allocate(Cnpp(1:nlon,1:nlat,1:ndays))
      allocate(Cgpp(1:nlon,1:nlat,1:ndays))
      allocate(tairk(1:nlon,1:nlat,1:ndays))
      allocate(ndep(1:nlon,1:nlat,1:ndays))
      allocate(tsoil(1:nlon,1:nlat,1:ms,1:ndays))
      allocate(moist(1:nlon,1:nlat,1:ms,1:ndays))
      allocate(frznmoist(1:nlon,1:nlat,1:ms,1:ndays))

      !Get variable ids, then data arrays

      ! clmgrid%lon1d(nlon) 

      start1 = (/ 1 /)
      count1 = (/ nlon /)

      if (verbose .ge. 1) print *, "  Reading longitude..."
      status = nf_inq_varid(ncid, "lon", varid)
      if (status /= nf_noerr) call handle_err(status, "clmgrid%lon1d")

      status = nf_get_var(ncid, varid, clmgrid%lon1d, start1, count1)
      if (status /= nf_noerr) call handle_err(status, "clmgrid%lon1d")


      ! clmgrid%lat1d(nlat) 

      start1 = (/ 1 /)
      count1 = (/ nlat /)

      if (verbose .ge. 1) print *, "  Reading latitude..."
      status = nf_inq_varid(ncid, "lat", varid)
      if (status /= nf_noerr) call handle_err(status, "lat1d")

      status = nf_get_var(ncid, varid, clmgrid%lat1d, start1, count1)
      if (status /= nf_noerr) call handle_err(status, "clmgrid%lat1d")


      start2 = (/ 1, 1 /)
      count2 = (/ nlat, nlon /)

      ! clmgrid%cellMissing(nlon, nlat) 

      if (verbose .ge. 1) print *, "  Reading cellMissing..."
      status = nf_inq_varid(ncid, "cellMissing", varid)
      if (status /= nf_noerr) call handle_err(status, "cellMissing")

      status = nf_get_var(ncid, varid, clmgrid%cellMissing, start2, count2)
      if (status /= nf_noerr) call handle_err(status, "clmgrid%cellMissing")


      ! clmgrid%cellid(nlon, nlat) 

      if (verbose .ge. 1) print *, "  Reading cellid..."
      status = nf_inq_varid(ncid, "cellid", varid)
      if (status /= nf_noerr) call handle_err(status, "cellid")

      status = nf_get_var(ncid, varid, clmgrid%cellid, start2, count2)
      if (status /= nf_noerr) call handle_err(status, "clmgrid%cellid")


      start3 = (/ 1, 1 ,1 /)
      count3 = (/ ndays, nlat, nlon /)

      ! lai(nlon, nlat, ndays) 

!     if (verbose .ge. 1) print *, "  Reading xlai..."
!     status = nf_inq_varid(ncid, "xlai", varid)
!     if (status /= nf_noerr) call handle_err(status, "xlai")
!
!     status = nf_get_var(ncid, varid, lai, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "xlai")


      ! Cnpp(nlon, nlat, ndays) 
 
!     if (verbose .ge. 1) print *, "  Reading xcnpp..."
!     status = nf_inq_varid(ncid, "xcnpp", varid)
!     if (status /= nf_noerr) call handle_err(status, "xcnpp")
!
!     status = nf_get_var(ncid, varid, Cnpp, start3, count3)
!     if (status /= nf_noerr) call handle_err(status, "xcnpp")
 

      ! cgpp(nlon, nlat, ndays) 
 
      if (verbose .ge. 1) print *, "  Reading xcgpp..."
      status = nf_inq_varid(ncid, "xcgpp", varid)
      if (status /= nf_noerr) call handle_err(status, "xcgpp")
 
      status = nf_get_var(ncid, varid, cgpp, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "xcgpp")

 
      ! tairk(nlon, nlat, ndays) 
 
      if (verbose .ge. 1) print *, "  Reading xtairk..."
      status = nf_inq_varid(ncid, "xtairk", varid)
      if (status /= nf_noerr) call handle_err(status, "xtairk")
 
      status = nf_get_var(ncid, varid, tairk, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "xtairk")


      ! ndep(nlon, nlat, ndays) 
 
      if (verbose .ge. 1) print *, "  Reading ndep..."
      status = nf_inq_varid(ncid, "ndep", varid)
      if (status /= nf_noerr) call handle_err(status, "ndep")
 
      status = nf_get_var(ncid, varid, ndep, start3, count3)
      if (status /= nf_noerr) call handle_err(status, "ndep")
 
 
      start4 = (/ 1, 1 ,1, 1 /)
      count4 = (/ ndays, ms, nlat, nlon /)
 
 
      ! tsoil(nlon, nlat, ms, ndays) volumetric soil water content (mm3/mm3)
 
      if (verbose .ge. 1) print *, "  Reading xtsoil..."
      status = nf_inq_varid(ncid, "xtsoil", varid)
      if (status /= nf_noerr) call handle_err(status, "xtsoil")
 
      status = nf_get_var(ncid, varid, tsoil, start4, count4)
      if (status /= nf_noerr) call handle_err(status, "xtsoil")
 
 
      ! moist(nlon, nlat, ms, ndays) volumetric soil liquid water content (mm3/mm3)
 
      if (verbose .ge. 1) print *, "  Reading xmoist..."
      status = nf_inq_varid(ncid, "xmoist", varid)
      if (status /= nf_noerr) call handle_err(status, "xmoist")
 
      status = nf_get_var(ncid, varid, moist, start4, count4)
      if (status /= nf_noerr) call handle_err(status, "xmoist")

      ! frznmoist(nlon, nlat, ms, ndays) volumetric soil frozen water content (mm3/mm3)
      ! This variable will be reamin zero if it does not exist in the met.nc file. -mdh 3/13/2017
 
      if (verbose .ge. 1) print *, "  Reading xfrznmoist..."
      status = nf_inq_varid(ncid, "xfrznmoist", varid)
!     if (status /= nf_noerr) call handle_err(status, "xfrznmoist")
      if (status == nf_noerr) then
          status = nf_get_var(ncid, varid, frznmoist, start4, count4)
          if (status /= nf_noerr) call handle_err(status, "xfrznmoist")
      else
          write(*,*) '    WARNING: xfrznmoist does not exist in the file ', trim(filename_cnpmet)
          write(*,*) '    Soil frozen moisture will be zet to 0.0'
      endif
 
 
     ! Because dimensions in FORTRAN are in Column Major Order (the first 
     ! array index varies the most rapidly) dimensions of FORTRAN arrays 
     ! are in the opposite order that they appear in the NetCDF file with ncdump. 
 
      cellCnt = 0
      npt = 0
      !Count the number of grid cells with valid data
      !Optimize here. Use COUNT function. -mdh 7/22/2019
!     do i = 1, nlon
!        do j = 1, nlat
!           if (clmgrid%cellMissing(i,j) .eq. 0) then
!              npt = npt + 1
!           endif
!        enddo
!     enddo 
      npt = COUNT(clmgrid%cellMissing(1:nlon,1:nlat) == 0)

      !Make sure dimensions in netCDF file are consistent with casacnp
      if (mp .ne. npt) then
         print *, "Error in ReadMetNcFile: mp = ", mp, "; npt = ", npt
         STOP
      endif

      cellCnt = 0
      npt = 0
      !Fill xcgpp, xtairk, xndepDay, xtsoil, xmoist, and xfrznmoist with non-missing data only
      !Optimize here? Switch order of do i,j loops? Can't because of npt.
      do i = 1, nlon
         do j = 1, nlat
            cellCnt = cellCnt+1
            if (clmgrid%cellid(i,j) .ne. cellCnt) then
               print *, "Data alignment problem in ReadMetNcFile: "
               print *, "cellid = ", clmgrid%cellid(i,j), " cellCnt = ", cellCnt
            endif
            if (clmgrid%cellMissing(i,j) .eq. 0) then
               npt = npt + 1
               do iday = 1, ndays
                  ! Segmentation fault was occurring on next assignment until
                  ! I passed in the dimensions to these variables. -MDH 2/17/2014
                  !xlai(npt,iday) = lai(i,j,iday)
                  !xcnpp(npt,iday) = Cnpp(i,j,iday)
                  xcgpp(npt,iday) = Cgpp(i,j,iday)
                  xtairk(npt,iday) = tairk(i,j,iday)
                  xndepDay(npt,iday) = ndep(i,j,iday)
                  do ims = 1, nsoilyrs
                     xtsoil(npt,ims,iday) = tsoil(i,j,ims,iday)
                     xmoist(npt,ims,iday) = moist(i,j,ims,iday)
                     !write(*,*) 'xfrznmoist(',npt,ims,iday,') = frznmoist(',i,j,ims,iday,')'
                     xfrznmoist(npt,ims,iday) = frznmoist(i,j,ims,iday)
                  enddo
               enddo
            endif
         enddo
      enddo 

      status = nf_close(ncid)

      if (verbose .ge. 0) print *, "Done reading met file ", trim(filename_cnpmet), "..."

END SUBROUTINE ReadMetNcFile

!----------------------------------------------------------------------------------------------------
   SUBROUTINE handle_err(status, errmsg)
      include 'netcdf.inc'
      integer, intent(in) :: status
      character(len=*), intent(in)    :: errmsg  ! append error message
     
      if (status /= nf_noerr) then
         print *, trim(nf_strerror(status)), ": ", errmsg
         stop "Stopped"
      endif
   end SUBROUTINE handle_err

!----------------------------------------------------------------------------------------------------
   SUBROUTINE get_time_and_date(date_string, time_string)

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

   end SUBROUTINE get_time_and_date

!----------------------------------------------------------------------------------------------------
   SUBROUTINE PutVariableAttributeReal(ncid, varid, attr_name, attr_units, missing_value)

    implicit none
    include 'netcdf.inc'

!   ARGUMENTS
      integer, intent(in) :: ncid                     ! netcdf file id
      integer, intent(in) :: varid                    ! netcdf variable id
      character(len=100), intent(in) :: attr_name     ! String for assigning variable "long_name" attribute
      character(len=100), intent(in) :: attr_units    ! String for assigning variable "units" attribute
      real(4), intent(in) :: missing_value
      
!   LOCAL VARIABLES
      integer status    ! NetCDF error status

!     write(*,*)
!     write(*,*) 'PutVariableAttributeReal: ncid = ', ncid
!     write(*,*) 'PutVariableAttributeReal: varid = ', varid
!     write(*,*) 'PutVariableAttributeReal: attr_name = ', trim(attr_name)
!     write(*,*) 'PutVariableAttributeReal: attr_units = ', trim(attr_units)
!     write(*,*) 'PutVariableAttributeReal: missing_value = ', missing_value

      status = nf_put_att_text(ncid, varid, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
 
      status = nf_put_att_text(ncid, varid, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
 
      status = nf_put_att_real(ncid, varid, '_FillValue', NF_REAL, 1, missing_value)
      if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
      status = nf_put_att_real(ncid, varid, 'missing_value', NF_REAL, 1, missing_value)
      if (status /= nf_noerr) call handle_err(status, "missing_value")

   end SUBROUTINE PutVariableAttributeReal

!----------------------------------------------------------------------------------------------------
   SUBROUTINE PutVariableAttributeDouble(ncid, varid, attr_name, attr_units, missing_value)

    implicit none
    include 'netcdf.inc'

!   ARGUMENTS
      integer, intent(in) :: ncid                      ! netcdf file id
      integer, intent(in) :: varid                     ! netcdf variable id
      character(len=100), intent(in) :: attr_name      ! String for assigning variable "long_name" attribute
      character(len=100), intent(in) :: attr_units     ! String for assigning variable "units" attribute
      real(8), intent(in) :: missing_value
      
!   LOCAL VARIABLES
      integer status    ! NetCDF error status

      status = nf_put_att_text(ncid, varid, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
 
      status = nf_put_att_text(ncid, varid, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
 
      status = nf_put_att_double(ncid, varid, '_FillValue', NF_DOUBLE, 1, missing_value)
      if (status /= nf_noerr) call handle_err(status, "_FillValue")
 
      status = nf_put_att_double(ncid, varid, 'missing_value', NF_DOUBLE, 1, missing_value)
      if (status /= nf_noerr) call handle_err(status, "missing_value")
 
   END SUBROUTINE PutVariableAttributeDouble

!----------------------------------------------------------------------------------------------------
    SUBROUTINE WritePoolFluxNcFile_casacnp_annual(filename_ncOut_casa, mp, year)
!
!   DESCRIPTION
!   Create, define, and write mean annual output to netcdf file filename_ncOut_casa.  
!   The file contains mean annual output for one year.
!   Called by casacnpdriver.
!
!   Melannie Hartman. March 9, 2014

    USE casadimension
    USE clmgridvariable
    USE casavariable
    USE casaparm
    USE define_types
    implicit none 
    include 'netcdf.inc'
    real(4), parameter :: MISSING_VALUE = 1.e+36
    integer, parameter :: MISSING_INT = -9999

!   ARGUMENTS
      character(len=*), intent(in) :: filename_ncOut_casa    ! NetCDF output file name 
      integer, intent(in) :: mp                              ! number of points with valid data
      integer, intent(in) :: year                            ! output year
!   LOCAL VARIABLES
      !integer :: i
      integer :: ncid                  ! netcdf file ID
      integer :: status                ! function return status
!     integer :: dimid_mp              ! netcdf dimension id
      integer :: dimid_lat             ! netcdf dimension id
      integer :: dimid_lon             ! netcdf dimension id
      integer :: dimid_time            ! netcdf dimension id
      integer :: nlon, nlat, ntimes    ! Dimension sizes for NetCDf file
      integer :: dims(3)               ! Array of NetCDF dimension IDs for defining variables
      integer :: start1(1), count1(1)  ! start and count arrays for writing 1-D data from netcdf files
      !integer :: start2(2), count2(2) ! start and count arrays for writing 2-D data from netcdf files
      integer :: start3(3), count3(3)  ! start and count arrays for writing 3-D data from netcdf files
      integer :: varid_lon, varid_lat  ! NetCDF variable ID for latitude and longitude
      integer :: varid_time            ! NetCDF variable ID for time
!     integer :: varid_year            ! NetCDF variable ID for year
      integer :: varid_mask            ! NetCDF variable ID for cellMissing(nlon,nlat)
      integer :: varid_cellid          ! NetCDF variable ID for cellid(nlon,nlat)
      integer :: varid_igbp            ! NetCDF variable ID for IGBP_PFT(nlon,nlat)
      integer :: varid_landarea        ! NetCDF variable ID for land area
      integer :: varid_cleaf           ! NetCDF variable ID for leaf C
      integer :: varid_cwood           ! NetCDF variable ID for wood C
      integer :: varid_cfroot          ! NetCDF variable ID for fine root C
      integer :: varid_nleaf           ! NetCDF variable ID for leaf N
      integer :: varid_nwood           ! NetCDF variable ID for wood N
      integer :: varid_nfroot          ! NetCDF variable ID for fine root N
      integer :: varid_clitmetb        ! NetCDF variable ID for metablic litter C
      integer :: varid_clitstr         ! NetCDF variable ID for structural litter C
      integer :: varid_clitcwd         ! NetCDF variable ID for coarse weeody debris C
      integer :: varid_nlitmetb        ! NetCDF variable ID for metablic litter N
      integer :: varid_nlitstr         ! NetCDF variable ID for structural litter N
      integer :: varid_nlitcwd         ! NetCDF variable ID for coarse weeody debris N
      integer :: varid_csoilmic        ! NetCDF variable ID for microbial soil C
      integer :: varid_csoilslow       ! NetCDF variable ID for slow soil C
      integer :: varid_csoilpass       ! NetCDF variable ID for passive soil C
      integer :: varid_nsoilmic        ! NetCDF variable ID for microbial soil N
      integer :: varid_nsoilslow       ! NetCDF variable ID for slow soil N
      integer :: varid_nsoilpass       ! NetCDF variable ID for passive soil N
      integer :: varid_cnpp            ! NetCDF variable ID for NPP C
      integer :: varid_cgpp            ! NetCDF variable ID for GPP C
      integer :: varid_cresp           ! NetCDF variable ID for soil respiration C
      integer :: varid_tairC           ! NetCDF variable ID for air temperature (C)
      integer :: varid_tsoilC          ! NetCDF variable ID for soil temperature (C)
      integer :: varid_clitInptMet     ! NetCDF variable ID for metabolic litter C inputs
      integer :: varid_clitInptStruc   ! NetCDF variable ID for structural litter C inputs
      integer :: varid_cpassInpt       ! NetCDF variable ID for passive pool C inputs
      integer :: varid_fT, varid_fW    ! NetCDF variable ID soil temperature and moisture decomposition rate multipliers
      integer :: varid_thetaLiq        ! NetCDF variable ID liquid soil moisture fraction
      !! Added N fluxes and pools -mdh 11/25/2019
      integer :: varid_nlitInptMet     ! NetCDF variable ID for metabolic litter N inputs
      integer :: varid_nlitInptStruc   ! NetCDF variable ID for structural litter N inputs
      integer :: varid_nmindep         ! NetCDF variable ID for N deposition
      integer :: varid_nminfix         ! NetCDF variable ID for N fixation
      integer :: varid_nminuptake      ! NetCDF variable ID for plant N uptake
      integer :: varid_nminleach       ! NetCDF variable ID for mineral N leaching
      integer :: varid_nminloss        ! NetCDF variable ID for mineral N loss (N2O)
      integer :: varid_nlittermin      ! NetCDF variable ID for litter N mineralization
      integer :: varid_nsmin           ! NetCDF variable ID for soil N mineralization
      integer :: varid_nsimm           ! NetCDF variable ID for soil N immobilization
      integer :: varid_nsnet           ! NetCDF variable ID for net N mineralization
      integer :: varid_nminrl          ! NetCDF variable ID for mineral soil N
      character(len=100) :: attr_name   ! String for assigning global and variable attributes
      character(len=100) :: attr_units  ! String for assigning global and variable attributes
      character(len=10)  :: date_string ! String for assigning date to global attributes
      character(len=8)   :: time_string ! String for assigning time to global attributes
      integer :: verbose=0
      integer :: npt, ilon, ilat, itime
      integer, allocatable :: IGBP_PFT(:,:)  ! IGBP_PFT(nlon,nlat) IGBP PFT classification (1-18)
      real(4), allocatable :: landarea(:,:)  ! landarea(nlon,nlat) km^2
      real(4), allocatable :: var1(:,:,:)    ! gridded output variable
      real(4), allocatable :: var2(:,:,:)    ! gridded output variable
      real(4), allocatable :: var3(:,:,:)    ! gridded output variable
      real(4), allocatable :: var4(:,:,:)    ! gridded output variable
      real(4), allocatable :: var5(:,:,:)    ! gridded output variable
      real(4), allocatable :: var6(:,:,:)    ! gridded output variable
      real(4), allocatable :: var7(:,:,:)    ! gridded output variable
      real(4), allocatable :: var8(:,:,:)    ! gridded output variable
      real(4) :: time                        ! time to write to netcdf file
    

      if (verbose .ge. 0) print *, "Writing output to file ", trim(filename_ncOut_casa), "..."

!Output these variables in NetCDF file:
!       veg%iveg(npt)           - IGBP PFTs
!       casamet%lat(npt)        - latitudes    
!       casamet%lon(npt)        - longitudes
!       casamet%areacell(npt)*(1.0e-6) - landarea (km^2)
!       casabal%FCnppyear(npt)         - NPP (gC/m2/yr) 
!       casabal%FCrsyear(npt)          - soil respiration (gC/m2/yr)        
!       casabal%FCrpyear(npt)          - plant respiration (gC/m2/yr)    
!       casapool%cplantAn(npt,LEAF)    - plant leaf C (gC/m2)
!       casapool%cplantAn(npt,WOOD)    - plant wood C (gC/m2)
!       casapool%cplantAn(npt,FROOT)   - plant fine root C (gC/m2)
!       casapool%clitterAn(npt,METB)   - metabolic litter C (gC/m2)
!       casapool%clitterAn(npt,STR)    - structural litter C (gC/m2)
!       casapool%clitterAn(npt,CWD)    - coarse woody debris C (gC/m2)
!       casapool%csoilAn(npt,MIC)      - microbial soil  C (gC/m2)
!       casapool%csoilAn(npt,SLOW)     - slow soil C (gC/m2)
!       casapool%csoilAn(npt,PASS)     - passive soil C (gC/m2)
!       casapool%nplantAn(npt,LEAF)    - plant leaf N (gN/m2)
!       casapool%nplantAn(npt,WOOD)    - plant wood N (gN/m2)
!       casapool%nplantAn(npt,FROOT)   - plant fine root N (gN/m2)
!       casapool%nlitterAn(npt,METB)   - metabolic litter N (gN/m2)
!       casapool%nlitterAn(npt,STR)    - structural litter N (gN/m2)
!       casapool%nlitterAn(npt,CWD)    - coarse woody debris N (gN/m2)
!       casapool%nsoilAn(npt,MIC)      - microbial soil N (gN/m2)
!       casapool%nsoilAn(npt,SLOW)     - slow soil N (gN/m2)
!       casapool%nsoilAn(npt,PASS)     - passive soil N (gN/m2)

   dims(1) = 0
   dims(2) = 0
   ntimes = 1
   nlat = clmgrid%nlat
   nlon = clmgrid%nlon


   ! Create the netcdf file 

   status = nf_create(filename_ncOut_casa, NF_CLOBBER, ncid)
   if (status /= nf_noerr) call handle_err(status, trim(filename_ncOut_casa))


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

   status = nf_def_var(ncid, 'cresp', NF_REAL, 3, dims, varid_cresp)
   if (status /= nf_noerr) call handle_err(status, "cresp")

   status = nf_def_var(ncid, 'cnpp', NF_REAL, 3, dims, varid_cnpp)
   if (status /= nf_noerr) call handle_err(status, "cnpp")

   status = nf_def_var(ncid, 'cgpp', NF_REAL, 3, dims, varid_cgpp)
   if (status /= nf_noerr) call handle_err(status, "cgpp")

   status = nf_def_var(ncid, 'cleaf', NF_REAL, 3, dims, varid_cleaf)
   if (status /= nf_noerr) call handle_err(status, "cleaf")

   status = nf_def_var(ncid, 'nleaf', NF_REAL, 3, dims, varid_nleaf)
   if (status /= nf_noerr) call handle_err(status, "nleaf")


   status = nf_def_var(ncid, 'cwood', NF_REAL, 3, dims, varid_cwood)
   if (status /= nf_noerr) call handle_err(status, "cwood")

   status = nf_def_var(ncid, 'nwood', NF_REAL, 3, dims, varid_nwood)
   if (status /= nf_noerr) call handle_err(status, "nwood")


   status = nf_def_var(ncid, 'cfroot', NF_REAL, 3, dims, varid_cfroot)
   if (status /= nf_noerr) call handle_err(status, "cfroot")

   status = nf_def_var(ncid, 'nfroot', NF_REAL, 3, dims, varid_nfroot)
   if (status /= nf_noerr) call handle_err(status, "nfroot")


   status = nf_def_var(ncid, 'clitmetb', NF_REAL, 3, dims, varid_clitmetb)
   if (status /= nf_noerr) call handle_err(status, "clitmetb")

   status = nf_def_var(ncid, 'nlitmetb', NF_REAL, 3, dims, varid_nlitmetb)
   if (status /= nf_noerr) call handle_err(status, "nlitmetb")


   status = nf_def_var(ncid, 'clitstr', NF_REAL, 3, dims, varid_clitstr)
   if (status /= nf_noerr) call handle_err(status, "clitstr")

   status = nf_def_var(ncid, 'nlitstr', NF_REAL, 3, dims, varid_nlitstr)
   if (status /= nf_noerr) call handle_err(status, "nlitstr")


   status = nf_def_var(ncid, 'clitcwd', NF_REAL, 3, dims, varid_clitcwd)
   if (status /= nf_noerr) call handle_err(status, "clitcwd")

   status = nf_def_var(ncid, 'nlitcwd', NF_REAL, 3, dims, varid_nlitcwd)
   if (status /= nf_noerr) call handle_err(status, "nlitcwd")


   status = nf_def_var(ncid, 'csoilmic', NF_REAL, 3, dims, varid_csoilmic)
   if (status /= nf_noerr) call handle_err(status, "csoilmic")

   status = nf_def_var(ncid, 'nsoilmic', NF_REAL, 3, dims, varid_nsoilmic)
   if (status /= nf_noerr) call handle_err(status, "nsoilmic")


   status = nf_def_var(ncid, 'csoilslow', NF_REAL, 3, dims, varid_csoilslow)
   if (status /= nf_noerr) call handle_err(status, "csoilslow")

   status = nf_def_var(ncid, 'nsoilslow', NF_REAL, 3, dims, varid_nsoilslow)
   if (status /= nf_noerr) call handle_err(status, "nsoilslow")


   status = nf_def_var(ncid, 'csoilpass', NF_REAL, 3, dims, varid_csoilpass)
   if (status /= nf_noerr) call handle_err(status, "csoilpass")

   status = nf_def_var(ncid, 'nsoilpass', NF_REAL, 3, dims, varid_nsoilpass)
   if (status /= nf_noerr) call handle_err(status, "nsoilpass")

   status = nf_def_var(ncid, 'tairC', NF_REAL, 3, dims, varid_tairC)
   if (status /= nf_noerr) call handle_err(status, "tairC")

   status = nf_def_var(ncid, 'tsoilC', NF_REAL, 3, dims, varid_tsoilC)
   if (status /= nf_noerr) call handle_err(status, "tsoilC")

   status = nf_def_var(ncid, 'cLitInptMet', NF_REAL, 3, dims, varid_clitInptMet)
   if (status /= nf_noerr) call handle_err(status, "cLitInptMet")

   status = nf_def_var(ncid, 'cLitInptStruc', NF_REAL, 3, dims, varid_clitInptStruc)
   if (status /= nf_noerr) call handle_err(status, "cLitInptStruc")

   status = nf_def_var(ncid, 'cpassInpt', NF_REAL, 3, dims, varid_cpassInpt)
   if (status /= nf_noerr) call handle_err(status, "cpassInpt")

   status = nf_def_var(ncid, 'thetaLiq', NF_REAL, 3, dims, varid_thetaLiq)
   if (status /= nf_noerr) call handle_err(status, "thetaLiq")

   status = nf_def_var(ncid, 'fT', NF_REAL, 3, dims, varid_fT)
   if (status /= nf_noerr) call handle_err(status, "f(T)")

   status = nf_def_var(ncid, 'fW', NF_REAL, 3, dims, varid_fW)
   if (status /= nf_noerr) call handle_err(status, "f(W)")

   !! Added N fluxes and pools -mdh 11/25/2019

   status = nf_def_var(ncid, 'nLitInptMet', NF_REAL, 3, dims, varid_nlitInptMet)
   if (status /= nf_noerr) call handle_err(status, "nLitInptMet")

   status = nf_def_var(ncid, 'nLitInptStruc', NF_REAL, 3, dims, varid_nlitInptStruc)
   if (status /= nf_noerr) call handle_err(status, "nLitInptStruc")

   status = nf_def_var(ncid, 'nMinDep', NF_REAL, 3, dims, varid_nmindep)
   if (status /= nf_noerr) call handle_err(status, "nMinDep")

   status = nf_def_var(ncid, 'nMinFix', NF_REAL, 3, dims, varid_nminfix)
   if (status /= nf_noerr) call handle_err(status, "nMinFix")

   status = nf_def_var(ncid, 'nMinUptake', NF_REAL, 3, dims, varid_nminuptake)
   if (status /= nf_noerr) call handle_err(status, "nMinUptake")

   status = nf_def_var(ncid, 'nMinLeach', NF_REAL, 3, dims, varid_nminleach)
   if (status /= nf_noerr) call handle_err(status, "nMinLeach")

   status = nf_def_var(ncid, 'nMinLoss', NF_REAL, 3, dims, varid_nminloss)
   if (status /= nf_noerr) call handle_err(status, "nMinLoss")

   status = nf_def_var(ncid, 'nLitMineralization', NF_REAL, 3, dims, varid_nlittermin)
   if (status /= nf_noerr) call handle_err(status, "nLitMineralization")

   status = nf_def_var(ncid, 'nSoilMineralization', NF_REAL, 3, dims, varid_nsmin)
   if (status /= nf_noerr) call handle_err(status, "nSoilMineralization")

   status = nf_def_var(ncid, 'nSoilImmob', NF_REAL, 3, dims, varid_nsimm)
   if (status /= nf_noerr) call handle_err(status, "nSoilImmob")

   status = nf_def_var(ncid, 'nNetMineralization', NF_REAL, 3, dims, varid_nsnet)
   if (status /= nf_noerr) call handle_err(status, "nNetMineralization")

   status = nf_def_var(ncid, 'nMineral', NF_REAL, 3, dims, varid_nminrl)
   if (status /= nf_noerr) call handle_err(status, "nMineral")


   ! Global attributes
   attr_name = 'CASACNP model output'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'title', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "title")
 
   attr_name = 'NOTE: None of the variables are weighted by land fraction!'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'comment', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "comment")
 
   call get_time_and_date(date_string, time_string)
   attr_name = 'created on ' // date_string // ' ' // time_string
   status = nf_put_att_text(ncid, NF_GLOBAL, 'history', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "history")
 
   attr_name = 'CASACNP Model'
   status = nf_put_att_text(ncid, NF_GLOBAL, 'source', len(trim(attr_name)), trim(attr_name))
   if (status /= nf_noerr) call handle_err(status, "source")

   attr_name = trim(filename_cnpbiome)
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
   !status = nf_put_att_real(ncid, varid_time, '_FillValue', NF_REAL, 1, missing_value)
   !if (status /= nf_noerr) call handle_err(status, "_FillValue")
   !status = nf_put_att_real(ncid, varid_time, 'missing_value', NF_REAL, 1, missing_value)
   !if (status /= nf_noerr) call handle_err(status, "missing_value")

!  ! Attributes of year variable
!  attr_name = 'calendar year'
!  status = nf_put_att_text(ncid, varid_year, 'long_name', len(trim(attr_name)), trim(attr_name))
!  if (status /= nf_noerr) call handle_err(status, "long_name")
!  attr_units = 'year'
!  status = nf_put_att_text(ncid, varid_year, 'units', len(trim(attr_units)), trim(attr_units))
!  if (status /= nf_noerr) call handle_err(status, "units")
!  !status = nf_put_att_int(ncid, varid_year, '_FillValue', NF_INT, 1, MISSING_INT)
!  !if (status /= nf_noerr) call handle_err(status, "_FillValue")
!  !status = nf_put_att_int(ncid, varid_year, 'missing_value', NF_INT, 1, MISSING_INT)
!  !if (status /= nf_noerr) call handle_err(status, "missing_value")

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

   ! Attributes of cnpp variable
   attr_name = 'net primary production'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cnpp, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cgpp variable
   attr_name = 'gross primary production'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cgpp, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cresp variable
   attr_name = 'soil heterotrophic respiration'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cresp, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cleaf variable
   attr_name = 'leaf carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cleaf, attr_name, attr_units, MISSING_VALUE)
   
   ! Attributes of nleaf variable
   attr_name = 'leaf nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nleaf, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cwood variable
   attr_name = 'wood carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cwood, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nwood variable
   attr_name = 'wood nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nwood, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cfroot variable
   attr_name = 'fine root carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_cfroot, attr_name, attr_units, MISSING_VALUE)


   ! Attributes of nfroot variable
   attr_name = 'fine root nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nfroot, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of clitmetb variable
   attr_name = 'metabolic litter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_clitmetb, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nlitmetb variable
   attr_name = 'metabolic litter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nlitmetb, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of clitstr variable
   attr_name = 'structural litter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_clitstr, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nlitstr variable
   attr_name = 'structural litter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nlitstr, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of clitcwd variable
   attr_name = 'coarse woody debris carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_clitcwd, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nlitcwd variable
   attr_name = 'coarse woody debris nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nlitcwd, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of csoilmic variable
   attr_name = 'microbial soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_csoilmic, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nsoilmic variable
   attr_name = 'microbial soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nsoilmic, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of csoilslow variable
   attr_name = 'slow soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_csoilslow, attr_name, attr_units, MISSING_VALUE)
 
  ! Attributes of nsoilslow variable
   attr_name = 'slow soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nsoilslow, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of csoilpass variable
   attr_name = 'passive soil organic matter carbon'
   attr_units = 'gC m-2'
   call PutVariableAttributeReal(ncid, varid_csoilpass, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nsoilpass variable
   attr_name = 'passive soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nsoilpass, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of tairC variable
   attr_name = 'mean annual air temperature'
   attr_units = 'degrees C'
   call PutVariableAttributeReal(ncid, varid_tairC, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of tsoilC variable
   attr_name = 'mean annual soil temperature in top 50 cm'
   attr_units = 'degrees C'
   call PutVariableAttributeReal(ncid, varid_tsoilC, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLitInptMet variable
   attr_name = 'Metabolic Litter C Inputs'
   attr_units = 'gC m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_clitInptMet, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLitInptStruc variable
   attr_name = 'Structural Litter C Inputs'
   attr_units = 'gC m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_clitInptStruc, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cpassInpt variable
   attr_name = 'Passive Pool C Inputs'
   attr_units = 'gC m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_cpassInpt, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of thetaLiq variable
   attr_name = 'Fraction of liquid soil water saturation (0.0-1.0)'
   attr_units = 'fraction'
   call PutVariableAttributeReal(ncid, varid_thetaLiq, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of f(T) variable
   attr_name = 'soil temperature multiplier on decomposition rate (0.0-1.0)'
   attr_units = 'multiplier'
   call PutVariableAttributeReal(ncid, varid_fT, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of f(W) variable
   attr_name = 'soil moisture multiplier on decomposition rate (0.0-1.0)'
   attr_units = 'multiplier'
   call PutVariableAttributeReal(ncid, varid_fW, attr_name, attr_units, MISSING_VALUE)

   !! Added annual N flux and pools -mdh 11/25/2019

   ! Attributes of nLitInptMet variable
   attr_name = 'Metabolic Litter N Inputs'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nlitInptMet, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nLitInptStruc variable
   attr_name = 'Structural Litter N Inputs'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nlitInptStruc, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMinDep variable
   attr_name = 'N deposition'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nmindep, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMinFix variable
   attr_name = 'N Fixation'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nminfix, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMinUptake variable
   attr_name = 'Plant Mineral N Uptake'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nminuptake, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMinLeach variable
   attr_name = 'Mineral N Leaching'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nminleach, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMinLoss variable
   attr_name = 'Mineral N Loss (N2O)'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nminloss, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nLitMineralization variable
   attr_name = 'Litter N Mineralization'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nlittermin, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nSoilMineralization variable
   attr_name = 'Soil N Mineralization'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nsmin, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nSoilImmob variable
   attr_name = 'Soil N Immobilization'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nsimm, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nNetMineralization variable
   attr_name = 'Net N mineralization'
   attr_units = 'gN m-2 yr-1'
   call PutVariableAttributeReal(ncid, varid_nsnet, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMineral variable
   attr_name = 'Mineral Soil N'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nminrl, attr_name, attr_units, MISSING_VALUE)


   ! --------------- End the definition phase so that variables can be written to the file ---------------
   status = nf_enddef(ncid)
   if (status /= nf_noerr) call handle_err(status, "enddef")


   ! ------------------------------  Write variable values to filename_ncOut_casa ------------------------------

 
   !! Write time using nf_put_vara_real instead of nf_put_var so
   !! that time gets written correctly. -mdh 12/28/2016

   time = real(year)
   start1 = (/ 1 /)
   count1 = (/ ntimes /)    ! ntimes = 1 in this subroutine

!! status =  nf_put_var(ncid, varid_time, time)
   status =  nf_put_vara_real(ncid, varid_time, start1, count1, time)
   if (status /= nf_noerr) call handle_err(status, "put_vara_real(time)")

!! status =  nf_put_var(ncid, varid_year, year)
!  status =  nf_put_vara_int(ncid, varid_year, start1, count1, year)
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


!  casabal%FCnppyear(npt)     - NPP (gC/m2/yr) 
!  casabal%FCgppyear(npt)     - GPP (gC/m2/yr) 
!  casabal%FCrsyear(npt)    - soil respiration (gC/m2/yr)?        
!  casabal%FCrpyear(npt)    - plant respiration (gC/m2/yr)?        
!  casapoolAn%cplantAn(npt,LEAF)  - average annual plant leaf C (gC/m2)
!  casapoolAn%cplantAn(npt,WOOD)  - average annual plant wood C (gC/m2)
!  casapoolAn%cplantAn(npt,FROOT) - average annual plant fine root C (gC/m2)

   IGBP_PFT(:,:) = MISSING_INT
   landarea(:,:) = MISSING_VALUE
   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      !write(*,*) 'WritePoolFluxNcFile_casacnp_annual: npt = ', npt, 'ilat =', ilat, 'ilon =', ilon
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
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
      var1(ilon,ilat,itime) = casabal%FCnppyear(npt)
      var2(ilon,ilat,itime) = casabal%FCrsyear(npt)
      var3(ilon,ilat,itime) = casapoolAn%cplantAn(npt,LEAF)
      var4(ilon,ilat,itime) = casapoolAn%cplantAn(npt,WOOD)
      var5(ilon,ilat,itime) = casapoolAn%cplantAn(npt,FROOT)
      var6(ilon,ilat,itime) = casabal%FCgppyear(npt)
   enddo


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
!! status =  nf_put_var(ncid, varid_cnpp, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cnpp, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(cnpp)")
 
   status =  nf_put_var(ncid, varid_cgpp, var6, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cgpp)")

   status =  nf_put_var(ncid, varid_cresp, var2, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cresp)")
 
   status =  nf_put_var(ncid, varid_cleaf, var3, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cleaf)")

   status =  nf_put_var(ncid, varid_cwood, var4, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cwood)")

   status =  nf_put_var(ncid, varid_cfroot, var5, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cfroot)")


!  casapoolAn%clitterAn(npt,METB)    - average annual metabolic litter C (gC/m2)
!  casapoolAn%clitterAn(npt,STR)    - average annual structural litter C (gC/m2)
!  casapoolAn%clitterAn(npt,CWD)    - average annual coarse woody debris C (gC/m2)
!  casapoolAn%csoilAn(npt,MIC)        - average annual microbial soil  C (gC/m2)
!  casapoolAn%csoilAn(npt,SLOW)        - average annual slow soil C (gC/m2)
!  casapoolAn%csoilAn(npt,PASS)        - average annual passive soil C (gC/m2)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapoolAn%clitterAn(npt,METB)
      var2(ilon,ilat,itime) = casapoolAn%clitterAn(npt,STR)
      var3(ilon,ilat,itime) = casapoolAn%clitterAn(npt,CWD)
      var4(ilon,ilat,itime) = casapoolAn%csoilAn(npt,MIC)
      var5(ilon,ilat,itime) = casapoolAn%csoilAn(npt,SLOW)
      var6(ilon,ilat,itime) = casapoolAn%csoilAn(npt,PASS)
   enddo

   status =  nf_put_var(ncid, varid_clitmetb, var1, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(clitmetb)")
 
   status =  nf_put_var(ncid, varid_clitstr, var2, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(clitstr)")
 
   status =  nf_put_var(ncid, varid_clitcwd, var3, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(clitcwd)")

   status =  nf_put_var(ncid, varid_csoilmic, var4, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(csoilmic)")

   status =  nf_put_var(ncid, varid_csoilslow, var5, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(csoilslow)")

   status =  nf_put_var(ncid, varid_csoilpass, var6, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(csoilpass)")


!  casapoolAn%nplantAn(npt,LEAF)    - average annual plant leaf N (gN/m2)
!  casapoolAn%nplantAn(npt,WOOD)    - average annual plant wood N (gN/m2)
!  casapoolAn%nplantAn(npt,FROOT)    - average annual plant fine root N (gN/m2)
!  casapoolAn%nlitterAn(npt,METB)    - average annual metabolic litter N (gN/m2)
!  casapoolAn%nlitterAn(npt,STR)    - average annual structural litter N (gN/m2)
!  casapoolAn%nlitterAn(npt,CWD)    - average annual coarse woody debris N (gN/m2)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapoolAn%nplantAn(npt,LEAF)
      var2(ilon,ilat,itime) = casapoolAn%nplantAn(npt,WOOD)
      var3(ilon,ilat,itime) = casapoolAn%nplantAn(npt,FROOT)
      var4(ilon,ilat,itime) = casapoolAn%nlitterAn(npt,METB)
      var5(ilon,ilat,itime) = casapoolAn%nlitterAn(npt,STR)
      var6(ilon,ilat,itime) = casapoolAn%nlitterAn(npt,CWD)
   enddo

   status =  nf_put_var(ncid, varid_nleaf, var1, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nleaf)")
 
   status =  nf_put_var(ncid, varid_nwood, var2, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nwood)")
 
   status =  nf_put_var(ncid, varid_nfroot, var3, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nfroot)")

   status =  nf_put_var(ncid, varid_nlitmetb, var4, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nlitmetb)")

   status =  nf_put_var(ncid, varid_nlitstr, var5, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nlitstr)")

   status =  nf_put_var(ncid, varid_nlitcwd, var6, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nlitcwd)")


!  casapoolAn%nsoilAn(npt,MIC)    - average annual microbial soil N (gN/m2)
!  casapoolAn%nsoilAn(npt,SLOW)   - average annual slow soil N (gN/m2)
!  casapoolAn%nsoilAn(npt,PASS)   - average annual passive soil N (gN/m2)
!  casapoolAn%tairAn(npt)         - average annual air temperature (C)
!  casapoolAn%tsoilAn(npt)        - average annual soil temperature (C)
!  casaflux%ClitInptMetAn(npt)    - metabolic litter inputs (gC/m2/yr)
!  casaflux%ClitInptStrucAn(npt)  - structural litter inputs (gC/m2/yr)
!  casaflux%CpassInptAn(npt)      - inputs to passive soil pool (gC/m2/yr)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE
   var8(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapoolAn%nsoilAn(npt,MIC)
      var2(ilon,ilat,itime) = casapoolAn%nsoilAn(npt,SLOW)
      var3(ilon,ilat,itime) = casapoolAn%nsoilAn(npt,PASS)
      var4(ilon,ilat,itime) = casapoolAn%tairAn(npt)
      var5(ilon,ilat,itime) = casapoolAn%tsoilAn(npt)
      var6(ilon,ilat,itime) = casaflux%ClitInptMetAn(npt)
      var7(ilon,ilat,itime) = casaflux%ClitInptStrucAn(npt)
      var8(ilon,ilat,itime) = casaflux%CpassInptAn(npt)
      !write(*,*) 'var7(',ilon,ilat,itime,')=', var7(ilon,ilat,itime)
   enddo

   status =  nf_put_var(ncid, varid_nsoilmic, var1, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nsoilmic)")
 
   status =  nf_put_var(ncid, varid_nsoilslow, var2, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nsoilslow)")
 
   status =  nf_put_var(ncid, varid_nsoilpass, var3, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nsoilpass)")

   status =  nf_put_var(ncid, varid_tairC, var4, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(tairC)")

   status =  nf_put_var(ncid, varid_tsoilC, var5, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(tsoilC)")

   status =  nf_put_var(ncid, varid_clitInptMet, var6, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cLitInptMet)")

   status =  nf_put_var(ncid, varid_clitInptStruc, var7, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cLitInptStruc)")

   status =  nf_put_var(ncid, varid_cpassInpt, var8, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cpassInpt)")


!  f(T) and f(W) output added 11/20/2017. ThetaLiq added 11/27/2017.
!  casaflux%fTAn(npt)  - mean annual soil temperature multiplier on decomposition rate (0.0-1.0+)
!  casaflux%fWAn(npt)  - mean annual soil moisure multiplier on decomposition rate (0.0-1.0)
!  casaflux%thetaLiqAn(npt)  - mean annual fraction of liquid water saturation (0.0-1.0)
!  casapoolAn%NsoilminAn(npt) - mean annual mineral soil N (gN/m2)
!  casaflux%NlitInptMetAn(npt) - total annual metabolic litter N inputs (gN/m2/yr)
!  casaflux%NlitInptStrucAn(npt) - total annual structural litter N inputs (gN/m2/yr)
!  casaflux%NmindepAn(npt) - total annual N deposition (gN/m2/yr)
!  casaflux%NminfixAn(npt) - total annual N fixation (gN/m2/yr)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE
   var8(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapoolAn%fTAn(npt)
      var2(ilon,ilat,itime) = casapoolAn%fWAn(npt)
      var3(ilon,ilat,itime) = casapoolAn%thetaLiqAn(npt)

      var4(ilon,ilat,itime) = casapoolAn%NsoilminAn(npt)
      var5(ilon,ilat,itime) = casaflux%NlitInptMetAn(npt)
      var6(ilon,ilat,itime) = casaflux%NlitInptStrucAn(npt)
      var7(ilon,ilat,itime) = casaflux%NmindepAn(npt)
      var8(ilon,ilat,itime) = casaflux%NminfixAn(npt)
   enddo

!! status =  nf_put_var(ncid, varid_fT, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_fT, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(fT)")
 
   status =  nf_put_vara_real(ncid, varid_fW, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(fW)")

   status =  nf_put_vara_real(ncid, varid_thetaLiq, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(thetaLiq)")

   status =  nf_put_vara_real(ncid, varid_nminrl, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMineral)")

   status =  nf_put_vara_real(ncid, varid_nlitinptmet, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(NlitInptMet)")

   status =  nf_put_vara_real(ncid, varid_nlitinptstruc, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(NlitInptStruc)")

   status =  nf_put_vara_real(ncid, varid_nmindep, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinDep)")

   status =  nf_put_vara_real(ncid, varid_nminfix, start3, count3, var8)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinFix)")


!  casaflux%NminuptakeAn(npt) - total annual Plant N uptake (gN/m2/yr)
!  casaflux%NminleachAn(npt) - total annual Mineral N Leaching (gN/m2/yr)
!  casaflux%NminlossAn(npt) - total annual Mineral N Loss (gN/m2/yr)
!  casaflux%NlitterminAn(npt) - total annual Litter N mineralization (gN/m2/yr)
!  casaflux%NsminAn(npt) - total annual Soil N mineralization (gN/m2/yr)
!  casaflux%NsimmAn(npt) - total annual Soil N Immobilization (gN/m2/yr)
!  casaflux%NsnetAn(npt) - total annual Net N mineralization (gN/m2/yr)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casaflux%NminuptakeAn(npt)
      var2(ilon,ilat,itime) = casaflux%NminleachAn(npt)
      var3(ilon,ilat,itime) = casaflux%NminlossAn(npt)
      var4(ilon,ilat,itime) = casaflux%NlitterminAn(npt)
      var5(ilon,ilat,itime) = casaflux%NsminAn(npt)
      var6(ilon,ilat,itime) = casaflux%NsimmAn(npt)
      var7(ilon,ilat,itime) = casaflux%NsnetAn(npt)
   enddo

   status =  nf_put_vara_real(ncid, varid_nminuptake, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinUptake)")

   status =  nf_put_vara_real(ncid, varid_nminleach, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinLeach)")

   status =  nf_put_vara_real(ncid, varid_nminloss, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinLoss)")

   status =  nf_put_vara_real(ncid, varid_nlittermin, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLitMineralization)")

   status =  nf_put_vara_real(ncid, varid_nsmin, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSoilMineralization)")

   status =  nf_put_vara_real(ncid, varid_nsimm, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSoilImmob)")

   status =  nf_put_vara_real(ncid, varid_nsnet, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "put_var(nNetMineralization)")


   deallocate(var1,var2,var3,var4,var5,var6,var7,var8,IGBP_PFT,landarea)

   status = nf_close(ncid)

   if (verbose .ge. 0) print *, "Done writing output to ", trim(filename_ncOut_casa), "..."

END SUBROUTINE WritePoolFluxNcFile_casacnp_annual

!----------------------------------------------------------------------------------------------------
    SUBROUTINE WritePoolFluxNcFile_casacnp_daily(filename_ncOut_casa, mp, year, iday)
!
!   DESCRIPTION
!   Write daily output to netcdf file filename_ncOut_casa.  
!   If iday=1, then create and define the file also.
!   The file contains 365 days of output for one year.
!   Called by casacnpdriver.
!
!   Melannie Hartman. June 6, 2016

    USE casadimension
    USE clmgridvariable
    USE casavariable
    USE casaparm
    USE define_types
    implicit none
    include 'netcdf.inc'
    real(4), parameter :: MISSING_VALUE = 1.e+36
    integer, parameter :: MISSING_INT = -9999

!   ARGUMENTS
      character(len=*), intent(in) :: filename_ncOut_casa    ! NetCDF output file name 
      integer, intent(in) :: mp                    ! number of points with valid data
      integer, intent(in) :: year            ! output year
      integer, intent(in) :: iday            ! output day
!   LOCAL VARIABLES
      integer :: i
      integer :: ncid               ! netcdf file ID
      integer :: status             ! function return status
!     integer :: dimid_mp           ! netcdf dimension id
      integer :: dimid_lat          ! netcdf dimension id
      integer :: dimid_lon          ! netcdf dimension id
      integer :: dimid_time         ! netcdf dimension id
      integer :: nlon, nlat, ntimes     ! Dimension sizes for NetCDf file
      integer :: nwrtimes               ! Number of times that will be written when this subroutine is called
      integer :: dims(3)                ! Array of NetCDF dimension IDs for defining variables
      integer :: start1(1), count1(1)   ! start and count arrays for writing 1-D data from netcdf files
      !!integer :: start2(2), count2(2) ! start and count arrays for writing 2-D data from netcdf files
      integer :: start3(3), count3(3)   ! start and count arrays for writing 3-D data from netcdf files
      integer :: varid_lon, varid_lat   ! NetCDF variable ID for latitude and longitude
      integer :: varid_time           ! NetCDF variable ID for time
      integer :: varid_day            ! NetCDF variable ID for day
      integer :: varid_mask           ! NetCDF variable ID for cellMissing(nlon,nlat)
      integer :: varid_cellid         ! NetCDF variable ID for cellid(nlon,nlat)
      integer :: varid_igbp           ! NetCDF variable ID for IGBP_PFT(nlon,nlat)
      integer :: varid_landarea       ! NetCDF variable ID for land area
      integer :: varid_cleaf          ! NetCDF variable ID for leaf C
      integer :: varid_cwood          ! NetCDF variable ID for wood C
      integer :: varid_cfroot         ! NetCDF variable ID for fine root C
      integer :: varid_nleaf          ! NetCDF variable ID for leaf N
      integer :: varid_nwood          ! NetCDF variable ID for wood N
      integer :: varid_nfroot         ! NetCDF variable ID for fine root N
      integer :: varid_clitmetb       ! NetCDF variable ID for metablic litter C
      integer :: varid_clitstr        ! NetCDF variable ID for structural litter C
      integer :: varid_clitcwd        ! NetCDF variable ID for coarse weeody debris C
      integer :: varid_nlitmetb       ! NetCDF variable ID for metablic litter N
      integer :: varid_nlitstr        ! NetCDF variable ID for structural litter N
      integer :: varid_nlitcwd        ! NetCDF variable ID for coarse weeody debris N
      integer :: varid_csoilmic       ! NetCDF variable ID for microbial soil C
      integer :: varid_csoilslow      ! NetCDF variable ID for slow soil C
      integer :: varid_csoilpass      ! NetCDF variable ID for passive soil C
      integer :: varid_nsoilmic       ! NetCDF variable ID for microbial soil N
      integer :: varid_nsoilslow      ! NetCDF variable ID for slow soil N
      integer :: varid_nsoilpass      ! NetCDF variable ID for passive soil N
      integer :: varid_cnpp           ! NetCDF variable ID for NPP C
      integer :: varid_cgpp           ! NetCDF variable ID for GPP C
      integer :: varid_cresp          ! NetCDF variable ID for soil respiration C
      integer :: varid_tairC          ! NetCDF variable ID for air temperature (C)
      integer :: varid_tsoilC         ! NetCDF variable ID for soil temperature (C)
      integer :: varid_clitInptMet    ! NetCDF variable ID for metabolic litter inputs
      integer :: varid_clitInptStruc  ! NetCDF variable ID for structural litter inputs
      integer :: varid_cpassInpt      ! NetCDF variable ID for passive pool C inputs
      integer :: varid_fT, varid_fW   ! NetCDF variable ID for soil temperature and moisture decomposition rate modifiers
      integer :: varid_thetaLiq       ! NetCDF variable ID for fraction of liquid soil water saturation
      !! Added N fluxes and pools -mdh 11/25/2019
      integer :: varid_nlitInptMet     ! NetCDF variable ID for metabolic litter N inputs
      integer :: varid_nlitInptStruc   ! NetCDF variable ID for structural litter N inputs
      integer :: varid_nmindep         ! NetCDF variable ID for N deposition
      integer :: varid_nminfix         ! NetCDF variable ID for N fixation
      integer :: varid_nminuptake      ! NetCDF variable ID for plant N uptake
      integer :: varid_nminleach       ! NetCDF variable ID for mineral N leaching
      integer :: varid_nminloss        ! NetCDF variable ID for mineral N loss (N2O)
      integer :: varid_nlittermin      ! NetCDF variable ID for litter N mineralization
      integer :: varid_nsmin           ! NetCDF variable ID for soil N mineralization
      integer :: varid_nsimm           ! NetCDF variable ID for soil N immobilization
      integer :: varid_nsnet           ! NetCDF variable ID for net N mineralization
      integer :: varid_nminrl          ! NetCDF variable ID for mineral soil N
      character*100 :: attr_name       ! String for assigning global and variable attributes
      character*100 :: attr_units      ! String for assigning global and variable attributes
      character*10 :: date_string      ! String for assigning date to global attributes
      character*8 :: time_string       ! String for assigning time to global attributes
      integer :: verbose=0
      integer :: npt, ilon, ilat, itime
      integer, allocatable :: IGBP_PFT(:,:)  ! IGBP_PFT(nlon,nlat) IGBP PFT classification (1-18)
      real(4), allocatable :: time(:)        ! time array (1..ntimes)
      integer, allocatable :: days(:)        ! day array (1..ntimes)
      real(4), allocatable :: landarea(:,:)  ! landarea(nlon,nlat) km^2
      real(4), allocatable :: var1(:,:,:)    ! gridded output variable
      real(4), allocatable :: var2(:,:,:)    ! gridded output variable
      real(4), allocatable :: var3(:,:,:)    ! gridded output variable
      real(4), allocatable :: var4(:,:,:)    ! gridded output variable
      real(4), allocatable :: var5(:,:,:)    ! gridded output variable
      real(4), allocatable :: var6(:,:,:)    ! gridded output variable
      real(4), allocatable :: var7(:,:,:)    ! gridded output variable
      real(4), allocatable :: var8(:,:,:)    ! gridded output variable
    

      if (verbose .ge. 1) print *, "Writing output to file ", trim(filename_ncOut_casa), " on day", iday, "..."

!Output these variables in NetCDF file:
!       veg%iveg(npt)           - IGBP PFTs
!       casamet%lat(npt)        - latitudes    
!       casamet%lon(npt)        - longitudes
!       casamet%areacell(npt)*(1.0e-6)     - landarea (km^2)
!       casaflux%Cnpp(npt)      - NPP (gC/m2/day) 
!       casaflux%Cgpp(npt)      - GPP (gC/m2/day) 
!       casaflux%Crsoil(npt)    - soil respiration (gC/m2/day)         
!       casaflux%Crp(npt)       - plant respiration (gC/m2/day) 
!       casapool%cplant(npt,LEAF)    - plant leaf C (gC/m2)
!       casapool%cplant(npt,WOOD)    - plant wood C (gC/m2)
!       casapool%cplant(npt,FROOT)   - plant fine root C (gC/m2)
!       casapool%clitter(npt,METB)   - metabolic litter C (gC/m2)
!       casapool%clitter(npt,STR)    - structural litter C (gC/m2)
!       casapool%clitter(npt,CWD)    - coarse woody debris C (gC/m2)
!       casapool%csoil(npt,MIC)      - microbial soil  C (gC/m2)
!       casapool%csoil(npt,SLOW)     - slow soil C (gC/m2)
!       casapool%csoil(npt,PASS)     - passive soil C (gC/m2)
!       casapool%nplant(npt,LEAF)    - plant leaf N (gN/m2)
!       casapool%nplant(npt,WOOD)    - plant wood N (gN/m2)
!       casapool%nplant(npt,FROOT)   - plant fine root N (gN/m2)
!       casapool%nlitter(npt,METB)   - metabolic litter N (gN/m2)
!       casapool%nlitter(npt,STR)    - structural litter N (gN/m2)
!       casapool%nlitter(npt,CWD)    - coarse woody debris N (gN/m2)
!       casapool%nsoil(npt,MIC)      - microbial soil N (gN/m2)
!       casapool%nsoil(npt,SLOW)     - slow soil N (gN/m2)
!       casapool%nsoil(npt,PASS)     - passive soil N (gN/m2)
!       casaflux%fT(npt)             - temperature effect on soil decomposition rate 
!       casaflux%fW(npt)             - moisture effect on soil decomposition rate 
!       casaflux%thetaLiq(npt)       - fraction of liquid soil water saturation (0.0-1.0)

   dims(1) = 0
   dims(2) = 0
   ntimes = 365
   nlat = clmgrid%nlat
   nlon = clmgrid%nlon

   if (iday == 1) then 

      ! Create the netcdf file 
   
      status = nf_create(filename_ncOut_casa, NF_CLOBBER, ncid)
      if (status /= nf_noerr) call handle_err(status, trim(filename_ncOut_casa))
   
   
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
   
      status = nf_def_var(ncid, 'cresp', NF_REAL, 3, dims, varid_cresp)
      if (status /= nf_noerr) call handle_err(status, "cresp")
   
      status = nf_def_var(ncid, 'cnpp', NF_REAL, 3, dims, varid_cnpp)
      if (status /= nf_noerr) call handle_err(status, "cnpp")
   
      status = nf_def_var(ncid, 'cgpp', NF_REAL, 3, dims, varid_cgpp)
      if (status /= nf_noerr) call handle_err(status, "cgpp")
   
      status = nf_def_var(ncid, 'cleaf', NF_REAL, 3, dims, varid_cleaf)
      if (status /= nf_noerr) call handle_err(status, "cleaf")
   
      status = nf_def_var(ncid, 'nleaf', NF_REAL, 3, dims, varid_nleaf)
      if (status /= nf_noerr) call handle_err(status, "nleaf")
   
   
      status = nf_def_var(ncid, 'cwood', NF_REAL, 3, dims, varid_cwood)
      if (status /= nf_noerr) call handle_err(status, "cwood")
   
      status = nf_def_var(ncid, 'nwood', NF_REAL, 3, dims, varid_nwood)
      if (status /= nf_noerr) call handle_err(status, "nwood")
   
   
      status = nf_def_var(ncid, 'cfroot', NF_REAL, 3, dims, varid_cfroot)
      if (status /= nf_noerr) call handle_err(status, "cfroot")
   
      status = nf_def_var(ncid, 'nfroot', NF_REAL, 3, dims, varid_nfroot)
      if (status /= nf_noerr) call handle_err(status, "nfroot")
   
   
      status = nf_def_var(ncid, 'clitmetb', NF_REAL, 3, dims, varid_clitmetb)
      if (status /= nf_noerr) call handle_err(status, "clitmetb")
   
      status = nf_def_var(ncid, 'nlitmetb', NF_REAL, 3, dims, varid_nlitmetb)
      if (status /= nf_noerr) call handle_err(status, "nlitmetb")
   
   
      status = nf_def_var(ncid, 'clitstr', NF_REAL, 3, dims, varid_clitstr)
      if (status /= nf_noerr) call handle_err(status, "clitstr")
   
      status = nf_def_var(ncid, 'nlitstr', NF_REAL, 3, dims, varid_nlitstr)
      if (status /= nf_noerr) call handle_err(status, "nlitstr")
   
   
      status = nf_def_var(ncid, 'clitcwd', NF_REAL, 3, dims, varid_clitcwd)
      if (status /= nf_noerr) call handle_err(status, "clitcwd")
   
      status = nf_def_var(ncid, 'nlitcwd', NF_REAL, 3, dims, varid_nlitcwd)
      if (status /= nf_noerr) call handle_err(status, "nlitcwd")
   
   
      status = nf_def_var(ncid, 'csoilmic', NF_REAL, 3, dims, varid_csoilmic)
      if (status /= nf_noerr) call handle_err(status, "csoilmic")
   
      status = nf_def_var(ncid, 'nsoilmic', NF_REAL, 3, dims, varid_nsoilmic)
      if (status /= nf_noerr) call handle_err(status, "nsoilmic")
   
   
      status = nf_def_var(ncid, 'csoilslow', NF_REAL, 3, dims, varid_csoilslow)
      if (status /= nf_noerr) call handle_err(status, "csoilslow")
   
      status = nf_def_var(ncid, 'nsoilslow', NF_REAL, 3, dims, varid_nsoilslow)
      if (status /= nf_noerr) call handle_err(status, "nsoilslow")
   
   
      status = nf_def_var(ncid, 'csoilpass', NF_REAL, 3, dims, varid_csoilpass)
      if (status /= nf_noerr) call handle_err(status, "csoilpass")
   
      status = nf_def_var(ncid, 'nsoilpass', NF_REAL, 3, dims, varid_nsoilpass)
      if (status /= nf_noerr) call handle_err(status, "nsoilpass")
   
      status = nf_def_var(ncid, 'tairC', NF_REAL, 3, dims, varid_tairC)
      if (status /= nf_noerr) call handle_err(status, "tairC")
   
      status = nf_def_var(ncid, 'tsoilC', NF_REAL, 3, dims, varid_tsoilC)
      if (status /= nf_noerr) call handle_err(status, "tsoilC")
   
      status = nf_def_var(ncid, 'cLitInptMet', NF_REAL, 3, dims, varid_clitInptMet)
      if (status /= nf_noerr) call handle_err(status, "cLitInptMet")

      status = nf_def_var(ncid, 'cLitInptStruc', NF_REAL, 3, dims, varid_clitInptStruc)
      if (status /= nf_noerr) call handle_err(status, "cLitInptStruc")

      status = nf_def_var(ncid, 'cpassInpt', NF_REAL, 3, dims, varid_cpassInpt)
      if (status /= nf_noerr) call handle_err(status, "cpassInpt")

      status = nf_def_var(ncid, 'thetaLiq', NF_REAL, 3, dims, varid_thetaLiq)
      if (status /= nf_noerr) call handle_err(status, "thetaLiq")

      status = nf_def_var(ncid, 'fT', NF_REAL, 3, dims, varid_fT)
      if (status /= nf_noerr) call handle_err(status, "fT")

      status = nf_def_var(ncid, 'fW', NF_REAL, 3, dims, varid_fW)
      if (status /= nf_noerr) call handle_err(status, "fW")


      !! Added daily N fluxes and pools -mdh 11/25/2019
   
      status = nf_def_var(ncid, 'nLitInptMet', NF_REAL, 3, dims, varid_nlitInptMet)
      if (status /= nf_noerr) call handle_err(status, "nLitInptMet")
   
      status = nf_def_var(ncid, 'nLitInptStruc', NF_REAL, 3, dims, varid_nlitInptStruc)
      if (status /= nf_noerr) call handle_err(status, "nLitInptStruc")
   
      status = nf_def_var(ncid, 'nMinDep', NF_REAL, 3, dims, varid_nmindep)
      if (status /= nf_noerr) call handle_err(status, "nMinDep")
   
      status = nf_def_var(ncid, 'nMinFix', NF_REAL, 3, dims, varid_nminfix)
      if (status /= nf_noerr) call handle_err(status, "nMinFix")
   
      status = nf_def_var(ncid, 'nMinUptake', NF_REAL, 3, dims, varid_nminuptake)
      if (status /= nf_noerr) call handle_err(status, "nMinUptake")
   
      status = nf_def_var(ncid, 'nMinLeach', NF_REAL, 3, dims, varid_nminleach)
      if (status /= nf_noerr) call handle_err(status, "nMinLeach")
   
      status = nf_def_var(ncid, 'nMinLoss', NF_REAL, 3, dims, varid_nminloss)
      if (status /= nf_noerr) call handle_err(status, "nMinLoss")
   
      status = nf_def_var(ncid, 'nLitMineralization', NF_REAL, 3, dims, varid_nlittermin)
      if (status /= nf_noerr) call handle_err(status, "nLitMineralization")
   
      status = nf_def_var(ncid, 'nSoilMineralization', NF_REAL, 3, dims, varid_nsmin)
      if (status /= nf_noerr) call handle_err(status, "nSoilMineralization")
   
      status = nf_def_var(ncid, 'nSoilImmob', NF_REAL, 3, dims, varid_nsimm)
      if (status /= nf_noerr) call handle_err(status, "nSoilImmob")
   
      status = nf_def_var(ncid, 'nNetMineralization', NF_REAL, 3, dims, varid_nsnet)
      if (status /= nf_noerr) call handle_err(status, "nNetMineralization")
   
      status = nf_def_var(ncid, 'nMineral', NF_REAL, 3, dims, varid_nminrl)
      if (status /= nf_noerr) call handle_err(status, "nMineral")
   
      ! Global attributes
      attr_name = 'CASACNP model output'
      status = nf_put_att_text(ncid, NF_GLOBAL, 'title', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "title")
    
      attr_name = 'NOTE: None of the variables are weighted by land fraction!'
      status = nf_put_att_text(ncid, NF_GLOBAL, 'comment', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "comment")
    
      call get_time_and_date(date_string, time_string)
      attr_name = 'created on ' // date_string // ' ' // time_string
      status = nf_put_att_text(ncid, NF_GLOBAL, 'history', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "history")
    
      attr_name = 'CASACNP Model'
      status = nf_put_att_text(ncid, NF_GLOBAL, 'source', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "source")
   
      attr_name = trim(filename_cnpbiome)
      status = nf_put_att_text(ncid, NF_GLOBAL, 'parameters', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "parameters")

      ! ------------------------------ Attributes of the variables ------------------------------
   
      ! Attributes of coordinate time variable
      attr_name = 'coordinate time'
      status = nf_put_att_text(ncid, varid_time, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
      !!attr_units = 'simulation time in years'
      attr_units = '1..ntimes'
      status = nf_put_att_text(ncid, varid_time, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
      !status = nf_put_att_real(ncid, varid_time, '_FillValue', NF_REAL, 1, missing_value)
      !if (status /= nf_noerr) call handle_err(status, "_FillValue")
      !status = nf_put_att_real(ncid, varid_time, 'missing_value', NF_REAL, 1, missing_value)
      !if (status /= nf_noerr) call handle_err(status, "missing_value")

      ! Attributes of day variable
      attr_name = 'day of year'
      status = nf_put_att_text(ncid, varid_day, 'long_name', len(trim(attr_name)), trim(attr_name))
      if (status /= nf_noerr) call handle_err(status, "long_name")
      attr_units = '1..365'
      status = nf_put_att_text(ncid, varid_day, 'units', len(trim(attr_units)), trim(attr_units))
      if (status /= nf_noerr) call handle_err(status, "units")
      !status = nf_put_att_int(ncid, varid_day, '_FillValue', NF_INT, 1, MISSING_INT)
      !if (status /= nf_noerr) call handle_err(status, "_FillValue")
      !status = nf_put_att_int(ncid, varid_day, 'missing_value', NF_INT, 1, MISSING_INT)
      !if (status /= nf_noerr) call handle_err(status, "missing_value")

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
   
      ! Attributes of cnpp variable
      attr_name = 'net primary production'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cnpp, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cgpp variable
      attr_name = 'gross primary production'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cgpp, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cresp variable
      attr_name = 'soil heterotrophic respiration'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cresp, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cleaf variable
      attr_name = 'leaf carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cleaf, attr_name, attr_units, MISSING_VALUE)
      
      ! Attributes of nleaf variable
      attr_name = 'leaf nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nleaf, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cwood variable
      attr_name = 'wood carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cwood, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nwood variable
      attr_name = 'wood nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nwood, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cfroot variable
      attr_name = 'fine root carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_cfroot, attr_name, attr_units, MISSING_VALUE)
   
   
      ! Attributes of nfroot variable
      attr_name = 'fine root nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nfroot, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of clitmetb variable
      attr_name = 'metabolic litter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_clitmetb, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nlitmetb variable
      attr_name = 'metabolic litter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nlitmetb, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of clitstr variable
      attr_name = 'structural litter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_clitstr, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nlitstr variable
      attr_name = 'structural litter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nlitstr, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of clitcwd variable
      attr_name = 'coarse woody debris carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_clitcwd, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nlitcwd variable
      attr_name = 'coarse woody debris nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nlitcwd, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of csoilmic variable
      attr_name = 'microbial soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_csoilmic, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nsoilmic variable
      attr_name = 'microbial soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nsoilmic, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of csoilslow variable
      attr_name = 'slow soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_csoilslow, attr_name, attr_units, MISSING_VALUE)
    
     ! Attributes of nsoilslow variable
      attr_name = 'slow soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nsoilslow, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of csoilpass variable
      attr_name = 'passive soil organic matter carbon'
      attr_units = 'gC m-2'
      call PutVariableAttributeReal(ncid, varid_csoilpass, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nsoilpass variable
      attr_name = 'passive soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nsoilpass, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of tairC variable
      attr_name = 'mean daily air temperature'
      attr_units = 'degrees C'
      call PutVariableAttributeReal(ncid, varid_tairC, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of tsoilC variable
      attr_name = 'mean daily soil temperature in top 50 cm'
      attr_units = 'degrees C'
      call PutVariableAttributeReal(ncid, varid_tsoilC, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of cLitInptMet variable
      attr_name = 'Metabolic Litter C Inputs'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_clitInptMet, attr_name, attr_units, MISSING_VALUE)
 
      ! Attributes of cLitInptStruc variable
      attr_name = 'Structural Litter C Inputs'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_clitInptStruc, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of cpassInpt variable
      attr_name = 'Passive Pool C Inputs'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cpassInpt, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of thetaLiq variable
      attr_name = 'fraction of liquid soil water saturation (0.0-1.0)'
      attr_units = 'fraction'
      call PutVariableAttributeReal(ncid, varid_thetaLiq, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of f(T) variable
      attr_name = 'soil temperature multiplier on decomposition rate (0.0-1.0)'
      attr_units = 'multiplier'
      call PutVariableAttributeReal(ncid, varid_fT, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of f(W) variable
      attr_name = 'soil moisture multiplier on decomposition rate (0.0-1.0)'
      attr_units = 'multiplier'
      call PutVariableAttributeReal(ncid, varid_fW, attr_name, attr_units, MISSING_VALUE)


      !! Added daily N flux and pools -mdh 11/25/2019

      ! Attributes of nLitInptMet variable
      attr_name = 'Metabolic Litter N Inputs'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nlitInptMet, attr_name, attr_units, MISSING_VALUE)
 
      ! Attributes of nLitInptStruc variable
      attr_name = 'Structural Litter N Inputs'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nlitInptStruc, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nMinDep variable
      attr_name = 'N deposition'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nmindep, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nMinFix variable
      attr_name = 'N Fixation'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nminfix, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nMinUptake variable
      attr_name = 'Plant Mineral N Uptake'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nminuptake, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nMinLeach variable
      attr_name = 'Mineral N Leaching'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nminleach, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nMinLoss variable
      attr_name = 'Mineral N Loss (N2O)'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nminloss, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nLitMineralization variable
      attr_name = 'Litter N Mineralization'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nlittermin, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nSoilMineralization variable
      attr_name = 'Soil N Mineralization'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nsmin, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nSoilImmob variable
      attr_name = 'Soil N Immobilization'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nsimm, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nNetMineralization variable
      attr_name = 'Net N mineralization'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nsnet, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of nMineral variable
      attr_name = 'Mineral Soil N'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nminrl, attr_name, attr_units, MISSING_VALUE)


      ! --------------- End the definition phase so that variables can be written to the file ---------------
      status = nf_enddef(ncid)
      if (status /= nf_noerr) call handle_err(status, "enddef")
   
      ! ------------------------------  Write variable values to filename_ncOut_casa ------------------------------
   
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
            print *, 'WritePoolFluxNcFile_casacnp_daily: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
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
      enddo

      ! IGBP PFT
      status =  nf_put_var(ncid, varid_igbp, IGBP_PFT)
      if (status /= nf_noerr) call handle_err(status, "put_var(IGBP_PFT)")
   
      ! Land area 
      status =  nf_put_var(ncid, varid_landarea, landarea)
      if (status /= nf_noerr) call handle_err(status, "put_var(landarea)")

      deallocate(IGBP_PFT,landarea,time,days)

   else

      status = nf_open(filename_ncOut_casa, NF_WRITE, ncid)
      if (status /= nf_noerr) call handle_err(status, "open(filename_ncOut_casa)")

      ! get variable ids when iday > 1

      status = nf_inq_varid(ncid, 'cleaf', varid_cleaf)
      if (status /= nf_noerr) call handle_err(status, "cleaf")

      status = nf_inq_varid(ncid, 'cwood', varid_cwood)
      if (status /= nf_noerr) call handle_err(status, "cwood")

      status = nf_inq_varid(ncid, 'cfroot',varid_cfroot)
      if (status /= nf_noerr) call handle_err(status, "cfroot")

      status = nf_inq_varid(ncid, 'nleaf',varid_nleaf)
      if (status /= nf_noerr) call handle_err(status, "nleaf")

      status = nf_inq_varid(ncid, 'nwood',varid_nwood)
      if (status /= nf_noerr) call handle_err(status, "nwood")

      status = nf_inq_varid(ncid, 'nfroot',varid_nfroot)
      if (status /= nf_noerr) call handle_err(status, "nfroot")

      status = nf_inq_varid(ncid, 'clitmetb',varid_clitmetb)
      if (status /= nf_noerr) call handle_err(status, "clitmetb")

      status = nf_inq_varid(ncid, 'clitstr',varid_clitstr)
      if (status /= nf_noerr) call handle_err(status, "clitstr")

      status = nf_inq_varid(ncid, 'clitcwd',varid_clitcwd)
      if (status /= nf_noerr) call handle_err(status, "clitcwd")

      status = nf_inq_varid(ncid, 'nlitmetb',varid_nlitmetb)
      if (status /= nf_noerr) call handle_err(status, "nlitmetb")

      status = nf_inq_varid(ncid, 'nlitstr',varid_nlitstr)
      if (status /= nf_noerr) call handle_err(status, "nlitstr")

      status = nf_inq_varid(ncid, 'nlitcwd',varid_nlitcwd)
      if (status /= nf_noerr) call handle_err(status, "nlitcwd")

      status = nf_inq_varid(ncid, 'csoilmic',varid_csoilmic)
      if (status /= nf_noerr) call handle_err(status, "csoilmic")

      status = nf_inq_varid(ncid, 'csoilslow',varid_csoilslow)
      if (status /= nf_noerr) call handle_err(status, "csoilslow")

      status = nf_inq_varid(ncid, 'csoilpass',varid_csoilpass)
      if (status /= nf_noerr) call handle_err(status, "csoilpass")

      status = nf_inq_varid(ncid, 'nsoilmic',varid_nsoilmic)
      if (status /= nf_noerr) call handle_err(status, "nsoilmic")

      status = nf_inq_varid(ncid, 'nsoilslow',varid_nsoilslow)
      if (status /= nf_noerr) call handle_err(status, "nsoilslow")

      status = nf_inq_varid(ncid, 'nsoilpass',varid_nsoilpass)
      if (status /= nf_noerr) call handle_err(status, "nsoilpass")

      status = nf_inq_varid(ncid, 'cnpp',varid_cnpp)
      if (status /= nf_noerr) call handle_err(status, "cnpp")

      status = nf_inq_varid(ncid, 'cgpp',varid_cgpp)
      if (status /= nf_noerr) call handle_err(status, "cgpp")

      status = nf_inq_varid(ncid, 'cresp',varid_cresp)
      if (status /= nf_noerr) call handle_err(status, "cresp")

      status = nf_inq_varid(ncid, 'tairC',varid_tairC)
      if (status /= nf_noerr) call handle_err(status, "tairC")

      status = nf_inq_varid(ncid, 'tsoilC',varid_tsoilC)
      if (status /= nf_noerr) call handle_err(status, "tsoilC")

      status = nf_inq_varid(ncid, 'cLitInptMet',varid_clitInptMet)
      if (status /= nf_noerr) call handle_err(status, "cLitInptMet")

      status = nf_inq_varid(ncid, 'cLitInptStruc',varid_clitInptStruc)
      if (status /= nf_noerr) call handle_err(status, "cLitInptStruc")

      status = nf_inq_varid(ncid, 'cpassInpt',varid_cpassInpt)
      if (status /= nf_noerr) call handle_err(status, "cpassInpt")

      status = nf_inq_varid(ncid, 'thetaLiq',varid_thetaLiq)
      if (status /= nf_noerr) call handle_err(status, "thetaLiq")

      status = nf_inq_varid(ncid, 'fT',varid_fT)
      if (status /= nf_noerr) call handle_err(status, "f(T)")

      status = nf_inq_varid(ncid, 'fW',varid_fW)
      if (status /= nf_noerr) call handle_err(status, "f(W)")


      !! Added daily N fluxes and pools -mdh 11/14/2019

      status = nf_inq_varid(ncid, 'nLitInptMet',varid_nlitInptMet)
      if (status /= nf_noerr) call handle_err(status, "nLitInptMet")

      status = nf_inq_varid(ncid, 'nLitInptStruc',varid_nlitInptStruc)
      if (status /= nf_noerr) call handle_err(status, "nLitInptStruc")

      status = nf_inq_varid(ncid, 'nMinDep',varid_nmindep)
      if (status /= nf_noerr) call handle_err(status, "nMinDep")

      status = nf_inq_varid(ncid, 'nMinFix',varid_nminfix)
      if (status /= nf_noerr) call handle_err(status, "nMinFix")

      status = nf_inq_varid(ncid, 'nMinUptake',varid_nminuptake)
      if (status /= nf_noerr) call handle_err(status, "nMinUptake")

      status = nf_inq_varid(ncid, 'nMinLeach',varid_nminleach)
      if (status /= nf_noerr) call handle_err(status, "nMinLeach")

      status = nf_inq_varid(ncid, 'nMinLoss',varid_nminloss)
      if (status /= nf_noerr) call handle_err(status, "nMinLoss")

      status = nf_inq_varid(ncid, 'nLitMineralization',varid_nlittermin)
      if (status /= nf_noerr) call handle_err(status, "nLitMineralization")

      status = nf_inq_varid(ncid, 'nSoilMineralization',varid_nsmin)
      if (status /= nf_noerr) call handle_err(status, "nSoilMineralization")

      status = nf_inq_varid(ncid, 'nSoilImmob',varid_nsimm)
      if (status /= nf_noerr) call handle_err(status, "nSoilImmob")

      status = nf_inq_varid(ncid, 'nNetMineralization',varid_nsnet)
      if (status /= nf_noerr) call handle_err(status, "nNetMineralization")

      status = nf_inq_varid(ncid, 'nMineral',varid_nminrl)
      if (status /= nf_noerr) call handle_err(status, "nMineral")

   endif  !end iday==1

!  Define start3 and count3 for record variables (those with unlimited time dimension)
!  Write up to 6 output variables at a time

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


!  casaflux%Cnpp(npt)         - NPP (gC/m2/day) 
!  casaflux%Cgpp(npt)         - GPP (gC/m2/day) 
!  casaflux%Crsoil(npt)       - soil respiration (gC/m2/day)         
!  casaflux%Crp(npt)          - plant respiration (gC/m2/day)         
!  casapool%cplant(npt,LEAF)  - plant leaf C (gC/m2)
!  casapool%cplant(npt,WOOD)  - plant wood C (gC/m2)
!  casapool%cplant(npt,FROOT) - plant fine root C (gC/m2)


   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_daily: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casaflux%Cnpp(npt)
      var2(ilon,ilat,itime) = casaflux%Crsoil(npt)
      var3(ilon,ilat,itime) = casapool%cplant(npt,LEAF)
      var4(ilon,ilat,itime) = casapool%cplant(npt,WOOD)
      var5(ilon,ilat,itime) = casapool%cplant(npt,FROOT)
      var6(ilon,ilat,itime) = casaflux%Cgpp(npt)
   enddo


!! ATTENTION: I could not get the variables with a "time" dimension to write 
!! to the netcdf file when the time dimension was unlimited, UNLESS I substituted 
!! one "nf_put_vara_real" for a "nf_put_var".  I DON'T UNDERSTAND!
!! Otherwise nf_put_var seemed to ignore start3 and count3.
!! Melannie 6/3/2014
!! Had to use "nf_put_vara_real" for all variables to avoid segmentation violation
!! Melannie 6/6/2016

!! status =  nf_put_var(ncid, varid_cnpp, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cnpp, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(cnpp)")
 
!! status =  nf_put_var(ncid, varid_cgpp, var6, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cgpp, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(cgpp)")

!! status =  nf_put_var(ncid, varid_cresp, var2, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cresp, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(cresp)")
 
!! status =  nf_put_var(ncid, varid_cleaf, var3, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cleaf, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cleaf)")

!! status =  nf_put_var(ncid, varid_cwood, var4, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cwood, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(cwood)")

!! status =  nf_put_var(ncid, varid_cfroot, var5, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cfroot, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(cfroot)")


!  casapool%clitter(npt,METB) - metabolic litter C (gC/m2)
!  casapool%clitter(npt,STR)  - structural litter C (gC/m2)
!  casapool%clitter(npt,CWD)  - coarse woody debris C (gC/m2)
!  casapool%csoil(npt,MIC)    - microbial soil  C (gC/m2)
!  casapool%csoil(npt,SLOW)   - slow soil C (gC/m2)
!  casapool%csoil(npt,PASS)   - passive soil C (gC/m2)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_daily: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapool%clitter(npt,METB)
      var2(ilon,ilat,itime) = casapool%clitter(npt,STR)
      var3(ilon,ilat,itime) = casapool%clitter(npt,CWD)
      var4(ilon,ilat,itime) = casapool%csoil(npt,MIC)
      var5(ilon,ilat,itime) = casapool%csoil(npt,SLOW)
      var6(ilon,ilat,itime) = casapool%csoil(npt,PASS)
   enddo

!! status =  nf_put_var(ncid, varid_clitmetb, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_clitmetb, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(clitmetb)")
 
!! status =  nf_put_var(ncid, varid_clitstr, var2, start3, count3)
   status =  nf_put_vara_real(ncid, varid_clitstr, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(clitstr)")
 
!! status =  nf_put_var(ncid, varid_clitcwd, var3, start3, count3)
   status =  nf_put_vara_real(ncid, varid_clitcwd, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(clitcwd)")

!! status =  nf_put_var(ncid, varid_csoilmic, var4, start3, count3)
   status =  nf_put_vara_real(ncid, varid_csoilmic, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(csoilmic)")

!! status =  nf_put_var(ncid, varid_csoilslow, var5, start3, count3)
   status =  nf_put_vara_real(ncid, varid_csoilslow, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(csoilslow)")

!! status =  nf_put_var(ncid, varid_csoilpass, var6, start3, count3)
   status =  nf_put_vara_real(ncid, varid_csoilpass, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(csoilpass)")


!  casapool%nplant(npt,LEAF)  - plant leaf N (gN/m2)
!  casapool%nplant(npt,WOOD)  - plant wood N (gN/m2)
!  casapool%nplant(npt,FROOT) - plant fine root N (gN/m2)
!  casapool%nlitter(npt,METB) - metabolic litter N (gN/m2)
!  casapool%nlitter(npt,STR)  - structural litter N (gN/m2)
!  casapool%nlitter(npt,CWD)  - coarse woody debris N (gN/m2)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_daily: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapool%nplant(npt,LEAF)
      var2(ilon,ilat,itime) = casapool%nplant(npt,WOOD)
      var3(ilon,ilat,itime) = casapool%nplant(npt,FROOT)
      var4(ilon,ilat,itime) = casapool%nlitter(npt,METB)
      var5(ilon,ilat,itime) = casapool%nlitter(npt,STR)
      var6(ilon,ilat,itime) = casapool%nlitter(npt,CWD)
   enddo

!! status =  nf_put_var(ncid, varid_nleaf, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nleaf, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(nleaf)")
 
!! status =  nf_put_var(ncid, varid_nwood, var2, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nwood, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(nwood)")
 
!! status =  nf_put_var(ncid, varid_nfroot, var3, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nfroot, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nfroot)")

!! status =  nf_put_var(ncid, varid_nlitmetb, var4, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nlitmetb, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(nlitmetb)")

!! status =  nf_put_var(ncid, varid_nlitstr, var5, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nlitstr, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(nlitstr)")

!! status =  nf_put_var(ncid, varid_nlitcwd, var6, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nlitcwd, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(nlitcwd)")


!  casapool%nsoil(npt,MIC)    - microbial soil N (gN/m2)
!  casapool%nsoil(npt,SLOW)    - slow soil N (gN/m2)
!  casapool%nsoil(npt,PASS)    - passive soil N (gN/m2)
!  casamet%tairk(npt)           - air temperature (K)
!  casamet%tsoilavg(npt)        - soil temperature (K)
!  casaflux%ClitInptMet(npt)    - metabolic litter inputs (gC/m2/dy)
!  casaflux%ClitInptStruc(npt)  - structural litter inputs (gC/m2/dy)
!  casaflux%CpassInpt(npt)      - inputs to passive soil pool (gC/m2/dy)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE
   var8(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_daily: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapool%nsoil(npt,MIC)
      var2(ilon,ilat,itime) = casapool%nsoil(npt,SLOW)
      var3(ilon,ilat,itime) = casapool%nsoil(npt,PASS)
      var4(ilon,ilat,itime) = casamet%tairk(npt) - tkzeroc
      var5(ilon,ilat,itime) = casamet%tsoilavg(npt) - tkzeroc
      var6(ilon,ilat,itime) = casaflux%ClitInptMet(npt)
      var7(ilon,ilat,itime) = casaflux%ClitInptStruc(npt)
      var8(ilon,ilat,itime) = casaflux%CpassInpt(npt)
!     write(*,*) 'casaflux%ClitInptStruc(',npt,')=', casaflux%ClitInptStruc(npt)
!     write(*,*) 'var7(',ilon,ilat,itime,')=', var7(ilon,ilat,itime)
   enddo

!! status =  nf_put_var(ncid, varid_nsoilmic, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nsoilmic, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(nsoilmic)")
 
!! status =  nf_put_var(ncid, varid_nsoilslow, var2, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nsoilslow, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(nsoilslow)")
 
!! status =  nf_put_var(ncid, varid_nsoilpass, var3, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nsoilpass, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nsoilpass)")

!! status =  nf_put_var(ncid, varid_tairC, var4, start3, count3)
   status =  nf_put_vara_real(ncid, varid_tairC, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(tairC)")

!! status =  nf_put_var(ncid, varid_tsoilC, var5, start3, count3)
   status =  nf_put_vara_real(ncid, varid_tsoilC, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(tsoilC)")

   status =  nf_put_vara_real(ncid, varid_clitInptMet, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(cLitInptMet)")

   status =  nf_put_vara_real(ncid, varid_clitInptStruc, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "put_var(cLitInptStruc)")

   status =  nf_put_vara_real(ncid, varid_cpassInpt, start3, count3, var8)
   if (status /= nf_noerr) call handle_err(status, "put_var(cpassInpt)")


!  f(T) and f(W) output added 11/20/2017. ThetaLiq added 11/27/2017
!  casapool%fT(npt)  - daily soil temperature multiplier on decomposition rate (0.0-1.0+)
!  casapool%fW(npt)  - daily soil moisture multiplier on decomposition rate (0.0-1.0)
!  casapool%thetaLiq(npt)  - fraction of liquid soil water saturation (0.0-1.0)
!  casapoolAn%Nsoilmin(npt) - daily mineral soil N (gN/m2)
!  casaflux%NlitInptMet(npt) - daily metabolic litter N inputs (gN/m2/yr)
!  casaflux%NlitInptStruc(npt) - daily structural litter N inputs (gN/m2/yr)
!  casaflux%Nmindep(npt) - daily N deposition (gN/m2/yr)
!  casaflux%Nminfix(npt) - daily N fixation (gN/m2/yr)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE
   var8(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casapool%fT(npt)
      var2(ilon,ilat,itime) = casapool%fW(npt)
      var3(ilon,ilat,itime) = casapool%thetaLiq(npt)

      var4(ilon,ilat,itime) = casapool%Nsoilmin(npt)
      var5(ilon,ilat,itime) = casaflux%NlitInptMet(npt)
      var6(ilon,ilat,itime) = casaflux%NlitInptStruc(npt)
      var7(ilon,ilat,itime) = casaflux%Nmindep(npt)
      var8(ilon,ilat,itime) = casaflux%Nminfix(npt)
   enddo

!! status =  nf_put_var(ncid, varid_fT, var1, start3, count3)
   status =  nf_put_vara_real(ncid, varid_fT, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(fT)")
 
   status =  nf_put_vara_real(ncid, varid_fW, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(fW)")

   status =  nf_put_vara_real(ncid, varid_thetaLiq, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(thetaLiq)")

   status =  nf_put_vara_real(ncid, varid_nminrl, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMineral)")

   status =  nf_put_vara_real(ncid, varid_nlitinptmet, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(NlitInptMet)")

   status =  nf_put_vara_real(ncid, varid_nlitinptstruc, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(NlitInptStruc)")

   status =  nf_put_vara_real(ncid, varid_nmindep, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinDep)")

   status =  nf_put_vara_real(ncid, varid_nminfix, start3, count3, var8)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinFix)")


!  casaflux%Nminuptake(npt) - Plant N uptake (gN/m2/day)
!  casaflux%Nminleach(npt) - Mineral N Leaching (gN/m2/day)
!  casaflux%Nminloss(npt) - Mineral N Loss (gN/m2/day)
!  casaflux%Nlittermin(npt) - Litter N mineralization (gN/m2/day)
!  casaflux%Nsmin(npt) - Soil N mineralization (gN/m2/day)
!  casaflux%Nsimm(npt) - Soil N Immobilization (gN/m2/day)
!  casaflux%Nsnet(npt) - Net N mineralization (gN/m2/day)

   var1(:,:,:) = MISSING_VALUE
   var2(:,:,:) = MISSING_VALUE
   var3(:,:,:) = MISSING_VALUE
   var4(:,:,:) = MISSING_VALUE
   var5(:,:,:) = MISSING_VALUE
   var6(:,:,:) = MISSING_VALUE
   var7(:,:,:) = MISSING_VALUE

   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      if (casamet%ijgcm(npt) .ne. clmgrid%cellid(ilon,ilat)) then
         print *, 'WritePoolFluxNcFile_casacnp_annual: casamet%ijgcm(', npt, ')=', casamet%ijgcm(npt)
         print *, '   clmgrid%cellid(', ilon, ',', ilat, ')=', clmgrid%cellid(ilon,ilat)
         STOP
      endif
      var1(ilon,ilat,itime) = casaflux%Nminuptake(npt)
      var2(ilon,ilat,itime) = casaflux%Nminleach(npt)
      var3(ilon,ilat,itime) = casaflux%Nminloss(npt)
      var4(ilon,ilat,itime) = casaflux%Nlittermin(npt)
      var5(ilon,ilat,itime) = casaflux%Nsmin(npt)
      var6(ilon,ilat,itime) = casaflux%Nsimm(npt)
      var7(ilon,ilat,itime) = casaflux%Nsnet(npt)
   enddo

   status =  nf_put_vara_real(ncid, varid_nminuptake, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinUptake)")

   status =  nf_put_vara_real(ncid, varid_nminleach, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinLeach)")

   status =  nf_put_vara_real(ncid, varid_nminloss, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMinLoss)")

   status =  nf_put_vara_real(ncid, varid_nlittermin, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLitMineralization)")

   status =  nf_put_vara_real(ncid, varid_nsmin, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSoilMineralization)")

   status =  nf_put_vara_real(ncid, varid_nsimm, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSoilImmob)")

   status =  nf_put_vara_real(ncid, varid_nsnet, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "put_var(nNetMineralization)")


   deallocate(var1,var2,var3,var4,var5,var6,var7,var8)

   status = nf_close(ncid)

   if (verbose .ge. 1) print *, "Done writing output to ", trim(filename_ncOut_casa), " on day", iday, "..."

END SUBROUTINE WritePoolFluxNcFile_casacnp_daily

!----------------------------------------------------------------------------------------------------
SUBROUTINE casacnpdriver(filename_cnpmet, filename_cnpepool, filename_cnpflux, filename_ncOut_casa, &
                         filename_mimicsepool, filename_ncOut_mimics, &
                         filename_corpseepool, filename_ncOut_corpse, &
                         mloop, mdaily, co2air, deltsoil, deltair, deltYr, nppMult, calyr)
   use define_dimensions
   use define_types
   use casadimension
   use casaparm
   use casavariable
   use mimicsvariable
   use corpsevariable
   implicit none

   !Subroutine arguments
   character(len=100), INTENT(IN) :: filename_cnpmet
   character(len=100), INTENT(IN) :: filename_cnpepool
   character(len=100), INTENT(IN) :: filename_cnpflux
   character(len=100), INTENT(IN) :: filename_ncOut_casa
   character(len=100), INTENT(IN) :: filename_mimicsepool
   character(len=100), INTENT(IN) :: filename_ncOut_mimics
   character(len=100), INTENT(IN) :: filename_corpseepool
   character(len=100), INTENT(IN) :: filename_ncOut_corpse

   ! mloop = sumber of times to cycle through the current met.nc file.
   ! mdaily = 0 for annual output, =1 for daily output
   ! deltYr = the value of iYrCnt (defined below) to begin using deltsoil and deltair
   ! calyr = the calendar year (transient runs only)

   integer,   INTENT(IN) :: mloop, mdaily, deltYr, calyr

   ! co2air = atmospheric CO2 concentration (ppm)
   ! deltsoil = add this to soil temperature when iYrCnt >= deltYr
   ! deltair  = add this to air temperature when iYrCnt >= deltYr
   ! nppMult  = multiple NPP by this when iYrCnt >= deltYr

   real(r_2), INTENT(IN) :: co2air, deltsoil, deltair, nppMult

   !Local Variables
   integer :: nloop,npt,nyear,nday,ns,iyear,imon,iday,iday1,iday2,ndoy,idoy
   integer :: ndoy1,ndoy2
   integer :: myear          ! Number of years in met.nc file
   integer :: iYrCnt         ! Counter of years executed in this subroutine call (1..mloop*myear)
   integer :: wrtYr          ! Time variable value in the output NetCDF files associated with current simulation year
   integer :: linesWritten   ! For CORPSE, used to determine the the position to write the current year's daily values in NetCDF file.
   logical writeAnSpinNcOutput, writeToRestartCSVfile
   
   ! Names if daily and annual NetCDF output files for the 3 SOM models
   character(len=100) :: filename_ncOut_casa_spin_yr
   character(len=100) :: filename_ncOut_mimics_spin_yr
   character(len=100) :: filename_ncOut_corpse_spin_yr
   character(len=100) :: filename_ncOut_corpse_spin_day
   character(len=100) :: filename_ncOut_casa_day
   character(len=100) :: filename_ncOut_mimics_day
   character(len=100) :: filename_ncOut_corpse_day

   real(r_2), dimension(:,:,:), allocatable   :: xtsoil     ! daily soil temperature (K) (read from met.nc file) 
   real(r_2), dimension(:,:,:), allocatable   :: xmoist     ! daily volumetric liquid soil moisture (read from met.nc file) 
   real(r_2), dimension(:,:,:), allocatable   :: xfrznmoist ! daily volumetric frozen soil moisture (read from met.nc file) 
   real(r_2), dimension(:,:),   allocatable   :: xtairk     ! daily air temperature (K) (read from met.nc file)
   real(r_2), dimension(:,:),   allocatable   :: xndepDay   ! daily atmospheric N deposition (read from met.nc file)
   real(r_2), dimension(:,:),   allocatable   :: xcnpp      ! daily NPP for the current simulation year (calculated day by day)
   real(r_2), dimension(:,:),   allocatable   :: xcgpp      ! daily GPP for the current simulation year (read from met.nc)
!  real(r_2), dimension(:,:),   allocatable   :: xlai
!  character(len=100) fnpp,fwb,ftgg,ftairk
!  real(r_2), dimension(:),     allocatable   :: totynpp
   real(r_2), dimension(:,:),   allocatable   :: monlai
   integer mdoy(12),middoy(12)
!  real(r_2), dimension(mp,12)  :: FCgppmon,FCnppmon,FCrsmon,FCneemon 
   real(r_2), dimension(mp,3)   :: clitterinput,csoilinput
   real(r_2), dimension(mp)     :: xnlimit,xplimit,xnplimit,xnuptake,xpuptake,xnpuptake
   real(r_2), dimension(mp)     :: yavgtair,betaco2
   integer,   dimension(mp)     :: nyavgtair
   real(r_2)                    :: xratio,co2cp
   real(r_2)                    :: nppScalar     ! Multiplier on NPP: this is set to 1.0 or the value nppMult
   real(r_2), dimension(mp)     :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd,nwd2str

   data mdoy/31,59,90,120,151,181,212,243,273,304,334,365/
   data middoy/15,46,74,105,135,166,196,227,258,288,319,349/
!   data mdays/31,28,31,30,31,30,31,31,30,31,30,31/

101 format(a100)

!     open(98,file='laimonth_mod.csv')
      allocate(monlai(mp,12))
!     allocate(totynpp(mp))
!     totynpp=0.0

      !---------------------------------------------------------------------------------------------------
      ! Open and read all of the met.nc file
      !---------------------------------------------------------------------------------------------------
      
      ! Get the dimensions of variables for use in allocation
      call GetMetNcFileDim(filename_cnpmet, ms, myear)

      !Segmentation fault occurred on this allocation when it was in ReadMetNcFile.
      !Variables are allocated here instead.
!     allocate(xlai(1:mp,1:myear*mdyear))
      allocate(xcnpp(1:mp,1:myear*mdyear))
      allocate(xcgpp(1:mp,1:myear*mdyear))
      allocate(xtairk(1:mp,1:myear*mdyear))
      allocate(xndepDay(1:mp,1:myear*mdyear))
      allocate(xtsoil(1:mp,1:ms,1:myear*mdyear))
      allocate(xmoist(1:mp,1:ms,1:myear*mdyear))
      allocate(xfrznmoist(1:mp,1:ms,1:myear*mdyear))
      !xlai(:,:) = 0.0
      xcnpp(:,:) = 0.0
      xcgpp(:,:) = 0.0
      xtairk(:,:) = 0.0
      xndepDay(:,:) = 0.0
      xtsoil(:,:,:) = 0.0
      xmoist(:,:,:) = 0.0
      xfrznmoist(:,:,:) = 0.0

!     call ReadMetNcFile(filename_cnpmet, mp, ms, myear, xlai, xcnpp, xcgpp, &
!                        xtairk, xndepDay, xtsoil, xmoist, xfrznmoist)
      call ReadMetNcFile(filename_cnpmet, mp, ms, myear, xcgpp, &
                         xtairk, xndepDay, xtsoil, xmoist, xfrznmoist)
      xcnpp = xcgpp / 2.0  !Initialization only

      ! Compute the average air temperature over the entire simulation cycle.
      yavgtair(:) = 0.0
      nyavgtair(:) = 0
      do nyear=1,myear
         ndoy1=(nyear-1)*365+1
         ndoy2=nyear*365
         do npt=1,mp
            do iday=ndoy1,ndoy2
               if(xtairk(npt,iday) > 273.12) then
                  nyavgtair(npt) = nyavgtair(npt) + 1
                  yavgtair(npt)  = yavgtair(npt) + xtairk(npt,iday) - 273.12
               endif
            enddo
         enddo
      enddo   

      ! Set ice land point NPP to zero
      do npt=1,mp
         if(casamet%iveg2(npt)==icewater) then
            xcnpp(npt,:)=0.0
            xcgpp(npt,:)=0.0
            yavgtair(npt) = 0.0
         else
           yavgtair(npt) = deltair + yavgtair(npt) /real(nyavgtair(npt))
         endif
         co2cp=42.7+1.68*(yavgtair(npt)-25.0)+0.012*(yavgtair(npt)-25.0)*(yavgtair(npt)-25.0)
         xratio = 0.7
         betaco2(npt) = 3.0*xratio*co2air*co2cp &
                      /((xratio*co2air-co2cp)*(xratio*co2air+2.0*co2cp)) 

         ! Adjustment for savannas? (MDH)
         if(veg%iveg(npt)==8.or.veg%iveg(npt)==9) then
            betaco2(npt) = betaco2(npt)*0.5
         endif

!        write(*,717) npt,veg%iveg(npt),yavgtair(npt),casamet%lat(npt),casamet%lon(npt),betaco2(npt)
      enddo

717 format(i4,2x,i3,2x,4(f10.3,2x))

  xnplimit = 1.0
  xnpuptake = 1.0

  !!The monthly values are not used anywhere. -mdh 3/9/2020
  !FCgppmon=0.0;FCnppmon=0.0;FCrsmon=0.0;FCneemon=0.0

  iYrCnt = 0    ! Count the number of years (1..mloop*myear)
  DO nloop=1,mloop

     ! If writeAnSpinNcOutput == .true., write to netcdf files each year during a spinup
     ! Eventually I will need to make writeAnSpinNcOutput an option in the .lst file
     writeAnSpinNcOutput = .true.

     do iyear=1,myear

        !! Initialize average annual fluxes for the year.
        !! For output only. These are accumulated in casa_cnpflux.
        !! The annual means are computed in casa_fluxout.
        casabal%FCgppyear=0.0; casabal%FCrpyear=0.0
        casabal%FCnppyear=0;   casabal%FCrsyear=0.0;    casabal%FCneeyear=0.0
        casabal%FNdepyear=0.0; casabal%FNfixyear=0.0;   casabal%FNsnetyear=0.0
        casabal%FNupyear=0.0;  casabal%FNleachyear=0.0; casabal%FNlossyear=0.0
        casabal%FPweayear=0.0; casabal%FPdustyear=0.0;  casabal%FPsnetyear=0.0
        casabal%FPupyear=0.0;  casabal%FPleachyear=0.0; casabal%FPlossyear=0.0 
        clitterinput = 0.0;    csoilinput = 0.0    
   
        !! Initialize average daily pool values over the most recent myear years for 
        !! CASACNP. For output only. These daily values are accumulated in casa_cnppool 
        !! and the daily mean is computed in casa_poolout.  -MDH 9/29/2014
        casapoolAn%CsoilAn=0.0; casapoolAn%CplantAn=0.0; casapoolAn%ClitterAn=0.0
        casapoolAn%NsoilAn=0.0; casapoolAn%NplantAn=0.0; casapoolAn%NlitterAn=0.0
        casapoolAn%PsoilAn=0.0; casapoolAn%PplantAn=0.0; casapoolAn%PlitterAn=0.0
        casaflux%ClitInptMetAn = 0.0; casaflux%ClitInptStrucAn = 0.0; casaflux%CpassInptAn = 0.0
        casaflux%NlitInptMetAn = 0.0; casaflux%NlitInptStrucAn = 0.0
        casapoolAn%tairAn=0.0; casapoolAn%tsoilAn=0.0
        casapoolAn%fTAn=0.0; casapoolAn%fWAn=0.0; casapoolAn%thetaLiqAn=0.0;
        !! Added annual N pools and fluxes 11/25/2019 -mdh
        casapoolAn%NsoilminAn = 0.0; casaflux%NmindepAn = 0.0; casaflux%NminfixAn = 0.0
        casaflux%NminuptakeAn = 0.0; casaflux%NminlossAn = 0.0; casaflux%NminleachAn = 0.0
        casaflux%NlitterminAn = 0.0; casaflux%NsminAn = 0.0; casaflux%NsimmAn = 0.0; casaflux%NsnetAn = 0.0
   
        if (isomModel == MIMICS) then
            !! Initialize average annual fluxes and average daily pool values over the 
            !! most recent myear years for MIMICS. For output only. These fluxes and pools 
            !! are accumulated in mimics_caccum (mimics_cycle.f90). The means are computed 
            !! in mimics_poolfluxout (mimics_cycle.f90). -MDH 01/26/2015
            mimicsfluxAn%ClitInputAn=0.0; mimicsfluxAn%ChrespAn=0.0 ;  mimicsfluxAn%CSOMpInputAn=0.0
            mimicspoolAn%ClitterAn=0.0;   mimicspoolAn%CmicrobeAn=0.0; mimicspoolAn%CsoilAn=0.0
            mimicspoolAn%NlitterAn=0.0;   mimicspoolAn%NmicrobeAn=0.0; mimicspoolAn%NsoilAn=0.0; 
            mimicsfluxAn%NlitInputAn=0.0; mimicspoolAn%fTAn=0.0;       mimicspoolAn%fWAn=0.0; 
            mimicspoolAn%thetaLiqAn=0.0;  mimicspoolAn%thetaFrznAn=0.0
        endif
   
        iYrCnt = iYrCnt + 1
        if (initcasa < 2) then
           wrtYr = iYrCnt
        else
           wrtYr = calyr
        endif
        if (iYrCnt >= deltYr) then
           write(*,34) 'Soil Temperature Increment for year', iYrCnt, ' = ', deltsoil
           write(*,34) 'Air Temperature Increment for year', iYrCnt, ' = ', deltair
           write(*,34) 'NPP multiplier', iYrCnt, ' = ', nppMult
34 format(a30,1x,i4,a3,f6.2)
        endif
 
        if(mdaily == 1) then
           ! This is the beginning of another simulation year.
           ! Determine the names of the daily .nc output files for the current year.
           ! The daily names are derived from the annual names.
           if (initcasa < 2) then
              !Non-transint run.  Use iYrCnt for simulation year
              call CreateDailyNcFileName(filename_ncOut_casa_day, filename_ncOut_casa, iYrCnt)
              if (isomModel == MIMICS) then
                 call CreateDailyNcFileName(filename_ncOut_mimics_day, filename_ncOut_mimics, iYrCnt)
              endif
              if (isomModel == CORPSE) then
                 call CreateDailyNcFileName(filename_ncOut_corpse_day, filename_ncOut_corpse, iYrCnt)
              endif
           else
              ! Transient file names already include calendar year so do not replace year with iYrCnt.
              call CreateDailyNcFileName(filename_ncOut_casa_day, filename_ncOut_casa, -1)
              if (isomModel == MIMICS) then
                 call CreateDailyNcFileName(filename_ncOut_mimics_day, filename_ncOut_mimics, -1)
              endif
              if (isomModel == CORPSE) then
                 call CreateDailyNcFileName(filename_ncOut_corpse_day, filename_ncOut_corpse, -1)
              endif
           endif
        endif  ! if (mdaily==1)

        monlai = -1.0
  
        !Annual NPP from the previous year is needed for MIMICS and CORPSE models
        !Only do this calculation if CnppAn has not yet been initialized. -mdh 3/9/2020
        do npt=1,mp
           if (casaflux%CnppAn(npt) <= -99.0) then
              casaflux%CnppAn(npt) = SUM(xcnpp(npt,(iyear-1)*365+1:(iyear-1)*365+365))
           endif
        enddo

        do imon=1,12
           if(imon==1) then
              iday1=1
           else
              iday1=mdoy(imon-1)+1
           endif
           iday2=mdoy(imon)
  
           do iday=iday1,iday2

              ndoy=(iyear-1)*365+iday
              casaflux%Cgpp(:)  = xcgpp(:,ndoy)                ! Added -mdh 3/24/2014 
              casaflux%Cnpp(:)  = xcnpp(:,ndoy) * xnplimit(:)  ! casaflux%Cnpp(:) is reassigned in SUBROUTINE casa_rplant
              if (iYrCnt >= deltYr) then
                 casamet%tairk(:)  = xtairk(:,ndoy)    + deltair
                 casamet%tsoil(:,:)= xtsoil(:,:,ndoy)  + deltsoil
              else
                 casamet%tairk(:)  = xtairk(:,ndoy)
                 casamet%tsoil(:,:)= xtsoil(:,:,ndoy)
              endif
              casamet%moist(:,:)= xmoist(:,:,ndoy)
              casamet%frznmoist(:,:)= xfrznmoist(:,:,ndoy)  ! Added 3/13/2017
              casaflux%Nmindep(:) = xndepDay(:,ndoy) ! Added 8/22/2016 -mdh. Overwrites value read from gridinfo_igbpz.csv
    
              !! Updated the call to biogeochem (-MDH 6/9/2014)
              !!call biogeochem(iday,xnlimit,xplimit,xnplimit,xnuptake,xpuptake,xnpuptake)
              !!call biogeochem(iday,veg,soil,casabiome,casapool,casaflux, &
              !!                casamet,casabal,phen)
              if (iYrCnt >= deltYr) then
                 nppScalar = nppMult
              else
                 nppScalar = 1.0
              endif

              !call biogeochem(iYrCnt,iday,nppScalar)
              !When calling biogeochem, pass back plant litter fluxes for daily output. -mdh 5/14/2018
              call biogeochem(iYrCnt,iday,mdaily,nppScalar,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)

              if (casafile%iptToSaveIndx > 0) then
                 ! Write point-specific output. Write end-of-year values only when mdaily==0. -mdh 5/14/2018
                 if ((mdaily == 1) .or. (iday == 365)) then
                     call WritePointCASA(iYrCnt,iday,mp,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)
                 endif
              endif

              !! casa_cnpflux is called in biogeochem now (-MDH 6/30/2014)
              !!call casa_cnpflux(clitterinput,csoilinput)
              !casapool%Nsoilmin(:) = max(0.1,casapool%Nsoilmin(:))
              !casapool%Psoillab(:) = max(0.01,casapool%Psoillab(:))
              if(nloop==mloop) then
                 !! FCnppmon, FCrsmon, and FCneemon are not being used anywhere. -mdh 3/9/2020
                 !FCnppmon(:,imon) = FCnppmon(:,imon) + casaflux%Cnpp(:)   * deltpool
                 !FCrsmon(:,imon)  = FCrsmon(:,imon)  + casaflux%Crsoil(:) * deltpool
                 !FCneemon(:,imon) = FCneemon(:,imon) + (-casaflux%Cnpp(:)+casaflux%Crsoil(:)) * deltpool
                 if(iday==middoy(imon)) monlai(:,imon) = casamet%glai(:)
              endif 

              if(mdaily == 1) then
                 !Write current day's values to NetCDF output file
                 !Annual output occurs below, after the year loop.
                 if (initcasa > 1) then
                    ! Write daily output to NetCDF file each year for transient runs
                    call WritePoolFluxNcFile_casacnp_daily(filename_ncOut_casa_day, mp, wrtYr, iday)
                    if (isomModel == MIMICS) then
                       call WritePoolFluxNcFile_mimics_daily(filename_ncOut_mimics_day, mp, wrtYr, iday)
                    else if (isomModel == CORPSE) then
                       ! CORPSE daily NetCDF output occurs at the end year, all 365 days at once
                    endif
                !!else if (iyear == myear .and. writeAnSpinNcOutput) then
                else if (writeAnSpinNcOutput) then
                   ! Write daily spinup files each year, not just every myear years.  -mdh 9/19/2016
                   call WritePoolFluxNcFile_casacnp_daily(filename_ncOut_casa_day, mp, wrtYr, iday)
                   if (isomModel == MIMICS) then
                      call WritePoolFluxNcFile_mimics_daily(filename_ncOut_mimics_day, mp, wrtYr, iday)
                   else if (isomModel == CORPSE) then
                      ! CORPSE daily NetCDF output occurs at the end year, all 365 days at once
                   endif
                endif
             endif

          enddo  !end do iday
       enddo  ! end do imon

      !END OF SIMULATION YEAR

!      !Write to laimonth_mod.csv
!      if(nloop==mloop) then
!         do npt=1,mp
!            write(98,981) npt,iyear,veg%iveg(npt),soil%isoilm(npt), &
!                  casamet%lat(npt),casamet%lon(npt),casamet%areacell(npt)*(1.0e-9), &
!                  monlai(npt,:)
!         enddo
!      endif

       !!if (initcasa < 2 .and. writeAnSpinNcOutput .and. iyear .eq. myear) then
       ! Write annual spinup files each year, not just every myear years.  -mdh 9/19/2016
       ! Daily output occurs above, within the day loop.
       ! Output for transient runs occurs "at the end of the simulation" since transient 
       ! simulations are one year long.
       if (initcasa < 2 .and. writeAnSpinNcOutput) then

          filename_ncOut_casa_spin_yr = filename_ncOut_casa
          filename_ncOut_mimics_spin_yr = filename_ncOut_mimics
          filename_ncOut_corpse_spin_yr = filename_ncOut_corpse

          call InsertYearInNcFileName(filename_ncOut_casa_spin_yr, iYrCnt)
          call InsertYearInNcFileName(filename_ncOut_mimics_spin_yr, iYrCnt)
          call InsertYearInNcFileName(filename_ncOut_corpse_spin_yr, iYrCnt)
          call CreateDailyNcFileName(filename_ncOut_corpse_spin_day, filename_ncOut_corpse, iYrCnt)
  
          !! Compute annual means in casa_poolout and casa_fluxout but only write 
          !! to output CSV files every 100 years.
          if (MOD((nloop-1)*myear+iyear, 100) == 0) then
             writeToRestartCSVfile = .true.
             write(*,*) 'Writing to restart files year: ', (nloop-1)*myear+iyear
          else
             writeToRestartCSVfile = .false.
          endif
          call casa_poolout(filename_cnpepool,iYrCnt,myear,writeToRestartCSVfile)
          call casa_fluxout(filename_cnpflux,myear,clitterinput,csoilinput,writeToRestartCSVfile)

          if (mdaily /= 1) then
             !NOTE: The daily .nc output for CASACNP occurs above within the iday loop
             call WritePoolFluxNcFile_casacnp_annual(filename_ncOut_casa_spin_yr, mp, wrtYr)
          endif

          if (isomModel == MIMICS) then
             if (icycle == 1) then
                call mimics_poolfluxout(filename_mimicsepool,mp,iYrCnt,myear,writeToRestartCSVfile)
             else
                call mimics_poolfluxout_CN(filename_mimicsepool,mp,iYrCnt,myear,writeToRestartCSVfile)
             endif
             if (mdaily /= 1) then
                !NOTE: The daily .nc output for MIMICS occurs above within the iday loop
                call WritePoolFluxNcFile_mimics_annual(filename_ncOut_mimics_spin_yr, mp, wrtYr)
             endif
          else if (isomModel == CORPSE) then
             ! Output current year's CORPSE results for non-transient run(-mdh 5/16/2016)
             call corpse_poolfluxout(filename_corpseepool,mp,writeToRestartCSVfile)
             if (mdaily == 1) then
                linesWritten = pt(mp)%litterlayer_outputs%linesWritten  ! all pools and points have the same # of output lines
                call InitPoolFluxNcFile_corpse(filename_ncOut_corpse_spin_day, 365, mdaily)
                call WritePoolFluxNcFile_corpse(filename_ncOut_corpse_spin_day, pt(mp)%litterlayer_outputs, &
                                                pt(mp)%soil_outputs(1), mp, linesWritten-364, linesWritten, mdaily, wrtYr)
             else
                linesWritten = pt(mp)%litterlayer_outputs%nYear         ! all pools and points have the same # of output years
                call InitPoolFluxNcFile_corpse(filename_ncOut_corpse_spin_yr, 1, mdaily)
                call WritePoolFluxNcFile_corpse(filename_ncOut_corpse_spin_yr, pt(mp)%litterlayer_outputs, &
                                                pt(mp)%soil_outputs(1), mp, linesWritten, linesWritten, mdaily, wrtYr)
             endif
          endif
        endif  ! if writeAnSpinNcOutput

        !Annual NPP from the previous year is needed for MIMICS and CORPSE models.
        !Calculation added here. -mdh 3/9/2020
        casaflux%CnppAn(:) = casabal%FCnppyear(:)

     enddo     ! end do iyear

     !END OF MET.NC FILE (iyear=myear)

     !Write END-OF-SIMULATION-CYCLE pools and fluxes to NetCDF file
     if(nloop==mloop) then

        ! This is the end of the simulation cyle.  Always write to CSV files and netcdf files.
        writeToRestartCSVfile = .true.
        call casa_poolout(filename_cnpepool,iYrCnt,myear,writeToRestartCSVfile)
        call casa_fluxout(filename_cnpflux,myear,clitterinput,csoilinput,writeToRestartCSVfile)
        if (mdaily /= 1) then
           call WritePoolFluxNcFile_casacnp_annual(filename_ncOut_casa, mp, wrtYr)
        endif

        if (isomModel == MIMICS) then

           if (icycle == 1) then
              call mimics_poolfluxout(filename_mimicsepool,mp,iYrCnt,myear,writeToRestartCSVfile)
           else
              call mimics_poolfluxout_CN(filename_mimicsepool,mp,iYrCnt,myear,writeToRestartCSVfile)
           endif
           if (mdaily /= 1) then
              call WritePoolFluxNcFile_mimics_annual(filename_ncOut_mimics, mp, wrtYr)
           endif

        else if (isomModel == CORPSE) then

           call corpse_poolfluxout(filename_corpseepool,mp,writeToRestartCSVfile)
           if (casafile%iptToSaveIndx > 0) then
              !Write point-specific output to sPtFileNameCORPSE
              call WritePointCORPSE(sPtFileNameCORPSE,casafile%iptToSaveIndx)
           endif
           ! Output current year's CORPSE results for transient run (-mdh 5/16/2016)
           if (mdaily == 1) then
              linesWritten = pt(mp)%litterlayer_outputs%linesWritten  ! all pools and points have the same # of output lines
              call InitPoolFluxNcFile_corpse(filename_ncOut_corpse_day, 365, mdaily)
              call WritePoolFluxNcFile_corpse(filename_ncOut_corpse_day, pt(mp)%litterlayer_outputs, &
                                                  pt(mp)%soil_outputs(1), mp, linesWritten-364, linesWritten, mdaily, wrtYr)
           else
              linesWritten = pt(mp)%litterlayer_outputs%nYear          ! all pools and points have the same # of output years
              call InitPoolFluxNcFile_corpse(filename_ncOut_corpse, 1, mdaily)
              call WritePoolFluxNcFile_corpse(filename_ncOut_corpse, pt(mp)%litterlayer_outputs, &
                                              pt(mp)%soil_outputs(1), mp, linesWritten, linesWritten, mdaily, wrtYr)
           endif

        endif ! end isomModel==CORPSE

     endif !nloop==mloop

  enddo  ! end do nloop 

  !END OF SIMULATION CYCLE

! close(98)
!93    format(4(i6,',',2x),100(f15.6,',',2x))
!981   format(4(i6,',',2x),100(f15.6,',',2x))
END SUBROUTINE casacnpdriver

!----------------------------------------------------------------------------------------------------
! This subroutine is called to write file headers to .csv files when saving daily output 
! for a specific point.
!
SUBROUTINE WritePointFileHeaders(dirPtFile,mp)
  USE casavariable
  USE mimicsvariable
  USE corpsevariable
  implicit none
  character(len=100), intent(in) :: dirPtFile   !Directory where point file is to be written. If empty assum "./".
  integer, intent(in) :: mp                     ! Number of points in the simulation (indices 1..mp)

  !Local Variables
  character(len=10) ptstr
  integer :: npt

  write(ptstr, '(i10)') casafile%iptToSave
  ptstr = adjustl(ptstr)
  casafile%sPtFileNameCASA = trim(dirPtFile) // 'POINT_casa_' // trim(ptstr) // '.csv' 
  open(213,file=casafile%sPtFileNameCASA)
  !Create the header for the point file.  Append to the file later.
!!                         10        20        30        40        50        60        70        80        90        100
!!                1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
  write(213,703) 'npt,ijgcm,iYrCnt,idoy,iveg,tsoilavg,moistavg,casapool%Cplant(LEAF),casapool%Cplant(WOOD),', &
                 'casapool%Cplant(FROOT),casapool%Nplant(LEAF),casapool%Nplant(WOOD),casapool%Nplant(FROOT),', &
                 'casapool%Clitter(MET),casapool%Clitter(STR),casapool%Clitter(CWD),', &
                 'casapool%Nlitter(MET),casapool%Nlitter(STR),casapool%Nlitter(CWD),', &
                 'casapool%Csoil(MIC),casapool%Csoil(SLOW),casapool%Csoil(PASS),', & 
                 'casapool%Nsoil(MIC),casapool%Nsoil(SLOW),casapool%Nsoil(PASS),', & 
                 'casaflux%Crsoil,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd,', & 
                 'casaflux%Cgpp,casaflux%Cnpp,casaflux%Crmplant(LEAF),casaflux%Crmplant(WOOD),casaflux%Crmplant(FROOT),', &
                 'casaflux%Crgplant,phen%phase,casaflux%fracCalloc(LEAF),casaflux%fracCalloc(WOOD),casaflux%fracCalloc(FROOT),', &
                 'casapool%fT,casapool%fW,casaflux%CpassInpt,casapool%Nsoilmin,casaflux%Nminuptake,', &
                 'casaflux%Nlittermin,casaflux%Nsmin,casaflux%Nsimm,casaflux%Nminleach,casaflux%Nminloss,xkNlimiting,', &
                 'klitter(METB),klitter(STR),klitter(CWD),ksoil(MIC),ksoil(SLOW),ksoil(PASS),', &
                 'casaflux%NminDep,casaflux%NminFix,casaflux%NlitInptMet,casaflux%NlitInptStruc,', &
                 'casaflux%Cexudate,casaflux%Nexudate'
  close(213)
  703 format(a89,a90,a66,a66,a62,a62,a82,a101,a108,a81,a99,a75,a78,a35)

  ! The contents of this file are written each day in subroutine WritePointCASA which is called
  ! at the end of each day from subroutine biogeochem.

  if (isomModel == MIMICS) then
      iptToSave_mimics = casafile%iptToSave
      write(ptstr, '(i10)') iptToSave_mimics
      ptstr = adjustl(ptstr)
      sPtFileNameMIMICS = trim(dirPtFile) // 'POINT_mimics_' // trim(ptstr) // '.csv' 
      !Create the header for the point file.  Append to the file later.
      open(214,file=sPtFileNameMIMICS)

      if (icycle == 2) then
!!                             10        20        30        40        50        60        70        80        90        100
!!                    12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
      write(214,121) 'npt,ijgcm,iYrCnt,doy,tsoilavg,moistavg,hresp,inptMetC,inputStrC,LITm,LITs,MICr,MICk,SOMa,SOMc,SOMp,',   &
                     'dLITm,dLITs,dMICr,dMICk,dSOMa,dSOMc,dSOMp,', &
                     'LITmin(1),LITmin(2),LITmin(3),LITmin(4),',  &
                     'MICtrn(1),MICtrn(2),MICtrn(3),MICtrn(4),MICtrn(5),MICtrn(6),', &
                     'SOMmin(1),SOMmin(2),DEsorb,OXIDAT,', &
                     'fmet,tauMod,tauR,tauK,', &
                     'Vmax(R1),Vmax(R2),Vmax(R3),Vmax(K1),Vmax(K2),Vmax(K3),', &
                     'Km(R1),Km(R2),Km(R3),Km(K1),Km(K2),Km(K3),', &
                     'Vslope(R1),Vslope(R2),Vslope(R3),', &
                     'Vslope(K1),Vslope(K2),Vslope(K3),', &
                     'Kslope(R1),Kslope(R2),Kslope(R3),', &
                     'Kslope(K1),Kslope(K2),Kslope(K3),', &
                     'Vint,Kint,Tsoil,Tmoist,NPPan,CWD,ratLigN,', &
                     'cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd,', &
                     'fAVAL(1),fAVAL(2),fCHEM(1),fCHEM(2),fPHYS(1),fPHYS(2),', &
                     'Kmod(R1),Kmod(R2),Kmod(R3),Kmod(K1),Kmod(K2),Kmod(K3),cSOMpIn,', &

                     'inptMetN,inputStrN,LITmN,LITsN,MICrN,MICkN,SOMaN,SOMcN,SOMpN,DIN,', &
                     'dLITmN,dLITsN,dMICrN,dMICkN,dSOMaN,dSOMcN,dSOMpN,dDIN,', &
                     'LITminN(1),LITminN(2),LITminN(3),LITminN(4),', &
                     'MICtrnN(1),MICtrnN(2),MICtrnN(3),MICtrnN(4),MICtrnN(5),MICtrnN(6),', &
                     'SOMminN(1),SOMminN(2),DEsorbN,OXIDATN,', &
                     'DINup_r,DINup_k,upMICrC,upMICrN,upMICkC,upMICkN,', &
                     'Overflow_r,Overflow_k,Nspill_r,Nspill_k,Cbalance,Nbalance,', &
                     'nleaf2met,nleaf2str,nroot2met,nroot2str,nwd2str,nwood2cwd'

          ! The contents of this file are written each day in subroutine mimics_soil_reverseMM_CN

      else
 
!!                             10        20        30        40        50        60        70        80        90        100
!!                    12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
      write(214,122) 'npt,ijgcm,iYrCnt,doy,tsoilavg,moistavg,hresp,inptMetC,inputStrC,LITm,LITs,MICr,MICk,SOMa,SOMc,SOMp,',   &
                     'dLITm,dLITs,dMICr,dMICk,dSOMa,dSOMc,dSOMp,', &
                     'LITmin(1),LITmin(2),LITmin(3),LITmin(4),',  &
                     'MICtrn(1),MICtrn(2),MICtrn(3),MICtrn(4),MICtrn(5),MICtrn(6),', &
                     'SOMmin(1),SOMmin(2),DEsorb,OXIDAT,', &
                     'fmet,tauMod,tauR,tauK,', &
                     'Vmax(R1),Vmax(R2),Vmax(R3),Vmax(K1),Vmax(K2),Vmax(K3),', &
                     'Km(R1),Km(R2),Km(R3),Km(K1),Km(K2),Km(K3),', &
                     'Vslope(R1),Vslope(R2),Vslope(R3),', &
                     'Vslope(K1),Vslope(K2),Vslope(K3),', &
                     'Kslope(R1),Kslope(R2),Kslope(R3),', &
                     'Kslope(K1),Kslope(K2),Kslope(K3),', &
                     'Vint,Kint,Tsoil,Tmoist,NPPan,CWD,ratLigN,', &
                     'cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd,', &
                     'fAVAL(1),fAVAL(2),fCHEM(1),fCHEM(2),fPHYS(1),fPHYS(2),', &
                     'Kmod(R1),Kmod(R2),Kmod(R3),Kmod(K1),Kmod(K2),Kmod(K3),cSOMpIn,Cbalance,'

          ! The contents of this file are written each day in subroutine mimics_soil_forwardMM
          ! or mimics_soil_reverseMM.

      endif

      close(214)
      121 format(a99,a42,a40,a60,a34,a22,a54,a42,a33,a33,a33,a33,a41,a66,a54,a62, a65,a54,a44,a66,a38,a48,a58,a57)
      122 format(a99,a42,a40,a60,a34,a22,a54,a42,a33,a33,a33,a33,a41,a66,a54,a71)



  else if (isomModel == CORPSE) then
      iptToSave_corpse = casafile%iptToSave
      write(ptstr, '(i10)') iptToSave_corpse
      ptstr = adjustl(ptstr)
      sPtFileNameCORPSE = trim(dirPtFile) // 'POINT_corpse_' // trim(ptstr) // '.csv' 
      ! Subroutine WritePointCORPSE is called from subroutine casacnpdriver a the end of the
      ! simulation to write the header and contents to sPtFileNameCORPSE.
  endif

  ! Write to Cbalance.csv and Nbalance.csv. -mdh 3/2/2020
  open(220,file='Cbalance.csv')
  write(220,'(a62)') 'npt,ijgcm,iYrCnt,doy,Cplant,Clitter,Cmic,Csom,Clabile,Cin,Cout'
  close(220)
  if (icycle == 2) then
      open(221,file='Nbalance.csv')
      write(221,'(a65)') 'npt,ijgcm,iYrCnt,doy,Nplant,Nlitter,Nmic,Nsom,Nsoilminrl,Nin,Nout'
      close(221)
  endif


  ! Save the Index of the point to save 1..mp
! do npt = 1,mp
!     if (casamet%ijgcm(npt) .eq. casafile%iptToSave) then 
!         casafile%iptToSaveIndx = npt
!         exit
!     endif
! enddo
  casafile%iptToSaveIndx = casafile%iptToSave

END SUBROUTINE WritePointFileHeaders

!----------------------------------------------------------------------------------------------------
! This subroutine is called to write the contents to the daily point .csv file casafile%sPtFileNameCASA.
! It is called at the end of each day from subroutine biogeochem.
!
SUBROUTINE WritePointCASA(iYrCnt,idoy,mp,cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd)
  USE casaparm
  USE casavariable
  USE define_types
  USE phenvariable
  USE mimicsparam
  USE mimicsvariable

  implicit none
  integer, intent(in) :: iYrCnt, idoy, mp
  real(r_2), dimension(mp), intent(in)   :: cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd

  ! Local Variables
  integer :: npt
  real(r_2):: Cplant,Clitter,Cmic,Csom,Clabile,Cin,Cout
  real(r_2):: Nplant,Nlitter,Nmic,Nsom,Nminrl,Nin,Nout
  real(r_2):: unitConv

! if (casamet%ijgcm(npt) .eq. casafile%iptToSave)  
  npt = casafile%iptToSaveIndx
  open(213,file=casafile%sPtFileNameCASA,access='APPEND')
  write(213,701) npt,casamet%ijgcm(npt),iYrCnt,idoy,veg%iveg(npt),casamet%tsoilavg(npt),casamet%moistavg(npt), &
                 casapool%cplant(npt,LEAF),casapool%cplant(npt,WOOD),casapool%cplant(npt,FROOT), &
                 casapool%nplant(npt,LEAF),casapool%nplant(npt,WOOD),casapool%nplant(npt,FROOT), &
                 casapool%Clitter(npt,LEAF),casapool%Clitter(npt,WOOD),casapool%Clitter(npt,FROOT), &
                 casapool%Nlitter(npt,LEAF),casapool%Nlitter(npt,WOOD),casapool%Nlitter(npt,FROOT), &
                 casapool%Csoil(npt,MIC),casapool%Csoil(npt,SLOW),casapool%Csoil(npt,PASS), &
                 casapool%Nsoil(npt,MIC),casapool%Nsoil(npt,SLOW),casapool%Nsoil(npt,PASS), &
                 casaflux%Crsoil(npt),cleaf2met(npt),cleaf2str(npt),croot2met(npt),croot2str(npt), &
                 cwd2str(npt),cwd2co2(npt),cwood2cwd(npt), &
                 casaflux%Cgpp(npt),casaflux%Cnpp(npt), &
                 casaflux%Crmplant(npt,LEAF),casaflux%Crmplant(npt,WOOD),casaflux%Crmplant(npt,FROOT), &
                 casaflux%Crgplant(npt),real(phen%phase(npt)), &
                 casaflux%fracCalloc(npt,LEAF),casaflux%fracCalloc(npt,WOOD),casaflux%fracCalloc(npt,FROOT), &
                 casapool%fT(npt),casapool%fW(npt),casaflux%CpassInpt(npt),casapool%Nsoilmin(npt),casaflux%Nminuptake(npt), &
                 casaflux%Nlittermin(npt),casaflux%Nsmin(npt),casaflux%Nsimm(npt), &
                 casaflux%Nminleach(npt),casaflux%Nminloss(npt),casaflux%xkNlimiting(npt), &
                 casaflux%klitter(npt,METB),casaflux%klitter(npt,STR),casaflux%klitter(npt,CWD), &
                 casaflux%ksoil(npt,MIC),casaflux%ksoil(npt,SLOW),casaflux%ksoil(npt,PASS), &
                 casaflux%NminDep(npt),casaflux%NminFix(npt),casaflux%NlitInptMet(npt),casaflux%NlitInptStruc(npt), &
                 casaflux%Cexudate(npt),casaflux%Nexudate(npt)
  close(213)

  701 format(5(i6,','),54(f12.5,','),f12.8,4(',',f12.5),2(',',f12.8))

  ! Write to Cbalance.csv. All units are gC/m2 for pools and gC/m2/day for fluxes. -mdh 3/2/2020
  ! Write to Nbalance.csv. All units are gN/m2 for pools and gN/m2/day for fluxes. -mdh 3/2/2020


  if (isomModel == CASACNP) then

      open(220,file='Cbalance.csv',access='APPEND')
      Cplant = sum(casapool%Cplant(npt,:))
      Clitter = sum(casapool%Clitter(npt,:))
      Cmic = casapool%Csoil(npt,MIC)
      Csom = casapool%Csoil(npt,SLOW)+casapool%Csoil(npt,PASS)
      !Clabile comes from GPP, and was once part of leaf maintenance respiration, but that is commented out.
      Clabile = casapool%Clabile(npt)
      !Cin = casaflux%Cnpp(npt) 
      Cin = casaflux%Cgpp(npt) 
      Cout = casaflux%Crsoil(npt)+sum(casaflux%Crmplant(npt,:))+casaflux%Crgplant(npt)+casaflux%Clabloss(npt)
      write(220,702) npt,casamet%ijgcm(npt),iYrCnt,idoy,Cplant,Clitter,Cmic,Csom,Clabile,Cin,Cout
      close(220)
    
      if (icycle == 2) then
          open(221,file='Nbalance.csv',access='APPEND')
          Nplant = sum(casapool%Nplant(npt,:))
          Nlitter = sum(casapool%Nlitter(npt,:))
          Nmic = casapool%Nsoil(npt,MIC)
          Nsom = casapool%Nsoil(npt,SLOW)+casapool%Nsoil(npt,PASS)
          Nminrl = casapool%Nsoilmin(npt)
          Nin = casaflux%NminDep(npt) + casaflux%NminFix(npt)
          Nout = casaflux%Nminleach(npt)+casaflux%Nminloss(npt)
          write(221,703) npt,casamet%ijgcm(npt),iYrCnt,idoy,Nplant,Nlitter,Nmic,Nsom,Nminrl,Nin,Nout
          close(221)
      endif

  else if (isomModel == MIMICS) then

      ! (mg C/cm3)*(1/1000)(g/mg)*(10000)(cm2/m2)*depth(cm) = g C/m2
      unitConv = 10.0*mimicsbiome%depth(veg%iveg(npt))    ! Convert mgC/cm3 to gC/m2 by multipling by this factor

      open(220,file='Cbalance.csv',access='APPEND')
      Cplant = sum(casapool%Cplant(npt,:))
      Clitter = (mimicspool%LITm(npt) + mimicspool%LITs(npt))*unitConv
      Cmic = (mimicspool%MICr(npt) + mimicspool%MICk(npt))*unitConv
      Csom = (mimicspool%SOMa(npt) + mimicspool%SOMc(npt) + mimicspool%SOMp(npt))*unitConv
      !Clabile comes from GPP, and was once part of leaf maintenance respiration, but that is commented out.
      Clabile = casapool%Clabile(npt)
      !Cin = casaflux%Cnpp(npt) 
      Cin = casaflux%Cgpp(npt) 
      Cout = mimicsflux%Chresp(npt)*unitConv + sum(casaflux%Crmplant(npt,:))+casaflux%Crgplant(npt)+casaflux%Clabloss(npt)
      write(220,702) npt,casamet%ijgcm(npt),iYrCnt,idoy,Cplant,Clitter,Cmic,Csom,Clabile,Cin,Cout
      close(220)
    
      if (icycle == 2) then
          open(221,file='Nbalance.csv',access='APPEND')
          Nplant = sum(casapool%Nplant(npt,:))
          Nlitter = (mimicspool%LITmN(npt) + mimicspool%LITsN(npt))*unitConv
          Nmic = (mimicspool%MICrN(npt) + mimicspool%MICkN(npt))*unitConv
          Nsom = (mimicspool%SOMaN(npt) + mimicspool%SOMcN(npt) + mimicspool%SOMpN(npt))*unitConv
          Nminrl = casapool%Nsoilmin(npt)
          Nin = casaflux%NminDep(npt) + casaflux%NminFix(npt)
          Nout = casaflux%Nminleach(npt)+casaflux%Nminloss(npt)
          write(221,703) npt,casamet%ijgcm(npt),iYrCnt,idoy,Nplant,Nlitter,Nmic,Nsom,Nminrl,Nin,Nout
          close(221)
      endif

  endif

  702 format(4(i6,','),4(f16.8,','),3(f12.8,','))
  703 format(4(i6,','),4(f15.8,','),3(f12.8,','))

END SUBROUTINE WritePointCASA

!----------------------------------------------------------------------------------------------------
! Replace the last 4 characters in the filename (usally 'yyyy') with a 4-digit year

SUBROUTINE InsertYearInNcFileName(fileName_ncOut, year)

    character(len=100), intent(INOUT) :: filename_ncOut
    integer, INTENT(IN)               :: year

    !Local Variables
    integer :: idx
    character(len=4) syear
    logical :: back

    back = .true.

    write(syear,'(i4)') year    ! copy year to string syear
    if (year < 1000) then
        syear(1:1) = '0'
    endif
    if (year < 100) then
        syear(2:2) = '0'
    endif
    if (year < 10) then
        syear(3:3) = '0'
    endif

    idx = index(filename_ncOut, '.nc', back)    ! start position of '.nc' in file name (count starts at 1)

    if (idx == 0) then
        write(*,*) 'InsertYearInNcFileName: Error in NetCDF file name ', trim(filename_ncOut) 
        write (*,*)'The name does not end with .nc extension.'
        STOP
    endif

    filename_ncOut(idx-4:idx-1) = syear         ! replace 'yyyy' with syear

END SUBROUTINE InsertYearInNcFileName

!----------------------------------------------------------------------------------------------------
! Create a file name by inserting a year followed by the suffix '_daily' before the .nc extension

SUBROUTINE CreateDailyNcFileName(fileName_ncOut_day, fileName_ncOut, year)

    character(len=100), intent(INOUT) :: filename_ncOut_day
    character(len=100), intent(IN)    :: filename_ncOut
    integer, INTENT(IN)               :: year

    !Local Variables
    integer :: idx
    logical :: back

    back = .true.

    fileName_ncOut_day = fileName_ncOut
    !Don't insert year unless year >= 0
    if (year >= 0) then
        call InsertYearInNcFileName(fileName_ncOut_day, year)
    endif
    idx = index(fileName_ncOut_day, '.nc', back)    ! start position of '.nc' in file name (count starts at 1)
    if (idx == 0) then
        write(*,*) 'CreateDailyNcFileName: Error in NetCDF file name ', trim(fileName_ncOut_day) 
        write (*,*)'The name does not end with .nc extension.'
        STOP
    endif
    filename_ncOut_day(idx:idx+5) = '_daily'     
    filename_ncOut_day(idx+6:idx+8) = '.nc'

END SUBROUTINE CreateDailyNcFileName

!----------------------------------------------------------------------------------------------------
! Replace any non-space characers to the right of a ' ' or '!' with a ' ';
! these are assumed to be comments in the string

SUBROUTINE rmcomments(fileName)
    character(len=100), intent(INOUT) :: filename
   
    !Local variables
    integer :: idx, i

    filename = adjustl(filename)  ! remove leading spaces
    idx = index(filename, ' ')
    if (idx > 0) then
        do i = idx, len(filename)
            filename(i:i) = ' '
        enddo
    endif

    idx = index(filename, '!')
    if (idx > 0) then
        do i = idx, len(filename)
            filename(i:i) = ' '
        enddo
    endif

    !!write(*,*) 'Removed extra characters from ', trim(filename)

END SUBROUTINE rmcomments

