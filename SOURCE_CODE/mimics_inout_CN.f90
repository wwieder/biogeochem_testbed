!--------------------------------------------------------------------------------
! FILE: mimics_inout_CN.f90
!
! Purpose: 
!   Input/Output subroutines for the MIMICS model
!
!     SUBROUTINE mimics_readbiome - read pftlookup_mimics.csv (mimics parameters)
!     SUBROUTINE mimics_init - initialize mimics pools. If initcasa >= 1, read .csv file 
!     SUBROUTINE mimics_poolfluxout - write mimics C pools to restart .csv output file
!     SUBROUTINE mimics_poolfluxout_CN - write mimics C&N pools to restart .csv output file
!     SUBROUTINE WritePoolFluxNcFile_mimics_annual - write annual mimics pools and fluxes to netCDF file
!     SUBROUTINE WritePoolFluxNcFile_mimics_daily - write daily mimics pools and fluxes to netCDF file
!     SUBROUTINE WritePointMIMICS - write daily mimics pools and other quantities to a .csv file daily
!
! Contact: Melannie Hartman
!          melannie@ucar.edu
!
! History:
!   1/5/2015 - Created
!   11/7/2017 - Echo contents of parameter and restart files
!   12/3/2018 - Add inputs to SOMp to daily and annual NetCDF files, 
!               and to daily point output file
!   6/23/2019 - Add subroutine mimics_poolfluxout_CN
!   11/30/2019 - Add MIMICS N variables to netCDF output
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
  integer   :: ioErr
  integer   :: npt, nP, nv, nv1
  real(r_2) :: Pscalar, kmod_val(6)
  character(len=200) :: buffer 

  open(101,file=fname_mimicsbiome)

  !Skip past vegetation descriptions
  read(101,*) 
  do nv=1,mvtype
      read(101,*) 
  end do

  ! Fixed parameters
  ! TO DO - check input file variable name to assure corrert placement of parameter values
  write(*,*)
  write(*,*) "Reading MIMICS parameters from file ", trim(fname_mimicsbiome), "..."
  read(101,*) buffer
  write(*,*) trim(buffer)
 
  read(101,*) mimicsbiome%Vslope(R1) 
  read(101,*) mimicsbiome%Vslope(R2) 
  read(101,*) mimicsbiome%Vslope(R3)
  read(101,*) mimicsbiome%Vslope(K1)
  read(101,*) mimicsbiome%Vslope(K2)
  read(101,*) mimicsbiome%Vslope(K3) 

  write(*,'(2x,a13,2x,f10.6)') 'Vslope(R1)=', mimicsbiome%Vslope(R1) 
  write(*,'(2x,a13,2x,f10.6)') 'Vslope(R2)=', mimicsbiome%Vslope(R2) 
  write(*,'(2x,a13,2x,f10.6)') 'Vslope(R3)=', mimicsbiome%Vslope(R3)
  write(*,'(2x,a13,2x,f10.6)') 'Vslope(K1)=', mimicsbiome%Vslope(K1)
  write(*,'(2x,a13,2x,f10.6)') 'Vslope(K2)=', mimicsbiome%Vslope(K2)
  write(*,'(2x,a13,2x,f10.6)') 'Vslope(K3)=', mimicsbiome%Vslope(K3) 


  read(101,*) mimicsbiome%Vint(R1) 
  read(101,*) mimicsbiome%Vint(R2) 
  read(101,*) mimicsbiome%Vint(R3)
  read(101,*) mimicsbiome%Vint(K1)
  read(101,*) mimicsbiome%Vint(K2)
  read(101,*) mimicsbiome%Vint(K3) 

  write(*,'(2x,a13,2x,f10.6)') 'Vint(R1)=', mimicsbiome%Vint(R1) 
  write(*,'(2x,a13,2x,f10.6)') 'Vint(R2)=', mimicsbiome%Vint(R2) 
  write(*,'(2x,a13,2x,f10.6)') 'Vint(R3)=', mimicsbiome%Vint(R3)
  write(*,'(2x,a13,2x,f10.6)') 'Vint(K1)=', mimicsbiome%Vint(K1)
  write(*,'(2x,a13,2x,f10.6)') 'Vint(K2)=', mimicsbiome%Vint(K2)
  write(*,'(2x,a13,2x,f10.6)') 'Vint(K3)=', mimicsbiome%Vint(K3)


  read(101,*) mimicsbiome%av(R1) 
  read(101,*) mimicsbiome%av(R2) 
  read(101,*) mimicsbiome%av(R3)
  read(101,*) mimicsbiome%av(K1)
  read(101,*) mimicsbiome%av(K2)
  read(101,*) mimicsbiome%av(K3) 

  write(*,'(2x,a13,2x,f12.8)') 'av(R1)=', mimicsbiome%av(R1) 
  write(*,'(2x,a13,2x,f12.8)') 'av(R2)=', mimicsbiome%av(R2) 
  write(*,'(2x,a13,2x,f12.8)') 'av(R3)=', mimicsbiome%av(R3)
  write(*,'(2x,a13,2x,f12.8)') 'av(K1)=', mimicsbiome%av(K1)
  write(*,'(2x,a13,2x,f12.8)') 'av(K2)=', mimicsbiome%av(K2)
  write(*,'(2x,a13,2x,f12.8)') 'av(K3)=', mimicsbiome%av(K3) 


  read(101,*) mimicsbiome%Kslope(R1) 
  read(101,*) mimicsbiome%Kslope(R2) 
  read(101,*) mimicsbiome%Kslope(R3)  
  read(101,*) mimicsbiome%Kslope(K1)  
  read(101,*) mimicsbiome%Kslope(K2) 
  read(101,*) mimicsbiome%Kslope(K3) 

  write(*,'(2x,a13,2x,f10.6)') 'Kslope(R1)=', mimicsbiome%Kslope(R1) 
  write(*,'(2x,a13,2x,f10.6)') 'Kslope(R2)=', mimicsbiome%Kslope(R2) 
  write(*,'(2x,a13,2x,f10.6)') 'Kslope(R3)=', mimicsbiome%Kslope(R3)  
  write(*,'(2x,a13,2x,f10.6)') 'Kslope(K1)=', mimicsbiome%Kslope(K1)  
  write(*,'(2x,a13,2x,f10.6)') 'Kslope(K2)=', mimicsbiome%Kslope(K2) 
  write(*,'(2x,a13,2x,f10.6)') 'Kslope(K3)=', mimicsbiome%Kslope(K3) 


  read(101,*) mimicsbiome%Kint(R1) 
  read(101,*) mimicsbiome%Kint(R2) 
  read(101,*) mimicsbiome%Kint(R3)  
  read(101,*) mimicsbiome%Kint(K1)  
  read(101,*) mimicsbiome%Kint(K2) 
  read(101,*) mimicsbiome%Kint(K3)

  write(*,'(2x,a13,2x,f10.6)') 'Kint(R1)=', mimicsbiome%Kint(R1) 
  write(*,'(2x,a13,2x,f10.6)') 'Kint(R2)=', mimicsbiome%Kint(R2) 
  write(*,'(2x,a13,2x,f10.6)') 'Kint(R3)=', mimicsbiome%Kint(R3)  
  write(*,'(2x,a13,2x,f10.6)') 'Kint(K1)=', mimicsbiome%Kint(K1)  
  write(*,'(2x,a13,2x,f10.6)') 'Kint(K2)=', mimicsbiome%Kint(K2) 
  write(*,'(2x,a13,2x,f10.6)') 'Kint(K3)=', mimicsbiome%Kint(K3)


  read(101,*) mimicsbiome%ak(R1) 
  read(101,*) mimicsbiome%ak(R2) 
  read(101,*) mimicsbiome%ak(R3)
  read(101,*) mimicsbiome%ak(K1)
  read(101,*) mimicsbiome%ak(K2)
  read(101,*) mimicsbiome%ak(K3) 

  write(*,'(2x,a13,2x,f10.6)') 'ak(R1)=', mimicsbiome%ak(R1) 
  write(*,'(2x,a13,2x,f10.6)') 'ak(R2)=', mimicsbiome%ak(R2) 
  write(*,'(2x,a13,2x,f10.6)') 'ak(R3)=', mimicsbiome%ak(R3)
  write(*,'(2x,a13,2x,f10.6)') 'ak(K1)=', mimicsbiome%ak(K1)
  write(*,'(2x,a13,2x,f10.6)') 'ak(K2)=', mimicsbiome%ak(K2)
  write(*,'(2x,a13,2x,f10.6)') 'ak(K3)=', mimicsbiome%ak(K3) 


  read(101,*) mimicsbiome%Vmod(R1)
  read(101,*) mimicsbiome%Vmod(R2)
  read(101,*) mimicsbiome%Vmod(R3)
  read(101,*) mimicsbiome%Vmod(K1)
  read(101,*) mimicsbiome%Vmod(K2)
  read(101,*) mimicsbiome%Vmod(K3)

  write(*,'(2x,a13,2x,f10.6)') 'Vmod(R1)=', mimicsbiome%Vmod(R1)
  write(*,'(2x,a13,2x,f10.6)') 'Vmod(R2)=', mimicsbiome%Vmod(R2)
  write(*,'(2x,a13,2x,f10.6)') 'Vmod(R3)=', mimicsbiome%Vmod(R3)
  write(*,'(2x,a13,2x,f10.6)') 'Vmod(K1)=', mimicsbiome%Vmod(K1)
  write(*,'(2x,a13,2x,f10.6)') 'Vmod(K2)=', mimicsbiome%Vmod(K2)
  write(*,'(2x,a13,2x,f10.6)') 'Vmod(K3)=', mimicsbiome%Vmod(K3)


  read(101,*) kmod_val(R1)
  read(101,*) kmod_val(R2)
  read(101,*) kmod_val(R3)
  read(101,*) kmod_val(K1)
  read(101,*) kmod_val(K2)
  read(101,*) kmod_val(K3)

  write(*,'(2x,a13,2x,f10.6)') 'kmod_val(R1)=', kmod_val(R1)
  write(*,'(2x,a13,2x,f10.6)') 'kmod_val(R2)=', kmod_val(R2)
  write(*,'(2x,a13,2x,f10.6)') 'kmod_val(R3)=', kmod_val(R3)
  write(*,'(2x,a13,2x,f10.6)') 'kmod_val(K1)=', kmod_val(K1)
  write(*,'(2x,a13,2x,f10.6)') 'kmod_val(K2)=', kmod_val(K2)
  write(*,'(2x,a13,2x,f10.6)') 'kmod_val(K3)=', kmod_val(K3)


  read(101,*) mimicsbiome%KO(1)
  read(101,*) mimicsbiome%KO(2)

  write(*,'(2x,a13,2x,f10.6)') 'KO(1)=', mimicsbiome%KO(1)
  write(*,'(2x,a13,2x,f10.6)') 'KO(2)=', mimicsbiome%KO(2)


  read(101,*) mimicsbiome%MGE(1)
  read(101,*) mimicsbiome%MGE(2)
  read(101,*) mimicsbiome%MGE(3)
  read(101,*) mimicsbiome%MGE(4)

  write(*,'(2x,a13,2x,f10.6)') 'MGE(1)=', mimicsbiome%MGE(1)
  write(*,'(2x,a13,2x,f10.6)') 'MGE(2)=', mimicsbiome%MGE(2)
  write(*,'(2x,a13,2x,f10.6)') 'MGE(3)=', mimicsbiome%MGE(3)
  write(*,'(2x,a13,2x,f10.6)') 'MGE(4)=', mimicsbiome%MGE(4)


  read(101,*) mimicsbiome%tau_r(1)
  read(101,*) mimicsbiome%tau_r(2)

  write(*,'(2x,a13,2x,f10.6)') 'tau_r(1)=', mimicsbiome%tau_r(1)
  write(*,'(2x,a13,2x,f10.6)') 'tau_r(2)=', mimicsbiome%tau_r(2)


  read(101,*) mimicsbiome%tau_k(1)
  read(101,*) mimicsbiome%tau_k(2)

  write(*,'(2x,a13,2x,f10.6)') 'tau_k(1)=', mimicsbiome%tau_k(1)
  write(*,'(2x,a13,2x,f10.6)') 'tau_k(2)=', mimicsbiome%tau_k(2)


  read(101,*) mimicsbiome%tauModDenom
  read(101,*) mimicsbiome%tauMod_MIN
  read(101,*) mimicsbiome%tauMod_MAX

  write(*,'(2x,a13,2x,f10.6)') 'tauModDenom=', mimicsbiome%tauModDenom
  write(*,'(2x,a13,2x,f10.6)') 'tauMod_MIN=', mimicsbiome%tauMod_MIN
  write(*,'(2x,a13,2x,f10.6)') 'tauMod_MAX=', mimicsbiome%tauMod_MAX


  read(101,*) mimicsbiome%fPHYS_r(1)
  read(101,*) mimicsbiome%fPHYS_r(2)

  write(*,'(2x,a13,2x,f10.6)') 'fPHYS_r(1)=', mimicsbiome%fPHYS_r(1)
  write(*,'(2x,a13,2x,f10.6)') 'fPHYS_r(2)=', mimicsbiome%fPHYS_r(2)


  read(101,*) mimicsbiome%fPHYS_K(1)
  read(101,*) mimicsbiome%fPHYS_K(2)

  write(*,'(2x,a13,2x,f10.6)') 'fPHYS_K(1)=', mimicsbiome%fPHYS_K(1)
  write(*,'(2x,a13,2x,f10.6)') 'fPHYS_K(2)=', mimicsbiome%fPHYS_K(2)


  read(101,*) mimicsbiome%fCHEM_r(1)
  read(101,*) mimicsbiome%fCHEM_r(2)
  read(101,*) mimicsbiome%fCHEM_r(3)

  write(*,'(2x,a13,2x,f10.6)') 'fCHEM_r(1)=', mimicsbiome%fCHEM_r(1)
  write(*,'(2x,a13,2x,f10.6)') 'fCHEM_r(2)=', mimicsbiome%fCHEM_r(2)
  write(*,'(2x,a13,2x,f10.6)') 'fCHEM_r(3)=', mimicsbiome%fCHEM_r(3)


  read(101,*) mimicsbiome%fCHEM_K(1)
  read(101,*) mimicsbiome%fCHEM_K(2)
  read(101,*) mimicsbiome%fCHEM_K(3)

  write(*,'(2x,a13,2x,f10.6)') 'fCHEM_K(1)=', mimicsbiome%fCHEM_K(1)
  write(*,'(2x,a13,2x,f10.6)') 'fCHEM_K(2)=', mimicsbiome%fCHEM_K(2)
  write(*,'(2x,a13,2x,f10.6)') 'fCHEM_K(3)=', mimicsbiome%fCHEM_K(3)


  read(101,*) mimicsbiome%fSOM_p(1)
  read(101,*) mimicsbiome%fSOM_p(2)

  write(*,'(2x,a13,2x,f10.6)') 'fSOM_p(1)=', mimicsbiome%fSOM_p(1)
  write(*,'(2x,a13,2x,f10.6)') 'fSOM_p(2)=', mimicsbiome%fSOM_p(2)


  read(101,*) mimicsbiome%phys_scalar(1)
  read(101,*) mimicsbiome%phys_scalar(2)

  write(*,'(2x,a15,2x,f10.6)') 'phys_scalar(1)=', mimicsbiome%phys_scalar(1)
  write(*,'(2x,a15,2x,f10.6)') 'phys_scalar(2)=', mimicsbiome%phys_scalar(2)


  read(101,*) mimicsbiome%Fi(metbc)
  read(101,*) mimicsbiome%Fi(struc)

  write(*,'(2x,a13,2x,f10.6)') 'Fi(metbc)=', mimicsbiome%Fi(metbc)
  write(*,'(2x,a13,2x,f10.6)') 'Fi(struc)=',mimicsbiome%Fi(struc)


  read(101,*) mimicsbiome%fmet_p(1)
  read(101,*) mimicsbiome%fmet_p(2)
  read(101,*) mimicsbiome%fmet_p(3)

  write(*,'(2x,a13,2x,f10.6)') 'fmet_p(1)=', mimicsbiome%fmet_p(1)
  write(*,'(2x,a13,2x,f10.6)') 'fmet_p(2)=', mimicsbiome%fmet_p(2)
  write(*,'(2x,a13,2x,f10.6)') 'fmet_p(3)=', mimicsbiome%fmet_p(3)

  !! ----------------------------------------------------------------------
  !! Additional N-related parameters for MIMICS-CN. -mdh 6/21/2019
 
  read(101,*) mimicsbiome%densDep
  read(101,*)
  read(101,*)
  !read(101,*) mimicsbiome%fracNimport_r
  !read(101,*) mimicsbiome%fracNimport_k

  write(*,'(2x,a8,2x,f10.6)')  'densDep=', mimicsbiome%densDep
  !write(*,'(2x,a14,2x,f10.6)') 'fracNimport_r=', mimicsbiome%fracNimport_r
  !write(*,'(2x,a14,2x,f10.6)') 'fracNimport_k=', mimicsbiome%fracNimport_k

  read(101,*) mimicsbiome%NUE(1)
  read(101,*) mimicsbiome%NUE(2)
  read(101,*) mimicsbiome%NUE(3)
  read(101,*) mimicsbiome%NUE(4)

  write(*,'(2x,a7,2x,f10.6)') 'NUE(1)=', mimicsbiome%NUE(1)
  write(*,'(2x,a7,2x,f10.6)') 'NUE(2)=', mimicsbiome%NUE(2)
  write(*,'(2x,a7,2x,f10.6)') 'NUE(3)=', mimicsbiome%NUE(3)
  write(*,'(2x,a7,2x,f10.6)') 'NUE(4)=', mimicsbiome%NUE(4)

  read(101,*) mimicsbiome%CN_r
  read(101,*) mimicsbiome%CN_k
  read(101,*) mimicsbiome%fracDINavailMIC

  write(*,'(2x,a5,2x,f10.6)')  'CN_r=', mimicsbiome%CN_r
  write(*,'(2x,a5,2x,f10.6)')  'CN_k=', mimicsbiome%CN_k
  write(*,'(2x,a16,2x,f10.6)') 'fracDINavailMIC=', mimicsbiome%fracDINavailMIC
  !! ----------------------------------------------------------------------


  !Biome-specific parameters
  read(101,*) buffer
  write(*,*) trim(buffer)
  read(101,*) buffer
  write(*,*) trim(buffer)
  do nv=1,mvtype
     read(101,*) nv1,mimicsbiome%depth(nv)
     write(*,'(2x,i2,2x,f7.2)') nv1,mimicsbiome%depth(nv)
     if (mimicsbiome%depth(nv) <= 0.0) then
        write(*,*) 'Error in ', trim(filename_mimicsbiome),' (depth <= 0): depth(',nv,')=',mimicsbiome%depth(nv)
        STOP
     endif
  end do

  ! The fWFunction switch determines what water function to use (CASACNP=1,MIMICS=2,CORPSE=3)
  ! If this parameter was not inserted (unexpected EOF), then default to MIMICS f(W)=1.0
  read(101,*,IOSTAT=ioErr) mimicsbiome%fWFunction
  if (ioErr .eq. 0) THEN
    write(*,'(2x,a11,2x,i2)') 'fWFunction=', mimicsbiome%fWFunction
    if (mimicsbiome%fWFunction .eq. CASACNP) then
      write(*,*) '  Using f(W) function from CASACNP...'
    elseif (mimicsbiome%fWFunction .eq. CORPSE) then
      write(*,*) '  Using f(W) function from CORPSE...'
    else
      write(*,*) '  Using f(W)=1.0 from MIMICS...'
      mimicsbiome%fWFunction = MIMICS
    endif
  else
    ! Error or EOF. Default to MIMICS f(W)=1.0
    write(*,*) '  Using f(W)=1.0 from MIMICS...'
    mimicsbiome%fWFunction = MIMICS
  endif


  close(101)
  write(*,*) "Done reading MIMICS parameters from file ", trim(fname_mimicsbiome), "..."

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
          ! When icycle > 1, mimicsbiome%ligninNratio will be reset in mimics_coeffplant. 
          ! The value of mimicsbiome%ligninNratio below is used for C-only simulations. -mdh 1/20/2020
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
! This subroutine should be called after subroutine init_casa (not in place of)
! and after subroutine mimics_readbiome.

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
  !integer,           INTENT(IN) :: icycle       ! 1 = C only, 2 = C+N
  
  !Local Variables
  integer   :: np,npt,npz,nl,ns,nland,nlandz
  real(r_2) :: nyearz,ivtz,latz,lonz,areacellz
  character(len=200) :: buffer

  mimicspool%LITm(:) = 0.0
  mimicspool%LITs(:) = 0.0 
  mimicspool%MICr(:) = 0.0 
  mimicspool%MICk(:) = 0.0 
  mimicspool%SOMa(:) = 0.0 
  mimicspool%SOMc(:) = 0.0 
  mimicspool%SOMp(:) = 0.0 

  mimicspool%LITmN(:) = 0.0
  mimicspool%LITsN(:) = 0.0 
  mimicspool%MICrN(:) = 0.0 
  mimicspool%MICkN(:) = 0.0 
  mimicspool%SOMaN(:) = 0.0 
  mimicspool%SOMcN(:) = 0.0 
  mimicspool%SOMpN(:) = 0.0 

  WHERE(casamet%iveg2 /= icewater)
      mimicspool%LITm(:) = 1.0
      mimicspool%LITs(:) = 1.0
      mimicspool%MICr(:) = 0.015
      mimicspool%MICk(:) = 0.025
      mimicspool%SOMa(:) = 1.0
      mimicspool%SOMc(:) = 1.0
      mimicspool%SOMp(:) = 1.0

      mimicspool%LITmN(:) = mimicspool%LITm(:)/10.0
      mimicspool%LITsN(:) = mimicspool%LITs(:)/10.0
      mimicspool%MICrN(:) = mimicspool%MICr(:)/mimicsbiome%CN_r
      mimicspool%MICkN(:) = mimicspool%MICk(:)/mimicsbiome%CN_k
      mimicspool%SOMaN(:) = mimicspool%SOMa(:)/10.0
      mimicspool%SOMcN(:) = mimicspool%SOMc(:)/10.0
      mimicspool%SOMpN(:) = mimicspool%SOMp(:)/10.0
  END WHERE

  !If not a spinup run (initcasa .ne. 0) read initial pool values from file,
  !overwriting the C&N pool assignments above. 
  if (initcasa >= 1) then
      write(*,*)
      write(*,*) "Reading initial MIMICS pool file: ", trim(filename_mimicsipool), "..."
      open(105,file=filename_mimicsipool)
      read(105,*)  buffer ! Skip past header line
      !write(*,*) trim(buffer)
      do npt =1, mp
          if (icycle == 1) then
              read(105,*) nyearz,npz,ivtz,latz,lonz,areacellz,         &
                      mimicspool%LITm(npt), mimicspool%LITs(npt),  &
                      mimicspool%MICr(npt), mimicspool%MICk(npt),  &
                      mimicspool%SOMa(npt), mimicspool%SOMc(npt), mimicspool%SOMp(npt)
          else
              read(105,*) nyearz,npz,ivtz,latz,lonz,areacellz,         &
                      mimicspool%LITm(npt), mimicspool%LITs(npt),  &
                      mimicspool%MICr(npt), mimicspool%MICk(npt),  &
                      mimicspool%SOMa(npt), mimicspool%SOMc(npt), mimicspool%SOMp(npt), &
                      mimicspool%LITmN(npt), mimicspool%LITsN(npt),  &
                      mimicspool%MICrN(npt), mimicspool%MICkN(npt),  &
                      mimicspool%SOMaN(npt), mimicspool%SOMcN(npt), mimicspool%SOMpN(npt)
          endif

!         write(*,123) nyearz,npz,ivtz,latz,lonz,areacellz,        &
!                     mimicspool%LITm(npt), mimicspool%LITs(npt),  &
!                     mimicspool%MICr(npt), mimicspool%MICk(npt),  &
!                     mimicspool%SOMa(npt), mimicspool%SOMc(npt), mimicspool%SOMp(npt), &
!                     mimicspool%LITmN(npt), mimicspool%LITsN(npt),  &
!                     mimicspool%MICrN(npt), mimicspool%MICkN(npt),  &
!                     mimicspool%SOMaN(npt), mimicspool%SOMcN(npt), mimicspool%SOMpN(npt)
!
! 123 format(f4.0,',',i4,',',f4.0,',',17(f8.4,','))

          !ATTENTION: Check npz, ivtz, latz, lonz, areacellz against values read by casa_init
          !TO DO
      end do
      close(105)
      print *, "Done reading initial MIMICS pool file: ", filename_mimicsipool, "..."
  endif

  ! check for negative pool sizes
  mimicspool%LITm = max(0.0, mimicspool%LITm)
  mimicspool%LITs = max(0.0, mimicspool%LITs)
  mimicspool%MICr = max(0.0, mimicspool%MICr)
  mimicspool%MICk = max(0.0, mimicspool%MICk)
  mimicspool%SOMa = max(0.0, mimicspool%SOMa)
  mimicspool%SOMc = max(0.0, mimicspool%SOMc)
  mimicspool%SOMp = max(0.0, mimicspool%SOMp)

  mimicspool%LITmN = max(0.0, mimicspool%LITmN)
  mimicspool%LITsN = max(0.0, mimicspool%LITsN)
  mimicspool%MICrN = max(0.0, mimicspool%MICrN)
  mimicspool%MICkN = max(0.0, mimicspool%MICkN)
  mimicspool%SOMaN = max(0.0, mimicspool%SOMaN)
  mimicspool%SOMcN = max(0.0, mimicspool%SOMcN)
  mimicspool%SOMpN = max(0.0, mimicspool%SOMpN)

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
! casapool%Nlitter(:,:) = 0.0
! casapool%Nsoilmin(:)  = 0.0
! casapool%Psoil(:,:)   = 0.0
! casapool%Plitter(:,:) = 0.0

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
  mimicspoolAn%thetaLiqAn(:)   = mimicspoolAn%thetaLiqAn(:)   * xyear
  mimicspoolAn%thetaFrznAn(:)  = mimicspoolAn%thetaFrznAn(:)  * xyear
  mimicspoolAn%fTAn(:)         = mimicspoolAn%fTAn(:)         * xyear
  mimicspoolAn%fWAn(:)         = mimicspoolAn%fWAn(:)         * xyear

  !Divide by number of simulation years to get average annual flux
  mimicsfluxAn%ChrespAn(:)      = mimicsfluxAn%ChrespAn(:) * xyear2    
  mimicsfluxAn%CLitInputAn(:,:) = mimicsfluxAn%CLitInputAn(:,:) * xyear2   
  mimicsfluxAn%CSOMpInputAn(:)  = mimicsfluxAn%CSOMpInputAn(:) * xyear2    

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

!--------------------------------------------------------------------------------

SUBROUTINE mimics_poolfluxout_CN(filename_mimicsepool,mp,iYrCnt,myear,writeToRestartCSVfile)

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
  mimicspoolAn%NlitterAn(:,:)  = mimicspoolAn%NlitterAn(:,:)  * xyear
  mimicspoolAn%NmicrobeAn(:,:) = mimicspoolAn%NmicrobeAn(:,:) * xyear
  mimicspoolAn%NsoilAn(:,:)    = mimicspoolAn%NsoilAn(:,:)    * xyear
  mimicspoolAn%DINAn(:)        = mimicspoolAn%DINAn(:)        * xyear
  mimicspoolAn%thetaLiqAn(:)   = mimicspoolAn%thetaLiqAn(:)   * xyear
  mimicspoolAn%thetaFrznAn(:)  = mimicspoolAn%thetaFrznAn(:)  * xyear
  mimicspoolAn%fTAn(:)         = mimicspoolAn%fTAn(:)         * xyear
  mimicspoolAn%fWAn(:)         = mimicspoolAn%fWAn(:)         * xyear

  !Divide by number of simulation years to get average annual flux
  mimicsfluxAn%ChrespAn(:)      = mimicsfluxAn%ChrespAn(:) * xyear2    
  mimicsfluxAn%CLitInputAn(:,:) = mimicsfluxAn%CLitInputAn(:,:) * xyear2   
  mimicsfluxAn%CSOMpInputAn(:)  = mimicsfluxAn%CSOMpInputAn(:) * xyear2    
  mimicsfluxAn%NLitInputAn(:,:) = mimicsfluxAn%NLitInputAn(:,:) * xyear2   

!--------------------------------------------------------------------------

  if (writeToRestartCSVfile) then

      nout=106
      open(nout,file=filename_mimicsepool)

      ! mimicsbal%sumcbal=min(9999.0,max(-9999.0,mimicsbal%sumcbal))

      !Write header line
      write(nout,91) "iYrCnt,npt,veg%iveg,",               &
                     "casamet%lat,casamet%lon,",           &
                     "casamet%areacell,",                  &
                     "mimicspool%LITm,mimicspool%LITs,",   &
                     "mimicspool%MICr,mimicspool%MICk,",   &
                     "mimicspool%SOMa,mimicspool%SOMc,",   &
                     "mimicspool%SOMp,",                   &
                     "mimicspool%LITmN,mimicspool%LITsN,", &
                     "mimicspool%MICrN,mimicspool%MICkN,", &
                     "mimicspool%SOMaN,mimicspool%SOMcN,", &
                     "mimicspool%SOMpN"
    91 format(a20,a24,a17,a32,a32,a32,a16,a34,a34,a34,a16)

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
                         mimicspool%SOMa(npt), mimicspool%SOMc(npt), mimicspool%SOMp(npt), &
                         mimicspool%LITmN(npt), mimicspool%LITsN(npt),  &
                         mimicspool%MICrN(npt), mimicspool%MICkN(npt),  &
                         mimicspool%SOMaN(npt), mimicspool%SOMcN(npt), mimicspool%SOMpN(npt)
      end do

      CLOSE(nout)

  endif

92  format(3(i6,',',1x),17(f18.10,',',1x))

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

END SUBROUTINE mimics_poolfluxout_CN

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
      character(len=*), intent(in) :: filename_ncOut    ! NetCDF output file name 
      integer, intent(in) :: mp                         ! number of points with valid data
      integer, intent(in) :: year                       ! output year
!   LOCAL VARIABLES
      integer :: i
      integer :: ncid                           ! netcdf file ID
      integer :: status                         ! function return status
!     integer :: dimid_mp                       ! netcdf dimension id
      integer :: dimid_lat                      ! netcdf dimension id
      integer :: dimid_lon                      ! netcdf dimension id
      integer :: dimid_time                     ! netcdf dimension id
      integer :: nlon, nlat, ntimes             ! Dimension sizes for NetCDf file
      integer :: dims(3)                        ! Array of NetCDF dimension IDs for defining variables
      integer :: start1(1), count1(1)           ! start and count arrays for writing 1-D data from netcdf files
      integer :: start2(2), count2(2)           ! start and count arrays for writing 2-D data from netcdf files
      integer :: start3(3), count3(3)           ! start and count arrays for writing 3-D data from netcdf files
      integer :: varid_lon, varid_lat           ! NetCDF variable ID for latitude and longitude
      integer :: varid_time                     ! NetCDF variable ID for time
!     integer :: varid_year                     ! NetCDF variable ID for year
      integer :: varid_mask                     ! NetCDF variable ID for cellMissing(nlon,nlat)
      integer :: varid_cellid                   ! NetCDF variable ID for cellid(nlon,nlat)
      integer :: varid_igbp                     ! NetCDF variable ID for IGBP_PFT(nlon,nlat)
      integer :: varid_landarea                 ! NetCDF variable ID for land area
      integer :: varid_cLITm                    ! NetCDF variable ID for metablic litter C
      integer :: varid_cLITs                    ! NetCDF variable ID for structural litter C
      integer :: varid_cMICr                    ! NetCDF variable ID for microbe r-selected pool C
      integer :: varid_cMICk                    ! NetCDF variable ID for microbe K-selected pool C
      integer :: varid_cSOMa                    ! NetCDF variable ID for available soil C
      integer :: varid_cSOMc                    ! NetCDF variable ID for chemically protected soil C
      integer :: varid_cSOMp                    ! NetCDF variable ID for physically protected soil C
      integer :: varid_cHresp                   ! NetCDF variable ID for soil heterotrophic respiration C
      integer :: varid_cSOMpIn                  ! NetCDF variable ID for physically protected soil C inputs
      integer :: varid_cLitIn_m                 ! NetCDF variable ID for metabolic litter Input C
      integer :: varid_cLitIn_s                 ! NetCDF variable ID for structural litter Input C

      integer :: varid_nLITm                    ! NetCDF variable ID for metablic litter N
      integer :: varid_nLITs                    ! NetCDF variable ID for structural litter N
      integer :: varid_nMICr                    ! NetCDF variable ID for microbe r-selected pool N
      integer :: varid_nMICk                    ! NetCDF variable ID for microbe K-selected pool N
      integer :: varid_nSOMa                    ! NetCDF variable ID for available soil N
      integer :: varid_nSOMc                    ! NetCDF variable ID for chemically protected soil N
      integer :: varid_nSOMp                    ! NetCDF variable ID for physically protected soil N
      integer :: varid_DIN                      ! NetCDF variable ID for dissolved inorganic N
      integer :: varid_nLitIn_m                 ! NetCDF variable ID for metabolic litter Input N
      integer :: varid_nLitIn_s                 ! NetCDF variable ID for structural litter Input N

      integer :: varid_thetaLiq                 ! NetCDF variable ID for fraction of liquid soil water saturation
      integer :: varid_thetaFrzn                ! NetCDF variable ID for fraction of frozen soil water saturation
!     integer :: varid_fT                       ! NetCDF variable ID for soil temperature amultiplier on decompostion
      integer :: varid_fW                       ! NetCDF variable ID for soil moisture multiplier on decompostion
      character*100 :: attr_name                ! String for assigning global and variable attributes
      character*100 :: attr_units               ! String for assigning global and variable attributes
      character*10 :: date_string               ! String for assigning date to global attributes
      character*8 :: time_string                ! String for assigning time to global attributes
      integer :: verbose=0
      integer :: npt, ilon, ilat, itime
      integer, allocatable :: IGBP_PFT(:,:)     ! IGBP_PFT(nlon,nlat) IGBP PFT classification (1-18)
      real(4), allocatable :: landarea(:,:)     ! landarea(nlon,nlat) km^2
      real(4), allocatable :: var1(:,:,:)       ! gridded output variable
      real(4), allocatable :: var2(:,:,:)       ! gridded output variable
      real(4), allocatable :: var3(:,:,:)       ! gridded output variable
      real(4), allocatable :: var4(:,:,:)       ! gridded output variable
      real(4), allocatable :: var5(:,:,:)       ! gridded output variable
      real(4), allocatable :: var6(:,:,:)       ! gridded output variable
      real(4), allocatable :: var7(:,:,:)       ! gridded output variable
      real(4), allocatable :: var8(:,:,:)       ! gridded output variable
      real(4), allocatable :: var9(:,:,:)       ! gridded output variable
      real(4), allocatable :: var10(:,:,:)      ! gridded output variable
      real(4), allocatable :: var11(:,:,:)      ! gridded output variable
      real(4), allocatable :: var12(:,:,:)      ! gridded output variable
      real(4), allocatable :: var13(:,:,:)      ! gridded output variable
      real(4), allocatable :: var14(:,:,:)      ! gridded output variable
      real(4), allocatable :: var15(:,:,:)      ! gridded output variable
      real(4) :: time                           ! time in years (1..ntimes)

      if (verbose .ge. 0) print *, "Writing output to file ", trim(filename_ncOut), "..."

! Output these variables in NetCDF file:
!     veg%iveg(npt)                             - IGBP PFTs
!     casamet%lat(npt)                          - latitudes
!     casamet%lon(npt)                          - longitudes
!     casamet%areacell(npt)*(1.0e-6)            - landarea (km^2)
!   Average annual values for fluxes and pools      
!     mimicsfluxAn%CLitInputAn(npt,METBC)       - metabolic litter inputs (gC/m2/yr)
!     mimicsfluxAn%CLitInputAn(npt,STRUC)       - structural litter inputs (gC/m2/yr)
!     mimicsfluxAn%ChrespAn(npt)                - heterotrphic respiration flux (gC/m2/yr)
!     mimicsfluxAn%CSOMpInputAn(npt)            - inputs to SOMp (gC/m2/yr)
!     mimicspoolAn%ClitterAn(npt,METBC)         - metabolic litter pool (gC/m2)
!     mimicspoolAn%ClitterAn(npt,STRUC)         - structural litter pool (gC/m2)
!     mimicspoolAn%CmicrobeAn(npt,RSEL)         - r-selected microbe pool (gC/m2)
!     mimicspoolAn%CmicrobeAn(npt,KSEL)         - K-selected microbe pool (gC/m2)
!     mimicspoolAn%CsoilAn(npt,AVAL)            - available SOM pool (gC/m2)
!     mimicspoolAn%CsoilAn(npt,CHEM)            - chemically protected SOM pool (gC/m2)
!     mimicspoolAn%CsoilAn(npt,PHYS)            - physically protected SOM pool (gC/m2)

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

   status = nf_def_var(ncid, 'cSOMpIn', NF_REAL, 3, dims, varid_cSOMpIn)
   if (status /= nf_noerr) call handle_err(status, "cSOMpIn")

   status = nf_def_var(ncid, 'cLitInput_metb', NF_REAL, 3, dims, varid_cLitIn_m)
   if (status /= nf_noerr) call handle_err(status, "cLitInput_metb")

   status = nf_def_var(ncid, 'cLitInput_struc', NF_REAL, 3, dims, varid_cLitIn_s)
   if (status /= nf_noerr) call handle_err(status, "cLitInput_struc")


   ! Added Annual N variables. -mdh 11/30/2019

   status = nf_def_var(ncid, 'nLITm', NF_REAL, 3, dims, varid_nLITm)
   if (status /= nf_noerr) call handle_err(status, "nLITm")

   status = nf_def_var(ncid, 'nLITs', NF_REAL, 3, dims, varid_nLITs)
   if (status /= nf_noerr) call handle_err(status, "nLITs")

   status = nf_def_var(ncid, 'nMICr', NF_REAL, 3, dims, varid_nMICr)
   if (status /= nf_noerr) call handle_err(status, "nMICr")

   status = nf_def_var(ncid, 'nMICk', NF_REAL, 3, dims, varid_nMICk)
   if (status /= nf_noerr) call handle_err(status, "nMICk")

   status = nf_def_var(ncid, 'nSOMa', NF_REAL, 3, dims, varid_nSOMa)
   if (status /= nf_noerr) call handle_err(status, "nSOMa")

   status = nf_def_var(ncid, 'nSOMc', NF_REAL, 3, dims, varid_nSOMc)
   if (status /= nf_noerr) call handle_err(status, "nSOMc")

   status = nf_def_var(ncid, 'nSOMp', NF_REAL, 3, dims, varid_nSOMp)
   if (status /= nf_noerr) call handle_err(status, "nSOMp")

   status = nf_def_var(ncid, 'DIN', NF_REAL, 3, dims, varid_DIN)
   if (status /= nf_noerr) call handle_err(status, "DIN")

   status = nf_def_var(ncid, 'nLitInput_metb', NF_REAL, 3, dims, varid_nLitIn_m)
   if (status /= nf_noerr) call handle_err(status, "nLitInput_metb")

   status = nf_def_var(ncid, 'nLitInput_struc', NF_REAL, 3, dims, varid_nLitIn_s)
   if (status /= nf_noerr) call handle_err(status, "nLitInput_struc")


   status = nf_def_var(ncid, 'thetaLiq', NF_REAL, 3, dims, varid_thetaLiq)
   if (status /= nf_noerr) call handle_err(status, "thetaLiq")

   status = nf_def_var(ncid, 'thetaFrzn', NF_REAL, 3, dims, varid_thetaFrzn)
   if (status /= nf_noerr) call handle_err(status, "thetaFrzn")

!  status = nf_def_var(ncid, 'fT', NF_REAL, 3, dims, varid_fT)
!  if (status /= nf_noerr) call handle_err(status, "f(T)")

   status = nf_def_var(ncid, 'fW', NF_REAL, 3, dims, varid_fW)
   if (status /= nf_noerr) call handle_err(status, "f(W)")


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

   ! Attributes of cSOMpIn variable
   attr_name = 'inputs to physically protected soil pool'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cSOMpIn, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLitInput_metb variable
   attr_name = 'metabolic litter C inputs'
   attr_units = 'gC m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_cLitIn_m, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of cLitInput_struc variable
   attr_name = 'structural litter C inputs'
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
   attr_name = 'K-selected microbial soil organic matter carbon'
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


   ! Added N variables. -mdh 11/30/2019

   ! Attributes of nLITm variable
   attr_name = 'metabolic litter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nLITm, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nLITs variable
   attr_name = 'structural litter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nLITs, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMICr variable
   attr_name = 'r-selected microbial soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nMICr, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nMICk variable
   attr_name = 'K-selected microbial soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nMICk, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nSOMa variable
   attr_name = 'active soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nSOMa, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nSOMc variable
   attr_name = 'chemically protected soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nSOMc, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nSOMp variable
   attr_name = 'physically protected soil organic matter nitrogen'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_nSOMp, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of DIN variable
   attr_name = 'dissolved inorganic N available to microbes'
   attr_units = 'gN m-2'
   call PutVariableAttributeReal(ncid, varid_DIN, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nLitInput_metb variable
   attr_name = 'metabolic litter N inputs'
   attr_units = 'gN m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_nLitIn_m, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of nLitInput_struc variable
   attr_name = 'structural litter N inputs'
   attr_units = 'gN m-2 year-1'
   call PutVariableAttributeReal(ncid, varid_nLitIn_s, attr_name, attr_units, MISSING_VALUE)


   ! Attributes of thetaLiq variable
   attr_name = 'fraction of liquid soil water saturation (0.0-1.0)'
   attr_units = 'fraction'
   call PutVariableAttributeReal(ncid, varid_thetaLiq, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of thetaFrzn variable
   attr_name = 'fraction of frozen soil water saturation (0.0-1.0)'
   attr_units = 'fraction'
   call PutVariableAttributeReal(ncid, varid_thetaFrzn, attr_name, attr_units, MISSING_VALUE)

!  ! Attributes of fT variable
!  attr_name = 'soil temperature multiplier on decomposition (0.0-1.0)'
!  attr_units = 'multiplier'
!  call PutVariableAttributeReal(ncid, varid_fT, attr_name, attr_units, MISSING_VALUE)

   ! Attributes of fW variable
   attr_name = 'soil moisture multiplier on decomposition (0.0-1.0)'
   attr_units = 'multiplier'
   call PutVariableAttributeReal(ncid, varid_fW, attr_name, attr_units, MISSING_VALUE)



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
   allocate(var11(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var12(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var13(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var14(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))
   allocate(var15(1:clmgrid%nlon,1:clmgrid%nlat,1:ntimes))

!  Average annual values for fluxes and pools
!  mimicsfluxAn%CLitInputAn(npt,METBC)  - metabolic litter inputs (gC/m2/yr)
!  mimicsfluxAn%CLitInputAn(npt,STRUC)  - structural litter inputs (gC/m2/yr)
!  mimicsfluxAn%ChrespAn(npt)           - heterotrophic respiration flux (gC/m2/yr)
!  mimicsfluxAn%CSOMpInputAn(npt)       - inputs to SOMp (gC/m2/yr)
!  mimicspoolAn%ClitterAn(npt,METBC)    - metabolic litter pool (gC/m2)
!  mimicspoolAn%ClitterAn(npt,STRUC)    - structural litter pool (gC/m2)
!  mimicspoolAn%CmicrobeAn(npt,RSEL)    - r-selected microbe pool (gC/m2)
!  mimicspoolAn%CmicrobeAn(npt,KSEL)    - K-selected microbe pool (gC/m2)
!  mimicspoolAn%CsoilAn(npt,AVAL)       - available SOM pool (gC/m2)
!  mimicspoolAn%CsoilAn(npt,CHEM)       - chemically protected SOM pool (gC/m2)
!  mimicspoolAn%CsoilAn(npt,PHYS)       - physically protected SOM pool (gC/m2)
!  mimicspoolAn%thetaLiqAn(npt)         - MIMICS mean annual fraction of liquid water saturation
!  mimicspoolAn%thetaFrznAn(npt)        - MIMICS mean annual fraction of frozen water saturation
!  mimicspoolAn%fTAn(npt)               - MIMICS mean annual soil temperature function
!  mimicspoolAn%fWAn(npt)               - MIMICS mean annual soil moisture function

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
   var11(:,:,:) = MISSING_VALUE
   var12(:,:,:) = MISSING_VALUE
   var13(:,:,:) = MISSING_VALUE
   var14(:,:,:) = MISSING_VALUE
   var15(:,:,:) = MISSING_VALUE

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
      var11(ilon,ilat,itime) = mimicspoolAn%thetaLiqAn(npt)   
      var12(ilon,ilat,itime) = mimicspoolAn%thetaFrznAn(npt) 
!     var13(ilon,ilat,itime) = mimicspoolAn%fTAn(npt)   
      var14(ilon,ilat,itime) = mimicspoolAn%fWAn(npt) 
      var15(ilon,ilat,itime) = mimicsfluxAn%CSOMpInputAn(npt) 
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

   status =  nf_put_var(ncid, varid_thetaLiq, var11, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(thetaLiq)")

   status =  nf_put_var(ncid, varid_thetaFrzn, var12, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(thetaFrzn)")

!  status =  nf_put_var(ncid, varid_fT, var13, start3, count3)
!  if (status /= nf_noerr) call handle_err(status, "put_var(fT)")

   status =  nf_put_var(ncid, varid_fW, var14, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(fW)")

   status =  nf_put_var(ncid, varid_cSOMpIn, var15, start3, count3)
   if (status /= nf_noerr) call handle_err(status, "put_var(cSOMpIn)")


   ! Added annual N variables. -mdh 12/1/2019

!  Average annual values for fluxes and pools
!  mimicsfluxAn%NLitInputAn(npt,METBC)  - metabolic litter inputs (gN/m2/yr)
!  mimicsfluxAn%NLitInputAn(npt,STRUC)  - structural litter inputs (gN/m2/yr)
!  mimicspoolAn%NlitterAn(npt,METBC)    - metabolic litter pool (gN/m2)
!  mimicspoolAn%NlitterAn(npt,STRUC)    - structural litter pool (gN/m2)
!  mimicspoolAn%NmicrobeAn(npt,RSEL)    - r-selected microbe pool (gN/m2)
!  mimicspoolAn%NmicrobeAn(npt,KSEL)    - K-selected microbe pool (gN/m2)
!  mimicspoolAn%NsoilAn(npt,AVAL)       - available SOM pool (gN/m2)
!  mimicspoolAn%NsoilAn(npt,CHEM)       - chemically protected SOM pool (gN/m2)
!  mimicspoolAn%NsoilAn(npt,PHYS)       - physically protected SOM pool (gN/m2)
!  mimicspoolAn%DINAn(npt)              - dissolved inorganic N available to microbes (gN/m2)

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
      endif

      var1(ilon,ilat,itime)  = mimicsfluxAn%NLitInputAn(npt,METBC)
      var2(ilon,ilat,itime)  = mimicsfluxAn%NLitInputAn(npt,STRUC)
      var3(ilon,ilat,itime)  = mimicspoolAn%NlitterAn(npt,METBC)
      var4(ilon,ilat,itime)  = mimicspoolAn%NlitterAn(npt,STRUC)
      var5(ilon,ilat,itime)  = mimicspoolAn%NmicrobeAn(npt,RSEL)
      var6(ilon,ilat,itime)  = mimicspoolAn%NmicrobeAn(npt,KSEL)
      var7(ilon,ilat,itime)  = mimicspoolAn%NsoilAn(npt,AVAL)
      var8(ilon,ilat,itime)  = mimicspoolAn%NsoilAn(npt,CHEM)
      var9(ilon,ilat,itime)  = mimicspoolAn%NsoilAn(npt,PHYS)
      var10(ilon,ilat,itime) = mimicspoolAn%DINAn(npt)

   end do

   status =  nf_put_vara_real(ncid, varid_nLitIn_m, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_vara_real(nLitIn_m)")

   !!status =  nf_put_var(ncid, varid_nLitIn_s, var2, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nLitIn_s, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLitIn_s)")

   !!status =  nf_put_var(ncid, varid_nLITm, var3, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nLITm, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLITm)")

   !!status =  nf_put_var(ncid, varid_nLITs, var4, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nLITs, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLITs)")

   !!status =  nf_put_var(ncid, varid_nMICr, var5, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nMICr, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMICr)")

   !!status =  nf_put_var(ncid, varid_nMICk, var6, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nMICk, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMICk)")

   !!status =  nf_put_var(ncid, varid_nSOMa, var7, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nSOMa, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSOMa)")

   !!status =  nf_put_var(ncid, varid_nSOMc, var8, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nSOMc, start3, count3, var8)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSOMc)")

   !!status =  nf_put_var(ncid, varid_nSOMp, var9, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nSOMp, start3, count3, var9)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSOMp)")

   !!status =  nf_put_var(ncid, varid_DIN, var10, start3, count3)
   status =  nf_put_vara_real(ncid, varid_DIN, start3, count3, var10)
   if (status /= nf_noerr) call handle_err(status, "put_var(DIN)")


   deallocate(var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15,IGBP_PFT,landarea)

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
      character(len=*), intent(in) :: filename_ncOut    ! NetCDF output file name 
      integer, intent(in) :: mp                         ! number of points with valid data
      integer, intent(in) :: year                       ! output year
      integer, intent(in) :: iday                       ! output day (1..365)
!   LOCAL VARIABLES
      integer :: i
      integer :: ncid                           ! netcdf file ID
      integer :: status                         ! function return status
!     integer :: dimid_mp                       ! netcdf dimension id
      integer :: dimid_lat                      ! netcdf dimension id
      integer :: dimid_lon                      ! netcdf dimension id
      integer :: dimid_time                     ! netcdf dimension id
      integer :: nlon, nlat, ntimes             ! Dimension sizes for NetCDf file
      integer :: nwrtimes                       ! Number of times that will be written when this subroutine is called
      integer :: dims(3)                        ! Array of NetCDF dimension IDs for defining variables
      integer :: start1(1), count1(1)           ! start and count arrays for writing 1-D data from netcdf files
      integer :: start2(2), count2(2)           ! start and count arrays for writing 2-D data from netcdf files
      integer :: start3(3), count3(3)           ! start and count arrays for writing 3-D data from netcdf files
      integer :: varid_lon, varid_lat           ! NetCDF variable ID for latitude and longitude
      integer :: varid_time                     ! NetCDF variable ID for time
      integer :: varid_day                      ! NetCDF variable ID for day
      integer :: varid_mask                     ! NetCDF variable ID for cellMissing(nlon,nlat)
      integer :: varid_cellid                   ! NetCDF variable ID for cellid(nlon,nlat)
      integer :: varid_igbp                     ! NetCDF variable ID for IGBP_PFT(nlon,nlat)
      integer :: varid_landarea                 ! NetCDF variable ID for land area
      integer :: varid_cLITm                    ! NetCDF variable ID for metablic litter C
      integer :: varid_cLITs                    ! NetCDF variable ID for structural litter C
      integer :: varid_cMICr                    ! NetCDF variable ID for microbe r-selected pool C
      integer :: varid_cMICk                    ! NetCDF variable ID for microbe K-selected pool C
      integer :: varid_cSOMa                    ! NetCDF variable ID for available soil C
      integer :: varid_cSOMc                    ! NetCDF variable ID for chemically protected soil C
      integer :: varid_cSOMp                    ! NetCDF variable ID for physically protected soil C
      integer :: varid_cHresp                   ! NetCDF variable ID for soil heterotrophic respiration C
      integer :: varid_cSOMpIn                  ! NetCDF variable ID for physically protected soil C inputs
      integer :: varid_cLitIn_m                 ! NetCDF variable ID for metabolic litter Input C
      integer :: varid_cLitIn_s                 ! NetCDF variable ID for structural litter Input C

      integer :: varid_nLITm                    ! NetCDF variable ID for metablic litter N
      integer :: varid_nLITs                    ! NetCDF variable ID for structural litter N
      integer :: varid_nMICr                    ! NetCDF variable ID for microbe r-selected pool N
      integer :: varid_nMICk                    ! NetCDF variable ID for microbe K-selected pool N
      integer :: varid_nSOMa                    ! NetCDF variable ID for available soil N
      integer :: varid_nSOMc                    ! NetCDF variable ID for chemically protected soil N
      integer :: varid_nSOMp                    ! NetCDF variable ID for physically protected soil N
      integer :: varid_DIN                      ! NetCDF variable ID for dissolved inorganic N
      integer :: varid_nLitIn_m                 ! NetCDF variable ID for metabolic litter Input N
      integer :: varid_nLitIn_s                 ! NetCDF variable ID for structural litter Input N

      integer :: varid_thetaLiq                 ! NetCDF variable ID for fraction of liquid soil water saturation
      integer :: varid_thetaFrzn                ! NetCDF variable ID for fraction of frozen soil water saturation
!     integer :: varid_fT                       ! NetCDF variable ID for soil temperature multiplier on decompostion
      integer :: varid_fW                       ! NetCDF variable ID for soil moisture multiplier on decompostion
      character*100 :: attr_name                ! String for assigning global and variable attributes
      character*100 :: attr_units               ! String for assigning global and variable attributes
      character*10 :: date_string               ! String for assigning date to global attributes
      character*8 :: time_string                ! String for assigning time to global attributes
      integer :: verbose=0
      integer :: npt, ilon, ilat, itime
      integer, allocatable :: IGBP_PFT(:,:)     ! IGBP_PFT(nlon,nlat) IGBP PFT classification (1-18)
      integer, allocatable :: days(:)           ! day array (1..ntimes)
      real(4), allocatable :: time(:)           ! time array (time in years)
      real(4), allocatable :: landarea(:,:)     ! landarea(nlon,nlat) km^2
      real(4), allocatable :: var1(:,:,:)       ! gridded output variable
      real(4), allocatable :: var2(:,:,:)       ! gridded output variable
      real(4), allocatable :: var3(:,:,:)       ! gridded output variable
      real(4), allocatable :: var4(:,:,:)       ! gridded output variable
      real(4), allocatable :: var5(:,:,:)       ! gridded output variable
      real(4), allocatable :: var6(:,:,:)       ! gridded output variable
      real(4), allocatable :: var7(:,:,:)       ! gridded output variable
      real(4), allocatable :: var8(:,:,:)       ! gridded output variable
      real(4), allocatable :: var9(:,:,:)       ! gridded output variable
      real(4), allocatable :: var10(:,:,:)      ! gridded output variable
      real(4), allocatable :: var11(:,:,:)      ! gridded output variable
      real(4), allocatable :: var12(:,:,:)      ! gridded output variable
      real(4), allocatable :: var13(:,:,:)      ! gridded output variable
      real(4), allocatable :: var14(:,:,:)      ! gridded output variable
      real(4), allocatable :: var15(:,:,:)      ! gridded output variable
      real(r_2) :: unitConv                     ! mgC/cm3 * depth(cm)* (1g/10^3mg)*(10^4cm2)/m2 = gC/m2

      if (verbose .gt. 0) print *, iday, "Writing output to file ", trim(filename_ncOut), "..."

! Output these variables in NetCDF file:
!     veg%iveg(npt)                           - IGBP PFTs
!     casamet%lat(npt)                        - latitudes
!     casamet%lon(npt)                        - longitudes
!     casamet%areacell(npt)*(1.0e-6)          - landarea (km^2)
!   Daily values for fluxes and pools     
!     mimicsflux%ClitInput(npt,metbc)         - metabolic litter inputs (mgC/cm3/day)
!     mimicsflux%ClitInput(npt,struc)         - structural litter inputs (mgC/cm3/day)
!     mimicsflux%Chresp(npt)                  - heterotrphic respiration flux (mgC/cm3/day)
!     mimicsflux%CSOMpInput(npt)              - inputs to SOMp (mgC/cm3/day)
!     mimicspool%LITm(npt)                    - metabolic litter pool (mgC/cm3)
!     mimicspool%LITs(npt)                    - structural litter pool (mgC/cm3)
!     mimicspool%MICr(npt)                    - r-selected microbe pool (mgC/cm3)
!     mimicspool%MICk(npt)                    - K-selected microbe pool (mgC/cm3)
!     mimicspool%SOMa(npt)                    - available SOM pool (mgC/cm3)
!     mimicspool%SOMc(npt)                    - chemically protected SOM pool (mgC/cm3)
!     mimicspool%SOMp(npt)                    - physically protected SOM pool (mgC/cm3)

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
   
      status = nf_def_var(ncid, 'cSOMpIn', NF_REAL, 3, dims, varid_cSOMpIn)
      if (status /= nf_noerr) call handle_err(status, "cSOMpIn")

      status = nf_def_var(ncid, 'cLitInput_metb', NF_REAL, 3, dims, varid_cLitIn_m)
      if (status /= nf_noerr) call handle_err(status, "cLitInput_metb")
   
      status = nf_def_var(ncid, 'cLitInput_struc', NF_REAL, 3, dims, varid_cLitIn_s)
      if (status /= nf_noerr) call handle_err(status, "cLitInput_struc")


      ! Added daily N variables. -mdh 11/30/2019

      status = nf_def_var(ncid, 'nLITm', NF_REAL, 3, dims, varid_nLITm)
      if (status /= nf_noerr) call handle_err(status, "nLITm")
   
      status = nf_def_var(ncid, 'nLITs', NF_REAL, 3, dims, varid_nLITs)
      if (status /= nf_noerr) call handle_err(status, "nLITs")
   
      status = nf_def_var(ncid, 'nMICr', NF_REAL, 3, dims, varid_nMICr)
      if (status /= nf_noerr) call handle_err(status, "nMICr")
   
      status = nf_def_var(ncid, 'nMICk', NF_REAL, 3, dims, varid_nMICk)
      if (status /= nf_noerr) call handle_err(status, "nMICk")
   
      status = nf_def_var(ncid, 'nSOMa', NF_REAL, 3, dims, varid_nSOMa)
      if (status /= nf_noerr) call handle_err(status, "nSOMa")
   
      status = nf_def_var(ncid, 'nSOMc', NF_REAL, 3, dims, varid_nSOMc)
      if (status /= nf_noerr) call handle_err(status, "nSOMc")
   
      status = nf_def_var(ncid, 'nSOMp', NF_REAL, 3, dims, varid_nSOMp)
      if (status /= nf_noerr) call handle_err(status, "nSOMp")

      status = nf_def_var(ncid, 'DIN', NF_REAL, 3, dims, varid_DIN)
      if (status /= nf_noerr) call handle_err(status, "DIN")
   
      status = nf_def_var(ncid, 'nLitInput_metb', NF_REAL, 3, dims, varid_nLitIn_m)
      if (status /= nf_noerr) call handle_err(status, "nLitInput_metb")
   
      status = nf_def_var(ncid, 'nLitInput_struc', NF_REAL, 3, dims, varid_nLitIn_s)
      if (status /= nf_noerr) call handle_err(status, "nLitInput_struc")


      status = nf_def_var(ncid, 'thetaLiq', NF_REAL, 3, dims, varid_thetaLiq)
      if (status /= nf_noerr) call handle_err(status, "thetaLiq")

      status = nf_def_var(ncid, 'thetaFrzn', NF_REAL, 3, dims, varid_thetaFrzn)
      if (status /= nf_noerr) call handle_err(status, "thetaFrzn")

!     status = nf_def_var(ncid, 'fT', NF_REAL, 3, dims, varid_fT)
!     if (status /= nf_noerr) call handle_err(status, "f(T)")

      status = nf_def_var(ncid, 'fW', NF_REAL, 3, dims, varid_fW)
      if (status /= nf_noerr) call handle_err(status, "f(W)")
   
   
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

      ! Attributes of cSOMpIn variable
      attr_name = 'inputs to protected soil C pool'
      attr_units = 'gC m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_cSOMpIn, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of cLitInput_metb variable
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
      attr_name = 'K-selected microbial soil organic matter carbon'
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


      ! Added N variables. -mdh 11/30/2019
   
      ! Attributes of nLITm variable
      attr_name = 'metabolic litter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nLITm, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nLITs variable
      attr_name = 'structural litter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nLITs, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nMICr variable
      attr_name = 'r-selected microbial soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nMICr, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nMICk variable
      attr_name = 'K-selected microbial soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nMICk, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nSOMa variable
      attr_name = 'active soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nSOMa, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nSOMc variable
      attr_name = 'chemically protected soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nSOMc, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nSOMp variable
      attr_name = 'physically protected soil organic matter nitrogen'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_nSOMp, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of DIN variable
      attr_name = 'dissolved inorganic N available to microbes'
      attr_units = 'gN m-2'
      call PutVariableAttributeReal(ncid, varid_DIN, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nLitInput_metb variable
      attr_name = 'metabolic litter N inputs'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nLitIn_m, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of nLitInput_struc variable
      attr_name = 'structural litter N inputs'
      attr_units = 'gN m-2 day-1'
      call PutVariableAttributeReal(ncid, varid_nLitIn_s, attr_name, attr_units, MISSING_VALUE)

   
      ! Attributes of thetaLiq variable
      attr_name = 'fraction of liquid soil water saturation (0.0-1.0)'
      attr_units = 'fraction'
      call PutVariableAttributeReal(ncid, varid_thetaLiq, attr_name, attr_units, MISSING_VALUE)

      ! Attributes of thetaFrzn variable
      attr_name = 'fraction of frozen soil water saturation (0.0-1.0)'
      attr_units = 'fraction'
      call PutVariableAttributeReal(ncid, varid_thetaFrzn, attr_name, attr_units, MISSING_VALUE)

!     ! Attributes of fT variable
!     attr_name = 'soil temperature multiplier on decomposition (0.0-1.0)'
!     attr_units = 'multiplier'
!     call PutVariableAttributeReal(ncid, varid_fT, attr_name, attr_units, MISSING_VALUE)
   
      ! Attributes of fW variable
      attr_name = 'soil moisture multiplier on decomposition (0.0-1.0)'
      attr_units = 'multiplier'
      call PutVariableAttributeReal(ncid, varid_fW, attr_name, attr_units, MISSING_VALUE)

   
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
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cLITm)")
   
      status = nf_inq_varid(ncid, 'cLITs', varid_cLITs)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cLITs)")
   
      status = nf_inq_varid(ncid, 'cMICr', varid_cMICr)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cMICr)")
   
      status = nf_inq_varid(ncid, 'cMICk', varid_cMICk)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cMICk)")
   
      status = nf_inq_varid(ncid, 'cSOMa', varid_cSOMa)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cSOMa)")
   
      status = nf_inq_varid(ncid, 'cSOMc', varid_cSOMc)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cSOMc)")
   
      status = nf_inq_varid(ncid, 'cSOMp', varid_cSOMp)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cSOMp)")
   
      status = nf_inq_varid(ncid, 'cHresp', varid_cHresp)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cHresp)")
   
      status = nf_inq_varid(ncid, 'cSOMpIn', varid_cSOMpIn)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cSOMpIn)")
   
      status = nf_inq_varid(ncid, 'cLitInput_metb', varid_cLitIn_m)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cLitInput_metb)")
   
      status = nf_inq_varid(ncid, 'cLitInput_struc', varid_cLitIn_s)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(cLitInput_struc)")


      ! Added daily N variables. -mdh 11/30/2019

      status = nf_inq_varid(ncid, 'nLITm',varid_nLITm)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nLITm)")
   
      status = nf_inq_varid(ncid, 'nLITs', varid_nLITs)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nLITs)")
   
      status = nf_inq_varid(ncid, 'nMICr', varid_nMICr)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nMICr)")
   
      status = nf_inq_varid(ncid, 'nMICk', varid_nMICk)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nMICk)")
   
      status = nf_inq_varid(ncid, 'nSOMa', varid_nSOMa)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nSOMa)")
   
      status = nf_inq_varid(ncid, 'nSOMc', varid_nSOMc)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nSOMc)")
   
      status = nf_inq_varid(ncid, 'nSOMp', varid_nSOMp)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nSOMp)")
   
      status = nf_inq_varid(ncid, 'DIN', varid_DIN)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(DIN)")
   
      status = nf_inq_varid(ncid, 'nLitInput_metb', varid_nLitIn_m)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nLitInput_metb)")
   
      status = nf_inq_varid(ncid, 'nLitInput_struc', varid_nLitIn_s)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(nLitInput_struc)")


      status = nf_inq_varid(ncid, 'thetaLiq', varid_thetaLiq)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(thetaLiq)")

      status = nf_inq_varid(ncid, 'thetaFrzn', varid_thetaFrzn)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(thetaFrzn)")

!     status = nf_inq_varid(ncid, 'fT', varid_fT)
!     if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(fT)")
   
      status = nf_inq_varid(ncid, 'fW', varid_fW)
      if (status /= nf_noerr) call handle_err(status, "nf_inq_varid(fW)")

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
   allocate(var11(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var12(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var13(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var14(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))
   allocate(var15(1:clmgrid%nlon,1:clmgrid%nlat,1:nwrtimes))

!  Daily values for fluxes and pools
!  mimicsflux%ClitInput(npt,metbc)   - metabolic litter inputs (mgC/cm3/day)
!  mimicsflux%ClitInput(npt,struc)   - structural litter inputs (mgC/cm3/day)
!  mimicsflux%Chresp(npt)            - heterotrphic respiration flux (mgC/cm3/day)
!  mimicsflux%CSOMpInput(npt)        - inputs to SOMp (mgC/cm3/day)
!  mimicspool%LITm(npt)              - metabolic litter pool (mgC/cm3)
!  mimicspool%LITs(npt)              - structural litter pool (mgC/cm3)
!  mimicspool%MICr(npt)              - r-selected microbe pool (mgC/cm3)
!  mimicspool%MICk(npt)              - K-selected microbe pool (mgC/cm3)
!  mimicspool%SOMa(npt)              - available SOM pool (mgC/cm3)
!  mimicspool%SOMc(npt)              - chemically protected SOM pool (mgC/cm3)
!  mimicspool%SOMp(npt)              - physically protected SOM pool (mgC/cm3)
!  mimicspool%thetaLiq(npt)          - fraction of liquid soil water saturation (0.0-1.0)
!  mimicspool%thetaFrzn(npt)         - fraction of frozen soil water saturation (0.0-1.0)
!  mimicspool%fT(npt)                - MIMICS soil temperature function
!  mimicspool%fW(npt)                - MIMICS soil moisture function

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
   var11(:,:,:) = MISSING_VALUE
   var12(:,:,:) = MISSING_VALUE
   var13(:,:,:) = MISSING_VALUE
   var14(:,:,:) = MISSING_VALUE
   var15(:,:,:) = MISSING_VALUE



   itime = 1
   do npt = 1, mp
      ilon = casamet%ilon(npt)
      ilat = casamet%ilat(npt)
      unitConv = 10.0 * mimicsbiome%depth(veg%iveg(npt))     ! mgC/cm3 * depth(cm)* (1g/10^3mg)*(10^4cm2)/m2 = gC/m2
      !!unitConv = 1.0 !! FOR TESTING - comment out after mimics-cn variables have been checked. -mdh 12/2/2019
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
      var11(ilon,ilat,itime) = mimicspool%thetaLiq(npt)   
      var12(ilon,ilat,itime) = mimicspool%thetaFrzn(npt) 
!     var13(ilon,ilat,itime) = mimicspool%fT(npt)   
      var14(ilon,ilat,itime) = mimicspool%fW(npt) 
      var15(ilon,ilat,itime) = mimicsflux%CSOMpInput(npt) * unitConv
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

!! status =  nf_put_var(ncid, varid_thetaLiq, var11, start3, count3)
   status =  nf_put_vara_real(ncid, varid_thetaLiq, start3, count3, var11)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(thetaLiq)")

!! status =  nf_put_var(ncid, varid_thetaFrzn, var12, start3, count3)
   status =  nf_put_vara_real(ncid, varid_thetaFrzn, start3, count3, var12)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(thetaFrzn)")

!! status =  nf_put_var(ncid, varid_fT, var13, start3, count3)
!  status =  nf_put_vara_real(ncid, varid_fT, start3, count3, var13)
!  if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(fT)")

!! status =  nf_put_var(ncid, varid_fW, var14, start3, count3)
   status =  nf_put_vara_real(ncid, varid_fW, start3, count3, var14)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(fW)")

!! status =  nf_put_var(ncid, varid_cSOMpIn, var15, start3, count3)
   status =  nf_put_vara_real(ncid, varid_cSOMpIn, start3, count3, var15)
   if (status /= nf_noerr) call handle_err(status, "nf_put_vara_real(cSOMpIn)")


   ! Added daily N variables. -mdh 12/1/2019

!  Daily values for fluxes and pools
!  mimicsflux%NLitInput(npt,METBC)  - metabolic litter inputs (gN/m2/day)
!  mimicsflux%NLitInput(npt,STRUC)  - structural litter inputs (gN/m2/day)
!  mimicspool%Nlitter(npt,METBC)    - metabolic litter pool (gN/m2)
!  mimicspool%Nlitter(npt,STRUC)    - structural litter pool (gN/m2)
!  mimicspool%Nmicrobe(npt,RSEL)    - r-selected microbe pool (gN/m2)
!  mimicspool%Nmicrobe(npt,KSEL)    - K-selected microbe pool (gN/m2)
!  mimicspool%Nsoil(npt,AVAL)       - available SOM pool (gN/m2)
!  mimicspool%Nsoil(npt,CHEM)       - chemically protected SOM pool (gN/m2)
!  mimicspool%Nsoil(npt,PHYS)       - physically protected SOM pool (gN/m2)
!  mimicspool%DIN(npt)              - dissolved inorganic N available to microbes (gN/m2)

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
      endif

      unitConv = 10.0 * mimicsbiome%depth(veg%iveg(npt))     ! mgN/cm3 * depth(cm)* (1g/10^3mg)*(10^4cm2)/m2 = gN/m2
      !!unitConv = 1.0 !! FOR TESTING - comment out after mimics-cn variables have been checked. -mdh 12/2/2019
      var1(ilon,ilat,itime)  = mimicsflux%NLitInput(npt,METBC) * unitConv
      var2(ilon,ilat,itime)  = mimicsflux%NLitInput(npt,STRUC) * unitConv
      var3(ilon,ilat,itime)  = mimicspool%LITmN(npt) * unitConv
      var4(ilon,ilat,itime)  = mimicspool%LITsN(npt) * unitConv
      var5(ilon,ilat,itime)  = mimicspool%MICrN(npt) * unitConv
      var6(ilon,ilat,itime)  = mimicspool%MICkN(npt) * unitConv
      var7(ilon,ilat,itime)  = mimicspool%SOMaN(npt) * unitConv
      var8(ilon,ilat,itime)  = mimicspool%SOMcN(npt) * unitConv
      var9(ilon,ilat,itime)  = mimicspool%SOMpN(npt) * unitConv
      var10(ilon,ilat,itime) = mimicspool%DIN(npt)   * unitConv

   end do

   status =  nf_put_vara_real(ncid, varid_nLitIn_m, start3, count3, var1)
   if (status /= nf_noerr) call handle_err(status, "put_vara_real(nLitIn_m)")

   !!status =  nf_put_var(ncid, varid_nLitIn_s, var2, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nLitIn_s, start3, count3, var2)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLitIn_s)")

   !!status =  nf_put_var(ncid, varid_nLITm, var3, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nLITm, start3, count3, var3)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLITm)")

   !!status =  nf_put_var(ncid, varid_nLITs, var4, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nLITs, start3, count3, var4)
   if (status /= nf_noerr) call handle_err(status, "put_var(nLITs)")

   !!status =  nf_put_var(ncid, varid_nMICr, var5, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nMICr, start3, count3, var5)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMICr)")

   !!status =  nf_put_var(ncid, varid_nMICk, var6, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nMICk, start3, count3, var6)
   if (status /= nf_noerr) call handle_err(status, "put_var(nMICk)")

   !!status =  nf_put_var(ncid, varid_nSOMa, var7, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nSOMa, start3, count3, var7)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSOMa)")

   !!status =  nf_put_var(ncid, varid_nSOMc, var8, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nSOMc, start3, count3, var8)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSOMc)")

   !!status =  nf_put_var(ncid, varid_nSOMp, var9, start3, count3)
   status =  nf_put_vara_real(ncid, varid_nSOMp, start3, count3, var9)
   if (status /= nf_noerr) call handle_err(status, "put_var(nSOMp)")

   !!status =  nf_put_var(ncid, varid_DIN, var10, start3, count3)
   status =  nf_put_vara_real(ncid, varid_DIN, start3, count3, var10)
   if (status /= nf_noerr) call handle_err(status, "put_var(DIN)")




   deallocate(var1,var2,var3,var4,var5,var6,var7,var8,var9,var10,var11,var12,var13,var14,var15)

   status = nf_close(ncid)

   if (verbose .gt. 0) print *, iday, "Done writing output to ", trim(filename_ncOut), "..."

END SUBROUTINE WritePoolFluxNcFile_mimics_daily

!-------------------------------------------------------------------------------------
! Write MIMICS C only pools and other quantities to .csv output daily

SUBROUTINE WritePointMIMICS(unit1, sPtFileName, npt, mp, iYrCnt, idoy, &
    cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
    LITmin, MICtrn, SOMmin, DEsorb, OXIDAT, &
    dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil, Cbalance) 

    USE define_types
    USE casadimension
    USE casaparm
    USE casavariable
    USE clmgridvariable
    USE mimicsdimension
    USE mimicsparam
    USE mimicsvariable
    implicit none

!   ARGUMENTS
    integer, intent(in)            :: unit1         ! FORTRAN file unit
    character(len=*), intent(in)   :: sPtFileName   ! .csv output file name 
    integer, intent(IN)            :: mp, npt       ! # points, point index
    integer, intent(IN)            :: iYrCnt, idoy  ! simulation year count, day of year
    real(r_2), dimension(mp),intent(IN) :: cleaf2met,cleaf2str,croot2met,croot2str
    real(r_2), dimension(mp),intent(IN) :: cwd2str,cwd2co2,cwood2cwd
    real(r_2), intent(IN)          :: LITmin(4), MICtrn(6), SOMmin(2), DEsorb, OXIDAT
    real(r_2), intent(IN)          :: dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk
    real(r_2), intent(IN)          :: Tsoil, Cbalance


    open(unit1,file=sPtFileNameMIMICS, access='APPEND')

    write(unit1,102) npt,casamet%ijgcm(npt),iYrCnt,idoy,casamet%tsoilavg(npt),casamet%moistavg(npt),&
                   mimicsflux%Chresp(npt),mimicsflux%ClitInput(npt,metbc),mimicsflux%ClitInput(npt,struc), &
                   mimicspool%LITm(npt),mimicspool%LITs(npt),mimicspool%MICr(npt), &
                   mimicspool%MICk(npt),mimicspool%SOMa(npt),mimicspool%SOMc(npt),mimicspool%SOMp(npt), &
                   dLITm, dLITs, dMICr, dMICk, dSOMa, dSOMc, dSOMp,&

                   LITmin(1), LITmin(2), LITmin(3), LITmin(4), &
                   MICtrn(1), MICtrn(2), MICtrn(3), MICtrn(4), MICtrn(5), MICtrn(6), &
                   SOMmin(1), SOMmin(2), DEsorb, OXIDAT, mimicsbiome%fmet(npt), &

                   mimicsbiome%tauMod(npt), mimicsbiome%tauR(npt), mimicsbiome%tauK(npt), &
                   mimicsbiome%Vmax(npt,R1),mimicsbiome%Vmax(npt,R2),mimicsbiome%Vmax(npt,R3), &
                   mimicsbiome%Vmax(npt,K1),mimicsbiome%Vmax(npt,K2),mimicsbiome%Vmax(npt,K3), &
                   mimicsbiome%Km(npt,R1),mimicsbiome%Km(npt,R2),mimicsbiome%Km(npt,R3), &
                   mimicsbiome%Km(npt,K1),mimicsbiome%Km(npt,K2),mimicsbiome%Km(npt,K3), &
                   mimicsbiome%Vslope(R1),mimicsbiome%Vslope(R2),mimicsbiome%Vslope(R3), &
                   mimicsbiome%Vslope(K1),mimicsbiome%Vslope(K2),mimicsbiome%Vslope(K3), &
                   mimicsbiome%Kslope(R1),mimicsbiome%Kslope(R2),mimicsbiome%Kslope(R3), &
                   mimicsbiome%Kslope(K1),mimicsbiome%Kslope(K2),mimicsbiome%Kslope(K3), &
                   mimicsbiome%Vint(R1),mimicsbiome%Kint(R1),Tsoil,casamet%moistavg(npt), &

                   casaflux%CnppAn(npt),casapool%Clitter(npt,cwd),mimicsbiome%ligninNratioAvg(npt), &
                   cleaf2met(npt),cleaf2str(npt),croot2met(npt),croot2str(npt), &
                   cwd2str(npt),cwd2co2(npt),cwood2cwd(npt), &

                   mimicsbiome%fAVAL(npt,1), mimicsbiome%fAVAL(npt,2), mimicsbiome%fCHEM(npt,1), &
                   mimicsbiome%fCHEM(npt,2), mimicsbiome%fPHYS(npt,1), mimicsbiome%fPHYS(npt,2), &

                   mimicsbiome%Kmod(npt,R1),mimicsbiome%Kmod(npt,R2),mimicsbiome%Kmod(npt,R3), &
                   mimicsbiome%Kmod(npt,K1),mimicsbiome%Kmod(npt,K2),mimicsbiome%Kmod(npt,K3), &
                   mimicsflux%CSOMpInput(npt), Cbalance


    close(unit1)

102  format(4(i6,','),19(f18.8,','),15(f18.8,','),31(f18.8,','),24(f18.8,','))

END SUBROUTINE WritePointMIMICS

!-------------------------------------------------------------------------------------
! Write MIMICS CN pools and other quantities to .csv output daily

SUBROUTINE WritePointMIMICS_CN(unit1, sPtFileName, npt, mp, iYrCnt, idoy, &
    cleaf2met,cleaf2str,croot2met,croot2str,cwd2str,cwd2co2,cwood2cwd, &
    LITmin, MICtrn, SOMmin, DEsorb, OXIDAT, &
    dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil, &

    nleaf2met,nleaf2str,nroot2met,nroot2str,nwd2str,nwood2cwd, &
    LITminN, MICtrnN, SOMminN, DEsorbN, OXIDATN, &
    dLITmN, dLITsN, dMICrN, dMICkN, dSOMaN, dSOMcN, dSOMpN, dDIN, &
    DINup_r, DINup_k, upMICrC, upMICrN, upMICkC, upMICkN, &
    Overflow_r, Overflow_k, Nspill_r, Nspill_k, Cbalance, Nbalance) 

    USE define_types
    USE casadimension
    USE casaparm
    USE casavariable
    USE clmgridvariable
    USE mimicsdimension
    USE mimicsparam
    USE mimicsvariable
    implicit none

!   ARGUMENTS
    integer, intent(in)            :: unit1         ! FORTRAN file unit
    character(len=*), intent(in)   :: sPtFileName   ! .csv output file name 
    integer, intent(IN)            :: mp, npt       ! # points, point index
    integer, intent(IN)            :: iYrCnt, idoy  ! simulation year count, day of year
    real(r_2), dimension(mp),intent(IN) :: cleaf2met,cleaf2str,croot2met,croot2str
    real(r_2), dimension(mp),intent(IN) :: cwd2str,cwd2co2,cwood2cwd
    real(r_2), dimension(mp),intent(IN) :: nleaf2met,nleaf2str,nroot2met,nroot2str
    real(r_2), dimension(mp),intent(IN) :: nwd2str,nwood2cwd
    real(r_2), intent(IN)          :: LITmin(4), MICtrn(6), SOMmin(2), DEsorb, OXIDAT
    real(r_2), intent(IN)          :: LITminN(4), MICtrnN(6), SOMminN(2), DEsorbN, OXIDATN
    real(r_2), intent(IN)          :: dLITm, dLITs, dSOMa, dSOMc, dSOMp, dMICr, dMICk, Tsoil
    real(r_2), intent(IN)          :: dLITmN, dLITsN, dMICrN, dMICkN, dSOMaN, dSOMcN, dSOMpN, dDIN
    real(r_2), intent(IN)          :: DINup_r, DINup_k, upMICrC, upMICrN, upMICkC, upMICkN 
    real(r_2), intent(IN)          :: Overflow_r, Overflow_k, Nspill_r, Nspill_k, Cbalance, Nbalance


    open(unit1,file=sPtFileNameMIMICS, access='APPEND')

    write(unit1,102) npt,casamet%ijgcm(npt),iYrCnt,idoy,casamet%tsoilavg(npt),casamet%moistavg(npt),&
                   mimicsflux%Chresp(npt),mimicsflux%ClitInput(npt,metbc),mimicsflux%ClitInput(npt,struc), &
                   mimicspool%LITm(npt),mimicspool%LITs(npt),mimicspool%MICr(npt), &
                   mimicspool%MICk(npt),mimicspool%SOMa(npt),mimicspool%SOMc(npt),mimicspool%SOMp(npt), &
                   dLITm, dLITs, dMICr, dMICk, dSOMa, dSOMc, dSOMp,&

                   LITmin(1), LITmin(2), LITmin(3), LITmin(4), &
                   MICtrn(1), MICtrn(2), MICtrn(3), MICtrn(4), MICtrn(5), MICtrn(6), &
                   SOMmin(1), SOMmin(2), DEsorb, OXIDAT, mimicsbiome%fmet(npt), &

                   mimicsbiome%tauMod(npt), mimicsbiome%tauR(npt), mimicsbiome%tauK(npt), &
                   mimicsbiome%Vmax(npt,R1),mimicsbiome%Vmax(npt,R2),mimicsbiome%Vmax(npt,R3), &
                   mimicsbiome%Vmax(npt,K1),mimicsbiome%Vmax(npt,K2),mimicsbiome%Vmax(npt,K3), &
                   mimicsbiome%Km(npt,R1),mimicsbiome%Km(npt,R2),mimicsbiome%Km(npt,R3), &
                   mimicsbiome%Km(npt,K1),mimicsbiome%Km(npt,K2),mimicsbiome%Km(npt,K3), &
                   mimicsbiome%Vslope(R1),mimicsbiome%Vslope(R2),mimicsbiome%Vslope(R3), &
                   mimicsbiome%Vslope(K1),mimicsbiome%Vslope(K2),mimicsbiome%Vslope(K3), &
                   mimicsbiome%Kslope(R1),mimicsbiome%Kslope(R2),mimicsbiome%Kslope(R3), &
                   mimicsbiome%Kslope(K1),mimicsbiome%Kslope(K2),mimicsbiome%Kslope(K3), &
                   mimicsbiome%Vint(R1),mimicsbiome%Kint(R1),Tsoil,casamet%moistavg(npt), &

                   casaflux%CnppAn(npt),casapool%Clitter(npt,cwd),mimicsbiome%ligninNratioAvg(npt), &
                   cleaf2met(npt),cleaf2str(npt),croot2met(npt),croot2str(npt), &
                   cwd2str(npt),cwd2co2(npt),cwood2cwd(npt), &

                   mimicsbiome%fAVAL(npt,1), mimicsbiome%fAVAL(npt,2), mimicsbiome%fCHEM(npt,1), &
                   mimicsbiome%fCHEM(npt,2), mimicsbiome%fPHYS(npt,1), mimicsbiome%fPHYS(npt,2), &

                   mimicsbiome%Kmod(npt,R1),mimicsbiome%Kmod(npt,R2),mimicsbiome%Kmod(npt,R3), &
                   mimicsbiome%Kmod(npt,K1),mimicsbiome%Kmod(npt,K2),mimicsbiome%Kmod(npt,K3), &
                   mimicsflux%CSOMpInput(npt), &

                   mimicsflux%NlitInput(npt,metbc), mimicsflux%NlitInput(npt,struc), &
                   mimicspool%LITmN(npt), mimicspool%LITsN(npt), mimicspool%MICrN(npt), &
                   mimicspool%MICkN(npt), mimicspool%SOMaN(npt), mimicspool%SOMcN(npt), mimicspool%SOMpN(npt), &
                   mimicspool%DIN(npt), &
                   dLITmN, dLITsN, dMICrN, dMICkN, dSOMaN, dSOMcN, dSOMpN, dDIN, &

                   LITminN(1), LITminN(2), LITminN(3), LITminN(4), &
                   MICtrnN(1), MICtrnN(2), MICtrnN(3), MICtrnN(4), MICtrnN(5), MICtrnN(6), &
                   SOMminN(1), SOMminN(2), DEsorbN, OXIDATN, &
                   DINup_r, DINup_k, upMICrC, upMICrN, upMICkC, upMICkN, &
                   Overflow_r, Overflow_k, Nspill_r, Nspill_k, Cbalance, Nbalance, &
                   nleaf2met(npt), nleaf2str(npt), nroot2met(npt), nroot2str(npt), nwd2str(npt), nwood2cwd(npt)

                   

    close(unit1)

102  format(4(i6,','),19(f18.10,','),15(f18.10,','),31(f18.10,','),23(f10.4,','), 18(f12.6,','),32(f12.6,','))

END SUBROUTINE WritePointMIMICS_CN

!-------------------------------------------------------------------------------------
