MODULE define_dimensions
  INTEGER            :: mp      ! # points in this simulation
  INTEGER, PARAMETER :: mf = 2  ! # leaves (sunlit, shaded)
  INTEGER, PARAMETER :: nrb = 3 ! # radiation bands
  INTEGER, PARAMETER :: ms = 6  ! # soil layers
  INTEGER, PARAMETER :: ncp = 3 ! # vegetation carbon stores
  INTEGER, PARAMETER :: ncs = 2 ! # soil carbon stores
END MODULE define_dimensions

MODULE define_types
  ! i_d is default kind for representing integer values.
!  INTEGER, PARAMETER :: i_d = KIND(9)
  ! r_1 is default kind for representing REAL values (typically 32 or 64 bits).
!  INTEGER, PARAMETER :: r_1  = KIND(1.0)
  ! r_2 is kind for representing REAL values with at least 10-digit precision
  ! (typically 64 bits).
!  INTEGER, PARAMETER :: r_2  = SELECTED_REAL_KIND(12, 50)

 TYPE veg_parameter_type
      integer, dimension(:), pointer   :: iveg
      real,    dimension(:,:), pointer :: froot
 END type veg_parameter_type

 TYPE soil_parameter_type
      integer, dimension(:), pointer  :: isoilm
      real,    dimension(:), pointer  :: sfc, swilt
      real,    dimension(:), pointer  :: clay, silt, sand, ssat
      real,    dimension(:), pointer  :: dzsoil
 END type soil_parameter_type

 TYPE (veg_parameter_type)  :: veg
 TYPE (soil_parameter_type) :: soil

Contains

   SUBROUTINE alloc_casavegsoil(mp,ms)
     IMPLICIT none

! TYPE (veg_parameter_type)  :: veg
! TYPE (soil_parameter_type) :: soil

     integer, INTENT(IN) :: mp,ms
  
     ALLOCATE(veg%iveg(mp))
     ALLOCATE(veg%froot(mp,ms))
     ALLOCATE(soil%isoilm(mp))
     ALLOCATE(soil%sfc(mp), soil%swilt(mp), soil%ssat(mp))
     ALLOCATE(soil%clay(mp), soil%silt(mp), soil%sand(mp))
     ALLOCATE(soil%dzsoil(ms))
   END subroutine alloc_casavegsoil

END MODULE define_types

