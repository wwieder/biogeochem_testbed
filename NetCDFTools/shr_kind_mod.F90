module shr_kind_mod

   !===========================================================================!
   ! Define "kind" variable for double precision
   !===========================================================================!

   integer, parameter :: r8 = selected_real_kind(12) ! kind=8 (8-byte real, equivalent to double precision)

end module shr_kind_mod
