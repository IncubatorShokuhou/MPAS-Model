! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
! Column mask for hierarchical time stepping. column_mask_on stays false unless
! a multi-level step is active, so the its:ite loops keep their original range.

module mpas_atmphys_column_mask

   use mpas_kind_types, only: RKIND

   implicit none
   public

   logical :: column_mask_on = .false.
   logical, allocatable :: column_active(:)

   ! Multi-level HTS moves the dt-integrating surface updates onto the level
   ! step. One level and HTS-off leave this false, and the time tracker keeps
   ! the original calls.
   logical :: column_mask_defer_sfc_dt = .false.
   logical :: column_sfc_leap = .false.
   real(kind=RKIND) :: column_sfc_julday = 0.0_RKIND

   contains

   logical function column_is_active(i)

      ! True for every index when the mask is off, so the its:ite loops keep
      ! their original columns. The mask array is not read on that path.

      integer, intent(in) :: i

      if (.not. column_mask_on) then
         column_is_active = .true.
         return
      end if
      if (.not. allocated(column_active)) then
         column_is_active = .false.
         return
      end if
      if (i < lbound(column_active, 1) .or. i > ubound(column_active, 1)) then
         column_is_active = .false.
         return
      end if
      column_is_active = column_active(i)

   end function column_is_active

   subroutine column_mask_next_segment(cursor, ite, segStart, segEnd, found)

      ! Next contiguous run of columns in [cursor, ite]. With the mask off
      ! that run is the whole range, one call, the same its:ite as before.

      implicit none

      integer, intent(inout) :: cursor
      integer, intent(in) :: ite
      integer, intent(out) :: segStart, segEnd
      logical, intent(out) :: found

      integer :: i

      if (cursor > ite) then
         found = .false.
         return
      end if

      if (.not. column_mask_on) then
         segStart = cursor
         segEnd = ite
         cursor = ite + 1
         found = .true.
         return
      end if

      i = cursor
      do while (i <= ite)
         if (i >= lbound(column_active, 1) .and. i <= ubound(column_active, 1)) then
            if (column_active(i)) exit
         end if
         i = i + 1
      end do
      if (i > ite) then
         found = .false.
         cursor = ite + 1
         return
      end if

      segStart = i
      i = i + 1
      do while (i <= ite)
         if (i < lbound(column_active, 1) .or. i > ubound(column_active, 1)) exit
         if (.not. column_active(i)) exit
         i = i + 1
      end do
      segEnd = i - 1
      cursor = i
      found = .true.

   end subroutine column_mask_next_segment

end module mpas_atmphys_column_mask
