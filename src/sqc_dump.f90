! sqcalc - structure factors from LAMMPS dump trajectories
! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
!
! SPDX-License-Identifier: GPL-3.0-or-later

!> Streaming reader for LAMMPS dump files (default "atoms" dump style).
module sqc_dump
   use sqc_kinds
   use sqc_cell, only: cell_t
   implicit none
   private

   public :: frame_t, frame_reader_t, lammps_dump_reader_t, max_line_length

   integer, parameter :: max_line_length = 1024
   integer, parameter :: max_columns = 64
   !> Sentinel for "no unit is associated with this reader" (NEWUNIT may return
   !! negative unit numbers, so zero or a sign test would be wrong).
   integer, parameter :: unopened_unit = -huge(1)

   !> A single trajectory snapshot.
   type :: frame_t
      !> Number of atoms in the frame.
      integer(ik) :: natoms = 0
      !> Time step as written in the dump.
      integer(ik) :: timestep = 0
      !> Cartesian coordinates, pos(:, i) for atom i.
      real(rk), allocatable :: pos(:, :)
      !> LAMMPS atom id of each atom (only when the dump has an "id" column).
      integer(ik), allocatable :: id(:)
      !> LAMMPS atom type of each atom (1-based).
      integer(ik), allocatable :: type_id(:)
      !> Simulation cell of this frame.
      type(cell_t) :: cell
      !> Periodicity of the box as written in the dump (pp = periodic).
      logical :: pbc(3) = .true.
      !> True when pos holds unwrapped coordinates: either the dump itself
      !! provided xu/yu/zu, or the reader reconstructed them from wrapped
      !! x/y/z by tracking the image each atom crossed.
      logical :: unwrapped = .false.
      !> True when the unwrapped coordinates came straight from the dump.
      logical :: unwrapped_from_dump = .false.
   contains
      procedure :: destroy => frame_destroy
   end type frame_t

   !> Abstract interface for trajectory readers.
   type, abstract :: frame_reader_t
      !> Number of frames delivered so far.
      integer(ik) :: nframes = 0
   contains
      procedure(frame_next_iface), deferred :: next_frame
      procedure(frame_close_iface), deferred :: close
   end type frame_reader_t

   abstract interface
      !> Read the next frame into frame.  ierr is 0 on success and 1 at the end
      !! of the trajectory; larger values are fatal errors.
      subroutine frame_next_iface(self, frame, ierr)
         import :: frame_reader_t, frame_t
         class(frame_reader_t), intent(inout) :: self
         type(frame_t), intent(inout) :: frame
         integer, intent(out) :: ierr
      end subroutine frame_next_iface

      subroutine frame_close_iface(self)
         import :: frame_reader_t
         class(frame_reader_t), intent(inout) :: self
      end subroutine frame_close_iface
   end interface

   !> Reader for the default LAMMPS atoms dump style.
   type, extends(frame_reader_t) :: lammps_dump_reader_t
      character(len=:), allocatable :: path
      integer :: unit = unopened_unit
      !> Column positions of the quantities we need (0 = not present).
      integer :: idx_id = 0, idx_type = 0, idx_x = 0, idx_y = 0, idx_z = 0
      integer :: ncol = 0
      real(rk), allocatable :: buffer(:)
      !> The x/y/z columns were the unwrapped xu/yu/zu variants.
      logical :: from_unwrapped = .false.
      !> Unwrapping state: the wrapped fractional coordinates of the previous
      !! frame and the accumulated (unwrapped) fractional coordinates.
      real(rk), allocatable :: prev_frac(:, :), unwrapped_frac(:, :)
      logical :: have_prev = .false.
      !> Atom ids of the previous frame (used to detect a reordering that would
      !! make the index based unwrapping invalid).
      integer(ik), allocatable :: prev_id(:)
      !> Frames whose displacement was so large that unwrapping is ambiguous.
      integer :: n_ambiguous = 0
      !> Timestep bookkeeping: the increment must be constant, it fixes the
      !! physical time between frames together with the --dt time step.
      integer(ik) :: prev_timestep = 0, step_stride = 0
      logical :: have_timestep = .false., uniform_timestep = .true.
      !> True while the next expected line is "ITEM: TIMESTEP".
      logical :: expect_header = .true.
      !> Line number of the line last read (for error messages).
      integer :: lineno = 0
      !> Whether the current frame used the triclinic box convention.
      logical :: triclinic = .false.
   contains
      procedure :: open => lammps_open
      procedure :: next_frame => lammps_next_frame
      procedure :: close => lammps_close
      final :: lammps_final
   end type lammps_dump_reader_t

contains

   subroutine frame_destroy(self)
      class(frame_t), intent(inout) :: self
      if (allocated(self%pos)) deallocate(self%pos)
      if (allocated(self%id)) deallocate(self%id)
      if (allocated(self%type_id)) deallocate(self%type_id)
      self%natoms = 0
      self%timestep = 0
      self%unwrapped = .false.
      self%unwrapped_from_dump = .false.
   end subroutine frame_destroy

   !> Open a dump file and position the reader on the first frame header.
   subroutine lammps_open(self, path, ierr)
      class(lammps_dump_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: path
      integer, intent(out) :: ierr
      character(len=max_line_length) :: line
      integer :: ios

      self%path = trim(path)
      open(newunit=self%unit, file=self%path, status='old', action='read', &
           form='formatted', iostat=ios)
      if (ios /= 0) then
         ierr = 2
         return
      end if
      self%lineno = 0
      self%expect_header = .true.
      self%nframes = 0
      self%have_prev = .false.
      self%have_timestep = .false.
      self%step_stride = 0
      self%uniform_timestep = .true.
      self%n_ambiguous = 0
      ierr = 0
   end subroutine lammps_open

   subroutine lammps_close(self)
      class(lammps_dump_reader_t), intent(inout) :: self
      if (self%unit /= unopened_unit) then
         close(self%unit)
         self%unit = unopened_unit
      end if
   end subroutine lammps_close

   subroutine lammps_final(self)
      type(lammps_dump_reader_t), intent(inout) :: self
      call lammps_close(self)
      if (allocated(self%buffer)) deallocate(self%buffer)
      self%expect_header = .true.
   end subroutine lammps_final

   !> Read the next frame; ierr = 1 signals a clean end of file.
   subroutine lammps_next_frame(self, frame, ierr)
      class(lammps_dump_reader_t), intent(inout) :: self
      type(frame_t), intent(inout) :: frame
      integer, intent(out) :: ierr
      character(len=max_line_length) :: line, label
      character(len=16) :: tokens(max_columns)
      real(rk) :: bounds(3, 3), values(3)
      integer :: ios, ntok, n, i, irow, nvals
      logical :: have_atoms, have_box, triclinic

      ierr = 0
      if (self%unit == unopened_unit) then
         ierr = 3
         return
      end if

      ! Position the reader on the "ITEM: TIMESTEP" line of the next frame.
      if (self%expect_header) then
         call read_line(self, line, ios)
         if (ios /= 0) then
            ierr = 1
            return
         end if
      else
         do
            call read_line(self, line, ios)
            if (ios /= 0) then
               ierr = 1
               return
            end if
            if (is_item(line, 'TIMESTEP')) exit
         end do
      end if
      if (.not. is_item(line, 'TIMESTEP')) then
         ierr = 4
         return
      end if
      self%expect_header = .false.

      call read_line(self, line, ios)
      if (ios /= 0) then
         ierr = 4
         return
      end if
      read (line, *, iostat=ios) frame%timestep
      if (ios /= 0) then
         ierr = 4
         return
      end if

      ! Walk through the remaining header records of this frame.
      have_atoms = .false.
      have_box = .false.
      triclinic = .false.
      bounds = 0.0_rk
      do
         call read_line(self, line, ios)
         if (ios /= 0) then
            ierr = 4
            return
         end if
         if (is_item(line, 'NUMBER OF ATOMS')) then
            call read_line(self, line, ios)
            if (ios /= 0) then
               ierr = 4
               return
            end if
            read (line, *, iostat=ios) frame%natoms
            if (ios /= 0 .or. frame%natoms <= 0) then
               ierr = 5
               return
            end if
         else if (is_item(line, 'BOX BOUNDS')) then
            triclinic = has_tilt_header(line)
            call read_pbc(line, frame%pbc)
            do irow = 1, 3
               call read_line(self, line, ios)
               if (ios /= 0) then
                  ierr = 4
                  return
               end if
               values = 0.0_rk
               nvals = 3
               read (line, *, iostat=ios) values(1:nvals)
               if (ios /= 0) then
                  nvals = 2
                  read (line, *, iostat=ios) values(1:nvals)
                  if (ios /= 0) then
                     ierr = 6
                     return
                  end if
               end if
               bounds(:, irow) = values
            end do
            have_box = .true.
         else if (is_item(line, 'ATOMS')) then
            label = adjustl(line(len('ITEM: ATOMS') + 1:))
            call split_tokens(label, tokens, ntok)
            call resolve_columns(self, tokens, ntok, ierr)
            if (ierr /= 0) return
            have_atoms = .true.
            exit
         end if
      end do

      if (.not. have_atoms .or. .not. have_box) then
         ierr = 7
         return
      end if

      call frame_alloc_atoms(frame, frame%natoms)
      if (.not. allocated(self%buffer)) allocate(self%buffer(self%ncol))

      do i = 1, frame%natoms
         read (self%unit, *, iostat=ios) self%buffer
         if (ios /= 0) then
            ierr = 8
            return
         end if
         frame%pos(1, i) = self%buffer(self%idx_x)
         frame%pos(2, i) = self%buffer(self%idx_y)
         frame%pos(3, i) = self%buffer(self%idx_z)
         frame%type_id(i) = nint(self%buffer(self%idx_type))
         if (self%idx_id > 0) then
            frame%id(i) = nint(self%buffer(self%idx_id))
         else
            frame%id(i) = int(i, ik)
         end if
      end do

      call frame%cell%init_from_bounds(bounds(1, :), bounds(2, :), &
                                      [bounds(3, 1), bounds(3, 2), bounds(3, 3)])
      self%triclinic = triclinic

      call unwrap_positions(self, frame, ierr)
      if (ierr /= 0) return

      ! The time step between dumped frames fixes the physical time axis; the
      ! dynamics method needs it to be constant.
      if (.not. self%have_timestep) then
         self%have_timestep = .true.
      else if (self%step_stride == 0) then
         self%step_stride = frame%timestep - self%prev_timestep
      else if (frame%timestep - self%prev_timestep /= self%step_stride) then
         self%uniform_timestep = .false.
      end if
      self%prev_timestep = frame%timestep

      self%nframes = self%nframes + 1
      self%expect_header = .true.
      ierr = 0
   end subroutine lammps_next_frame

   !> Make the frame coordinates unwrapped.
   !!
   !! An arbitrary (off lattice) q picks up the phase of every box crossing, so
   !! the dynamic structure factor needs continuous coordinates.  When the dump
   !! already carries xu/yu/zu nothing has to be done; otherwise the image each
   !! atom crossed since the previous frame is accumulated from the
   !! minimum-image displacement in fractional coordinates (the same trick
   !! `sqc_neighbour_list` uses for its Verlet list).  The reconstruction is
   !! index based, so a reordering of the atoms (`ITEM: ATOMS id` changing
   !! between frames) is rejected.
   subroutine unwrap_positions(self, frame, ierr)
      type(lammps_dump_reader_t), intent(inout) :: self
      type(frame_t), intent(inout) :: frame
      integer, intent(out) :: ierr
      real(rk) :: s(3), ds(3)
      integer :: i, k
      logical :: ambiguous

      ierr = 0
      if (self%from_unwrapped) then
         frame%unwrapped = .true.
         frame%unwrapped_from_dump = .true.
         return
      end if

      ! An atom count change cannot be unwrapped across (the state is indexed
      ! by atom): report it so the caller can complain about the dump.
      if (self%have_prev) then
         if (size(self%prev_frac, 2) /= frame%natoms) then
            ierr = 11
            return
         end if
      end if

      ! Unwrapping is index based: the atom order must not change.
      if (self%idx_id > 0) then
         if (allocated(self%prev_id)) then
            if (size(self%prev_id) /= frame%natoms .or. any(self%prev_id /= frame%id)) then
               ierr = 10
               return
            end if
         else
            allocate (self%prev_id(frame%natoms))
         end if
         self%prev_id = frame%id
      end if

      if (.not. self%have_prev) then
         allocate (self%prev_frac(3, frame%natoms), self%unwrapped_frac(3, frame%natoms))
         do i = 1, frame%natoms
            s = wrap_periodic(frame%cell%fractional(frame%pos(:, i)), frame%pbc)
            self%prev_frac(:, i) = s
            self%unwrapped_frac(:, i) = s
         end do
         self%have_prev = .true.
      else
         do i = 1, frame%natoms
            s = wrap_periodic(frame%cell%fractional(frame%pos(:, i)), frame%pbc)
            ds = s - self%prev_frac(:, i)
            ambiguous = .false.
            do k = 1, 3
               if (.not. frame%pbc(k)) cycle
               ds(k) = ds(k) - anint(ds(k))
               if (abs(ds(k)) > 0.4_rk) ambiguous = .true.
            end do
            if (ambiguous) self%n_ambiguous = self%n_ambiguous + 1
            self%prev_frac(:, i) = s
            self%unwrapped_frac(:, i) = self%unwrapped_frac(:, i) + ds
         end do
      end if

      do i = 1, frame%natoms
         frame%pos(:, i) = matmul(frame%cell%a, self%unwrapped_frac(:, i))
      end do
      frame%unwrapped = .true.
      frame%unwrapped_from_dump = .false.
   end subroutine unwrap_positions

   !> Wrap fractional coordinates into [0, 1) along the periodic directions
   !! only.  A non-periodic direction has no lattice to wrap by, and the
   !! coordinate the dump writes for it (x, not xu) is already the true one.
   pure function wrap_periodic(s, pbc) result(wrapped)
      real(rk), intent(in) :: s(3)
      logical, intent(in) :: pbc(3)
      real(rk) :: wrapped(3)
      integer :: k

      wrapped = s
      do k = 1, 3
         if (pbc(k)) wrapped(k) = wrapped(k) - floor(wrapped(k))
      end do
   end function wrap_periodic

   !> (Re)allocate the per atom arrays of a frame for a given atom count.
   !!
   !! The methods require a constant atom count, but the reader has to survive
   !! a frame that breaks it: reallocating here keeps the atom loop in bounds
   !! (and the unwrapping state check below detects the change) so that the
   !! caller reports the mismatch instead of the reader overrunning memory.
   subroutine frame_alloc_atoms(frame, natoms)
      type(frame_t), intent(inout) :: frame
      integer(ik), intent(in) :: natoms

      if (allocated(frame%pos)) then
         if (size(frame%pos, 2) /= natoms) deallocate (frame%pos)
      end if
      if (.not. allocated(frame%pos)) allocate (frame%pos(3, natoms))

      if (allocated(frame%id)) then
         if (size(frame%id) /= natoms) deallocate (frame%id)
      end if
      if (.not. allocated(frame%id)) allocate (frame%id(natoms))

      if (allocated(frame%type_id)) then
         if (size(frame%type_id) /= natoms) deallocate (frame%type_id)
      end if
      if (.not. allocated(frame%type_id)) allocate (frame%type_id(natoms))
   end subroutine frame_alloc_atoms

   subroutine read_line(self, line, ios)
      type(lammps_dump_reader_t), intent(inout) :: self
      character(len=max_line_length), intent(out) :: line
      integer, intent(out) :: ios
      read (self%unit, '(A)', iostat=ios) line
      if (ios == 0) self%lineno = self%lineno + 1
   end subroutine read_line

   !> True if the record is the LAMMPS "ITEM: <label> ..." header for label.
   pure logical function is_item(line, label) result(found)
      character(len=*), intent(in) :: line, label
      character(len=:), allocatable :: trimmed
      trimmed = adjustl(line)
      found = index(trimmed, 'ITEM: '//label) == 1
   end function is_item

   !> Detect the "xy xz yz" tilt keywords on a BOX BOUNDS header line.
   pure logical function has_tilt_header(line) result(triclinic)
      character(len=*), intent(in) :: line
      triclinic = index(line, 'xy') > 0 .or. index(line, 'xy xz yz') > 0
   end function has_tilt_header

   !> Read the periodicity tokens at the end of a BOX BOUNDS header line
   !! ("pp", "ff" or "ss" per direction).
   subroutine read_pbc(line, pbc)
      character(len=*), intent(in) :: line
      logical, intent(out) :: pbc(3)
      character(len=16) :: tokens(max_columns)
      integer :: ntok, i

      pbc = .true.
      call split_tokens(adjustl(line(len('ITEM: BOX BOUNDS') + 1:)), tokens, ntok)
      if (ntok < 3) return
      do i = 1, 3
         pbc(i) = tokens(ntok - 3 + i)(1:1) == 'p'
      end do
   end subroutine read_pbc

   !> Split a record into whitespace separated tokens.
   subroutine split_tokens(line, tokens, ntok)
      character(len=*), intent(in) :: line
      character(len=*), intent(out) :: tokens(:)
      integer, intent(out) :: ntok
      integer :: pos, start, n

      ntok = 0
      n = len_trim(line)
      pos = 1
      do while (pos <= n .and. ntok < size(tokens))
         ! Fortran does not short-circuit .and., so index tests must be explicit
         do
            if (pos > n) exit
            if (line(pos:pos) /= ' ') exit
            pos = pos + 1
         end do
         if (pos > n) exit
         start = pos
         do
            if (pos > n) exit
            if (line(pos:pos) == ' ') exit
            pos = pos + 1
         end do
         ntok = ntok + 1
         tokens(ntok) = line(start:pos - 1)
      end do
   end subroutine split_tokens

   !> Locate the atom record columns we need.
   subroutine resolve_columns(self, tokens, ntok, ierr)
      type(lammps_dump_reader_t), intent(inout) :: self
      character(len=*), intent(in) :: tokens(:)
      integer, intent(in) :: ntok
      integer, intent(out) :: ierr
      integer :: i
      character(len=16) :: name

      self%ncol = ntok
      self%idx_id = 0
      self%idx_type = 0
      self%idx_x = 0
      self%idx_y = 0
      self%idx_z = 0
      self%from_unwrapped = .false.
      block
         logical :: unwrapped_x, unwrapped_y, unwrapped_z
         unwrapped_x = .false.
         unwrapped_y = .false.
         unwrapped_z = .false.
         do i = 1, ntok
            name = lowercase(trim(tokens(i)))
            select case (name)
            case ('id')
               self%idx_id = i
            case ('type')
               self%idx_type = i
            case ('x', 'xu')
               if (self%idx_x == 0) then
                  self%idx_x = i
                  unwrapped_x = name == 'xu'
               end if
            case ('y', 'yu')
               if (self%idx_y == 0) then
                  self%idx_y = i
                  unwrapped_y = name == 'yu'
               end if
            case ('z', 'zu')
               if (self%idx_z == 0) then
                  self%idx_z = i
                  unwrapped_z = name == 'zu'
               end if
            end select
         end do
         self%from_unwrapped = unwrapped_x .and. unwrapped_y .and. unwrapped_z
      end block

      ierr = 0
      if (self%idx_type == 0 .or. self%idx_x == 0 .or. self%idx_y == 0 &
          .or. self%idx_z == 0) ierr = 9
   end subroutine resolve_columns

   pure function lowercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, code
      out = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) out(i:i) = achar(code + 32)
      end do
   end function lowercase

end module sqc_dump
