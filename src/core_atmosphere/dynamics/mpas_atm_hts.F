! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
! Hierarchical time stepping for the MPAS-Atmosphere split-explicit integrator.
! See docs/hts/DESIGN.md. config_hts = false does not touch the prognostic state.

module mpas_atm_hts

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use mpas_constants
   use mpas_timer, only: mpas_timer_start, mpas_timer_stop

   implicit none

   public :: atm_hts_setup, atm_hts_finalize
   public :: atm_hts_begin_intervals, atm_hts_apply_interpolation, atm_hts_restore_above
   public :: atm_hts_sync_finer, atm_hts_accumulate_flux, atm_hts_prepare_targets
   public :: atm_hts_clear_consumed, atm_hts_correct_mass
   public :: hts_restrict_physics_tendencies
   public :: hts_exch_layer1
   public :: hts_diagnostics_from_state
   public :: atm_hts_mark_acoustic_cells
   public :: hts_cell_acoustic

   logical, public :: hts_enabled = .false.
   logical, public :: hts_masking = .false.
   logical, public :: hts_impose = .false.
   integer, public :: hts_active_level = 0
   integer, public :: hts_tick = 0
   integer, public :: hts_max_level_used = 0

   integer, public, allocatable :: hts_cell_level(:)
   integer, public, allocatable :: hts_edge_level(:)
   logical, public, allocatable :: hts_edge_interface(:)
   ! True when a cell has an edge finer than itself. The cell stays on the
   ! interval snapshot during the finer step. Advancing only its acoustic
   ! rho_pp and rtheta_pp (bb575a1b2) NaNed u at 00:18 on the eds land VR
   ! (job 199953), earlier than the pure hold.
   logical, public, allocatable :: hts_cell_open_iface(:)
   ! Acoustic columns of the active level. Open-interface parents are not
   ! included. hts_cell_acoustic does not read this mask.
   logical, public, allocatable :: hts_cell_live(:)
   logical, public, allocatable :: hts_impose_edge(:)
   real (kind=RKIND), public, allocatable :: hts_ru_target(:,:)

   integer, parameter :: hts_level_cap = 8

   ! Index lists. With one thread the hot loops iterate these instead of every
   ! cell. Work lists are the active level plus a 2-cell ring, because the
   ! tendency and hyperdiffusion stencils read neighbors. Acoustic state updates
   ! use the strict active list. hts_ac_edge also contains the interface edges
   ! a coarse level imposes.
   integer, public :: hts_nloop = 0
   integer, public :: hts_ii = 0
   ! Tendency and dissipation index lists. hts_lists_work false walks the
   ! original index range. Diagnostics use a wider list (hts_diag_rings),
   ! two rings past the tendency work list, so the APVM and kinetic-energy
   ! values a level reads are the full-domain values.
   logical, public :: hts_lists_work = .true.
   integer, parameter :: hts_work_rings = 2
   integer, parameter :: hts_diag_rings = 4
   integer, public :: hts_nCellsSolve = 0
   integer, public :: hts_nEdgesSolve = 0
   integer, public :: hts_nVertices = 0
   integer, public :: hts_nVerticesSolve = 0

   integer, public :: hts_n_work_own(0:hts_level_cap) = 0
   integer, public :: hts_n_work_all(0:hts_level_cap) = 0
   integer, public :: hts_n_active_own(0:hts_level_cap) = 0
   integer, public :: hts_n_active_all(0:hts_level_cap) = 0
   integer, public :: hts_n_work_edge_own(0:hts_level_cap) = 0
   integer, public :: hts_n_work_edge_all(0:hts_level_cap) = 0
   integer, public :: hts_n_ac_edge(0:hts_level_cap) = 0
   integer, public :: hts_n_off_own(0:hts_level_cap) = 0
   integer, public :: hts_n_vert_all(0:hts_level_cap) = 0
   integer, public :: hts_n_vert_own(0:hts_level_cap) = 0
   integer, public :: hts_n_diag_all(0:hts_level_cap) = 0
   integer, public :: hts_n_diag_edge_all(0:hts_level_cap) = 0
   integer, public :: hts_n_diag_vert_all(0:hts_level_cap) = 0
   integer, public :: hts_n_level_edge(0:hts_level_cap) = 0
   integer, public :: hts_n_ac_cell(0:hts_level_cap) = 0
   integer, public :: hts_n_recon_own(0:hts_level_cap) = 0
   integer, public :: hts_n_live_own(0:hts_level_cap) = 0
   integer, public :: hts_n_live_all(0:hts_level_cap) = 0
   ! True once a level's state arrays already hold its snapshot. The next
   ! fine step does not rewrite that level. A level step clears the flag so
   ! the following fine step copies the new snapshot into both time levels.
   logical, public :: hts_hold_valid(0:hts_level_cap) = .false.
   logical, public :: hts_hold_rewrote(0:hts_level_cap) = .false.
   integer, public, allocatable :: hts_work_own(:,:)
   integer, public, allocatable :: hts_work_all(:,:)
   integer, public, allocatable :: hts_active_own(:,:)
   integer, public, allocatable :: hts_active_all(:,:)
   integer, public, allocatable :: hts_work_edge_own(:,:)
   integer, public, allocatable :: hts_work_edge_all(:,:)
   integer, public, allocatable :: hts_ac_edge(:,:)
   integer, public, allocatable :: hts_off_own(:,:)
   integer, public, allocatable :: hts_vert_all(:,:)
   integer, public, allocatable :: hts_vert_own(:,:)
   integer, public, allocatable :: hts_diag_all(:,:)
   integer, public, allocatable :: hts_diag_edge_all(:,:)
   integer, public, allocatable :: hts_diag_vert_all(:,:)
   integer, public, allocatable :: hts_level_edge(:,:)
   integer, public, allocatable :: hts_ac_cell(:,:)
   integer, public, allocatable :: hts_recon_own(:,:)
   integer, public, allocatable :: hts_live_own(:,:)
   integer, public, allocatable :: hts_live_all(:,:)

   ! Per-substep acoustic halo. rho_pp and rtheta_pp are read only on the
   ! two cells of an active edge, so the message carries cells that touch the
   ! active level. The non-HTS path does not call this exchange.
   logical, public :: hts_ah_ready = .false.
   integer :: hts_ah_nvert = 0
   integer :: hts_ah_nsend(0:hts_level_cap) = 0
   integer :: hts_ah_nrecv(0:hts_level_cap) = 0
   integer, allocatable :: hts_ah_send_proc(:,:)
   integer, allocatable :: hts_ah_send_n(:,:)
   integer, allocatable :: hts_ah_send_off(:,:)
   integer, allocatable :: hts_ah_send_cell(:)
   integer, allocatable :: hts_ah_recv_proc(:,:)
   integer, allocatable :: hts_ah_recv_n(:,:)
   integer, allocatable :: hts_ah_recv_off(:,:)
   integer, allocatable :: hts_ah_recv_cell(:)
   real (kind=RKIND), allocatable :: hts_ah_sbuf(:)
   real (kind=RKIND), allocatable :: hts_ah_rbuf(:)

   real (kind=R8KIND), allocatable :: hts_flux_sum(:,:)

   real (kind=RKIND), allocatable :: snap_rho(:,:,:), delta_rho(:,:,:)
   real (kind=RKIND), allocatable :: snap_theta(:,:,:), delta_theta(:,:,:)
   real (kind=RKIND), allocatable :: snap_w(:,:,:), delta_w(:,:,:)
   real (kind=RKIND), allocatable :: snap_rw(:,:,:), delta_rw(:,:,:)
   real (kind=RKIND), allocatable :: snap_u(:,:,:), delta_u(:,:,:)

   logical, allocatable :: hts_level_started(:)
   integer :: hts_nCells, hts_nEdges, hts_nVert
   integer :: hts_correct_logs
   real (kind=RKIND) :: hts_max_flux_resid

   contains


   subroutine atm_hts_setup(domain)

      implicit none

      type (domain_type), intent(inout) :: domain

      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: mesh
      logical, pointer :: config_hts
      integer, pointer :: config_hts_max_level
      integer, pointer :: nCells, nCellsSolve, nEdges, nEdgesSolve, nVertLevels
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      real (kind=RKIND), dimension(:), pointer :: meshDensity
      integer :: iCell, iEdge, iNbr, cell1, cell2, L, nbr
      integer :: nBlocks
      integer, dimension(0:hts_level_cap) :: localCount, globalCount
      integer :: localIface, globalIface, globalOwned, globalOpen
      real (kind=RKIND) :: mdMin, mdMax, mdMinG, mdMaxG, spacing, logtwo
      logical :: uniform

      hts_enabled = .false.
      hts_masking = .false.
      hts_impose = .false.
      hts_active_level = 0
      hts_tick = 0
      hts_max_level_used = 0
      hts_correct_logs = 0
      hts_max_flux_resid = 0.0_RKIND
      hts_hold_valid(:) = .false.
      hts_hold_rewrote(:) = .false.

      block => domain % blocklist
      nBlocks = 0
      do while (associated(block))
         nBlocks = nBlocks + 1
         block => block % next
      end do
      block => domain % blocklist
      if (nBlocks /= 1) then
         call mpas_log_write('HTS expects one block per MPI task; found $i.', &
                             intArgs=(/nBlocks/), messageType=MPAS_LOG_WARN)
      end if

      call mpas_pool_get_config(block % configs, 'config_hts', config_hts)
      call mpas_pool_get_config(block % configs, 'config_hts_max_level', config_hts_max_level)
      hts_enabled = config_hts

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(mesh, 'meshDensity', meshDensity)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)

      hts_nCells = nCells
      hts_nEdges = nEdges
      hts_nVert = nVertLevels
      hts_nCellsSolve = nCellsSolve
      hts_nEdgesSolve = nEdgesSolve

      if (allocated(hts_cell_level)) deallocate(hts_cell_level)
      if (allocated(hts_edge_level)) deallocate(hts_edge_level)
      if (allocated(hts_edge_interface)) deallocate(hts_edge_interface)
      if (allocated(hts_cell_open_iface)) deallocate(hts_cell_open_iface)
      if (allocated(hts_cell_live)) deallocate(hts_cell_live)
      if (allocated(hts_impose_edge)) deallocate(hts_impose_edge)
      if (allocated(hts_ru_target)) deallocate(hts_ru_target)
      if (allocated(hts_flux_sum)) deallocate(hts_flux_sum)

      allocate(hts_cell_level(nCells+1))
      allocate(hts_edge_level(nEdges+1))
      allocate(hts_edge_interface(nEdges+1))
      allocate(hts_cell_open_iface(nCells+1))
      allocate(hts_cell_live(nCells+1))
      allocate(hts_impose_edge(nEdges+1))
      hts_cell_level(:) = 0
      hts_edge_level(:) = 0
      hts_edge_interface(:) = .false.
      hts_cell_open_iface(:) = .false.
      hts_cell_live(:) = .false.
      hts_impose_edge(:) = .false.

      mdMin = meshDensity(1)
      mdMax = meshDensity(1)
      do iCell = 1, nCellsSolve
         mdMin = min(mdMin, meshDensity(iCell))
         mdMax = max(mdMax, meshDensity(iCell))
      end do
      call mpas_dmpar_min_real(domain % dminfo, mdMin, mdMinG)
      call mpas_dmpar_max_real(domain % dminfo, mdMax, mdMaxG)
      mdMin = mdMinG
      mdMax = mdMaxG
      uniform = (mdMin <= 0.0_RKIND) .or. (mdMax <= mdMin * (1.0_RKIND + 1.0e-5_RKIND))

      logtwo = log(2.0_RKIND)
      L = config_hts_max_level
      if (L < 0) L = 0
      if (L > hts_level_cap) L = hts_level_cap

      do iCell = 1, nCells
         if (uniform .or. meshDensity(iCell) <= 0.0_RKIND) then
            hts_cell_level(iCell) = 0
         else
            spacing = meshDensity(iCell) ** (-0.25_RKIND)
            if (spacing < 1.0_RKIND) spacing = 1.0_RKIND
            hts_cell_level(iCell) = nint(log(spacing) / logtwo)
            if (hts_cell_level(iCell) < 0) hts_cell_level(iCell) = 0
            if (hts_cell_level(iCell) > L) hts_cell_level(iCell) = L
         end if
      end do
      hts_cell_level(nCells+1) = 0

      ! A cell that touches a finer neighbor is pulled down to that neighbor's
      ! level only at the edge. Cell levels stay as binned. Edges take the finer cell.
      do iEdge = 1, nEdges
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         if (cell1 < 1 .or. cell1 > nCells) cell1 = nCells + 1
         if (cell2 < 1 .or. cell2 > nCells) cell2 = nCells + 1
         hts_edge_level(iEdge) = min(hts_cell_level(cell1), hts_cell_level(cell2))
         hts_edge_interface(iEdge) = hts_cell_level(cell1) /= hts_cell_level(cell2)
         if (cell1 >= 1 .and. cell1 <= nCells) then
            if (hts_edge_level(iEdge) < hts_cell_level(cell1)) hts_cell_open_iface(cell1) = .true.
         end if
         if (cell2 >= 1 .and. cell2 <= nCells) then
            if (hts_edge_level(iEdge) < hts_cell_level(cell2)) hts_cell_open_iface(cell2) = .true.
         end if
      end do

      localCount(:) = 0
      do iCell = 1, nCellsSolve
         localCount(hts_cell_level(iCell)) = localCount(hts_cell_level(iCell)) + 1
      end do
      do L = 0, hts_level_cap
         call mpas_dmpar_sum_int(domain % dminfo, localCount(L), globalCount(L))
      end do
      localIface = 0
      do iEdge = 1, nEdgesSolve
         if (hts_edge_interface(iEdge)) localIface = localIface + 1
      end do
      call mpas_dmpar_sum_int(domain % dminfo, localIface, globalIface)
      localIface = 0
      do iCell = 1, nCellsSolve
         if (hts_cell_open_iface(iCell)) localIface = localIface + 1
      end do
      call mpas_dmpar_sum_int(domain % dminfo, localIface, globalOpen)
      call mpas_dmpar_sum_int(domain % dminfo, nCellsSolve, globalOwned)

      hts_max_level_used = 0
      do L = 0, hts_level_cap
         if (globalCount(L) > 0) hts_max_level_used = L
      end do

      call mpas_log_write('')
      call mpas_log_write('--- hierarchical time stepping ---')
      if (hts_enabled) then
         call mpas_log_write(' config_hts = true, config_hts_max_level = $i', intArgs=(/config_hts_max_level/))
      else
         call mpas_log_write(' config_hts = false (integrator is unchanged)')
      end if
      call mpas_log_write(' owned cells = $i', intArgs=(/globalOwned/))
      do L = 0, hts_max_level_used
         call mpas_log_write(' level $i cells = $i', intArgs=(/L, globalCount(L)/))
      end do
      call mpas_log_write(' interface edges = $i', intArgs=(/globalIface/))
      call mpas_log_write(' open-interface cells = $i', intArgs=(/globalOpen/))
      if (uniform) then
         call mpas_log_write(' meshDensity does not vary; every cell is level 0')
      end if
      call mpas_log_write('----------------------------------')

      if (.not. hts_enabled) return

#ifdef MPAS_OPENACC
      call mpas_log_write('config_hts is not supported in an OpenACC build.', messageType=MPAS_LOG_CRIT)
#endif

      if (hts_max_level_used == 0) then
         call mpas_log_write('HTS enabled with a single level; masked integrator, no subcycling')
         return
      end if

      call mpas_log_write(' column physics runs on the active level; microphysics skips the other columns')

      call atm_hts_build_lists(block)
      call atm_hts_build_acoustic_halo(block)

      allocate(hts_ru_target(nVertLevels, nEdges+1))
      allocate(hts_flux_sum(nVertLevels, nEdges+1))
      hts_ru_target(:, :) = 0.0_RKIND
      hts_flux_sum(:, :) = 0.0_R8KIND

      allocate(snap_rho(nVertLevels, nCells+1, hts_max_level_used))
      allocate(delta_rho(nVertLevels, nCells+1, hts_max_level_used))
      allocate(snap_theta(nVertLevels, nCells+1, hts_max_level_used))
      allocate(delta_theta(nVertLevels, nCells+1, hts_max_level_used))
      allocate(snap_w(nVertLevels+1, nCells+1, hts_max_level_used))
      allocate(delta_w(nVertLevels+1, nCells+1, hts_max_level_used))
      allocate(snap_rw(nVertLevels+1, nCells+1, hts_max_level_used))
      allocate(delta_rw(nVertLevels+1, nCells+1, hts_max_level_used))
      allocate(snap_u(nVertLevels, nEdges+1, hts_max_level_used))
      allocate(delta_u(nVertLevels, nEdges+1, hts_max_level_used))
      allocate(hts_level_started(hts_max_level_used))
      snap_rho = 0.0_RKIND
      delta_rho = 0.0_RKIND
      snap_theta = 0.0_RKIND
      delta_theta = 0.0_RKIND
      snap_w = 0.0_RKIND
      delta_w = 0.0_RKIND
      snap_rw = 0.0_RKIND
      delta_rw = 0.0_RKIND
      snap_u = 0.0_RKIND
      delta_u = 0.0_RKIND
      hts_level_started(:) = .false.

      ! Silence unused-variable warnings on builds that do not walk neighbors here.
      iNbr = nEdgesOnCell(1)
      nbr = cellsOnCell(1, 1)
      if (iNbr < 0 .or. nbr < 0) call mpas_log_write('HTS neighbor tables are empty', messageType=MPAS_LOG_WARN)

   end subroutine atm_hts_setup


   subroutine atm_hts_finalize(domain)

      implicit none

      type (domain_type), intent(inout) :: domain

      ! domain is unused; arrays are process-local
      if (associated(domain % blocklist)) then
         continue
      end if

      if (allocated(hts_cell_level)) deallocate(hts_cell_level)
      if (allocated(hts_edge_level)) deallocate(hts_edge_level)
      if (allocated(hts_edge_interface)) deallocate(hts_edge_interface)
      if (allocated(hts_cell_open_iface)) deallocate(hts_cell_open_iface)
      if (allocated(hts_cell_live)) deallocate(hts_cell_live)
      if (allocated(hts_impose_edge)) deallocate(hts_impose_edge)
      if (allocated(hts_ru_target)) deallocate(hts_ru_target)
      if (allocated(hts_flux_sum)) deallocate(hts_flux_sum)
      if (allocated(snap_rho)) deallocate(snap_rho)
      if (allocated(delta_rho)) deallocate(delta_rho)
      if (allocated(snap_theta)) deallocate(snap_theta)
      if (allocated(delta_theta)) deallocate(delta_theta)
      if (allocated(snap_w)) deallocate(snap_w)
      if (allocated(delta_w)) deallocate(delta_w)
      if (allocated(snap_rw)) deallocate(snap_rw)
      if (allocated(delta_rw)) deallocate(delta_rw)
      if (allocated(snap_u)) deallocate(snap_u)
      if (allocated(delta_u)) deallocate(delta_u)
      if (allocated(hts_level_started)) deallocate(hts_level_started)
      if (allocated(hts_work_own)) deallocate(hts_work_own)
      if (allocated(hts_work_all)) deallocate(hts_work_all)
      if (allocated(hts_active_own)) deallocate(hts_active_own)
      if (allocated(hts_active_all)) deallocate(hts_active_all)
      if (allocated(hts_work_edge_own)) deallocate(hts_work_edge_own)
      if (allocated(hts_work_edge_all)) deallocate(hts_work_edge_all)
      if (allocated(hts_ac_edge)) deallocate(hts_ac_edge)
      if (allocated(hts_off_own)) deallocate(hts_off_own)
      if (allocated(hts_vert_all)) deallocate(hts_vert_all)
      if (allocated(hts_vert_own)) deallocate(hts_vert_own)
      if (allocated(hts_diag_all)) deallocate(hts_diag_all)
      if (allocated(hts_diag_edge_all)) deallocate(hts_diag_edge_all)
      if (allocated(hts_diag_vert_all)) deallocate(hts_diag_vert_all)
      if (allocated(hts_level_edge)) deallocate(hts_level_edge)
      if (allocated(hts_ac_cell)) deallocate(hts_ac_cell)
      if (allocated(hts_recon_own)) deallocate(hts_recon_own)
      if (allocated(hts_live_own)) deallocate(hts_live_own)
      if (allocated(hts_live_all)) deallocate(hts_live_all)
      if (allocated(hts_ah_send_proc)) deallocate(hts_ah_send_proc)
      if (allocated(hts_ah_send_n)) deallocate(hts_ah_send_n)
      if (allocated(hts_ah_send_off)) deallocate(hts_ah_send_off)
      if (allocated(hts_ah_send_cell)) deallocate(hts_ah_send_cell)
      if (allocated(hts_ah_recv_proc)) deallocate(hts_ah_recv_proc)
      if (allocated(hts_ah_recv_n)) deallocate(hts_ah_recv_n)
      if (allocated(hts_ah_recv_off)) deallocate(hts_ah_recv_off)
      if (allocated(hts_ah_recv_cell)) deallocate(hts_ah_recv_cell)
      if (allocated(hts_ah_sbuf)) deallocate(hts_ah_sbuf)
      if (allocated(hts_ah_rbuf)) deallocate(hts_ah_rbuf)
      hts_ah_ready = .false.

   end subroutine atm_hts_finalize


   subroutine atm_hts_build_lists(block)

      ! Owned and halo index lists for each level. Cell work lists add
      ! hts_work_rings neighbor rings so a masked tendency still sees its stencil.

      implicit none

      type (block_type), intent(in) :: block

      type (mpas_pool_type), pointer :: mesh
      integer, pointer :: nCells, nCellsSolve, nEdges, nEdgesSolve
      integer, pointer :: nVertices, nVerticesSolve, vertexDegree
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, verticesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell
      integer :: iCell, iEdge, i, iPass, iVertex, L, nbr, cell1, cell2, v1, v2, coarse
      integer :: nOwn, nAll, nEOwn, nEAll, nAc, nOff, nVAll, nVOwn, nMiss, nExtra
      logical :: edgeUpdated
      logical, allocatable :: mark(:), mark2(:), edgeMark(:), vertMark(:)

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_dimension(mesh, 'nVerticesSolve', nVerticesSolve)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)

      hts_nVertices = nVertices
      hts_nVerticesSolve = nVerticesSolve

      if (allocated(hts_work_own)) deallocate(hts_work_own)
      if (allocated(hts_work_all)) deallocate(hts_work_all)
      if (allocated(hts_active_own)) deallocate(hts_active_own)
      if (allocated(hts_active_all)) deallocate(hts_active_all)
      if (allocated(hts_work_edge_own)) deallocate(hts_work_edge_own)
      if (allocated(hts_work_edge_all)) deallocate(hts_work_edge_all)
      if (allocated(hts_ac_edge)) deallocate(hts_ac_edge)
      if (allocated(hts_off_own)) deallocate(hts_off_own)
      if (allocated(hts_vert_all)) deallocate(hts_vert_all)
      if (allocated(hts_vert_own)) deallocate(hts_vert_own)
      if (allocated(hts_diag_all)) deallocate(hts_diag_all)
      if (allocated(hts_diag_edge_all)) deallocate(hts_diag_edge_all)
      if (allocated(hts_diag_vert_all)) deallocate(hts_diag_vert_all)
      if (allocated(hts_level_edge)) deallocate(hts_level_edge)
      if (allocated(hts_ac_cell)) deallocate(hts_ac_cell)
      if (allocated(hts_recon_own)) deallocate(hts_recon_own)
      if (allocated(hts_live_own)) deallocate(hts_live_own)
      if (allocated(hts_live_all)) deallocate(hts_live_all)

      allocate(hts_work_own(nCells, 0:hts_level_cap))
      allocate(hts_work_all(nCells, 0:hts_level_cap))
      allocate(hts_active_own(nCells, 0:hts_level_cap))
      allocate(hts_active_all(nCells, 0:hts_level_cap))
      allocate(hts_work_edge_own(nEdges, 0:hts_level_cap))
      allocate(hts_work_edge_all(nEdges, 0:hts_level_cap))
      allocate(hts_ac_edge(nEdges, 0:hts_level_cap))
      allocate(hts_off_own(nCells, 0:hts_level_cap))
      allocate(hts_vert_all(nVertices, 0:hts_level_cap))
      allocate(hts_vert_own(nVertices, 0:hts_level_cap))
      allocate(hts_diag_all(nCells, 0:hts_level_cap))
      allocate(hts_diag_edge_all(nEdges, 0:hts_level_cap))
      allocate(hts_diag_vert_all(nVertices, 0:hts_level_cap))
      allocate(hts_level_edge(nEdges, 0:hts_level_cap))
      allocate(hts_ac_cell(nCells, 0:hts_level_cap))
      allocate(hts_recon_own(nCells, 0:hts_level_cap))
      allocate(hts_live_own(nCells, 0:hts_level_cap))
      allocate(hts_live_all(nCells, 0:hts_level_cap))
      hts_n_work_own = 0
      hts_n_work_all = 0
      hts_n_active_own = 0
      hts_n_active_all = 0
      hts_n_work_edge_own = 0
      hts_n_work_edge_all = 0
      hts_n_ac_edge = 0
      hts_n_off_own = 0
      hts_n_vert_all = 0
      hts_n_vert_own = 0
      hts_n_diag_all = 0
      hts_n_diag_edge_all = 0
      hts_n_diag_vert_all = 0
      hts_n_level_edge = 0
      hts_n_ac_cell = 0
      hts_n_recon_own = 0
      hts_n_live_own = 0
      hts_n_live_all = 0

      allocate(mark(nCells), mark2(nCells), edgeMark(nEdges), vertMark(nVertices))

      call mpas_log_write(' HTS work rings = $i, diagnostic rings = $i', intArgs=(/hts_work_rings, hts_diag_rings/))

      do L = 0, hts_max_level_used
         mark(:) = .false.
         do iCell = 1, nCells
            if (hts_cell_level(iCell) == L) mark(iCell) = .true.
         end do
         if (hts_work_rings <= 0) then
            do iCell = 1, nCells
               mark(iCell) = .true.
            end do
         else
            do iPass = 1, hts_work_rings
               mark2(:) = mark(:)
               do iCell = 1, nCells
                  if (.not. mark2(iCell)) cycle
                  do i = 1, nEdgesOnCell(iCell)
                     nbr = cellsOnCell(i, iCell)
                     if (nbr >= 1 .and. nbr <= nCells) mark(nbr) = .true.
                  end do
               end do
            end do
         end if

         nOwn = 0
         nAll = 0
         nOff = 0
         do iCell = 1, nCells
            if (hts_cell_level(iCell) == L) then
               nAll = nAll + 1
               hts_active_all(nAll, L) = iCell
               if (iCell <= nCellsSolve) then
                  nOwn = nOwn + 1
                  hts_active_own(nOwn, L) = iCell
               end if
            else if (iCell <= nCellsSolve) then
               nOff = nOff + 1
               hts_off_own(nOff, L) = iCell
            end if
         end do
         hts_n_active_all(L) = nAll
         hts_n_active_own(L) = nOwn
         hts_n_off_own(L) = nOff

         ! Acoustic columns are this level only. A coarser neighbor that
         ! shares an edge stays on the snapshot, including rho_pp and
         ! rtheta_pp. Job 199953 updated those perturbations alone and
         ! NaNed u at 00:18.
         hts_n_live_all(L) = hts_n_active_all(L)
         hts_n_live_own(L) = hts_n_active_own(L)
         do i = 1, hts_n_active_all(L)
            hts_live_all(i, L) = hts_active_all(i, L)
         end do
         do i = 1, hts_n_active_own(L)
            hts_live_own(i, L) = hts_active_own(i, L)
         end do

         nOwn = 0
         nAll = 0
         do iCell = 1, nCells
            if (.not. mark(iCell)) cycle
            nAll = nAll + 1
            hts_work_all(nAll, L) = iCell
            if (iCell <= nCellsSolve) then
               nOwn = nOwn + 1
               hts_work_own(nOwn, L) = iCell
            end if
         end do
         hts_n_work_all(L) = nAll
         hts_n_work_own(L) = nOwn

         edgeMark(:) = .false.
         nEOwn = 0
         nEAll = 0
         nAc = 0
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (cell1 < 1 .or. cell1 > nCells) cell1 = 0
            if (cell2 < 1 .or. cell2 > nCells) cell2 = 0
            if ((cell1 > 0 .and. mark(cell1)) .or. (cell2 > 0 .and. mark(cell2))) then
               edgeMark(iEdge) = .true.
               nEAll = nEAll + 1
               hts_work_edge_all(nEAll, L) = iEdge
               if (iEdge <= nEdgesSolve) then
                  nEOwn = nEOwn + 1
                  hts_work_edge_own(nEOwn, L) = iEdge
               end if
            end if
            ! Same edges the acoustic loop updates: any edge touching an owned cell.
            ! A missing neighbor (the garbage cell) still counts.
            if (.not. ((cell1 >= 1 .and. cell1 <= nCellsSolve) .or. &
                       (cell2 >= 1 .and. cell2 <= nCellsSolve))) cycle
            if (hts_edge_level(iEdge) == L) then
               nAc = nAc + 1
               hts_ac_edge(nAc, L) = iEdge
            else if (hts_edge_interface(iEdge) .and. cell1 >= 1 .and. cell2 >= 1) then
               coarse = max(hts_cell_level(cell1), hts_cell_level(cell2))
               if (coarse == L) then
                  nAc = nAc + 1
                  hts_ac_edge(nAc, L) = iEdge
               end if
            end if
         end do
         hts_n_work_edge_all(L) = nEAll
         hts_n_work_edge_own(L) = nEOwn
         hts_n_ac_edge(L) = nAc

         vertMark(:) = .false.
         do iEdge = 1, nEdges
            if (.not. edgeMark(iEdge)) cycle
            v1 = verticesOnEdge(1, iEdge)
            v2 = verticesOnEdge(2, iEdge)
            if (v1 >= 1 .and. v1 <= nVertices) vertMark(v1) = .true.
            if (v2 >= 1 .and. v2 <= nVertices) vertMark(v2) = .true.
         end do
         nVAll = 0
         nVOwn = 0
         do i = 1, nVertices
            if (.not. vertMark(i)) cycle
            nVAll = nVAll + 1
            hts_vert_all(nVAll, L) = i
            if (i <= nVerticesSolve) then
               nVOwn = nVOwn + 1
               hts_vert_own(nVOwn, L) = i
            end if
         end do
         hts_n_vert_all(L) = nVAll
         hts_n_vert_own(L) = nVOwn

         ! Four rings cover the APVM and kinetic-energy reads of the two-ring
         ! tendency, including the edges a vertex kinetic energy needs. The
         ! diagnostic write uses this set. There is no extra pass before RK1.
         do iPass = hts_work_rings + 1, hts_diag_rings
            mark2(:) = mark(:)
            do iCell = 1, nCells
               if (.not. mark2(iCell)) cycle
               do i = 1, nEdgesOnCell(iCell)
                  nbr = cellsOnCell(i, iCell)
                  if (nbr >= 1 .and. nbr <= nCells) mark(nbr) = .true.
               end do
            end do
         end do

         nAll = 0
         do iCell = 1, nCells
            if (.not. mark(iCell)) cycle
            nAll = nAll + 1
            hts_diag_all(nAll, L) = iCell
         end do
         hts_n_diag_all(L) = nAll

         edgeMark(:) = .false.
         nEAll = 0
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (cell1 < 1 .or. cell1 > nCells) cell1 = 0
            if (cell2 < 1 .or. cell2 > nCells) cell2 = 0
            if ((cell1 > 0 .and. mark(cell1)) .or. (cell2 > 0 .and. mark(cell2))) then
               edgeMark(iEdge) = .true.
               nEAll = nEAll + 1
               hts_diag_edge_all(nEAll, L) = iEdge
            end if
         end do

         vertMark(:) = .false.
         do iEdge = 1, nEdges
            if (.not. edgeMark(iEdge)) cycle
            v1 = verticesOnEdge(1, iEdge)
            v2 = verticesOnEdge(2, iEdge)
            if (v1 >= 1 .and. v1 <= nVertices) vertMark(v1) = .true.
            if (v2 >= 1 .and. v2 <= nVertices) vertMark(v2) = .true.
         end do
         do iVertex = 1, nVertices
            if (.not. vertMark(iVertex)) cycle
            do i = 1, vertexDegree
               iEdge = edgesOnVertex(i, iVertex)
               if (iEdge < 1 .or. iEdge > nEdges) cycle
               if (edgeMark(iEdge)) cycle
               edgeMark(iEdge) = .true.
               nEAll = nEAll + 1
               hts_diag_edge_all(nEAll, L) = iEdge
            end do
         end do
         hts_n_diag_edge_all(L) = nEAll

         nVAll = 0
         do i = 1, nVertices
            if (.not. vertMark(i)) cycle
            nVAll = nVAll + 1
            hts_diag_vert_all(nVAll, L) = i
         end do
         hts_n_diag_vert_all(L) = nVAll

         call mpas_log_write(' level $i index lists: active owned $i, work owned $i, edges $i', &
                             intArgs=(/L, hts_n_active_own(L), hts_n_work_own(L), hts_n_work_edge_own(L)/))
         call mpas_log_write(' level $i diagnostic list: cells $i, edges $i, vertices $i', &
                             intArgs=(/L, hts_n_diag_all(L), hts_n_diag_edge_all(L), hts_n_diag_vert_all(L)/))

         ! The acoustic list must be the same edges the unmasked loop would update.
         ! A raw cell id of 0 passes the original "<= nCellsSolve" test.
         nMiss = 0
         nExtra = 0
         do iEdge = 1, nEdges
            edgeMark(iEdge) = .false.
         end do
         do i = 1, nAc
            edgeMark(hts_ac_edge(i, L)) = .true.
         end do
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
               if (hts_edge_level(iEdge) == L .and. .not. edgeMark(iEdge)) nMiss = nMiss + 1
            end if
            if (edgeMark(iEdge) .and. hts_edge_level(iEdge) == L) then
               if (.not. (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve)) nExtra = nExtra + 1
            end if
         end do
         if (nMiss > 0 .or. nExtra > 0) then
            call mpas_log_write(' level $i acoustic list mismatch: missed $i extra $i', &
                                intArgs=(/L, nMiss, nExtra/), messageType=MPAS_LOG_WARN)
         end if

         nEAll = 0
         do iEdge = 1, nEdges
            if (hts_edge_level(iEdge) /= L) cycle
            nEAll = nEAll + 1
            hts_level_edge(nEAll, L) = iEdge
         end do
         hts_n_level_edge(L) = nEAll

         ! Cells an acoustic edge reads. The per-substep rtheta_pp_old copy and
         ! the inactive perturbation clear use this set instead of every cell.
         mark(:) = .false.
         nAll = 0
         do i = 1, nAc
            iEdge = hts_ac_edge(i, L)
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (cell1 >= 1 .and. cell1 <= nCells) then
               if (.not. mark(cell1)) then
                  mark(cell1) = .true.
                  nAll = nAll + 1
                  hts_ac_cell(nAll, L) = cell1
               end if
            end if
            if (cell2 >= 1 .and. cell2 <= nCells) then
               if (.not. mark(cell2)) then
                  mark(cell2) = .true.
                  nAll = nAll + 1
                  hts_ac_cell(nAll, L) = cell2
               end if
            end if
         end do
         hts_n_ac_cell(L) = nAll

         ! Owned cells whose reconstructed wind changes when this level updates
         ! an edge. Halo edges of this level are included: the exchange has
         ! already filled them before the reconstruction.
         mark(:) = .false.
         nOwn = 0
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (.not. ((cell1 >= 1 .and. cell1 <= nCellsSolve) .or. &
                       (cell2 >= 1 .and. cell2 <= nCellsSolve))) cycle
            edgeUpdated = hts_edge_level(iEdge) == L
            if (hts_edge_interface(iEdge) .and. cell1 >= 1 .and. cell1 <= nCells .and. &
                cell2 >= 1 .and. cell2 <= nCells) then
               if (max(hts_cell_level(cell1), hts_cell_level(cell2)) == L .and. &
                   hts_edge_level(iEdge) < L) edgeUpdated = .true.
            end if
            if (.not. edgeUpdated) cycle
            if (cell1 >= 1 .and. cell1 <= nCellsSolve) then
               if (.not. mark(cell1)) then
                  mark(cell1) = .true.
                  nOwn = nOwn + 1
                  hts_recon_own(nOwn, L) = cell1
               end if
            end if
            if (cell2 >= 1 .and. cell2 <= nCellsSolve) then
               if (.not. mark(cell2)) then
                  mark(cell2) = .true.
                  nOwn = nOwn + 1
                  hts_recon_own(nOwn, L) = cell2
               end if
            end if
         end do
         hts_n_recon_own(L) = nOwn

         call mpas_log_write(' level $i hold lists: edges $i, acoustic cells $i, recon owned $i', &
                             intArgs=(/L, hts_n_level_edge(L), hts_n_ac_cell(L), hts_n_recon_own(L)/))
         call mpas_log_write(' level $i acoustic columns: owned $i, parents excluded', &
                             intArgs=(/L, hts_n_live_own(L)/))
      end do

      deallocate(mark, mark2, edgeMark, vertMark)

   end subroutine atm_hts_build_lists


   subroutine atm_hts_mark_acoustic_cells()

      ! Mark the active level. Open-interface parents are left unmarked.
      ! Job 199953 included them in the acoustic pressure update and NaNed
      ! u at 00:18 on the eds land VR. The acoustic step uses the cell
      ! level directly and does not read this mask.

      implicit none

      integer :: i, iCell

      if (.not. allocated(hts_cell_live)) return
      hts_cell_live(:) = .false.
      if (hts_active_level < 0 .or. hts_active_level > hts_level_cap) return
      if (.not. allocated(hts_live_all)) return
      do i = 1, hts_n_live_all(hts_active_level)
         iCell = hts_live_all(i, hts_active_level)
         hts_cell_live(iCell) = .true.
      end do

   end subroutine atm_hts_mark_acoustic_cells


   logical function hts_cell_acoustic(iCell)

      ! True only for the active level. A coarser open-interface parent is
      ! not advanced. Eds job 199953 updated that parent's rho_pp and
      ! rtheta_pp while u and w stayed on the snapshot: at 00:17:30 u was
      ! about ±7.2e5 m/s and w about -1e4 / +8e4 m/s, and u was NaN at
      ! 00:18. Ramping held-ring rho and theta (cf1dc652d, job 199977)
      ! NaNed u at 00:12:30. The pure hold on d3fe6fb75 stayed finite until
      ! about 01:57.

      implicit none

      integer, intent(in) :: iCell

      hts_cell_acoustic = hts_cell_level(iCell) == hts_active_level

   end function hts_cell_acoustic


   subroutine hts_fetch_fields(block, rho1, rho2, th1, th2, w1, w2, u1, u2, &
                               rho_p, rho_p_save, rtheta_p, rtheta_p_save, &
                               pressure_p, exner, exner_base, rtheta_base, rho_base, zz, &
                               ru, ru_save, rw, rw_save, rho_edge)

      implicit none

      type (block_type), intent(in) :: block
      real (kind=RKIND), dimension(:,:), pointer :: rho1, rho2, th1, th2, w1, w2, u1, u2
      real (kind=RKIND), dimension(:,:), pointer :: rho_p, rho_p_save, rtheta_p, rtheta_p_save
      real (kind=RKIND), dimension(:,:), pointer :: pressure_p, exner, exner_base, rtheta_base, rho_base, zz
      real (kind=RKIND), dimension(:,:), pointer :: ru, ru_save, rw, rw_save, rho_edge

      type (mpas_pool_type), pointer :: state, diag, mesh

      nullify(state, diag, mesh)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)

      call mpas_pool_get_array(state, 'rho_zz', rho1, 1)
      call mpas_pool_get_array(state, 'rho_zz', rho2, 2)
      call mpas_pool_get_array(state, 'theta_m', th1, 1)
      call mpas_pool_get_array(state, 'theta_m', th2, 2)
      call mpas_pool_get_array(state, 'w', w1, 1)
      call mpas_pool_get_array(state, 'w', w2, 2)
      call mpas_pool_get_array(state, 'u', u1, 1)
      call mpas_pool_get_array(state, 'u', u2, 2)

      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rho_p_save', rho_p_save)
      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(diag, 'rtheta_p_save', rtheta_p_save)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'exner_base', exner_base)
      call mpas_pool_get_array(diag, 'rtheta_base', rtheta_base)
      call mpas_pool_get_array(diag, 'rho_base', rho_base)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'ru_save', ru_save)
      call mpas_pool_get_array(diag, 'rw', rw)
      call mpas_pool_get_array(diag, 'rw_save', rw_save)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)
      call mpas_pool_get_array(mesh, 'zz', zz)

   end subroutine hts_fetch_fields


   subroutine hts_diagnostics_from_state(k, iCell, rho_zz, theta_m, rho_base, rtheta_base, zz, &
                                         exner_base, rho_p, rtheta_p, exner, pressure_p)

      implicit none

      integer, intent(in) :: k, iCell
      real (kind=RKIND), intent(in) :: rho_zz, theta_m, rho_base, rtheta_base, zz, exner_base
      real (kind=RKIND), intent(out) :: rho_p, rtheta_p, exner, pressure_p
      real (kind=RKIND) :: rcv, p0, rtheta_tot

      rcv = rgas / (cp - rgas)
      p0 = 1.0e5_RKIND
      rho_p = rho_zz - rho_base
      rtheta_tot = theta_m * rho_zz
      rtheta_p = rtheta_tot - rtheta_base
      exner = (zz * (rgas / p0) * rtheta_tot) ** rcv
      pressure_p = zz * rgas * (exner * rtheta_p + rtheta_base * (exner - exner_base))

      ! iCell is part of the signature so callers can pass the column index
      ! without a dummy local. The column is horizontally local.
      if (iCell < 0) pressure_p = pressure_p

   end subroutine hts_diagnostics_from_state


   subroutine atm_hts_begin_intervals(block)

      implicit none

      type (block_type), intent(in) :: block

      real (kind=RKIND), dimension(:,:), pointer :: rho1, rho2, th1, th2, w1, w2, u1, u2
      real (kind=RKIND), dimension(:,:), pointer :: rho_p, rho_p_save, rtheta_p, rtheta_p_save
      real (kind=RKIND), dimension(:,:), pointer :: pressure_p, exner, exner_base, rtheta_base, rho_base, zz
      real (kind=RKIND), dimension(:,:), pointer :: ru, ru_save, rw, rw_save, rho_edge
      integer :: L, i, iCell, iEdge, k, ratio

      if (hts_max_level_used < 1) return

      call mpas_timer_start('hts_begin_intervals')

      call hts_fetch_fields(block, rho1, rho2, th1, th2, w1, w2, u1, u2, &
                            rho_p, rho_p_save, rtheta_p, rtheta_p_save, &
                            pressure_p, exner, exner_base, rtheta_base, rho_base, zz, &
                            ru, ru_save, rw, rw_save, rho_edge)

      do L = 1, hts_max_level_used
         ratio = ishft(1, L)
         if (mod(hts_tick, ratio) /= 0) cycle

         if (.not. hts_level_started(L)) then
            do i = 1, hts_n_active_all(L)
               iCell = hts_active_all(i, L)
               do k = 1, hts_nVert
                  snap_rho(k, iCell, L) = rho1(k, iCell)
                  snap_theta(k, iCell, L) = th1(k, iCell)
                  delta_rho(k, iCell, L) = 0.0_RKIND
                  delta_theta(k, iCell, L) = 0.0_RKIND
               end do
               do k = 1, hts_nVert + 1
                  snap_w(k, iCell, L) = w1(k, iCell)
                  snap_rw(k, iCell, L) = rw(k, iCell)
                  delta_w(k, iCell, L) = 0.0_RKIND
                  delta_rw(k, iCell, L) = 0.0_RKIND
               end do
            end do
            do i = 1, hts_n_level_edge(L)
               iEdge = hts_level_edge(i, L)
               do k = 1, hts_nVert
                  snap_u(k, iEdge, L) = u1(k, iEdge)
                  delta_u(k, iEdge, L) = 0.0_RKIND
               end do
            end do
            hts_level_started(L) = .true.
         else
            do i = 1, hts_n_active_all(L)
               iCell = hts_active_all(i, L)
               do k = 1, hts_nVert
                  delta_rho(k, iCell, L) = rho1(k, iCell) - snap_rho(k, iCell, L)
                  delta_theta(k, iCell, L) = th1(k, iCell) - snap_theta(k, iCell, L)
                  snap_rho(k, iCell, L) = rho1(k, iCell)
                  snap_theta(k, iCell, L) = th1(k, iCell)
               end do
               do k = 1, hts_nVert + 1
                  delta_w(k, iCell, L) = w1(k, iCell) - snap_w(k, iCell, L)
                  delta_rw(k, iCell, L) = rw(k, iCell) - snap_rw(k, iCell, L)
                  snap_w(k, iCell, L) = w1(k, iCell)
                  snap_rw(k, iCell, L) = rw(k, iCell)
               end do
            end do
            do i = 1, hts_n_level_edge(L)
               iEdge = hts_level_edge(i, L)
               do k = 1, hts_nVert
                  delta_u(k, iEdge, L) = u1(k, iEdge) - snap_u(k, iEdge, L)
                  snap_u(k, iEdge, L) = u1(k, iEdge)
               end do
            end do
         end if
      end do

      call mpas_timer_stop('hts_begin_intervals')

   end subroutine atm_hts_begin_intervals


   subroutine atm_hts_apply_interpolation(block, activeLevel)

      implicit none

      type (block_type), intent(in) :: block
      integer, intent(in) :: activeLevel

      real (kind=RKIND), dimension(:,:), pointer :: rho1, rho2, th1, th2, w1, w2, u1, u2
      real (kind=RKIND), dimension(:,:), pointer :: rho_p, rho_p_save, rtheta_p, rtheta_p_save
      real (kind=RKIND), dimension(:,:), pointer :: pressure_p, exner, exner_base, rtheta_base, rho_base, zz
      real (kind=RKIND), dimension(:,:), pointer :: ru, ru_save, rw, rw_save, rho_edge
      integer, dimension(:,:), pointer :: cellsOnEdge
      type (mpas_pool_type), pointer :: mesh
      integer :: L, i, iCell, iEdge, k, cell1, cell2
      real (kind=RKIND) :: alpha, rho_here, th_here, u_here
      logical :: any

      if (hts_max_level_used <= activeLevel) return

      call mpas_timer_start('hts_apply_interpolation')

      ! A level whose state already matches its snapshot is left alone. The
      ! fine step does not write that level's prognostic fields, and the time
      ! levels stay equal across the shift until the level itself steps.
      any = .false.
      do L = activeLevel + 1, hts_max_level_used
         if (.not. hts_hold_valid(L)) any = .true.
      end do
      if (.not. any) then
         call mpas_timer_stop('hts_apply_interpolation')
         return
      end if

      call hts_fetch_fields(block, rho1, rho2, th1, th2, w1, w2, u1, u2, &
                            rho_p, rho_p_save, rtheta_p, rtheta_p_save, &
                            pressure_p, exner, exner_base, rtheta_base, rho_base, zz, &
                            ru, ru_save, rw, rw_save, rho_edge)
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      ! Hold coarser levels at the start-of-interval snapshot (alpha = 0)
      ! while a finer level steps. A diagnostic-ring rho/theta ramp of the
      ! previous increment (cf1dc652d) NaNed u at 00:12:30 on eds job
      ! 199977, earlier than this hold. Extrapolating u and w fed the coarse
      ! acoustic transient into the fine flux. Cells are written before
      ! edges so an edge density sees both columns.
      alpha = 0.0_RKIND
      do L = activeLevel + 1, hts_max_level_used
         if (hts_hold_valid(L)) cycle
         do i = 1, hts_n_active_all(L)
            iCell = hts_active_all(i, L)
            do k = 1, hts_nVert
               rho_here = snap_rho(k, iCell, L) + alpha * delta_rho(k, iCell, L)
               th_here = snap_theta(k, iCell, L) + alpha * delta_theta(k, iCell, L)
               rho1(k, iCell) = rho_here
               rho2(k, iCell) = rho_here
               th1(k, iCell) = th_here
               th2(k, iCell) = th_here
               call hts_diagnostics_from_state(k, iCell, rho_here, th_here, rho_base(k, iCell), &
                    rtheta_base(k, iCell), zz(k, iCell), exner_base(k, iCell), &
                    rho_p(k, iCell), rtheta_p(k, iCell), exner(k, iCell), pressure_p(k, iCell))
               rho_p_save(k, iCell) = rho_p(k, iCell)
               rtheta_p_save(k, iCell) = rtheta_p(k, iCell)
            end do
            do k = 1, hts_nVert + 1
               w1(k, iCell) = snap_w(k, iCell, L) + alpha * delta_w(k, iCell, L)
               w2(k, iCell) = w1(k, iCell)
               rw(k, iCell) = snap_rw(k, iCell, L) + alpha * delta_rw(k, iCell, L)
               rw_save(k, iCell) = rw(k, iCell)
            end do
         end do
      end do
      do L = activeLevel + 1, hts_max_level_used
         if (hts_hold_valid(L)) cycle
         do i = 1, hts_n_level_edge(L)
            iEdge = hts_level_edge(i, L)
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            do k = 1, hts_nVert
               u_here = snap_u(k, iEdge, L) + alpha * delta_u(k, iEdge, L)
               u1(k, iEdge) = u_here
               u2(k, iEdge) = u_here
               ru(k, iEdge) = 0.5_RKIND * (rho2(k, cell1) + rho2(k, cell2)) * u_here
               ru_save(k, iEdge) = ru(k, iEdge)
               rho_edge(k, iEdge) = 0.5_RKIND * (rho2(k, cell1) + rho2(k, cell2))
            end do
         end do
         hts_hold_valid(L) = .true.
         hts_hold_rewrote(L) = .true.
      end do

      call mpas_timer_stop('hts_apply_interpolation')

   end subroutine atm_hts_apply_interpolation


   subroutine atm_hts_restore_above(block, activeLevel)

      implicit none

      type (block_type), intent(in) :: block
      integer, intent(in) :: activeLevel

      real (kind=RKIND), dimension(:,:), pointer :: rho1, rho2, th1, th2, w1, w2, u1, u2
      real (kind=RKIND), dimension(:,:), pointer :: rho_p, rho_p_save, rtheta_p, rtheta_p_save
      real (kind=RKIND), dimension(:,:), pointer :: pressure_p, exner, exner_base, rtheta_base, rho_base, zz
      real (kind=RKIND), dimension(:,:), pointer :: ru, ru_save, rw, rw_save, rho_edge
      integer, dimension(:,:), pointer :: cellsOnEdge
      type (mpas_pool_type), pointer :: mesh
      integer :: L, i, iCell, iEdge, k, cell1, cell2
      logical :: any

      if (hts_max_level_used <= activeLevel) return

      call mpas_timer_start('hts_restore_above')

      ! The active step writes its own columns. Levels that still match their
      ! snapshot were not written, so putting the snapshot back is a no-op.
      any = .false.
      do L = activeLevel + 1, hts_max_level_used
         if (.not. hts_hold_valid(L)) any = .true.
      end do
      if (.not. any) then
         call mpas_timer_stop('hts_restore_above')
         return
      end if

      call hts_fetch_fields(block, rho1, rho2, th1, th2, w1, w2, u1, u2, &
                            rho_p, rho_p_save, rtheta_p, rtheta_p_save, &
                            pressure_p, exner, exner_base, rtheta_base, rho_base, zz, &
                            ru, ru_save, rw, rw_save, rho_edge)
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      do L = activeLevel + 1, hts_max_level_used
         if (hts_hold_valid(L)) cycle
         do i = 1, hts_n_active_all(L)
            iCell = hts_active_all(i, L)
            do k = 1, hts_nVert
               rho1(k, iCell) = snap_rho(k, iCell, L)
               rho2(k, iCell) = snap_rho(k, iCell, L)
               th1(k, iCell) = snap_theta(k, iCell, L)
               th2(k, iCell) = snap_theta(k, iCell, L)
               call hts_diagnostics_from_state(k, iCell, snap_rho(k, iCell, L), snap_theta(k, iCell, L), &
                    rho_base(k, iCell), rtheta_base(k, iCell), zz(k, iCell), exner_base(k, iCell), &
                    rho_p(k, iCell), rtheta_p(k, iCell), exner(k, iCell), pressure_p(k, iCell))
               rho_p_save(k, iCell) = rho_p(k, iCell)
               rtheta_p_save(k, iCell) = rtheta_p(k, iCell)
            end do
            do k = 1, hts_nVert + 1
               w1(k, iCell) = snap_w(k, iCell, L)
               w2(k, iCell) = snap_w(k, iCell, L)
               rw(k, iCell) = snap_rw(k, iCell, L)
               rw_save(k, iCell) = snap_rw(k, iCell, L)
            end do
         end do
      end do
      do L = activeLevel + 1, hts_max_level_used
         if (hts_hold_valid(L)) cycle
         do i = 1, hts_n_level_edge(L)
            iEdge = hts_level_edge(i, L)
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            do k = 1, hts_nVert
               u1(k, iEdge) = snap_u(k, iEdge, L)
               u2(k, iEdge) = snap_u(k, iEdge, L)
               ru(k, iEdge) = 0.5_RKIND * (rho2(k, cell1) + rho2(k, cell2)) * snap_u(k, iEdge, L)
               ru_save(k, iEdge) = ru(k, iEdge)
               rho_edge(k, iEdge) = 0.5_RKIND * (rho2(k, cell1) + rho2(k, cell2))
            end do
         end do
         hts_hold_valid(L) = .true.
         hts_hold_rewrote(L) = .true.
      end do

      call mpas_timer_stop('hts_restore_above')

   end subroutine atm_hts_restore_above


   subroutine atm_hts_sync_finer(block, activeLevel)

      ! Make levels that already advanced this tick look like a consistent
      ! current state to the coarser solve (time level 1 matches time level 2).

      implicit none

      type (block_type), intent(in) :: block
      integer, intent(in) :: activeLevel

      real (kind=RKIND), dimension(:,:), pointer :: rho1, rho2, th1, th2, w1, w2, u1, u2
      real (kind=RKIND), dimension(:,:), pointer :: rho_p, rho_p_save, rtheta_p, rtheta_p_save
      real (kind=RKIND), dimension(:,:), pointer :: pressure_p, exner, exner_base, rtheta_base, rho_base, zz
      real (kind=RKIND), dimension(:,:), pointer :: ru, ru_save, rw, rw_save, rho_edge
      integer :: iCell, iEdge, k

      call hts_fetch_fields(block, rho1, rho2, th1, th2, w1, w2, u1, u2, &
                            rho_p, rho_p_save, rtheta_p, rtheta_p_save, &
                            pressure_p, exner, exner_base, rtheta_base, rho_base, zz, &
                            ru, ru_save, rw, rw_save, rho_edge)

      do iCell = 1, hts_nCells
         if (hts_cell_level(iCell) >= activeLevel) cycle
         do k = 1, hts_nVert
            rho1(k, iCell) = rho2(k, iCell)
            th1(k, iCell) = th2(k, iCell)
            rho_p_save(k, iCell) = rho_p(k, iCell)
            rtheta_p_save(k, iCell) = rtheta_p(k, iCell)
         end do
         do k = 1, hts_nVert + 1
            w1(k, iCell) = w2(k, iCell)
            rw_save(k, iCell) = rw(k, iCell)
         end do
      end do
      do iEdge = 1, hts_nEdges
         if (hts_edge_level(iEdge) >= activeLevel) cycle
         do k = 1, hts_nVert
            u1(k, iEdge) = u2(k, iEdge)
            ru_save(k, iEdge) = ru(k, iEdge)
         end do
      end do

   end subroutine atm_hts_sync_finer


   subroutine atm_hts_accumulate_flux(block, level, dt)

      implicit none

      type (block_type), intent(in) :: block
      integer, intent(in) :: level
      real (kind=RKIND), intent(in) :: dt

      type (mpas_pool_type), pointer :: diag
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg
      integer :: iEdge, k

      if (.not. allocated(hts_flux_sum)) return

      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)

      do iEdge = 1, hts_nEdges
         if (hts_edge_level(iEdge) /= level) cycle
         if (.not. hts_edge_interface(iEdge)) cycle
         do k = 1, hts_nVert
            hts_flux_sum(k, iEdge) = hts_flux_sum(k, iEdge) + real(ruAvg(k, iEdge), R8KIND) * real(dt, R8KIND)
         end do
      end do

   end subroutine atm_hts_accumulate_flux


   subroutine atm_hts_prepare_targets(block, level, dt_level)

      implicit none

      type (block_type), intent(in) :: block
      integer, intent(in) :: level
      real (kind=RKIND), intent(in) :: dt_level

      integer, dimension(:,:), pointer :: cellsOnEdge
      type (mpas_pool_type), pointer :: mesh
      integer :: iEdge, k, cell1, cell2, coarse

      hts_impose_edge(:) = .false.
      if (.not. allocated(hts_ru_target)) return
      hts_ru_target(:, :) = 0.0_RKIND

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)

      do iEdge = 1, hts_nEdges
         if (.not. hts_edge_interface(iEdge)) cycle
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         if (cell1 < 1 .or. cell1 > hts_nCells .or. cell2 < 1 .or. cell2 > hts_nCells) cycle
         coarse = max(hts_cell_level(cell1), hts_cell_level(cell2))
         if (coarse /= level) cycle
         if (hts_edge_level(iEdge) >= level) cycle
         hts_impose_edge(iEdge) = .true.
         do k = 1, hts_nVert
            hts_ru_target(k, iEdge) = real(hts_flux_sum(k, iEdge) / real(dt_level, R8KIND), RKIND)
         end do
      end do

      ! block is used through mesh; silence unused if a compiler cannot see it
      if (.false.) call mpas_log_write('prepare targets on level $i', intArgs=(/level/))

   end subroutine atm_hts_prepare_targets


   subroutine atm_hts_clear_consumed(level)

      implicit none

      integer, intent(in) :: level

      integer :: iEdge, k

      if (.not. allocated(hts_flux_sum)) return

      ! Impose flags mark exactly the edges this level consumes.
      ! `level` selects the caller; the flags were built for that level.
      if (level < 0) return
      do iEdge = 1, hts_nEdges
         if (.not. hts_impose_edge(iEdge)) cycle
         do k = 1, hts_nVert
            hts_flux_sum(k, iEdge) = 0.0_R8KIND
         end do
         hts_impose_edge(iEdge) = .false.
      end do

   end subroutine atm_hts_clear_consumed


   subroutine atm_hts_correct_mass(block, level, dt_level)

      ! If the acoustic impose did not land on the accumulated fine flux,
      ! patch the coarse column so the committed horizontal flux matches.

      implicit none

      type (block_type), intent(in) :: block
      integer, intent(in) :: level
      real (kind=RKIND), intent(in) :: dt_level

      type (mpas_pool_type), pointer :: mesh, diag, state
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, pointer :: nCellsSolve
      real (kind=RKIND), dimension(:,:), pointer :: edgesOnCell_sign
      real (kind=RKIND), dimension(:), pointer :: dvEdge, invAreaCell
      real (kind=RKIND), dimension(:,:), pointer :: ruAvg, ru, rho_p, rtheta_p, pressure_p, exner, exner_base
      real (kind=RKIND), dimension(:,:), pointer :: rtheta_base, rho_base, zz, rho2, th2, u2
      integer :: iCell, iEdge, i, k, cell1, cell2
      real (kind=RKIND) :: drho, theta_edge, resid, residMax
      real (kind=RKIND) :: rcv, p0, rtheta_tot
      real (kind=RKIND) :: residOwned

      if (.not. allocated(hts_ru_target)) return

      residMax = 0.0_RKIND
      residOwned = 0.0_RKIND
      rcv = rgas / (cp - rgas)
      p0 = 1.0e5_RKIND

      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'state', state)
      call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'dvEdge', dvEdge)
      call mpas_pool_get_array(mesh, 'invAreaCell', invAreaCell)
      call mpas_pool_get_array(mesh, 'zz', zz)
      call mpas_pool_get_array(diag, 'ruAvg', ruAvg)
      call mpas_pool_get_array(diag, 'ru', ru)
      call mpas_pool_get_array(diag, 'rho_p', rho_p)
      call mpas_pool_get_array(diag, 'rtheta_p', rtheta_p)
      call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
      call mpas_pool_get_array(diag, 'exner', exner)
      call mpas_pool_get_array(diag, 'exner_base', exner_base)
      call mpas_pool_get_array(diag, 'rtheta_base', rtheta_base)
      call mpas_pool_get_array(diag, 'rho_base', rho_base)
      call mpas_pool_get_array(state, 'rho_zz', rho2, 2)
      call mpas_pool_get_array(state, 'theta_m', th2, 2)
      call mpas_pool_get_array(state, 'u', u2, 2)

      ! Re-mark impose edges; clear_consumed has not run yet, but prepare
      ! may have been followed by the solve which does not clear the flags.
      do iEdge = 1, hts_nEdges
         if (.not. hts_impose_edge(iEdge)) cycle
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         do k = 1, hts_nVert
            resid = abs(ruAvg(k, iEdge) - hts_ru_target(k, iEdge))
            if (resid > residMax) residMax = resid
            if ((cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) .and. resid > residOwned) then
               residOwned = resid
            end if
         end do
      end do
      if (residMax > hts_max_flux_resid) hts_max_flux_resid = residMax

      do iCell = 1, nCellsSolve
         if (hts_cell_level(iCell) /= level) cycle
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (.not. hts_impose_edge(iEdge)) cycle
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            do k = 1, hts_nVert
               drho = -edgesOnCell_sign(i, iCell) * dvEdge(iEdge) * invAreaCell(iCell) * dt_level &
                      * (hts_ru_target(k, iEdge) - ruAvg(k, iEdge))
               theta_edge = 0.5_RKIND * (th2(k, cell1) + th2(k, cell2))
               rho2(k, iCell) = rho2(k, iCell) + drho
               rho_p(k, iCell) = rho_p(k, iCell) + drho
               rtheta_p(k, iCell) = rtheta_p(k, iCell) + theta_edge * drho
               rtheta_tot = rtheta_p(k, iCell) + rtheta_base(k, iCell)
               if (rho2(k, iCell) /= 0.0_RKIND) then
                  th2(k, iCell) = rtheta_tot / rho2(k, iCell)
               end if
               exner(k, iCell) = (zz(k, iCell) * (rgas / p0) * rtheta_tot) ** rcv
               pressure_p(k, iCell) = zz(k, iCell) * rgas * (exner(k, iCell) * rtheta_p(k, iCell) &
                    + rtheta_base(k, iCell) * (exner(k, iCell) - exner_base(k, iCell)))
               ruAvg(k, iEdge) = hts_ru_target(k, iEdge)
               ru(k, iEdge) = hts_ru_target(k, iEdge)
               u2(k, iEdge) = 2.0_RKIND * ru(k, iEdge) / (rho2(k, cell1) + rho2(k, cell2))
            end do
         end do
      end do

      hts_correct_logs = hts_correct_logs + 1
      if (hts_correct_logs <= 3 .or. residOwned > 1.0e-3_RKIND) then
         call mpas_log_write('HTS level $i owned flux residual max $r (halo residual $r)', &
                             intArgs=(/level/), realArgs=(/residOwned, residMax/))
      end if

   end subroutine atm_hts_correct_mass


   subroutine hts_restrict_physics_tendencies(tend_ru, tend_rtheta, tend_rho)

      ! Column physics runs on the active level only and writes per-second
      ! rates there. physics_get_tend copies those rates, including rates
      ! left from the previous step on columns this call did not touch, and
      ! the RK multiplies by this step's dt. Keep the copy only on the
      ! active level so a coarse column is not forced on a fine tick and a
      ! fine column is not forced again over the long step. Single-level
      ! and all-fine return without writing.

      implicit none

      real (kind=RKIND), intent(inout) :: tend_ru(:,:)
      real (kind=RKIND), intent(inout) :: tend_rtheta(:,:)
      real (kind=RKIND), intent(inout) :: tend_rho(:,:)

      integer :: iCell, iEdge, k, nVert, nCells, nEdges

      if (.not. hts_masking) return
      if (hts_max_level_used <= 0) return
      if (.not. allocated(hts_cell_level)) return
      if (.not. allocated(hts_edge_level)) return

      nVert = size(tend_rho, 1)
      nCells = min(size(tend_rho, 2), size(hts_cell_level))
      nEdges = min(size(tend_ru, 2), size(hts_edge_level))

      do iCell = 1, nCells
         if (hts_cell_level(iCell) == hts_active_level) cycle
         do k = 1, nVert
            tend_rtheta(k, iCell) = 0.0_RKIND
            tend_rho(k, iCell) = 0.0_RKIND
         end do
      end do

      nVert = min(nVert, size(tend_ru, 1))
      do iEdge = 1, nEdges
         if (hts_edge_level(iEdge) == hts_active_level) cycle
         do k = 1, nVert
            tend_ru(k, iEdge) = 0.0_RKIND
         end do
      end do

   end subroutine hts_restrict_physics_tendencies


   subroutine atm_hts_build_acoustic_halo(block)

      ! Compact layer-1 lists for cells that touch the active level. A cell is
      ! in the list when its own level matches or one of its edges does, which
      ! is the same answer on every rank that holds the cell. Send and receive
      ! counts are checked with an all-to-all. A mismatch keeps the full
      ! dynamics:rho_pp and dynamics:rtheta_pp groups.

#ifdef MPAS_USE_MPI_F08
      use mpi_f08, only : MPI_Comm, MPI_INTEGER, MPI_SUM, MPI_Alltoall, MPI_Allreduce
#else
      use mpi, only : MPI_INTEGER, MPI_SUM, MPI_Alltoall, MPI_Allreduce
#endif

      implicit none

      type (block_type), intent(in) :: block

      type (mpas_pool_type), pointer :: diag, mesh
      type (field2DReal), pointer :: fld
      type (mpas_exchange_list), pointer :: exch
      integer, dimension(:,:), pointer :: cellsOnEdge
      logical, allocatable :: touch(:,:)
      integer, allocatable :: sendCounts(:), recvExpect(:), myRecv(:)
      integer :: L, iEdge, iCell, cell1, cell2, j, n, nmsg, sendOff, recvOff, fullN, sent
      integer :: nprocs, maxMsg, maxBuf, localBad, globalBad, ierr
#ifdef MPAS_USE_MPI_F08
      type (MPI_Comm) :: comm
#else
      integer :: comm
#endif

      hts_ah_ready = .false.
      hts_ah_nsend(:) = 0
      hts_ah_nrecv(:) = 0
      if (hts_max_level_used < 1) return
      if (.not. allocated(hts_cell_level)) return

      nullify(diag, mesh, fld)
      call mpas_pool_get_subpool(block % structs, 'diag', diag)
      call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
      call mpas_pool_get_field(diag, 'rho_pp', fld)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      if (.not. associated(fld)) return
      if (.not. associated(fld % sendList)) return
      if (.not. associated(fld % recvList)) return
      if (.not. associated(fld % sendList % halos)) return
      if (size(fld % sendList % halos) < 1) return
      if (size(fld % recvList % halos) < 1) return

      nprocs = block % domain % dminfo % nprocs
      comm = block % domain % dminfo % comm

      allocate(touch(hts_nCells, 0:hts_max_level_used))
      touch(:, :) = .false.
      do iCell = 1, hts_nCells
         L = hts_cell_level(iCell)
         if (L >= 0 .and. L <= hts_max_level_used) touch(iCell, L) = .true.
      end do
      do iEdge = 1, hts_nEdges
         L = hts_edge_level(iEdge)
         if (L < 0 .or. L > hts_max_level_used) cycle
         cell1 = cellsOnEdge(1, iEdge)
         cell2 = cellsOnEdge(2, iEdge)
         if (cell1 >= 1 .and. cell1 <= hts_nCells) touch(cell1, L) = .true.
         if (cell2 >= 1 .and. cell2 <= hts_nCells) touch(cell2, L) = .true.
      end do

      maxMsg = 0
      exch => fld % sendList % halos(1) % exchList
      do while (associated(exch))
         maxMsg = maxMsg + 1
         exch => exch % next
      end do
      exch => fld % recvList % halos(1) % exchList
      n = 0
      do while (associated(exch))
         n = n + 1
         exch => exch % next
      end do
      maxMsg = max(maxMsg, n, 1)

      allocate(hts_ah_send_proc(maxMsg, 0:hts_max_level_used))
      allocate(hts_ah_send_n(maxMsg, 0:hts_max_level_used))
      allocate(hts_ah_send_off(maxMsg, 0:hts_max_level_used))
      allocate(hts_ah_recv_proc(maxMsg, 0:hts_max_level_used))
      allocate(hts_ah_recv_n(maxMsg, 0:hts_max_level_used))
      allocate(hts_ah_recv_off(maxMsg, 0:hts_max_level_used))
      allocate(hts_ah_send_cell(max(hts_nCells, 1) * (hts_max_level_used + 1)))
      allocate(hts_ah_recv_cell(max(hts_nCells, 1) * (hts_max_level_used + 1)))
      hts_ah_send_n(:, :) = 0
      hts_ah_recv_n(:, :) = 0

      allocate(sendCounts(0:nprocs-1), recvExpect(0:nprocs-1), myRecv(0:nprocs-1))

      localBad = 0
      maxBuf = 0
      sendOff = 0
      recvOff = 0
      do L = 0, hts_max_level_used
         nmsg = 0
         fullN = 0
         sent = 0
         exch => fld % sendList % halos(1) % exchList
         do while (associated(exch))
            fullN = fullN + exch % nList
            n = 0
            do iCell = 1, exch % nList
               cell1 = exch % srcList(iCell)
               if (cell1 >= 1 .and. cell1 <= hts_nCells) then
                  if (touch(cell1, L)) then
                     if (sendOff + sent + n + 1 > size(hts_ah_send_cell)) then
                        localBad = 1
                     else
                        hts_ah_send_cell(sendOff + sent + n + 1) = cell1
                     end if
                     n = n + 1
                  end if
               end if
            end do
            if (n > 0) then
               nmsg = nmsg + 1
               if (nmsg > maxMsg) then
                  localBad = 1
               else
                  hts_ah_send_proc(nmsg, L) = exch % endPointID
                  hts_ah_send_n(nmsg, L) = n
                  hts_ah_send_off(nmsg, L) = sendOff + sent
               end if
               sent = sent + n
            end if
            exch => exch % next
         end do
         hts_ah_nsend(L) = nmsg

         nmsg = 0
         n = 0
         exch => fld % recvList % halos(1) % exchList
         do while (associated(exch))
            j = 0
            do iCell = 1, exch % nList
               cell1 = exch % destList(iCell)
               if (cell1 >= 1 .and. cell1 <= hts_nCells) then
                  if (touch(cell1, L)) then
                     if (recvOff + n + j + 1 > size(hts_ah_recv_cell)) then
                        localBad = 1
                     else
                        hts_ah_recv_cell(recvOff + n + j + 1) = cell1
                     end if
                     j = j + 1
                  end if
               end if
            end do
            if (j > 0) then
               nmsg = nmsg + 1
               if (nmsg > maxMsg) then
                  localBad = 1
               else
                  hts_ah_recv_proc(nmsg, L) = exch % endPointID
                  hts_ah_recv_n(nmsg, L) = j
                  hts_ah_recv_off(nmsg, L) = recvOff + n
               end if
               n = n + j
            end if
            exch => exch % next
         end do
         hts_ah_nrecv(L) = nmsg

         sendCounts(:) = 0
         myRecv(:) = 0
         do j = 1, hts_ah_nsend(L)
            cell1 = hts_ah_send_proc(j, L)
            if (cell1 < 0 .or. cell1 >= nprocs) then
               localBad = 1
            else
               sendCounts(cell1) = hts_ah_send_n(j, L)
            end if
         end do
         do j = 1, hts_ah_nrecv(L)
            cell1 = hts_ah_recv_proc(j, L)
            if (cell1 < 0 .or. cell1 >= nprocs) then
               localBad = 1
            else
               myRecv(cell1) = hts_ah_recv_n(j, L)
            end if
         end do
         call MPI_Alltoall(sendCounts, 1, MPI_INTEGER, recvExpect, 1, MPI_INTEGER, comm, ierr)
         do cell1 = 0, nprocs - 1
            if (recvExpect(cell1) /= myRecv(cell1)) localBad = 1
         end do

         maxBuf = max(maxBuf, sent, n)
         call mpas_log_write(' level $i acoustic halo send cells $i of $i', intArgs=(/L, sent, fullN/))
         sendOff = sendOff + sent
         recvOff = recvOff + n
      end do

      deallocate(touch, sendCounts, recvExpect, myRecv)

      call MPI_Allreduce(localBad, globalBad, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      if (globalBad /= 0) then
         hts_ah_ready = .false.
         call mpas_log_write(' acoustic halo lists do not match; keeping the full layer-1 groups')
         return
      end if

      hts_ah_nvert = hts_nVert
      allocate(hts_ah_sbuf(max(maxBuf, 1) * hts_nVert))
      allocate(hts_ah_rbuf(max(maxBuf, 1) * hts_nVert))
      hts_ah_ready = .true.
      call mpas_log_write(' acoustic halo lists match; layer-1 rho_pp and rtheta_pp use them')

   end subroutine atm_hts_build_acoustic_halo


   subroutine hts_exch_layer1(domain, field)

      ! Exchange one layer-1 cell field on the active level's compact list.
      ! Every rank enters the calls. The caller has already checked hts_ah_ready.

#ifdef MPAS_USE_MPI_F08
      use mpi_f08, only : MPI_Comm, MPI_Request, MPI_Datatype, MPI_STATUSES_IGNORE, &
                          MPI_REAL, MPI_DOUBLE_PRECISION, MPI_Irecv, MPI_Isend, MPI_Waitall
#else
      use mpi, only : MPI_STATUSES_IGNORE, MPI_REAL, MPI_DOUBLE_PRECISION, &
                      MPI_Irecv, MPI_Isend, MPI_Waitall
#endif

      implicit none

      type (domain_type), intent(in) :: domain
      real (kind=RKIND), intent(inout) :: field(:,:)

#ifdef MPAS_USE_MPI_F08
#ifdef SINGLE_PRECISION
      type (MPI_Datatype), parameter :: MPI_REALKIND = MPI_REAL
#else
      type (MPI_Datatype), parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
#endif
      type (MPI_Comm) :: comm
      type (MPI_Request), allocatable :: sreq(:), rreq(:)
#else
#ifdef SINGLE_PRECISION
      integer, parameter :: MPI_REALKIND = MPI_REAL
#else
      integer, parameter :: MPI_REALKIND = MPI_DOUBLE_PRECISION
#endif
      integer :: comm
      integer, allocatable :: sreq(:), rreq(:)
#endif

      integer :: L, i, j, k, n, iCell, base, rank, ierr
      integer :: nVert

      L = hts_active_level
      nVert = hts_ah_nvert
      if (L < 0 .or. L > hts_max_level_used) return
      if (nVert < 1) return

      call mpas_timer_start('hts_acoustic_halo')

      comm = domain % dminfo % comm
      rank = domain % dminfo % my_proc_id
      allocate(sreq(max(hts_ah_nsend(L), 1)))
      allocate(rreq(max(hts_ah_nrecv(L), 1)))

      base = 0
      do i = 1, hts_ah_nrecv(L)
         n = hts_ah_recv_n(i, L)
         call MPI_Irecv(hts_ah_rbuf(base + 1), n * nVert, MPI_REALKIND, hts_ah_recv_proc(i, L), &
                        hts_ah_recv_proc(i, L), comm, rreq(i), ierr)
         base = base + n * nVert
      end do

      base = 0
      do i = 1, hts_ah_nsend(L)
         n = hts_ah_send_n(i, L)
         do j = 1, n
            iCell = hts_ah_send_cell(hts_ah_send_off(i, L) + j)
            do k = 1, nVert
               hts_ah_sbuf(base + (j - 1) * nVert + k) = field(k, iCell)
            end do
         end do
         call MPI_Isend(hts_ah_sbuf(base + 1), n * nVert, MPI_REALKIND, hts_ah_send_proc(i, L), &
                        rank, comm, sreq(i), ierr)
         base = base + n * nVert
      end do

      if (hts_ah_nrecv(L) > 0) then
         call MPI_Waitall(hts_ah_nrecv(L), rreq, MPI_STATUSES_IGNORE, ierr)
      end if

      base = 0
      do i = 1, hts_ah_nrecv(L)
         n = hts_ah_recv_n(i, L)
         do j = 1, n
            iCell = hts_ah_recv_cell(hts_ah_recv_off(i, L) + j)
            do k = 1, nVert
               field(k, iCell) = hts_ah_rbuf(base + (j - 1) * nVert + k)
            end do
         end do
         base = base + n * nVert
      end do

      if (hts_ah_nsend(L) > 0) then
         call MPI_Waitall(hts_ah_nsend(L), sreq, MPI_STATUSES_IGNORE, ierr)
      end if
      deallocate(sreq, rreq)

      call mpas_timer_stop('hts_acoustic_halo')

   end subroutine hts_exch_layer1


end module mpas_atm_hts
