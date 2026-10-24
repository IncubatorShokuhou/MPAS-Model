! Index-list loop headers for hierarchical time stepping.
! Each macro opens one loop; the existing "end do" closes it.
! hts_masking false, a single level, or hts_lists_work false keeps the
! original index order. Acoustic and damping lists ignore hts_lists_work.

#define HTS_WORK_CELL_SOLVE_DO \
   hts_nloop = cellSolveEnd - cellSolveStart + 1; \
   if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. cellSolveStart == 1 .and. cellSolveEnd == hts_nCellsSolve) hts_nloop = hts_n_work_own(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. cellSolveStart == 1 .and. cellSolveEnd == hts_nCellsSolve) then; \
         iCell = hts_work_own(hts_ii, hts_active_level); \
      else; \
         iCell = cellSolveStart + hts_ii - 1; \
      end if

#define HTS_WORK_CELL_ALL_DO \
   hts_nloop = cellEnd - cellStart + 1; \
   if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. cellStart == 1 .and. cellEnd == hts_nCells) hts_nloop = hts_n_work_all(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. cellStart == 1 .and. cellEnd == hts_nCells) then; \
         iCell = hts_work_all(hts_ii, hts_active_level); \
      else; \
         iCell = cellStart + hts_ii - 1; \
      end if

#define HTS_ACTIVE_CELL_SOLVE_DO \
   hts_nloop = cellSolveEnd - cellSolveStart + 1; \
   if (hts_masking .and. hts_max_level_used > 0 .and. cellSolveStart == 1 .and. cellSolveEnd == hts_nCellsSolve) hts_nloop = hts_n_active_own(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_masking .and. hts_max_level_used > 0 .and. cellSolveStart == 1 .and. cellSolveEnd == hts_nCellsSolve) then; \
         iCell = hts_active_own(hts_ii, hts_active_level); \
      else; \
         iCell = cellSolveStart + hts_ii - 1; \
      end if

#define HTS_ACTIVE_CELL_ALL_DO \
   hts_nloop = cellEnd - cellStart + 1; \
   if (hts_masking .and. hts_max_level_used > 0 .and. cellStart == 1 .and. cellEnd == hts_nCells) hts_nloop = hts_n_active_all(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_masking .and. hts_max_level_used > 0 .and. cellStart == 1 .and. cellEnd == hts_nCells) then; \
         iCell = hts_active_all(hts_ii, hts_active_level); \
      else; \
         iCell = cellStart + hts_ii - 1; \
      end if

#define HTS_WORK_EDGE_SOLVE_DO \
   hts_nloop = edgeSolveEnd - edgeSolveStart + 1; \
   if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. edgeSolveStart == 1 .and. edgeSolveEnd == hts_nEdgesSolve) hts_nloop = hts_n_work_edge_own(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. edgeSolveStart == 1 .and. edgeSolveEnd == hts_nEdgesSolve) then; \
         iEdge = hts_work_edge_own(hts_ii, hts_active_level); \
      else; \
         iEdge = edgeSolveStart + hts_ii - 1; \
      end if

#define HTS_WORK_EDGE_ALL_DO \
   hts_nloop = edgeEnd - edgeStart + 1; \
   if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. edgeStart == 1 .and. edgeEnd == hts_nEdges) hts_nloop = hts_n_work_edge_all(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. edgeStart == 1 .and. edgeEnd == hts_nEdges) then; \
         iEdge = hts_work_edge_all(hts_ii, hts_active_level); \
      else; \
         iEdge = edgeStart + hts_ii - 1; \
      end if

#define HTS_AC_EDGE_DO \
   hts_nloop = edgeEnd - edgeStart + 1; \
   if (hts_masking .and. hts_max_level_used > 0 .and. edgeStart == 1 .and. edgeEnd == hts_nEdges) hts_nloop = hts_n_ac_edge(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_masking .and. hts_max_level_used > 0 .and. edgeStart == 1 .and. edgeEnd == hts_nEdges) then; \
         iEdge = hts_ac_edge(hts_ii, hts_active_level); \
      else; \
         iEdge = edgeStart + hts_ii - 1; \
      end if

#define HTS_VERT_ALL_DO \
   hts_nloop = vertexEnd - vertexStart + 1; \
   if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. vertexStart == 1 .and. vertexEnd == hts_nVertices) hts_nloop = hts_n_vert_all(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. vertexStart == 1 .and. vertexEnd == hts_nVertices) then; \
         iVertex = hts_vert_all(hts_ii, hts_active_level); \
      else; \
         iVertex = vertexStart + hts_ii - 1; \
      end if

#define HTS_DIAG_CELL_ALL_DO \
   hts_nloop = cellEnd - cellStart + 1; \
   if (hts_masking .and. hts_max_level_used > 0 .and. cellStart == 1 .and. cellEnd == hts_nCells) hts_nloop = hts_n_diag_all(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_masking .and. hts_max_level_used > 0 .and. cellStart == 1 .and. cellEnd == hts_nCells) then; \
         iCell = hts_diag_all(hts_ii, hts_active_level); \
      else; \
         iCell = cellStart + hts_ii - 1; \
      end if

#define HTS_DIAG_EDGE_ALL_DO \
   hts_nloop = edgeEnd - edgeStart + 1; \
   if (hts_masking .and. hts_max_level_used > 0 .and. edgeStart == 1 .and. edgeEnd == hts_nEdges) hts_nloop = hts_n_diag_edge_all(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_masking .and. hts_max_level_used > 0 .and. edgeStart == 1 .and. edgeEnd == hts_nEdges) then; \
         iEdge = hts_diag_edge_all(hts_ii, hts_active_level); \
      else; \
         iEdge = edgeStart + hts_ii - 1; \
      end if

#define HTS_DIAG_VERT_ALL_DO \
   hts_nloop = vertexEnd - vertexStart + 1; \
   if (hts_masking .and. hts_max_level_used > 0 .and. vertexStart == 1 .and. vertexEnd == hts_nVertices) hts_nloop = hts_n_diag_vert_all(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_masking .and. hts_max_level_used > 0 .and. vertexStart == 1 .and. vertexEnd == hts_nVertices) then; \
         iVertex = hts_diag_vert_all(hts_ii, hts_active_level); \
      else; \
         iVertex = vertexStart + hts_ii - 1; \
      end if

#define HTS_VERT_SOLVE_DO \
   hts_nloop = vertexSolveEnd - vertexSolveStart + 1; \
   if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. vertexSolveStart == 1 .and. vertexSolveEnd == hts_nVerticesSolve) hts_nloop = hts_n_vert_own(hts_active_level); \
   hts_ii = 0; \
   do; \
      hts_ii = hts_ii + 1; \
      if (hts_ii > hts_nloop) exit; \
      if (hts_lists_work .and. hts_masking .and. hts_max_level_used > 0 .and. vertexSolveStart == 1 .and. vertexSolveEnd == hts_nVerticesSolve) then; \
         iVertex = hts_vert_own(hts_ii, hts_active_level); \
      else; \
         iVertex = vertexSolveStart + hts_ii - 1; \
      end if
