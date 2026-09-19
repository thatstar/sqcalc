find_package(Python3 COMPONENTS Interpreter)

# Helper that dumps an HDF5 result file as text (used by the python tests).
if(SQC_HAVE_HDF5)
  add_executable(h5read h5read.f90)
  target_link_libraries(h5read PRIVATE HDF5::HDF5)
endif()

# Layout check for the FINUFFT options mirror (fails loudly on an ABI change).
add_executable(test_finufft_opts test_finufft_opts.f90)
target_link_libraries(test_finufft_opts PRIVATE sqcalc_lib)
add_test(NAME finufft_opts COMMAND test_finufft_opts)

# Cell list: pair list versus brute force enumeration.
add_executable(test_neighbour_list test_neighbour_list.f90)
target_link_libraries(test_neighbour_list PRIVATE sqcalc_lib)
add_test(NAME neighbour_list COMMAND test_neighbour_list)

if(Python3_Interpreter_FOUND)
  set(_sqc_h5read "")
  if(SQC_HAVE_HDF5)
    set(_sqc_h5read "--h5read=$<TARGET_FILE:h5read>")
  endif()
  add_test(
    NAME sqcalc_physics
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_sqcalc.py
            --executable $<TARGET_FILE:sqcalc> --workdir ${CMAKE_CURRENT_BINARY_DIR}/work
            ${_sqc_h5read}
  )
  set_tests_properties(sqcalc_physics PROPERTIES LABELS "integration")

  if(SQC_ENABLE_CUDA)
    # Compares the cufinufft backend against the CPU results; skips itself when
    # no CUDA device is visible (so it can also run on CPU-only build hosts).
    add_test(
      NAME sqcalc_gpu
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_sqcalc.py
              --executable $<TARGET_FILE:sqcalc>
              --workdir ${CMAKE_CURRENT_BINARY_DIR}/work_gpu --gpu
              ${_sqc_h5read}
    )
    set_tests_properties(sqcalc_gpu PROPERTIES LABELS "integration;gpu")
  endif()
endif()

# The distributed skill documents the command line flags and the CMake options;
# keep those documents in step with the parser and this build.
set(_sqc_skill_check "${CMAKE_SOURCE_DIR}/skills/sq-calc/scripts/check_doc_sync.py")
if(Python3_Interpreter_FOUND AND EXISTS "${_sqc_skill_check}")
  add_test(
    NAME skills_doc_sync
    COMMAND ${Python3_EXECUTABLE} "${_sqc_skill_check}" "${CMAKE_SOURCE_DIR}"
  )
endif()
