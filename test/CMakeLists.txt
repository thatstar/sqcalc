find_package(Python3 COMPONENTS Interpreter)

# Layout check for the FINUFFT options mirror (fails loudly on an ABI change).
add_executable(test_finufft_opts test_finufft_opts.f90)
target_link_libraries(test_finufft_opts PRIVATE sqcalc_lib)
add_test(NAME finufft_opts COMMAND test_finufft_opts)

if(Python3_Interpreter_FOUND)
  add_test(
    NAME sqcalc_physics
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_sqcalc.py
            --executable $<TARGET_FILE:sqcalc> --workdir ${CMAKE_CURRENT_BINARY_DIR}/work
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
    )
    set_tests_properties(sqcalc_gpu PROPERTIES LABELS "integration;gpu")
  endif()
endif()
