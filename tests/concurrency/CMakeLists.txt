set(UNIT_TESTS
    test_task_search_parameter
    test_task_ranking_parameter
    test_buffer
    test_concurrency_manager
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    if(COVERAGE)
        target_compile_options(${TEST} PUBLIC ${COVERAGE_COMPILER_FLAGS})
    endif()    
    target_link_libraries(${TEST} ariadne)
    add_test(${TEST} ${TEST})    
    set_property(TEST ${TEST} PROPERTY LABELS "concurrency;primary")
endforeach()

add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
