set(UNIT_TESTS
    test_expression
)

foreach(TEST ${UNIT_TESTS})
    add_executable(${TEST} ${TEST}.cpp)
    if(COVERAGE)
        target_compile_options(${TEST} PUBLIC ${COVERAGE_COMPILER_FLAGS})
    endif()    
    target_link_libraries(${TEST} ariadne-kernel)  
    add_test(${TEST} ${TEST})
    set_property(TEST ${TEST} PROPERTY LABELS "symbolic;primary;kernel")           
endforeach()

add_dependencies(tests-kernel ${UNIT_TESTS})
add_dependencies(tests-cpp ${UNIT_TESTS})
add_dependencies(tests ${UNIT_TESTS})
