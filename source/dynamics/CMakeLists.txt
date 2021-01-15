set(LIBRARY_NAME ariadne-dynamics)

add_library(${LIBRARY_NAME} OBJECT
    enclosure.cpp
    orbit.cpp
    inclusion_vector_field.cpp
    inclusion_evolver.cpp
    first_order_pde.cpp
    map_evolver.cpp
    vector_field_evolver.cpp
    simulator.cpp
    reachability_analyser.cpp
)

cotire(${LIBRARY_NAME} CONFIGURATIONS None)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
