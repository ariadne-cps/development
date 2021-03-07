set(LIBRARY_NAME ariadne-configuration)

add_library(${LIBRARY_NAME} OBJECT
        searchable_configuration.cpp
        configuration_search_space.cpp
        configuration_property_path.cpp
        configuration_search_parameter.cpp
        configuration_search_point.cpp
        configuration_property.cpp
        )

cotire(${LIBRARY_NAME} CONFIGURATIONS None)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
