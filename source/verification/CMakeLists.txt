set(LIBRARY_NAME ariadne-verification)

add_library(${LIBRARY_NAME} OBJECT
    verification_manager.cpp
)

cotire(${LIBRARY_NAME} CONFIGURATIONS None)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()