set(LIBRARY_NAME ariadne-function)

add_library(${LIBRARY_NAME} OBJECT
    affine.cpp
    polynomial.cpp
    chebyshev_polynomial.cpp
    formula.cpp
    procedure.cpp
    function.cpp
    measurable_function.cpp
    multifunction.cpp
    function_patch.cpp
    function_model.cpp
    affine_model.cpp
    taylor_model.cpp
    taylor_function.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
