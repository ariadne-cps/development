set(LIBRARY_NAME ariadne-numeric)

add_library(${LIBRARY_NAME} OBJECT
    dyadic.cpp
    logical.cpp
    operators.cpp
    integer.cpp
    rational.cpp
    decimal.cpp
    number.cpp
    real.cpp
    complex.cpp
    floatdp.cpp
    floatmp.cpp
    float_approximation.cpp
    float_lower_bound.cpp
    float_upper_bound.cpp
    float_bounds.cpp
    float_ball.cpp
    float_value.cpp
    float_error.cpp
    float_literals.cpp

    concepts.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
     target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
     target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()

