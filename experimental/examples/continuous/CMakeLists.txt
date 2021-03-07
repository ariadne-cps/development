set(EXAMPLES_EXPERIMENTAL_CONTINUOUS
    harmonic_oscillator
    laub-loomis-nonoise
	lotka-volterra-nonoise
	higgins-selkov-nonoise
)

foreach(example ${EXAMPLES_EXPERIMENTAL_CONTINUOUS})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()

add_subdirectory(noisy)