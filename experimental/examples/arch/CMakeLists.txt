set(EXAMPLES_EXPERIMENTAL_ARCH
    vanderpol_arch
    vanderpol_2coupled_arch
	laub_loomis_arch
	quadrotor_arch
    lotka_volterra_arch
    production_destruction_arch
    space_rendezvous_arch
)

foreach(example ${EXAMPLES_EXPERIMENTAL_ARCH})
    add_executable(${example} ${example}.cpp)
    target_link_libraries(${example} ariadne)
endforeach()