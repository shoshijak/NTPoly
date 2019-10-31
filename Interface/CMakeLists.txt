################################################################################
add_executable(localdriv local.f90)
target_link_libraries(localdriv NTPoly ${MPI_Fortran_LIBRARIES}
                      ${OpenMP_Fortran_LIBRARIES} ${BLAS_LIBRARIES}
                      ${TOOLCHAIN_LIBS} ${DBCSR_LIBS})
target_include_directories(localdriv PRIVATE ${MPI_INCLUDE_PATH}
                           ${DBCSR_INCLUDE_PATH})
install(TARGETS localdriv DESTINATION bin)

configure_file(input.mtx input.mtx COPYONLY)
configure_file(output.mtx output.mtx COPYONLY)
configure_file(blocks.inp blocks.inp COPYONLY)
