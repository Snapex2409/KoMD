# no option here - we must use KOKKOS
add_subdirectory(${CMAKE_SOURCE_DIR}/dependencies-external/kokkos-4.3.01)

include_directories(${Kokkos_INCLUDE_DIRS_RET})

set(KOKKOS_LIB kokkos)