cmake_minimum_required(VERSION 3.28)
project(KoMD)

set(CMAKE_CXX_STANDARD 20)

#find_package(OpenMP REQUIRED)

#list(APPEND CMAKE_PREFIX_PATH "~/libs")
#find_package(Kokkos REQUIRED)

# modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)
# include vtk module
include(vtk)
# include kokkos module
include(kokkos)
# include OpenMP module
#include(OpenMP)

file(GLOB_RECURSE files src/*.cpp)
if(NOT ENABLE_VTK)
    list(FILTER files EXCLUDE REGEX ".*/IO/vtk/.*")
endif ()

add_library(KoMD_lib STATIC ${files})
target_include_directories(KoMD_lib PUBLIC src)
target_link_libraries(KoMD_lib
        ${KOKKOS_LIB}           # for Kokkos
)

if(ENABLE_VTK)
    target_link_libraries(KoMD_lib
            ${VTK_LIB}              # for VTK/xerces
    )
    target_compile_definitions(KoMD_lib PUBLIC -DUSE_VTK)
endif ()

add_executable(KoMD "src/main.cpp")
target_include_directories(KoMD PUBLIC src)
target_link_libraries(KoMD PRIVATE KoMD_lib)
target_compile_options(KoMD PRIVATE -Wall -Wextra)
