cmake_minimum_required(VERSION 3.28)
project(KoMD)

set(CMAKE_CXX_STANDARD 23)

find_package(OpenMP REQUIRED)
#list(APPEND CMAKE_PREFIX_PATH "~/libs")
#find_package(Kokkos REQUIRED)

file(GLOB_RECURSE files src/*.cpp)

add_library(KoMD_lib STATIC ${files})
target_include_directories(KoMD_lib PUBLIC src)
target_link_libraries(KoMD_lib
        PRIVATE OpenMP::OpenMP_CXX
#        PUBLIC Kokkos::kokkos
)


add_executable(KoMD "src/main.cpp")
target_include_directories(KoMD PUBLIC src)
target_link_libraries(KoMD PRIVATE KoMD_lib)
target_compile_options(KoMD PRIVATE -Wall -Wextra)
