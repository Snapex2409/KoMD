# KoMD
Small MD framework


The intention is to provide a full MD framework for non-rigid multi-site molecular models. 
Portability is provided by using Kokkos. Currently, OpenMP and Cuda backends are supported.

## Prequisites
* vtk Library: ```sudo apt install libvtk<version>-qt-dev``` minimum version: 7
* For CUDA support: Install CUDA Toolkit, follow the instructions on the official CUDA Toolkit page.
* For OpenMP support: Install OpenMP lib and runtime

## Compilation
* Create build directory ``` mkdir build && cd build``` 
* Configure CMake ``` cmake .. <options>``` or ```ccmake ..```
* CMake options:
  * ENABLE_VTK:BOOL=ON|OFF
  * CUDA options:
    * CUDA_TOOLKIT_ROOT_DIR:FILEPATH=```<path to cuda installation>```
    * CUDAToolkit_INCLUDE_DIRECTORIES:FILEPATH=```<path to cuda installation>/targets/x86_64-linux/include```
    * CUDA_CUDART:FILEPATH=```<path to cuda installation>/targets/x86_64-linux/lib```
  * Kokkos options:
    * Kokkos_ENABLE_```<backend>```:BOOL=ON backend options are CUDA and OpenMP
    * Kokkos_ARCH_NATIVE:BOOL=ON
    * Kokkos_NVCC_WRAPPER:FILEPATH=```<path to KoMD repo>```/dependencies-external/kokkos-4.3.01/bin/nvcc_wrapper
    * Kokkos_ARCH_```<arch>```:BOOL=ON set arch to used GPU architecture
* Compile with ```make all```

## Troubleshooting
* When installing CUDA Toolkit, make sure it is fully installed and sourced. If it is not, then many CUDA tools and libs need to be symlinked into PATH directories.
* Recommendation: Use ccmake to configure CMake. There all options will be visible.

## Input File format
* see ```examples/inp.dat```

## Running KoMD
* ```KoMD <path to input file>```
* When running with OpenMP, set these env variables: OMP_PROC_BIND=spread;OMP_NUM_THREADS=8
