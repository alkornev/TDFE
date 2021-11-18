cmake_minimum_required(VERSION 3.16)
project(TDSE_Solver VERSION 1.0 DESCRIPTION "TDSE Solver" LANGUAGES CXX)


set(INCLUDE_FILES "./include/")
set(LIBRARY_FILES "./lib/")
set(WRAPPER_FOLDER "./wrappers")

#include_directories(include/basis/)
#add_library(S32 ${LIBRARY_FILES}/S32.cpp)
#target_include_directories(S32 PUBLIC ${INCLUDE_FILES})

#include_directories(${INCLUDE_FILES})
set(BLA_VENDOR Intel10_64lp)
set(MKL_INTERFACE_FULL intel_lp64)
find_package(MKL)
find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)


add_subdirectory(thirdparty/pybind11)
include_directories(thirdparty/eigen)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_compile_options(-qopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
endif()
add_compile_options(-O3 -qmkl -mavx -mfma -qopenmp)


pybind11_add_module(hamiltonian SHARED ${LIBRARY_FILES}/s32.cpp ${LIBRARY_FILES}/bornoppenheimer.cpp ${LIBRARY_FILES}/hamiltonian2d.cpp ${LIBRARY_FILES}/hamiltonian3d.cpp ${WRAPPER_FOLDER}/wrapper.cpp)
target_include_directories(hamiltonian PUBLIC ${INCLUDE_FILES} ${MKLROOT}/include)
target_link_libraries(hamiltonian PUBLIC OpenMP::OpenMP_CXX -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
#target_link_libraries(hamiltonian PUBLIC OpenMP::OpenMP_CXX -lmkl_sequential -lmkl_core -lpthread -lm -ldl)# )