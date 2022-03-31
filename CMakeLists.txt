cmake_minimum_required(VERSION 3.16)
project(TDSE_Solver VERSION 1.0 DESCRIPTION "TDSE Solver" LANGUAGES CXX)


set(INCLUDE_FILES "./cpp/include/")
set(LIBRARY_FILES "./cpp/lib/")
set(WRAPPER_FOLDER "./cpp/wrappers")


add_subdirectory(thirdparty/pybind11)
include_directories("thirdparty/eigen" "thirdparty/spectra/include")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -fopenmp)
    add_compile_options(-O3 -msse4 -DNDEBUG -g -std=c++11 -parallel -mavx -mfma)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-fopenmp)
    add_compile_options(-O3 -msse4 -DNDEBUG -g -std=c++11 -parallel -qmkl -mavx -mfma)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(BLA_VENDOR Intel10_64lp)
    set(MKL_INTERFACE_FULL intel_lp64)
    find_package(MKL)
    find_package(OpenMP REQUIRED)
    find_package(LAPACK REQUIRED)
    add_compile_options(-qopenmp)
    add_compile_options(-O3 -msse4 -DNDEBUG -g -std=c++11 -parallel -qmkl -mavx -mfma -qopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  # using Visual Studio C++
endif()


pybind11_add_module(hamiltonian SHARED 
    ${WRAPPER_FOLDER}/wrapper.cpp 
    ${LIBRARY_FILES}/s32.cpp
   # ${LIBRARY_FILES}/bornoppenheimer.cpp
   # ${LIBRARY_FILES}/hamiltonian2d.cpp 
    #${LIBRARY_FILES}/hamiltonian3d.cpp
)
#target_include_directories(hamiltonian SHARED thirdparty/eigen)
target_include_directories(hamiltonian PUBLIC ${INCLUDE_FILES})


if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    target_link_libraries(hamiltonian PUBLIC -Wl, -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl)
else()
    target_link_libraries(hamiltonian PUBLIC)
endif()

    #target_link_libraries(hamiltonian PUBLIC OpenMP::OpenMP_CXX -lmkl_sequential -lmkl_core -lpthread -lm -ldl)# )