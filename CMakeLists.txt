cmake_minimum_required(VERSION 2.8)
project(Fast_RNRR)

#--- CMake configuration
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#--- Compiler configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

#--- OpenMP (cannot work in clang)
find_package(OpenMP QUIET)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

#--- Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# enable Intel-MKL (only for Intel compilers)
# NB requires correct specification of BLA_VENDOR, e.g.
#   export BLA_VENDOR=Intel10_64lp_seq
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  find_package(BLAS)
  include_directories(${BLAS_INCLUDE_DIRS})
  set(PROJECT_LIBS ${BLAS_LIBRARIES})
  add_definitions(-DEIGEN_USE_MKL_ALL)
endif()

#--- Nanoflann
find_package(NanoFlann REQUIRED)
include_directories(${NANOFLANN_INCLUDE_DIR})

#--- Build example
include_directories("./src/")
file(GLOB SOURCES "./src/*.cpp" "./src/geodesic/*.cpp" "./src/tools/*.cpp")
file(GLOB HEADERS "./src/*.h" "./src/geodesic/*.h" "./src/tools/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

#--- openmesh libs
target_link_libraries(${PROJECT_NAME}
    OpenMeshCore OpenMeshTools)


if(APPLE OR UNIX)
    #--- Deploy data folder link
    execute_process(COMMAND ln -f -s ${CMAKE_SOURCE_DIR}/data WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    #--- "Dummy" target to have data appear in QtCreator
    add_custom_target(data SOURCES ${CMAKE_SOURCE_DIR}/data)
else()
    file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})
endif()
