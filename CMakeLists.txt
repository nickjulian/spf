cmake_minimum_required(VERSION 3.1.0)
project(spf_3d.x VERSION 0.1.0 LANGUAGES C CXX)

# set( variableName "variableValue")
# message(STATUS "message to print")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_BUILD_TYPE Debug) # might not work, maybe supposed to send -g to compiler

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # creates compile_commands.json
set(CMAKE_VERBOSE_MAKEFILE TRUE)

message("CMAKE_SYSTEM_NAME IS: ${CMAKE_SYSTEM_NAME}")

set(CMAKE_CXX_FLAGS "-std=gnu++17 -Wall -Wextra -march=native -Ofast -g")

set(SPF_INCLUDE_Dir "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SPF_SRC_Dir "${CMAKE_CURRENT_SOURCE_DIR}/src")
message("PROJECT_NAME: ${PROJECT_NAME}")
message("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")
message("SPF_SRC_Dir: ${SPF_SRC_Dir}")
message("SPF_INCLUDE_Dir: ${SPF_INCLUDE_Dir}")

set(SPF ${PROJECT_NAME})

###############
find_package(MPI REQUIRED COMPONENTS C)
if(NOT MPI_FOUND)
   message(FATAL_ERROR "MPI library not found")
endif()
message(STATUS "MPI_VERSION:" ${MPI_VERSION})
message(STATUS "MPI_C_COMPILE_FLAGS:" ${MPI_C_COMPILE_FLAGS})
message(STATUS "MPI_INCLUDE_PATH:" ${MPI_INCLUDE_PATH})
message(STATUS "MPI_C_HEADER_DIR:" ${MPI_C_HEADER_DIR})
###############
set(HDF5_PREFER_PARALLEL TRUE)
message(STATUS "HDF5_PREFER_PARALLEL: ${HDF5_PREFER_PARALLEL}")
find_package(HDF5 REQUIRED) # TODO: switch to something that will download it if it's not found
if(NOT ${HDF5_FOUND})
   message(FATAL_ERROR "hdf5 library not found")
endif()
message(STATUS "HDF5_DEFINITIONS: ${HDF5_DEFINITIONS}")
message(STATUS "HDF5_C_DEFINITIONS: ${HDF5_C_DEFINITIONS}")
message(STATUS "HDF5_IS_PARALLEL: ${HDF5_IS_PARALLEL}")
if( NOT ${HDF5_IS_PARALLEL})
   message(FATAL_ERROR "HDF5_IS_PARALLEL: ${HDF5_IS_PARALLEL}. This program requires a version of HDF5 with MPI capabilities enabled.")
endif()
message(STATUS "SPF_LINK_LIBRARIES_LIST: ${SPF_LINK_LIBRARIES_LIST}")
###############
list(APPEND SPF_LINK_LIBRARIES_LIST ${MPI_C_LIBRARIES})
list(APPEND SPF_LINK_LIBRARIES_LIST HDF5::HDF5)

message(STATUS "SPF_LINK_LIBRARIES_LIST: ${SPF_LINK_LIBRARIES_LIST}")

add_executable(${SPF} ${SPF_SRC_Dir}/spf_3d.cpp)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/ src)
#
target_compile_definitions(${SPF} PUBLIC ${HDF5_DEFINITIONS})
target_compile_definitions(${SPF} PUBLIC ${HDF5_C_DEFINITIONS})
#
target_include_directories(${SPF} PUBLIC ${SPF_INCLUDE_Dir}/)
target_include_directories(${SPF} PUBLIC ${SPF_INCLUDE_Dir}/Random123/)
target_include_directories(${SPF} PUBLIC ${HDF5_C_INCLUDE_DIRS})
target_include_directories(${SPF} PUBLIC ${HDF5_INCLUDE_DIRS})
target_include_directories(${SPF} PUBLIC ${MPI_C_INCLUDE_DIRS})
#
target_link_libraries(${SPF} PRIVATE ${SPF_LINK_LIBRARIES_LIST})
