cmake_minimum_required(VERSION 3.0.0)
project(infer VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_STANDARD 11)

message(STATUS "start running cmake...")
find_package(Boost 1.58.0 COMPONENTS iostreams program_options REQUIRED)
if(Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_executable(infer src_o/BaseGADA.cc src_o/read_para.cpp src_o/format.cc src_o/infer.cpp src_o/prob.cpp src_o/model_selection.cpp)
add_executable(GADA src_o/BaseGADA.cc src_o/read_para.cpp src_o/format.cc src_o/GADA.cc)
if(Boost_FOUND)
    target_link_libraries(infer ${Boost_LIBRARIES})
    target_link_libraries(GADA ${Boost_LIBRARIES})
endif()

find_package(GSL REQUIRED)
target_link_libraries(infer GSL::gsl GSL::gslcblas)

find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(infer  ${HDF5_CXX_LIBRARIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
