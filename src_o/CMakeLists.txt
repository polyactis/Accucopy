cmake_minimum_required(VERSION 3.0.0)
project(infer_test VERSION 0.1.0)

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

add_executable(infer BaseGADA.cc read_para.cpp format.cc infer.cpp prob.cpp)
add_executable(GADA BaseGADA.cc read_para.cpp format.cc GADA.cc)
if(Boost_FOUND)
    target_link_libraries(infer ${Boost_LIBRARIES})
    target_link_libraries(GADA ${Boost_LIBRARIES})
endif()

find_package(GSL REQUIRED)
target_link_libraries(infer GSL::gsl GSL::gslcblas)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
