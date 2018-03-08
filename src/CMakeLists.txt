if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

cmake_minimum_required (VERSION 2.6)
project (CFDV)

add_subdirectory(adap)
# add_subdirectory(chem)
add_subdirectory(Utility)
add_subdirectory(FiniteElement)
add_subdirectory(FDV)
# add_subdirectory(gridgen)

find_package(MPI REQUIRED)
find_package(Boost 1.5.8 COMPONENTS program_options REQUIRED)

include_directories(${MPI_INCLUDE_PATH}
${Boost_INCLUDE_DIR}
"${PROJECT_SOURCE_DIR}/adap"
"${PROJECT_SOURCE_DIR}/Utility"
"${CMAKE_SOURCE_DIR}/Utility/IO"
"${CMAKE_SOURCE_DIR}/Utility/Solvers"
"${PROJECT_SOURCE_DIR}/FiniteElement"
"${PROJECT_SOURCE_DIR}/FDV"
"${CMAKE_SOURCE_DIR}/FDV/Function")

add_executable(CFDV main.cpp)
add_executable(CONVERT ConvertMain.cpp)
add_executable(EXTRACT ExtractMain.cpp)

target_link_libraries(CFDV ${MPI_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(CONVERT ${MPI_LIBRARIES} ${Boost_LIBRARIES})
target_link_libraries(EXTRACT ${MPI_LIBRARIES} ${Boost_LIBRARIES})


if(MPI_COMPILE_FLAGS)
  set_target_properties(CFDV PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(CFDV PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# TODO clean up code to compile without -fpermissive
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
