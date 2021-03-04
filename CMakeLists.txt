cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(virtual_scanner)

find_package(PCL REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenCV REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(${OpenCV_LIBRARY_DIRS})
add_definitions(${OpenCV_DEFINITIONS})

add_executable (virtual_scanner virtual_scanner.cpp)
target_link_libraries (virtual_scanner ${PCL_LIBRARIES} ${Boost_LIBRARIES} ${OpenCV_LIBS})
