project(segmentation)
cmake_minimum_required(VERSION 3.0)

find_package(Boost 1.58 COMPONENTS heap)
include_directories(Boost_INCLUDE_DIRS)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
find_package(Qt5Gui)

add_library(${PROJECT_NAME}_core baatzalgorithm.cpp segment.cpp)

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_core Qt5::Gui)

set_property(TARGET ${PROJECT_NAME}_core ${PROJECT_NAME}
             PROPERTY CXX_STANDARD 11)
