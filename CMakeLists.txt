cmake_minimum_required(VERSION 3.12)
project(mes)

set(CMAKE_CXX_STANDARD 17)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

add_executable(mes main.cpp Node.h Element.h Grid.h Grid.cpp Point.h JacobianMatrix.cpp JacobianMatrix.h Constants.h)
target_link_libraries(mes ${ARMADILLO_LIBRARIES})