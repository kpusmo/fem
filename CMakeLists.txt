cmake_minimum_required(VERSION 3.12)
project(mes)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOUIC_SEARCH_PATHS  UI)

find_package(Qt5Widgets REQUIRED)

set(CMAKE_CXX_STANDARD 17)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS} Headers)

add_executable(mes main.cpp Headers/Node.h Headers/Element.h Headers/Grid.h Sources/Grid.cpp Headers/Point.h Sources/JacobianMatrix.cpp Headers/JacobianMatrix.h Headers/Constants.h Sources/Element.cpp Sources/MainWindow.cpp Headers/MainWindow.h Sources/Solver.cpp Headers/Solver.h Headers/Helpers.h Sources/Result.cpp Headers/Result.h Sources/ResultsModel.cpp Headers/ResultsModel.h)
target_link_libraries(mes ${ARMADILLO_LIBRARIES} Qt5::Widgets)