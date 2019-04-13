#include <iostream>
#include "Grid.h"

int main() {
    Grid grid("indata/test");
//    grid.drawGrid(std::cout);
    grid.evaluate(500, 50);
}