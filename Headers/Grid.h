#ifndef MES_GRID_H
#define MES_GRID_H

#include "Node.h"
#include "Element.h"
#include <vector>
#include <string>
#include <ostream>

class Grid {
public:
    Grid(double width, double height, int columns, int rows);

    ~Grid();

    void init(double elementInitialTemperature, double elementKFactor);

    void drawGrid(std::ostream &out);

    std::vector<Element *> getElements() const;

    Node **getNodes() const;

    int getRows() const;

    int getColumns() const;

    static const int GRID_DIMENSIONS = 2;
protected:
    std::vector<Element *> elements;
    Node **nodes{};
    double height;
    double width;
    int rows;
    int columns;
};

#endif //MES_GRID_H
