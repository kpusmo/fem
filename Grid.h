#ifndef MES_GRID_H
#define MES_GRID_H

#include "Node.h"
#include "Element.h"
#include "JacobianMatrix.h"
#include <vector>
#include <string>
#include <ostream>

class Grid {
public:
    explicit Grid(std::string filename);

    ~Grid();

    void drawGrid(std::ostream &out);

    void evaluate(const double &time, const double &tau);

    static const int GRID_DIMENSIONS = 2;
protected:
    void init(double **temperatures, double **ks);

    double **calculateLocalHMatrix(const Element *element, const std::vector<JacobianMatrix> &jacobianMatrices);

    double **calculateLocalCMatrix(const std::vector<JacobianMatrix> &jacobianMatrices);

    double *calculateLocalPVector(const Element *element);

    double **initializeTwoDimensionalArray(const int &n, const int &m);

    template<typename T>
    void deleteTwoDimensionalArray(T **array, const int &rows);

    double *pointShapeFunctions(const Point &point);

    double *calculateSegmentLengths(const Element *element);

    double *solveEquation(double **H, double *p);

    std::vector<Element *> elements;
    Node **nodes{};
    double height;
    double width;
    int rows;
    int columns;
    double heatCapacity;
    double density;
    double alpha;
    double ambientTemperature;
    bool boundaryConditions[Constants::BOUNDARY_COUNT];
};

#endif //MES_GRID_H
