#ifndef MES_SOLVER_H
#define MES_SOLVER_H

#include <vector>
#include "Element.h"
#include "JacobianMatrix.h"
#include "Grid.h"
#include "Result.h"
#include <armadillo>

class Solver {
public:
    Solver() = default;

    ~Solver();

    void setParameters(double heatCapacity, double density, double alpha, double ambientTemperature);

    void setBoundaryConditions(bool bottom, bool right, bool top, bool left);

    void setGrid(Grid *grid);

    std::vector<Result> evaluate(const double &time, const double &tau);

protected:
    double **calculateLocalHMatrix(const Element *element, const std::vector<JacobianMatrix> &jacobianMatrices);

    double **calculateLocalCMatrix(const std::vector<JacobianMatrix> &jacobianMatrices);

    double *calculateLocalPVector(const Element *element);

    double **initializeTwoDimensionalArray(const int &n, const int &m);

    double *pointShapeFunctions(const Point &point);

    double *calculateSegmentLengths(const Element *element);

    arma::vec solveEquation(double **H, const double *p);

    Grid *grid{nullptr};
    double heatCapacity{};
    double density{};
    double alpha{};
    double ambientTemperature{};
    bool boundaryConditions[Constants::BOUNDARY_COUNT]{};
};


#endif //MES_SOLVER_H
