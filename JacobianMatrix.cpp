#include "JacobianMatrix.h"
#include "Grid.h"

JacobianMatrix::JacobianMatrix(const Element *e, const int &integrationPointIndex) : element(e) {
    matrix = new double *[Grid::GRID_DIMENSIONS];
    inversedMatrix = new double *[Grid::GRID_DIMENSIONS];
    for (int i = 0; i < Grid::GRID_DIMENSIONS; ++i) {
        matrix[i] = new double[Grid::GRID_DIMENSIONS];
        inversedMatrix[i] = new double[Grid::GRID_DIMENSIONS];
    }
    evaluate(integrationPointIndex);
}

JacobianMatrix::~JacobianMatrix() {
    for (int i = 0; i < Grid::GRID_DIMENSIONS; ++i) {
        delete[] matrix[i];
        delete[] inversedMatrix[i];
    }
    delete[] matrix;
    delete[] inversedMatrix;
}

void JacobianMatrix::evaluate(const int &integrationPointIndex) {
    evaluateMatrix(integrationPointIndex);
    evaluateDeterminant();
    evaluateInversedMatrix();
}

void JacobianMatrix::evaluateMatrix(const int &integrationPointIndex) {
    double dXdKsi = 0;
    double dYdEta = 0;
    double dXdEta = 0;
    double dYdKsi = 0;
    for (int i = 0; i < Constants::ELEMENT_NODE_COUNT; ++i) {
        dXdKsi += Constants::dNdKsi[integrationPointIndex][i] * element->getNode(i)->getX();
        dYdEta += Constants::dNdEta[integrationPointIndex][i] * element->getNode(i)->getY();
        dXdEta += Constants::dNdEta[integrationPointIndex][i] * element->getNode(i)->getX();
        dYdKsi += Constants::dNdKsi[integrationPointIndex][i] * element->getNode(i)->getY();
    }
    matrix[0][0] = dYdEta;
    matrix[0][1] = -dYdKsi;
    matrix[1][0] = -dXdEta;
    matrix[1][1] = dXdKsi;
}

void JacobianMatrix::evaluateDeterminant() {
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}

void JacobianMatrix::evaluateInversedMatrix() {
    inversedMatrix[0][0] = matrix[1][1] / determinant;
    inversedMatrix[1][1] = matrix[0][0] / determinant;
    inversedMatrix[0][1] = -matrix[0][1] / determinant;
    inversedMatrix[1][0] = -matrix[1][0] / determinant;
}
