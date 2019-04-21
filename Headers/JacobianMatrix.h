#ifndef MES_JACOBIANMATRIX_H
#define MES_JACOBIANMATRIX_H


#include "Element.h"

class JacobianMatrix {
public:
    JacobianMatrix(const Element *e, const int &integrationPointIndex);

    ~JacobianMatrix();

    const double *operator[](const int &i) const {
        return matrix[i];
    }

    double **getMatrix() const {
        return matrix;
    }

    double **getInversedMatrix() const {
        return inversedMatrix;
    }

    double getDeterminant() const {
        return determinant;
    }

protected:
    void evaluate(const int &integrationPointIndex);

    void evaluateMatrix(const int &integrationPointIndex);

    void evaluateDeterminant();

    void evaluateInversedMatrix();

    double **matrix;
    double **inversedMatrix;
    double determinant;
    const Element *element;
};


#endif //MES_JACOBIANMATRIX_H
